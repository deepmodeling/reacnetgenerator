# cython: language_level=3
# cython: linetrace=True
"""Provide utils for ReacNetGenerator."""


import os
import shutil
import itertools
import logging
import pickle
import hashlib
import asyncio
from multiprocessing import Pool, Semaphore

import lz4.frame
import numpy as np
import pybase64
import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm


class WriteBuffer:
    """Store a buffer for writing files.

    It is expensive to write to a file, so we need to make a buffer.
    
    Parameters
    ----------
    f: fileObject
        The file object to write.
    linenumber: int, default: 1200
        The number of contents to store in the buffer. The buffer will be flushed
        if it exceeds the set number.
    sep: str or bytes, default: None
        The separator for contents. If None (default), there will be no separator.
    """
    def __init__(self, f, linenumber=1200, sep=None):
        self.f = f
        if sep is not None:
            self.sep = sep
        elif f.mode == 'w':
            self.sep = ''
        elif f.mode == 'wb':
            self.sep = b''
        else:
            raise RuntimeError("File mode should be w or wb!")
        self.linenumber = linenumber
        self.buff = []
        self.name = self.f.name

    def append(self, text):
        """Append a text.
        
        Parameters
        ----------
        text: str
            The text to be appended.
        """
        self.buff.append(text)
        self.check()

    def extend(self, text):
        """Extend texts.

        Paramenters
        -----------
        text: list of strs
            Texts to be extended.
        """
        self.buff.extend(text)
        self.check()

    def check(self):
        """Check if the number of stored contents exceeds.
        
        If so, the buffer will be flushed.
        """
        if len(self.buff) > self.linenumber:
            self.flush()

    def flush(self):
        """Flush the buffer."""
        if self.buff:
            self.f.writelines([self.sep.join(self.buff), self.sep])
            self.buff[:] = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.flush()
        self.f.__exit__(exc_type, exc_value, traceback)


def appendIfNotNone(f, wbytes):
    if wbytes is not None:
        f.append(wbytes)


def produce(semaphore, plist, parameter):
    """Item producer with a semaphore.
    
    Prevent large memory usage due to slow IO.
    
    Parameters
    ----------
    semaphore: multiprocessing.Semaphore
        The semaphore to acquire.
    plist: list of objects
        The list of items to be passed.
    parameter: object
        The parameter yielded with each item.

    Yield
    -----
    item: object
        The item to be yielded.
    parameter: object
        The parameter yielded with each item.
    """
    for item in plist:
        semaphore.acquire()
        if parameter is not None:
            item = (item, parameter)
        yield item


def compress(x, isbytes=False):
    """Compress the line.

    This function reduces IO overhead to speed up the program. The functions will
    use lz4 to compress and base64 to encode, since lz4 has better performance
    that any others.

    Parameters
    ----------
    x: str or bytes
        The line to compress.
    isbytes: bool, optional, default: False
        If `x` is bytes. If not, `x` will be converted to bytes first.
    
    Returns
    -------
    bytes
        The compressed line, with a linebreak in the end.
    """
    if isbytes:
        return pybase64.b64encode(lz4.frame.compress(x, compression_level=0))+b'\n'
    return pybase64.b64encode(lz4.frame.compress(x.encode(), compression_level=-1))+b'\n'


def decompress(x, isbytes=False):
    """Decompress the line.
    
    Parameters
    ----------
    x: bytes
        The line to decompress.
    isbytes: bool, optional, default: False
        If the decompressed content is bytes. If not, the line will be decoded.
    
    Returns
    -------
    str or bytes
        The decompressed line.
    """
    if isbytes:
        return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True))
    return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True)).decode()


def listtobytes(x):
    """Convert an object to a compressed line.
    
    Parameters
    ----------
    x: object
        The object to convert, such as numpy.ndarray.
    
    Returns
    -------
    bytes
        The compressed line.
    """
    return compress(pickle.dumps(x), isbytes=True)


def bytestolist(x):
    """Convert a compressed line to an object.
    
    Parameters
    ----------
    x: bytes
        The compressed line.
    
    Returns
    -------
    object
        The decompressed object.
    """
    return pickle.loads(decompress(x, isbytes=True))


def listtostirng(l, sep):
    """Convert a list to string, that is easier to store.

    Parameters
    ----------
    l: list of strs or lists
        The list to convert, which can contain any number of dimensions.
    sep: list of strs
        The seperators for each dimension.

    Returns
    -------
    str
        The converted string.
    """
    if isinstance(l, str):
        return l
    if isinstance(l, (list, tuple, np.ndarray)):
        return sep[0].join(map(lambda x: listtostirng(x, sep[1:]), l))
    return str(l)


def multiopen(pool, func, l, semaphore=None, nlines=None, unordered=True, return_num=False, start=0, extra=None, interval=None, bar=True, desc=None, unit="it", total=None):
    """Returns an interated object for process a file with multiple processors.

    Parameters
    ----------
    pool: multiprocessing.Pool
        The pool for multiprocessing.
    func: function
        The function to process lines.
    l: File object
        The file object.
    semaphore: multiprocessing.Semaphore, optional, default: None
        The semaphore to acquire. If None (default), the object will be passed
        without control.
    nlines: int, optional, default: None
        The number of lines to pass to the function each time. If None (default),
        only one line will be passed to the function.
    unordered: bool, optional, default: True
        Whether the process can be unordered.
    return_num: bool, optional, default: False
        If True, adds a counter to an iterable.
    start: int, optional, default: 0
        The start number of the counter.
    extra: object, optional, default: None
        The extra object passed to the item.
    interval: obj, optional, default: None
        The interval of items that will be passed to the function. For example,
        if set to 10, a item will be passed once every 10 items and others will
        be dropped.
    bar: bool, optional, default: True
        If True, show a tqdm bar for the iteration.
    desc: str, optional, default: None
        The description of the iteration shown in the bar.
    unit: str, optional, default: it
        The unit of the iteration shown in the bar.
    total: int, optional, default: None
        The total number of the iteration shown in the bar.
    
    Returns
    -------
    object
        An object that can be iterated.
    """
    obj = l
    if nlines:
        obj = itertools.zip_longest(*[obj] * nlines)
    if interval:
        obj = itertools.islice(obj, 0, None, interval)
    if return_num:
        obj = enumerate(obj, start)
    if semaphore:
        obj = produce(semaphore, obj, extra)
    if unordered:
        obj = pool.imap_unordered(func, obj, 100)
    else:
        obj = pool.imap(func, obj, 100)
    if bar:
        obj = tqdm(obj, desc=desc, unit=unit, total=total)
    return obj


class SCOUROPTIONS:
    """Scour (SVG optimization) options."""
    strip_xml_prolog = True
    remove_titles = True
    remove_descriptions = True
    remove_metadata = True
    remove_descriptive_elements = True
    strip_comments = True
    enable_viewboxing = True
    strip_xml_space_attribute = True
    strip_ids = True
    shorten_ids = True
    newlines = False


class SharedRNGData:
    """Share ReacNetGenerator data with a submodule.

    Parameters
    ----------
    rng: reacnetgenerator.ReacNetGenerator
        The centered ReacNetGenerator class.
    usedRNGKeys: list of strs
        Keys that needs to pass from ReacNetGenerator class to the submodule.
    returnedRNGKeys: list of strs
        Keys that needs to pass from the submodule to ReacNetGenerator class.
    extraNoneKeys: list of strs, optional, default: None
        Set keys to None, which will be used in the submodule.
    """
    def __init__(self, rng, usedRNGKeys, returnedRNGKeys, extraNoneKeys = None):
        self.rng = rng
        self.returnedRNGKeys = returnedRNGKeys
        for key in usedRNGKeys:
            setattr(self, key, getattr(self.rng, key))
        for key in returnedRNGKeys:
            setattr(self, key, None)
        if extraNoneKeys is not None:
            for key in extraNoneKeys:
                setattr(self, key, None)

    def returnkeys(self):
        """Return back keys to ReacNetGenerator class."""
        for key in self.returnedRNGKeys:
            setattr(self.rng, key, getattr(self, key))


def checksha256(filename, sha256_check):
    """Check sha256 of a file is correct.
    
    Parameters
    ----------
    filename: str
        The filename.
    sha256_check: str or list of strs
        The sha256 to be checked.
    
    Returns
    -------
    bool
        Indicate whether sha256 is correct.
    """
    if not os.path.isfile(filename):
        return
    h = hashlib.sha256()
    b = bytearray(128*1024)
    mv = memoryview(b)
    with open(filename, 'rb', buffering=0) as f:
        for n in iter(lambda: f.readinto(mv), 0):
            h.update(mv[:n])
    sha256 = h.hexdigest()
    logging.info(f"SHA256 of {filename}: {sha256}")
    if sha256 in must_be_list(sha256_check):
        return True
    logging.warning("SHA256 is not correct.")
    logging.warning(open(filename).read())
    return False


async def download_file(urls, pathfilename, sha256):
    """Download files from remote urls if not exists.
    
    Parameters
    ----------
    urls: str or list of strs
        The url(s) that is available to download.
    pathfilename: str
        The downloading path of the file.
    sha256: str
        Sha256 of the file. If not None and match the file, the download will be skiped.
    
    Returns
    -------
    pathfilename: str
        The downloading path of the file.
    """
    s = requests.Session()
    s.mount('http://', HTTPAdapter(max_retries=3))
    s.mount('https://', HTTPAdapter(max_retries=3))
    # download if not exists
    if os.path.isfile(pathfilename) and (sha256 is None or checksha256(pathfilename, sha256)):
        return pathfilename

    # from https://stackoverflow.com/questions/16694907
    for url in must_be_list(urls):
        logging.info(f"Try to download {pathfilename} from {url}")
        with s.get(url, stream=True) as r, open(pathfilename, 'wb') as f:
            try:
                shutil.copyfileobj(r.raw, f)
                break
            except requests.exceptions.RequestException as e:
                logging.warning(f"Request {pathfilename} Error.", exc_info=e)
    else:
        raise RuntimeError(f"Cannot download {pathfilename}.")

    return pathfilename


async def gather_download_files(urls):
    """See download_multifiles function for details."""
    await asyncio.gather(*[download_file(jdata["url"], jdata["fn"], jdata.get("sha256", None)) for jdata in urls])


def download_multifiles(urls):
    """Download multiple files from dicts.

    Parameters
    ----------
    urls: list of dicts
        The information of download files. Each dict should contain the following key:
            - url: str or list of strs
                The url(s) that is available to download.
            - pathfilename: str
                The downloading path of the file.
            - sha256: str, optional, default: None
                Sha256 of the file. If not None and match the file, the download will be skiped.
    """
    asyncio.run(gather_download_files(urls))


def run_mp(nproc, **arg):
    """Process a file with multiple processors.

    Parameters
    ----------
    nproc: int
        The number of processors to be used.
    Other parameters can be found in the `multiopen` function.
    """
    pool = Pool(nproc, maxtasksperchild=1000)
    semaphore = Semaphore(nproc*150)
    try:
        results = multiopen(pool=pool, semaphore=semaphore, **arg)
        for item in results:
            yield item
            semaphore.release()
    except:
        logging.exception("run_mp failed")
        pool.terminate()
        raise
    else:
        pool.close()
    finally:
        pool.join()


def must_be_list(obj):
    """Convert a object to a list if the object is not a list.
    
    Parameters
    ----------
    obj: Object
        The object to convert.

    Returns
    -------
    obj: list
        If the input object is not a list, returns a list that only contains that
        object. Otherwise, returns that object.
    """
    if isinstance(obj, list):
        return obj
    return [obj]
