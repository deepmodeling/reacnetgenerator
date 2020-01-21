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
        self.buff.append(text)
        self.check()

    def extend(self, text):
        self.buff.extend(text)
        self.check()

    def check(self):
        if len(self.buff) > self.linenumber:
            self.flush()

    def flush(self):
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
    """Prevent large memory usage due to slow IO."""
    for item in plist:
        semaphore.acquire()
        if parameter is not None:
            item = (item, parameter)
        yield item


def compress(x, isbytes=False):
    """Compress the line.

    This function reduces IO overhead to speed up the program.
    """
    if isbytes:
        return pybase64.b64encode(lz4.frame.compress(x, compression_level=0))+b'\n'
    return pybase64.b64encode(lz4.frame.compress(x.encode(), compression_level=-1))+b'\n'


def decompress(x, isbytes=False):
    """Decompress the line."""
    if isbytes:
        return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True))
    return lz4.frame.decompress(pybase64.b64decode(x.strip(), validate=True)).decode()


def listtobytes(x):
    return compress(pickle.dumps(x), isbytes=True)


def bytestolist(x):
    return pickle.loads(decompress(x, isbytes=True))


def listtostirng(l, sep):
    if isinstance(l, str):
        return l
    if isinstance(l, (list, tuple, np.ndarray)):
        return sep[0].join(map(lambda x: listtostirng(x, sep[1:]), l))
    return str(l)


def multiopen(pool, func, l, semaphore=None, nlines=None, unordered=True, return_num=False, start=0, extra=None, interval=None, bar=True, desc=None, unit="it", total=None):
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
        for key in self.returnedRNGKeys:
            setattr(self.rng, key, getattr(self, key))


def checksha256(filename, sha256_check):
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
    s = requests.Session()
    s.mount('http://', HTTPAdapter(max_retries=3))
    s.mount('https://', HTTPAdapter(max_retries=3))
    # download if not exists
    if os.path.isfile(pathfilename) and (checksha256(pathfilename, sha256) or sha256 is None):
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
    await asyncio.gather(*[download_file(jdata["url"], jdata["fn"], jdata.get("sha256", None)) for jdata in urls])


def download_multifiles(urls):
    asyncio.run(gather_download_files(urls))


def run_mp(nproc, **arg):
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
    if isinstance(obj, list):
        return obj
    return [obj]
