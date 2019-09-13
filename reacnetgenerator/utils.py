"""Provide utils for ReacNetGenerator."""


import lz4.frame
import pybase64
import pickle


class WriteBuffer:
    def __init__(self, f, linenumber=1200, sep=None):
        self.f = f
        if sep is not None:
            self.sep = sep
        elif f.mode=='w':
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
            self.f.write(self.sep.join(self.buff))
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
        yield item, parameter

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

def setparam(x, default):
    return x if x is not None else default

def listtobytes(x):
    return compress(pickle.dumps(x), isbytes=True)

def bytestolist(x):
    return pickle.loads(decompress(x, isbytes=True))

def listtostirng(l, sep):
    if isinstance(l, str):
        return l
    if isinstance(l, (int, float)):
        return str(l)
    if isinstance(l, (list, tuple)):
        return sep[0].join(map(lambda x: listtostirng(x, sep[1:]), l))
    raise RuntimeError("Unsupported type")


class SCOUROPTIONS:
    strip_xml_prolog = True
    remove_titles= True
    remove_descriptions= True
    remove_metadata= True
    remove_descriptive_elements= True
    strip_comments= True
    enable_viewboxing= True
    strip_xml_space_attribute= True
    strip_ids= True
    shorten_ids= True
    newlines= False
