# SPDX-License-Identifier: LGPL-3.0-or-later
# cython: language_level=3
# cython: linetrace=True
"""Provide utils for ReacNetGenerator."""

import asyncio
import hashlib
import itertools
import os
import pickle
import shutil
import tempfile
import time
from collections.abc import Callable, Generator, Iterable
from contextlib import ExitStack
from functools import partial
from multiprocessing import Event, Pool, Semaphore, SimpleQueue
from multiprocessing.pool import ExceptionWithTraceback  # type: ignore[attr-defined]
from multiprocessing.reduction import ForkingPickler
from queue import Empty, Queue
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    AnyStr,
    BinaryIO,
    Generic,
    Optional,
    cast,
    overload,
)

import lz4.frame
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from tqdm.auto import tqdm

from ._logging import logger

if TYPE_CHECKING:
    import multiprocessing.pool
    import multiprocessing.synchronize

    import reacnetgenerator


class WriteBuffer(Generic[AnyStr]):
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

    def __init__(
        self, f: IO[AnyStr], linenumber: int = 1200, sep: AnyStr | None = None
    ) -> None:
        self.f = f
        if sep is not None:
            self.sep = sep
        elif f.mode == "w":
            self.sep = cast(AnyStr, "")
        elif f.mode == "wb":
            self.sep = cast(AnyStr, b"")
        else:
            raise RuntimeError("File mode should be w or wb!")
        self.linenumber = linenumber
        self.buff: list[AnyStr] = []
        self.name = self.f.name

    def append(self, text: AnyStr) -> None:
        """Append a text.

        Parameters
        ----------
        text : str or bytes
            The text to be appended.
        """
        self.buff.append(text)
        self.check()

    def extend(self, text: Iterable[AnyStr]) -> None:
        """Extend texts.

        Parameters
        ----------
        text : list of strs or bytes
            Texts to be extended.
        """
        self.buff.extend(text)
        self.check()

    def check(self) -> None:
        """Check if the number of stored contents exceeds.

        If so, the buffer will be flushed.
        """
        if len(self.buff) > self.linenumber:
            self.flush()

    def flush(self) -> None:
        """Flush the buffer."""
        if self.buff:
            self.f.writelines([cast(Any, self.sep).join(self.buff), self.sep])
            self.buff[:] = []

    def __enter__(self) -> "WriteBuffer[AnyStr]":
        """Enter the context."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Exit the context."""
        self.flush()
        self.f.__exit__(exc_type, exc_value, traceback)


@overload
def appendIfNotNone(f: WriteBuffer[str] | ExitStack, wbytes: str | None) -> None: ...
@overload
def appendIfNotNone(
    f: WriteBuffer[bytes] | ExitStack, wbytes: bytes | None
) -> None: ...
def appendIfNotNone(f: WriteBuffer[AnyStr] | ExitStack, wbytes: AnyStr | None) -> None:
    """Append a line to a file if the line is not None.

    Parameters
    ----------
    f : WriteBuffer
        The file to write.
    wbytes : str or bytes
        The line to write.
    """
    if wbytes is not None:
        assert not isinstance(f, ExitStack)
        f.append(wbytes)


def produce(
    semaphore: "multiprocessing.synchronize.Semaphore",
    plist: Iterable[Any],
    parameter: Any,
) -> Generator[tuple[Any, Any], None, None]:
    """Item producer with a semaphore.

    Prevent large memory usage due to slow IO.

    Parameters
    ----------
    semaphore : multiprocessing.Semaphore
        The semaphore to acquire.
    plist : list of objects
        The list of items to be passed.
    parameter : object
        The parameter yielded with each item.

    Yields
    ------
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


def compress(x: str | bytes) -> bytes:
    """Compress the line.

    This function reduces IO overhead to speed up the program. The functions will
    use lz4 to compress, since lz4 has better performance
    that any others.

    The compressed format is size + data + size + data + ..., where size is a 64-bit
    little-endian integer.

    Parameters
    ----------
    x : str or bytes
        The line to compress.

    Returns
    -------
    bytes
        The compressed line, with a linebreak in the end.
    """
    if isinstance(x, str):
        x = x.encode()
    compress_block = lz4.frame.compress(x, compression_level=0)
    length_bytes = len(compress_block).to_bytes(64, byteorder="little")
    return length_bytes + compress_block


def decompress(x: bytes, isbytes: bool = False) -> str | bytes:
    """Decompress the line.

    Parameters
    ----------
    x : bytes
        The line to decompress.
    isbytes : bool, optional, default: False
        If the decompressed content is bytes. If not, the line will be decoded.

    Returns
    -------
    str or bytes
        The decompressed line.
    """
    x = lz4.frame.decompress(x[64:])
    if isbytes:
        return x
    return x.decode()


def listtobytes(x: Any) -> bytes:
    """Convert an object to a compressed line.

    Parameters
    ----------
    x : object
        The object to convert, such as numpy.ndarray.

    Returns
    -------
    bytes
        The compressed line.
    """
    return compress(pickle.dumps(x))


def read_compressed_block(f: BinaryIO) -> Generator[bytes, None, None]:
    """Read compressed binary file, assuming the format is size + data + size + data + ...

    Parameters
    ----------
    f : fileObject
        The file object to read.

    Yields
    ------
    data: bytes
        The compressed block.
    """
    while True:
        sizeb = f.read(64)
        if not sizeb:
            break
        size = int.from_bytes(sizeb, byteorder="little")
        yield sizeb + f.read(size)


def bytestolist(x: bytes) -> Any:
    """Convert a compressed line to an object.

    Parameters
    ----------
    x : bytes
        The compressed line.

    Returns
    -------
    object
        The decompressed object.
    """
    data = decompress(x, isbytes=True)
    assert isinstance(data, bytes)
    return pickle.loads(data)


class _DiskOrderedResultSpool:
    """Store completed out-of-order results on disk until they can be yielded."""

    def __init__(self, directory: str | None = None) -> None:
        self.directory = directory
        self._temporary_directory: tempfile.TemporaryDirectory | None = None
        self._paths: dict[int, str] = {}

    def __enter__(self) -> "_DiskOrderedResultSpool":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()

    @property
    def pending_count(self) -> int:
        """Return the number of results waiting for ordered consumption."""
        return len(self._paths)

    def has(self, index: int) -> bool:
        """Return whether a result is available for an input index."""
        return index in self._paths

    def put(self, index: int, result: tuple[bool, bytes]) -> None:
        """Write one serialized out-of-order result to the spool."""
        if self.has(index):
            raise RuntimeError("Ordered result contains a duplicate index")
        if self._temporary_directory is None:
            self._temporary_directory = tempfile.TemporaryDirectory(
                prefix="reacnetgenerator-ordered-results-",
                dir=self.directory,
            )
        path = os.path.join(self._temporary_directory.name, f"{index}.pickle")
        try:
            with open(path, "wb") as f:
                succeeded, payload = result
                # The one-byte status header is parent-owned metadata; the
                # remaining bytes are the worker's untouched ForkingPickler
                # payload and are decoded only when their index is consumed.
                f.write(b"\x01" if succeeded else b"\x00")
                f.write(payload)
        except BaseException:
            if os.path.exists(path):
                os.unlink(path)
            raise
        self._paths[index] = path

    def pop(self, index: int) -> Any:
        """Read and remove one result from the spool."""
        try:
            path = self._paths.pop(index)
        except KeyError as e:
            raise RuntimeError("Ordered result is not available") from e
        with open(path, "rb") as f:
            data = f.read()
        if not data:
            raise RuntimeError("Ordered result spool entry is empty")
        succeeded = data[0] == 1
        value = (succeeded, data[1:])
        os.unlink(path)
        return value

    def close(self) -> None:
        """Remove every temporary result and its containing directory."""
        self._paths.clear()
        if self._temporary_directory is not None:
            self._temporary_directory.cleanup()
            self._temporary_directory = None


def _record_pool_result(
    completed: Queue[tuple[int, bool, Any]],
    index: int,
    succeeded: bool,
    value: Any,
) -> None:
    """Send an asynchronous pool result to the bounded consumer."""
    completed.put((index, succeeded, value))


def _record_serialized_pool_result(
    completed: Queue[tuple[int, bool, Any]],
    index: int,
    result: tuple[bool, bytes],
) -> None:
    """Forward safe worker bytes so ordered mode can spool before decoding."""
    completed.put((index, True, result))


def _decode_pool_result(
    callback_succeeded: bool,
    value: Any,
) -> tuple[bool, Any]:
    """Decode a delivered worker result while preserving callback failures."""
    if not callback_succeeded:
        return False, value
    worker_succeeded, payload = value
    try:
        decoded = pickle.loads(payload)
    except BaseException as error:
        return False, error
    return worker_succeeded, decoded


_POOL_TASK_EVENTS: Any | None = None


def _init_bounded_pool_worker(task_events: Any) -> None:
    """Give bounded-pool workers a synchronous task-lifecycle channel."""
    global _POOL_TASK_EVENTS
    _POOL_TASK_EVENTS = task_events


def _bounded_pool_initializer(
    task_events: Any, start_gate: Any, startup_events: Any
) -> None:
    """Hold workers until the parent records their initial process objects."""
    start_gate.wait()
    _init_bounded_pool_worker(task_events)
    startup_events.put(os.getpid())


def _forking_pickle_dumps(value: Any) -> bytes:
    """Serialize with the same registered reducers used by multiprocessing."""
    return bytes(ForkingPickler.dumps(value, protocol=pickle.HIGHEST_PROTOCOL))


def _run_bounded_pool_item(index: int, call_payload: bytes) -> tuple[bool, bytes]:
    """Serialize one result before declaring its worker task complete."""
    if _POOL_TASK_EVENTS is None:
        raise RuntimeError("Bounded pool worker was not initialized")
    pid = os.getpid()
    _POOL_TASK_EVENTS.put(("started", pid, index))
    try:
        func, item = pickle.loads(call_payload)
        try:
            value = func(item)
        except Exception as error:
            succeeded = False
            # Match Pool's legacy error propagation, including the worker-side
            # traceback restored as the error cause in the parent.
            value = ExceptionWithTraceback(error, error.__traceback__)
        else:
            succeeded = True
        try:
            payload = _forking_pickle_dumps(value)
        except Exception as error:
            # Match Pool's behavior for values or exceptions that cannot be
            # serialized by reporting the serialization failure to the parent.
            succeeded = False
            payload = _forking_pickle_dumps(error)
    except BaseException:
        # This includes SystemExit raised by either the task or an object's
        # reduction method. os._exit bypasses this handler, leaving the task
        # active so the parent can identify the exited worker by PID.
        _POOL_TASK_EVENTS.put(("fatal", pid, index))
        raise

    # Pool now only needs to serialize a bool and bytes. User-controlled
    # reduction code cannot run after this lifecycle event.
    _POOL_TASK_EVENTS.put(("finished", pid, index))
    return succeeded, payload


def _drain_pool_task_events(task_events: Any, active_tasks: dict[int, int]) -> None:
    """Apply all task lifecycle events written synchronously by workers."""
    while task_events._reader.poll():
        state, pid, index = task_events.get()
        if state == "started":
            active_tasks[pid] = index
        elif state == "finished":
            if active_tasks.get(pid) == index:
                del active_tasks[pid]
        elif state == "fatal":
            raise RuntimeError(
                f"run_mp worker {pid} exited while processing item {index} "
                "without reporting a result"
            )
        else:
            raise RuntimeError(f"run_mp received an unknown worker event {state!r}")


def _check_pool_workers(
    pool,
    known_workers: dict[int, Any],
    task_events: Any,
    active_tasks: dict[int, int],
) -> None:
    """Raise if a worker exits while it still owns an unreported task."""
    _drain_pool_task_events(task_events, active_tasks)
    for worker in getattr(pool, "_pool", ()):
        known_workers[id(worker)] = worker

    workers_by_pid = {
        worker.pid: worker
        for worker in known_workers.values()
        if worker.pid is not None
    }
    for pid, index in active_tasks.items():
        worker = workers_by_pid.get(pid)
        if worker is None:
            raise RuntimeError(
                f"run_mp worker {pid} exited while processing item {index} "
                "without reporting a result"
            )
        if worker.exitcode is not None:
            raise RuntimeError(
                f"run_mp worker {pid} exited unexpectedly with code "
                f"{worker.exitcode} while processing item {index}"
            )

    for worker_id, worker in tuple(known_workers.items()):
        if worker.exitcode is not None:
            raise RuntimeError(
                f"run_mp worker {worker.pid} exited unexpectedly "
                f"with code {worker.exitcode}"
            )


def _wait_for_pool_initialization(
    known_workers: dict[int, Any], startup_events: Any
) -> None:
    """Wait until every initially captured worker has completed initialization."""
    pending = {
        worker.pid for worker in known_workers.values() if worker.pid is not None
    }
    while pending:
        while startup_events._reader.poll():
            pending.discard(startup_events.get())
        for worker in known_workers.values():
            if worker.pid in pending and worker.exitcode is not None:
                raise RuntimeError(
                    f"run_mp worker {worker.pid} exited unexpectedly "
                    f"with code {worker.exitcode} during initialization"
                )
        if pending:
            time.sleep(0.01)


def _wait_pool_result(
    pool,
    completed: Queue[tuple[int, bool, Any]],
    known_workers: dict[int, Any],
    task_events: Any,
    active_tasks: dict[int, int],
) -> tuple[int, bool, Any]:
    """Wait for a result while detecting workers that cannot report an error."""
    while True:
        _check_pool_workers(
            pool,
            known_workers,
            task_events,
            active_tasks,
        )
        try:
            return completed.get(timeout=0.05)
        except Empty:
            pass


def _submit_pool_item(
    pool,
    completed,
    pending,
    func,
    index,
    item,
) -> None:
    """Submit one item and arrange for success or failure to reach the consumer."""
    # Pool's task feeder normally serializes func/item before the worker can
    # identify which task it owns. Pre-serialize them so the worker receives
    # only safe builtins and can mark the task active before deserialization.
    call_payload = _forking_pickle_dumps((func, item))
    pending[index] = pool.apply_async(
        _run_bounded_pool_item,
        (index, call_payload),
        callback=partial(_record_serialized_pool_result, completed, index),
        error_callback=partial(_record_pool_result, completed, index, False),
    )


def _bounded_pool_results(
    pool,
    func: Callable,
    items: Iterable[Any],
    max_inflight: int,
    *,
    ordered_total: int | None = None,
    spool: _DiskOrderedResultSpool | None = None,
    task_events: Any,
    known_workers: dict[int, Any] | None = None,
) -> Generator[Any, None, None]:
    """Yield pool results with an explicit submitted-but-not-consumed bound."""
    source = iter(items)
    completed: Queue[tuple[int, bool, Any]] = Queue(maxsize=max_inflight)
    pending: dict[int, Any] = {}
    if known_workers is None:
        known_workers = {}
    active_tasks: dict[int, int] = {}
    submitted = 0
    # Capture the initial pool members before task submission can let the
    # reaper replace a worker that fails during startup or module import.
    for worker in getattr(pool, "_pool", ()):
        known_workers[id(worker)] = worker

    if spool is None:
        outstanding = 0
        source_exhausted = False
        while not source_exhausted or outstanding:
            while not source_exhausted and outstanding < max_inflight:
                try:
                    item = next(source)
                except StopIteration:
                    source_exhausted = True
                    break
                _submit_pool_item(
                    pool,
                    completed,
                    pending,
                    func,
                    submitted,
                    item,
                )
                submitted += 1
                outstanding += 1
            if not outstanding:
                break
            index, callback_succeeded, value = _wait_pool_result(
                pool,
                completed,
                known_workers,
                task_events,
                active_tasks,
            )
            if pending.pop(index, None) is None:
                raise RuntimeError("run_mp received an unknown result index")
            succeeded, value = _decode_pool_result(callback_succeeded, value)
            if not succeeded:
                if isinstance(value, BaseException):
                    raise value
                raise RuntimeError("run_mp worker failed without an exception")
            yield value
            outstanding -= 1
        return

    assert ordered_total is not None
    next_index = 0
    source_exhausted = False
    while next_index < ordered_total:
        while (
            not source_exhausted
            and submitted < ordered_total
            and submitted - next_index < max_inflight
        ):
            try:
                item = next(source)
            except StopIteration:
                source_exhausted = True
                break
            _submit_pool_item(
                pool,
                completed,
                pending,
                func,
                submitted,
                item,
            )
            submitted += 1

        if not pending:
            raise RuntimeError("Ordered result count does not match declared total")

        result_index, callback_succeeded, value = _wait_pool_result(
            pool,
            completed,
            known_workers,
            task_events,
            active_tasks,
        )
        if pending.pop(result_index, None) is None:
            raise RuntimeError("run_mp received an unknown result index")
        if not callback_succeeded:
            if isinstance(value, BaseException):
                raise value
            raise RuntimeError("run_mp worker failed without an exception")
        if result_index < next_index or spool.has(result_index):
            raise RuntimeError("Ordered result contains a duplicate index")
        if result_index == next_index:
            succeeded, value = _decode_pool_result(True, value)
            if not succeeded:
                if isinstance(value, BaseException):
                    raise value
                raise RuntimeError("run_mp worker failed without an exception")
            yield value
            next_index += 1
            while spool.has(next_index):
                succeeded, payload = spool.pop(next_index)
                succeeded, value = _decode_pool_result(True, (succeeded, payload))
                if not succeeded:
                    if isinstance(value, BaseException):
                        raise value
                    raise RuntimeError("run_mp worker failed without an exception")
                yield value
                next_index += 1
        else:
            spool.put(result_index, value)

    try:
        next(source)
    except StopIteration:
        pass
    else:
        raise RuntimeError("Ordered result count does not match declared total")
    if pending or spool.pending_count:
        raise RuntimeError("Ordered result count does not match declared total")


def _bounded_multiopen(
    pool,
    func: Callable,
    l: IO,
    *,
    max_inflight: int,
    spool: _DiskOrderedResultSpool | None,
    task_events: Any,
    known_workers: dict[int, Any] | None = None,
    nlines: int | None = None,
    unordered: bool = True,
    return_num: bool = False,
    start: int = 0,
    extra: Any | None = None,
    interval: int | None = None,
    bar: bool = True,
    desc: str | None = None,
    unit: str = "it",
    total: int | None = None,
    chunksize: int = 1,
) -> Iterable[Any]:
    """Prepare inputs for the opt-in bounded ``run_mp`` execution path."""
    if chunksize != 1:
        raise ValueError("bounded run_mp requires chunksize=1")
    obj = l
    if nlines:
        obj = itertools.zip_longest(*[obj] * nlines)
    if interval:
        obj = itertools.islice(obj, 0, None, interval)
    if return_num:
        obj = enumerate(obj, start)
    if extra is not None:
        obj = ((item, extra) for item in obj)

    ordered_total = None
    if spool is not None:
        if unordered:
            raise ValueError("disk_ordered requires unordered=False")
        if total is None:
            raise ValueError("disk_ordered requires an exact total")
        ordered_total = int(total)
        if ordered_total < 0:
            raise ValueError("disk_ordered total must be non-negative")
    elif not unordered:
        raise ValueError("bounded ordered run_mp requires disk_ordered=True")

    results = _bounded_pool_results(
        pool,
        func,
        obj,
        max_inflight,
        ordered_total=ordered_total,
        spool=spool,
        task_events=task_events,
        known_workers=known_workers,
    )
    if bar:
        results = tqdm(results, desc=desc, unit=unit, total=total, disable=None)
    return results


def listtostirng(
    l: str | list | tuple | np.ndarray, sep: list[str] | tuple[str, ...]
) -> str:
    """Convert a list to string, that is easier to store.

    Parameters
    ----------
    l : str or array-like
        The list to convert, which can contain any number of dimensions.
    sep : list of strs
        The seperators for each dimension.

    Returns
    -------
    str
        The converted string.
    """
    if isinstance(l, str):
        return l
    if isinstance(l, (list, tuple, np.ndarray)):
        return sep[0].join(listtostirng(x, sep[1:]) for x in l)
    return str(l)


def multiopen(
    pool: "multiprocessing.pool.Pool",
    func: Callable,
    l: IO,
    semaphore: Optional["multiprocessing.synchronize.Semaphore"] = None,
    nlines: int | None = None,
    unordered: bool = True,
    return_num: bool = False,
    start: int = 0,
    extra: Any | None = None,
    interval: int | None = None,
    bar: bool = True,
    desc: str | None = None,
    unit: str = "it",
    total: int | None = None,
) -> Iterable:
    """Return an interated object for process a file with multiple processors.

    Parameters
    ----------
    pool : multiprocessing.Pool
        The pool for multiprocessing.
    func : function
        The function to process lines.
    l : File object
        The file object.
    semaphore : multiprocessing.Semaphore, optional, default: None
        The semaphore to acquire. If None (default), the object will be passed
        without control.
    nlines : int, optional, default: None
        The number of lines to pass to the function each time. If None (default),
        only one line will be passed to the function.
    unordered : bool, optional, default: True
        Whether the process can be unordered.
    return_num : bool, optional, default: False
        If True, adds a counter to an iterable.
    start : int, optional, default: 0
        The start number of the counter.
    extra : object, optional, default: None
        The extra object passed to the item.
    interval : int, optional, default: None
        The interval of items that will be passed to the function. For example,
        if set to 10, a item will be passed once every 10 items and others will
        be dropped.
    bar : bool, optional, default: True
        If True, show a tqdm bar for the iteration.
    desc : str, optional, default: None
        The description of the iteration shown in the bar.
    unit : str, optional, default: it
        The unit of the iteration shown in the bar.
    total : int, optional, default: None
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
        obj = tqdm(obj, desc=desc, unit=unit, total=total, disable=None)
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
    """Share ReacNetGenerator data with a class of the submodule.

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

    def __init__(
        self,
        rng: "reacnetgenerator.ReacNetGenerator",
        usedRNGKeys: list[str],
        returnedRNGKeys: list[str],
        extraNoneKeys: list[str] | None = None,
    ) -> None:
        self.rng = rng
        self.returnedRNGKeys = returnedRNGKeys
        for key in usedRNGKeys:
            setattr(self, key, getattr(self.rng, key))
        for key in returnedRNGKeys:
            setattr(self, key, None)
        if extraNoneKeys is not None:
            for key in extraNoneKeys:
                setattr(self, key, None)

    def returnkeys(self) -> None:
        """Return back keys to ReacNetGenerator class."""
        for key in self.returnedRNGKeys:
            setattr(self.rng, key, getattr(self, key))


def checksha256(filename: str, sha256_check: str | list[str]):
    """Check sha256 of a file is correct.

    Parameters
    ----------
    filename : str
        The filename.
    sha256_check : str or list of strs
        The sha256 to be checked.

    Returns
    -------
    bool
        Indicate whether sha256 is correct.
    """
    if not os.path.isfile(filename):
        return
    h = hashlib.sha256()
    b = bytearray(128 * 1024)
    mv = memoryview(b)
    with open(filename, "rb", buffering=0) as f:
        for n in iter(lambda: f.readinto(mv), 0):
            h.update(mv[:n])
    sha256 = h.hexdigest()
    logger.info(f"SHA256 of {filename}: {sha256}")
    if sha256 in must_be_list(sha256_check):
        return True
    logger.warning("SHA256 is not correct.")
    logger.warning(open(filename).read())
    return False


async def download_file(
    urls: str | list[str], pathfilename: str, sha256: str | None
) -> str:
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
    s.mount("http://", HTTPAdapter(max_retries=3))
    s.mount("https://", HTTPAdapter(max_retries=3))
    # download if not exists
    if os.path.isfile(pathfilename) and (
        sha256 is None or checksha256(pathfilename, sha256)
    ):
        return pathfilename

    # from https://stackoverflow.com/questions/16694907
    for url in must_be_list(urls):
        logger.info(f"Try to download {pathfilename} from {url}")
        with s.get(url, stream=True) as r, open(pathfilename, "wb") as f:
            try:
                shutil.copyfileobj(r.raw, f)
                break
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request {pathfilename} Error.", exc_info=e)
    else:
        raise RuntimeError(f"Cannot download {pathfilename}.")

    return pathfilename


async def gather_download_files(urls: list[dict]) -> None:
    """Asynchronously download files from remote urls if not exists.

    See download_multifiles function for details.

    See Also
    --------
    download_multifiles
    """
    await asyncio.gather(
        *[
            download_file(jdata["url"], jdata["fn"], jdata.get("sha256", None))
            for jdata in urls
        ]
    )


def download_multifiles(urls: list[dict]) -> None:
    """Download multiple files from dicts.

    Parameters
    ----------
    urls : list of dicts
        The information of download files. Each dict should contain the following key:
            - url: str or list of strs
                The url(s) that is available to download.
            - pathfilename: str
                The downloading path of the file.
            - sha256: str, optional, default: None
                Sha256 of the file. If not None and match the file, the download will be skiped.
    """
    asyncio.run(gather_download_files(urls))


def run_mp(
    nproc: int,
    *,
    max_inflight: int | None = None,
    disk_ordered: bool = False,
    ordered_spool_dir: str | None = None,
    **kwargs: Any,
) -> Iterable[Any]:
    """Process a file with multiple processors.

    Parameters
    ----------
    nproc : int
        The number of processors to be used.
    max_inflight : int, optional
        Maximum number of submitted inputs not yet delivered to the consumer.
        If omitted, the existing ``nproc * 150`` semaphore path is unchanged.
    disk_ordered : bool, optional, default: False
        Run workers to completion out of order, spool early results to disk, and
        yield them in input order. This requires ``unordered=False`` and an exact
        ``total``.
    ordered_spool_dir : str, optional
        Parent directory for temporary ordered-result files.
    **kwargs : dict, optional
        Other parameters can be found in the `multiopen` method.

    Yields
    ------
    object
        The yielded object from the `multiopen` method.

    See Also
    --------
    multiopen
    """
    if max_inflight is not None:
        max_inflight = int(max_inflight)
        if max_inflight <= 0:
            raise ValueError("max_inflight must be a positive integer")
    elif disk_ordered:
        max_inflight = nproc * 150

    task_events = SimpleQueue() if max_inflight is not None else None
    start_gate = Event() if task_events is not None else None
    startup_events = SimpleQueue() if task_events is not None else None
    if task_events is None:
        pool = Pool(nproc, maxtasksperchild=1000)
    else:
        # ForkingPickler reducers such as socket/DupFd can depend on a worker's
        # resource_sharer until the parent reconstructs the delivered result.
        # Keep bounded workers alive through consumption instead of recycling
        # them immediately after the 1000th task.
        pool = Pool(
            nproc,
            initializer=_bounded_pool_initializer,
            initargs=(task_events, start_gate, startup_events),
        )
    known_workers = (
        {id(worker): worker for worker in getattr(pool, "_pool", ())}
        if task_events is not None
        else None
    )
    if start_gate is not None:
        # Workers wait in the initializer so no startup failure can be reaped
        # before this snapshot preserves its Process object for monitoring.
        start_gate.set()
    semaphore = Semaphore(nproc * 150) if max_inflight is None else None
    try:
        if task_events is not None:
            assert startup_events is not None
            assert known_workers is not None
            _wait_for_pool_initialization(known_workers, startup_events)
        if max_inflight is None:
            results = multiopen(pool=pool, semaphore=semaphore, **kwargs)
        elif disk_ordered:
            assert task_events is not None
            with _DiskOrderedResultSpool(ordered_spool_dir) as spool:
                results = _bounded_multiopen(
                    pool=pool,
                    max_inflight=max_inflight,
                    spool=spool,
                    task_events=task_events,
                    known_workers=known_workers,
                    **kwargs,
                )
                for item in results:
                    yield item
        else:
            assert task_events is not None
            results = _bounded_multiopen(
                pool=pool,
                max_inflight=max_inflight,
                spool=None,
                task_events=task_events,
                known_workers=known_workers,
                **kwargs,
            )
            for item in results:
                yield item
        if max_inflight is None:
            for item in results:
                yield item
                if semaphore is not None:
                    semaphore.release()
    except:
        logger.exception("run_mp failed")
        pool.terminate()
        raise
    else:
        pool.close()
    finally:
        pool.join()
        if task_events is not None:
            task_events.close()
        if startup_events is not None:
            startup_events.close()


def must_be_list(obj: Any | list[Any]) -> list[Any]:
    """Convert a object to a list if the object is not a list.

    Parameters
    ----------
    obj : Object
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


def get_timestep_value(timestep: Any) -> Any:
    """Normalize stored timestep metadata to the timestep value."""
    if isinstance(timestep, tuple):
        timestep = timestep[-1]
    if isinstance(timestep, np.generic):
        return timestep.item()
    return timestep


def check_zero_signal(signal: np.ndarray) -> bool:
    """Check if the given signal contains only zeros.

    Parameters
    ----------
    signal : 1D array of bool
        The signal to check. The dtype should be bool.

    Returns
    -------
    bool
        False if the signal contains only zeros, True otherwise.
    """
    # Benchmark
    # one_million_ones = np.ones(10**6, dtype=bool)
    # %timeit reacnetgenerator.utils_np.check_zero_signal(one_million_ones)
    # 808 ns ± 197 ns per loop (mean ± std. dev. of 7 runs, 1,000,000 loops each)
    # %timeit one_million_ones.any()
    # 17.9 µs ± 1.01 µs per loop (mean ± std. dev. of 7 runs, 100,000 loops each)
    # %timeit one_million_ones[one_million_ones.argmax()]
    # 561 ns ± 135 ns per loop (mean ± std. dev. of 7 runs, 1,000,000 loops each)
    #
    # any() doesn't have short-circuits, but argmax() does for bool.
    # See https://stackoverflow.com/a/45774536/9567349
    # ty 0.0.15 cannot resolve NumPy's no-argument argmax overload.
    return signal[signal.argmax()].item()  # ty: ignore[no-matching-overload]


def idx_to_signal(idx: np.ndarray, step: int):
    """Convert an index array to a signal array.

    Parameters
    ----------
    idx : array_like
        Index array.
    step : int
        Step size.

    Returns
    -------
    signal : ndarray
        Signal array in int8.
    """
    # Cython implementation is 1-2x faster than Python implementation.
    # However, with it, it's hard to use the limited Python API:
    # (1) https://github.com/cython/cython/issues/5697
    # (2) memory view limited API only available since python 3.11
    # when 1 is resolved, we may add back the Cython implementation
    # for Python 3.11+ only (build two wheels).

    signal = np.zeros(step, dtype=np.int8)
    signal[idx] = 1
    return signal.reshape((step, 1))
