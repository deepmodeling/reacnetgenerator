# SPDX-License-Identifier: LGPL-3.0-or-later
"""Per-record SMILES state, real process contexts and exceptional cleanup."""

import multiprocessing
from pathlib import Path
from types import MethodType

import pytest

from reacnetgenerator import _path, utils
from reacnetgenerator import _smilesworker as workers

from .test_smiles_batching import (
    _assert_closed,
    _collector,
    _exited_task,
    _failed_initializer,
    _failed_task,
    _fallback_initializer,
    _records,
    _selected_failure,
    _serialization_initializer,
)


@pytest.fixture(params=multiprocessing.get_all_start_methods())
def process_context(request, monkeypatch):
    """Use actual context-specific pools and synchronization primitives."""
    context = multiprocessing.get_context(request.param)
    for name in ("Pool", "SimpleQueue", "Event", "Semaphore"):
        monkeypatch.setattr(utils, name, getattr(context, name))
    monkeypatch.setattr(workers, "get_start_method", lambda: request.param)
    return request.param


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_record_context_matches_local_output(tmp_path, process_context, miso):
    """Minimal workers preserve names, representatives, frames and timeline order."""
    local = _collector(tmp_path / "local", _records(20), miso, 1)
    parallel = _collector(tmp_path / "parallel", _records(20), miso)
    # A real collector owns an RNG. Make unrelated state impossible to pickle.
    parallel.rng = {"unrelated": lambda: None}
    children = {p.pid for p in multiprocessing.active_children()}
    for collector in (local, parallel):
        collector._printmoleculename()
    assert list(parallel.mname) == list(local.mname)
    for attribute in ("moleculefilename", "moleculetimelinefilename"):
        assert (
            Path(getattr(parallel, attribute)).read_bytes()
            == Path(getattr(local, attribute)).read_bytes()
        )
    _assert_closed(tmp_path, children)


@pytest.mark.parametrize("filter_kind", ["frames", "timesteps"])
def test_record_filters_stay_in_parent(tmp_path, process_context, filter_kind):
    """Workers return frames when a parent-only filter requires a timeline."""
    local = _collector(tmp_path / "local", _records(8), nproc=1)
    parallel = _collector(tmp_path / "parallel", _records(8))
    for collector in (local, parallel):
        collector.printmoleculetime = False
        if filter_kind == "frames":
            collector._moleculeframefilter = {0, 2, 6}
        else:
            collector._moleculetimestepfilter = {0, 20, 60}
        collector._printmoleculename()
    assert (
        Path(parallel.moleculetimelinefilename).read_bytes()
        == Path(local.moleculetimelinefilename).read_bytes()
    )


@pytest.mark.parametrize("miso", [0, 1, 2])
def test_record_fallback_preserves_frequency(
    tmp_path, monkeypatch, process_context, miso
):
    """ValueError candidates still use parent VF2 and global occurrence counts."""
    records = [
        ([0, 1], [(0, 1, order)], list(range(count)))
        for order, count in ((1, 4), (2, 4), (3, 7))
    ]
    local = _collector(tmp_path / "local", records, miso, 1)
    local.convertSMILES = MethodType(_selected_failure, local)
    parallel = _collector(tmp_path / "parallel", records, miso)

    def inject(nproc, **kwargs):
        kwargs["initializer"] = _fallback_initializer
        return utils.run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", inject)
    for collector in (local, parallel):
        collector._printmoleculename()
    assert list(parallel.mname) == list(local.mname)
    assert (
        Path(parallel.moleculefilename).read_bytes()
        == Path(local.moleculefilename).read_bytes()
    )


class _CustomCollector(_path._CollectSMILESPaths):
    def convertSMILES(self, atoms, bonds):
        return "custom-name"


def _custom_name(atoms, bonds):
    return "custom-name"


@pytest.mark.parametrize("custom", ["subclass", "instance"])
def test_custom_record_converter_keeps_bound_state(
    tmp_path, monkeypatch, process_context, custom
):
    """Existing extension points retain the original bound-method runner."""
    collector = _collector(tmp_path, _records(4))
    if custom == "subclass":
        collector.__class__ = _CustomCollector
    else:
        collector.convertSMILES = _custom_name
    calls = []

    def checked(nproc, **kwargs):
        assert kwargs["func"].__self__ is collector
        assert "initializer" not in kwargs
        calls.append(nproc)
        return utils.run_mp(nproc, **kwargs)

    monkeypatch.setattr(_path, "run_mp", checked)
    collector._printmoleculename()
    assert calls == [2]
    assert list(collector.mname) == ["custom-name"] * 4


@pytest.mark.parametrize(
    "failure", ["initializer", "task", "exit", "serialize", "consumer"]
)
def test_record_failures_close_workers_with_live_traceback(
    tmp_path, monkeypatch, process_context, failure
):
    """Close the iterator even when a retained exception keeps its frames alive."""
    collector = _collector(tmp_path, _records(8))

    def inject(nproc, **kwargs):
        if failure == "initializer":
            kwargs["initializer"] = _failed_initializer
        elif failure == "serialize":
            kwargs["initializer"] = _serialization_initializer
        elif failure in ("task", "exit"):
            kwargs["func"] = _failed_task if failure == "task" else _exited_task
        return utils.run_mp(nproc, **kwargs)

    monkeypatch.setattr(workers, "run_mp", inject)
    if failure == "consumer":

        def fail_consumer(*args, **kwargs):
            raise RuntimeError("intentional consumer failure")

        monkeypatch.setattr(
            _path._CollectSMILESPaths, "_formatmoleculename", fail_consumer
        )
    children = {p.pid for p in multiprocessing.active_children()}
    with pytest.raises(
        (RuntimeError, AttributeError), match=r"intentional|exited|pickl|local object"
    ) as caught:
        collector._printmoleculename()
    assert caught.value.__traceback__ is not None
    assert Path(collector.moleculetemp2filename).exists()
    _assert_closed(tmp_path, children)
