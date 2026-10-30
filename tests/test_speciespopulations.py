# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare every population layout with the original per-frame Counters."""

from collections import Counter

import numpy as np
import pytest

from reacnetgenerator import _speciespopulations as populations
from reacnetgenerator._matrix import _GenerateMatrix
from reacnetgenerator._step3state import _MoleculeNameTable
from reacnetgenerator.utils import listtobytes


def _source(path, frames):
    """Encode four-field molecule records with the requested observed frame indices."""
    with path.open("wb") as handle:
        for values in frames:
            for value in ([0], [], [], np.asarray(values, dtype=np.int64)):
                handle.write(listtobytes(value))


def _legacy(names, records, timesteps):
    """Build an independent byte-for-byte oracle from the original ordered Counters."""
    rows = [Counter() for _ in timesteps]
    for name, frames in zip(names, records):
        for frame in frames:
            rows[frame][name] += 1
    return "".join(
        f"Timestep {timestep}:"
        + "".join(f" {name} {count}" for name, count in row.items())
        + "\n"
        for timestep, row in zip(timesteps, rows, strict=True)
    ).encode()


def _run(tmp_path, names, records, timesteps):
    """Check matrix-entry output against the oracle and verify scratch-file cleanup."""
    source = tmp_path / "records.bin"
    output = tmp_path / "result.species"
    _source(source, records)
    matrix = object.__new__(_GenerateMatrix)
    matrix.moleculetemp2filename = str(source)
    matrix.speciesfilename = str(output)
    matrix.mname, matrix.timestep = names, timesteps
    stats = matrix._printspecies()
    assert output.read_bytes() == _legacy(names, records, timesteps)
    assert source.exists()
    assert not list(tmp_path.glob("reacnetgenerator-species-*"))
    return stats


@pytest.mark.parametrize("compact_names", [False, True])
@pytest.mark.parametrize(
    "layout",
    [
        "small",
        "dense",
        "sparse",
        "spool-cap",
        "spool-failure",
        "spool-write-failure",
        "rescan",
    ],
)
def test_layouts_preserve_exact_counts_and_first_appearance(
    tmp_path, monkeypatch, compact_names, layout
):
    """Shared names, negative indices, duplicate visits and local order survive."""
    monkeypatch.setattr(populations.tempfile, "tempdir", str(tmp_path))
    names = np.array(["氢", "B", "氢", "C"], dtype=object)
    records = [[0, 1, 3, 3, 399], [1, 5, 6, -1], [5, 4, 0, 1], [2, 5]]
    if layout == "dense":
        records = [record + list(range(400)) * 3 for record in records]
    if compact_names:
        names = _MoleculeNameTable.from_names(names)
    if layout != "small":
        monkeypatch.setattr(populations, "_SMALL_OBSERVATION_ROWS", 0)
    if layout not in ("small", "dense"):
        monkeypatch.setattr(populations, "_SPECIES_COUNT_BLOCK_BYTES", 4096)
    if layout == "spool-cap":
        monkeypatch.setattr(populations, "_SPECIES_OBSERVATION_SPOOL_BYTES", 1)
    original = populations.tempfile.mkstemp

    def create(**kwargs):
        """Inject selected scratch-creation failures without replacing unrelated files."""
        if layout == "rescan" or (
            layout == "spool-failure" and "observed" in kwargs["prefix"]
        ):
            raise OSError("scratch storage unavailable")
        return original(**kwargs)

    monkeypatch.setattr(populations.tempfile, "mkstemp", create)
    if layout == "spool-write-failure":

        class FailedWriter:
            def __init__(self, wrapped):
                """Retain the real file so only observation writes are replaced."""
                self.wrapped = wrapped

            def write(self, data):
                """Simulate a storage error while buffered observations are being flushed."""
                raise OSError("observation scratch write failed")

            def __getattr__(self, name):
                """Delegate closing and other file operations to the owned real handle."""
                return getattr(self.wrapped, name)

        fdopen = populations.os.fdopen
        monkeypatch.setattr(
            populations.os, "fdopen", lambda *args: FailedWriter(fdopen(*args))
        )
        # Trigger an actual buffered write during observation; the incomplete
        # spool must be discarded and observations reread from the source.
        monkeypatch.setattr(populations, "_SPECIES_EVENT_CHUNK_ROWS", 3)
    stats = _run(tmp_path, names, records, np.arange(400) * 10 + 100)
    expected = "sparse" if layout.startswith("spool-") else layout
    assert stats["layout"] == expected
    if layout.startswith("spool-"):
        assert not stats["observation_reused"]


@pytest.mark.parametrize("frame_count,records", [(0, []), (5, []), (5, [[], []])])
def test_empty_timelines_need_no_scratch(tmp_path, monkeypatch, frame_count, records):
    """Empty output and empty frames must not allocate temporary layouts."""

    def fail(**kwargs):
        """Reject any scratch allocation on the small or empty input path."""
        raise AssertionError("Small/empty input created scratch storage")

    monkeypatch.setattr(populations.tempfile, "mkstemp", fail)
    names = ["A"] * len(records)
    assert _run(tmp_path, names, records, list(range(frame_count)))["layout"] == "small"


def test_many_species_in_one_frame_exceed_small_memory_budget(tmp_path, monkeypatch):
    """A short timeline alone must not qualify for the dictionary fast path."""
    monkeypatch.setattr(populations, "_SMALL_COUNT_BYTES", 2048)
    names = [f"S{i}" for i in range(100, 0, -1)]
    assert _run(tmp_path, names, [[0]] * 100, [77])["layout"] == "dense"


@pytest.mark.parametrize("spool_limit", [0, 2**20])
def test_sparse_refinement_and_single_frame_overflow(
    tmp_path, monkeypatch, spool_limit
):
    """Duplicate-heavy frames use counts, not an arbitrarily large event sort."""
    monkeypatch.setattr(populations, "_SMALL_OBSERVATION_ROWS", 0)
    monkeypatch.setattr(populations, "_SPECIES_COUNT_BLOCK_BYTES", 1024)
    monkeypatch.setattr(populations, "_SPECIES_OBSERVATION_SPOOL_BYTES", spool_limit)
    monkeypatch.setattr(populations.tempfile, "tempdir", str(tmp_path))
    stats = _run(tmp_path, ["A", "B"], [[0] * 600, [0] * 300], range(10000))
    assert stats["layout"] == "dense"
    assert stats["working_bytes"] <= 1024


@pytest.mark.parametrize("layout", ["small", "dense", "sparse"])
def test_random_populations_match_counter_order(tmp_path, monkeypatch, layout):
    """Byte equality checks per-frame order as well as duplicate multiplicity."""
    if layout != "small":
        monkeypatch.setattr(populations, "_SMALL_OBSERVATION_ROWS", 0)
    if layout == "sparse":
        monkeypatch.setattr(populations, "_SPECIES_COUNT_BLOCK_BYTES", 512)
    random = np.random.default_rng(892)
    for _ in range(15):
        names = random.choice(["A", "B", "C", "D"], 12)
        records = [
            random.integers(-400, 400, random.integers(0, 10)).tolist() for _ in names
        ]
        if layout == "dense":
            records = [record + list(range(400)) for record in records]
        stats = _run(tmp_path, names, records, range(400))
        assert stats["layout"] == layout


def test_scratch_is_removed_when_output_fails(tmp_path, monkeypatch):
    """An output error must propagate after closing both temporary files."""
    monkeypatch.setattr(populations, "_SMALL_OBSERVATION_ROWS", 0)
    monkeypatch.setattr(populations, "_SPECIES_COUNT_BLOCK_BYTES", 4096)
    monkeypatch.setattr(populations.tempfile, "tempdir", str(tmp_path))
    source = tmp_path / "records.bin"
    _source(source, [[0, 1], [2]])
    with pytest.raises(IsADirectoryError):
        populations._write_species_populations(source, ["A", "B"], range(400), tmp_path)
    assert not list(tmp_path.glob("reacnetgenerator-species-*"))


def test_truncated_observation_record_is_rejected(tmp_path):
    """Skipping structure payloads must still detect a truncated last record."""
    source = tmp_path / "records.bin"
    _source(source, [[0]])
    source.write_bytes(source.read_bytes()[:-1])
    with pytest.raises(ValueError, match="Incomplete molecule record"):
        populations._write_species_populations(
            source, ["A"], range(4), tmp_path / "out"
        )


@pytest.mark.parametrize("spool_limit", [0, 2**20])
@pytest.mark.parametrize("overflow", [False, True])
def test_sparse_refinement_uses_logarithmic_passes(
    tmp_path, monkeypatch, spool_limit, overflow
):
    """A partition just over budget must not repeatedly reread almost all input."""
    monkeypatch.setattr(populations, "_SMALL_OBSERVATION_ROWS", 0)
    monkeypatch.setattr(populations, "_SPECIES_COUNT_BLOCK_BYTES", 1024)
    monkeypatch.setattr(populations, "_SPECIES_OBSERVATION_SPOOL_BYTES", spool_limit)
    monkeypatch.setattr(populations.tempfile, "tempdir", str(tmp_path))
    frame_count = 4096
    with populations._SpeciesTimelineStore(frame_count, 2, 3) as store:
        row_limit = (
            populations._SPECIES_COUNT_BLOCK_BYTES
            // store._sparse_working_bytes_per_row
        )
        refinement_limit = (store.block_frame_count - 1).bit_length()

    reads = {"source": 0, "spool": 0}
    source_records = populations._frame_records
    spool_records = populations._SpeciesTimelineStore._iter_observations

    def read_source(*args, **kwargs):
        """Count full input passes, including small-input probing and population."""
        reads["source"] += 1
        assert reads["source"] <= refinement_limit + 3
        yield from source_records(*args, **kwargs)

    def read_spool(self):
        """Count spool passes, allowing one final pass to populate the layout."""
        reads["spool"] += 1
        assert reads["spool"] <= refinement_limit + 1
        yield from spool_records(self)

    monkeypatch.setattr(populations, "_frame_records", read_source)
    monkeypatch.setattr(
        populations._SpeciesTimelineStore, "_iter_observations", read_spool
    )
    crowded_frame = frame_count // 2
    # Either one frame is just over budget or adjacent frames can be separated
    # into valid sparse partitions. Both preserve the same near-limit workload.
    records = [
        [crowded_frame] * row_limit,
        [crowded_frame if overflow else crowded_frame + 1],
        [0, frame_count - 1],
    ]
    stats = _run(tmp_path, ["B", "A", "B"], records, range(frame_count))
    assert stats["layout"] == ("dense" if overflow else "sparse")
    assert stats["observation_reused"] == bool(spool_limit)
    assert reads["spool"] > 1 if spool_limit else reads["source"] > 3
