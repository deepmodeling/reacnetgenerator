# SPDX-License-Identifier: LGPL-3.0-or-later
"""Timeline contracts, end-to-end semantics, and bounded resource regressions."""

import csv
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator._detect import _Detect
from reacnetgenerator._hmmfilter import _HMMFilter
from reacnetgenerator._path import _CollectPaths
from reacnetgenerator._reaction import ReactionsFinder
from reacnetgenerator._timedoutput import _Column, _TimedOutputWriter
from reacnetgenerator.commandline import main_parser, parm2cmd
from reacnetgenerator.timedoutput import (
    iter_frames,
    iter_molecule_ranges,
    iter_molecules,
    iter_reaction_events,
    iter_reaction_types,
    iter_species,
    read_metadata,
)
from reacnetgenerator.utils import bytestolist, read_compressed_block

INPUT = Path(__file__).parent / "inputs" / "reaction.bond"


def generator(tmp_path, **kwargs):
    """Construct a small real bond-input analysis with isolated outputs."""
    return ReacNetGenerator(
        inputfilename=kwargs.pop("inputfilename", str(INPUT)),
        inputfiletype="bond",
        atomname=["H", "He"],
        nproc=kwargs.pop("nproc", 1),
        runHMM=kwargs.pop("runHMM", False),
        output_dir=tmp_path,
        **kwargs,
    )


def prepare(rng):
    """Run native detection and filtering without removing the signal oracle."""
    _Detect.gettype(rng).detect()
    _HMMFilter(rng).filter()


@pytest.mark.parametrize("run_hmm", [False, True])
@pytest.mark.parametrize("nproc", [1, 2])
def test_pipeline(tmp_path, run_hmm, nproc):
    """Ranges match actual matrix signals; aggregate events match legacy CSV."""
    path = tmp_path / "timeline.h5"
    rng = generator(
        tmp_path,
        timed_output=path,
        runHMM=run_hmm,
        nproc=nproc,
        printreactionevent=True,
    )
    prepare(rng)
    expected = []
    with open(rng.hmmfilename if run_hmm else rng.originfilename, "rb") as file:
        for molecule_id, blob in enumerate(read_compressed_block(file), 1):
            expected.extend(
                (molecule_id, frame) for frame in np.flatnonzero(bytestolist(blob))
            )
    _CollectPaths.getstype(rng).collect()
    actual = [
        (r.molecule_id, frame)
        for r in iter_molecule_ranges(path, block_rows=2)
        for frame in range(r.start_frame, r.end_frame + 1)
    ]
    assert actual == expected
    assert len(list(iter_molecules(path, block_rows=2))) == rng.hmmit
    with open(rng.reactioneventfilename, newline="") as file:
        csv_events = Counter(
            (int(row["Timestep_Index"]), row["Reactant"], row["Product"])
            for row in csv.DictReader(file)
        )
    types = {i: (left, right) for i, left, right, _ in iter_reaction_types(path)}
    h5_events = Counter(
        {
            (e.transition, *types[e.reaction_type_id]): e.count
            for e in iter_reaction_events(path, block_rows=1)
        }
    )
    assert h5_events == csv_events
    totals = Counter()
    for (_, left, right), count in h5_events.items():
        totals[left, right] += count
    assert {
        (left, right): count for _, left, right, count in iter_reaction_types(path)
    } == totals
    metadata = read_metadata(path)
    assert metadata["schema_version"] == "1.0"
    assert metadata["configuration"]["parameters"]["runHMM"] is run_hmm
    assert not list(tmp_path.glob("*.incomplete"))


def test_optional_and_existing_outputs(tmp_path):
    """Opt-in preserves text results and does not implicitly create timed CSV."""
    outputs = []
    for enabled in [False, True]:
        directory = tmp_path / str(enabled)
        path = directory / "timeline.h5"
        rng = generator(directory, timed_output=path if enabled else None)
        prepare(rng)
        _CollectPaths.getstype(rng).collect()
        outputs.append(
            tuple(
                Path(p).read_bytes()
                for p in [
                    rng.moleculefilename,
                    rng.atomroutefilename,
                    rng.reactionabcdfilename,
                ]
            )
        )
        assert path.exists() is enabled
        assert not Path(rng.reactioneventfilename).exists()
        assert not Path(rng.moleculetimelinefilename).exists()
    assert outputs[0] == outputs[1]


def test_source_occurrences_and_stride(tmp_path):
    """Map repeated input occurrences and global stride to source-local frames."""
    path = tmp_path / "timeline.h5"
    rng = generator(
        tmp_path,
        inputfilename=[str(INPUT), str(INPUT)],
        timed_output=path,
        stepinterval=3,
    )
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    assert rng.source_frame_counts == [8, 8]
    assert [(x.source_id, x.source_frame) for x in iter_frames(path, block_rows=2)] == [
        (0, 0),
        (0, 3),
        (0, 6),
        (1, 1),
        (1, 4),
        (1, 7),
    ]
    assert [x.timestep for x in iter_frames(path)] == [0, 3, 6, 1, 4, 7]


def test_empty_source_occurrences(tmp_path):
    """Preserve empty occurrences before, between, and after repeated inputs."""
    path = tmp_path / "timeline.h5"
    empty = tmp_path / "empty.bond"
    empty.write_text("")
    rng = generator(
        tmp_path,
        inputfilename=[str(empty), str(INPUT), str(empty), str(INPUT), str(empty)],
        timed_output=path,
        stepinterval=3,
    )
    prepare(rng)
    _CollectPaths.getstype(rng).collect()

    assert rng.source_frame_counts == [0, 8, 0, 8, 0]
    assert [(x.source_id, x.source_frame) for x in iter_frames(path)] == [
        (1, 0),
        (1, 3),
        (1, 6),
        (3, 1),
        (3, 4),
        (3, 7),
    ]


def test_atomic_failure(tmp_path, monkeypatch):
    """A PATH failure never replaces the previous artifact and closes handles."""
    path = tmp_path / "timeline.h5"
    path.write_bytes(b"previous result")
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)

    def fail(*args, **kwargs):
        raise RuntimeError("reaction failed")

    monkeypatch.setattr(ReactionsFinder, "findreactions", fail)
    with pytest.raises(RuntimeError, match="reaction failed"):
        _CollectPaths.getstype(rng).collect()
    assert path.read_bytes() == b"previous result"
    (partial,) = tmp_path.glob("*.incomplete")
    with h5py.File(partial, "r+") as file:
        assert file.attrs["status"] == "incomplete"
    with pytest.raises(ValueError, match="incomplete"):
        read_metadata(partial)


def test_missing_frame_mapping(tmp_path):
    """Missing meaning-changing metadata fails instead of fabricating a frame."""
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    rng.source_frame_counts = None
    with pytest.raises(ValueError, match="source frame mapping"):
        with _TimedOutputWriter(rng):
            pass
    assert not (tmp_path / "timeline.h5").exists()


def test_compact_counts_and_blocks(tmp_path, monkeypatch):
    """Read the first block without loading the full table or expanding counts."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        writer.write_events(
            [{"Timestep_Index": 0, "Reactant": "A", "Product": "B"}] * 5
        )
        for i in range(50):
            writer._append("molecule_ranges", molecule_id=1, start_frame=i, end_frame=i)
    calls = []
    getitem = h5py.Dataset.__getitem__

    def bounded(dataset, key, *args, **kwargs):
        if dataset.name.startswith("/molecule_ranges/"):
            assert isinstance(key, slice) and key.stop - key.start <= 3
            calls.append(key)
        return getitem(dataset, key, *args, **kwargs)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", bounded)
    ranges = iter_molecule_ranges(path, block_rows=3)
    assert next(ranges).start_frame == 0
    assert len(calls) == 3
    ranges.close()
    assert len(list(iter_reaction_events(path))) == 1
    assert next(iter_reaction_events(path)).count == 5
    with h5py.File(path, "r+"):
        pass


@pytest.mark.parametrize("block_rows", [0, -1, 1.5])
def test_invalid_block_size(tmp_path, block_rows):
    """Reject invalid block sizes before reading result arrays."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    with _TimedOutputWriter(rng):
        pass
    with pytest.raises((ValueError, TypeError)):
        next(iter_frames(path, block_rows=block_rows))


def test_column_byte_bound(tmp_path, monkeypatch):
    """Flush before the byte budget even when the row limit has not been met."""
    monkeypatch.setattr("reacnetgenerator._timedoutput._BLOCK_BYTES", 12)
    with h5py.File(tmp_path / "batch.h5", "w") as file:
        column = _Column(file, "name", h5py.string_dtype("utf-8"))
        column.append("abcdefgh")
        column.append("ijklmnop")
        assert len(column.dataset) == 1
        assert column.bytes == 8
        column.flush()
        assert list(column.dataset.asstr()[:]) == ["abcdefgh", "ijklmnop"]


def test_cli_and_aliases(tmp_path):
    """Round-trip the opt-in switch and reject paths that could destroy inputs."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    args = main_parser().parse_args(parm2cmd(rng.parameters)[1:])
    assert args.timed_output == str(path)
    with pytest.raises(ValueError, match="alias"):
        generator(tmp_path, timed_output=INPUT)
    alias = tmp_path / "alias"
    alias.symlink_to(INPUT)
    with pytest.raises(ValueError, match="alias"):
        generator(tmp_path, timed_output=alias)
    with pytest.raises(ValueError, match="executes PATH"):
        rng.run_items(["species"])


def test_empty_tables(tmp_path):
    """Complete empty tables remain readable; dictionaries need no sentinel row."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    with _TimedOutputWriter(rng):
        pass
    assert list(iter_molecule_ranges(path)) == []
    assert list(iter_reaction_events(path)) == []
    assert list(iter_species(path)) == []
    assert list(iter_reaction_types(path)) == []
    with h5py.File(path) as file:
        assert json.loads(file.attrs["capabilities"]) == [
            "molecule_ranges",
            "reaction_events",
        ]


@pytest.mark.parametrize("filetype", ["dump", "xyz", "extxyz"])
def test_coordinate_sources(tmp_path, filetype):
    """Coordinate parsers map their real input frames without changing detection."""
    path = tmp_path / "timeline.h5"
    rng = ReacNetGenerator(
        inputfilename=str(INPUT.with_name(f"water.{filetype}")),
        inputfiletype=filetype,
        atomname=["H", "O"],
        nproc=1,
        runHMM=False,
        pbc=False,
        output_dir=tmp_path,
        timed_output=path,
    )
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    assert rng.source_frame_counts == [1]
    assert [(f.source_id, f.source_frame) for f in iter_frames(path)] == [(0, 0)]
    assert list(iter_molecules(path))


def test_source_boundary_errors(tmp_path):
    """Do not silently assign a frame assembled from two partial files."""
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.write_text("a\nb\nc\n")
    second.write_text("d\n")
    rng = generator(tmp_path, inputfilename=[str(first), str(second)])
    detector = _Detect.gettype(rng)
    with pytest.raises(ValueError, match="boundary"):
        list(detector._source_lines(2))
    first.write_text("a\n")
    rng.inputfilename = [str(first)]
    detector = _Detect.gettype(rng)
    with pytest.raises(ValueError, match="final"):
        list(detector._source_lines(2))


def test_reader_headers_alignment_and_close(tmp_path):
    """Fail on unknown schemas and mismatched columns; close on early exit."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    with _TimedOutputWriter(rng):
        pass
    before = h5py.h5f.get_obj_count(types=h5py.h5f.OBJ_FILE)
    reader = iter_frames(path, block_rows=1)
    next(reader)
    assert h5py.h5f.get_obj_count(types=h5py.h5f.OBJ_FILE) == before + 1
    reader.close()
    assert h5py.h5f.get_obj_count(types=h5py.h5f.OBJ_FILE) == before
    with h5py.File(path, "r+") as file:
        file.attrs["schema_version"] = "2.0"
    with pytest.raises(ValueError, match="Unsupported"):
        read_metadata(path)
    with h5py.File(path, "r+") as file:
        file.attrs["schema_version"] = "1.0"
        file["frames/source_id"].resize((1,))
    with pytest.raises(ValueError, match="Misaligned"):
        next(iter_frames(path))


@pytest.mark.parametrize(
    "dataset",
    [
        "molecules/atom_offsets",
        "molecules/bond_offsets",
        "molecules/atom_index",
        "molecules/bond_atom_index_1",
        "molecules/bond_atom_index_2",
        "molecules/bond_order",
    ],
)
def test_molecule_reader_rejects_non_integer_columns(tmp_path, dataset):
    """Reject malformed direct numeric columns before converting their values."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    with h5py.File(path, "r+") as file:
        values = file[dataset][:].astype(float)
        del file[dataset]
        file.create_dataset(dataset, data=values)

    with pytest.raises(ValueError, match="Invalid numeric column"):
        next(iter_molecules(path))


@pytest.mark.parametrize(
    "dataset", ["molecules/atom_offsets", "molecules/bond_offsets"]
)
def test_molecule_reader_rejects_offset_length(tmp_path, dataset):
    """Require one more atom and bond offset than molecule definitions."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    with h5py.File(path, "r+") as file:
        values = file[dataset][:-1]
        del file[dataset]
        file.create_dataset(dataset, data=values)

    with pytest.raises(ValueError, match="offset length"):
        next(iter_molecules(path))


def test_molecule_reader_rejects_global_bond_misalignment(tmp_path):
    """Reject an unused tail that per-molecule bond slices would not observe."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    with h5py.File(path, "r+") as file:
        dataset = file["molecules/bond_order"]
        dataset.resize((len(dataset) + 1,))
        dataset[-1] = 1

    with pytest.raises(ValueError, match="bond columns"):
        next(iter_molecules(path))


def test_reaction_type_reader_rejects_non_integer_total(tmp_path):
    """Do not truncate malformed floating-point reaction totals."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    with h5py.File(path, "r+") as file:
        dataset = "reaction_types/total_count"
        values = file[dataset][:].astype(float) + 0.5
        del file[dataset]
        file.create_dataset(dataset, data=values)

    with pytest.raises(ValueError, match="Invalid numeric column"):
        next(iter_reaction_types(path))


def test_replace_failure(tmp_path, monkeypatch):
    """A publish failure keeps the previous destination and closes the sibling."""
    path = tmp_path / "timeline.h5"
    path.write_bytes(b"old")
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)

    def fail(*args):
        raise OSError("publish failed")

    monkeypatch.setattr("reacnetgenerator._timedoutput.os.replace", fail)
    with pytest.raises(OSError, match="publish failed"):
        with _TimedOutputWriter(rng):
            pass
    assert path.read_bytes() == b"old"
    (partial,) = tmp_path.glob("*.incomplete")
    with h5py.File(partial, "r+") as file:
        assert file.attrs["status"] == "complete"


@pytest.mark.parametrize(
    "signal,expected",
    [
        ([], []),
        ([0, 0, 0], []),
        ([1, 1, 1, 1, 1], [(0, 4)]),
        ([0, 1, 1, 1, 0, 1, 0, 1], [(1, 3), (5, 5), (7, 7)]),
    ],
)
@pytest.mark.parametrize("column", [False, True])
def test_signal_block_boundaries(signal, expected, column):
    """Keep runs maximal when they cross scan blocks or reach the final frame."""
    from reacnetgenerator._timedoutput import _signal_ranges

    values = np.array(signal)
    if column:
        values = values.reshape(-1, 1)
    assert list(_signal_ranges(values, block_rows=2)) == expected
