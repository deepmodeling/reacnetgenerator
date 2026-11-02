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
    BondChange,
    TransitionEvidence,
    TransitionParticipant,
    iter_frames,
    iter_molecule_ranges,
    iter_molecules,
    iter_reaction_events,
    iter_reaction_types,
    iter_species,
    iter_transition_evidence,
    read_metadata,
    validate_timed_output,
)
from reacnetgenerator.utils import bytestolist, idx_to_signal, read_compressed_block

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
    with open(rng.hmmfilename if run_hmm else rng.moleculetemp2filename, "rb") as file:
        blocks = read_compressed_block(file)
        signals = (
            (bytestolist(blob) for blob in blocks)
            if run_hmm
            else (
                idx_to_signal(bytestolist(record[-1]), rng.step)
                for record in zip(*[blocks] * 4, strict=True)
            )
        )
        for molecule_id, signal in enumerate(signals, 1):
            expected.extend((molecule_id, frame) for frame in np.flatnonzero(signal))
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
    evidence = list(iter_transition_evidence(path, block_rows=1))
    assert (
        Counter((item.transition, item.reactant, item.product) for item in evidence)
        == csv_events
    )
    assert all(isinstance(item, TransitionEvidence) for item in evidence)
    assert all(
        isinstance(participant, TransitionParticipant)
        for item in evidence
        for participant in item.participants
    )
    assert all(
        isinstance(change, BondChange)
        for item in evidence
        for change in item.bond_changes
    )
    metadata = read_metadata(path)
    assert metadata["schema_version"] == "1.1"
    assert metadata["capabilities"] == [
        "molecule_ranges",
        "reaction_events",
        "transition_evidence",
    ]
    assert metadata["configuration"]["parameters"]["runHMM"] is run_hmm
    assert validate_timed_output(path, block_rows=2).frames == rng.step
    assert not list(tmp_path.glob("*.incomplete"))


def test_optional_and_existing_outputs(tmp_path):
    """An explicit timeline preserves text results and does not enable timed CSV."""
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


@pytest.mark.parametrize(
    ("options", "legacy_output"),
    [
        ({"printmoleculetime": True}, "moleculetimelinefilename"),
        ({"printreactionevent": True}, "reactioneventfilename"),
        ({"moleculeframes": [0]}, "moleculetimelinefilename"),
        ({"moleculetimesteps": [0]}, "moleculetimelinefilename"),
    ],
)
@pytest.mark.parametrize("explicit_path", [False, True])
def test_timed_flags_add_default_hdf5_companion(
    tmp_path, options, legacy_output, explicit_path
):
    """Existing timed switches also publish the normalized HDF5 timeline."""
    directory = tmp_path / "outputs"
    path = tmp_path / "custom.h5" if explicit_path else directory / "timeline.h5"
    rng = generator(directory, timed_output=path if explicit_path else None, **options)

    assert rng.timed_output == str(path)
    assert rng.artifacts["timeline"] == str(path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()

    assert path.exists()
    legacy_path = Path(getattr(rng, legacy_output))
    assert legacy_path.exists()
    if "moleculeframes" in options or "moleculetimesteps" in options:
        with legacy_path.open(newline="") as file:
            assert {row["Timestep"] for row in csv.DictReader(file)} == {"0"}
    assert validate_timed_output(path, block_rows=2).frames == rng.step
    if explicit_path:
        assert not (directory / "timeline.h5").exists()


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
    assert validate_timed_output(path, block_rows=2).frames == rng.step


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
    assert validate_timed_output(path, block_rows=2).frames == rng.step


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


def test_cached_writer_still_rejects_absent_participant(tmp_path):
    """A warmed cache must reject a gap and preserve the previous complete file."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    path.write_bytes(b"previous complete result")
    with pytest.raises(ValueError, match="participant is absent"):
        with _TimedOutputWriter(rng) as writer:
            for sid, name in enumerate(("A", "B")):
                writer._append("species", name=name)
                writer.species[name] = sid
                writer._append(
                    "molecules",
                    species_id=sid,
                    atom_index=0,
                    atom_offsets=sid + 1,
                    bond_offsets=0,
                )
                writer._append(
                    "molecule_ranges",
                    molecule_id=sid + 1,
                    start_frame=sid,
                    end_frame=sid,
                )
            event = {
                "Timestep_Index": 0,
                "Reactant": "A",
                "Product": "B",
                "ReactantMoleculeIDs": (1,),
                "ProductMoleculeIDs": (2,),
            }
            writer.write_events([event])
            assert writer._range_lookup is not None
            writer.write_events([dict(event, Timestep_Index=1)])
    assert writer._range_lookup is None
    assert writer._species_names is None
    assert not writer._molecule_cache
    assert writer._molecule_cache_bytes == 0
    assert path.read_bytes() == b"previous complete result"
    (partial,) = tmp_path.glob("*.incomplete")
    with h5py.File(partial, "r+") as file:
        assert file.attrs["status"] == "incomplete"


def test_species_id_lookup_does_not_rescan_the_dictionary(tmp_path):
    """Dictionary order is not ID order, and warm events must not scan names."""

    class Once(dict):
        scans = 0

        def items(self):
            self.scans += 1
            if self.scans > 1:
                raise AssertionError("species dictionary scanned for every event")
            return super().items()

    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        writer.species = Once({"B": 1, "A": 0})
        writer._append("molecules", species_id=0)
        writer._append("molecules", species_id=1)
        assert writer._participant_species([1, 2, 1]) == {"A": 2, "B": 1}
        assert writer._participant_species([2]) == {"B": 1}
        assert writer.species.scans == 1
    assert writer._species_names is None


@pytest.mark.parametrize("species", [{"A": -1}, {"A": 1}, {"A": 0, "B": 0}])
def test_species_index_rejects_invalid_dictionary_without_caching_it(tmp_path, species):
    """Reject out-of-range or duplicate IDs before publishing a name index."""
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        writer.species = species
        writer._append("molecules", species_id=0)
        with pytest.raises(ValueError, match="invalid species_id"):
            writer._participant_species([1])
        assert writer._species_names is None
        writer.species = {"A": 0}
        assert writer._participant_species([1]) == {"A": 1}


@pytest.mark.parametrize("invalid_id", [-1, 2])
def test_warm_species_index_rejects_invalid_ids(tmp_path, invalid_id):
    """A warm name index must still reject invalid stored species IDs."""
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        writer.species = {"A": 0}
        writer._append("molecules", species_id=0)
        assert writer._participant_species([1]) == {"A": 1}
        writer.columns["molecules/species_id"].dataset[0] = invalid_id
        with pytest.raises(ValueError, match="invalid species_id"):
            writer._participant_species([1])


def _append_metadata_molecules(writer, *, atom_start=0, empty=False):
    for row in range(3):
        writer._append("molecules", species_id=0, bond_offsets=row + 1)
        writer._append("molecules", atom_offsets=0 if empty else row + 1)
        if not empty:
            writer._append("molecules", atom_index=atom_start + row)
        writer._append(
            "molecules", bond_atom_index_1=0, bond_atom_index_2=1, bond_order=row + 1
        )


def test_warm_metadata_avoids_reads_and_rechecks_overlap_and_bonds(
    tmp_path, monkeypatch
):
    """Reuse stored definitions while preserving per-event consistency checks."""
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        _append_metadata_molecules(writer)
        assert writer._participant_atoms([1, 2]) == {0, 1}
        assert writer._participant_bonds([1]) == {(0, 1): 1}
        assert writer._participant_bonds([2]) == {(0, 1): 2}

        def no_read(*args):
            raise AssertionError("warm molecule metadata reread HDF5")

        with monkeypatch.context() as context:
            context.setattr(h5py.Dataset, "__getitem__", no_read)
            assert writer._participant_atoms([2, 1]) == {0, 1}
            assert writer._participant_bonds([1]) == {(0, 1): 1}
            with pytest.raises(ValueError, match="overlap"):
                writer._participant_atoms([1, 1])
            with pytest.raises(ValueError, match="conflicting bonds"):
                writer._participant_bonds([1, 2])
    assert not writer._molecule_cache
    assert writer._molecule_cache_bytes == 0


@pytest.mark.parametrize("empty", [False, True])
def test_metadata_cache_is_bounded_and_oversized_values_bypass_it(
    tmp_path, monkeypatch, empty
):
    """Bound bytes and empty entries without retaining an oversized molecule."""
    monkeypatch.setattr("reacnetgenerator._timedoutput._MOLECULE_CACHE_BYTES", 16)
    monkeypatch.setattr("reacnetgenerator._timedoutput._MOLECULE_CACHE_ENTRIES", 2)
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        _append_metadata_molecules(writer, empty=empty)
        writer._flush_molecule_tables()
        for mid in (1, 2, 1, 3):
            values = writer._molecule_values(mid, "atoms")
            assert not values.flags.writeable
            assert writer._molecule_cache_bytes <= 16
            assert len(writer._molecule_cache) <= 2
        assert list(writer._molecule_cache) == [("atoms", 1), ("atoms", 3)]
        assert writer._molecule_values(1, "bonds").tolist() == [[0, 1, 1]]
        assert ("bonds", 1) not in writer._molecule_cache
        assert writer._molecule_cache_bytes == sum(
            v.nbytes for v in writer._molecule_cache.values()
        )


def test_metadata_read_error_is_retryable_and_cache_is_file_local(
    tmp_path, monkeypatch
):
    """Failed reads do not poison later requests or leak into another file."""
    rng = generator(tmp_path, timed_output=tmp_path / "timeline.h5")
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        _append_metadata_molecules(writer)
        assert writer._participant_atoms([1]) == {0}
        original = h5py.Dataset.__getitem__

        def fail_atoms(dataset, key):
            if dataset.name.endswith("/atom_index"):
                raise OSError("atom read failed")
            return original(dataset, key)

        with monkeypatch.context() as context:
            context.setattr(h5py.Dataset, "__getitem__", fail_atoms)
            with pytest.raises(OSError, match="atom read failed"):
                writer._participant_atoms([2])
        assert ("atoms", 2) not in writer._molecule_cache
        assert writer._participant_atoms([2]) == {1}
    rng.timed_output = str(tmp_path / "other.h5")
    with _TimedOutputWriter(rng) as writer:
        _append_metadata_molecules(writer, atom_start=10)
        assert writer._participant_atoms([1]) == {10}


def test_compact_counts_and_blocks(tmp_path, monkeypatch):
    """Read the first block without loading the full table or expanding counts."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    with _TimedOutputWriter(rng) as writer:
        writer._append("species", name="A")
        writer._append("species", name="B")
        writer.species.update(A=0, B=1)
        for species_id in range(2):
            writer._append(
                "molecules",
                species_id=species_id,
                atom_index=0,
                atom_offsets=species_id + 1,
                bond_offsets=0,
            )
        for i in range(50):
            writer._append("molecule_ranges", molecule_id=1, start_frame=i, end_frame=i)
        writer._append("molecule_ranges", molecule_id=2, start_frame=1, end_frame=1)
        event = {
            "Timestep_Index": 0,
            "Reactant": "A",
            "Product": "B",
            "ReactantMoleculeIDs": (1,),
            "ProductMoleculeIDs": (2,),
        }
        writer.write_events([event] * 5)
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
    """Round-trip an explicit timeline and reject paths that could destroy inputs."""
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
    path.symlink_to(INPUT)
    with pytest.raises(ValueError, match="alias"):
        generator(tmp_path, printreactionevent=True)


def test_default_timeline_paths(tmp_path):
    """Timed switches and filters select a predictable default HDF5 path."""
    directory = tmp_path / "artifacts"
    assert generator(directory).timed_output is None
    assert (
        generator(directory, moleculeframes=[], moleculetimesteps=[]).timed_output
        is None
    )
    filtered = generator(directory, moleculeframes=[0])
    expected = str(directory / "timeline.h5")
    assert filtered.timed_output == expected
    assert filtered.artifacts["timeline"] == expected
    with pytest.raises(ValueError, match="executes PATH"):
        filtered.run_items(["species"])

    source = tmp_path / "trajectory.bond"
    rng = ReacNetGenerator(
        inputfilename=str(source),
        inputfiletype="bond",
        atomname=["H"],
        printreactionevent=True,
    )
    assert rng.timed_output == f"{source}.timeline.h5"


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
            "transition_evidence",
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


def test_species_reader_rejects_multidimensional_names(tmp_path):
    """Keep the public species name field scalar by rejecting array rows."""
    path = tmp_path / "timeline.h5"
    rng = generator(tmp_path, timed_output=path)
    prepare(rng)
    _CollectPaths.getstype(rng).collect()
    with h5py.File(path, "r+") as file:
        dataset = "species/name"
        values = file[dataset].asstr()[:].reshape(-1, 1)
        del file[dataset]
        file.create_dataset(dataset, data=values, dtype=h5py.string_dtype("utf-8"))

    with pytest.raises(ValueError, match="Invalid column shape"):
        next(iter_species(path))


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
