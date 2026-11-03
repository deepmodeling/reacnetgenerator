# SPDX-License-Identifier: LGPL-3.0-or-later
"""Classify reaction instances even when their formatted species sides coincide."""

import json
import shutil

import h5py
import pytest

from reacnetgenerator import ReacNetGenerator
from reacnetgenerator.reactioncounts import main
from reacnetgenerator.timedoutput import (
    _order_only_instances,
    iter_reaction_types,
    iter_transition_evidence,
    validate_timed_output,
)


@pytest.fixture(scope="module")
def oxygen_timeline(tmp_path_factory):
    """Three isolated O pairs, two of which exchange partners at each switch."""
    directory = tmp_path_factory.mktemp("oxygen")
    source = directory / "oxygen.xyz"
    initial = (0, 1.23529, 5, 6.23529, 10, 11.23529)
    switched = (0, 1.22136, 5, 11.23529, 10, 6.22136)
    with source.open("w") as stream:
        for positions in [initial] * 60 + [switched] * 60 + [initial] * 60:
            stream.write("6\noxygen pairs\n")
            for x in positions:
                stream.write(f"O {x} 0 0\n")
    path = directory / "timeline.h5"
    ReacNetGenerator(
        inputfilename=source,
        inputfiletype="xyz",
        pbc=False,
        atomname=["O"],
        output_dir=directory,
        timed_output=path,
        nproc=1,
    ).run()
    assert validate_timed_output(path).frames == 180
    return path


def test_select_instances_before_counting(oxygen_timeline, capsys):
    """At each transition a pure order change and a partner swap share a type."""
    totals = list(iter_reaction_types(oxygen_timeline))
    assert [row[-1] for row in totals] == [2, 2]
    for selection in ("only", "exclude"):
        evidence = list(
            iter_transition_evidence(
                oxygen_timeline, inferred_order=selection, block_rows=1
            )
        )
        assert [event.transition for event in evidence] == [59, 119]
        for event in evidence:
            kinds = {change.kind for change in event.bond_changes}
            if selection == "only":
                assert kinds == {"order_changed"}
                assert event.participants[0].atom_index == (0, 1)
            else:
                assert kinds == {"formed", "broken"}
        counts = list(
            iter_reaction_types(oxygen_timeline, inferred_order=selection, block_rows=1)
        )
        assert counts == [(*row[:3], 1) for row in totals]
        main([str(oxygen_timeline), "--inferred-order", selection])
        assert capsys.readouterr().out.splitlines() == [
            f"1 {left}->{right}" for _, left, right, _ in totals
        ]


@pytest.mark.parametrize("input_type", ["bond", "lammpsbondfile", "custom"])
def test_filter_requires_coordinate_provenance(oxygen_timeline, tmp_path, input_type):
    """Interpret the same recorded bond changes according to their declared source."""
    path = tmp_path / "copy.h5"
    shutil.copyfile(oxygen_timeline, path)
    with h5py.File(path, "r+") as file:
        config = json.loads(file.attrs["configuration"])
        config["parameters"]["inputfiletype"] = input_type
        file.attrs["configuration"] = json.dumps(config)
    if input_type == "custom":
        with pytest.raises(ValueError, match="bond-order source"):
            list(iter_reaction_types(path, inferred_order="only"))
        return
    assert list(iter_reaction_types(path, inferred_order="only")) == []
    assert list(iter_reaction_types(path, inferred_order="exclude")) == list(
        iter_reaction_types(path)
    )


def test_classify_complete_change_list(tmp_path):
    """Empty evidence and an instance with both an order change and a break stay in exclude."""
    with h5py.File(tmp_path / "changes.h5", "w") as file:
        evidence = file.create_group("transition_evidence")
        evidence["bond_change_offsets"] = [0, 0, 2, 4]
        evidence["before_order"] = [1, 2, 1, 3]
        evidence["after_order"] = [2, 0, 2, 1]
        assert list(_order_only_instances(file, block_rows=1)) == [False, False, True]
