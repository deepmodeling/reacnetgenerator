# SPDX-License-Identifier: LGPL-3.0-or-later
# Copyright 2018-2022, East China Normal University
"""ReacNetGenerator is an automatic reaction network generator for
reactive molecular dynamics simulation.[1]_.

References
----------
.. [1] Jinzhe Zeng, Liqun Cao, Chih-Hao Chin, Haisheng Ren, John Z. H.
   Zhang, Tong Zhu, ReacNetGenerator: an automatic reaction network
   generator for reactive molecular dynamic simulations, Phys. Chem.
   Chem. Phys., 2020, 22 (2): 683-691, doi: 10.1039/C9CP05091D.
"""

__date__ = "2018-03-11"
__author__ = "Jinzhe Zeng"
__email__ = "jinzhe.zeng@ustc.edu.cn"
__credits__ = ["Jinzhe Zeng", "Tong Zhu", "Liqun Cao", "Chih-Hao Chin", "John ZH Zhang"]
__copyright__ = (
    "Copyright 2018-2024, East China Normal University; Copyright 2024, DeepModeling"
)

from typing import TYPE_CHECKING

from ._version import __version__


class ReacNetGenerator:
    """Factory class for :class:`reacnetgenerator.reacnetgen.ReacNetGenerator`."""

    def __new__(cls, *args, **kwargs):
        """Create a new ReacNetGenerator instance."""
        from .reacnetgen import ReacNetGenerator as RealRNG

        return RealRNG(*args, **kwargs)


if TYPE_CHECKING:
    from .reacnetgen import ReacNetGenerator

def run(
    *,
    input_path,
    output_dir,
    input_type,
    atomname,
    items=("species", "reactions", "network", "report"),
    **kwargs,
):
    """Run ReacNetGenerator and return semantic artifact paths.

    Parameters
    ----------
    input_path : str or pathlib.Path or sequence
        Input trajectory or bond file(s).
    output_dir : str or pathlib.Path
        Directory receiving all default-generated artifacts.
    input_type : str
        ReacNetGenerator input type, such as "dump" or "bond".
    atomname : sequence of str
        Element names in the input trajectory.
    items : sequence of str, optional
        Requested stages. "species" and "reactions" are produced by
        the core run; "network" and "report" control the optional
        drawing/report stages.
    **kwargs
        Additional arguments forwarded to ReacNetGenerator.
    """
    if isinstance(items, str):
        items = (items,)
    requested = set(items)
    allowed = {"species", "reactions", "network", "report"}
    unknown = requested - allowed
    if unknown:
        raise ValueError(f"Unsupported output items: {sorted(unknown)}")
    if not requested:
        raise ValueError("items must contain at least one output stage")

    generator = ReacNetGenerator(
        inputfilename=input_path,
        output_dir=output_dir,
        inputfiletype=input_type,
        atomname=atomname,
        **kwargs,
    )
    return generator.runanddraw(
        run=True,
        draw="network" in requested,
        report="report" in requested,
    )


__all__ = ["ReacNetGenerator", "__version__", "run"]
