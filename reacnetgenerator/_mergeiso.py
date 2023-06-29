# SPDX-License-Identifier: LGPL-3.0-or-later
import itertools

import numpy as np
from tqdm.auto import tqdm

from .utils import SharedRNGData, bytestolist, listtobytes, read_compressed_block


class _mergeISO(SharedRNGData):
    miso: int
    moleculetempfilename: str
    temp1it: int

    def __init__(self, rng):
        SharedRNGData.__init__(self, rng, ["miso", "moleculetempfilename"], ["temp1it"])

    def merge(self):
        if self.miso > 0:
            self._mergeISO()
        self.returnkeys()

    def _mergeISO(self):
        items = []
        with open(self.moleculetempfilename, "rb") as ft:
            for item in itertools.zip_longest(*[read_compressed_block(ft)] * 4):
                items.append(item)
        new_items = []
        _oldbitem = b""
        _oldbbond = b"0"
        _oldindex = []
        _oldfreq = 0
        for _bitem, _bbond0, _bbond1, _bindex in tqdm(
            sorted(items, key=lambda x: (x[0], x[1])),
            desc="Merge isomers:",
            disable=None,
        ):
            _index = bytestolist(_bindex)
            _freq = len(_index)
            if (_bitem == _oldbitem) and ((_bbond0 == _oldbbond[0]) or (self.miso > 1)):
                _oldindex = np.hstack((_oldindex, _index))
            else:
                if _oldbitem:
                    new_items.append([_oldbitem, *_oldbbond, listtobytes(_oldindex)])
                _oldbitem = _bitem
                _oldindex = _index
                if _freq > _oldfreq:
                    _oldbbond = (_bbond0, _bbond1)
        new_items.append([_oldbitem, *_oldbbond, listtobytes(_oldindex)])
        new_items.sort(key=lambda x: len(x[0]))
        self.temp1it = len(new_items)
        with open(self.moleculetempfilename, "wb") as ft:
            for item in new_items:
                [ft.write(i) for i in item]
