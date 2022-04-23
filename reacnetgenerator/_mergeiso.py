from tqdm import tqdm
import numpy as np

from .utils import bytestolist, listtobytes,  SharedRNGData


class _mergeISO(SharedRNGData):
    def __init__(self, rng):
        SharedRNGData.__init__(
            self, rng, ['miso', 'moleculetempfilename'], ['temp1it'])

    def merge(self):
        if self.miso > 0:
            self._mergeISO()
        self.returnkeys()

    def _mergeISO(self):
        with open(self.moleculetempfilename, 'rb') as ft:
            items = ft.readlines()
        items = [[items[i], items[i+1], items[i+2]]
                 for i in range(0, len(items), 3)]
        new_items = []
        _oldbitem = b''
        _oldbbond = b'0'
        _oldindex = []
        _oldfreq = 0
        for _bitem, _bbond, _bindex in tqdm(sorted(items, key=lambda x: (x[0], x[1].split()[0])), desc='Merge isomers:'):
            _index = bytestolist(_bindex)
            _freq = len(_index)
            if (_bitem == _oldbitem) and ((_bbond.split()[0] == _oldbbond.split()[0]) or (self.miso > 1)):
                _oldindex = np.hstack((_oldindex, _index))
            else:
                if _oldbitem:
                    new_items.append(
                        [_oldbitem, _oldbbond, listtobytes(_oldindex)])
                _oldbitem = _bitem
                _oldindex = _index
                if _freq > _oldfreq:
                    _oldbbond = _bbond
        new_items.append([_oldbitem, _oldbbond, listtobytes(_oldindex)])
        new_items.sort(key=lambda x: len(x[0]))
        self.temp1it = len(new_items)
        with open(self.moleculetempfilename, 'wb') as ft:
            for item in new_items:
                [ft.write(i) for i in item]
