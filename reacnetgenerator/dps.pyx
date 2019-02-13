# distutils: language = c++
# cython: language_level=3
"""Connect molecule with Depth-First Search."""
import numpy as np
cimport numpy as np

from .stack import Stack, StackEmpty

def dps(bonds, levels):
    molecule = []
    bondlist = []
    cdef int _N = len(bonds)
    cdef np.ndarray[np.npy_bool, ndim=1, cast=True] visited
    visited = np.zeros((_N), dtype=bool)
    cdef int i, s, b, l
    st = Stack()

    for i in range(_N):
        if not visited[i]:
            mol = []
            bond = []
            st.push(i)
            visited[i] = True
            try:
                while True:
                    s = st.pop()
                    mol.append(s)
                    for b, l in zip(bonds[s], levels[s]):
                        if not visited[b]:
                            bond.append((s, b, l) if i < b else (b, s, l))
                            st.push(b)
                            visited[b]=True
            except StackEmpty:
                pass
            mol.sort()
            bond.sort()
            molecule.append(mol)
            bondlist.append(bond)
    return molecule, bondlist