# distutils: language = c++
# cython: language_level=3
# cython: linetrace=True
# cython: infer_types=True
"""Connect molecule with Depth-First Search."""
from libc.stdlib cimport malloc, free

cdef extern from 'c_stack.h':
    # This function is copied from https://zhuanlan.zhihu.com/p/38212302
    cdef cppclass C_Stack:
        void push(int val);
        int pop();

def dps(bonds, levels):
    molecule = []
    bondlist = []
    cdef int _N = len(bonds)
    cdef int *visited = <int *> malloc(_N * sizeof(int))
    cdef int i, s, b_c, l
    cdef C_Stack st
    for i in range(_N):
        visited[i]=0

    for i in range(_N):
        if visited[i]==0:
            mol = []
            bond = []
            st.push(i)
            while True:
                s = st.pop()
                if s < 0:
                    break
                elif visited[s]==1:
                    continue
                mol.append(s)
                for b, l in zip(bonds[s], levels[s]):
                    b_c = b
                    if visited[b_c]==0:
                        bond.append((s, b, l) if i < b else (b, s, l))
                        st.push(b_c)
                visited[s]=1
            mol.sort()
            bond.sort()
            molecule.append(mol)
            bondlist.append(bond)
    free(visited)
    return molecule, bondlist
