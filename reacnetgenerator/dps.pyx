# distutils: language = c++
# cython: language_level=3
# cython: linetrace=True
# cython: infer_types=True
"""Connect molecule with Depth-First Search."""
from libc.stdlib cimport malloc, free
from collections import defaultdict

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

def dps_reaction(reactdict):
    """A+B->C+D"""
    # left in index 0 and right in index 1
    reactions = []
    cdef set visited_left = set()
    cdef set visited_right = set()
    visited = [visited_left, visited_right]
    cdef C_Stack st

    for init_mol in reactdict[0]:
        if init_mol not in visited[0]:
            reaction = [[], []]
            st.push(init_mol)
            st.push(0)
            while True:
                side, mol = st.pop(), st.pop()
                if mol < 0:
                    break
                elif mol in visited[side]:
                    continue
                reaction[side].append(mol)
                for r in reactdict[side][mol]:
                    if r < 0:
                        if r not in reaction[1-side]:
                            reaction[1-side].append(r)
                    elif r not in visited[1-side]:
                        st.push(r)
                        st.push(1-side)
                visited[side].add(mol)
            reactions.append(reaction)
    return reactions
