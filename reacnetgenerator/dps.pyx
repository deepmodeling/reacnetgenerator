# distutils: language = c++
# cython: language_level=3
"""Connect molecule with Depth-First Search."""
def dps(bond, level):
    molecule = []
    bondlist = []
    _N = len(bond)
    done = [False]*_N
    def _mo(i):    
        molecule[-1].append(i)
        done[i-1] = True
        for b, l in zip(bond[i-1], level[i-1]):
            if not done[b-1]:
                bondlist[-1].append((i, b, l) if i < b else (b, i, l))
                _mo(b)
    
    for i in range(_N):
        if not done[i]:
            molecule.append([])
            bondlist.append([])
            _mo(i+1)
            molecule[-1].sort()
            bondlist[-1].sort()
    return molecule, bondlist