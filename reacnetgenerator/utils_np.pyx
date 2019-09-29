# distutils: language = c
# cython: language_level=3
# cython: linetrace=True
# cython: infer_types=True

import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.int
ctypedef np.int_t DTYPE_t
DTYPE8 = np.int8
ctypedef np.int8_t DTYPE8_t

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef idx_to_signal(DTYPE_t[:] idx, int step):
    """Benchmark for 1,000,000 loops (step=250000):
    Cython: 18.61 s
    Python: 31.30 s
    """
    signal = np.zeros((step, 1), dtype=DTYPE8, order='F')
    cdef int n = idx.size
    cdef DTYPE8_t[::1, :] signal_view = signal
    cdef DTYPE_t[:] idx_view = idx
    cdef int i
    with nogil:
        for i in range(n):
            signal_view[idx_view[i]][0] = 1
    return signal

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef check_zero_signal(DTYPE_t[:] signal):
    """Benchmark for 1,000,000 loops (1000/2000):
    Cython: 1.45 s
    Python: 3.67 s
    """
    cdef int n = signal.size
    cdef int i
    for i in range(n):
        if signal[i] == 1:
            return True
    return False
    
