import numpy as np
cimport numpy as np

DTYPE = np.int16
ctypedef np.int16_t DTYPE_t
cimport cython
@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True)
def revert_preconditioning_inner(np.ndarray[DTYPE_t, ndim=1] d):
    cdef int N = d.shape[0]
    cdef int i
    for i in range(2, N):
        d[i] += (d[i-1] + d[i-2]) / 2
    return d