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

# declare the interface to the C code

cdef extern void c_multiply (double* array, double value, int m, int n)


@cython.boundscheck(False)
@cython.wraparound(False)
def multiply(np.ndarray[double, ndim=2, mode="c"] input not None, double value):
    """
    multiply (arr, value)

    Takes a numpy arry as input, and multiplies each elemetn by value, in place

    param: array -- a 2-d numpy array of np.float64
    param: value -- a number that will be multiplied by each element in the array

    """
    cdef int m, n

    m, n = input.shape[0], input.shape[1]

    c_multiply (&input[0,0], value, m, n)

    return None

cdef extern void Decode_dom (const np.uint8_t *bufin, size_t bufinlen, np.int16_t *bufout, size_t bufoutlen)

@cython.boundscheck(False)
@cython.wraparound(False)
def huff_decode(np.ndarray[np.uint8_t, ndim=1] bufin not None, np.ndarray[np.int16_t, ndim=1] bufout not None):
    cdef int bufinlen = bufin.shape[0]
    cdef int bufoutlen = bufout.shape[0]
    Decode_dom(&bufin[0], bufinlen, &bufout[0], bufoutlen)
    return None
