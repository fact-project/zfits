# distutils: language = c++
# distutils: sources = remove_spikes_source.cpp
import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def remove_spikes_4(
        np.ndarray[np.float32_t, ndim=2] calib_data not None,
        ):

    cdef int number_of_pixel = calib_data.shape[0]
    cdef unsigned int roi = calib_data.shape[1]

    remove_spikes_4_dom(&calib_data[0, 0], number_of_pixel, roi)
    return None
