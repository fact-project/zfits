# distutils: language = c++
from fitsio import FITS
import numpy as np
cimport numpy as np
cimport cython
from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector

# maybe nice to know ... not needed at the moment.
fits_to_np_map = {
    "L": ("Logical", 1),
    "B": ("Unsigned byte", 1, 'u1'),
    "I": ("16-bit integer", 2, 'i2'),
    "J": ("32-bit integer", 4, 'i4'),
    "K": ("64-bit integer", 8, 'i8'),
    "A": ("Character", 1, 'c'),
    "E": ("Single precision floating point", 4, 'f4'),
    "D": ("Double precision floating point", 8, 'f8'),
    "C": ("Single precision complex", 8),
    "M": ("Double precision complex", 16),
    "P": ("Array Descriptor (32-bit)", 8),
    "Q": ("Array Descriptor (64-bit)", 16),
}


cdef extern from "factfits.h":
    cdef cppclass factfits:
        cppclass Table:
            vector[string] GetColumnNames()
            vector[char] GetColumnTypes()
            vector[size_t] GetColumnWidth()

        Table fTable

        factfits(
            const string fname,
            const string tableName,
            bool_t force
        ) except +

        size_t GetNumRows() except +

        bool_t SetPtrAddress[T](
            const string name,
            T* ptr,
            size_t cnt)

        bool_t GetRow(size_t row, bool_t check)

cdef class Pyfactfits:
    cdef factfits* c_factfits

    def __cinit__(self, fname, tablename="", force=False):
        self.c_factfits = new factfits(
            bytes(fname, 'ascii'),
            bytes(tablename, 'ascii'),
            force)

    def __dealloc__(self):
        del self.c_factfits

    def GetRow(self, row, check=True):
        return self.c_factfits.GetRow(row, check)

    def GetNumRows(self):
        return self.c_factfits.GetNumRows()

    @property
    def cols_dtypes(self):

        numpy_type_map = {
            'I': np.int16,
            'J': np.int32,
        }

        dtypes = {}
        for name, type_code, width in zip(
            self.c_factfits.fTable.GetColumnNames(),
            list(map(chr, self.c_factfits.fTable.GetColumnTypes())),
            self.c_factfits.fTable.GetColumnWidth()
        ):
            dtypes[name] = numpy_type_map[type_code], width

        return dtypes

    def SetPtrAddress_int16(self, name):
        dtype, width = self.cols_dtypes[name]
        assert dtype == np.int16, "Must be int16"

        cdef np.ndarray[np.int16_t] _array = np.zeros(width, dtype=np.int16)

        self.c_factfits.SetPtrAddress(
            name,
            <np.int16_t*>_array.data,
            _array.shape[0]
        )

        return _array

    def SetPtrAddress_int32(self, name):
        dtype, width = self.cols_dtypes[name]
        assert dtype == np.int32, "Must be int32"

        cdef np.ndarray[np.int32_t] _array = np.zeros(width, dtype=np.int32)

        self.c_factfits.SetPtrAddress(
            name,
            <np.int32_t*>_array.data,
            _array.shape[0]
        )

        return _array


class FactFits:

    def __init__(self, fname):
        self.f = Pyfactfits(fname)
        self.fitsio = FITS(fname)
        self.row = 0
        self.rows = self.f.GetNumRows()

        self.data = {}
        set_ptr_address = {
            np.int16: self.f.SetPtrAddress_int16,
            np.int32: self.f.SetPtrAddress_int32,
        }
        for name, (dtype, width) in self.f.cols_dtypes.items():
            self.data[name] = set_ptr_address[dtype](name)

    def header(self):
        return self.fitsio[2].read_header()

    def __next__(self):
        if self.row < self.rows:
            self.f.GetRow(self.row)
            self.row += 1
            return {k: v.copy() for k, v in self.data.items()}
        else:
            raise StopIteration

    def __iter__(self):
        return self

@cython.boundscheck(False)
@cython.wraparound(False)
def remove_spikes_4(
        np.ndarray[np.float32_t, ndim=2] calib_data not None,
        ):

    cdef int number_of_pixel = calib_data.shape[0]
    cdef unsigned int roi = calib_data.shape[1]

    remove_spikes_4_dom(&calib_data[0,0], number_of_pixel, roi)
    return None
