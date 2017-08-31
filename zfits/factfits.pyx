# distutils: language = c++
import numpy as np
cimport numpy as np

from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from libcpp.map cimport map as _map


cdef extern from "factfits.h":
    cdef cppclass factfits:
        cppclass Table:
            cppclass Column:
                size_t offset, num, size, bytes
                char type
                string unit
                int comp

            vector[string] GetColumnNames()
            vector[char] GetColumnTypes()
            vector[size_t] GetColumnWidth()

            ctypedef _map[string, Column] Columns

        Table fTable

        factfits(
            const string fname,
            const string tableName,
            bool_t force
        ) except +

        size_t GetNumRows() except +
        bool_t HasKey(const string key)
        bool_t HasColumn(const string col)

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


    def GetRow(self, row, check=True):
        return self.c_factfits.GetRow(row, check)

    def __dealloc__(self):
        del self.c_factfits

    def GetNumRows(self):
        return self.c_factfits.GetNumRows()

    def HasKey(self, key):
        return self.c_factfits.HasKey(bytes(key, 'ascii'))

    def HasColumn(self, col):
        return self.c_factfits.HasColumn(bytes(col, 'ascii'))

    def SetPtrAddress(self, name, array):
        print(array, array, array.dtype, array.data, array.shape)

    def ColumnNames(self):
        return self.c_factfits.fTable.GetColumnNames()

    def ColumnTypes(self):
        return list(map(chr, self.c_factfits.fTable.GetColumnTypes()))

    def ColumnWidth(self):
        return self.c_factfits.fTable.GetColumnWidth()

    @property
    def cols_dtypes(self):

        numpy_type_map = {
            'I': np.int16,
            'J': np.int32,
        }

        dtypes = {}
        for name, type_code, width in zip(
            self.ColumnNames(),
            self.ColumnTypes(),
            self.ColumnWidth()
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


class AnotherFoo:

    def __init__(self, fname, columns):
        self.f = Pyfactfits(fname)
        self.row = 0
        self.rows = self.f.GetNumRows()

        self.data = {}
        set_ptr_address = {
            np.int16: self.f.SetPtrAddress_int16,
            np.int32: self.f.SetPtrAddress_int32,
        }
        for name, (dtype, width) in self.f.cols_dtypes.items():
            if name not in columns:
                continue
            self.data[name] = set_ptr_address[dtype](name)

    def __next__(self):
        if self.row < self.rows:
            self.f.GetRow(self.row)
            self.row += 1
            return {k: v.copy() for k, v in self.data.items()}
        else:
            raise StopIteration

    def __iter__(self):
        return self

DTYPE = np.int16
ctypedef np.int16_t DTYPE_t
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def revert_preconditioning_inner(np.ndarray[DTYPE_t, ndim=1] d):
    cdef int N = d.shape[0]
    cdef int i
    for i in range(2, N):
        d[i] += (d[i-1] + d[i-2]) / 2
    return d

cdef extern void Decode_dom (const np.uint8_t *bufin, size_t bufinlen, np.int16_t *bufout, size_t bufoutlen)

@cython.boundscheck(False)
@cython.wraparound(False)
def huff_decode(np.ndarray[np.uint8_t, ndim=1] bufin not None, np.ndarray[np.int16_t, ndim=1] bufout not None):
    cdef int bufinlen = bufin.shape[0]
    cdef int bufoutlen = bufout.shape[0]
    Decode_dom(&bufin[0], bufinlen, &bufout[0], bufoutlen)
    return None


cdef extern void remove_spikes_4_dom (
        np.float32_t *calib_data,
        size_t number_of_pixel,
        np.uint32_t roi
    )

@cython.boundscheck(False)
@cython.wraparound(False)
def remove_spikes_4(
        np.ndarray[np.float32_t, ndim=2] calib_data not None,
        ):

    cdef int number_of_pixel = calib_data.shape[0]
    cdef unsigned int roi = calib_data.shape[1]

    remove_spikes_4_dom(&calib_data[0,0], number_of_pixel, roi)
    return None
