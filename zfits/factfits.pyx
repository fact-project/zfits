# distutils: language = c++
from fitsio import FITS
import numpy as np
cimport numpy as np
from libcpp.string cimport string
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector
from collections import namedtuple

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
            'B': np.uint8,
        }

        dtypes = {}
        for name, type_code, width in zip(
            self.c_factfits.fTable.GetColumnNames(),
            list(map(chr, self.c_factfits.fTable.GetColumnTypes())),
            self.c_factfits.fTable.GetColumnWidth()
        ):
            dtypes[name] = numpy_type_map[type_code], width

        return dtypes

    def SetPtrAddress_uint8(self, name):
        dtype, width = self.cols_dtypes[name]
        assert dtype == np.uint8, "Must be uint8"

        cdef np.ndarray[np.uint8_t] _array = np.zeros(width, dtype=np.uint8)

        self.c_factfits.SetPtrAddress(
            name,
            <np.uint8_t*>_array.data,
            _array.shape[0]
        )

        return _array


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
        self.fits = FITS(fname)

        header = self.header()
        if 'ZTABLE' in header and header['ZTABLE']:
            self.zfits = True
        else:
            self.zfits = False

        if self.zfits:
            self.fact_fits = Pyfactfits(fname)

        self.row = 0

        if self.zfits:
            self.rows = self.fact_fits.GetNumRows()
        else:
            self.rows = self.fits['Events'].get_nrows()

        if self.zfits:
            set_ptr_address = {
                np.int16: self.fact_fits.SetPtrAddress_int16,
                np.int32: self.fact_fits.SetPtrAddress_int32,
                np.uint8: self.fact_fits.SetPtrAddress_uint8,
            }

            self.data = {}
            for name, (dtype, width) in self.fact_fits.cols_dtypes.items():
                self.data[name] = set_ptr_address[dtype](name)

            colnames = [
                c.decode('utf-8')
                for c in self.fact_fits.cols_dtypes.keys()
            ]

        else:
            colnames = self.fits['Events'].get_colnames()

    def header(self):
        return self.fits['Events'].read_header()

    def __next__(self):
        if self.row >= self.rows:
            raise StopIteration


        evt_dict = {}

        if self.zfits:
            self.fact_fits.GetRow(self.row)
            for k, v in self.data.items():
                key = k.decode('utf-8')
                value = v.copy()
                evt_dict[key] = value
        else:
            data = self.fits['Events'][self.row]
            for column in data.dtype.names:
                evt_dict[column] = data[column][0]

        self.row += 1

        evt_dict['Data'].shape = (1440, -1)
        return evt_dict

    def __iter__(self):
        return self

