import io
from . import tools
from .tools import unpack
from fitsio import FITS
import functools
import numpy as np
import time

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

class ZFits(FITS):

    def __init__(self, filename):
        self.file_path = filename
        super().__init__(self.file_path, mode='r')

    def _read_dtype(self, extension, colname):
        colnum = self[extension].get_colnames().index(colname)
        h = self[extension].read_header()
        x = str(h["ZFORM{0:d}".format(colnum+1)])
        length = int(x[:-1])
        np_type = fits_to_np_map[x[-1].upper()][2]

        return np_type

    def _uncompress_block(self, array, colname, dtype):
        stream = io.BytesIO(array.tobytes())
        length = unpack(stream, 'q')[0]
        ordering = unpack(stream, 'c')[0]
        num_proc = unpack(stream, 'B')[0]
        proc = unpack(stream, "%dH"%num_proc)
        for proc_key in proc[::-1]:
            if proc_key == 0:
                stream = tools.convert(stream, dtype)
            elif proc_key == 1:
                stream = tools.revert_preconditioning(stream)
            elif proc_key == 2:
                stream = tools.uncompress_huffman(stream, array)
        return stream

    def get(self, extension, colname, rownum):
        dtype = self._read_dtype(extension, colname)
        array = self[extension][colname][rownum][0]
        return self._uncompress_block(array, colname, dtype)

    @property
    @functools.lru_cache()
    def z_drs_offset(self):
        z_drs_offset = self.get("ZDrsCellOffsets","OffsetCalibration",0)
        z_drs_offset = z_drs_offset.reshape(1440, -1)
        z_drs_offset = np.concatenate((z_drs_offset, z_drs_offset), axis=1)
        return z_drs_offset

    def get_raw_data(self, row, rolled=False):

        data = self.get("Events", "Data", row)
        data = data.reshape(1440, -1)
        if not (self.z_drs_offset == 0).all():
            for i in range(1440):
                data[i] += self.z_drs_offset[i, sc[i]:sc[i]+len(data[i])]

        if rolled:
            sc = self.get("Events", "StartCellData", row)
            for i, s in enumerate(sc):
                data[i] = np.roll(data[i], s)
            """
            assert data.shape[1] == 1024

            sc_per_chip = sc[::9]

            x = np.copy(data)
            for cid, scpc in enumerate(sc_per_chip):
                data[cid*9:(cid+1)*9, scpc:] = x[cid*9:(cid+1)*9, :1024 - scpc]
                data[cid*9:(cid+1)*9, :scpc] = x[cid*9:(cid+1)*9, 1024 - scpc:]
            """
        return data
