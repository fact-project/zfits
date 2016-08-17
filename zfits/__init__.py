import os
import io
import shutil
import tempfile
import struct
from math import ceil
import numpy as np
from fitsio import FITS

def unpack(stream, fmt):
    size = struct.calcsize(fmt)
    buf = stream.read(size)
    return struct.unpack(fmt, buf)

def create_temporary_copy(path):
    temp_dir = tempfile.gettempdir()
    temp_path = os.path.join(temp_dir, os.path.basename(path))
    shutil.copy2(path, temp_path)
    return temp_path

def modify_copy_THEAP(path):
    temp_path = create_temporary_copy(path)
    f = FITS(temp_path, "rw")
    f[2].write_key("THEAP", f[2].read_header()["ZHEAPPTR"])
    f[1].write_key("THEAP", f[1].read_header()["ZHEAPPTR"])
    f.close()
    return temp_path

def read_hufftree(stream):
    number_of_symbols = unpack(stream, "Q")[0]

    hufftree = {}
    for symbol_id in range(number_of_symbols):
        sym = unpack(stream, "h")[0]
        numbits = unpack(stream, "B")[0]
        numbytes = ceil(numbits/8)
        numbits = numbits % 8 if not numbits % 8 == 0 else 8
        code = unpack(stream, "{0:d}B".format(numbytes))

        sub_tree = hufftree
        for byte in code[:-1]:
            sub_tree = sub_tree.setdefault(byte, {})

        for i in np.arange(2**(8-numbits), dtype=np.uint8) << numbits:
            sub_tree[code[-1] | i] = sym, numbits

    return hufftree

def uncompress_huffman(stream, *args):
    compressedSizes = unpack(stream, "I")[0]
    data_count = unpack(stream, "Q")[0]
    
    hufftree = read_hufftree(stream)
    
    cur_tree = hufftree
    
    reservoir = 0
    fill = 0
    nbits = 8
    found_symbols = np.zeros(data_count, dtype=np.int16)
    symbol_id = 0
    while symbol_id < data_count:
        if fill < 8:
            reservoir |= unpack(stream, "B")[0] << fill
            fill += 8

        result = cur_tree[0xff & reservoir]
        if not isinstance(result, tuple):
            cur_tree = result
            nbits = 8
        else:
            sym, nbits = result
            found_symbols[symbol_id] = sym
            symbol_id += 1
            cur_tree = hufftree
        reservoir >>= nbits
        fill -= nbits
    return found_symbols

def revert_preconditioning(stream, *args):
    d = stream
    for i in range(2, len(d)):
        d[i] += int((d[i-1] + d[i-2])/2)

    return d

def convert(stream, dtype):
    return np.frombuffer(stream.read(), dtype)


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

process_raw_data = {
    0: convert,
    1: revert_preconditioning,
    2: uncompress_huffman,
}

class ZFits(FITS):

    def __init__(self, filename):
        self.orig_path = filename
        self.temp_path = modify_copy_THEAP(filename)
        super().__init__(self.temp_path, mode='r')

    def __del__(self):
        import os
        os.unlink(self.temp_path)

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
            stream = process_raw_data[proc_key](stream, dtype)
        return stream

    def get(self, extension, colname, rownum):
        dtype = self._read_dtype(extension, colname)
        array = self[extension][colname][rownum]
        return self._uncompress_block(array, colname, dtype)

class FactFits:

    def __init__(self, data_path, calib_path):
        self.data_file = ZFits(data_path)
        self.drs_file = FITS(calib_path)

        z_drs_offset = self.data_file.get("ZDrsCellOffsets","OffsetCalibration",0)
        z_drs_offset = z_drs_offset.reshape(1440, -1)
        z_drs_offset = np.concatenate((z_drs_offset, z_drs_offset), axis=1)
        self.z_drs_offset = z_drs_offset

        bsl = self.drs_file[1]["BaselineMean"][0]
        bsl = bsl.reshape(1440, -1)
        bsl = np.concatenate((bsl, bsl), axis=1)
        bsl *= 4096 / 2000
        self.bsl = bsl

        self.off = self.z_drs_offset - self.bsl

        gain = self.drs_file[1]["GainMean"][0]
        gain = gain.reshape(1440, -1)
        gain = np.concatenate((gain, gain), axis=1)
        gain *= 1 # some ominous factor missing here 1.7 or so.
        self.gain = gain

        
        trg = self.drs_file[1]["TriggerOffsetMean"][0]
        trg = trg.reshape(1440, -1)
        trg *= 4096 / 2000
        self.trg = trg


    def get(self, colname, row):
        data = self.data_file.get("Events", colname, row)

        return data 

    def get_data_calibrated(self, row):

        data = self.data_file.get("Events", "Data", row)        
        sc = self.data_file.get("Events", "StartCellData", row)
        data = data.reshape(1440, -1)

        calib_data = np.zeros_like(data, np.float32)
        
        for i in range(1440):
            calib_data[i] = data[i] + self.off[i, sc[i]:sc[i]+len(calib_data[i])] - self.trg[i]
        
        return calib_data


    def __repr__(self):
        return repr(self.data_file[2]) + repr(self.drs_file[1])
