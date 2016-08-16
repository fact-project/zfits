from fitsio import FITS

import tempfile, shutil, os
import numpy as np
import io
import bitarray as ba
import struct


def create_temporary_copy(path):
    temp_dir = tempfile.gettempdir()
    temp_path = os.path.join(temp_dir, os.path.basename(path))
    shutil.copy2(path, temp_path)
    return temp_path

def modify_copy_THEAP(path):
    temp_path = create_temporary_copy(path)
    f = FITS(temp_path, "rw")
    f[2].write_key("THEAP", f[2].read_header()["ZHEAPPTR"])
    f.close()
    return temp_path

def read_hufftree(odr):
    count = unpack(odr, "Q")[0]
    syms = {}
    for sym_id in range(count):
        sym = unpack(odr, "h")[0]
        numbits = unpack(odr, "B")[0]
        numbytes = np.ceil(numbits/8).astype(np.uint8)
        bits = unpack(odr, "{0:d}B".format(numbytes))

        numbits = numbits%8
        if numbits == 0:
            numbits = 8

        syms[sym] = bits, numbits

    hufftree = {}
    for sym in syms:
        code, numbits = syms[sym]

        sub_tree = hufftree
        for byte in code[:-1]:
            sub_tree = sub_tree.setdefault(byte, {})

        for i in np.arange(2**(8-numbits), dtype=np.uint8) << numbits:
            sub_tree[code[-1] | i] = sym, numbits

    return hufftree

def uncompress_huffman(stream):
    compressedSizes = unpack(stream, "I")[0]
    data_count = unpack(stream, "Q")[0]
    
    t = read_hufftree(stream)
    
    a = ba.bitarray()
    cur_tree = t
    found_syms = []
    while len(found_syms) < data_count:
        if len(a) < 8:
            for jj in range(4):
                xr = stream.read(1)
                aa = ba.bitarray()
                aa.frombytes(xr)
                a = aa + a

        N = (8 - len(a) % 8)
        aa = ba.bitarray('0'*N)
        b = (aa +a).tobytes()[-1]

        result = cur_tree[b]
        if not isinstance(result, tuple):
            cur_tree = result
            for i in range(8):
                a.pop()
            continue

        sym, nbits = result
        cur_tree = t
        found_syms.append(sym)

        for i in range(nbits):
            a.pop()
    return found_syms

def unpack(stream, fmt):
    size = struct.calcsize(fmt)
    buf = stream.read(size)
    return struct.unpack(fmt, buf)

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
    0: lambda x : x,
    1: uncompress_huffman,
    2: lambda x : x,
}

class ZFits(FITS):

    def __init__(self, filename):
        self.orig_path = filename
        self.temp_path = modify_copy_THEAP(filename)
        super().__init__(self.temp_path, mode='r')

        self.dtypes = self.learn()
        self.dtypes["BoardTime"] = (self.dtypes["BoardTime"][0], 'u4')
        self.dtypes["Data"] = (self.dtypes["Data"][0], 'B')

    def __del__(self):
        import os
        os.unlink(self.temp_path)

    def learn(self):
        colnames = self[2].get_colnames()
        h = self[2].read_header()
        uncompressed_coltypes = {}
        for colnum, colname in enumerate(colnames):
            x = str(h["ZFORM{0:d}".format(colnum+1)])
            length = int(x[:-1])
            np_type = fits_to_np_map[x[-1].upper()][2]
            uncompressed_coltypes[colname] = (length, np_type)
        return uncompressed_coltypes

    def get(self, colname, arg):
        x = self[2][colname][arg]
        return self.compression_block(x, colname)


    def compression_block(self, x, colname):
        stream = io.BytesIO(x.tobytes())
        length = unpack(stream, 'q')[0]
        ordering = unpack(stream, 'c')[0]
        num_proc = unpack(stream, 'B')[0]
        proc = unpack(stream, "%dH"%num_proc)
        for p in proc:
            stream = process_raw_data[p](stream)
        return self.convert(stream, colname)

    def convert(self, stream, colname):
        if colname != "Data":
            dtype = self.dtypes[colname][1]
            return np.frombuffer(stream.read(), dtype)
        else:
            return stream



def main():
    z = ZFits("20160809_121.fits.fz")
    x = z.get("Data", 0)
    return x
    
if __name__ == "__main__":
    x = main()

    assert x[:10] == [-72, -10, -70, 37, 68, 9, 1, -3, 2, 3]
    assert x[295:305] == [0, -5, 3, -2, -2, -79, 5, -38, 40, 53]
    assert x[-10:] == [-1, -2, 0, 1, 2, 4, 3, 1, -2, -2]