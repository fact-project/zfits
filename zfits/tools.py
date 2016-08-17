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



