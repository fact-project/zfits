from fitsio import FITS

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
