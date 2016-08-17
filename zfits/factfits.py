import os
import io
import shutil
import tempfile
import struct
from math import ceil
import numpy as np
from fitsio import FITS
from fitsio import FITS
from .zfits import ZFits

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
