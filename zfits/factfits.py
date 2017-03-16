import numpy as np
from fitsio import FITS
from .zfits import ZFits

from .cython_tools import remove_spikes_4

class FactFits:

    def __init__(self, data_path, calib_path):
        self.data_file = ZFits(data_path)
        self.drs_file = FITS(calib_path)

        z_drs_offset = self.data_file.get(
            "ZDrsCellOffsets",
            "OffsetCalibration",
            0
        )
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
        gain /= 1907.35
        self.gain = gain


        trg = self.drs_file[1]["TriggerOffsetMean"][0]
        trg = trg.reshape(1440, -1)
        trg *= 4096 / 2000
        self.trg = trg

        self.previous_start_cells = []
        self.fMaxNumPrevEvents = 5

        self.current_row = None
    def get(self, colname, row):
        data = self.data_file.get("Events", colname, row)

        return data

    def get_data_calibrated(self, row):
        if self.current_row == row:
            return self.calib_data

        data = self.data_file.get("Events", "Data", row)
        sc = self.data_file.get("Events", "StartCellData", row)
        data = data.reshape(1440, -1)

        calib_data = np.zeros_like(data, np.float32)

        for i in range(1440):
            calib_data[i] = data[i] + self.off[i, sc[i]:sc[i]+len(calib_data[i])] - self.trg[i]
            calib_data[i] *= self.gain[i, sc[i]:sc[i]+len(calib_data[i])]

        calib_data = self._remove_jumps(calib_data, sc)
        self._remove_spikes_in_place(calib_data)

        self.calib_data = calib_data
        self.current_row = row
        return calib_data

    def _remove_jumps(self, calib_data, sc):
        roi = calib_data.shape[1]

        for old_sc in self.previous_start_cells:
            correct_step(
                calib_data,
                dists=(old_sc - sc + roi+10 + 1024)%1024
            )
            correct_step(
                calib_data,
                dists=(old_sc - sc + 3 + 1024)%1024
            )

        self.previous_start_cells.append(np.copy(sc))
        self.previous_start_cells = self.previous_start_cells[-self.fMaxNumPrevEvents:]
        return calib_data

    def _remove_spikes_in_place(self, calib_data):
        remove_spikes_4(calib_data)

    def __repr__(self):
        return repr(self.data_file[2]) + repr(self.drs_file[1])


def correct_step(calib_data, dists):
    roi = calib_data.shape[1]
    dists[dists >= roi] = 0
    steps = find_steps(calib_data, dists)
    patch_steps = steps.reshape(-1 , 9)[:, :8].mean(axis=1)

    if np.isnan(patch_steps).all():
        return
    average_step = np.nanmean(patch_steps)
    if average_step == 0.:
        return

    if np.nanstd(patch_steps) > 5:
        # truncated mean
        patch_steps = np.sort(patch_steps)[10:-10]
        if np.isnan(patch_steps).all():
            return
        average_step = np.nanmean(patch_steps)

    if average_step > 0:
        mask = dists[:,None] <= np.arange(calib_data.shape[1])
    else:
        mask = dists[:,None] > np.arange(calib_data.shape[1])
    calib_data[mask] -= np.abs(average_step)
    return average_step


def find_steps(data, dists):
    diff = np.diff(
        data[
            np.arange(data.shape[0])[:, None],
            dists[:, None] + [-1, 0]
        ],
        axis=1
    )
    #treat special cases
    diff[dists == 0] = np.nan
    diff[dists == data.shape[1]] = np.nan
    return diff