import numpy as np
from fitsio import FITS
from .factfits import FactFits
from .remove_spikes import remove_spikes_4


class FactFitsCalib:
    def __init__(self, data_path, calib_path, pixel_ids=None):
        self.data_file = FactFits(data_path)
        self.drs_file = FITS(calib_path)

        bsl = self.drs_file[1]["BaselineMean"][0]
        bsl = bsl.reshape(1440, -1)
        bsl = np.concatenate((bsl, bsl), axis=1)
        self.bsl = bsl

        gain = self.drs_file[1]["GainMean"][0]
        gain = gain.reshape(1440, -1)
        gain = np.concatenate((gain, gain), axis=1)
        self.gain = gain

        trg = self.drs_file[1]["TriggerOffsetMean"][0]
        trg = trg.reshape(1440, -1)
        self.trg = trg

        self.previous_start_cells = []
        self.fMaxNumPrevEvents = 5

        self.rows = self.data_file.rows
        if pixel_ids is None:
            self.pixel_ids = np.arange(1440)
        else:
            self.pixel_ids = np.array(pixel_ids)

    @property
    def row(self):
        return self.data_file.row

    def __iter__(self):
        return self

    def __next__(self):
        if self.row < self.rows:
            event = next(self.data_file)
            calib_data = self.get_data_calibrated(event)

            event["CalibData"] = calib_data
            return event
        else:
            raise StopIteration

    def get_data_calibrated(self, event):
        data = event["Data"]
        sc = event["StartCellData"]

        calib_data = np.empty((len(self.pixel_ids), data.shape[1]), np.float32)
        roi = calib_data.shape[1]

        for i, pix in enumerate(self.pixel_ids):
            if sc[pix] == -1:
                continue
            sl = slice(sc[pix], sc[pix] + roi)
            calib_data[i] = data[pix] * 2000.0 / 4096.0
            calib_data[i] -= self.bsl[pix, sl]
            calib_data[i] -= self.trg[pix]
            calib_data[i] /= self.gain[pix, sl]
            calib_data[i] *= 1907.35

        calib_data = self._remove_jumps(calib_data, sc)
        self._remove_spikes_in_place(calib_data)

        self.calib_data = calib_data
        return calib_data

    def _remove_jumps(self, calib_data, sc):
        roi = calib_data.shape[1]

        for old_sc in self.previous_start_cells:
            correct_step(
                calib_data,
                dists=((old_sc - sc + roi + 10 + 1024) % 1024)[self.pixel_ids],
            )
            correct_step(
                calib_data, dists=((old_sc - sc + 3 + 1024) % 1024)[self.pixel_ids]
            )

        self.previous_start_cells.append(np.copy(sc))
        self.previous_start_cells = self.previous_start_cells[-self.fMaxNumPrevEvents :]
        return calib_data

    def _remove_spikes_in_place(self, calib_data):
        remove_spikes_4(calib_data)


def correct_step(calib_data, dists):
    roi = calib_data.shape[1]
    dists[dists >= roi] = 0
    steps = find_steps(calib_data, dists)
    patch_steps = steps.reshape(-1, 9)[:, :8].mean(axis=1)

    if np.isnan(patch_steps).all():
        return
    average_step = np.nanmean(patch_steps)
    if average_step == 0.0:
        return

    if np.nanstd(patch_steps) > 5:
        # truncated mean
        patch_steps = np.sort(patch_steps)[10:-10]
        if np.isnan(patch_steps).all():
            return
        average_step = np.nanmean(patch_steps)

    if average_step > 0:
        mask = dists[:, None] <= np.arange(calib_data.shape[1])
    else:
        mask = dists[:, None] > np.arange(calib_data.shape[1])
    calib_data[mask] -= np.abs(average_step)
    return average_step


def find_steps(data, dists):
    diff = np.diff(
        data[np.arange(data.shape[0])[:, None], dists[:, None] + [-1, 0]], axis=1
    )
    # treat special cases
    diff[dists == 0] = np.nan
    diff[dists == data.shape[1]] = np.nan
    return diff
