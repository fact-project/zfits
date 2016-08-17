# zfits
a pure python zfits/factfits reader

# Beware .. this has no setup.py yet 

you need to `pip install fitsio` (currently version 0.9.10) before using it.

# Example? 

    In [1]: from zfits import FactFits

    In [2]: f = FactFits("20160809_121.fits.fz", "20160809_120.drs.fits")

    In [3]: f
    Out[3]: 

      file: /tmp/20160809_121.fits.fz
      extension: 2
      type: BINARY_TBL
      extname: Events
      rows: 1000
      column info:
        EventNum            u1  varray[16]
        TriggerNum          u1  varray[16]
        TriggerType         u1  varray[14]
        NumBoards           u1  varray[16]
        UnixTimeUTC         u1  varray[20]
        BoardTime           u1  varray[172]
        StartCellData       u1  varray[2892]
        StartCellTimeMarker
                            u1  varray[332]
        Data                u1  varray[223270]
        TimeMarker          u1  varray[0]
      file: 20160809_120.drs.fits
      extension: 1
      type: BINARY_TBL
      extname: DrsCalibration
      rows: 1
      column info:
        RunNumberBaseline
                            i4  
        RunNumberGain       i4  
        RunNumberTriggerOffset
                            i4  
        BaselineMean        f4  array[1474560]
        BaselineRms         f4  array[1474560]
        GainMean            f4  array[1474560]
        GainRms             f4  array[1474560]
        TriggerOffsetMean
                            f4  array[432000]
        TriggerOffsetRms
                            f4  array[432000]
        TriggerOffsetTMMean
                            f4  
        TriggerOffsetTMRms
                            f4  

    In [5]: f.get("BoardTime", 10)
    Out[5]: 
    array([1472590205, 1472602002, 1472626666, 1472581883, 1472601218,
           1472610432, 1472625263, 1472607044, 1472593864, 1472585330,
            966863783,  966861390,  966853879,  966867758,  966853486,
            966871886,  966866640,  966844337,  966858949,  966867093,
           1472625978, 1472601586, 1472563220, 1472589956, 1472637092,
           1472604001, 1472605704, 1472618909, 1472588973, 1472630187,
            -36434850,  -36409086,  -36393419,  -36385907,  -36436477,
            -36401814,  -36406230,  -36422998,  -36413420,  -36383791], dtype=int32)

    In [6]: %time d = f.get_data_calibrated(0)
    CPU times: user 1.65 s, sys: 4 ms, total: 1.66 s
    Wall time: 1.66 s
