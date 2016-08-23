# zfits
a pure python zfits/factfits reader

# Beware !

In order to use this, you'll need to use a fork of the original fitsio:

    pip install git+https://github.com/dneise/fistio

# Example? 

    In [1]: from zfits import FactFits

    In [2]: f = FactFits("20160809_121.fits.fz", "20160809_120.drs.fits")

    In [3]: f.get("BoardTime", 10)
    Out[3]: 
    array([1472590205, 1472602002, 1472626666, 1472581883, 1472601218,
           1472610432, 1472625263, 1472607044, 1472593864, 1472585330,
            966863783,  966861390,  966853879,  966867758,  966853486,
            966871886,  966866640,  966844337,  966858949,  966867093,
           1472625978, 1472601586, 1472563220, 1472589956, 1472637092,
           1472604001, 1472605704, 1472618909, 1472588973, 1472630187,
            -36434850,  -36409086,  -36393419,  -36385907,  -36436477,
            -36401814,  -36406230,  -36422998,  -36413420,  -36383791], dtype=int32)

    In [4]: %time d = f.get_data_calibrated(0)
    CPU times: user 1.65 s, sys: 4 ms, total: 1.66 s
    Wall time: 1.66 s
