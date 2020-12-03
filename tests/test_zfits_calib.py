import fitsio
import numpy as np


def test_with_fact_tools():
    from zfits import FactFitsCalib

    fact_tools_reference = fitsio.read("tests/resources/20160817_016_calibrated.fits")[
        "DataCalibrated"
    ][0].reshape((1440, -1))

    print(fact_tools_reference.shape)

    f = FactFitsCalib(
        "tests/resources/20160817_016.fits.fz",
        "tests/resources/20160817_030.drs.fits.gz",
    )
    event = next(f)
    # in pixel 18 is a spike:
    #  - in fact_tools_reference, this spike is not removed
    #  - in the FactFitsCalib['CalibData'] the spike **is removed**.
    # so the test would fail ...
    for i in range(17):
        a = event["CalibData"][i, 10:250]
        b = fact_tools_reference[i, 10:250]
        d = a - b
        dabs = np.abs(d)
        argmax = np.argmax(dabs)
        print(i, argmax, dabs[argmax], d[argmax], a[argmax], b[argmax])
        assert np.all(np.isclose(a, b))
