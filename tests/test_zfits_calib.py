import fitsio
import numpy as np


def test_with_fact_tools():
    from zfits import FactFits

    fact_tools_reference = fitsio.read(
        'tests/resources/20160817_016_calibrated.fits'
    )['DataCalibrated'][0].reshape((1440, -1))

    print(fact_tools_reference.shape)

    f = FactFits(
        'tests/resources/20160817_016.fits.fz',
        'tests/resources/20160817_030.drs.fits.gz',
    )

    assert np.all(np.isclose(
        f.get_data_calibrated(0)[0, 10:250], fact_tools_reference[0, 10:250]
    ))
