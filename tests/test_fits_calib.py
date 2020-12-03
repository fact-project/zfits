def test_calib():
    from zfits import FactFitsCalib

    f = FactFitsCalib(
        "tests/resources/testDataFile.fits.gz",
        "tests/resources/testDrsFile.drs.fits.gz",
    )
    next(f)
