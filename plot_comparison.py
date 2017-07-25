# coding: utf-8
from zfits import FactFits
from astropy.io import fits
import matplotlib.pyplot as plt

f = FactFits('zfits/test_data/20160817_016.fits.fz', 'zfits/test_data/20160817_030.drs.fits.gz')
facttools = fits.open('zfits/test_data/20160817_016_calibrated.fits')

e = f.get_data_calibrated(0)
ft_e = facttools[1].data['DataCalibrated'][0].reshape((1440, -1))

plt.subplot(2, 1, 1)
plt.plot(ft_e[0], label='FACT-Tools')
plt.plot(e[0], label='zfits')
plt.legend()

plt.subplot(2, 1, 2)
plt.axhline(2000 / 4096, color='C1')
plt.plot(e[0] / ft_e[0])
plt.ylabel('zfits / facttools')
plt.ylim(0.8, 1.2)
plt.savefig('comp2.png', dpi=300)
