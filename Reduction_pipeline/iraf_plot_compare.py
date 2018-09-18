from astropy.io import fits
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt

# read current reference solar spectrum
solar_data = fits.open('solar.fits')
flx_solar = solar_data[0].data
wvl_solar = solar_data[0].header['CRVAL1'] + np.arange(len(flx_solar))*solar_data[0].header['CDELT1']
plt.plot(wvl_solar, flx_solar, label='Reference')

# read 1D reduced Echelle spectrum
for fits_name in ['EC60966','EC60970','EC60974','EC60972','EC60976']:
	spec_data = fits.open(fits_name+'_Echelle_RVS.fits')
	flx_spec = spec_data[0].data
	wvl_spec = spec_data[0].header['CRVAL1'] + np.arange(len(flx_spec))*spec_data[0].header['CDELT1']
	plt.plot(wvl_spec, flx_spec, label='Observed')

plt.legend()
plt.tight_layout()
plt.show()
plt.close()




