import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join

import imp
imp.load_source('helper_functions', '../../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize
imp.load_source('spectra_collection_functions', '../../Carbon-Spectra/spectra_collection_functions.py')
from spectra_collection_functions import read_pkl_spectra, CollectionParameters

data_dir = '/shared/ebla/cotar/'
obs_dir = data_dir + 'Asiago_reduced_data/GREGOR_TEST_2/'

sola_data = np.genfromtxt(data_dir + 'Solar_data_dr53/solar_spectra_conv.txt', delimiter=' ')
sola_flx = sola_data[:, 1]
sola_wvl = sola_data[:, 0]

tell_data = np.genfromtxt(data_dir + 'telluric_spectra_conv.dat', delimiter=' ')
tell_flx = tell_data[:, 1]
tell_wvl = tell_data[:, 0]

obs_file = 'EC59444'
obs_data = fits.open(obs_dir + obs_file + '_1D_vh_norm.0001.fits')

obs_flx_orig = obs_data[0].data
obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']
# remove VHELIO correction shift
# obs_wvl_rv = obs_wvl_orig * (1. - obs_data[0].header['VHELIO'] / 299792.458)
# add actual determined RV shift
rv_star = -10.5
obs_wvl_rv = obs_wvl_orig * (1. - rv_star / 299792.458)

plt.plot(obs_wvl_rv, obs_flx_orig, lw=1, c='black')
plt.plot(sola_wvl, sola_flx, lw=1, c='C0')
plt.plot(tell_wvl, 1. - (1. - tell_flx)*0.5, lw=1, c='C2')
plt.ylim(0, 1.2)
plt.show()
plt.close()

# resample to the same wavelength scale
sola_flx_res = np.interp(obs_wvl_rv, sola_wvl, sola_flx)

plt.plot(obs_wvl_rv, sola_flx_res-obs_flx_orig)
plt.ylim(-0.5, 0.5)
plt.show()
plt.close()