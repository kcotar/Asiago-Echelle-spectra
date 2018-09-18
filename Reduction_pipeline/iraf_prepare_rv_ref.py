from astropy.io import fits
import numpy as np
import pandas as pd
from glob import glob

# R20000 wavelengths
wvl = np.loadtxt('LAMBDA_R20.DAT')
d_wvl = 0.07
wvl_t = np.arange(3200., 7800., d_wvl)
print wvl

# read current reference solar spectrum
ref_fits_file = fits.open('solar.fits')
ref_fits_file[0].header['CRVAL1'] = wvl_t[0]
ref_fits_file[0].header['CDELT1'] = d_wvl

"""
# Read R20000 spectra and save them in a format appropriate for the RV analysis
for asc_f in glob('T*.ASC'):
	print asc_f

	flx = np.loadtxt(asc_f)
	flx_new = np.interp(wvl_t, wvl, flx)  # linear intepolation

	# save everything
	ref_fits_file[0].data = flx_new
	ref_fits_file.writeto(asc_f.split('.')[0]+'.fits', overwrite=True)
"""

# Prepare convolved GALAH Solar spectrum to be used as a reference in fxcor
sol_ref_data = np.loadtxt('solar_spectra_conv.txt')
print sol_ref_data
wvl_ref = sol_ref_data[:,0]
flx_ref = sol_ref_data[:,1]
print wvl_ref, flx_ref
ref_fits_file[0].header['CRVAL1'] = wvl_ref[0]
ref_fits_file[0].header['CDELT1'] = 0.005
ref_fits_file[0].data = flx_ref
ref_fits_file.writeto('solar_spectra_conv.fits', overwrite=True)

# Prepare observed Echelle spectra to be used as a reference in fxcor
s_dir = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/observations/rv_stand2/'
in_spectra = ['EC60966', 'EC60970', 'EC60972', 'EC60974', 'EC60976']
rv_spectra = [8.627, -21.748, 13.495, -46.883, -83.98]
# rv_spectra = [2.108, -26.821, 6.904, -53.250, -90.766]  # incorrect test values
for i_s in range(len(in_spectra)):
	fits_data = fits.open(s_dir+in_spectra[i_s]+'_1D_vh_norm.0001.fits')
	flx_o = fits_data[0].data
	wvl_o = fits_data[0].header['CRVAL1'] + np.arange(len(flx_o)) * fits_data[0].header['CDELT1']
	wvl_o *= (1 - rv_spectra[i_s] / 299792.458)
	d_wvl_new = fits_data[0].header['CDELT1']/10.
	wvl_t = np.arange(wvl_o[0], wvl_o[-1], d_wvl_new)
	flx_new = np.interp(wvl_t, wvl_o, flx_o)  # linear intepolation
	print flx_new
	print wvl_t
	# remove spikes in data
	flx_new[flx_new>2.]=2.
	flx_new[flx_new<0.]=0.
	# Save new data back
	fits_data[0].header['CDELT1'] = d_wvl_new
	fits_data[0].header['CD1_1'] = d_wvl_new
	fits_data[0].header['CRVAL1'] = wvl_t[0]
	fits_data[0].data = flx_new
	fits_data.writeto(in_spectra[i_s]+'_Echelle_RVS.fits', overwrite=True)


