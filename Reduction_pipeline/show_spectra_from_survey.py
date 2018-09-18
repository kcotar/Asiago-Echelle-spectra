from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

survey = 'NGC6940'
survey_dir = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/observations/'+survey+'/'
survey_fits = ''+survey+'_reduced.fits'

survey_data = Table.read(survey_dir + survey_fits)
print survey_data.colnames

for obj in survey_data:
	print obj['Asiago_id']
	f_name = obj['Asiago_id'].split('_')[0]
	fits_data = fits.open(survey_dir + f_name + '_1D_vh_comb_norm.0001.fits')
	flx_o = fits_data[0].data
	wvl_o = fits_data[0].header['CRVAL1'] + np.arange(len(flx_o)) * fits_data[0].header['CDELT1']
	wvl_o *= (1 - obj['rv'] / 299792.458)
	plt.plot(wvl_o, flx_o, label=f_name)
plt.ylim(0, 1.5)
plt.xlim(4100, 6900)
plt.legend()
plt.tight_layout()
plt.show()
plt.close()

