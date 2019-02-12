import numpy as np
import matplotlib.pyplot as plt
from sklearn.externals import joblib
from glob import glob
from astropy.io import fits
from astropy.table import Table, join
from copy import deepcopy
import os

import imp
imp.load_source('helper_functions', '../../Carbon-Spectra/helper_functions.py')
from helper_functions import get_spectra_dr52, spectra_normalize

import varconvolve as varcon
from scipy import mgrid
def kernel(s):
    """
    Constructs a normalized discrete 1D gaussian kernel
    """
    size_grid = int(s*4)
    x = mgrid[-size_grid:size_grid+1]
    g = np.exp(-(x**2/float(s**2)/2.))
    return g / np.sum(g)

min_wvl = np.array([4720, 5660, 6480, 7700])
max_wvl = np.array([4900, 5870, 6730, 7880])

data_dir = '/data4/cotar/'
obs_dir = data_dir + 'Asiago_reduced_data/GAIA_BENCH/'
ben_dir = data_dir + 'Gaia_benchmark/'
d53_dir = data_dir + 'dr5.3/'

obs_data = Table.read(obs_dir + 'GAIA_BENCH_reduced.fits')
galah_data = Table.read(ben_dir + 'Gaia_benchmark_stars_GALAH.csv')

for obs_object in obs_data:
    # load spectrum
    obs_file = obs_object['Asiago_id'].split('_')[0]
    print obs_object['obj_name'], obs_file

    try:
        if 'comb' in obs_object['Asiago_id']:
            obs_data = fits.open(obs_dir + obs_file + '_1D_vh_comb_norm.0001.fits')
        else:
            obs_data = fits.open(obs_dir + obs_file + '_1D_vh_norm.0001.fits')
    except:
        print ' Asiago spectrum not fund'
        continue
    obs_flx_orig = obs_data[0].data
    obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']
    # RV shift
    obs_wvl_rv = obs_wvl_orig * (1. - obs_object['rv'] / 299792.458)

    # plt.plot(obs_wvl_rv, obs_flx_orig, lw=2, c='black', alpha=0.75)

    for i_b in [0, 1, 2]:
        idx_renorm = np.where(np.logical_and(obs_wvl_rv >= min_wvl[i_b], obs_wvl_rv <= max_wvl[i_b]))[0]
        obs_flx_orig[idx_renorm] = spectra_normalize(obs_wvl_rv[idx_renorm] - np.mean(obs_wvl_rv[idx_renorm]),
                                                     obs_flx_orig[idx_renorm],
                                                     steps=32, sigma_low=1.5, sigma_high=3., order=11, n_min_perc=5.,
                                                     return_fit=False, func='poly')

    # use only subset
    idx_sub_asiago = np.logical_and(obs_wvl_rv >= 4500, obs_wvl_rv <= 7000)
    obs_flx_orig = obs_flx_orig[idx_sub_asiago]
    obs_wvl_rv = obs_wvl_rv[idx_sub_asiago]

    # find GALAH equivalents
    galah_get = galah_data[galah_data['std_name'] == obs_object['obj_name']]

    if len(galah_get) == 0:
        print ' No GALAH matched objects'
        continue

    for galah_object in galah_get:
        object_id = str(galah_object['sobject_id'])
        try:
            flx, wvl = get_spectra_dr52(object_id, root=d53_dir, bands=[1,2,3], extension=4)
        except:
           print ' Problem reading GALAH spectrum '+object_id
           continue

        for i_b in range(len(flx)):
            R_source = 28000.
            R_target = 19000.
            s_orig = wvl[i_b] / (2.3548 * R_source)
            s_targ = wvl[i_b] / (2.3548 * R_target)
            kernel_widths = np.sqrt(s_targ ** 2 - s_orig ** 2)
            flx[i_b] = varcon.varconvolve(wvl[i_b], flx[i_b], kernel, kernel_widths)
            # plot convolved spectrum
            plt.plot(wvl[i_b], flx[i_b], lw=2, c='red', alpha=0.5)

    plt.plot(obs_wvl_rv, obs_flx_orig, lw=2, c='blue', alpha=0.75)
    plt.ylim(0.2, 1.2)
    plt.title(obs_object['obj_name'])
    plt.tight_layout()
    plt.show()
    plt.close()