from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.modeling import models, fitting
from copy import deepcopy
from joblib import Parallel, delayed
import multiprocessing
from time import time
from scipy.spatial.distance import canberra, euclidean


# TODO : automatically find those files or supply them as function inputs
ref_dir = '/home/gregor/public_html/Astro_predmeti/normalized_spectra_R20000/'
res_list = glob(ref_dir+'*/T*V000K2SNWNVR20N.ASC')
res_list = [rf.split('.')[0] for rf in res_list]
print 'Number of ref spectra:', len(res_list)

ref_wvl_orig = np.loadtxt(ref_dir+'LAMBDA_R20.DAT')
wvl_range_param = [4600., 6800.]
idx_ref_use = np.logical_and(ref_wvl_orig >= wvl_range_param[0], ref_wvl_orig <= wvl_range_param[1])

# res_list = ['solar']
obs_list = ['EC60966', 'EC60968', 'EC60970', 'EC60972', 'EC60974', 'EC60976', 'EC61147']
obs_rv = [8.7, 1.4, -19.8, 14.9, -45.1, -74.7, -4.7]

for i_o, obs_file in enumerate(obs_list):
    rv_shifts_all = []
    rv_shifts_std_all = []

    print obs_file
    obs_data = fits.open(obs_file + '_1D_vh_norm.0001.fits')
    obs_flx_orig = obs_data[0].data
    obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']

    # RV shift
    obs_wvl_new = obs_wvl_orig * (1. - obs_rv[i_o] / 299792.458)

    def compute_similarity(ref_file):
        # print ' ', ref_file
        # ref_data = fits.open(ref_file+'.fits')
        # ref_flx_orig = ref_data[0].data
        # ref_wvl_orig = ref_data[0].header['CRVAL1'] + np.arange(len(ref_flx_orig)) * ref_data[0].header['CDELT1']
        ref_flx_orig = np.loadtxt(ref_file+'.ASC')

        obs_flx_use = np.interp(ref_wvl_orig[idx_ref_use], obs_wvl_new, obs_flx_orig)
        flx_dist = euclidean(obs_flx_use, ref_flx_orig[idx_ref_use])/len(obs_flx_use)

        return flx_dist

    ts = time()
    dist_vals = Parallel(n_jobs=10)(delayed(compute_similarity)(r_l) for r_l in res_list)
    print '   ', 'Processing time {:.2} min'.format((time()-ts)/60.)

    idx_dist_vals = np.argsort(dist_vals)[:4]
    best_ref = np.array(res_list)[idx_dist_vals]
    best_ref = [br.split('/')[-1] for br in best_ref]
    print '  Best synthetics:', best_ref
