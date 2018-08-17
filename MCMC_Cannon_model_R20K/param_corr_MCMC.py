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
import emcee, corner

# load Cannon model
import thecannon as tc

cannon_model = tc.CannonModel.read('model_cannon_R20000_4000K_7000K.dat')
thetas = cannon_model.theta
vectorizer = cannon_model.vectorizer
fid = cannon_model._fiducials
sca = cannon_model._scales
ref_wvl = cannon_model.dispersion
wvl_range_param = [4600., 6800.]
idx_ref_use = np.logical_and(ref_wvl >= wvl_range_param[0], ref_wvl <= wvl_range_param[1])
ref_wvl_use = ref_wvl[idx_ref_use]

def get_cannon(teff, logg, feh):
    sint = thetas[:, 0] * 0.0
    labs = (np.array([teff, logg, feh]) - fid) / sca
    vec = vectorizer(labs)
    for i, j in enumerate(vec):
        sint += thetas[:, i] * j
    return sint


def eval_params(param):
    t, l, f = param
    if not 3500. < t < 7000.:
        return False
    if not 0. < l < 5.5:
        return False
    if not -3. < f < 1.:
        return False
    return True


def lnprob_flx_fit(param, obs_flx, idx_flx):
    if eval_params(param):
        teff, logg, feh = param
        # get cannon model
        ref_flx = get_cannon(teff, logg, feh)[idx_flx]
        lnprob_flux = -0.5 * (np.sum((ref_flx - obs_flx) ** 2))  # / flx_s ** 2 + np.log(2 * np.pi * flx_s ** 2)))
        return lnprob_flux
    else:
        return -np.inf


def run_mcmc_params(p0, obs_flx, idx_flx, n_s):
    t_1 = time()
    sampler = emcee.EnsembleSampler(len(p0), len(p0[0]), lnprob_flx_fit,
                                    threads=10,
                                    args=(obs_flx, idx_flx))
    p0, lnp, _ = sampler.run_mcmc(p0, n_s)
    sampler.pool.close()
    sampler.pool = None
    print '   - took {:.2f} min'.format((time() - t_1) / 60.)
    return sampler


def _p0_generate(guess_low, guess_hig, nwalkers):
    # uniform distribution of initial values
    ndim = len(guess_low)
    p0 = list([])
    for i_w in range(nwalkers):
        p0_new = guess_low + np.random.rand(ndim) * (guess_hig - guess_low)
        p0.append(p0_new)
    return p0


ref_dir = '/home/gregor/public_html/Astro_predmeti/normalized_spectra_R20000/'
res_list = glob(ref_dir+'*/T*V000K2SNWNVR20N.ASC')
ref_wvl_orig = np.loadtxt(ref_dir+'LAMBDA_R20.DAT')
# read observed spectra or synthetic spectra (only for a test)
for i_o, ref_file in enumerate(res_list):
    ref_name = ref_file.split('/')[-1].split('.')[0]
    print ref_name
    ref_flx_orig = np.loadtxt(ref_file)
    ref_flx_use = np.interp(ref_wvl_use, ref_wvl_orig, ref_flx_orig)

    # perform MCMC fitting
    p0_0 = _p0_generate(np.array([3500., 0., -2.]), np.array([7500., 5., 0.5]), 100)
    sampler = run_mcmc_params(p0_0, ref_flx_use, idx_ref_use, 1200)

    kernel_fit = np.nanmedian(sampler.flatchain, axis=0)
    print '  Params are:', kernel_fit
    c_fig = corner.corner(sampler.flatchain,
                          truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                          labels=['teff', 'logg', 'feh'], bins=75)
    c_fig.savefig(ref_name + '_corner.png', dpi=200)
    plt.close(c_fig)


obs_list = ['EC60966', 'EC60968', 'EC60970', 'EC60972', 'EC60974', 'EC60976', 'EC61147']
obs_rv = [8.7, 1.4, -19.8, 14.9, -45.1, -74.7, -4.7]
# read observed spectra
for i_o, obs_file in enumerate(obs_list):
    print obs_file
    obs_data = fits.open(obs_file + '_1D_vh_norm.0001.fits')
    obs_flx_orig = obs_data[0].data
    obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']

    # RV shift
    obs_wvl_new = obs_wvl_orig * (1. - obs_rv[i_o] / 299792.458)
    obs_flx_use = np.interp(ref_wvl_use, obs_wvl_new, obs_flx_orig)

    # perform MCMC fitting
    p0_0 = _p0_generate(np.array([3500., 0., -2.]), np.array([7500., 5., 0.5]), 100)
    sampler = run_mcmc_params(p0_0, obs_flx_use, idx_ref_use, 1200)

    # sampler.flatchain
    # sampler.flatlnprobability
    kernel_fit = np.nanmedian(sampler.flatchain, axis=0)
    print '  Params are:', kernel_fit
    c_fig = corner.corner(sampler.flatchain,
                          truths=kernel_fit, quantiles=[0.16, 0.5, 0.84],
                          labels=['teff', 'logg', 'feh'], bins=75)
    c_fig.savefig(obs_file + '_corner.png', dpi=200)
    plt.close(c_fig)


