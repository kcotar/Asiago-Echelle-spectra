import numpy as np
import matplotlib.pyplot as plt
from sklearn.externals import joblib
from glob import glob
from astropy.io import fits
from astropy.table import Table, join

import imp
imp.load_source('helper_functions', '../../Carbon-Spectra/helper_functions.py')
from helper_functions import spectra_normalize
imp.load_source('spectra_collection_functions', '../../Carbon-Spectra/spectra_collection_functions.py')
from spectra_collection_functions import read_pkl_spectra, CollectionParameters

galah_data_dir = '/home/klemen/data4_mount/'

# R20000 synthetic reference spectra
# ref_dir = '/home/gregor/public_html/Astro_predmeti/normalized_spectra_R20000/'
# res_list = glob(ref_dir+'*/T*M05*V000K2SNWNVR20N.ASC')
# w_a = np.loadtxt(ref_dir+'LAMBDA_R20.DAT')
# ref_wvl_orig = np.loadtxt(ref_dir + 'LAMBDA_R20.DAT')

ref_cannon = Table.read(galah_data_dir+'sobject_iraf_53_reduced_20180327.fits')['sobject_id', 'galah_id', 'ra', 'dec']
ref_cannon = join(ref_cannon, Table.read(galah_data_dir+'sobject_iraf_iDR2_180325_cannon.fits'), keys='sobject_id', join_type='left')

date_string = '20180327'
spectra_file_list = ['galah_dr53_ccd1_4710_4910_wvlstep_0.040_ext4_R19000_'+date_string+'.pkl',
                     'galah_dr53_ccd2_5640_5880_wvlstep_0.050_ext4_R19000_'+date_string+'.pkl',
                     'galah_dr53_ccd3_6475_6745_wvlstep_0.060_ext4_R19000_'+date_string+'.pkl',
                     'galah_dr53_ccd4_7700_7895_wvlstep_0.070_ext4_R19000_'+date_string+'.pkl']

min_wvl = np.array([4710, 5640, 6475, 7700])
max_wvl = np.array([4910, 5880, 6775, 7895])

f_g_use = list([])
w_g_use = list([])
read_bands = [1,2]
ccd_str = '_ccd'+''.join([str(c+1) for c in read_bands])
for i_b in read_bands:  # [0, 1, 2, 3]:
    # parse interpolation and averaging settings from filename
    csv_param = CollectionParameters(spectra_file_list[i_b])
    ccd = csv_param.get_ccd()
    wvl_values = csv_param.get_wvl_values()

    # determine wvls that will be read from the spectra
    idx_wvl_read = np.where(np.logical_and(wvl_values >= min_wvl[i_b], wvl_values <= max_wvl[i_b]))[0]
    wvl_values = wvl_values[idx_wvl_read]

    # read limited number of columns instead of full spectral dataset
    print 'Reading resampled/interpolated GALAH spectra from ccd', i_b+1
    w_g_use.append(wvl_values)
    f_g_use.append(read_pkl_spectra(galah_data_dir + spectra_file_list[i_b],
                                    read_cols=idx_wvl_read, read_rows=None))

w_g_use = np.hstack(w_g_use)
f_g_use = np.hstack(f_g_use)
print f_g_use.shape

# nan values handling
idx_nan = ~ np.isfinite(spectral_data)
n_nan = np.sum(idx_nan)
if n_nan > 0:
    print 'Correcting '+str(n_nan)+' nan values'
    f_g_use[idx_nan] = 1.

# negative values handling
idx_neg = spectral_data < 0.
if np.sum(idx_neg) > 0:
    f_g_use[idx_neg] = 0.
# large positive values handling
idx_gros = f_g_use > 1.3
if np.sum(idx_gros) > 0:
    f_g_use[idx_gros] = 1.3


obs_dir = '/home/klemen/Asiago-Echelle-spectra/MCMC_Cannon_model_R20K/'
obs_list = ['EC60966', 'EC60968', 'EC60970', 'EC60972', 'EC60974', 'EC60976', 'EC61147']
obs_rv = [8.7, 1.4, -19.8, 14.9, -45.1, -74.7, -4.7]

# read observed spectra
for i_o, obs_file in enumerate(obs_list):
    print obs_file
    obs_data = fits.open(obs_dir + obs_file + '_1D_vh_norm.0001.fits')
    obs_flx_orig = obs_data[0].data
    obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']

    # RV shift
    obs_wvl_new = obs_wvl_orig * (1. - obs_rv[i_o] / 299792.458)
    obs_flx_use = np.interp(w_g_use, obs_wvl_new, obs_flx_orig)

    # re-normalize observed spectrum
    obs_flx_use = spectra_normalize(w_g_use - np.mean(w_g_use), obs_flx_use,
                                steps=11, sigma_low=2., sigma_high=3., order=9, n_min_perc=5.,
                                return_fit=False, func='poly')

    # quick spectral distance
    # idx_sim_mask = np.where(obs_flx_use < 0.975)[0]
    idx_sim_mask = np.where(obs_flx_use > 0.)[0]
    ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use)**2)[:, idx_sim_mask], axis=1))
    # ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use) ** 2), axis=1))
    idx_ref_dist = np.argsort(ref_dist)[:250]

    for i in idx_ref_dist:
        plt.plot(w_g_use, f_g_use[i,:], c='blue', alpha=0.01)
    plt.plot(w_g_use, obs_flx_use, c='black', alpha=1)
    plt.tight_layout()
    plt.savefig(obs_file + '_spectra_close_1.png', dpi=250)
    # plt.show()
    plt.close()

    # histogram of selected neighbours parameters
    print ' Ploting parameters histogram'
    fig, ax = plt.subplots(1, 3, figsize=(9, 5))
    val_median_all = list([])
    val_std_all = list([])
    for i_c, col in enumerate(['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']):
        vals = ref_cannon[col][idx_ref_dist]
        idx_flag_ok = ref_cannon['flag_cannon'][idx_ref_dist] == 0
        range = (np.nanmin(vals), np.nanmax(vals))
        val_median = np.nanmedian(vals[idx_flag_ok])
        val_std = np.nanstd(vals[idx_flag_ok])
        val_median_all.append(val_median)
        val_std_all.append(val_std)
        print val_median, val_std
        ax[i_c].hist(vals, range=range, bins=25, label='All Cannon')
        ax[i_c].hist(vals[idx_flag_ok], range=range, bins=25, label='Flag 0 Cannon')
        ax[i_c].axvline(val_median, color='black', ls='--')
        ax[i_c].set(title='{:.2f}'.format(val_median))
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(obs_file+'_params_cannon_1.png', dpi=250)
    plt.close(fig)

    # search the closest within the selected parameter space
    idx_param_space = np.abs(ref_cannon['Teff_cannon'] - val_median_all[0]) < 2*val_std_all[0]
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Logg_cannon'] - val_median_all[1]) < 2*val_std_all[1])
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Fe_H_cannon'] - val_median_all[2]) < 2*val_std_all[2])
    idx_param_space = np.logical_and(idx_param_space, ref_cannon['flag_cannon'] == 0)
    idx_param_space = np.where(idx_param_space)[0]
    idx_param_space_sel = np.argsort(ref_dist[idx_param_space])[:150]
    idx_ref_dist = idx_param_space[idx_param_space_sel]

    for i in idx_ref_dist:
        plt.plot(w_g_use, f_g_use[i,:], c='blue', alpha=0.01)
    plt.plot(w_g_use, obs_flx_use, c='black', alpha=1)
    plt.tight_layout()
    plt.savefig(obs_file + '_spectra_close_2.png', dpi=250)
    # plt.show()
    plt.close()

    # histogram of selected neighbours parameters
    print ' Ploting parameters histogram'
    fig, ax = plt.subplots(1, 3, figsize=(9, 5))
    val_median_all = list([])
    val_std_all = list([])
    for i_c, col in enumerate(['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']):
        vals = ref_cannon[col][idx_ref_dist]
        idx_flag_ok = ref_cannon['flag_cannon'][idx_ref_dist] == 0
        range = (np.nanmin(vals), np.nanmax(vals))
        val_median = np.nanmedian(vals[idx_flag_ok])
        val_std = np.nanstd(vals[idx_flag_ok])
        val_median_all.append(val_median)
        val_std_all.append(val_std)
        print val_median, val_std
        ax[i_c].hist(vals, range=range, bins=25, label='All Cannon')
        ax[i_c].hist(vals[idx_flag_ok], range=range, bins=25, label='Flag 0 Cannon')
        ax[i_c].axvline(val_median, color='black', ls='--')
        ax[i_c].set(title='{:.2f}'.format(val_median))
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(obs_file+'_params_cannon_2.png', dpi=250)
    plt.close(fig)

    # abundance determination in a similar ways
