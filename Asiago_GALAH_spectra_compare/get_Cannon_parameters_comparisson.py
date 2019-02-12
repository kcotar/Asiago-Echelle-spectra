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
from helper_functions import spectra_normalize
imp.load_source('spectra_collection_functions', '../../Carbon-Spectra/spectra_collection_functions.py')
from spectra_collection_functions import read_pkl_spectra, CollectionParameters


def get_linelist_mask(wvl_values_in, d_wvl=0., element=None):
    idx_lines_mask = wvl_values_in < 0.

    if element is None:
        galah_linelist_use = deepcopy(galah_linelist)
    else:
        galah_linelist_use = galah_linelist[galah_linelist['Element'] == element]

    for line in galah_linelist_use:
        idx_lines_mask[np.logical_and(wvl_values_in >= line['line_start'] - d_wvl, wvl_values_in <= line['line_end'] + d_wvl)] = True

    return idx_lines_mask


galah_data_dir = '/shared/ebla/cotar/'
asiago_data_dir = galah_data_dir+'Asiago_reduced_data/'

ref_cannon = Table.read(galah_data_dir+'sobject_iraf_53_reduced_20180327.fits')['sobject_id', 'galah_id', 'ra', 'dec']
ref_cannon = join(ref_cannon, Table.read(galah_data_dir+'sobject_iraf_iDR2_180325_cannon.fits'), keys='sobject_id', join_type='left')
galah_linelist = Table.read(galah_data_dir + 'GALAH_Cannon_linelist_newer.csv')

abundances = [c for c in ref_cannon.colnames if '_abund_' in c and len(c.split('_')) == 3]
elements = [ab.split('_')[0] for ab in abundances]

R_galah = '19000'
date_string = '20180327'
spectra_file_list = ['galah_dr53_ccd1_4710_4910_wvlstep_0.040_ext4_R'+R_galah+'_'+date_string+'.pkl',
                     'galah_dr53_ccd2_5640_5880_wvlstep_0.050_ext4_R'+R_galah+'_'+date_string+'.pkl',
                     'galah_dr53_ccd3_6475_6745_wvlstep_0.060_ext4_R'+R_galah+'_'+date_string+'.pkl',
                     'galah_dr53_ccd4_7700_7895_wvlstep_0.070_ext4_R'+R_galah+'_'+date_string+'.pkl']

min_wvl = np.array([4720, 5660, 6480, 7700])
max_wvl = np.array([4900, 5870, 6730, 7880])

f_g_use = list([])
w_g_use = list([])
read_bands = [0, 1, 2]
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
idx_nan = ~ np.isfinite(f_g_use)
n_nan = np.sum(idx_nan)
if n_nan > 0:
    print 'Correcting '+str(n_nan)+' nan values'
    f_g_use[idx_nan] = 1.

# negative values handling
idx_neg = f_g_use < 0.
if np.sum(idx_neg) > 0:
    print 'Correcting neg values'
    f_g_use[idx_neg] = 0.
# large positive values handling
idx_gros = f_g_use > 1.3
if np.sum(idx_gros) > 0:
    print 'Correcting large values'
    f_g_use[idx_gros] = 1.3

survey_name = 'ORION_CLUSTERS'
obs_dir = asiago_data_dir+survey_name+'/'
obs_fits = survey_name+'_reduced.fits'
subdir_out = asiago_data_dir+survey_name+'_params_R'+R_galah

obs_fits_out = obs_fits[:-5]+'_params.fits'
survey_data = Table.read(obs_dir + obs_fits)

# add additional cols to the table if not present
for elem in elements:
    if elem not in survey_data.colnames:
        survey_data[elem] = np.nan

os.system('mkdir '+subdir_out)
os.chdir(subdir_out)

# read observed spectra
for i_o, obs_spec_data in enumerate(survey_data):
    obs_file = obs_spec_data['Asiago_id'].split('_')[0]
    print obs_file
    if 'comb' in obs_spec_data['Asiago_id']:
        obs_data = fits.open(obs_dir + obs_file + '_1D_vh_comb_norm.0001.fits')
    else:
        obs_data = fits.open(obs_dir + obs_file + '_1D_vh_norm.0001.fits')
    obs_flx_orig = obs_data[0].data
    obs_wvl_orig = obs_data[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * obs_data[0].header['CDELT1']

    # re-normalize parts of the observed spectrum that will be used in the comparison
    print 'Renormalizing observed Asiago spectrum'
    for i_b in read_bands:
        idx_renorm = np.where(np.logical_and(obs_wvl_orig >= min_wvl[i_b], obs_wvl_orig <= max_wvl[i_b]))[0]
        obs_flx_orig[idx_renorm] = spectra_normalize(obs_wvl_orig[idx_renorm] - np.mean(obs_wvl_orig[idx_renorm]),
                                                     obs_flx_orig[idx_renorm],
                                                     steps=32, sigma_low=1.5, sigma_high=3., order=11, n_min_perc=5.,
                                                     return_fit=False, func='poly')

    # RV shift
    obs_wvl_new = obs_wvl_orig * (1. - obs_spec_data['rv'] / 299792.458)
    obs_flx_use = np.interp(w_g_use, obs_wvl_new, obs_flx_orig)

    # quick spectral distance
    # idx_sim_mask = np.where(obs_flx_use < 0.975)[0]
    idx_sim_mask = np.where(obs_flx_use > 0.)[0]
    ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use)**2)[:, idx_sim_mask], axis=1))
    # ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use) ** 2), axis=1))
    idx_ref_dist = np.argsort(ref_dist)[:250]

    plt.figure(figsize=(45, 5))
    for i in idx_ref_dist:
        # plt.plot(w_g_use, f_g_use[i,:], c='blue', alpha=0.01, lw=0.8)
        plt.plot(f_g_use[i,:], c='blue', alpha=0.01, lw=0.8)
    # plt.plot(w_g_use, obs_flx_use, c='black', alpha=1, lw=0.8)
    plt.plot(obs_flx_use, c='black', alpha=1, lw=0.8)
    plt.ylim(0.3, 1.07)
    plt.xlim(-1, len(obs_flx_use))
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
        range = (np.nanmin(vals), np.nanmax(vals))
        ax[i_c].hist(vals, range=range, bins=25, label='All Cannon')
        idx_flag_ok = ref_cannon['flag_cannon'][idx_ref_dist] == 0
        if np.sum(idx_flag_ok) > 0:
            val_median = np.nanmedian(vals[idx_flag_ok])
            val_std = np.nanstd(vals[idx_flag_ok])
            val_median_all.append(val_median)
            val_std_all.append(val_std)
            print val_median, val_std
            ax[i_c].hist(vals[idx_flag_ok], range=range, bins=25, label='Flag 0 Cannon')
            ax[i_c].axvline(val_median, color='black', ls='--')
            ax[i_c].set(title='{:.2f}'.format(val_median))
        else:
            val_median_all.append(np.nan)
            val_std_all.append(np.nan)
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(obs_file+'_params_cannon_1.png', dpi=250)
    plt.close(fig)

    print val_median_all
    # search the closest within the selected parameter space
    idx_param_space = np.abs(ref_cannon['Teff_cannon'] - val_median_all[0]) < 2*val_std_all[0]
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Logg_cannon'] - val_median_all[1]) < 2*val_std_all[1])
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Fe_H_cannon'] - val_median_all[2]) < 2*val_std_all[2])
    idx_param_space = np.logical_and(idx_param_space, ref_cannon['flag_cannon'] == 0)
    idx_param_space = np.where(idx_param_space)[0]
    idx_param_space_sel = np.argsort(ref_dist[idx_param_space])[:150]
    idx_ref_dist = idx_param_space[idx_param_space_sel]

    plt.figure(figsize=(45, 5))
    for i in idx_ref_dist:
        # plt.plot(w_g_use, f_g_use[i,:], c='blue', alpha=0.01, lw=0.8)
        plt.plot(f_g_use[i,:], c='blue', alpha=0.01, lw=0.8)
    # plt.plot(w_g_use, obs_flx_use, c='black', alpha=1, lw=0.8)
    plt.plot(obs_flx_use, c='black', alpha=1, lw=0.8)
    plt.ylim(0.3, 1.07)
    plt.xlim(-1, len(obs_flx_use))
    plt.tight_layout()
    plt.savefig(obs_file + '_spectra_close_2.png', dpi=250)
    # plt.show()
    plt.close()

    if np.sum(idx_ref_dist) <= 0:
        continue

    # histogram of selected neighbours parameters
    print ' Ploting parameters histogram'
    fig, ax = plt.subplots(1, 3, figsize=(9, 5))
    val_median_all = list([])
    val_std_all = list([])
    for i_c, col in enumerate(['Teff_cannon', 'Logg_cannon', 'Fe_H_cannon']):
        vals = ref_cannon[col][idx_ref_dist]
        range = (np.nanmin(vals), np.nanmax(vals))
        ax[i_c].hist(vals, range=range, bins=25, label='All Cannon')
        idx_flag_ok = ref_cannon['flag_cannon'][idx_ref_dist] == 0
        if np.sum(idx_flag_ok) > 0:
            val_median = np.nanmedian(vals[idx_flag_ok])
            val_std = np.nanstd(vals[idx_flag_ok])
            val_median_all.append(val_median)
            val_std_all.append(val_std)
            print val_median, val_std
            ax[i_c].hist(vals[idx_flag_ok], range=range, bins=25, label='Flag 0 Cannon')
            ax[i_c].axvline(val_median, color='black', ls='--')
            ax[i_c].set(title='{:.2f}'.format(val_median))
        else:
            val_median_all.append(np.nan)
            val_std_all.append(np.nan)
    ax[1].legend()
    plt.tight_layout()
    plt.savefig(obs_file+'_params_cannon_2.png', dpi=250)
    plt.close(fig)
    
    # add results to the tale
    survey_data[i_o]['teff'] = val_median_all[0]
    survey_data[i_o]['logg'] = val_median_all[1]
    survey_data[i_o]['feh'] = val_median_all[2]

    # abundance determination in a similar ways
    idx_param_space = np.abs(ref_cannon['Teff_cannon'] - val_median_all[0]) < 2*val_std_all[0]
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Logg_cannon'] - val_median_all[1]) < 2*val_std_all[1])
    idx_param_space = np.logical_and(idx_param_space, np.abs(ref_cannon['Fe_H_cannon'] - val_median_all[2]) < 2*val_std_all[2])
    idx_param_space = np.logical_and(idx_param_space, ref_cannon['flag_cannon'] == 0)
    idx_param_space = np.where(idx_param_space)[0]

    if np.sum(idx_param_space) <= 0:
        continue

    for i_a, abund in enumerate(abundances):
        elem = abund.split('_')[0]
        print abund
        idx_sim_mask = np.where(get_linelist_mask(w_g_use, d_wvl=0., element=elem))[0]

        if len(idx_sim_mask) == 0:
            print '  No lines found for element'
            continue

        ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use) ** 2)[:, idx_sim_mask], axis=1))
        idx_param_space_sel = np.argsort(ref_dist[idx_param_space])[:75]
        idx_ref_dist = idx_param_space[idx_param_space_sel]

        fig, ax = plt.subplots(1, 2, figsize=(9, 5))
        for i in idx_ref_dist:
            ax[0].plot(f_g_use[i, idx_sim_mask], c='blue', alpha=0.05)
        ax[0].plot(obs_flx_use[idx_sim_mask], c='black', alpha=1)

        vals = ref_cannon[abund][idx_ref_dist]
        idx_flag_ok = ref_cannon['flag_'+abund][idx_ref_dist] == 0
        if np.sum(idx_flag_ok) > 0:
            abund_median = np.nanmedian(vals[idx_flag_ok])
        else:
            abund_median = np.nan
        survey_data[i_o][elem] = abund_median

        ax[1].hist(vals, range=(-2.5, 2), bins=100, label='All abund')
        if np.sum(idx_flag_ok) > 0:
            ax[1].hist(vals[idx_flag_ok], range=(-2.5, 2), bins=100, label='Flag 0 abund')
            ax[1].axvline(abund_median, color='black', ls='--')
        ax[1].grid(color='black', ls='--', alpha=0.25)

        plt.tight_layout()
        plt.savefig(obs_file + '_params_cannon_3_'+abund+'.png', dpi=250)
        plt.close()

    survey_data.write(obs_fits_out, overwrite=True)
