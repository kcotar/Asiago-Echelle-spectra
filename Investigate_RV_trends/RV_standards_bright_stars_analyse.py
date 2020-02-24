import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as un
import astropy.constants as const
from astropy.io import fits
from astropy.table import Table, join
from glob import glob
from os import chdir, system, path

asiago_data_dir = '/shared/ebla/cotar/Asiago_reduced_data/'

# # ---- BRIGHT RV STANDARDS ----
subfolder_name = 'RVS_COMP'
out_dir = 'RV_standards_bright'
# # -----------------------------

# # ---- GAIA RV STANDARDS ------
# subfolder_name = 'GAIA_RV_STAND'
# out_dir = 'RV_standards_Gaia'
# # -----------------------------

# ---- ACQARIOUS STREAM CANDIDATES
# subfolder_name = 'GREGOR_TEST'
# out_dir = 'AQU_stream_candidates'
# -----------------------------

in_dir = asiago_data_dir + subfolder_name + '/'
# read data about reduced spectra
obs_data = Table.read(in_dir + subfolder_name + '_reduced.fits')

# order by MJD - from old to new
obs_data = obs_data[np.argsort(obs_data['MJD'])]

# go to a output subfolder
system('mkdir ' + out_dir)
chdir(out_dir)

c_value = const.c.to(un.km/un.s).value

# first plot RV values from different exposures
for star_name in np.unique(obs_data['obj_name']):
    star_data = obs_data[obs_data['obj_name'] == star_name]
    star_name_use = star_name.replace(' ', '_')
    print('Plotting RVs of', star_name)

    rvs_all = list()
    for ic, spectrum in enumerate(star_data):
        rvs_txt = in_dir + spectrum['Asiago_id'] + '_solar_rv.txt'
        rvs_data = np.genfromtxt(rvs_txt, delimiter=',')
        spectrum_fits = fits.open(in_dir + spectrum['Asiago_id'] + '.fits')
        obs_date = spectrum_fits[0].header['DATE-OBS'].split('T')[0]
        spectrum_fits.close()

        # print info about spectrum
        s_f_h = spectrum_fits[0].header
        if 'TEMP-EC' in s_f_h:
            print(spectrum['Asiago_id'].split('_')[0], '  ', s_f_h['DATE-OBS'][:-5], '  ', s_f_h['ALTITUDE'], '  ', s_f_h['AZIMUTH'], '  ', s_f_h['TEMP-EC'], '  ', s_f_h['HUM-INT'])
        else:
            print(spectrum['Asiago_id'].split('_')[0], '  ', s_f_h['DATE-OBS'][:-5], '  ', s_f_h['ALTITUDE'], '  ', s_f_h['AZIMUTH'])

        rvs_all.append(rvs_data[:, 1])
        plt.scatter(rvs_data[:, 0], rvs_data[:, 1], c='C'+str(ic), alpha=0.6, s=20, label=obs_date)
        plt.axhline(spectrum['rv'], ls='--', c='C'+str(ic), alpha=0.3, label='')
    plt.title(star_name)
    plt.xlabel('Echelle extracted order number')
    plt.ylabel('Radial velocity per order [km/s]')
    rvs_median = np.nanmedian(rvs_all)
    plt.ylim(rvs_median - 6., rvs_median + 6.)
    plt.xlim(0, 31)
    plt.grid(ls='--', c='black', alpha=0.3)
    plt.axvspan(-1, 5, alpha=0.3, color='black')
    plt.axvspan(25, 33, alpha=0.3, color='black')
    plt.tight_layout()
    plt.legend(loc=2)
    plt.savefig(star_name_use + '_rvs_202001.png', dpi=300)
    plt.close()

    print(' ')

raise SystemExit

# plot sections of individual spectra and arcs that were used in the analysis
# also add reference spectrum the plot
ref_fits = fits.open(in_dir + 'solar.fits')
ref_flx_orig = ref_fits[0].data
ref_wvl_orig = ref_fits[0].header['CRVAL1'] + np.arange(len(ref_flx_orig)) * ref_fits[0].header['CDELT1']
ref_fits.close()

for star_name in np.unique(obs_data['obj_name']):
    star_data = obs_data[obs_data['obj_name'] == star_name]
    star_name_use = star_name.replace(' ', '_')
    print('Plotting spectra of', star_name)

    # plot individual sections of spectra
    wvl_d = 30.
    wvl_s = 3700.
    wvl_e = 7300.

    for wvl_s in np.arange(wvl_s, wvl_e, wvl_d):

        fig1, ax1 = plt.subplots(1, 1, figsize=(16, 5))
        fig2, ax2 = plt.subplots(1, 1, figsize=(16, 5))

        for ic, spectrum in enumerate(star_data):
            spectrum_fits = fits.open(in_dir + spectrum['Asiago_id'].split('_')[0] + '_1D_vh_norm.0001.fits')
            arc_fits = fits.open(in_dir + spectrum['Asiago_arc_id'].split('_')[0] + '_1D_vh_norm.0001.fits')

            # extract data from the files
            obs_date = spectrum_fits[0].header['DATE-OBS'].split('T')[0]

            obs_flx_orig = spectrum_fits[0].data
            obs_wvl_orig = spectrum_fits[0].header['CRVAL1'] + np.arange(len(obs_flx_orig)) * spectrum_fits[0].header['CDELT1']
            obs_wvl_rv = obs_wvl_orig * (1. - spectrum['rv'] / c_value)

            # additional arc normalization/scaling
            flx_norm_val = np.percentile(obs_flx_orig[np.logical_and(obs_wvl_orig >= wvl_s, obs_wvl_orig <= wvl_s + wvl_d)], [1., 99.])
            obs_flx_orig /= flx_norm_val[1]

            arc_flx_orig = arc_fits[0].data
            arc_wvl_orig = arc_fits[0].header['CRVAL1'] + np.arange(len(arc_flx_orig)) * arc_fits[0].header['CDELT1']
            # additional arc normalization/scaling
            arc_norm_val = np.percentile(arc_flx_orig[np.logical_and(arc_wvl_orig >= wvl_s, arc_wvl_orig <= wvl_s + wvl_d)], [0.5, 99.5])
            arc_flx_orig = (arc_flx_orig-arc_norm_val[0])/arc_norm_val[1]

            spectrum_fits.close()
            arc_fits.close()

            ax1.plot(obs_wvl_rv, obs_flx_orig, c='C'+str(ic), alpha=0.7, label=obs_date, lw=0.6)
            ax2.plot(arc_wvl_orig, arc_flx_orig, c='C'+str(ic), alpha=0.7, label=obs_date, lw=0.6)

        # add RV reference spectrum to the spectra plot
        ax1.plot(ref_wvl_orig, ref_flx_orig, c='black', alpha=0.7, label='RV ref', lw=0.6)

        ax1.set(xlim=(wvl_s, wvl_s + wvl_d), ylim=(0.0, 1.1), title=star_name,
                ylabel='Normalized spectrum', xlabel='Determined wvl corrected for VHELIO and median RV')
        ax2.set(xlim=(wvl_s, wvl_s + wvl_d), title=star_name,
                ylim=(0.0, 1.0),
                ylabel='Scaled and normalized arc', xlabel='Determined wvl')

        ax1.grid(ls='--', c='black', alpha=0.2)
        ax2.grid(ls='--', c='black', alpha=0.2)
        ax1.legend(loc=3)
        ax2.legend(loc=2)
        fig1.tight_layout()
        fig2.tight_layout()

        fig1.savefig(star_name_use + '_' + str(int(wvl_s)) + '_spec.png', dpi=250)
        fig2.savefig(star_name_use + '_' + str(int(wvl_s)) + '_arc.png', dpi=250)
        plt.close('all')
