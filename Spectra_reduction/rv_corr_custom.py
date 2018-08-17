from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.signal import correlate
from lmfit.models import GaussianModel, VoigtModel
from astropy.modeling import models, fitting
from copy import deepcopy
from joblib import Parallel, delayed
import multiprocessing
from time import time


def spectra_logspace(flx, wvl):
    """

    :param flx:
    :param wvl:
    :return:
    """
    wvl_new = np.logspace(np.log10(wvl[0]), np.log10(wvl[-1]), num=len(wvl))
    return np.interp(wvl_new, wvl, flx), wvl_new


def correlate_spectra(obs_flx, obs_wvl, ref_flx, ref_wvl):

    # convert spectra sampling to logspace
    obs_flux_res_log, _ = spectra_logspace(obs_flx, obs_wvl)
    ref_flux_sub_log, wvl_log = spectra_logspace(ref_flx, ref_wvl)
    wvl_step = ref_wvl[1] - ref_wvl[0]

    # correlate the two spectra
    min_flux = 0.95
    ref_flux_sub_log[ref_flux_sub_log > min_flux] = 0.
    obs_flux_res_log[obs_flux_res_log > min_flux] = 0.
    corr_res = correlate(ref_flux_sub_log, obs_flux_res_log, mode='same', method='fft')

    # plt.plot(corr_res)
    # plt.show()
    # plt.close()

    # create a correlation subset that will actually be analysed
    corr_w_size = 100
    corr_c_off = np.int64(len(corr_res) / 2.)
    corr_pos_min = corr_c_off - corr_w_size
    corr_pos_max = corr_c_off + corr_w_size
    # print corr_pos_min, corr_pos_max
    corr_res_sub = corr_res[corr_pos_min:corr_pos_max]
    corr_res_sub -= np.median(corr_res_sub)
    corr_res_sub_x = np.arange(len(corr_res_sub))

    # analyze correlation function by fitting gaussian/voigt/lorentzian distribution to it
    fit_model = VoigtModel()
    parameters = fit_model.guess(corr_res_sub, x=corr_res_sub_x)
    corr_fit_res = fit_model.fit(corr_res_sub, parameters, x=corr_res_sub_x)
    corr_center = corr_fit_res.params['center'].value

    # plt.plot(corr_res_sub)
    # plt.axvline(corr_center)
    # plt.show()
    # plt.close()

    # determine the actual shift
    idx_no_shift = np.int32(len(corr_res) / 2.)
    idx_center = corr_c_off - corr_w_size + corr_center
    log_shift_px = idx_no_shift - idx_center
    log_shift_wvl = log_shift_px * wvl_step

    wvl_log_new = wvl_log - log_shift_wvl
    rv_shifts = (wvl_log_new[1:] - wvl_log_new[:-1]) / wvl_log_new[:-1] * 299792.458 * log_shift_px

    if log_shift_wvl < 2:
        return np.nanmedian(rv_shifts)
    else:
        # something went wrong
        return np.nan

# for file in ['EC60968','EC61147']:
#     ref_data = fits.open(file + '_1D_vh_norm.0001.fits')
#     ref_flx_orig = ref_data[0].data
#     ref_wvl_orig = ref_data[0].header['CRVAL1'] + np.arange(len(ref_flx_orig)) * ref_data[0].header['CDELT1']
#     plt.plot(ref_wvl_orig, ref_flx_orig)
# plt.ylim(0, 1.3)
# plt.show()
# plt.close()

# TODO : automatically find those files or supply them as function inputs
res_list = glob('T*K2SNWNVR20N.fits')
res_list = [rf.split('.')[0] for rf in res_list]
# res_list = ['solar']
obs_list = ['EC60966','EC60968','EC60970','EC60972','EC60974','EC60976','EC61147']

n_subsections = 31

for obs_file in obs_list:
    rv_shifts_all = []
    rv_shifts_std_all = []

    print obs_file

    for ref_file in res_list:
        print ' ', ref_file
        ref_data = fits.open(ref_file+'.fits')
        ref_flx_orig = ref_data[0].data
        ref_wvl_orig = ref_data[0].header['CRVAL1'] + np.arange(len(ref_flx_orig)) * ref_data[0].header['CDELT1']

        ref_wvl_d = ref_data[0].header['CDELT1']/5.
        ref_wvl = np.arange(ref_wvl_orig[0], ref_wvl_orig[-1], ref_wvl_d)
        ref_flx = np.interp(ref_wvl, ref_wvl_orig, ref_flx_orig)

        def correlate_order(i_s):
            # load the correct Echelle order
            try:
                obs_data = np.loadtxt(obs_file + '_vh_norm_order{:02.0f}.txt'.format(i_s))
                # print obs_data
            except:
                # print '  ', 'File for order {:.0f} not found'.format(i_s)
                return np.nan
            obs_flx = obs_data[:, 1]
            obs_wvl = obs_data[:, 0]

            # print i_s
            wvl_len = len(obs_wvl)
            wvl_beg = obs_wvl[int(wvl_len / 4.)]
            wvl_end = obs_wvl[int(wvl_len * 3. / 4.)]

            # resample data to the same wavelength step
            idx_ref_use = np.logical_and(ref_wvl >= wvl_beg, ref_wvl <= wvl_end)
            if np.sum(idx_ref_use) < 20:
                return np.nan

            ref_wvl_use = ref_wvl[idx_ref_use]
            ref_flx_use = ref_flx[idx_ref_use]

            obs_flx_use = np.interp(ref_wvl_use, obs_wvl, obs_flx)

            # add re-normalization step if needed

            # perform correlation between the datasets
            return correlate_spectra(obs_flx_use, ref_wvl_use, ref_flx_use, ref_wvl_use)

        ts = time()
        rv_shifts = Parallel(n_jobs=20)(delayed(correlate_order)(i) for i in range(1, n_subsections+1))
        print '   ', 'Processing time {:.2} s'.format(time()-ts)
        # for i_s in range(1, n_subsections+1):
        #     # load the correct Echelle order
        #     try:
        #         obs_data = np.loadtxt(obs_file+'_vh_norm_order{:02.0f}.txt'.format(i_s))
        #         # print obs_data
        #     except:
        #         print '  ', 'File for order {:.0f} not found'.format(i_s)
        #         rv_shifts.append(np.nan)
        #         continue
        #     obs_flx = obs_data[:, 1]
        #     obs_wvl = obs_data[:, 0]
        #
        #     # print i_s
        #     wvl_len = len(obs_wvl)
        #     wvl_beg = obs_wvl[int(wvl_len/4.)]
        #     wvl_end = obs_wvl[int(wvl_len*3./4.)]
        #
        #     # resample data to the same wavelength step
        #     idx_ref_use = np.logical_and(ref_wvl >= wvl_beg, ref_wvl <= wvl_end)
        #     if np.sum(idx_ref_use) < 20:
        #         rv_shifts.append([np.nan])
        #         continue
        #
        #     ref_wvl_use = ref_wvl[idx_ref_use]
        #     ref_flx_use = ref_flx[idx_ref_use]
        #
        #     obs_flx_use = np.interp(ref_wvl_use, obs_wvl, obs_flx)
        #
        #     # add re-normalization step if needed
        #
        #     # perform correlation between the datasets
        #     rv_shifts.append(correlate_spectra(obs_flx_use, ref_wvl_use, ref_flx_use, ref_wvl_use))

        rv_shifts = np.array(rv_shifts)
        rv_median = np.nanmedian(rv_shifts)
        rv_shifts[np.abs(rv_shifts - rv_median) > 4] = np.nan
        if np.sum(np.isfinite(rv_shifts)) < n_subsections/2:
            rv_std = np.nan
        else:
            rv_std = np.nanstd(rv_shifts)

        rv_shifts_all.append(rv_shifts)
        rv_shifts_std_all.append(rv_std)

    idx_best = np.argsort(rv_shifts_std_all)[0]
    rv_shifts = rv_shifts_all[idx_best]
    rv_median = np.nanmedian(rv_shifts)
    rv_x = np.arange(len(rv_shifts))
    plt.scatter(rv_x, rv_shifts)
    plt.axhline(rv_median)
    plt.ylim(rv_median-4, rv_median+4)
    plt.xlim(rv_x[0]-1, rv_x[-1]+1)
    plt.xlabel('Echelle order')
    plt.ylabel('Radial velocity')
    plt.title('Spectrum: '+obs_file+'     Reference: '+res_list[idx_best])
    plt.tight_layout()
    # plt.show()
    plt.savefig(obs_file+'.png', dpi=300)
    plt.close()

