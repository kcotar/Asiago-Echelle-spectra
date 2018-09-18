from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.signal import correlate
from lmfit.models import GaussianModel, VoigtModel
from astropy.modeling import models, fitting
from copy import deepcopy
from time import time



def get_julia_dates_header(filename):
    file_data = fits.open(filename+'.fits')
    jd = file_data[0].header['JD']
    mjd = file_data[0].header['MJD']
    # check if those dates are correct in the case when spectra are merged together
    file_data.close()
    return mjd, jd



def spectra_logspace(flx, wvl):
    """

    :param flx:
    :param wvl:
    :return:
    """
    wvl_new = np.logspace(np.log10(wvl[0]), np.log10(wvl[-1]), num=len(wvl))
    return np.interp(wvl_new, wvl, flx), wvl_new


def correlate_spectra(obs_flx, obs_wvl, ref_flx, ref_wvl, plot=None):

    # convert spectra sampling to logspace
    obs_flux_res_log, _ = spectra_logspace(obs_flx, obs_wvl)
    ref_flux_sub_log, wvl_log = spectra_logspace(ref_flx, ref_wvl)
    wvl_step = ref_wvl[1] - ref_wvl[0]

    # correlate the two spectra
    min_flux = 0.95
    ref_flux_sub_log[ref_flux_sub_log > min_flux] = 1.
    obs_flux_res_log[obs_flux_res_log > min_flux] = 1.
    corr_res = correlate(1.-ref_flux_sub_log, 1.-obs_flux_res_log, mode='same', method='fft')

    # create a correlation subset that will actually be analysed
    corr_w_size = 100
    corr_c_off = np.int64(len(corr_res) / 2.)
    corr_pos_min = corr_c_off - corr_w_size
    corr_pos_max = corr_c_off + corr_w_size

    if plot is not None:
      plt.plot(corr_res, lw=1)
      plt.axvline(corr_pos_min, color='black')
      plt.axvline(corr_pos_max, color='black')
      plt.savefig(plot+'_1.png', dpi=300)
      plt.close()

    # print corr_pos_min, corr_pos_max
    corr_res_sub = corr_res[corr_pos_min:corr_pos_max]
    corr_res_sub -= np.median(corr_res_sub)
    corr_res_sub_x = np.arange(len(corr_res_sub))

    # analyze correlation function by fitting gaussian/voigt/lorentzian distribution to it
    fit_model = VoigtModel()
    parameters = fit_model.guess(corr_res_sub, x=corr_res_sub_x)
    corr_fit_res = fit_model.fit(corr_res_sub, parameters, x=corr_res_sub_x)
    corr_center = corr_fit_res.params['center'].value

    if plot is not None:
      plt.plot(corr_res_sub, lw=1)
      plt.axvline(corr_center, color='black')
      plt.savefig(plot+'_2.png', dpi=200)
      plt.close()

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


def correlate_order(order_txt_file, ref_flx, ref_wvl, plot=False):
# load the correct Echelle order
  try:
    obs_data = np.loadtxt(order_txt_file)
    # print obs_data
  except:
    # print '  ', 'File not found',order_txt_file
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
  if plot:
    return correlate_spectra(obs_flx_use, ref_wvl_use, ref_flx_use, ref_wvl_use, plot=order_txt_file[:-4])
  else:
    return correlate_spectra(obs_flx_use, ref_wvl_use, ref_flx_use, ref_wvl_use)

def get_RV_ref_spectrum(spectrum_full_path):
  ref_data = fits.open(spectrum_full_path)
  ref_flx_orig = ref_data[0].data
  ref_wvl_orig = ref_data[0].header['CRVAL1'] + np.arange(len(ref_flx_orig)) * ref_data[0].header['CDELT1']
  # resample data to a finer resolution
  ref_wvl_d = ref_data[0].header['CDELT1']/5.
  ref_wvl = np.arange(ref_wvl_orig[0], ref_wvl_orig[-1], ref_wvl_d)
  ref_flx = np.interp(ref_wvl, ref_wvl_orig, ref_flx_orig)
  # output processed spectrum
  return ref_flx, ref_wvl


# res_list = glob('T*K2SNWNVR20N.fits')
# res_list = [rf.split('.')[0] for rf in res_list]
# res_list = ['solar']
# obs_list = ['EC60966','EC60968','EC60970','EC60972','EC60974','EC60976','EC61147']

