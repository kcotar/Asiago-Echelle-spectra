from astropy.table import Table
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.signal import correlate
from lmfit.models import GaussianModel, VoigtModel, ConstantModel
from astropy.modeling import models, fitting, polynomial
from scipy.signal import argrelextrema
from copy import deepcopy
from joblib import Parallel, delayed
import multiprocessing
from time import time


def fit_gaussians_to_spectrum(obs_flx, obs_wvl, idx_peaks):
    print '  Fitting multiple ({:.0f}) gaussinas to the extracted arc spectrum'.format(len(idx_peaks))
    #
    median_val = np.median(obs_flx)
    mean_wvl = 0  # np.mean(obs_wvl)
    fit_model = models.Const1D(amplitude=np.median(obs_flx)) + polynomial.Polynomial1D(4)
    for i_p, idx_p in enumerate(idx_peaks):
        peak_val = obs_flx[idx_p] - median_val
        # dela boljs brez nastavljenih boundov
        fit_model += models.Gaussian1D(amplitude=peak_val, mean=obs_wvl[idx_p]-mean_wvl, stddev=0.1)#,
                                       # bounds={'mean': (obs_wvl[idx_p]-0.5, obs_wvl[idx_p]+0.5),
                                       #         'amplitude': (peak_val*0.8, peak_val*1.2)})
    fit_t = fitting.LevMarLSQFitter()
    fitted_model = fit_t(fit_model, obs_wvl-mean_wvl, obs_flx)
    return fitted_model


def fit_gaussians_to_spectrum_linebyline(obs_flx, obs_wvl, idx_peaks, d_wvl=2):
    print '  Fitting {:.0f} individual gaussinas to the extracted arc spectrum'.format(len(idx_peaks))
    fitted_vals = []  # [mean, std] combinations
    for i_p, idx_p in enumerate(idx_peaks):
        # data subset
        idx_data_use = np.logical_and(obs_wvl > obs_wvl[idx_p]-d_wvl, obs_wvl < obs_wvl[idx_p]+d_wvl)
        use_flx = obs_flx[idx_data_use]
        use_wvl = obs_wvl[idx_data_use]
        #
        median_val = np.median(use_flx)
        fit_model = models.Const1D(amplitude=np.median(use_flx))
        peak_val = obs_flx[idx_p] - median_val
        # dela boljs brez nastavljenih boundov
        fit_model += models.Gaussian1D(amplitude=peak_val, mean=obs_wvl[idx_p], stddev=0.1)#,
                                       # bounds={'mean': (obs_wvl[idx_p]-0.5, obs_wvl[idx_p]+0.5),
                                       #         'amplitude': (peak_val*0.8, peak_val*1.2)})
        fit_t = fitting.LevMarLSQFitter()
        fitted_model = fit_t(fit_model, use_wvl, use_flx)
        # print fitted_model.param_names
        fitted_vals.append(fitted_model.parameters[2:])
    return np.array(fitted_vals)


# TODO : automatically find those files or supply them as function inputs
res_list = glob('T*K2SNWNVR20N.fits')
res_list = [rf.split('.')[0] for rf in res_list]
# res_list = ['solar']
obs_list = ['EC60967','EC60969','EC60971','EC60973','EC60975','EC60977']

n_subsections = 31

# create a dictionary that will hold all results
R_dict_all = {}
for i_s in range(1, n_subsections+1):
    R_dict_all[str(i_s)] = {'R_val':list([]), 'R_wvl':list([])}
print R_dict_all

for obs_file in obs_list:
    R_line = []
    R_line_wvl = []

    print obs_file
    for i_s in range(1, n_subsections+1):
        print ' order', i_s
        try:
            print obs_file + '.ec_wvl_order{:02.0f}.txt'.format(i_s)
            obs_data = np.loadtxt(obs_file + '.ec_wvl_order{:02.0f}.txt'.format(i_s))
            # print obs_data
        except:
            print '  ', 'File for order {:.0f} not found'.format(i_s)
            continue
        obs_flx = obs_data[:, 1]
        obs_wvl = obs_data[:, 0]

        # determine obvious peaks
        med_flx = np.median(obs_flx)
        # print ' Median flx', med_flx
        idx_peaks = argrelextrema(obs_flx, np.greater, order=2)[0]
        # first filtering
        n_max_strongest = 30
        min_peak_flx = np.sort(obs_flx[idx_peaks])[::-1][n_max_strongest]
        idx_peaks = idx_peaks[np.logical_and(obs_flx[idx_peaks] > 1.5*med_flx, obs_flx[idx_peaks] >= min_peak_flx)]

        arc_fit = fit_gaussians_to_spectrum(obs_flx, obs_wvl, idx_peaks)
        # gaus_params = fit_gaussians_to_spectrum_linebyline(obs_flx, obs_wvl, idx_peaks, d_wvl=1.5)
        # print arc_fit

        # plt.plot(obs_wvl, obs_flx, color='black')
        # wvl_oversample = np.linspace(obs_wvl[0], obs_wvl[-1], len(obs_wvl)*5)
        # for idx_p in idx_peaks:
        #     plt.scatter(obs_wvl[idx_p], obs_flx[idx_p], color='red')
        # plt.plot(wvl_oversample, arc_fit(wvl_oversample), color='blue')
        # plt.show()
        # plt.close()

        # analyze lines and determine R values
        arc_fit_params = arc_fit.param_names
        R_wvl = []
        R_val = []
        for i_l in np.arange(len(idx_peaks))+2:
            wvl = arc_fit.parameters[arc_fit_params.index('mean_'+str(i_l))]
            if wvl < np.nanmin(obs_wvl) or wvl > np.nanmax(obs_wvl):
                continue
            std = arc_fit.parameters[arc_fit_params.index('stddev_'+str(i_l))]
            fwhm = 2.3548 * std
            # print wvl/fwhm
            R_wvl.append(wvl)
            R_val.append(wvl/fwhm)

        # R_val2 = gaus_params[:,0] / (2.3548 * gaus_params[:,1])

        plt.scatter(R_wvl, R_val, label='combined')
        # plt.scatter(gaus_params[:,0], R_val2, label='individual')
        plt.ylim(5000, 25000)
        plt.legend()
        # plt.show()
        plt.savefig(obs_file + '_order{:02.0f}.png'.format(i_s), dpi=250)
        plt.close()

        R_val_std = np.std(R_val)
        if R_val_std < 300:
            # something strange here - almost no variation in R values
            print '  Variation in R strangely low'
            continue

        R_dict_all[str(i_s)]['R_val'].append(R_val)
        R_dict_all[str(i_s)]['R_wvl'].append(R_wvl)

for i_s in range(1, n_subsections+1):
    x_vals = np.hstack(R_dict_all[str(i_s)]['R_wvl'])
    plt.scatter(x_vals, np.hstack(R_dict_all[str(i_s)]['R_val']), label='Multi arc', lw=0, s=2)
    # plt.scatter(gaus_params[:,0], R_val2, label='individual')
    plt.ylim(5000, 25000)
    plt.xlim(np.percentile(x_vals, 1.), x_vals(x_vals, 99.))
    plt.legend()
    # plt.show()
    plt.savefig('Combined_R_order{:02.0f}.png'.format(i_s), dpi=250)
    plt.close()
