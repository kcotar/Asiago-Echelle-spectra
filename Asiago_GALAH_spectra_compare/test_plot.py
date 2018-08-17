import numpy as np
import matplotlib.pyplot as plt
from sklearn.externals import joblib
from glob import glob
from astropy.io import fits

ref_dir = '/home/gregor/public_html/Astro_predmeti/normalized_spectra_R20000/'
res_list = glob(ref_dir+'*/T*M05*V000K2SNWNVR20N.ASC')
w_a = np.loadtxt(ref_dir+'LAMBDA_R20.DAT')

# f_g = joblib.load('/home/klemen/data4_mount/galah_dr53_ccd1_4710_4910_wvlstep_0.040_ext4_R20000_20180327.pkl')
# w_g = np.arange(4710., 4910., 0.04)

# f_g = joblib.load('/home/klemen/data4_mount/galah_dr53_ccd2_5640_5880_wvlstep_0.050_ext4_R20000_20180327.pkl')
# w_g = np.arange(5640., 5880., 0.05)
# idx_use = np.where(np.logical_and(w_g > 5680., w_g < 5850.))[0]

f_g = joblib.load('/home/klemen/data4_mount/galah_dr53_ccd3_6475_6745_wvlstep_0.060_ext4_R20000_20180327_sig.pkl')
w_g = np.arange(6475., 6745., 0.06)
idx_use = np.where(np.logical_and(w_g > 6490., w_g < 6720.))[0]

# f_g = joblib.load('/home/klemen/data4_mount/galah_dr53_ccd4_7700_7895_wvlstep_0.070_ext4_R20000_20180327.pkl')
# w_g = np.arange(7700., 7895., 0.07)
# idx_use = np.where(np.logical_and(w_g > 7700., w_g < 7850.))[0]

f_g_use = f_g[:,idx_use]
w_g_use = w_g[idx_use]

ref_wvl_orig = np.loadtxt(ref_dir + 'LAMBDA_R20.DAT')

obs_dir = '/home/klemen/Asiago_prepare/MCMC_Cannon_model_R20K/'
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

    # quick spectral distance
    # idx_sim_mask = np.where(obs_flx_use < 0.975)[0]
    # ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use)**2)[:, idx_sim_mask], axis=1))
    ref_dist = np.sqrt(np.sum(((f_g_use - obs_flx_use) ** 2), axis=1))
    idx_ref_dist = np.argsort(ref_dist)[:100]

    for i in idx_ref_dist:
        plt.plot(w_g_use, f_g_use[i,:], c='blue', alpha=0.05)
    plt.plot(w_g_use, obs_flx_use, c='black', alpha=1)
    plt.tight_layout()
    plt.show()
    plt.close()
