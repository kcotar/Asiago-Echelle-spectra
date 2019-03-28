import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table, join
import glob

data_dir = '/shared/ebla/cotar/'
obs_dir = data_dir + 'Asiago_reduced_data/GREGOR_TEST_2/'

arcs = pd.read_csv('arc_lines.txt', delim_whitespace=True, header=None)#, format='ascii.fixed_width')
arcs_man = arcs[4].values
arcs_det = arcs[3].values

# print arcs.keys()
# print arcs[4].values

for txt_f in glob.glob(obs_dir+'EC59445.ec_wvl_order2*.txt'):
    print txt_f
    # arc_data = pd.read_csv(txt_f, delimiter='  ')
    arc_data = np.genfromtxt(txt_f, delimiter='  ')
    if len(arc_data) <= 0:
        continue
    # print arc_data
    plt.plot(arc_data[:, 0], arc_data[:, 1], color='black')
    for a1, a2 in zip(arcs_man, arcs_det):
        try:
            plt.axvline(np.float(a1), ls='--', color='C3')
            plt.axvline(np.float(a2), ls='--', color='C2')
        except:
            pass
    plt.xlim(np.min(arc_data[:, 0]), np.max(arc_data[:, 0]))
    plt.show()
    plt.close()
