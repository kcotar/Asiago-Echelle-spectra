import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Lasso
from matplotlib.path import Path
from astropy.table import Table, hstack, vstack
from os import path
from glob import glob


class LassoManager(object):
    def __init__(self, ax, d, file_out):
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.o_p = file_out
        self.d = d
        self.p = d['bp_rp','G_abs'].to_pandas().values
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def callback(self, verts):
        polyg = Path(verts)
        ind = polyg.contains_points(self.p)
        print d[ind]
        d[ind]['source_id','ra','dec','phot_g_mean_mag', 'radial_velocity'].write(self.o_p, overwrite=True, format='ascii.tab')
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso

    def onpress(self, event):
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes,
                           (event.xdata, event.ydata),
                           self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

data_dir = '/shared/ebla/cotar/'
asiago_data_dir = data_dir + 'Asiago_reduced_data/ORION_CLUSTERS/'
observed = Table.read(asiago_data_dir + 'ORION_CLUSTERS_reduced.fits')

all_observed = list([])
all_to_be = list([])
for i in range(1, 6):
    print i
    fits = 'Orion_group_'+str(i)+'.fits'
    f_out = 'Orion_group_' + str(i) + '_sel.txt'

    d = Table.read(fits)

    mag_abs = 5. + 5. * np.log10(d['parallax'] / 1e3) + d['phot_g_mean_mag']
    d['G_abs'] = mag_abs

    idx_rv = np.isfinite(d['radial_velocity'])
    plt.scatter(d['bp_rp'], d['G_abs'], lw=0, s=7, label='', c='black')
    plt.scatter(d['bp_rp'][idx_rv], d['G_abs'][idx_rv], lw=0, s=7, label='RV', c=d['radial_velocity'][idx_rv])
    for iso_csv in glob('gaia_*M.dat'):
        iso_data = Table.read(iso_csv, format='ascii.csv', delimiter=' ')
        idx_plot = (iso_data['Mini']-iso_data['Mass']) / iso_data['Mass'] < 0.1
        plt.plot((iso_data['G_BPmag']-iso_data['G_RPmag'])[idx_plot], iso_data['Gmag'][idx_plot], lw=1, label=iso_csv.split('.')[0].split('_')[1])
    # plt.plot(zams['Bmag']-zams['Vmag'], zams['Vmag'], lw=1, c='C3')
    plt.colorbar()
    plt.gca().invert_yaxis()
    plt.xlabel('Bp - Rp')
    plt.ylabel('G mag')
    plt.tight_layout()
    plt.legend()
    # lman = LassoManager(plt.gca(), d, f_out)
    plt.savefig('Orion_group_'+str(i)+'_rv.png', dpi=300)
    # plt.show()
    plt.close()

    plt.scatter(d['bp_rp'], d['G_abs'], lw=0, s=7, label='')
    if path.isfile(f_out):
        sel_list = Table.read(f_out, format='ascii.tab')
        ind_sel = np.in1d(d['source_id'], sel_list['source_id'])
        plt.scatter(d['bp_rp'][ind_sel], d['G_abs'][ind_sel], lw=0, s=7, label='Selected')
    plt.gca().invert_yaxis()
    plt.xlabel('Bp - Rp')
    plt.ylabel('Absolute G mag')
    plt.tight_layout()
    plt.legend()
    plt.savefig(f_out[:-3]+'png', dpi=300)
    plt.close()

    f_out_img = 'Orion_group_'+str(i)+'_obs.txt'
    ind_obs = np.in1d(d['source_id'], observed['source_id'])
    d[ind_obs]['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'radial_velocity'].write(f_out_img, overwrite=True, format='ascii.tab')
    plt.scatter(d['bp_rp'], d['phot_g_mean_mag'], lw=0, s=7, label='')
    plt.scatter(d['bp_rp'][ind_obs], d['phot_g_mean_mag'][ind_obs], lw=0, s=7, label='Observed')
    plt.gca().invert_yaxis()
    plt.xlabel('Bp - Rp')
    plt.ylabel('G mag')
    plt.tight_layout()
    plt.legend()
    plt.savefig(f_out_img[:-3]+'png', dpi=300)
    plt.close()

    all_observed.append(d[np.logical_and(ind_sel, ind_obs)])
    all_to_be.append(d[np.logical_and(ind_sel, ~ind_obs)])

all_observed = vstack(all_observed)
all_to_be = vstack(all_to_be)
all_observed['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'radial_velocity'].write('Orion_groups_together_observed.txt', overwrite=True, format='ascii.tab')
all_to_be['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'radial_velocity'].write('Orion_groups_together_to_be_observed.txt', overwrite=True, format='ascii.tab')

