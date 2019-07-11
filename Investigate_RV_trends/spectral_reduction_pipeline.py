import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from os import path
from astropy.io import fits
from astropy.table import Table, join
from scipy.interpolate import interp2d, RectBivariateSpline
from astropy.modeling import models, fitting

def get_line_val(line_str):
    return np.float64(line_str.split('\t')[2])


def norm_variable(x, xmin, xmax):
    return (2. * x - (xmax + xmin)) / (xmax - xmin)


def get_curve(trace_lines, i_s):

  n_read = get_line_val(trace_lines[i_s])
  curve_id = get_line_val(trace_lines[i_s+1])
  curve_order = get_line_val(trace_lines[i_s+2])
  def_min = get_line_val(trace_lines[i_s+3])
  def_max = get_line_val(trace_lines[i_s+4])
  #print curve_order, def_min, def_max
  x_vals = np.linspace(def_min, def_max, 4000)

  # curve params
  curve_vals = list([])
  for i_c in range(int(curve_order)):
    curve_vals.append(get_line_val(trace_lines[i_s+5+i_c]))

  # construct the curve
  if curve_id == 1:  # cheb polynomial
    y_vals = np.full_like(x_vals, 0.)
    z_vals = list([])  # inital z values
    z_vals.append(np.full_like(x_vals, 1.))
    z_vals.append(norm_variable(x_vals, def_min, def_max))
    for i_p, c_p in enumerate(curve_vals):
      if i_p < 2:
        z_curr = z_vals[i_p]
      else:
        #print z_vals[-1], z_vals[-2]
        z_curr = 2. * norm_variable(x_vals, def_min, def_max) * z_vals[i_p-1] - z_vals[i_p-2]
        z_vals.append(z_curr)
      y_vals += c_p * z_curr
      #print y_vals#, i_p, c_p

  return x_vals, y_vals

# ----------------
# -----TESTS------
# ----------------
data_dir = '/shared/ebla/cotar/'
obs_dir = data_dir + 'Asiago_reduced_data/GREGOR_TEST_2/'
apall_file_txt = obs_dir + '/database/apEC59444'
fits_file = obs_dir + 'EC59444.fits'
arc_file = obs_dir + 'EC59445.fits'

txt = open(apall_file_txt, 'r')
apall_lines = txt.read().split('\n')
txt.close()

fits_img = fits.open(fits_file)
arc_img = fits.open(arc_file)

img_vals = fits_img[0].data
arc_vals = arc_img[0].data

img_size = img_vals.shape
print img_vals.shape, arc_vals.shape
f_img = interp2d(np.arange(img_size[0])+1, np.arange(img_size[1])+1, img_vals, kind='cubic')
f_arc = interp2d(np.arange(img_size[0])+1, np.arange(img_size[1])+1, arc_vals, kind='cubic')

img_vals_interp = f_img(np.arange(img_size[0])+1, np.arange(img_size[1])+1)

i_order = 26
for i_l, line in enumerate(apall_lines):
    if i_l < 923:
        continue
    if 'center' in line:
        n_center_line = i_l
        center_vals = line.split('\t')[-1].split(' ')
        x_off = np.float64(center_vals[-2])
        y_off = np.float64(center_vals[-1])
    if 'curve' in line:
        print 'Found curve at pos', i_l
        x_curve, y_curve = get_curve(apall_lines, i_l)
        y_curve += y_off
        x_curve += 0.

        # print x_curve
        # print y_curve
        # print f_img(x_curve, y_curve)
        # print f_img(x_curve, y_curve).shape

        # simple extraction of spectrum and arc line
        spec_line = [f_img(XX, YY) for XX, YY in zip(x_curve, y_curve)]
        arc_line = [f_arc(XX, YY) for XX, YY in zip(x_curve, y_curve)]
        # export arc and spec lines
        f_out = 'spec_extracted_{:04.0f}.csv'.format(i_l)
        csv_file = open(f_out, 'w')
        for vals in zip(spec_line, arc_line, x_curve, y_curve):
            csv_file.write(''.join([str(vv) for vv in vals])+'\n')
        csv_file.close()

        # show both plots
        plt.plot(x_curve, spec_line, label='spectrum')
        plt.plot(x_curve, arc_line, label='arc')
        plt.show()
        plt.close()

        # identify extracted lines
        arc_ident = 'arc_iden_{:04.0f}.csv'.format(i_l)
        if not path.isfile(arc_ident):

            all_x_wvl_pos = []
            # interactive matplotlib plot
            fig, ax = plt.subplots()
            ax.plot(x_curve, arc_line)

            def onclick(event):
                if event.dblclick:
                    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f', event.x, event.y, event.xdata, event.ydata)
                    # fit peak
                    idx_fit_use = np.abs(x_curve - event.xdata) < 3.5
                    # print x_curve[idx_fit_use], y_curve[idx_fit_use]
                    gg = models.Gaussian1D(amplitude=np.max(y_curve[idx_fit_use]),
                                           mean=event.xdata)
                    fitter = fitting.LevMarLSQFitter()
                    gg_fit = fitter(gg, y_curve[idx_fit_use], x_curve[idx_fit_use])
                    x_pos = gg_fit.mean.value
                    ax.axvline(x_pos)
                    fig.canvas.draw()
                    # time.sleep(2.)
                    # print x_pos
                    wvl_peak = np.float(raw_input('Wavelength at {:.2f}:'.format(event.x)))
                    # print wvl_peak
                    all_x_wvl_pos.append([x_pos, wvl_peak])
                return True

            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            plt.tight_layout()
            plt.show()
            fig.canvas.mpl_disconnect(cid)
            plt.close()

            # save inputs to a txt file
            csv_out = open(arc_ident, 'w')
            for xy_wvl in all_x_wvl_pos:
                csv_out.write(','.join([str(xx) for xx in xy_wvl])+'\n')
            csv_out.close()

            all_x_wvl_pos = np.array(all_x_wvl_pos)
        else:
            # read predetermined arc point
            all_x_wvl_pos = np.genfromtxt(arc_ident, delimiter=',')
            print all_x_wvl_pos

        if len(all_x_wvl_pos) <= 0:
            continue
            
        # fit polynomial to the arc points
        gg = models.Polynomial1D(degree=4)
        fitter = fitting.LevMarLSQFitter()
        gg_fit = fitter(gg, all_x_wvl_pos[:, 0], all_x_wvl_pos[:, 1])
        print gg_fit
        diff_wvl = all_x_wvl_pos[:, 1] - gg_fit(all_x_wvl_pos[:, 0])
        print 'Diff: {:.2f} +/- {:.4f}'.format(np.nanmean(diff_wvl), np.nanstd(diff_wvl))

        # print arc wvl points
        x_vals = np.linspace(np.min(all_x_wvl_pos[:, 0]), np.max(all_x_wvl_pos[:, 0]), 250)
        plt.plot(x_vals, gg_fit(x_vals))
        plt.scatter(all_x_wvl_pos[:, 0], all_x_wvl_pos[:, 1])
        plt.title('Order '+str(i_order)+' calibration')
        plt.xlabel('Pixel')
        plt.ylabel('Wavelength')
	plt.tight_layout()
        plt.show()
        plt.close()

        # read the correct spectrum
        flx_file = obs_dir + 'EC59444.ec_order{:02.0f}.txt'.format(i_order)
        flx_iraf = np.genfromtxt(flx_file, delimiter='  ')

        plt.plot(flx_iraf[:, 0], flx_iraf[:, 1]/np.max(flx_iraf[:, 1]), label='IRAF complex reduction')
        plt.plot(gg_fit(x_curve), spec_line/np.max(spec_line), label='Custom fast per order reduction')
        plt.title('Order '+str(i_order))
        plt.xlabel('Wavelength')
        plt.ylabel('Max norm flux')
	plt.tight_layout()
        plt.legend()
        plt.show()
        plt.close()

        i_order += 1

