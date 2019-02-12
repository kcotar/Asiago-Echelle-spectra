import numpy as np
import matplotlib.pyplot as plt
import pyfits
from os import chdir
from scipy.interpolate import interp2d, RectBivariateSpline
from astropy.modeling import models, fitting
from lmfit.models import ConstantModel, GaussianModel

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
  x_vals = np.arange(def_min, def_max+1)

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

apall_file_txt = 'observations/TEST/database/apEC61329'
fits_file = 'observations/TEST/EC61329.fits'

#apall_file_txt = 'observations/GAIA_RV_STAND/database/apEC60970'
#fits_file = 'observations/GAIA_RV_STAND/EC60970.fits'

#apall_file_txt = 'observations/NGC6940/database/apEC61259'
#fits_file = 'observations/NGC6940/EC61259.fits'

txt = open(apall_file_txt, 'r')
apall_lines = txt.read().split('\n')
txt.close()

fits_img = pyfits.open(fits_file)
img_vals = fits_img[0].data

img_size = img_vals.shape
print img_size
f_img = interp2d(np.arange(img_size[0])+1, np.arange(img_size[1])+1, img_vals, kind='cubic')
img_vals_interp = f_img(np.arange(img_size[0])+1, np.arange(img_size[1])+1)
#f_img = RectBivariateSpline(np.arange(img_size[0]), np.arange(img_size[1]), img_vals, kx=3, ky=3)
print f_img(1264., 1068.)

plt.imshow(img_vals, origin='lower', vmin=10, vmax=np.percentile(img_vals, 90), interpolation='none', extent=(0.5, img_size[0]+0.5,0.5, img_size[0]+0.5))
#plt.show()
#plt.close()

for i_l, line in enumerate(apall_lines):
  if 'center' in line:
    n_center_line = i_l
    center_vals = line.split('\t')[-1].split(' ')
    x_off = np.float64(center_vals[-2])
    y_off = np.float64(center_vals[-1])
    #print x_off, y_off
  if 'curve' in line:
    #if y_off <= 1800:
    #  continue
    print 'Found curve at pos', i_l
    x_curve, y_curve = get_curve(apall_lines, i_l)
    y_curve += y_off
    x_curve += 0.
    
    #img_interp = f_img(np.linspace(1,2048,5000), np.linspace(np.min(y_curve),np.max(y_curve),250))
    #plt.imshow(img_interp, origin='lower', vmin=20, vmax=np.percentile(img_interp, 95), interpolation='none')
    #plt.show()
    #plt.close()

    '''
    x_multi_pos = np.linspace(np.min(x_curve), np.max(x_curve), 25)
    for x_p in [920]:#x_multi_pos:
      idx_x = np.where(x_curve == int(x_p))[0][0]
      print x_p, idx_x
      y_pos = np.int(y_curve[idx_x])
      #y_pos = np.linspace(y_pos-10, y_pos+10, 120) 
      y_pos = np.arange(y_pos-10, y_pos+10,)
      img_vals_interp = [f_img(np.int(x_curve[idx_x]),YY) for YY in y_pos]

      gg = models.Gaussian1D(amplitude=np.max(img_vals_interp), mean=y_curve[idx_x], stddev=1.5, bounds={"stddev_0": (0.1, 3.)}) + models.Const1D(amplitude=50, bounds={"amplitude_1": (0, 100)})
      fitter = fitting.LevMarLSQFitter()
      gg_fit = fitter(gg, y_pos, img_vals_interp)
      print gg_fit


      mod = GaussianModel()
      pars = mod.guess(img_vals_interp, x=y_pos)
      out = mod.fit(img_vals_interp, pars, x=y_pos)
      print out.fit_report()

      plt.plot(y_pos, img_vals_interp)
      plt.plot(y_pos, gg_fit(y_pos))
      plt.plot(y_pos, out.best_fit)
      plt.axvline(y_curve[idx_x])
      plt.title(str(x_curve[idx_x]))
      plt.show()
      plt.close()
    '''

    '''
    y_shifts = np.linspace(-2., 2., 91)
    x_multi_ranges = np.linspace(1,len(x_curve), int(len(x_curve)/50.))
    for i_x_r in range(len(x_multi_ranges)-1):
      #print ' X range:', i_x_r+1
      idx_x = np.logical_and(x_curve >= x_multi_ranges[i_x_r], x_curve <= x_multi_ranges[i_x_r+1])
      sum_shifts = list([])
      for y_shift in y_shifts:
        img_vals_interp = [f_img(XX,YY) for XX,YY in zip(x_curve[idx_x], y_curve[idx_x]+y_shift)]
        sum_shifts.append(np.sum(np.float64(img_vals_interp)))
        #plt.plot(img_vals_interp, label=str(y_shift))
      # change argmax with some kind of fit
      idx_max = np.argmax(sum_shifts)
      y_shift_best = y_shifts[idx_max]
      print '  Best Y corr', y_shift_best, sum_shifts[idx_max]
      plt.plot(x_curve[idx_x], y_curve[idx_x]+y_shift_best, c='red')
    #plt.axvline(y_shift_best, color='C1')
    #plt.plot(y_shifts, sum_shifts)
    #plt.legend()
    #plt.tight_layout()
    #plt.show()
    #plt.close()
    '''   

    '''
    sum_shifts = list([])
    interp_img_new = list([])
    y_shifts = np.linspace(-1., 1., 81)
    for y_shift in y_shifts:
      img_vals_interp = [f_img(XX,YY) for XX,YY in zip(x_curve, y_curve+y_shift)]
      sum_shifts.append(np.sum(img_vals_interp))
      interp_img_new.append(img_vals_interp)
      #plt.plot(img_vals_interp, label=str(y_shift))
    # change argmax with some kind of fit
    y_shift_best = y_shifts[np.argmax(sum_shifts)]
    print ' y corr', y_shift_best
    #plt.axvline(y_shift_best, color='C1')
    #plt.plot(y_shifts, sum_shifts)
    #plt.legend()
    #plt.imshow(np.array(interp_img_new)[:,:,0], origin='lower', vmin=10, vmax=np.percentile(interp_img_new, 95), interpolation='none')
    #plt.tight_layout()
    #plt.show()
    #plt.close()   
    '''

    plt.plot(x_curve, y_curve, c='black')
    #plt.plot(x_curve, y_curve+y_shift_best, c='red')
plt.xlim(0, 2048)
plt.ylim(0, 2048)
plt.tight_layout()
plt.show()
plt.close()
