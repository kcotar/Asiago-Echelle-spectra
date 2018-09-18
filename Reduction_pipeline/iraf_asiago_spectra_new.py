# -*- encoding: utf-8 -*-
#!/usr/bin/env python
import sys
import os
import glob
import pyfits
from pyraf import iraf
import shutil
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import operator
import time
import traceback
from lmfit import minimize, Parameters, report_fit

#import cosmics

iraf.noao(_doprint=0,Stdout="/dev/null")
iraf.rv(_doprint=0,Stdout="/dev/null")
iraf.imred(_doprint=0,Stdout="/dev/null")
iraf.ccdred(_doprint=0,Stdout="/dev/null")
iraf.images(_doprint=0,Stdout="/dev/null")
iraf.immatch(_doprint=0,Stdout="/dev/null")
iraf.onedspec(_doprint=0,Stdout="/dev/null")
iraf.twodspec(_doprint=0,Stdout="/dev/null")
iraf.apextract(_doprint=0,Stdout="/dev/null")
iraf.imutil(_doprint=0,Stdout="/dev/null")
iraf.echelle(_doprint=0,Stdout="/dev/null")
iraf.astutil(_doprint=0,Stdout="/dev/null")
iraf.apextract.dispaxi=1
iraf.echelle.dispaxi=1
#fixes a bug with latest versions of iraf
iraf.ccdred.instrum='blank.txt'


os.environ['PYRAF_BETA_STATUS'] = '1'

# REFS NEW: ecEC59550 ecEC59758 ecEC59795 ecEC59842 ecEC59844 ecEC59864 ecEC59866 ecEC59881 ecEC59883 ecEC59885

observations = {
  'cemp_cand': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/CEMP_obs/201712/',
    'biases':['EC60040', 'EC60041', 'EC60042', 'EC60043', 'EC60044'], 
    'flats':['EC60046', 'EC60047', 'EC60048'], 
    'objects': ['EC60089'], 
    'calibs':  ['EC60090'], 
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC59842.ec', 
    'REF_AP':'', 
    'ap_position':'left'
  },
  'bright': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Bright_stars/201804/',
    'biases':['bias1_1', 'bias1_2', 'bias1_3', 'bias1_4', 'bias1_5', 'bias2_1', 'bias2_2', 'bias2_3', 'bias2_4', 'bias2_5'], 
    'flats':['flat2_1', 'flat2_2', 'flat2_3', 'flat2_4', 'flat2_5'], 
    'objects': ['EC60835','EC60837','EC60839','EC60869','EC60871','EC60896','EC60898','EC60900'], 
    'calibs':  ['EC60836','EC60838','EC60840','EC60870','EC60872','EC60897','EC60899','EC60901'], 
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC59881.ec', 
    'REF_AP':'', 
    'ap_position':'center'
  },
  'rv_stand2': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/RV_standards/',
    'biases':['EC60944', 'EC60945', 'EC60946', 'EC60947', 'EC60948'], 
    'flats':['EC60949', 'EC60950', 'EC60951', 'EC60952', 'EC60953'], 
    'objects': ['EC60966','EC60968','EC60970','EC60972','EC60974','EC60976'], 
    'calibs': ['EC60967','EC60969','EC60971','EC60973','EC60975','EC60977'], 
#    'objects': ['EC61147'], 
#    'calibs': ['EC61148'], 
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  },
  'rv_stand3': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/RV_standards/',
    'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
    'flats':['EC61160', 'EC61161', 'EC61162'], 
    'objects': ['EC61147'], 
    'calibs': ['EC61148'], 
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  },
  'maxi': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/MAXI/',
    'biases':['EC60944', 'EC60945', 'EC60946', 'EC60947', 'EC60948'], 
    'flats':['EC60949', 'EC60950', 'EC60951', 'EC60952', 'EC60953'], 
    'objects': ['EC60982'], 
    'calibs':  ['EC60983'], 
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  },
# Solar twins candidates from here on
  'twin_1': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Twins_obs/201804/',
    'biases':['EC60791', 'EC60792', 'EC60793', 'EC60794', 'EC60795'], 
    'flats':['EC60820', 'EC60821', 'EC60822'], 
    'objects': ['EC60815','EC60817','EC60819'], 
    'calibs':  ['EC60816','EC60818','EC60818'], 
    'sobject_id': ['160524006601258'],
    'obj_name': ['TYC 502-985-1'],
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP':'', 
    'ap_position':'center'
  },
# Solar twins multiples only
  'twin_2': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Twin_obs/201806/',
    'biases':['EC60944', 'EC60945', 'EC60946', 'EC60947', 'EC60948'], 
    'flats':['EC60949', 'EC60950', 'EC60951', 'EC60952', 'EC60953'], 
    'objects': ['EC60984'], 
    'calibs':  ['EC60985'], 
    'sobject_id': ['160531006101153'],
    'obj_name': ['TYC 489-2258-1'],
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  },
  'twin_3': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Twin_obs/201807/',
    'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
    'flats':['EC61160', 'EC61161', 'EC61162'], 
    'objects': ['EC61131', 'EC61133', 'EC61135'], 
    'calibs':  ['EC61132', 'EC61134', 'EC61134'], 
    'sobject_id': ['160531006101153'],
    'obj_name': ['TYC 489-2258-1'],
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  },
  'twin_4': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Twin_obs/201807/',
    'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
    'flats':['EC61160', 'EC61161', 'EC61162'], 
    'objects': ['EC61136'], 
    'calibs':  ['EC61137'], 
    'sobject_id': ['160531005101256'],
    'obj_name': ['TYC 5805-1373-1'],
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  }
  'twin_5': {'ORIG_DIR':'/home/nandir/IRAF_Echelle_Asiago/Twin_obs/201807/',
    'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
    'flats':['EC61160', 'EC61161', 'EC61162'], 
    'objects': ['EC61138'], 
    'calibs':  ['EC61139'], 
    'sobject_id': ['170912002901303'],
    'obj_name': ['2MASS J23520336+1617418'],
    'REF_ARC':'/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/ecEC60975.ec', 
    'REF_AP': False, 
    'ap_position':'center'
  }

  # after and including april 2017, EC59550 is the new cal ref file because the orders are arranged differently on the ccd after upgrades 
}

# batches that will be processed
curr_work_dirs = ['rv_stand2']

c = 299792.458
AP_LIM = {'left': [-20, 10, '-19:-11'], 'center': [-15, 15, '-15:-10,10:15']}


def main():

  for workdir in curr_work_dirs:
    work_path = 'observations/' + workdir
    
    if not os.path.exists(work_path): 
      os.makedirs(work_path)

    if not os.path.exists(work_path+'/database'): 
      os.makedirs(work_path+'/database')

    if not os.path.exists(work_path+'/cosmics'): 
      os.makedirs(work_path+'/cosmics')

    # don't ask, we just need this
    os.system("touch %s/blank.txt" % work_path)
    
    # enter working folder, all the way to the end
    iraf.cd(work_path)
    
    # maybe useful, maybe not, not sure
    iraf.prcacheOff() 
    
    """
    # main steps, should work individually
    #check_spectra(workdir, '.ec.vh') 
    create_masterbias(workdir)    
    create_masterflat(workdir)    
    #do_ref_apall(workdir)    
    #normalize_masterflat(workdir)  
    clean_folder_write_header(workdir)      
    bias_correct(workdir)  
    #flat_correct(workdir)  
    #remove_cosmics(workdir)  
    
    do_apall_for_real(workdir)   
    apply_apall(workdir)  
    wavelength_solution(workdir) 
    #check_wav_solution(workdir)
    apply_wav_solution(workdir)   
    vhelio_correct(workdir)   
    #check_spectra(workdir, '.ec.vh') 
    combine_normalize_images(workdir, combine=False)
    """

    #get_RV(workdir, multi_ref=False, multi_sample=False)
    apply_wav_solution(workdir, include_cal=True)
    export_spectrum_to_txt(workdir, suffix='.ec_wvl', export_cal=True)
    export_spectrum_to_txt(workdir, suffix='_vh_norm')

    #check_spectra(workdir) 
    #new_flatfield_lens_problems(workdir)

    # back to mother folder
    iraf.cd('../')
    iraf.cd('../')  

def create_masterbias(WORK_DIR):
  print "\n + Creating masterbias\n"

  try: os.remove('masterbias.fits')
  except: pass

  iraf.zerocombine(input=','.join([observations[WORK_DIR]['ORIG_DIR']+i for i in observations[WORK_DIR]['biases']]), output='masterbias.fits', combine='median', reject='none', ccdtype='', process='no', Stdout="/dev/null")

def create_masterflat(WORK_DIR):
  print '\n + Creating masterflat\n'

  try: os.remove('masterflat.fits')
  except: pass

  # subtract masterbias from flats
  for flat in observations[WORK_DIR]['flats']:
    try: os.remove(flat+'.fits')
    except: pass
    iraf.imarith(operand1=observations[WORK_DIR]['ORIG_DIR']+flat, op='-', operand2='masterbias.fits', result=flat, Stdout="/dev/null")
  
  iraf.zerocombine(input=','.join([flat for flat in observations[WORK_DIR]['flats']]), output='masterflat.fits', combine='median', reject='none', ccdtype='', process='no')

  for flat in observations[WORK_DIR]['flats']:
    os.remove(flat+'.fits')

def do_ref_apall(WORK_DIR):
  print '\n + Doing apall for a reference to the masterflat\n'

  #do quick apall
  #This will serve as the reference where the apertures are in the process of creating a flat. 30 apertures are found automatically, including the red-most aperture that goes off the image. Change the position of apertures if needed. It is assumed that the object is positioned at roughly 1/3 from the left edge of the slit (1/3 from the right edge of the aperture).
  # d - briši aperturo
  # m - dodaj aperturo
  # o - posortiraj aperture po vrsti na novo
  # z - odstrani interval fitanja
  # s - levo in desno, dodaj interval fitanja

  iraf.unlearn('apall')

  try: os.remove(observations[WORK_DIR]['objects'][0]+'.fits')
  except: pass

  iraf.imcopy(input=observations[WORK_DIR]['ORIG_DIR']+observations[WORK_DIR]['objects'][0], output=observations[WORK_DIR]['objects'][0])
  
  iraf.apall(input=observations[WORK_DIR]['objects'][0], referen='', format='echelle', interac='yes', find='yes', recente='yes', resize='no', edit='yes', trace='yes', fittrac='yes',extract='no', extras='no', review='no', line=550, nsum=10, lower=-6, upper=6, width=11, radius=11, thresho=0.1, nfind=32, minsep=20, maxsep=155, t_func='chebyshev', t_order=4, t_sampl='51:2098', t_niter=5, backgro='none', ylevel='INDEF', llimit=AP_LIM[observations[WORK_DIR]['ap_position']][0], ulimit=AP_LIM[observations[WORK_DIR]['ap_position']][1])
  iraf.apresize(input=observations[WORK_DIR]['objects'][0], interac='no', find='no', recente='no', resize='yes', edit='yes', ylevel='INDEF', llimit=AP_LIM[observations[WORK_DIR]['ap_position']][0], ulimit=AP_LIM[observations[WORK_DIR]['ap_position']][1])

def normalize_masterflat(WORK_DIR):
  print '\n + Normalizing masterflat.\n'

  try: os.remove('masterflat_norm.fits')
  except: pass

  iraf.apflatten(input='masterflat.fits', output='masterflat_norm.fits', referen=observations[WORK_DIR]['objects'][0], interac='no', find='no', recenter='no', resize='no', edit='no', trace='no', fittrac='no', flatten='yes', fitspec='yes', functio='spline3', order=13, sample='51:2098', niterat=5)

def clean_folder_write_header(WORK_DIR):
  print '\n + Cleaning folder, truncating edges, fixing bad pixels\n'

  for obj in observations[WORK_DIR]['objects']:
    try: os.remove(obj+'.fits')
    except: pass
    try: os.remove(obj+'.ec.fits')
    except: pass
    iraf.imcopy(input=observations[WORK_DIR]['ORIG_DIR']+obj, output=obj) 

  for cal in observations[WORK_DIR]['calibs']:
    try: os.remove(cal+'.fits')
    except: pass  
    iraf.imcopy(input=observations[WORK_DIR]['ORIG_DIR']+cal, output=cal)

  with open('badpix.txt', 'w+') as file:
    file.write('704 704 1262  2048\n703 705 1262  1277')  

  iraf.hedit(images=','.join(observations[WORK_DIR]['objects']+observations[WORK_DIR]['calibs']), fields='CCDSEC', value='[51:2098,1:2048]', add='yes', verify='no')
  iraf.hedit(images=','.join(observations[WORK_DIR]['objects']+observations[WORK_DIR]['calibs']), fields='DATASEC', value='[51:2098,1:2048]', add='yes', verify='no') 

  # FIX badpix
  iraf.ccdproc(images=','.join(observations[WORK_DIR]['objects']+observations[WORK_DIR]['calibs']), ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', fixfile='badpix.txt')

def bias_correct(WORK_DIR):
  print '\n + Correcting images and arcs for bias\n'

  iraf.ccdproc(images=','.join(observations[WORK_DIR]['objects']+observations[WORK_DIR]['calibs']), ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='yes', darkcor='no', flatcor='no', zero='masterbias.fits', fixfile='badpix.txt')

def flat_correct(WORK_DIR):
  print '\n + Correcting images and arcs for flat\n'

  iraf.ccdproc(images=','.join(observations[WORK_DIR]['objects']+observations[WORK_DIR]['calibs']), ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='yes', flat='masterflat_norm.fits', fixfile='badpix.txt')

def remove_cosmics(WORK_DIR):
  print '\n + Removing cosmics\n'
  
  overdo = 5
  
  for obj in observations[WORK_DIR]['objects']+['mastercalib']:
    obj = obj+'.fits'

    rogain = iraf.hselect(images=obj, fields="GAIN", exp='yes', Stdout=1)
    ronoise = iraf.hselect(images=obj, fields="RDNOISE", exp='yes', Stdout=1)
    array, header = cosmics.fromfits(obj)
    exptime = float(iraf.hselect(images=obj, fields='EXPTIME', Stdout=1, expr='yes')[0])

    if exptime <= 600:
      c = cosmics.cosmicsimage(array, gain=float(rogain[0]), readnoise=float(ronoise[0]), sigclip=30.0+overdo, sigfrac=0.5, objlim=3.0)
    elif exptime <= 1200:
      c = cosmics.cosmicsimage(array, gain=float(rogain[0]), readnoise=float(ronoise[0]), sigclip=20.0+overdo, sigfrac=0.5, objlim=3.0)
    else:
      c = cosmics.cosmicsimage(array, gain=float(rogain[0]), readnoise=float(ronoise[0]), sigclip=15.0+overdo, sigfrac=0.5, objlim=3.0)
    
    c.run(maxiter = 4, xb=5.0, yb=3.0)
    cosmics.tofits(obj, c.cleanarray, header)
    cosmics.tofits("cosmics/" + obj, c.mask, header)

  iraf.imsum(input='cosmics/*.fits', output='cosmics/sum.fits', option='sum')
  os.system('ds9 cosmics/sum.fits')

def do_apall_for_real(WORK_DIR):
  print '\n + Doing apall for real\n'

  #do apall for real with background and all
  #This tracing will be used from now on. Pay attention to the background subtraction. You will have to verify each tracing on each reference image (first images in the log.txt).   
  try:
    shutil.copy(observations[WORK_DIR]['REF_AP'], 'database/' + observations[WORK_DIR]['REF_AP'].split('/')[-1])  
  except:
    pass

  for obj in observations[WORK_DIR]['objects']:
    if observations[WORK_DIR]['REF_AP']:      
      ref = observations[WORK_DIR]['REF_AP'].split('/')[-1].replace('ap', '')
    else:
      ref = obj

    iraf.apall(input=obj, referen=ref, format='echelle', interac='yes', find='yes', recente='yes', resize='yes', edit='yes', trace='yes', fittrac='yes', extract='no', extras='no', review='no', line=600, nsum=30, lower=-6, upper=6, width=11, radius=11, thresho=0.15, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro='median', b_order=1, b_sampl=AP_LIM[observations[WORK_DIR]['ap_position']][2], ylevel='INDEF', llimit=-5, ulimit=5, weights='none', clean='yes', order='increasing')#, apidtable='apert.txt')

def apply_apall(WORK_DIR):
  print '\n + Applying apall to rest of the objects and arcs\n'

  for obj in observations[WORK_DIR]['objects']:
    print ' - Processing image %s \n' % obj

    try: 
      os.remove(obj+'.ec.fits')
      os.remove(obj+'.ec.nosky.fits')
    except: pass 
    
    for j in [['.ec', 'median'], ['.ec.nosky', 'none']]:
      iraf.apall(input=obj, output=obj+j[0], format='echelle', referen=obj, interac='no', find='yes', recente='no', resize='no', edit='yes', trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro=j[1], b_order=1, b_sampl=AP_LIM[observations[WORK_DIR]['ap_position']][2], ylevel='INDEF', llimit=-5, ulimit=5, weights='none', clean='yes', order='increasing')#, apidtable='apert.txt')  
 
  for i, cal in enumerate(observations[WORK_DIR]['calibs']):
    print ' - Processing arc %s \n' % cal

    try: os.remove(cal+'.ec.fits')
    except: pass 
    print '   input:',cal,'    referen:',observations[WORK_DIR]['objects'][i]
    iraf.apall(input=cal, format='echelle', referen=observations[WORK_DIR]['objects'][i], interac='no', find='no', recente='no', resize='no', edit='no',trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro='none', ylevel='INDEF', llimit=-5, ulimit=5, order='increasing')#, apidtable='apert.txt')    
    
def check_spectra(WORK_DIR, ext):
  print "\n + Visually check spectra before proceeding\n"

  for obj in observations[WORK_DIR]['objects']:
    iraf.splot(obj+ext)

def wavelength_solution(WORK_DIR, skip_exist=False):
  print '\n + Finding the wavelength solution\n'

  # calibrate wavelength
  # This step will require you to calibrate first arc by hand. The rest will be done automatically.
  # after l (automatic find other lines) change x and y order to 4
  # x - choose plots of the fit (per order, pixel, ...)

  iraf.unlearn('ecreidentify')
  iraf.unlearn('ecidentify')

  if observations[WORK_DIR]['REF_ARC']:
    try:
      print 'shutil.copy'
      shutil.copy(observations[WORK_DIR]['REF_ARC'], 'database/' + observations[WORK_DIR]['REF_ARC'].split('/')[-1])  
    except:
      pass
  
    for cal in observations[WORK_DIR]['calibs']:
      print 'iraf.ecreident'
      iraf.ecreident(images=cal+'.ec', referenc=observations[WORK_DIR]['REF_ARC'].replace('/ec', '/').split('/')[-1], refit='yes', shift=0, cradius=5, thresho=10)
      pass
    
  for cal in observations[WORK_DIR]['calibs']:
    print 'iraf.ecident'
    iraf.ecident(images=cal+'.ec', coordli='linelists$thar.dat', match=1, maxfeat=1800, ftype='emission', fwidth=4, cradius=5, thresho=2, minsep=2, functio='chebyshev', xorder=3, yorder=3, niterat=5, lowreje=3, highreje=3, autowri='yes')
  
def check_wav_solution(WORK_DIR):
  print '\n + Check lines in the wavelength solution\n'

  for cal in observations[WORK_DIR]['calibs']:
    with open('database/ec%s.ec' % cal, 'rb') as f:
      lines = f.readlines()
      for i, j in enumerate(lines[::-1]):
        if 'features' in j:
          break
      diff = []
      for line in lines[len(lines)-i:-1]:
        if 'offset' not in line:
          try:
            diff.append(float(line.split()[3])-float(line.split()[4]))
          except:
            print_exception()
      
      fig = plt.figure()
      plt.scatter(range(len(diff)), diff, marker='x', s=15, c='black')
      plt.axhline(np.std(diff), c='red')
      plt.axhline(-np.std(diff), c='red')
      plt.title(cal)
      plt.show()
      #fig.savefig("residuals_%s.pdf" % cal)

def apply_wav_solution(WORK_DIR, include_cal=False):
  print '\n + Applying the wavelength solution\n'
  
  iraf.unlearn('dispcor')

  for i, obj in enumerate(observations[WORK_DIR]['objects']):
    cal = observations[WORK_DIR]['calibs'][i]+'.ec'
    # export wav calibrated spectra
    print 'Using cal:', cal
    for j in ['.ec', '.ec.nosky']:
      #iraf.hedit(images=obj, fields='refspec1', value=observations[WORK_DIR]['calibs'][0]+'.ec', add='yes', verify='no')
      iraf.refspectra(input=obj+j, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=obj+j, output=obj+j, lineari='no', verbose='yes')
    if include_cal:
      os.remove(cal+'_wvl.fits')
      # export wav calibrated cal image
      iraf.refspectra(input=cal, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=cal, output=cal+'_wvl', lineari='no', verbose='yes')
      
def vhelio_correct(WORK_DIR):
  print '\n + VHELIO correction\n'
    
  for obj in observations[WORK_DIR]['objects']:
    obj = obj+'.ec'

    utm = iraf.hselect(images=obj, fields='UTMIDDLE', Stdout=1, expr='yes')
    year = utm[0][:4]
    month = utm[0][5:7]
    day = utm[0][8:10]
    h = int(utm[0][11:13])
    m = int(utm[0][14:16])
    s = int(utm[0][17:19])

    ra = iraf.hselect(images=obj, fields='RA', Stdout=1, expr='yes')[0].split(':')
    dec = iraf.hselect(images=obj, fields='DEC', Stdout=1, expr='yes')[0].split(':')
    ra = float(ra[0])+int(ra[1])/60.+float(ra[2])/3600.
    dec = float(dec[0])+int(dec[1])/60.+float(dec[2])/3600.
    
    shutil.copy(obj+'.fits', obj+'.vh.fits')
    iraf.hedit(images=obj+'.vh', fields='UT', value=h+m/60.+s/3600., add='yes', verify='no')
    iraf.hedit(images=obj+'.vh', fields='EPOCH', value=year, add='yes', verify='no')
    iraf.rvcorrect(images=obj+'.vh', imupdat='yes', epoch=year, observa='ekar', year=year, month=month, day=day, ut=h+m/60.+s/3600., ra=ra, dec=dec)
    iraf.dopcor(input=obj+'.vh', redshift='-vhelio', isveloc='yes', dispers='yes')

    obj = obj+'.nosky'
    shutil.copy(obj+'.fits', obj+'.vh.fits')
    iraf.hedit(images=obj+'.vh', fields='UT', value=h+m/60.+s/3600., add='yes', verify='no')
    iraf.hedit(images=obj+'.vh', fields='EPOCH', value=year, add='yes', verify='no')
    iraf.rvcorrect(images=obj+'.vh', imupdat='yes', epoch=year, observa='ekar', year=year, month=month, day=day, ut=h+m/60.+s/3600., ra=ra, dec=dec)
    iraf.dopcor(input=obj+'.vh', redshift='-vhelio', isveloc='yes', dispers='yes')

def combine_normalize_images(WORK_DIR, combine=False):
  
  if combine:
    print '\n + Combine images\n'
    try: 
        os.remove('combined_sum.fits')       
    except: pass 
    iraf.scombine(input=','.join([obj+'.ec.vh' for obj in observations[WORK_DIR]['objects']]), output='combined_sum', group='apertures', combine='sum', reject='none', Stdout="/dev/null")

  # normalize spectra
  # This step will require you to manually normalize all the spectra.
  print '\n + Normalize spectra\n'

  if combine:
    try: 
        os.remove('combined_sum_cone_rv_echellet.fits')       
    except: pass 
    iraf.continuum(input='combined_sum', output='combined_sum_cont', type='fit', replace='no', listonly='no', functio='cheb', order=13, low_rej=2, high_rej=3, naverag=-3, niter=9, interac='no', markrej='no')
  else:  
    for obj in observations[WORK_DIR]['objects']:
      try: 
        os.remove(obj+'_cont.fits')
        os.remove(obj+'_vh_norm.fits')
      except: pass 
      iraf.continuum(input=obj+'.ec.vh', output=obj+'_cont', type='fit', replace='no', listonly='no', functio='cheb', order=13, low_rej=2, high_rej=3, naverag=-3, niter=9, interac='no', markrej='no', ask='yes')
      iraf.sarith(input1=obj+'.ec.vh', op='/', input2=obj+'_cont', output=obj+'_vh_norm', format='multispec', Stdout="/dev/null")

  #combine apertures
  print '\n + Combine apertures\n'

  if combine:
    try: 
        os.remove('combined_sum_data.fits')       
        os.remove('combined_sum_norm.fits')    
        os.remove('combined_final.0001.fits')       
    except: pass
    iraf.scombine(input='combined_sum', output='combined_sum_data', group='all', combine='sum', reject='none', Stdout="/dev/null")
    iraf.scombine(input='combined_sum_cont', output='combined_sum_norm', group='all', combine='sum', reject='none', Stdout="/dev/null")

    iraf.sarith(input1='combined_sum_data', op='/', input2='combined_sum_norm', output='combined_final', format='onedspec', Stdout="/dev/null")
  else:  
    for obj in observations[WORK_DIR]['objects']:
      try: 
        os.remove(obj+'_data1D.fits')
        os.remove(obj+'_cont1D.fits')
        os.remove(obj+'_1D_vh_norm.0001.fits')   
        os.remove(obj+'_1D_vh_norm.fits')        
      except: pass 
      iraf.scombine(input=obj+'.ec.vh', output=obj+'_data1D', group='all', combine='sum', reject='none', Stdout="/dev/null")
      iraf.scombine(input=obj+'_cont', output=obj+'_cont1D', group='all', combine='sum', reject='none', Stdout="/dev/null")

      iraf.sarith(input1=obj+'_data1D', op='/', input2=obj+'_cont1D', output=obj+'_1D_vh_norm', format='onedspec', Stdout="/dev/null")

def get_RV(WORK_DIR, multi_ref=False, multi_sample=False):
  # list RV ref files
  if multi_ref:
    rv_ref_list = glob.glob('../../*_Echelle_RVS.fits')
    rv_ref_list = [rf.split('/')[-1].split('.')[0] for rf in rv_ref_list]
  else:
    rv_ref_list = ['solar']
    # rv_ref_list = ['solar_spectra_conv']
  # print '  RV ref files list:', rv_ref_list

  for obj in observations[WORK_DIR]['objects']:
    for rv_ref in rv_ref_list:
      ts = time.time()
      rvs = []
      e_rvs = []
      # determine barycentric RV
      with open('RVS_'+obj+'_'+rv_ref+'.txt', 'w+') as file:

        for order in np.arange(2, 31, 1):
          try:
            iraf.wspectext(input=obj+'_vh_norm.fits[*,'+str(order)+',1]', output='template.txt', header='no')
            template = np.loadtxt('template.txt')
            
            wav = template[:, 0]        
            mid1 = int(len(wav)/4.0)
            mid2 = int(len(wav)/12.0)

            wav1 = wav[mid1*1]
            wav2 = wav[mid1*3]
            #print 'WAV:', wav, mid1, wav1,wav2

            # create multiple wav1, wav2 ranges if requested
            if multi_sample:
              wav_ranges = [[wav[mid2*2], wav[mid2*7]], [wav[mid2*4], wav[mid2*8]], [wav[mid2*5], wav[mid2*10]]]
            else:
              wav_ranges = [[wav1, wav2]]

            for wav_r in wav_ranges:

              iraf.unlearn('fxcor')
              try:
                os.remove('test.txt')
                os.remove('test.log')
                os.remove('test.gki')
              except: pass
              try:
                os.remove('template.txt')
              except: pass
            
              iraf.fxcor(obj+'.ec.vh', '../../'+rv_ref, apertures=order, pixcorr='no', continuum='both', 
                         osample='a %s %s' % (wav_r[0], wav_r[1]), rsample='a %s %s' % (wav_r[0], wav_r[1]), 
                         interactive='no', output='test', function = 'gaussian',
                         apodize = 0.0, rebin = "smallest")
              fields = np.loadtxt('test.txt', comments='#', dtype='str')
              # print wav_r
              # print fields
              vrel, verr = fields[-3], fields[-1]
              # vrel, verr = fields[-3], 0.

              if vrel == 'INDEF' or np.float(verr) > 10. or np.abs(np.float(vrel)) > 350.:
                vrel = np.nan
                verr = np.nan

              print order, vrel, verr              
              file.write("%s %s %s\n" % (order, vrel, verr))
          
              rvs.append(float(vrel))
              e_rvs.append(float(verr))

          except:
            print_exception()
      
      # export all values in a list
      # file.write('\n')
      # file.write("'rv_echelle':[" + ','.join([str(v) for v in rvs]) + "]\n")
      # file.write("'e_rv_echelle':[" + ','.join([str(v) for v in e_rvs]) + "]\n")
        
      print '%s: RV median: %s, RV mean: %s, RV sigma: %s' % (obj, np.nanmedian(rvs), np.nanmean(rvs), np.nanstd(rvs))
      print 'RV fit time: '+str((time.time()-ts)/60.)+' min'
      print ''


def export_spectrum_to_txt(WORK_DIR, suffix='_vh_norm', export_cal=False):
  if export_cal:
    objects = observations[WORK_DIR]['calibs']
  else:
    objects = observations[WORK_DIR]['objects']
  for obj in objects:
    print ' Exporting file:', obj+suffix
    for order in np.arange(1, 31, 1):
      try:
        iraf.wspectext(input=obj+suffix+'.fits[*,'+str(order)+',1]', output=obj+suffix+'_order{:02.0f}.txt'.format(order), header='no')
      except:
        pass


def new_flatfield_lens_problems(WORK_DIR):
  
  def mini(params, x, y, yerr):    
    model = params['a'] + params['b']*x + params['c']*(x-params['d'])**2
    LL = (y - model)**2 / (yerr**2)
    return LL

  
  
  fig = plt.figure()

  for i, obj in enumerate(observations[WORK_DIR]['objects']):
    print obj
    if obj == 'EC59466':
      continue
    data = np.loadtxt('test_april_RVS_%s.txt' % obj)    
    mask = (data[:, 1] < 1000) & (data[:, 1] > -1000) & (data[:, 0] > 3)
    x = data[:, 0][mask]
    y = data[:, 1][mask] - np.median(data[:, 1][mask])
    yerr = data[:, 2][mask]
    
    params = Parameters()  
    params.add('a', 0)
    params.add('b', 0, False)
    params.add('c', 0, False)
    params.add('d', 0, False)
    results = minimize(mini, params, args = (x, y, yerr))
    params = results.params
    params.pretty_print()
    #report_fit(params)
    pars = [h.value for h in params.values()]
    
    medi = np.median(data[:, 1][mask][:14])
    y = data[:, 1][mask] - medi
    #plt.plot(x, i*10+pars[0]+pars[1]*x+pars[2]*(x-pars[3])**2)
    plt.axhline(-i*15, c='black', lw=0.5)
    plt.errorbar(x, y-i*15, yerr=yerr, fmt='x', ms=4, label=obj + ('  %.1f' % medi))

  plt.legend(bbox_to_anchor=(1.0, 1.0))
  plt.xlabel('aperture')
  plt.title('April spectra')
  plt.show()
  fig.savefig("RVS.pdf")

def print_exception():
  e = sys.exc_info()
  exc_type, exc_value, exc_traceback = e
  a, j = (traceback.extract_tb(exc_traceback, 1))[0][0:2]
  k = (traceback.format_exception_only(exc_type, exc_value))[0]
  print a, j, k

if __name__ == "__main__":
  main()

