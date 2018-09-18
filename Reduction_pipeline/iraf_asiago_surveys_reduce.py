# -*- encoding: utf-8 -*-
#!/usr/bin/env python
import sys
import os
from glob import glob
import pyfits
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')  # removes conflicts with pyraf and tkinter
import matplotlib.pyplot as plt
import operator
import time
import traceback
from pyraf import iraf
from datetime import date
from astropy.table import Table
from lmfit import minimize, Parameters, report_fit
from rv_corr_custom import correlate_spectra, correlate_order, get_RV_ref_spectrum, get_julia_dates_header

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
iraf.ccdred.instrum='blank.txt'
os.environ['PYRAF_BETA_STATUS'] = '1'

surveys = {
# -----------------------------------------
# -----------------------------------------
# ----------- TEST ------------------------
# -----------------------------------------
# -----------------------------------------
'WAV_REF_GEN':[
{   
	'DATE_DIR':'201801',
	'ID_field': 'source_id',
	'biases':['EC60718', 'EC60719', 'EC60720', 'EC60721', 'EC60722'], 
        'flats':[], 
	'objects':
		[
		{
			'spectra': ['EC60707'], 
			'calibs':  ['EC60708'],
			'ID': -1,
			'obj_name': 'N/A',
		},
		],
	'REF_ARC': None, 
	'REF_AP': None, 
	'ap_position':'center',
        'combine_exp': False
}
],
# -----------------------------------------
# -----------------------------------------
# ----------- BRIGHT STARS ----------------
# -----------------------------------------
# -----------------------------------------
'BRIGHT_STARS':[
{   

}
],
# -----------------------------------------
# -----------------------------------------
# ----------- CEMP ------------------------
# -----------------------------------------
# -----------------------------------------
'CEMP':[
{   
	'DATE_DIR':'201712',
	'ID_field': 'sobject_id',
	'biases':['EC60040', 'EC60041', 'EC60042', 'EC60043', 'EC60044'], 
	'flats':['EC60046', 'EC60047', 'EC60048'], 
	'objects':
		[
		{
			'spectra': ['EC60089'], 
			'calibs':  ['EC60090'], 
			'ID': 150409005101291,
			'obj_name': '2MASS J11333341-0043060'
		}		
		],
	'REF_ARC': None,  # 'ecEC60975.ec', 
	'REF_AP': '', 
	'ap_position':'left',
        'combine_exp': True
}
],
# -----------------------------------------
# -----------------------------------------
# ----------- NGC 6940 --------------------
# -----------------------------------------
# -----------------------------------------
'NGC6940':[
{   
	'DATE_DIR':'201807',
	'ID_field': 'source_id',
	'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
	'flats':['EC61160', 'EC61161', 'EC61162'], 
	'objects':
		[
		{
			'spectra': ['EC61140', 'EC61142'], 
			'calibs':  ['EC61141', 'EC61143'], 
			'ID': 0,
			'obj_name': 'HD 341007',
		},
		{
			'spectra': ['EC61144'], 
			'calibs':  ['EC61145'], 
			'ID': 0,
			'obj_name': 'HD 334233',
		},
		{
			'spectra': ['EC61147'], 
			'calibs':  ['EC61148'], 
			'ID': 0,
			'obj_name': 'HD 193664',
		},
		{
			'spectra': ['EC61149'], 
			'calibs':  ['EC61150'], 
			'ID': 0,
			'obj_name': 'HD 334320',
		},
		{
			'spectra': ['EC61151'], 
			'calibs':  ['EC61152'], 
			'ID': 0,
			'obj_name': 'HD 340749',
		},
		{
			'spectra': ['EC61153'], 
			'calibs':  ['EC61154'], 
			'ID': 0,
			'obj_name': 'TYC 2160-224-1',
		}
		],
	'REF_ARC':'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201808',
	'ID_field': 'source_id',
	'biases':['EC61277','EC61278','EC61279','EC61280','EC61281'], 
	'flats':['EC61289','EC61290','EC61291','EC61292','EC61293'], 
	'objects':
		[
		{
			'spectra': ['EC61199'], 
			'calibs':  ['EC61200'], 
			'ID': 0,
			'obj_name': 'HD 339984',
		},
		{
			'spectra': ['EC61201'], 
			'calibs':  ['EC61202'], 
			'ID': 0,
			'obj_name': 'HD 340848',
		},
		{
			'spectra': ['EC61203'], 
			'calibs':  ['EC61204'], 
			'ID': 0,
			'obj_name': 'TYC 3686-1067-1',
		},
		{
			'spectra': ['EC61224','EC61226'], 
			'calibs':  ['EC61225','EC61225'], 
			'ID': 0,
			'obj_name': 'TYC 2182-704-1',
		},
		{
			'spectra': ['EC61227','EC61229'], 
			'calibs':  ['EC61228','EC61228'], 
			'ID': 0,
			'obj_name': 'HD 335192',
		},
		{
			'spectra': ['EC61232','EC61234'], 
			'calibs':  ['EC61233','EC61233'], 
			'ID': 0,
			'obj_name': 'TYC 2164-786-1',
		},
		{
			'spectra': ['EC61235','EC61237'], 
			'calibs':  ['EC61236','EC61236'], 
			'ID': 0,
			'obj_name': 'HD 334808',
		},
		{
			'spectra': ['EC61253','EC61255'], 
			'calibs':  ['EC61254','EC61254'], 
			'ID': 0,
			'obj_name': 'HD 341012',
		},
		{
			'spectra': ['EC61256','EC61258'], 
			'calibs':  ['EC61257','EC61257'], 
			'ID': 0,
			'obj_name': 'HD 340994',
		},
		{
			'spectra': ['EC61259','EC61261'], 
			'calibs':  ['EC61260','EC61260'], 
			'ID': 0,
			'obj_name': 'HD 340386',
		},
		{
			'spectra': ['EC61262','EC61264'], 
			'calibs':  ['EC61263','EC61263'], 
			'ID': 0,
			'obj_name': 'TYC 2168-953-1',
		},
		{
			'spectra': ['EC61265','EC61267'], 
			'calibs':  ['EC61266','EC61266'], 
			'ID': 0,
			'obj_name': 'HD 334663',
		},
		{
			'spectra': ['EC61268','EC61270'], 
			'calibs':  ['EC61269','EC61269'], 
			'ID': 0,
			'obj_name': 'TYC 2163-844-1',
		},
		{
			'spectra': ['EC61287','EC61294'], 
			'calibs':  ['EC61288','EC61295'], 
			'ID': 0,
			'obj_name': 'HD 335040',
		}
		],
	'REF_ARC':'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},

],
# -----------------------------------------
# -----------------------------------------
# ----------- TRIPLE TWINS ----------------
# -----------------------------------------
# -----------------------------------------
'TRIPLE':[
{   
	'DATE_DIR':'201806',
	'ID_field': 'sobject_id',
	'biases':['EC60944', 'EC60945', 'EC60946', 'EC60947', 'EC60948'], 
	'flats':['EC60949', 'EC60950', 'EC60951', 'EC60952', 'EC60953'], 
	'objects':
		[
		{
			'spectra': ['EC60984'], 
			'calibs':  ['EC60985'], 
			'ID': 160531006101153,
			'obj_name': 'TYC 489-2258-1',
		}		
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201807',
	'ID_field': 'sobject_id',
	'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
	'flats':['EC61160', 'EC61161', 'EC61162'], 
	'objects':
		[
		{
			'spectra': ['EC61131', 'EC61133', 'EC61135'], 
			'calibs':  ['EC61132', 'EC61134', 'EC61134'], 
			'ID': 160531006101153,
			'obj_name': 'TYC 489-2258-1',
		},
		{
			'spectra': ['EC61136'], 
			'calibs':  ['EC61137'], 
			'ID': 160531005101256,
			'obj_name': 'TYC 5805-1373-1',
		},
		{
			'spectra': ['EC61138'], 
			'calibs':  ['EC61139'], 
			'ID': 170912002901303,
			'obj_name': '2MASS J23520336+1617418',
		}
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},

{   
	'DATE_DIR':'201808',
	'ID_field': 'sobject_id',
	'biases':['EC61277','EC61278','EC61279','EC61280','EC61281'], 
	'flats':['EC61289','EC61290','EC61291','EC61292','EC61293'], 
	'objects':
		[
		{
			'spectra': ['EC61196','EC61198'], 
			'calibs':  ['EC61197','EC61197'], 
			'ID': 160531005101256,
			'obj_name': 'TYC 5805-1373-1',
		},
		{
			'spectra': ['EC61210'], 
			'calibs':  ['EC61211'], 
			'ID': 170121002801292,
			'obj_name': '2MASS J05174448+2022042',
		},
		{
			'spectra': ['EC61230'], 
			'calibs':  ['EC61231'], 
			'ID': 150703005601062,
			'obj_name': '2MASS J22212296-0937544',
		}
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
],
# -----------------------------------------
# -----------------------------------------
# ----------- ORION CLUSTERS FOR JANEZ ----
# -----------------------------------------
# -----------------------------------------
'ORION_CLUSTERS':[
{   
	'DATE_DIR':'201808',
	'ID_field': 'source_id',
	'biases':['EC61277','EC61278','EC61279','EC61280','EC61281'], 
	'flats':['EC61289','EC61290','EC61291','EC61292','EC61293'], 
	'objects':
		[
		{
			'spectra': ['EC61240'], 
			'calibs':  ['EC61241'], 
			'ID': 0,
			'obj_name': 'TYC 105-1958-1',
		},
                {
			'spectra': ['EC61242'], 
			'calibs':  ['EC61243'], 
			'ID': 0,
			'obj_name': 'TYC 105-15-1',
		},
                {
			'spectra': ['EC61244'], 
			'calibs':  ['EC61245'], 
			'ID': 0,
			'obj_name': 'HD 287842',
		},
                {
			'spectra': ['EC61273'], 
			'calibs':  ['EC61274'], 
			'ID': 0,
			'obj_name': 'HD 287845',
		},
                {
			'spectra': ['EC61275'], 
			'calibs':  ['EC61276'], 
			'ID': 0,
			'obj_name': 'HD 290380',
		}		
		],
	'REF_ARC':'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
}
]
,
# -----------------------------------------
# -----------------------------------------
# ----------- GAIA RV STANDARDS -----------
# -----------------------------------------
# -----------------------------------------
'GAIA_RV_STAND':[
{   
	'DATE_DIR':'201806',
	'ID_field': 'source_id',
	'biases':['EC60944', 'EC60945', 'EC60946', 'EC60947', 'EC60948'], 
        'flats':['EC60949', 'EC60950', 'EC60951', 'EC60952', 'EC60953'], 
	'objects':
		[
		{
			'spectra': ['EC60966'], 
			'calibs':  ['EC60967'], 
			'ID': 0,
			'obj_name': 'HIP110341',
		},
                {
			'spectra': ['EC60968'], 
			'calibs':  ['EC60969'], 
			'ID': 0,
			'obj_name': 'HIP100017',
		},
                {
			'spectra': ['EC60970'], 
			'calibs':  ['EC60971'], 
			'ID': 0,
			'obj_name': 'HIP078424',
		},
                {
			'spectra': ['EC60972'], 
			'calibs':  ['EC60973'], 
			'ID': 0,
			'obj_name': 'HIP085268',
		},
                {
			'spectra': ['EC60974'], 
			'calibs':  ['EC60975'], 
			'ID': 0,
			'obj_name': 'HIP083389',
		},
                {
			'spectra': ['EC60976'], 
			'calibs':  ['EC60977'], 
			'ID': 0,
			'obj_name': 'HIP085653',
		}		
		],
	'REF_ARC':None,  # 'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{   
	'DATE_DIR':'201807',
	'ID_field': 'source_id',
	'biases':['EC61155', 'EC61156', 'EC61157', 'EC61158', 'EC61159'], 
	'flats':['EC61160', 'EC61161', 'EC61162'], 
	'objects':
		[
		{
			'spectra': ['EC61147'], 
			'calibs':  ['EC61148'], 
			'ID': 0,
			'obj_name': 'HIP100017',
		}		
		],
	'REF_ARC':None,  # 'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
}
]

}

AP_LIM = {'left': [-20, 10, '-19:-11'], 'center': [-15, 15, '-15:-10,10:15']}

root_dir_data = '/home/nandir/gigli/home/janez/ftp/asiago_observations/'
root_dir_arcs = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files_manual/'
root_dir_rvref = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_rv_files/'

process_survey = 'TRIPLE'
REWRITE = False  # cheks if result is already present and skips processing of that stage accordingly
REANALYSE = True  # does RV and parameter determination have to be run again if results already exist

def main():

  work_path = 'observations/' + process_survey
    
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
  
  # get objects that have to be reduduced and/or analysed
  survey_obs = surveys[process_survey]

  # prepare Table object that will hold processing results
  res_tab_out_file = process_survey+'_reduced.fits'
  if os.path.isfile(res_tab_out_file):
    res_tab = Table.read(res_tab_out_file)
  else:
    res_tab = Table(names=['Asiago_id', 'dir', 'obj_name', survey_obs[0]['ID_field'], 'MJD', 'JD', 'rv', 'e_rv', 'teff', 'e_teff', 'feh', 'e_feh', 'logg', 'e_logg'], 
                    dtype=['S35', 'S6', 'S35', 'int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64'])

  # Process individually every month/date/day/etc
  for i_d in range(len(survey_obs)):
     survey_cur_date = survey_obs[i_d]
     
     root_dir_cur_date = root_dir_data + survey_cur_date['DATE_DIR'] + '/'
     # preapre data that are needed for reduction of every spectra in selected period
     mbias = create_masterbias(survey_cur_date, root_dir_cur_date, rewrite=REWRITE)    
     mflat = create_masterflat(survey_cur_date, root_dir_cur_date, rewrite=REWRITE)    
     # mflat_norm = normalize_masterflat(survey_cur_date, rewrite=REWRITE)  

     # Process with actual spectra reduction
     # Individually reduce every object observed on this date
     for i_o in range(len(survey_cur_date['objects'])):
       obj_spectra_data = survey_cur_date['objects'][i_o]
       folder_cleaned = clean_folder_write_header(obj_spectra_data, root_dir_cur_date, rewrite=REWRITE)
       if folder_cleaned:
         bias_correct(obj_spectra_data, mbias)  
         # flat_correct(obj_spectra_data, mflat_norm)  
         # remove_cosmics(obj_spectra_data) 

       # spectrum extraction
       AP_LIM_use = AP_LIM[survey_cur_date['ap_position']]
       do_apall_for_real(obj_spectra_data, AP_LIM_use, rewrite=REWRITE)  
       apply_apall(obj_spectra_data, AP_LIM_use, rewrite=REWRITE)  

       # determine wavelength solution
       if process_survey == 'WAV_REF_GEN':
         wavelength_solution_ref_create(obj_spectra_data)
         continue
 
       wavelength_solution(obj_spectra_data, root_dir_arcs, survey_cur_date['REF_ARC']) 
       apply_wav_solution(obj_spectra_data, include_cal=False)
       vhelio_correct(obj_spectra_data)   

       # check_spectra(obj_spectra_data, '.ec.vh') 
       final_norm_orders, final_norm_1d = combine_normalize_images(obj_spectra_data, combine=survey_cur_date['combine_exp'])

       # TODO: check if spectra was already analyzed
       #rvref_spec_list = glob(root_dir_rvref + 'T*K2SNWNVR20N.fits')
       #rvref_spec_list = [rf.split('/')[-1].split('.')[0] for rf in rvref_spec_list]
       rvref_spec_list = ['solar']
       #rvref_spec_list = ['EC60966','EC60968','EC60970','EC60972','EC60974','EC60976','EC61147']

       rv_med, rv_std = get_RV_custom_corr_perorder(final_norm_orders, root_dir_rvref, rvref_spec_list, plot_rv=True)
       print '    Final RV values:', rv_med, rv_std

       # add results to table
       for i_f, final_img in enumerate(final_norm_orders):
         idx_row_exist = np.where(res_tab['Asiago_id'] == final_img)[0]
         # first remove old results
         if len(idx_row_exist) > 0:
           res_tab.remove_rows(idx_row_exist)
         # add new reults
         mjd, jd = get_julia_dates_header(final_img)
         res_list = [final_img, survey_cur_date['DATE_DIR'], obj_spectra_data['obj_name'], obj_spectra_data['ID'],
                     mjd, jd, rv_med[i_f], rv_std[i_f], 0., 0., 0., 0., 0., 0.]
         print '    Final results:', res_list
         res_tab.add_row(res_list)
         res_tab.write(res_tab_out_file, overwrite=True)

def create_masterbias(obs_data, root_data, rewrite=True):
  print "\n + Creating masterbias "+obs_data['DATE_DIR']+"\n"
  mb_filename = 'masterbias_'+obs_data['DATE_DIR']+'.fits'

  if not rewrite:
    if os.path.isfile(mb_filename):
      print '  -- Masterbias already created.'
      return mb_filename

  try: os.remove(mb_filename)
  except: pass

  iraf.zerocombine(input=','.join([root_data+i for i in obs_data['biases']]), output=mb_filename, combine='median', reject='none', ccdtype='', process='no', Stdout="/dev/null")

  return mb_filename

def create_masterflat(obs_data, root_data, rewrite=True):
  print "\n + Creating masterflat "+obs_data['DATE_DIR']+"\n"

  mb_filename = 'masterbias_'+obs_data['DATE_DIR']+'.fits'
  mf_filename = 'masterflat_'+obs_data['DATE_DIR']+'.fits'

  if not rewrite:
    if os.path.isfile(mb_filename):
      print '  -- Masterflat already created.'
      return mf_filename

  try: os.remove(mf_filename)
  except: pass

  # subtract masterbias from flats
  for flat in obs_data['flats']:
    try: os.remove(flat+'.fits')
    except: pass
    iraf.imarith(operand1=root_data+flat, op='-', operand2=mb_filename, result=flat, Stdout="/dev/null")
  
  iraf.zerocombine(input=','.join([flat for flat in obs_data['flats']]), output=mf_filename, combine='median', reject='none', ccdtype='', process='no')

  for flat in obs_data['flats']:
    os.remove(flat+'.fits')

  return mf_filename

def normalize_masterflat(obs_data, rewrite=True):
  print "\n + Normalizing masterflat "+obs_data['DATE_DIR']+"\n"

  mf_filename_orig = 'masterflat_'+obs_data['DATE_DIR']+'.fits'
  mf_filename_norm = 'masterflat_'+obs_data['DATE_DIR']+'_norm.fits'

  if not rewrite:
    if os.path.isfile(mb_filename):
      print '  -- Normalized masterflat already created.'
      return mf_filename_norm

  try: os.remove(mf_filename_norm)
  except: pass

  iraf.apflatten(input=mf_filename_orig, output=mf_filename_norm, referen=obs_data['objects'][0], interac='no', find='no', recenter='no', resize='no', edit='no', trace='no', fittrac='no', flatten='yes', fitspec='yes', functio='spline3', order=13, sample='51:2098', niterat=5)

  return mf_filename_norm

def check_spectra(spec_data, ext):
  print "\n + Visually check spectra before proceeding\n"

  for obj in spec_data['spectra']:
    iraf.splot(obj+ext)

def clean_folder_write_header(spec_data, root_data, rewrite=True):
  print '\n + Cleaning folder, truncating edges, fixing bad pixels\n'

  if not rewrite:
    # check if all needed files already exist
    n_s = len(spec_data['spectra'])
    n_miss = 0
    for i_s in range(n_s):
      for chk_file in [spec_data['spectra'][i_s]+'.fits', spec_data['spectra'][i_s]+'.ec.fits', spec_data['calibs'][i_s]+'.fits']:
        if not os.path.isfile(chk_file): 
          n_miss += 1
    if n_miss == 0:
      print '  -- Folder already prepared.'
      return False

  for obj in spec_data['spectra']:
    try: os.remove(obj+'.fits')
    except: pass
    try: os.remove(obj+'.ec.fits')
    except: pass
    iraf.imcopy(input=root_data+obj, output=obj) 

  for cal in spec_data['calibs']:
    try: os.remove(cal+'.fits')
    except: pass  
    iraf.imcopy(input=root_data+cal, output=cal)

  with open('badpix.txt', 'w+') as file:
    file.write('704 704 1262  2048\n703 705 1262  1277')  

  edit_fits = ','.join(spec_data['spectra']+spec_data['calibs'])
  iraf.hedit(images=edit_fits, fields='CCDSEC', value='[51:2098,1:2048]', add='yes', verify='no')
  iraf.hedit(images=edit_fits, fields='DATASEC', value='[51:2098,1:2048]', add='yes', verify='no') 

  # FIX badpix
  iraf.ccdproc(images=edit_fits, ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='no', fixfile='badpix.txt')

  return True

def bias_correct(obs_data, master_bias):
  print '\n + Correcting images and arcs for bias\n'

  iraf.ccdproc(images=','.join(obs_data['spectra']+obs_data['calibs']), ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='yes', darkcor='no', flatcor='no', zero=master_bias, fixfile='badpix.txt')

def flat_correct(obs_data, master_flat):
  print '\n + Correcting images and arcs for flat\n'

  iraf.ccdproc(images=','.join(obs_data['spectra']+obs_data['calibs']), ccdtype='', fixpix='yes', oversca='no', trim='yes', zerocor='no', darkcor='no', flatcor='yes', flat=master_flat, fixfile='badpix.txt')

def do_apall_for_real(obs_data, AP_LIM, ref_ap=None, rewrite=True):
  print '\n + Doing apall for real\n'

  if not rewrite:
    n_miss = 0
    for obj in obs_data['spectra']:
      if not os.path.isfile('database/ap'+obj):
        n_miss += 1
    if n_miss == 0:
      print '  -- Appal tracing already performed.'
      return False

  #do apall for real with background and all
  #This tracing will be used from now on. Pay attention to the background subtraction. You will have to verify each tracing on each reference image (first images in the log.txt).   
  try:
    shutil.copy(ref_ap, 'database/' + ref_ap.split('/')[-1])  
  except:
    pass

  for obj in obs_data['spectra']:
    if ref_ap is not None:      
      ref = observations[WORK_DIR]['REF_AP'].split('/')[-1].replace('ap', '')
    else:
      ref = obj
    print 'APALL',obj, ref
    iraf.apall(input=obj, referen=ref, format='echelle', interac='yes', find='yes', recente='yes', resize='yes', edit='yes', trace='yes', fittrac='yes', extract='no', extras='no', review='no', line=600, nsum=30, lower=-6, upper=6, width=11, radius=11, thresho=0.15, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro='median', b_order=1, b_sampl=AP_LIM[2], ylevel='INDEF', llimit=-5, ulimit=5, weights='none', clean='yes', order='increasing')#, apidtable='apert.txt')

def apply_apall(obs_data, AP_LIM, rewrite=True):
  print '\n + Applying apall to rest of the objects and arcs\n'

  if not rewrite:
    n_s = len(obs_data['spectra'])
    n_miss = 0
    for i_s in range(n_s):
      for chk_file in [obs_data['spectra'][i_s]+'.ec.fits', obs_data['spectra'][i_s]+'.ec.nosky.fits', obs_data['calibs'][i_s]+'.ec.fits']:
        if not os.path.isfile(chk_file): 
          n_miss += 1
    if n_miss == 0:
      print '  -- Appal already applied to all spectra and calib.'
      return False

  for obj in obs_data['spectra']:
    print ' - Processing image %s \n' % obj

    try: 
      os.remove(obj+'.ec.fits')
      os.remove(obj+'.ec.nosky.fits')
    except: pass 
    
    for j in [['.ec', 'median'], ['.ec.nosky', 'none']]:
      iraf.apall(input=obj, output=obj+j[0], format='echelle', referen=obj, interac='no', find='yes', recente='no', resize='no', edit='yes', trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro=j[1], b_order=1, b_sampl=AP_LIM[2], ylevel='INDEF', llimit=-5, ulimit=5, weights='none', clean='yes', order='increasing')#, apidtable='apert.txt')  

  for i, cal in enumerate(obs_data['calibs']):
    print ' - Processing arc %s \n' % cal

    try: os.remove(cal+'.ec.fits')
    except: pass 
    print '   input:',cal,'    referen:', obs_data['spectra'][i]
    iraf.apall(input=cal, format='echelle', referen=obs_data['spectra'][i], interac='no', find='no', recente='no', resize='no', edit='no',trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro='none', ylevel='INDEF', llimit=-5, ulimit=5, order='increasing')#, apidtable='apert.txt')    

def wavelength_solution_ref_create(obs_data, rewrite=True):
  print 'Creating new reference arc wavelength solution'
  iraf.unlearn('ecidentify')
  iraf.unlearn('ecreidentify')
  for cal in obs_data['calibs']:
    # lower number of maxfeat for less work while removing bad emissions in arc, incresed minsep
    iraf.ecident(images=cal+'.ec', coordli='linelists$thar.dat', match=1, maxfeat=1250, ftype='emission', fwidth=4, cradius=5, thresho=2, minsep=3, functio='chebyshev', xorder=6, yorder=6, niterat=5, lowreje=3, highreje=3, autowri='yes')


def wavelength_solution(obs_data, arc_full_path, arc_file, rewrite=True):
  print '\n + Finding the wavelength solution\n'

  # calibrate wavelength
  # This step will require you to calibrate first arc by hand. The rest will be done automatically.
  # after l (automatic find other lines) change x and y order to 4
  # x - choose plots of the fit (per order, pixel, ...)

  if arc_file is not None:
    # ref arc was given
    arc_full_path_list = list([arc_full_path + arc_file])
  else:
    # ref arc was not given
    arc_full_path_list = glob(arc_full_path+'ec*.ec')

  for cal in obs_data['calibs']:
    iraf.unlearn('ecidentify')
    print ' Searching for the best arc' 
    arc_rms_all = list([])
    for arc_full_path in arc_full_path_list:
      # print arc_full_path
      if arc_full_path:
        try:
          # print 'shutil.copy'
          shutil.copy(arc_full_path, 'database/' + arc_full_path.split('/')[-1])  
        except:
          pass
      
        # print 'iraf.ecreident'
        iraf.unlearn('ecreidentify')
        log_arc = 'log_arc.txt'
        try: os.remove(log_arc)
        except: pass
        iraf.ecreident(images=cal+'.ec', referenc=arc_full_path.replace('/ec', '/').split('/')[-1], refit='yes', shift='INDEF', cradius=5, thresho=10, logfiles=log_arc)
        # parse output log gile to get RMS reading
        with open(log_arc) as log_h:
          log_lines = log_h.read()
        arc_rms = log_lines.split('\n')[-2].split(' ')[-1]
        arc_rms_all.append(arc_rms)
        #print arc_rms

    # select the best matching arc with the lowest RMS
    iraf.unlearn('ecreidentify')
    idx_arc = np.nanargmin(arc_rms_all)
    arc_full_use = arc_full_path_list[idx_arc]
    print ' All arc RMS:', arc_rms_all
    print ' Best:', arc_full_use.split('/')[-1], arc_rms_all[idx_arc]
    iraf.ecreident(images=cal+'.ec', referenc=arc_full_use.replace('/ec', '/').split('/')[-1], refit='yes', shift='INDEF', cradius=5, thresho=10)
    
    print 'iraf.ecident'
    # lower number of maxfeat for less work while removing bad emissions in arc, incresed minsep
    iraf.ecident(images=cal+'.ec', coordli='linelists$thar.dat', match=1, maxfeat=1250, ftype='emission', fwidth=4, cradius=5, thresho=2, minsep=3, functio='chebyshev', xorder=6, yorder=6, niterat=5, lowreje=3, highreje=3, autowri='yes')

def apply_wav_solution(obs_data, include_cal=False, rewrite=True):
  print '\n + Applying the wavelength solution\n'
  
  iraf.unlearn('dispcor')

  for i, obj in enumerate(obs_data['spectra']):
    cal = obs_data['calibs'][i]+'.ec'
    # export wav calibrated spectra
    print 'Using cal:', cal
    for j in ['.ec', '.ec.nosky']:
      iraf.refspectra(input=obj+j, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=obj+j, output=obj+j, lineari='no', verbose='yes')
    if include_cal:
      os.remove(cal+'_wvl.fits')
      # export wav calibrated cal image
      iraf.refspectra(input=cal, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=cal, output=cal+'_wvl', lineari='no', verbose='yes')

def vhelio_correct(obs_data):
  print '\n + VHELIO correction\n'
    
  for obj in obs_data['spectra']:
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

def combine_normalize_images(obs_data, combine=False):
  
  if combine:
    print '\n + Combine images\n'
    comb_sum_file = obs_data['spectra'][0]+'_vh_comb'
    comb_sum_norm_file = obs_data['spectra'][0]+'_vh_norm_comb'
    comb_sum_cont_file = obs_data['spectra'][0]+'_cont_comb'
    for sum_file in [comb_sum_file, comb_sum_norm_file, comb_sum_cont_file]:
      try: os.remove(sum_file+'.fits')       
      except: pass
 
    iraf.scombine(input=','.join([obj+'.ec.vh' for obj in obs_data['spectra']]), output=comb_sum_file, group='apertures', combine='sum', reject='none', Stdout="/dev/null")

  # normalize spectra
  # This step will require you to manually normalize all the spectra.
  print '\n + Normalize spectra\n'

  if combine:
    try: 
        os.remove('combined_sum_cone_rv_echellet.fits')       
    except: pass 
    iraf.continuum(input=comb_sum_file, output=comb_sum_cont_file, type='fit', replace='no', listonly='no', functio='cheb', order=13, low_rej=2, high_rej=3, naverag=-3, niter=9, interac='no', markrej='no')
    iraf.sarith(input1=comb_sum_file, op='/', input2=comb_sum_cont_file, output=comb_sum_norm_file, format='multispec', Stdout="/dev/null")
  else: 
    final_norm_files = list([]) 
    for obj in obs_data['spectra']:
      try: 
        os.remove(obj+'_cont.fits')
        os.remove(obj+'_vh_norm.fits')
      except: pass 
      iraf.continuum(input=obj+'.ec.vh', output=obj+'_cont', type='fit', replace='no', listonly='no', functio='cheb', order=13, low_rej=2, high_rej=3, naverag=-3, niter=9, interac='no', markrej='no', ask='yes')
      iraf.sarith(input1=obj+'.ec.vh', op='/', input2=obj+'_cont', output=obj+'_vh_norm', format='multispec', Stdout="/dev/null")
      final_norm_files.append(obj+'_vh_norm')

  #combine apertures
  print '\n + Combine apertures\n'

  if combine:
    try: 
        data_1d_file = obs_data['spectra'][0]+'_data1D_comb'
        cont_1d_file = obs_data['spectra'][0]+'_cont1D_comb'
	final_1d_file = obs_data['spectra'][0]+'_1D_vh_comb_norm'
        os.remove(data_1d_file+'.fits')       
        os.remove(cont_1d_file+'.fits')    
        os.remove(final_1d_file+'.0001.fits')
        os.remove(final_1d_file+'.fits')       
    except: pass
    iraf.scombine(input=comb_sum_file, output=data_1d_file, group='all', combine='sum', reject='none', Stdout="/dev/null")
    iraf.scombine(input=comb_sum_cont_file, output=cont_1d_file, group='all', combine='sum', reject='none', Stdout="/dev/null")

    iraf.sarith(input1=data_1d_file, op='/', input2=cont_1d_file, output=final_1d_file, format='onedspec', Stdout="/dev/null")

    # return name of the final file with combined and normalized spectra
    return [comb_sum_norm_file], [final_1d_file+'.0001']
  else:  
    final_1d_files = list([])
    for obj in obs_data['spectra']:
      try: 
        os.remove(obj+'_data1D.fits')
        os.remove(obj+'_cont1D.fits')
        os.remove(obj+'_1D_vh_norm.0001.fits')   
        os.remove(obj+'_1D_vh_norm.fits')        
      except: pass 
      iraf.scombine(input=obj+'.ec.vh', output=obj+'_data1D', group='all', combine='sum', reject='none', Stdout="/dev/null")
      iraf.scombine(input=obj+'_cont', output=obj+'_cont1D', group='all', combine='sum', reject='none', Stdout="/dev/null")

      iraf.sarith(input1=obj+'_data1D', op='/', input2=obj+'_cont1D', output=obj+'_1D_vh_norm', format='onedspec', Stdout="/dev/null")

      final_1d_files.append(obj+'_1D_vh_norm.0001')
    # return name of the final file with combined and normalized spectra
    return final_norm_files, final_1d_files

def export_observation_to_txt(obs_data, suffix='_vh_norm', export_cal=False):
  if export_cal:
    objects = obs_data['calibs']
  else:
    objects = obs_data['spectra']
  for obj in objects:
    print ' Exporting file:', obj+suffix
    for order in np.arange(1, 31, 1):
      try:
        iraf.wspectext(input=obj+suffix+'.fits[*,'+str(order)+',1]', output=obj+suffix+'_order{:02.0f}.txt'.format(order), header='no')
      except:
        pass

def export_spectrum_order_to_txt(exp_file, order=1):
    try:
      out_txt_file = exp_file+'_order{:02.0f}.txt'.format(order)
      iraf.wspectext(input=exp_file+'.fits[*,'+str(order)+',1]', output=out_txt_file, header='no')
      return out_txt_file
    except:
      return None
      pass

def get_RV_custom_corr_1D(obs_spectra, plot_rv=False, n_subsections = 31):
  pass

def get_RV_custom_corr_perorder(obs_spectra, rvref_dir, rvref_list, plot_rv=False, n_min_ord=5, n_max_ord=25):

  # prepare reference file
  rv_obs_mean_final = list([])
  rv_obs_std_final = list([])
  for obs_spectrum in obs_spectra:
    
    rv_obs_mean = list([])
    rv_obs_std = list([])
    for use_revref in rvref_list:
      ref_flx, ref_wvl = get_RV_ref_spectrum(rvref_dir+use_revref+'.fits')

      print '\n + Determinig RV for input per order spectrum', obs_spectrum, 'and reference', use_revref, '\n'    
      rv_shifts = list([])
      for i_order in range(n_min_ord, n_max_ord+1):
        order_txt_file = export_spectrum_order_to_txt(obs_spectrum, order=i_order)
        # print order_txt_file
        if order_txt_file is not None:
          rv_shifts.append(correlate_order(order_txt_file, ref_flx, ref_wvl, plot=True))
        else:
          rv_shifts.append(np.nan)
        # delete temp txt order file
        try: os.remove(order_txt_file)
        except: pass

      # analyse shifts
      rv_shifts = np.array(rv_shifts)
      rv_median = np.nanmedian(rv_shifts)
      # remove gross outlyers before computing rv std value
      rv_shifts[np.abs(rv_shifts - rv_median) > 10] = np.nan
      n_fin_rv = np.sum(np.isfinite(rv_shifts))
      if np.sum(np.isfinite(rv_shifts)) < (n_max_ord-n_min_ord)/5:
        rv_std = np.nan
        rv_median = np.nan
      elif n_fin_rv >= 3:
        rv_std = np.nanstd(rv_shifts)
      else:
        rv_std = np.nan
      rv_obs_mean.append(rv_median)
      rv_obs_std.append(rv_std)
      print rv_median, rv_std, n_fin_rv

      if plot_rv and np.isfinite(rv_median):
        rv_x = n_min_ord+np.arange(len(rv_shifts))
        plt.scatter(rv_x, rv_shifts)
        plt.axhline(rv_median)
        plt.ylim(rv_median-4., rv_median+4.)
        plt.xlim(rv_x[0]-1., rv_x[-1]+1.)
        plt.xlabel('Echelle order')
        plt.ylabel('Radial velocity')
        plt.title('Spectrum: '+obs_spectrum+'     Reference: '+use_revref)
        plt.tight_layout()
        plt.savefig(obs_spectrum+'_'+use_revref+'_rv.png', dpi=300)
        plt.close()

    idx_best_rv = np.nanargmin(rv_obs_std)
    rv_obs_mean_final.append(rv_obs_mean[idx_best_rv])
    rv_obs_std_final.append(rv_obs_std[idx_best_rv])

  # return all derived velocities and its std values
  return rv_obs_mean_final, rv_obs_std_final

if __name__ == "__main__":
  main()

