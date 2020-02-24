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
# ----------- TEST and REFERENCE ----------
# -----------------------------------------
# -----------------------------------------
'WAV_REF_GEN':[
#'TEST':[
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC53435'], 
			'calibs':  ['EC53436'], 
			'ID': 4,
			'obj_name': 'HR 6056',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
],
# -----------------------------------------
# -----------------------------------------
# ----------- BRIGHT STARS ----------------
# -----------------------------------------
# -----------------------------------------
'BRIGHT_STARS':[
{   
# TODO: fill with observed stars - mostly (or only) in 201804 and in 201907
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
	'REF_ARC': None, 
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
			'ID': 1856255973096189440,
			'obj_name': 'HD 341007',
		},
		{
			'spectra': ['EC61144'], 
			'calibs':  ['EC61145'], 
			'ID': 1861778716930466176,
			'obj_name': 'HD 334233',
		},
		{
			'spectra': ['EC61147'], 
			'calibs':  ['EC61148'], 
			'ID': 2246581329639103872,
			'obj_name': 'HD 193664',
		},
		{
			'spectra': ['EC61149'], 
			'calibs':  ['EC61150'], 
			'ID': 1860281216451967488,
			'obj_name': 'HD 334320',
		},
		{
			'spectra': ['EC61151'], 
			'calibs':  ['EC61152'], 
			'ID': 1855559741716879616,
			'obj_name': 'HD 340749',
		},
		{
			'spectra': ['EC61153'], 
			'calibs':  ['EC61154'], 
			'ID': 1832696256419483648,
			'obj_name': 'TYC 2160-224-1',
		}
		],
	'REF_ARC':None, 
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
			'ID': 1836225035909645696,
			'obj_name': 'HD 339984',
		},
		{
			'spectra': ['EC61201'], 
			'calibs':  ['EC61202'], 
			'ID': 1855632275121248896,
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
			'ID': 1858912702421124480,
			'obj_name': 'TYC 2182-704-1',
		},
		{
			'spectra': ['EC61227','EC61229'], 
			'calibs':  ['EC61228','EC61228'], 
			'ID': 1858044122605656704,
			'obj_name': 'HD 335192',
		},
		{
			'spectra': ['EC61232','EC61234'], 
			'calibs':  ['EC61233','EC61233'], 
			'ID': 1836001903772424832,
			'obj_name': 'TYC 2164-786-1',
		},
		{
			'spectra': ['EC61235','EC61237'], 
			'calibs':  ['EC61236','EC61236'], 
			'ID': 1859219809771537024,
			'obj_name': 'HD 334808',
		},
		{
			'spectra': ['EC61253','EC61255'], 
			'calibs':  ['EC61254','EC61254'], 
			'ID': 1857739557878767232,
			'obj_name': 'HD 341012',
		},
		{
			'spectra': ['EC61256','EC61258'], 
			'calibs':  ['EC61257','EC61257'], 
			'ID': 1856514942447936384,
			'obj_name': 'HD 340994',
		},
		{
			'spectra': ['EC61259','EC61261'], 
			'calibs':  ['EC61260','EC61260'], 
			'ID': 1832568335110921984,
			'obj_name': 'HD 340386',
		},
		{
			'spectra': ['EC61262','EC61264'], 
			'calibs':  ['EC61263','EC61263'], 
			'ID': 1860142540544973056,
			'obj_name': 'TYC 2168-953-1',
		},
		{
			'spectra': ['EC61265','EC61267'], 
			'calibs':  ['EC61266','EC61266'], 
			'ID': 1862525560201870592,
			'obj_name': 'HD 334663',
		},
		{
			'spectra': ['EC61268','EC61270'], 
			'calibs':  ['EC61269','EC61269'], 
			'ID': 1836217442407561216,
			'obj_name': 'TYC 2163-844-1',
		},
		{
			'spectra': ['EC61287','EC61294'], 
			'calibs':  ['EC61288','EC61295'], 
			'ID': 1859284577875371264,
			'obj_name': 'HD 335040',
		}
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201809',
	'ID_field': 'source_id',
	'biases':['EC61355','EC61356','EC61357','EC61358','EC61359','EC61360','EC61361'], 
	'flats':['EC61362','EC61363','EC61364','EC61365'], 
	'objects':
		[
		{
			'spectra': ['EC61320', 'EC61322'], 
			'calibs':  ['EC61321', 'EC61321'], 
			'ID': 1857223062294236416,
			'obj_name': 'HD 340770',
		},
		{
			'spectra': ['EC61323', 'EC61325'], 
			'calibs':  ['EC61324', 'EC61324'], 
			'ID': 1857461553227085056,
			'obj_name': 'TYC 2169-953-1',
		},
		{
			'spectra': ['EC61326', 'EC61328'], 
			'calibs':  ['EC61327', 'EC61327'], 
			'ID': 1860010530433894784,
			'obj_name': 'HD 340182',
		},
		{
			'spectra': ['EC61380', 'EC61382'],  # poor conditions for both exposures
			'calibs':  ['EC61381', 'EC61383'], 
			'ID': 1857535701556727680,
			'obj_name': 'HD 334739',
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
{   
	'DATE_DIR':'201812',
	'ID_field': 'sobject_id',
	'biases':['EC62044','EC62045','EC62046','EC62047','EC62048','EC62049','EC62050'], 
	'flats':['EC62051','EC62052','EC62053','EC62054','EC62055'], 
	'objects':
		[
		{
			'spectra': ['EC62040','EC62042'], 
			'calibs':  ['EC62041','EC62043'], 
			'ID': 150411004101331,
			'obj_name': 'TYC 4923-896-1',
		},
                {
			'spectra': ['EC62069','EC62071'], 
			'calibs':  ['EC62070','EC62072'], 
			'ID': 150409002601317,
			'obj_name': '2MASS J11432902-0430241',
		},
                {
			'spectra': ['EC62107','EC62109'], 
			'calibs':  ['EC62108','EC62110'], 
			'ID': 160401003901215,
			'obj_name': 'TYC 5546-59-1',
		}			
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201901',
	'ID_field': 'sobject_id',
	'biases':['EC62399','EC62400','EC62401','EC62402','EC62403'], 
	'flats':['EC62396','EC62397','EC62398'], 
	'objects':
		[
		{
			'spectra': ['EC62293','EC62295'], 
			'calibs':  ['EC62294','EC62294'], 
			'ID': 151111002101170,
			'obj_name': 'TYC 20-218-1',
		},
                {
			'spectra': ['EC62301','EC62303'], 
			'calibs':  ['EC62302','EC62304'], 
			'ID': 170102001901154,
			'obj_name': 'TYC 1387-334-1',
		},
                {
			'spectra': ['EC62305','EC62307'], 
			'calibs':  ['EC62306','EC62308'], 
			'ID': 170412002401082,
			'obj_name': 'TYC 6014-1231-1',
		},
                {
			'spectra': ['EC62309','EC62311'], 
			'calibs':  ['EC62310','EC62312'], 
			'ID': 160130005201082,
			'obj_name': 'TYC 6056-1107-1',
		},
                {
			'spectra': ['EC62323'], 
			'calibs':  ['EC62324'], 
			'ID': 170508004801312,
			'obj_name': 'TYC 6186-1546-1',
		},
                {
			'spectra': ['EC62328','EC62330'], 
			'calibs':  ['EC62329','EC62331'], 
			'ID': 151111002101059,
			'obj_name': 'TYC 21-93-1',
		},
                {
			'spectra': ['EC62332','EC62334'], 
			'calibs':  ['EC62333','EC62335'], 
			'ID': 151111002101116,
			'obj_name': 'TYC 4682-1181-1',
		},
                {
			'spectra': ['EC62345','EC62347'], 
			'calibs':  ['EC62346','EC62348'], 
			'ID': 161228002501030,
			'obj_name': 'tri0816',
		},
                {
			'spectra': ['EC62359','EC62361'], 
			'calibs':  ['EC62360','EC62362'], 
			'ID': 150409005601390,
			'obj_name': 'tri1101',
		},
                {
			'spectra': ['EC62369'], 
			'calibs':  ['EC62370'], 
			'ID': 161006004401018,
			'obj_name': 'TYC 12-345-1',
		}		
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201902',
	'ID_field': 'sobject_id',
	'biases':['EC62462','EC62463','EC62464','EC62465','EC62466','EC62467','EC62468'], 
	'flats':['EC62469','EC62470','EC62471'], 
	'objects':
		[
		{
			'spectra': ['EC62445','EC62447'], 
			'calibs':  ['EC62446','EC62448'], 
			'ID': 160403003601392,
			'obj_name': 'TYC 5531-230-1',
		},
		{
			'spectra': ['EC62498','EC62500','EC62543'], 
			'calibs':  ['EC62499','EC62501','EC62544'], 
			'ID': 170416004801356,
			'obj_name': 'tri1501',
		},
                {
			'spectra': ['EC62514'],  # bad focus in the middle of some orders/traces, too low for more exposures 
			'calibs':  ['EC62515'], 
			'ID': 151111002101116,
			'obj_name': 'TYC 4682-1181-1',
		},
                {
			'spectra': ['EC62516'],  # low signal, focus ok
			'calibs':  ['EC62517'], 
			'ID': 171102004501327,
			'obj_name': 'tri0406',
		}		
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
},
{   
	'DATE_DIR':'201907',
	'ID_field': 'sobject_id',
	'biases':['EC63188', 'EC63189', 'EC63190', 'EC63191', 'EC63192'], 
	'flats':['EC63193', 'EC63194', 'EC63195'], 
	'objects':
		[
		{
			'spectra': ['EC63273'], 
			'calibs':  ['EC63274'], 
			'ID': 170806004701033,
			'obj_name': '2MASS J22002743-0129328',
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
			'ID': 3236079328631677824,
			'obj_name': 'TYC 105-1958-1',
		},
                {
			'spectra': ['EC61242'], 
			'calibs':  ['EC61243'], 
			'ID': 3222282278769474304,
			'obj_name': 'TYC 101-15-1',
		},
                {
			'spectra': ['EC61244'], 
			'calibs':  ['EC61245'], 
			'ID': 3222261594207038464,
			'obj_name': 'HD 287842',
		},
                {
			'spectra': ['EC61273'], 
			'calibs':  ['EC61274'], 
			'ID': 3222164901606950784,
			'obj_name': 'HD 287845',
		},
                {
			'spectra': ['EC61275'], 
			'calibs':  ['EC61276'], 
			'ID': 3220462655745525632,
			'obj_name': 'HD 290380',
		}		
		],
	'REF_ARC':None, 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{   
	'DATE_DIR':'201809',
	'ID_field': 'source_id',
	'biases':[' ','EC61356','EC61357','EC61358','EC61359','EC61360','EC61361'], 
	'flats':['EC61362','EC61363','EC61364','EC61365'], 
	'objects':
		[
		{
			'spectra': ['EC61340'], 
			'calibs':  ['EC61341'], 
			'ID': 3223817128282174976,
			'obj_name': 'HD 35911',
		},
		{
			'spectra': ['EC61388'], 
			'calibs':  ['EC61389'], 
			'ID': 3223831971689138688,
			'obj_name': 'HD 287888',
		},
		{
			'spectra': ['EC61390'], 
			'calibs':  ['EC61391'], 
			'ID': 3223708826387363072,
			'obj_name': 'HD 288061',
		},
		{
			'spectra': ['EC61392'], 
			'calibs':  ['EC61393'], 
			'ID': 3223540841626655360,
			'obj_name': 'HD 287974',
		},
		{
			'spectra': ['EC61394'], 
			'calibs':  ['EC61395'], 
			'ID': 3222286127060186368,
			'obj_name': 'TYC 105-146-1',
		}		
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{   
	'DATE_DIR':'201812',
	'ID_field': 'source_id',
	'biases':['EC62044','EC62045','EC62046','EC62047','EC62048','EC62049','EC62050'], 
	'flats':['EC62051','EC62052','EC62053','EC62054','EC62055'], 
	'objects':
		[
		{
			'spectra': ['EC62086'], 
			'calibs':  ['EC62087'], 
			'ID': 3220620813620479488,
			'obj_name': 'HD 290356',
		},
                {
			'spectra': ['EC62088'], 
			'calibs':  ['EC62089'], 
			'ID': 3221694418005734656,
			'obj_name': 'HD 290198',
		},
                {
			'spectra': ['EC62090'], 
			'calibs':  ['EC62091'], 
			'ID': 3222164905902845952,
			'obj_name': 'HD 287845',
		},
                {
			'spectra': ['EC62092'], 
			'calibs':  ['EC62093'], 
			'ID': 3222191221167193216,
			'obj_name': 'HD 287854',
		},
                {
			'spectra': ['EC62094'], 
			'calibs':  ['EC62095'], 
			'ID': 3222068419461962240,
			'obj_name': 'HD 287802',
		},
                {
			'spectra': ['EC62096'], 
			'calibs':  ['EC62097'], 
			'ID': 3223747854754812544,
			'obj_name': '2MASS J05292353+0207150',
		},
                {
			'spectra': ['EC62098'], 
			'calibs':  ['EC62099'], 
			'ID': 3222243344890328576,
			'obj_name': 'TYC 101-171-1',
		},
                {
			'spectra': ['EC62100'], 
			'calibs':  ['EC62101'], 
			'ID': 3222306193147280000,
			'obj_name': 'TYC 105-350-1',
		}		
		],
	'REF_ARC':None,  
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
                        'RV_ref': 8.627
		},
                {
			'spectra': ['EC60968'], 
			'calibs':  ['EC60969'], 
			'ID': 0,
			'obj_name': 'HIP100017',
                        'RV_ref': -4.455
		},
                {
			'spectra': ['EC60970'], 
			'calibs':  ['EC60971'], 
			'ID': 0,
			'obj_name': 'HIP078424',
                        'RV_ref': -21.748
		},
                {
			'spectra': ['EC60972'], 
			'calibs':  ['EC60973'], 
			'ID': 0,
			'obj_name': 'HIP085268',
                        'RV_ref': 13.495
		},
                {
			'spectra': ['EC60974'], 
			'calibs':  ['EC60975'], 
			'ID': 0,
			'obj_name': 'HIP083389',
                        'RV_ref': -46.883
		},
                {
			'spectra': ['EC60976'], 
			'calibs':  ['EC60977'], 
			'ID': 0,
			'obj_name': 'HIP085653',
                        'RV_ref': -83.98
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
                        'RV_ref': -4.455
		}		
		],
	'REF_ARC':None,  # 'ecEC60975.ec', 
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
}
],
# -----------------------------------------
# -----------------------------------------
# ----------- GAIA BENCHMARK --------------
# -----------------------------------------
# -----------------------------------------
'GAIA_BENCH':[
{   
	'DATE_DIR':'201809',
	'ID_field': 'source_id',
	'biases':['EC61355','EC61356','EC61357','EC61358','EC61359','EC61360','EC61361'], 
	'flats':['EC61362','EC61363','EC61364','EC61365'], 
	'objects':
		[
		{
			'spectra': ['EC61329','EC61331'], 
			'calibs':  ['EC61330','EC61330'], 
			'ID': 2661005953843811456,
			'obj_name': 'HD 220009',
		},		
		{
			'spectra': ['EC61332','EC61334'], 
			'calibs':  ['EC61333','EC61333'], 
			'ID': 3250489115708824064,
			'obj_name': 'HD  22879',
		},		
		{
			'spectra': ['EC61342','EC61344'], 
			'calibs':  ['EC61343','EC61343'], 
			'ID': 1823067317695767552,
			'obj_name': 'gam Sge',
		},		
		{
			'spectra': ['EC61368','EC61370'], #,'EC61366'], - saturated
			'calibs':  ['EC61367','EC61369'], 
			'ID': -1,
			'obj_name': 'alf Cet',
		},		
		{
			'spectra': ['EC61371'], #,'EC61373'], - saturated
			'calibs':  ['EC61372'], #,'EC61372'], - saturated
			'ID': -1,
			'obj_name': 'alf Tau',
		}		
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': True
}
],

# -----------------------------------------
# -----------------------------------------
# ----------- GREGOR TEST OBJECTS ---------
# -----------------------------------------
# -----------------------------------------
'GREGOR_TEST':[
{	
	'DATE_DIR':'201704',
	'ID_field': 'source_id',
	'biases':['EC59521','EC59522','EC59523','EC59524','EC59525'], 
	'flats':['EC59518','EC59519','EC59520'], 
	'objects':
		[
		{
			'spectra': ['EC59433'], 
			'calibs':  ['EC59434'], 
			'ID': 0,
			'obj_name': 'EC59433',
		},
		{
			'spectra': ['EC59435'], 
			'calibs':  ['EC59436'], 
			'ID': 0,
			'obj_name': 'EC59435',
		},
		{
			'spectra': ['EC59437'], 
			'calibs':  ['EC59438'], 
			'ID': 0,
			'obj_name': 'EC59437',
		},
		{
			'spectra': ['EC59440'], 
			'calibs':  ['EC59441'], 
			'ID': 0,
			'obj_name': 'EC59440',
		},
		{
			'spectra': ['EC59442'], 
			'calibs':  ['EC59443'], 
			'ID': 0,
			'obj_name': 'EC59442',
		},
		{
			'spectra': ['EC59444'], 
			'calibs':  ['EC59445'], 
			'ID': 0,
			'obj_name': 'EC59444',
		},
		{
			'spectra': ['EC59448'], 
			'calibs':  ['EC59449'], 
			'ID': 0,
			'obj_name': 'EC59448',
		},
		{
			'spectra': ['EC59462'], 
			'calibs':  ['EC59463'], 
			'ID': 0,
			'obj_name': 'EC59462',
		},
		{
			'spectra': ['EC59464'], 
			'calibs':  ['EC59465'], 
			'ID': 0,
			'obj_name': 'EC59464',
		},
		{
			'spectra': ['EC59466'], 
			'calibs':  ['EC59467'], 
			'ID': 0,
			'obj_name': 'EC59466',
		},
		{
			'spectra': ['EC59493'], 
			'calibs':  ['EC59494'], 
			'ID': 0,
			'obj_name': 'EC59493',
		},
		{
			'spectra': ['EC59495'], 
			'calibs':  ['EC59496'], 
			'ID': 0,
			'obj_name': 'EC59495',
		},
		{
			'spectra': ['EC59497'], 
			'calibs':  ['EC59498'], 
			'ID': 0,
			'obj_name': 'EC59497',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'left',
        'combine_exp': False
},
],
'RVS_COMP':[
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC53824', 'EC54732', 'EC63245', 'EC63247', 'EC63263'], 
			'calibs':  ['EC53825', 'EC54733', 'EC63246', 'EC63248', 'EC63264'], 
			'ID': 0,
			'obj_name': 'Alpha Boo',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC57055', 'EC63255', 'EC63267'], 
			'calibs':  ['EC57057', 'EC63256', 'EC63268'], 
			'ID': 1,
			'obj_name': 'Alpha Cas',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC57289'], 
			'calibs':  ['EC57288'], 
			'ID': 2,
			'obj_name': 'Alpha Cet',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC53438', 'EC53828', 'EC63241', 'EC63261'], 
			'calibs':  ['EC53437', 'EC53830', 'EC63242', 'EC63262'], 
			'ID': 3,
			'obj_name': 'HR 5694',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
{	
	'DATE_DIR': 'RVS_STAND',
	'ID_field': 'no_id',
	'biases': [], 
	'flats': [], 
	'objects':
		[
		{
			'spectra': ['EC53435', 'EC63215', 'EC63237', 'EC63239'], 
			'calibs':  ['EC53436', 'EC63216', 'EC63238', 'EC63240'], 
			'ID': 4,
			'obj_name': 'HR 6056',
		},
		],
	'REF_ARC':None,  
	'REF_AP':'', 
	'ap_position':'center',
        'combine_exp': False
},
],
}

'''

'''

AP_LIM = {'left': [-20, 10, '-19:-11'], 'center': [-15, 15, '-15:-10,10:15']}

root_dir_data = '/home/nandir/gigli2/media/hdd/home2/janez/ftp/asiago_observations/'
root_dir_arcs = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files_manual/'
#root_dir_arcs = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_wav_files/'
root_dir_apref = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_ap_files/'
root_dir_rvref = '/home/nandir/IRAF_Echelle_Asiago/Reduction_pipeline/ref_rv_files/'

process_survey = 'RVS_COMP'
REWRITE = False  # cheks if result is already present and skips processing of that stage accordingly
REREDUCE = False  # run reduction procedure again?
REANALYSE = True  # does RV and parameter determination have to be run again if results already exist


def main():
  work_path = 'observations/' + process_survey
    
  if not os.path.exists(work_path): 
    os.makedirs(work_path)

  if not os.path.exists(work_path+'/database'): 
    os.makedirs(work_path+'/database')

  if not os.path.exists(work_path+'/cosmics'): 
    os.makedirs(work_path+'/cosmics')

  #iraf.cd(work_path)
  #for i_o_w in range(31):
  #  export_spectrum_order_to_txt('EC59444.ec', order=i_o_w)
  #raise SystemExit

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
    res_tab = Table(names=['Asiago_id', 'dir', 'obj_name', 'Asiago_arc_id', survey_obs[0]['ID_field'], 'MJD', 'JD', 'rv', 'e_rv', 'teff', 'e_teff', 'feh', 'e_feh', 'logg', 'e_logg'], 
                    dtype=['S35', 'S6', 'S35', 'S35', 'int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64'])

  # Process individually every month/date/day/etc
  for i_d in range(len(survey_obs)):
     survey_cur_date = survey_obs[i_d]
     
     root_dir_cur_date = root_dir_data + survey_cur_date['DATE_DIR'] + '/'
     # preapre data that are needed for reduction of every spectra in selected period
     if len(survey_cur_date['biases']) > 0: 
       mbias = create_masterbias(survey_cur_date, root_dir_cur_date, rewrite=REWRITE)    
     if len(survey_cur_date['flats']) > 0:
       mflat = create_masterflat(survey_cur_date, root_dir_cur_date, rewrite=REWRITE)    
       # mflat_norm = normalize_masterflat(survey_cur_date, rewrite=REWRITE)  

     # Process with actual spectra reduction
     # Individually reduce every object observed on this date
     for i_o in range(len(survey_cur_date['objects'])):
       obj_spectra_data = survey_cur_date['objects'][i_o]

       folder_cleaned = clean_folder_write_header(obj_spectra_data, root_dir_cur_date, rewrite=REWRITE)
       if folder_cleaned:
         if len(survey_cur_date['biases']) > 0:
           bias_correct(obj_spectra_data, mbias) 
         # if len(survey_cur_date['flats']) > 0: 
         #   flat_correct(obj_spectra_data, mflat_norm)  
         # remove_cosmics(obj_spectra_data) 

       # spectrum extraction
       if survey_cur_date['REF_AP'] is not '':
         ref_ap_use = root_dir_apref + survey_cur_date['REF_AP']
       else:
         ref_ap_use = None
       AP_LIM_use = AP_LIM[survey_cur_date['ap_position']]
       do_apall_for_real(obj_spectra_data, AP_LIM_use, rewrite=REWRITE, ref_ap=ref_ap_use)  
       apply_apall(obj_spectra_data, AP_LIM_use, rewrite=REWRITE)  

       # determine wavelength solution
       if process_survey == 'WAV_REF_GEN':
         wavelength_solution_ref_create(obj_spectra_data)
         continue
 
       wavelength_solution(obj_spectra_data, root_dir_arcs, survey_cur_date['REF_ARC'], rewrite=REWRITE) 
       apply_wav_solution(obj_spectra_data, include_cal=True, rewrite=REWRITE)
       vhelio_correct(obj_spectra_data, rewrite=REWRITE)   

       # check_spectra(obj_spectra_data, '.ec.vh')
       final_norm_arc_1d = combine_normalize_arcs(obj_spectra_data, normalize=True, rewrite=REWRITE)
       final_norm_orders, final_norm_1d = combine_normalize_images(obj_spectra_data, combine=survey_cur_date['combine_exp'], rewrite=REWRITE)

       # TODO: check if spectra was already analyzed
       if 'RV_ref' in obj_spectra_data:
         rv_ref_val = obj_spectra_data['RV_ref']
       else:
         rv_ref_val = np.nan
       #rvref_spec_list = glob(root_dir_rvref+'*M05V000K2SNWNVR20N.fits')
       #rvref_spec_list = [rf.split('/')[-1].split('.')[0] for rf in rvref_spec_list]
       rvref_spec_list = ['solar']
       rv_med, rv_std = get_RV_custom_corr_perorder(final_norm_orders, root_dir_rvref, rvref_spec_list, plot_rv=True, ref_val=rv_ref_val)
       print '    Final RV values:', rv_med, rv_std
       #get_RV(obj_spectra_data, root_dir_rvref, rvref_spec_list, multi_ref=False, multi_sample=False, ref_val=rv_ref_val)

       # add results to table
       for i_f, final_img in enumerate(final_norm_orders):
         idx_row_exist = np.where(res_tab['Asiago_id'] == final_img)[0]
         # first remove old results
         if len(idx_row_exist) > 0:
           res_tab.remove_rows(idx_row_exist)
         # add new reults
         mjd, jd = get_julia_dates_header(final_img)
         res_list = [final_img, survey_cur_date['DATE_DIR'], obj_spectra_data['obj_name'], final_norm_arc_1d[i_f], obj_spectra_data['ID'],
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
      ref = ref_ap.split('/')[-1].replace('ap', '')
    else:
      ref = obj
    print 'APALL',obj, ref
    iraf.apall(input=obj, referen=ref, format='echelle', interac='yes', find='yes', recente='yes', resize='yes', edit='yes', trace='yes', fittrac='yes', extract='no', extras='no', review='no', line=500, nsum=50, lower=-5, upper=5, width=5, radius=9, thresho=0.2, nfind=31, minsep=25, maxsep=255, t_func='chebyshev', t_order=11, t_niter=5, t_low_reject=2.5, t_high_reject=2.5, t_nsum=20, t_step=10, backgro='none', ylevel='INDEF', llimit=-5, ulimit=5, order='increasing')#, apidtable='apert.txt')

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
      iraf.apall(input=obj, output=obj+j[0], format='echelle', referen=obj, interac='no', find='no', recente='no', resize='no', edit='no', trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro=j[1], b_order=1, b_sampl=AP_LIM[2], ylevel='INDEF', llimit=-5, ulimit=5, weights='variance', clean='yes', pfit='fit1d', order='increasing')#, apidtable='apert.txt')  

  for i, cal in enumerate(obs_data['calibs']):
    print ' - Processing arc %s \n' % cal

    try: os.remove(cal+'.ec.fits')
    except: pass 
    print '   input:',cal,'    referen:', obs_data['spectra'][i]
    iraf.apall(input=cal, format='echelle', referen=obs_data['spectra'][i], interac='no', find='no', recente='no', resize='no', edit='no',trace='no', fittrac='no', extract='yes', extras='no', review='no', line=600, nsum=30, lower=-5, upper=5, width=8, radius=10, thresho=0.1, nfind=35, minsep=10, maxsep=255, t_func='chebyshev', t_order=9, t_niter=5, backgro='none', ylevel='INDEF', llimit=-5, ulimit=5, order='increasing', weights='variance', pfit='fit1d')#, apidtable='apert.txt')    

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

  if not rewrite:
    n_miss = 0
    for obj in obs_data['calibs']:
      if not os.path.isfile('database/ec'+obj+'.ec'):
        n_miss += 1
    if n_miss == 0:
      print '  -- Wavelength solutions already found.'
      return False

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
    # remove left-overs from trying to fit multiple arc sollutions
    try:
      os.remove('database/ec' + cal + '.ec')
    except:
      pass

    # select the best matching arc with the lowest RMS
    iraf.unlearn('ecreidentify')
    idx_arc = np.nanargmin(arc_rms_all)
    arc_full_use = arc_full_path_list[idx_arc]
    print ' All arc RMS:', arc_rms_all
    print ' Best:', arc_full_use.split('/')[-1], arc_rms_all[idx_arc]
    iraf.ecreident(images=cal+'.ec', referenc=arc_full_use.replace('/ec', '/').split('/')[-1], refit='yes', shift='INDEF', cradius=5, thresho=10)
    
    print 'iraf.ecident'
    # lower number of maxfeat for less work while removing bad emissions in arc, incresed minsep
    iraf.ecident(images=cal+'.ec', coordli='linelists$thar.dat', match=0.025, maxfeat=1200, ftype='emission', fwidth=4, cradius=5, thresho=2, minsep=3, functio='chebyshev', xorder=6, yorder=6, niterat=7, lowreje=2.5, highreje=2.5, autowri='yes')

def apply_wav_solution(obs_data, include_cal=False, rewrite=True):
  print '\n + Applying the wavelength solution\n'

  if not rewrite:
    # TODO: proper check for this
    # look into fits header
    pass
  
  iraf.unlearn('dispcor')

  for i, obj in enumerate(obs_data['spectra']):
    cal = obs_data['calibs'][i]+'.ec'
    # export wav calibrated spectra
    print 'Using cal:', cal
    for j in ['.ec', '.ec.nosky']:
      iraf.refspectra(input=obj+j, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=obj+j, output=obj+j, lineari='no', verbose='yes')
    if include_cal:
      try:
        os.remove(cal+'_wvl.fits')
      except:
        pass
      # export wav calibrated cal image
      iraf.refspectra(input=cal, referen=cal, sort='', group='', confirm='no', Stdout="/dev/null")
      iraf.dispcor(input=cal, output=cal+'_wvl', lineari='no', verbose='yes')
      # export wvl calibratied arc to txt file
      # print 'Exporting wvl calibrated arc spectrum'
      # for i_o_w in range(31):
      #   export_spectrum_order_to_txt(cal+'_wvl', order=i_o_w)

def vhelio_correct(obs_data, rewrite=True):
  print '\n + VHELIO correction\n'

  if not rewrite:
    n_miss = 0
    for obj in obs_data['spectra']:
      obj = obj+'.ec'
      if not os.path.isfile(obj+'.vh.fits'):
        n_miss += 1
    if n_miss == 0:
      print '  -- VHELIO correction already performed.'
      return False
    
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
    iraf.dopcor(input=obj+'.vh', output=None, redshift='-vhelio', isveloc='yes', dispers='yes', verbose='yes')

    obj = obj+'.nosky'
    shutil.copy(obj+'.fits', obj+'.vh.fits')
    iraf.hedit(images=obj+'.vh', fields='UT', value=h+m/60.+s/3600., add='yes', verify='no')
    iraf.hedit(images=obj+'.vh', fields='EPOCH', value=year, add='yes', verify='no')
    iraf.rvcorrect(images=obj+'.vh', imupdat='yes', epoch=year, observa='ekar', year=year, month=month, day=day, ut=h+m/60.+s/3600., ra=ra, dec=dec)
    iraf.dopcor(input=obj+'.vh', output=None, redshift='-vhelio', isveloc='yes', dispers='yes', verbose='yes')

def combine_normalize_arcs(obs_data, normalize=True, rewrite=True):

      print '\n + Normalization of arc spectra\n'

      final_norm_files = list([])
      for obj in obs_data['calibs']:
        final_norm_files.append(obj+'_vh_norm')

        if not rewrite:
          n_done = 0
          for chk_file in [obj+'_cont.fits', obj+'_vh_norm.fits', obj+'_1D_vh_norm.0001.fits', obj+'_data1D.fits', obj+'_cont1D.fits']:
            if os.path.isfile(chk_file):
              n_done += 1
          if n_done == 5:
            print '  -- Normalization already applied to arc spectrum.'
            continue

        try: 
          os.remove(obj+'_cont.fits')
          os.remove(obj+'_vh_norm.fits')
          os.remove(obj+'_data1D.fits')
          os.remove(obj+'_cont1D.fits')
          os.remove(obj+'_1D_vh_norm.0001.fits') 
        except: 
          pass 

        iraf.continuum(input=obj+'.ec_wvl', output=obj+'_cont', type='fit', replace='no', listonly='no', functio='cheb', order=7, low_rej=4, high_rej=2, naverag=-3, niter=13, interac='no', markrej='no', ask='yes')
        iraf.sarith(input1=obj+'.ec_wvl', op='/', input2=obj+'_cont', output=obj+'_vh_norm', format='multispec', Stdout="/dev/null")

        iraf.scombine(input=obj+'.ec_wvl', output=obj+'_data1D', group='all', combine='sum', reject='none', Stdout="/dev/null")
        iraf.scombine(input=obj+'_cont', output=obj+'_cont1D', group='all', combine='sum', reject='none', Stdout="/dev/null")

        iraf.sarith(input1=obj+'_data1D', op='/', input2=obj+'_cont1D', output=obj+'_1D_vh_norm', format='onedspec', Stdout="/dev/null")

      return final_norm_files

def combine_normalize_images(obs_data, combine=False, normalize=True, rewrite=True):
  
  if combine:
    print '\n + Combine images\n'
    comb_sum_file = obs_data['spectra'][0]+'_vh_comb'
    comb_sum_norm_file = obs_data['spectra'][0]+'_vh_norm_comb'
    comb_sum_cont_file = obs_data['spectra'][0]+'_cont_comb'
    
    n_done = 0
    if not rewrite:
      for chk_file in [comb_sum_file, comb_sum_norm_file, comb_sum_cont_file]:
        if os.path.isfile(chk_file):
          n_done += 1
    if n_done < 3:
      for sum_file in [comb_sum_file, comb_sum_norm_file, comb_sum_cont_file]:
        try: os.remove(sum_file+'.fits')       
        except: pass
      iraf.scombine(input=','.join([obj+'.ec.vh' for obj in obs_data['spectra']]), output=comb_sum_file, group='apertures', combine='sum', reject='none', Stdout="/dev/null")
    else:
      normalize = False
      print '  -- Normalization already applied to combined spectra.'

  else:
    print '\n + Imdividual spectra will not be combined\n'
    # if spectra are not to be combined
    final_norm_files = list([])
    n_miss = 0
    for obj in obs_data['spectra']:
      final_norm_files.append(obj+'_vh_norm')
      for chk_file in [obj+'_cont.fits', obj+'_vh_norm.fits']:
        if not os.path.isfile(chk_file):
          n_miss += 1
      if not rewrite and n_miss == 0:
        normalize = False
        print '  -- Normalization already applied to all individual spectra.'      
     
  # normalize spectra
  if normalize:
    # This step will require you to manually normalize all the spectra.
    print '\n + Normalize spectra\n'

    if combine:
      try: 
          os.remove('combined_sum_cone_rv_echellet.fits')       
      except: 
          pass
      iraf.continuum(input=comb_sum_file, output=comb_sum_cont_file, type='fit', replace='no', listonly='no', functio='cheb', order=13, low_rej=2, high_rej=3, naverag=-3, niter=9, interac='no', markrej='no')
      iraf.sarith(input1=comb_sum_file, op='/', input2=comb_sum_cont_file, output=comb_sum_norm_file, format='multispec', Stdout="/dev/null")
    else: 
      for obj in obs_data['spectra']:
        try: 
          os.remove(obj+'_cont.fits')
          os.remove(obj+'_vh_norm.fits')
        except: 
          pass 
        iraf.continuum(input=obj+'.ec.vh', output=obj+'_cont', type='fit', replace='no', listonly='no', functio='cheb', order=11, low_rej=2, high_rej=3, naverag=-3, niter=13, interac='no', markrej='no', ask='yes')
        iraf.sarith(input1=obj+'.ec.vh', op='/', input2=obj+'_cont', output=obj+'_vh_norm', format='multispec', Stdout="/dev/null")

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

def get_RV(obs_data, rvref_dir, rv_ref_list, multi_ref=False, multi_sample=False, ref_val=np.nan):

  for obj in obs_data['spectra']:
    for rv_ref in rv_ref_list:
      ts = time.time()
      rvs = []
      e_rvs = []
      # determine barycentric RV
      with open(obj+'_'+rv_ref+'_rvs.txt', 'w+') as txt_file:

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
            
              iraf.fxcor(obj+'.ec.vh', rvref_dir+rv_ref, apertures=order, pixcorr='no', continuum='both', 
                         osample='a %s %s' % (wav_r[0], wav_r[1]), rsample='a %s %s' % (wav_r[0], wav_r[1]), 
                         interactive='no', output='test', function = 'gaussian',
                         apodize = 0.0, rebin = "smallest")
              fields = np.loadtxt('test.txt', comments='#', dtype='str')
              # print wav_r
              # print fields
              vrel, verr = fields[-3], fields[-1]

              if vrel == 'INDEF' or np.float(verr) > 10. or np.abs(np.float(vrel)) > 350.:
                vrel = np.nan
                verr = np.nan

              print order, vrel, verr              
              txt_file.write("%s %s %s\n" % (order, vrel, verr))
          
              rvs.append(float(vrel))
              e_rvs.append(float(verr))

          except:
            print_exception()

      rv_median = np.nanmedian(rvs)
      plt.scatter(np.arange(len(rvs)),rvs)
      plt.axhline(rv_median)
      plt.ylim(rv_median-4., rv_median+4.)
      plt.xlim(0-1., len(rvs)+1.)
      if np.isfinite(ref_val):
        plt.axhline(ref_val, label='Ref RV', c='black', ls='--')
      plt.xlabel('Echelle order')
      plt.ylabel('Radial velocity')
      plt.title('Spectrum: '+obj+'     Reference: '+rv_ref)
      plt.tight_layout()
      plt.savefig(obj+'_'+rv_ref+'_rv_fxcor.png', dpi=250)
      plt.close()
        
      print '%s: RV median: %s, RV mean: %s, RV sigma: %s' % (obj, np.nanmedian(rvs), np.nanmean(rvs), np.nanstd(rvs))
      print 'RV fit time: '+str((time.time()-ts)/60.)+' min'
      print ''

def get_RV_custom_corr_perorder(obs_spectra, rvref_dir, rvref_list, plot_rv=False, n_min_ord=1, n_max_ord=31, ref_val=np.nan):

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
      with open(obs_spectrum+'_'+use_revref+'_rv.txt', 'w+') as txt_file:
        for i_order in range(n_min_ord, n_max_ord+1):
          order_txt_file = export_spectrum_order_to_txt(obs_spectrum, order=i_order)
          # print order_txt_file
          if order_txt_file is not None:
            rv_order_val, cent_wvl = correlate_order(order_txt_file, ref_flx, ref_wvl, plot=False)
          else:
            rv_order_val, cent_wvl = np.nan, np.nan
          rv_shifts.append(rv_order_val)
          txt_file.write("%s,%s,%s\n" % (i_order, rv_order_val, cent_wvl))
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
        if np.isfinite(ref_val):
          plt.axhline(ref_val, label='Ref RV', c='black', ls='--')
        plt.xlabel('Echelle order')
        plt.ylabel('Radial velocity')
        plt.title('Spectrum: '+obs_spectrum+'     Reference: '+use_revref)
        plt.tight_layout()
        plt.savefig(obs_spectrum+'_'+use_revref+'_rv.png', dpi=250)
        plt.close()

    if np.sum(np.isfinite(rv_obs_std)==0):
      # return nan if no valid RV was determined
      rv_obs_mean_final.append(np.nan)
      rv_obs_std_final.append(np.nan)
    else:
      # returns the best (lovest std between on orders) RV of all retrieved
      idx_best_rv = np.nanargmin(rv_obs_std)
      rv_obs_mean_final.append(rv_obs_mean[idx_best_rv])
      rv_obs_std_final.append(rv_obs_std[idx_best_rv])

  # return all derived velocities and its std values
  return rv_obs_mean_final, rv_obs_std_final

def print_exception():
  e = sys.exc_info()
  exc_type, exc_value, exc_traceback = e
  a, j = (traceback.extract_tb(exc_traceback, 1))[0][0:2]
  k = (traceback.format_exception_only(exc_type, exc_value))[0]
  print a, j, k

if __name__ == "__main__":
  main()

