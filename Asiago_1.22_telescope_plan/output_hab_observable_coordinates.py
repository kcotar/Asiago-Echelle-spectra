from astropy.table import Table, join
import numpy as np
import astropy.coordinates as coord
import astropy.units as un
from astropy.time import Time
from astroplan import Observer

# update IERS Bulletin A table
#from astroplan import download_IERS_A
#download_IERS_A()

st_data = Table.read('molec_new_selection.csv')

print(st_data.colnames)
print('N objects:', len(st_data))

# asiago hard location limit
st_data = st_data[st_data['dec'] > -40.]
# magnitude limit to get only brighter stars
st_data = st_data[st_data['phot_g_mean_mag'] < 10.]

# filter by ra
st_data = st_data[st_data['ra'] > 0.]
st_data = st_data[st_data['ra'] < 150.]

out_data = st_data['sobject_id','ra','dec','phot_g_mean_mag']
out_data = out_data[np.argsort(out_data['ra'])]
out_data.write('galileo_telescope_targets_1912.txt', format='ascii.fixed_width_two_line', overwrite=True)

# assiago coordinate
site = 'Mt. Ekar 182 cm. Telescope'

# requested day
day = 12
month = 8
year = 2019
time_day_start = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day), scale='utc')  # begin the same day
time_day_end = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day+1), scale='utc')  # end next day

prefix = '{:04.0f}_{:02.0f}_{:02.0f}_'.format(year, month, day)

out_file = prefix+'galileo_telescope_targets_1912.txt'

ekar_site = Observer.at_site(site, timezone="UTC")
ekar_location = coord.EarthLocation.of_site(site)
obs_start = ekar_site.twilight_evening_civil(time_day_start, which=u'nearest')
obs_end = ekar_site.twilight_morning_civil(time_day_end, which=u'nearest')

hour_begin = int(obs_start.utc.iso.split(':')[0].split(' ')[-1])
hour_end = int(obs_end.utc.iso.split(':')[0].split(' ')[-1]) + 1

print('Civil twilight start:', obs_start.utc.iso)
print('Civil twilight end:  ', obs_end.utc.iso)
print(hour_begin, hour_end)

# determine time stamps to be considered
chect_at_timestamps = list([])
utc_hours_list = list([])
for hour in range(hour_begin, hour_end+24+1):
    if hour < 24:
        chect_at_timestamps.append('{:04.0f}-{:02.0f}-{:02.0f} {:02.0f}:00:00'.format(year, month, day, hour))
        utc_hours_list.append(hour)
    else:
        chect_at_timestamps.append('{:04.0f}-{:02.0f}-{:02.0f} {:02.0f}:00:00'.format(year, month, day+1, hour-24))
        utc_hours_list.append(hour-24)
print(chect_at_timestamps)

# determine which clusters are observable in this time intervals
alt_min = 20
nmin_above = 2
txt_clusters = open(our_file_2, 'w')
txt_visibi = open(out_file, 'w')
# add some header information
txt_clusters.write('Observable star at Galileo and their coordinates \n')
txt_clusters.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_clusters.write('=============================================================================\n\n\n')
txt_clusters.write('GAIA DR2 source_id     sobject_id         ra          dec         G_mag       Abs G_mag\n\n')

txt_visibi.write('Observable star at Galileo visibility and altitude in brackets at given UT times  (hh:00:00).\n')
txt_visibi.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_visibi.write('=============================================================================\n\n\n')
# txt_visibi.write('Civil twilight start:'+ obs_start.utc.iso+'\n')
# txt_visibi.write('Civil twilight end:  '+ obs_end.utc.iso+'\n\n\n')

all_observable = list([])
for st_obj in st_data:
    print(st_obj['sobject_id'])
    star_pos = coord.ICRS(ra=st_obj['ra'] * un.deg, dec=st_obj['dec'] * un.deg)
    print(ekar_site.target_meridian_antitransit_time(obs_start.utc, star_pos, which='next'))



raise SystemExit
'''
observable_flag = list([])
alt_at_hour = list([])
for cur_time in chect_at_timestamps:
    cluster_pos = coord.ICRS(ra=st_obj['ra']*un.deg, dec=st_obj['dec']*un.deg)
    cluster_altaz = cluster_pos.transform_to(coord.AltAz(obstime=Time(cur_time, scale='utc'), location=ekar_location))
    az = cluster_altaz.az.value
    alt = cluster_altaz.alt.value
    observable_flag.append(alt >= alt_min)
    alt_at_hour.append('{:.0f}'.format(alt))

# determine if it is observable
idx_obs = np.array(observable_flag)
if np.sum(idx_obs) >= nmin_above:
    all_observable.append(str(st_obj['sobject_id']))
    # check if star was already observed in Asiago
    n_asiago_obs = np.sum(a_data['sobject_id'] == st_obj['sobject_id'])
    if n_asiago_obs > 0:
        asiago_obs_str = ' - {:.0f} Asiago spectra'.format(n_asiago_obs)
    else:
        asiago_obs_str = ''
    # export visibility string
    txt_visibi.write(str(st_obj['sobject_id']) + asiago_obs_str + '\n')
    # export data string
    txt_clusters.write(str(st_obj['source_id']) + '  ' + str(st_obj['sobject_id']) + '  ' + '{:.9f}'.format(st_obj['ra']) + '  ' + '{:.9f}'.format(st_obj['dec']) + '    ' + '{:.2f}'.format(st_obj['phot_g_mean_mag']) + '    ' + '{:.2f}'.format(st_obj['abs_g_mean_mag']) + asiago_obs_str + '\n')
    # export it to the list
    obs_time_alt_str = ''
    for h_p in np.where(idx_obs)[0]:
        obs_time_alt_str += str(utc_hours_list[h_p])+'h('+str(alt_at_hour[h_p])+') '
    print obs_time_alt_str
    txt_visibi.write(obs_time_alt_str+'\n\n')

# txt_clusters.write('\n\n')
# txt_clusters.write('All: ' + ','.join(all_observable)+'\n')
txt_clusters.close()
txt_visibi.close()
'''
