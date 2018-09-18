from astropy.table import Table, join
import numpy as np
import astropy.coordinates as coord
import astropy.units as un
from astropy.time import Time
from astroplan import Observer

c_data = Table.read('/home/klemen/data4_mount/sobject_iraf_iDR2_180325_cannon.fits')
g_data = Table.read('/home/klemen/data4_mount/sobject_iraf_53_gaia.fits')['sobject_id','source_id','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']
print g_data.colnames

c_data = join(c_data, g_data, keys='sobject_id', join_type='left')

# all investigated as possible multiple star
# list_objects = [140314004401277,140608002501303,140805003101303,140808002701338,140808003701104,141102002701267,141103003101379,150408005301116,150409002601317,150411004101331,150413003601344,150413005101096,150606005901339,150703005601062,151227005201142,160107003101157,160125004501038,160125004501256,160325002701048,160327006101355,160330002601095,160401003901215,160401004401123,160521004801353,160530005501077,160531005101256,160531006101153,160602001601307,160916001801263,161008002501018,161009002601018,161009005901171,161117004601245,161118002601376,161210004201315,161211003101387,161212002101397,161217002601138,161217004101075,161217004101234,161219005101228,170117003101044,170121002801292,170408004501048,170507007801271,170508002601001,170508004801312,170514002401099,170514003001180,170514003301001,170906003601357,170909002601291,170911003101383,170912002901303,171205004101175,171207003601278,171227003601367]
# all selected to have >2 members
# list_objects = [140608002501303,140808002701338,140808003701104,141102002701267,150408005301116,150413003601344,150413005101096,160107003101157,160125004501038,160325002701048,160401003901215,160531006101153,161009005901171,161118002601376,161212002101397,161217004101075,161217004101234,170121002801292,170507007801271,170508002601001,170508004801312,170514002401099,170514003301001,170911003101383]
# selected to have 3 members
list_objects = [140608002501303,140808003701104,150413003601344,160125004501038,160401003901215,161009005901171,161118002601376,161217004101234,170507007801271,170514002401099,170911003101383]
print 'N objects:', len(list_objects)

st_data = c_data[np.in1d(c_data['sobject_id'], list_objects)]

# asiago filter
st_data = st_data[st_data['dec'] > -40]

out_data = st_data['sobject_id','galah_id','ra','dec','phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']

print st_data['sobject_id','ra','dec']
print out_data
out_data.write('solatwins-multiples_asiago_avgust18.csv', format='ascii.fixed_width_two_line', overwrite=True)

# assiago coordinate
site = 'Mt. Ekar 182 cm. Telescope'

# requested day
day = 28
month = 8
year = 2018
time_day_start = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day), scale='utc')  # begin the same day
time_day_end = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day+1), scale='utc')  # end next day

prefix = '{:04.0f}_{:02.0f}_{:02.0f}_'.format(year, month, day)

out_file = prefix+'visibility_table_Solar_multiples.txt'
our_file_2 = prefix+'observable_Solar_multiples.txt'


ekar_site = Observer.at_site(site, timezone="UTC")
ekar_location = coord.EarthLocation.of_site(site)
obs_start = ekar_site.twilight_evening_civil(time_day_start, which=u'nearest')
obs_end = ekar_site.twilight_morning_civil(time_day_end, which=u'nearest')

hour_begin = int(obs_start.utc.iso.split(':')[0].split(' ')[-1])
hour_end = int(obs_end.utc.iso.split(':')[0].split(' ')[-1]) + 1

print 'Civil twilight start:', obs_start.utc.iso
print 'Civil twilight end:  ', obs_end.utc.iso
print hour_begin, hour_end

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
print chect_at_timestamps

# determine which clusters are observable in this time intervals
alt_min = 20
nmin_above = 2
txt_clusters = open(our_file_2, 'w')
txt_visibi = open(out_file, 'w')
# add some header information
txt_clusters.write('Observable Solar Multiples and their coordinates \n')
txt_clusters.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_clusters.write('=============================================================================\n\n\n')
txt_clusters.write('GAIA DR2 source_id     sobject_id         ra          dec         G_mag\n\n')

txt_visibi.write('Solar Multiples visibility and altitude in brackets at given UT times  (hh:00:00).\n')
txt_visibi.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_visibi.write('=============================================================================\n\n\n')
# txt_visibi.write('Civil twilight start:'+ obs_start.utc.iso+'\n')
# txt_visibi.write('Civil twilight end:  '+ obs_end.utc.iso+'\n\n\n')

all_observable = list([])
for st_obj in st_data:
    print st_obj['sobject_id']

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
        txt_visibi.write(str(st_obj['sobject_id']) + '\n')
        txt_clusters.write(str(st_obj['source_id']) + '  ' + str(st_obj['sobject_id']) + '  ' + str(st_obj['ra']) + '  ' + str(st_obj['dec']) + '  ' + str(st_obj['phot_g_mean_mag']) + '\n')
        # export it to the list
        obs_time_alt_str = ''
        for h_p in np.where(idx_obs)[0]:
            obs_time_alt_str += str(utc_hours_list[h_p])+'h('+str(alt_at_hour[h_p])+') '
        print obs_time_alt_str
        txt_visibi.write(obs_time_alt_str+'\n\n')

txt_clusters.write('\n\n')
txt_clusters.write('All: ' + ','.join(all_observable)+'\n')
txt_clusters.close()
txt_visibi.close()

