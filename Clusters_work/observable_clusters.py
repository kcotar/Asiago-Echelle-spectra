from astropy.table import Table
import numpy as np
import astropy.coordinates as coord
import astropy.units as un
from astropy.time import Time
from astroplan import Observer

# update IERS Bulletin A table
from astroplan import download_IERS_A
download_IERS_A()

c_data = Table.read('/data4/cotar/clusters/Kharchenko_2013/catalog.csv')

# assiago coordinate
site = 'Mt. Ekar 182 cm. Telescope'

# requested day
day = 25
month = 11
year = 2018
time_day_start = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day), scale='utc')  # begin the same day
time_day_end = Time('{:04.0f}-{:02.0f}-{:02.0f} 12:00:00'.format(year, month, day+1), scale='utc')  # end next day

prefix = '{:04.0f}_{:02.0f}_{:02.0f}_'.format(year, month, day)

out_file = prefix+'visibility_table.txt'
our_file_2 = prefix+'observable_clusters.txt'

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
alt_min = 30
nmin_above = 3
txt_clusters = open(our_file_2, 'w')
txt_visibi = open(out_file, 'w')
# add some header information
txt_clusters.write('Observable clusters and their coordinates \n')
txt_clusters.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_clusters.write('=============================================================================\n\n\n')
txt_visibi.write('Cluster visibility and altitude in brackets at given UT times.\n')
txt_visibi.write('Night between '+obs_start.utc.iso+' and '+obs_end.utc.iso+' UTC (civil twilight)\n')
txt_visibi.write('=============================================================================\n\n\n')

use_cluster_from = ['NGC', 'IC', 'ASCC', 'Melotte', 'Ruprecht', 'Platais', 'Blanco']

all_observable = list([])
all_interesting_clusters = list([])
for cluster in c_data:
    print cluster['Cluster']

    # check if the cluster is accepted
    n_ok = 0
    for c_list in use_cluster_from:
        if c_list in cluster['Cluster']:
            n_ok += 1
    if n_ok == 0:
        continue

    all_interesting_clusters.append(cluster['Cluster'])

    if cluster['d'] > 1000:
        continue

    observable_flag = list([])
    alt_at_hour = list([])
    for cur_time in chect_at_timestamps:
        cluster_pos = coord.ICRS(ra=cluster['RAdeg']*un.deg, dec=cluster['DEdeg']*un.deg)
        cluster_altaz = cluster_pos.transform_to(coord.AltAz(obstime=Time(cur_time, scale='utc'), location=ekar_location))
        az = cluster_altaz.az.value
        alt = cluster_altaz.alt.value
        observable_flag.append(alt >= alt_min)
        alt_at_hour.append('{:.0f}'.format(alt))

    # determine if it is observable
    idx_obs = np.array(observable_flag)
    if np.sum(idx_obs) >= nmin_above:
        all_observable.append(cluster['Cluster'])
        txt_visibi.write(cluster['Cluster'] + '\n')
        txt_clusters.write(cluster['Cluster'] + '  ' + str(cluster['RAdeg']) + '  ' + str(cluster['DEdeg']) + '\n')
        # export it to the list
        obs_time_alt_str = ''
        for h_p in np.where(idx_obs)[0]:
            obs_time_alt_str += str(utc_hours_list[h_p])+'h('+str(alt_at_hour[h_p])+') '
        print obs_time_alt_str
        txt_visibi.write(obs_time_alt_str+'\n\n')

txt_clusters.write('All: \n"' + '","'.join(all_observable)+'"\n\n')
txt_clusters.write('All interesting: \n"' + '","'.join(all_interesting_clusters)+'"\n\n')
txt_clusters.close()
txt_visibi.close()
