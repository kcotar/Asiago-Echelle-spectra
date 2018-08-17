from astropy.table import Table
import numpy as np

gaia_rv = Table.read('rvstdcat.csv', format='ascii.no_header', delimiter=',')
table_cols = ['ID','RV','eRV','s_RV','e_RV','Tbase','JDm','N','Source','Flag','RVS','e_RVS','o_RVS','2MASS',
              'Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag','Qflag','RAdeg','DEdeg','Epoch','Bmag','Vmag','SpType','otype']
for i_n in range(len(gaia_rv.colnames)):
    gaia_rv[gaia_rv.colnames[i_n]].name = table_cols[i_n]

# use flags if needed
# CAL1 = calibrator for DR2 and DR3
# gaia_rv = gaia_rv[gaia_rv['Flag']=='CAL1']

# filter by obsrvability
asiago_rv = gaia_rv[gaia_rv['DEdeg'] > -20.]
asiago_rv = asiago_rv[asiago_rv['RAdeg'] > 200.]
asiago_rv = asiago_rv[asiago_rv['RAdeg'] < 350.]
asiago_rv = asiago_rv[asiago_rv['Vmag'] < 11.]


# sort by number of observations
asiago_rv = asiago_rv[np.argsort(asiago_rv['N'])[::-1]]
data_out = asiago_rv['ID','RV','s_RV','e_RV','N','Flag','RVS','e_RVS','o_RVS','RAdeg','DEdeg','Bmag','Vmag','SpType','otype'][:150]

print data_out
data_out.write('rvstdcat_asiago_june18.csv', format='ascii.fixed_width_two_line', overwrite=True)
