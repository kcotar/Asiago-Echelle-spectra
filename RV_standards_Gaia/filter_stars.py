from astropy.table import Table
import numpy as np

# read catalog of RV standard stars as used by the Gaia RVS analysis and reduction team
gaia_rv = Table.read('rvstdcat.csv', format='ascii.no_header', delimiter=',')
table_cols = ['ID','RV','eRV','s_RV','e_RV','Tbase','JDm','N','Source','Flag','RVS','e_RVS','o_RVS','2MASS',
              'Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag','Qflag','RAdeg','DEdeg','Epoch','Bmag','Vmag','SpType','otype']
for i_n in range(len(gaia_rv.colnames)):
    gaia_rv[gaia_rv.colnames[i_n]].name = table_cols[i_n]

# use flags if needed
# CAL1 = calibrator for DR2 and DR3
# gaia_rv = gaia_rv[gaia_rv['Flag']=='CAL1']

# filter by observability
asiago_rv = gaia_rv[gaia_rv['DEdeg'] > -10.]
asiago_rv = asiago_rv[asiago_rv['RAdeg'] > 240.]
asiago_rv = asiago_rv[asiago_rv['RAdeg'] < 360.]
# prefer only bright stars as short exposures are needed, also have higher SNR
asiago_rv = asiago_rv[asiago_rv['Vmag'] < 10.]


# sort by number of observations
asiago_rv = asiago_rv[np.argsort(asiago_rv['N'])[::-1]]
data_out = asiago_rv['ID','RV','s_RV','e_RV','N','Flag','RVS','e_RVS','o_RVS','RAdeg','DEdeg','Bmag','Vmag','SpType','otype'][:150]

print data_out
# export printable version of the resulting table
data_out.write('rvstdcat_asiago_july19.csv', format='ascii.fixed_width_two_line', overwrite=True)
