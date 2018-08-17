import os

in_dir = '/home/gregor/public_html/Astro_predmeti/normalized_spectra_R20000/'
suffix = 'V000K2SNWNVR20N.ASC'

for feh in ['M10', 'M05', 'P00', 'P05']:
    for teff in ['03500', '04500', '05500', '06500']:
        for logg in ['10','20','30','40','50']:
            # example: T03500G00M05V000K2SNWNVR20N.ASC
            path = in_dir + 'T_' + teff + '/'
            filename_get = 'T'+teff+'G'+logg+feh+suffix
            os.system('cp '+path+filename_get+' '+filename_get)