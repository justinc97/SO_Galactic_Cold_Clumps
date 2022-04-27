import healpy as hp
import matplotlib.pyplot as plt
from add_planck_galccs.generate_galcc_maps import galcc_mapper

galccmap = galcc_mapper().galcc_map(catalogue = "PGCC.csv", 
                           freq_out = 93, 
                           nside = 2048, 
                           store_maps = True)

hp.gnomview(galccmap)
plt.show()
hp.mollview(galccmap)
plt.show()