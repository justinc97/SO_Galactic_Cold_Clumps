import numpy as np
import healpy as hp
import pandas as pd
import pysm3.units as u

from .galcc_utils import (modBB, 
                        arclength_sphere, 
                        gaussian_source
)


class galcc_mapper(object):
    def __init__(self):
        self.nside = 0
        self.npix = 0
        self.map = None
        self.per_pix_steradian = 0 
        
    def galcc_map(self, catalogue, freq_out, nside, store_maps = False, outdir = "", output_units = u.uK_CMB):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.per_pix_steradian = 1 / hp.nside2pixarea(self.nside) # 1/sr
        # Build empty map
        self.map = np.zeros(self.npix, dtype = np.float64)
        
        # Check catalogue exists in working directory
        try:
            df = pd.read_csv(catalogue)
        except:
            print("Catalogue not in working directory: ", os.getcwd() + catalogue)
            return
        # Restrict sources to flux_quality = 1 so T and beta are available
        df = df[df.flux_quality == 1]
        df = df.reset_index(drop = True)
        nsources = len(df)
        
        # Original flux values at 353 GHz
        planck_flux = np.array(df['flux_353_clump'])
        
        # Scale flux to desired frequency out
        scaled_flux = modBB(planck_flux, np.array(df['temp_clump']), freq_out, 353, np.array(df['beta_clump']))
        
        # Add sources to map at given locations with appropriate size and Gaussian flux distribution
        vecs = hp.ang2vec(df['glon'], df['glat'], lonlat = True)
        radii = (np.array(df['gau_major_axis']) * u.arcmin).to_value(u.rad)
        for i in range(nsources):
            pix_circ = hp.query_disc(nside = self.nside, vec = vecs[i], radius = 3 * radii[i])
            theta_pix, phi_pix = hp.pix2ang(self.nside, pix_circ)
            theta_cent, phi_cent = hp.vec2ang(vecs[i])
            pix_dists = arclength_sphere(theta_pix, phi_pix, theta_cent, phi_cent)
            flux_profile = scaled_flux[i] * gaussian_source(pix_dists, radii[i])
            self.map[pix_circ] += flux_profile # Jy
            
        # Convert map from Jy -> Jy/sr with appropriate pix area
        self.map *= self.per_pix_steradian # Jy/sr
        
        if output_units == u.uK_CMB or output_units == u.K_CMB:
            self.map = (self.map * u.Jy/u.sr).to_value(output_units, equivalencies = u.cmb_equivalencies(freq_out * u.GHz))
        else:
            self.map = (self.map * u.Jy/u.sr).to_value(output_units)
        
        if store_maps == True:
            hp.write_map(outdir + str(freq_out) + "_GHz_GCC_Map.fits", m = self.map, coord = "G", column_units = str(output_units))
            return self.map
        else:
            return self.map