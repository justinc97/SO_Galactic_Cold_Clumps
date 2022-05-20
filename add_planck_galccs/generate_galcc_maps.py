import numpy as np
import healpy as hp
import pandas as pd
import pysm3.units as u
import os
from .galcc_utils import (modBB, 
                          arclength_sphere, 
                          gaussian_source,
                          gaussian_source_circ,
                          map_unit_conversion
)


class build_galactic_clump_map(object):
    def __init__(self, nside, catalogue, b_range = None):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.map = np.zeros(self.npix, dtype = np.float64)
        self.per_pix_steradian = 1 / (hp.nside2pixarea(self.nside))

        try:
            df = pd.read_csv(catalogue)
        except:
            print("Catalogue not in working directory: ", os.getcwd() + catalogue)
            return
        # Restrict PGCC sources to flux_quality = 1 ensuring they have a measured
        # T and beta value
        df = df[df.flux_quality == 1]
        
        if b_range is not None:
            b_min = b_range[0]
            b_max = b_range[1]
            # Restrict sources under |b| > 30 for SO Visible
            df = df[df['glat'] < b_max]
            df = df[df['glat'] > b_min]

        
        df = df.reset_index(drop = True)
        
        self.Nsources = len(df)
        self.df = df
        
        # Read spectral and flux properties
        self.temperatures = np.array(df['temp_clump'])
        self.betas        = np.array(df['beta_clump'])
        self.flux_353     = np.array(df['flux_353_clump'])
        self.flux_545     = np.array(df['flux_545_clump'])
        self.flux_857     = np.array(df['flux_857_clump'])
        
        # Read source position, shape and orientation
        vecs = hp.ang2vec(df['glon'], df['glat'], lonlat = True)
        FWHM_maj = (np.array(df['gau_major_axis']) * u.arcmin).to_value(u.rad)
        FWHM_min = (np.array(df['gau_minor_axis']) * u.arcmin).to_value(u.rad)
        pos_ang = df['gau_position_angle']
        
        profiles = []
        sources = []
        for i in range(self.Nsources):
            pix_circ = hp.query_disc(nside = self.nside, vec = vecs[i],
                                     radius = 3 * FWHM_maj[i])
            theta_pix, phi_pix = hp.pix2ang(self.nside, pix_circ)
            theta_cent, phi_cent = hp.vec2ang(vecs[i])

            profile = gaussian_source((theta_pix, theta_cent),
                                          (phi_pix, phi_cent),
                                          FWHM_maj[i],
                                          FWHM_min[i],
                                          pos_ang[i])
            profiles.append(profile)
            sources.append(pix_circ)
        self.profiles = profiles
        self.sources = sources
        
    def mask_sources(self, threshold = 1e-2, store_maps = False, outdir = ""):
        m = self.map.copy() + 1
        
        flatprofiles = []
        for i in range(self.Nsources):
            maskval = []
            for j in range(len(self.profiles[i])):
                val = self.profiles[i][j]
                if val <= threshold:
                    val = 1
                elif val >= threshold:
                    val = 0
                maskval.append(val)
            flatprofiles.append(maskval)
        
        for i in range(self.Nsources):
            m[self.sources[i]] *= flatprofiles[i]
        
        if store_maps == True:
            hp.write_map(outdir + "PGCC_Source_Mask.fits", m = m, coord = "G",
                         overwrite = True)
        return m

    def cold_clumps_spectral(self, maptype, store_maps = False, outdir = ""):

        if maptype == 'temp':
            spectral_property = self.temperatures
            output_units = u.K
        elif maptype == 'beta':
            spectral_property = self.betas
            output_units = u.dimensionless_unscaled

        m = self.map.copy()
        for i in range(self.Nsources):
            source_profile = spectral_property[i] * self.profiles[i]

            m[self.sources[i]] += source_profile 
            
        if store_maps == True:
            hp.write_map(outdir + str(freq_out) + "_GHz_GCC_" + str(maptype) + "_Map.fits", m = m, coord = "G", column_units = str(output_units), overwrite = True)
            
        return m
        
    def cold_clumps_flux(self, freq_out, output_units = u.uK_CMB,
                  store_maps = False, outdir = ""):
        
        # Original flux values at 353 GHz
        if freq_out == 353:
            scaled_flux = self.flux_353
        elif freq_out == 545:
            scaled_flux = self.flux_545
        elif freq_out == 857:
            scaled_flux = self.flux_857
        else: 
            scaled_flux = modBB(np.array(self.flux_353), 
                                np.array(self.temperatures), 
                                freq_out, 
                                353, 
                                np.array(self.betas))
        m = self.map.copy()
        for i in range(self.Nsources):
            source_profile = scaled_flux[i] * self.profiles[i]
            m[self.sources[i]] += source_profile # Jy

            
        # Convert map from Jy -> Jy/sr with appropriate pix area
        m *= self.per_pix_steradian # Jy/sr
        
        m = map_unit_conversion(m * u.Jy / u.sr, freq_out, output_units)

        if store_maps == True:
            hp.write_map(outdir + str(freq_out) + "_GHz_GCC_Map.fits", m = m, coord = "G", column_units = str(output_units), overwrite = True)

        return m
