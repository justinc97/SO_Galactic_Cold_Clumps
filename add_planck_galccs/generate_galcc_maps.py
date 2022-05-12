import numpy as np
import healpy as hp
import pandas as pd
import pysm3.units as u
import os
from .galcc_utils import (modBB, 
                          arclength_sphere, 
                          gaussian_source,
                          gaussian_source_circ
)


class galcc_mapper(object):
    def __init__(self):
        self.nside = 0
        self.npix = 0
        self.map = np.zeros(self.npix, dtype = np.float64)
        self.per_pix_steradian = 0 

    def spectral_map(self, catalogue, maptype, nside,
                     shape_circ = False, shape_dist = 'gaussian', 
                     store_maps = False, outdir = ""):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.per_pix_steradian = 1 / (hp.nside2pixarea(self.nside))
        self.map = np.zeros(self.npix, dtype = np.float64)
        try:
            df = pd.read_csv(catalogue)
        except:
            print("Catalogue not in working directory: ", os.getcwd() + catalogue)
            return
        # Restrict sources to flux_quality = 1 so T and beta are available
        df = df[df.flux_quality == 1]
        df = df[df['glat'] < 30] # Restrict to sources under |b| > 30 for SO visible region
        df = df.reset_index(drop = True)
                
        nsources = len(df)
        temp = df['temp_clump']
        beta = df['beta_clump']
        
        if maptype == 'temp':
            spectral_property = temp
            output_units = u.K
        elif maptype == 'beta':
            spectral_property = temp
            output_units = u.dimensionless_unscaled


        # Add sources to map at given locations with appropriate size and Gaussian property distribution
        vecs = hp.ang2vec(df['glon'], df['glat'], lonlat = True)
        FWHM_maj = (np.array(df['gau_major_axis']) * u.arcmin).to_value(u.rad)
        FWHM_min = (np.array(df['gau_minor_axis']) * u.arcmin).to_value(u.rad)
        pos_ang = df['gau_position_angle'] # rad
        
        for i in range(nsources):
            pix_circ = hp.query_disc(nside = self.nside, vec = vecs[i],
                                     radius = 3 * FWHM_maj[i])
            theta_pix, phi_pix = hp.pix2ang(self.nside, pix_circ)
            theta_cent, phi_cent = hp.vec2ang(vecs[i])
            
            if shape_circ == True:
                pix_dists = arclength_sphere(theta_pix, 
                                             phi_pix, 
                                             theta_cent, 
                                             phi_cent)
                profile = gaussian_source_circ(pix_dists, FWHM_maj[i])
            else:
                profile = gaussian_source((theta_pix, theta_cent),
                                          (phi_pix, phi_cent),
                                          FWHM_maj[i],
                                          FWHM_min[i],
                                          pos_ang[i])
            
            spectral_profile = spectral_property[i] * profile
            self.map[pix_circ] += spectral_profile # Jy
            
        if store_maps == True:
            hp.write_map(outdir + str(freq_out) + "_GHz_GCC_Map.fits", m = self.map, coord = "G", column_units = str(output_units), overwrite = True)
            return self.map
        else:
            return self.map
        
    def galcc_map(self, catalogue, freq_out, nside, output_units = u.uK_CMB,
                  shape_circ = False, store_maps = False, outdir = ""):
        self.nside = nside
        self.npix = hp.nside2npix(self.nside)
        self.per_pix_steradian = 1 / (hp.nside2pixarea(self.nside)) # 1/sr
        self.map = np.zeros(self.npix, dtype = np.float64)
        # Check catalogue exists in working directory
        try:
            df = pd.read_csv(catalogue)
        except:
            print("Catalogue not in working directory: ", os.getcwd() + catalogue)
            return
        # Restrict sources to flux_quality = 1 so T and beta are available
        df = df[df.flux_quality == 1]
        df = df[df['glat'] < 30]
        df = df.reset_index(drop = True)
        
        # Restrict to sources under |b| > 30 for SO visible region
        nsources = len(df)
        temp = np.array(df['temp_clump'])
        beta = np.array(df['beta_clump'])
        # Scale flux to desired frequency out
        
        # Original flux values at 353 GHz
        if freq_out == 353:
            scaled_flux = df['flux_353_clump']
        elif freq_out == 545:
            scaled_flux = df['flux_545_clump']
        elif freq_out == 857:
            scaled_flux = df['flux_857_clump']
        else: 
            scaled_flux = modBB(np.array(df['flux_353_clump']), 
                                temp, 
                                freq_out, 353, 
                                beta)

        
        # Add sources to map at given locations with appropriate size and Gaussian flux distribution
        vecs = hp.ang2vec(df['glon'], df['glat'], lonlat = True)
        FWHM_maj = (np.array(df['gau_major_axis']) * u.arcmin).to_value(u.rad)
        FWHM_min = (np.array(df['gau_minor_axis']) * u.arcmin).to_value(u.rad)
        pos_ang = df['gau_position_angle'] # rad
        
        for i in range(nsources):
            pix_circ = hp.query_disc(nside = self.nside, vec = vecs[i],
                                     radius = 3 * FWHM_maj[i])
            theta_pix, phi_pix = hp.pix2ang(self.nside, pix_circ)
            theta_cent, phi_cent = hp.vec2ang(vecs[i])
            
            if shape_circ == True:
                pix_dists = arclength_sphere(theta_pix, 
                                             phi_pix, 
                                             theta_cent, 
                                             phi_cent)
                profile = gaussian_source_circ(pix_dists, FWHM_maj[i])
            else:
                profile = gaussian_source((theta_pix, theta_cent),
                                          (phi_pix, phi_cent),
                                          FWHM_maj[i],
                                          FWHM_min[i],
                                          pos_ang[i])
            flux_profile = scaled_flux[i] * profile
            self.map[pix_circ] += flux_profile # Jy
            
        # Convert map from Jy -> Jy/sr with appropriate pix area
        self.map *= self.per_pix_steradian # Jy/sr
        
        if output_units == u.uK_CMB or output_units == u.K_CMB:
            self.map = (self.map * u.Jy/u.sr).to_value(output_units, equivalencies = u.cmb_equivalencies(freq_out * u.GHz))
        else:
            self.map = (self.map * u.Jy/u.sr).to_value(output_units)

        if store_maps == True:
            hp.write_map(outdir + str(freq_out) + "_GHz_GCC_Map.fits", m = self.map, coord = "G", column_units = str(output_units), overwrite = True)
            return self.map
        else:
            return self.map
