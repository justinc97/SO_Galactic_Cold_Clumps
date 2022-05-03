import numpy as np
import healpy as hp
import pysm3.units as u
from astropy.modeling.models import BlackBody
import warnings
warnings.filterwarnings('ignore')

def modBB(F0, T, v, v0, b):
    """
        Modified black body equation to scale flux density measured at reference 
        frequency to a new frequency using measured spectral properties.
        -----Inputs-----
        F0   : Flux density fit at reference wavelength v0
        T    : Colour temperature
        v    : Frequency out
        v0   : Reference frequency in
        b    : emissivity spectral index
        -----Outputs-----
        F    : Flux scaled to frequency out
    """
    per_pixel_steradian = 1 / hp.nside2pixarea(2048) # From Planck map sizes
    bb = BlackBody(temperature = T * u.K)
    Bnu = bb(v) 
    Bnu0 = bb(v0)
    vvb = (int(v) / int(v0)) ** b
    F = F0 * (Bnu / Bnu0) * vvb
    
    return F.value

def arclength_sphere(theta, phi, theta_cent, phi_cent):
    """
        Compute arclength of a given pixel from the centre of a unitary sphere. 
        Angle between these points obtained with the dot product.
    """
    # Center position
    ccos = np.cos(theta_cent)
    csin = np.sin(theta_cent)
    # Pixel position
    pcos = np.cos(theta)
    psin = np.sin(theta)
    #Cos(center phi - pixel phi) -> For dot product
    cosphi = np.cos(phi_cent - phi)
    
    psi = np.arccos(ccos * pcos + csin * psin * cosphi)
    
    nanmask = np.isnan(psi) # Account for arccos ~>1 = NaN instead of 0
    psi[nanmask] = 0
    
    return psi

#def gaussian_source(pixel, fwhm_maj, fwhm_min):
    """
        Create a gaussian to represent spreading/smoothing of flux across source size
    """
    #sig_maj = (fwhm_maj) / np.sqrt(8 * np.log(2))
    #sig_min = (fwhm_min) / np.sqrt(8 * np.log(2))

    #exp = (pixel / (np.sqrt(2) * sig_maj)) ** 2 + (pixel / (np.sqrt(2) * sig_min)) ** 2
    
    #return np.exp(-exp)

def gaussian_source(theta_set, phi_set, a, b, rot):
    theta_pix, theta_cent = theta_set
    phi_pix, phi_cent = phi_set
    siga = a / np.sqrt(8 * np.log(2))
    sigb = b / np.sqrt(8 * np.log(2))
    
    """
    x = ((np.cos(rot)**2) / (2 * siga**2)) + ((np.sin(rot)**2) / (2 * sigb**2))
    y = -(np.sin(2 * rot)) / (4 * siga**2)  + (np.sin(2 * rot)) / (4 * sigb**2)
    z = (np.sin(rot)**2) / (2 * siga**2) + (np.cos(rot)**2) / (2 * sigb**2)
    g = np.exp(-(x * ((theta_pix - theta_cent)**2) + 2 * y * (theta_pix - theta_cent) * (phi_pix - phi_cent) + z * ((phi_pix - phi_cent)**2)))
    return g
  """
    A = ((theta_pix - theta_cent) * np.cos(rot) - (phi_pix - phi_cent) * np.sin(rot)) / (siga)
    B = ((theta_pix - theta_cent) * np.sin(rot) + (phi_pix - phi_cent) * np.cos(rot)) / (sigb)
    return np.exp(-0.5 * (A**2 + B**2))
"""
def gaussian_source(theta_set, phi_set, a, b, pa):
    theta_pix, theta_c = theta_set
    phi_pix, phi_c = phi_set
    sigma_maj = a / np.sqrt(8 * np.log(2))
    sigma_min = b / np.sqrt(8 * np.log(2))
    
    A = ((theta_pix - theta_c)**2) / (2 * sigma_maj**2)
    B = (pa * (theta_pix - theta_c) * (phi_pix - phi_c)) / (sigma_maj * sigma_min)
    C = ((phi_pix - phi_c)**2) / (2 * sigma_min**2)
    G = np.exp(-A-B-C)
    return G

def gaussian_source(theta_set, phi_set, a, b, pa):
    theta_pix, theta_c = theta_set
    phi_pix, phi_c = phi_set
    sigma_maj = a / np.sqrt(8 * np.log(2))
    sigma_min = b / np.sqrt(8 * np.log(2))
    
    A = (theta_pix - theta_c) * np.cos(pa) + (phi_pix - phi_c) * np.cos(pa)
    B = (theta_pix - theta_c) * np.sin(pa) + (phi_pix - phi_c) * np.cos(pa) 
    G = (A / sigma_maj) ** 2 + (B / sigma_min) ** 2
    return np.exp(-0.5 * G)
    
    """