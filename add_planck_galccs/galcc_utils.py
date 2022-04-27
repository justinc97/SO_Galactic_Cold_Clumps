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

def gaussian_source(pixel, size):
    """
        Create a gaussian to represent spreading/smoothing of flux across source size
    """
    sig = (size) / np.sqrt(8 * np.log(2))
    exp = (pixel / (np.sqrt(2) * sig)) ** 2
    
    return np.exp(-exp)