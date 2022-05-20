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
    bb = BlackBody(temperature = T * u.K)
    Bnu = bb(v * u.GHz) 
    Bnu0 = bb(v0 * u.GHz)
    vvb = ((v * u.GHz) / (v0 * u.GHz)) ** b
    F = (F0 * u.Jy) * (Bnu / Bnu0) * vvb
    
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

def gaussian_source_circ(pixel, radius):
    """
        Create a gaussian to represent spreading/smoothing of flux across source size
    """
    sig = (radius) / np.sqrt(8 * np.log(2))

    exp = (pixel / (np.sqrt(2) * sig)) ** 2
    
    return np.exp(-exp)

def gaussian_source(theta_set, phi_set, a, b, rot):
    phi_pix, phi_cent = phi_set
    theta_pix, theta_cent = theta_set
    siga = a / np.sqrt(8 * np.log(2))
    sigb = b / np.sqrt(8 * np.log(2))

    A = (phi_pix - phi_cent) * np.cos(rot) - (theta_pix - theta_cent) * np.sin(rot)
    B = (phi_pix - phi_cent) * np.sin(rot) + (theta_pix - theta_cent) * np.cos(rot) 
    return np.exp(-0.5 * ((A**2 / (siga**2)) + (B**2 / (sigb**2))))


def map_unit_conversion(m, freq_out, output_units):
    if output_units == u.uK_CMB or output_units == u.K_CMB:
        m = (m).to_value(output_units, equivalencies = u.cmb_equivalencies(freq_out * u.GHz))
    else:
        m = (m).to_value(output_units)
    return m