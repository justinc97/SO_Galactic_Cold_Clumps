# Credit Giuseppe Puglisi

import pylab as pl 
import scipy as sp 
from astropy.modeling import models 
from astropy import constants as const, units as u 
from scipy import integrate as integr
import numpy as np 
BB  = models.BlackBody(temperature=2.725 *u.K    ) 

def two_powerlaws(k ,k0=1e-1 , gamma1= -1, gamma2= -3 ) : 
    A2 =  k0 **(gamma1 - gamma2 )
    if k<=k0 :
        return  k **gamma1 
    else: 
        return A2 *k **gamma2 
    
def exp(k  ) : 
    return k/ (pl.exp(k) +1) 


def build_sampler_functional( ydata , xdata , debug =False  ): 
    
    
    #normalize it to make it a PDF 
    A=pl.trapz( ydata , xdata)
    
    #Interpolate the  function 
    #pdf  = sp.interpolate.interp1d(x=xdata, y=ydata /A ) 
    pdf = ydata /A
    #estimate the CDF 
    dx = pl.diff(xdata)
    cumulative   = pl.array([ pl.trapz(   pdf [ :i]  , xdata[:i]   ) for i in range(1,pdf.shape[0])] ) 
    xv =  xdata[1: ] 
    #interpolate the inverse CDF
    inverse_cdf = sp.interpolate.interp1d(cumulative  ,xdata[:-1] +dx , 
                                         bounds_error=False, fill_value="extrapolate" )   
    
    if debug : 
        print(f"normalization:{A:.3f}")
        fig = pl .figure(figsize=(10,10))
        pl.subplot(311)
        pl.loglog(xdata, pdf,  )
        pl.ylabel(r'$p(x)dx $', fontsize=15)
        pl.xlabel(r'$x$', fontsize=15)
        
        pl.yticks(fontsize=15)
        pl.xticks(fontsize=15)
        pl.subplot(312)
        pl.loglog(xv, cumulative )
        pl.ylabel(r'$F(>x) $', fontsize=15)
        pl.xlabel(r'$x$', fontsize=15)
        pl.yticks(fontsize=15)
        pl.xticks(fontsize=15)
        pl.subplot(313)
        pl.semilogy( cumulative  ,xdata[:-1] +dx)
        pl.ylim(xdata.min() , xdata.max() ) 
        pl.ylabel(r'$F^{-1}(u) $', fontsize=15)
        pl.xlabel(r'$u $', fontsize=15)
        pl.yticks(fontsize=15)
        pl.xticks(fontsize=15)
        pl.tight_layout() 
        

    return inverse_cdf 

def hist2sample(xdata, ydata, bstrap_size, debug = False):
    BB = models.BlackBody(temperature = 2.725 * u.K)
    bin_cents = 0.5*(xdata[1:] + xdata[:-1])
    pdf, cdf, sampler = build_sampler_functional(ydata = ydata, xdata = bin_cents, debug = debug)
    norm = pl.trapz(ydata, bin_cents)
    for i in range(10): 
        z = np.random.uniform(size = bstrap_size)
        bootstrapped_values = sampler(z)
        bins = pl.logspace(np.log10(bin_cents.min()), np.log10(bin_cents.max()), 30)
        h, edg = np.histogram(bootstrapped_values, bins = bins, density=True) 
        xb = [(edg[i] + edg[i+1])*.5 for i in range(edg.shape[0]-1)]
        pl.plot (xb, h*norm, 'k.', alpha=.2)
    pl.plot(bin_cents,ydata)
    pl.loglog
    pl.show()

    return pdf, cdf, [xb, h*norm]