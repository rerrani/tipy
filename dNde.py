#!/usr/bin/env python3

################# Construction of tidally stripped dN/de of mass-less stellar tracers 
################# embedded in NFW dark matter halos using the procedure described
################# in E+22 Appendix G.
################# 
################# Energies e=1-E/Phi0 are referred to the potential minimum Phi0 = Phi(r=0)
################# 
################# (i)   The initial dN/de is parametrised through E+22 Eq. 13, has a 
#################       power-law slope towards the most-bound energies, and an exponential
#################       truncation beyond some scale energy
################# (ii)  We use E+22 Eq. 9 to compute the selection of particles in the 
#################       initial conditions which will form the bound remnant
################# (iii) The relaxation process is modelled using the energy mapping of
#################       E+22 Eq. 12

import numpy as np
from scipy.integrate import quad
from scipy.optimize  import fsolve
import logging
logging.basicConfig(filename='warnings.log',level=logging.DEBUG)
logging.captureWarnings(True)


# Empirical energy distribution in the initial conditions e = 1 - E/Phi0
# E+22 Eq. 13
def dNde(e,es,alpha,beta): return e**alpha * np.exp(- (e/es)**beta )

# Filter function selecting those particles in the initial dN/dE which will remain bound
# applicable for Mmx/Mmx0 < 1/10 
# E+22 Eq. 9
def f(e,emxt): return 1./(1. + (0.85 * e/emxt)**12.)  

# Energy mapping initial ei to final ef
# applicable for Mmx/Mmx0 < 1/10 
# E+22 Eq. 12
def ef(ei,emxt): return (1. + (0.8 * ei/emxt)**(-3.)  )**(-1./3.) 

# Approximation for tidal truncation energy as a function of remnant mass fraction Mmx/Mmx0
# applicable for Mmx/Mmx0 < 1/3 
# E+22 Eq. 10
def emxt_fit(MmxMmx0): return  0.77 * MmxMmx0**0.43  

# Lognormal distribution around 'xbar' with standard deviation 'dex'
def Log10normal(x,xbar,dex): return np.log10(np.exp(1.))/(dex * np.sqrt(2*np.pi)*x) * np.exp(-0.5 * ( np.log10(x / xbar)/dex )**2. )

# tidal track Vmx/Vmx0 = VV0_track(rmx/rmx0), from EN21
def VV0_track(rr0): return (2.)**0.4 * rr0**0.65 / (1. + rr0**2.)**0.4

# numerically invert tidal track , rmx/rmx0 = rr0_track_from_M(MmxMmx0)
def rr0_track_from_M(MmxMmx0): return fsolve(lambda rr0:  rr0 * VV0_track(rr0) - MmxMmx0  , 1e-3)[0]


################# ################## ################# 
################# Worked out example ################# 
################# ################## ################# 
# if dNdE is called as a program
if __name__ == "__main__":

  # Integration accuracy 
  min_ei   = 10**-5. # lower integration limit (upper limit always log10(e) = 0)
  epsabs   = 1e-30   # absolute integration accuracy of quad routine
  epsrel   = 1e-30   # relative integration accuracy of quad routine

  # Array of energies considered (1e3 in total)
  e_array = np.logspace(np.log10(min_ei),0,num=int(1e3))

  # Initial dN/de , similar to exponential profile
  # scale energy es chosen so that Rh0/rmx0 = 1/4
  # see E+22 Eq. 13, Fig. 10 and Fig. C1
  alpha = 3.
  beta  = 3.
  es    = 10**-0.32

  MmxMmx0 = 1/100.

  # Tidal truncation energy corresponding to remnant mass fraction of Mmx/Mmx0 = 1/100
  emxt = emxt_fit( MmxMmx0 ) 

  # Array of sampled initial dN/dE
  dNde_i  = dNde(e_array ,es,alpha,beta)
  norm_i  = quad( lambda this_ei: np.interp(this_ei, e_array, dNde_i)  , min_ei , 1, epsrel=epsrel, epsabs=epsabs)[0]


  # Array of the truncated dN/dE in the initial conditions 
  dNde_it =  f(e_array,emxt)  * dNde_i
  norm_it = quad( lambda this_ei: np.interp(this_ei, e_array, dNde_it)  , min_ei , 1, epsrel=epsrel, epsabs=epsabs)[0] 

  # Relaxation process
  # energies mapped through E+22 Eq. 12
  # convolution with Log-normal to account for scatter around mean energy mapping
  convolution = np.vectorize (  lambda this_ef: quad( lambda ei: f(ei,emxt)  * dNde(ei ,es,alpha,beta) *   Log10normal(this_ef,ef(ei,emxt),0.03)  , min_ei , 1, epsrel=epsrel, epsabs=epsabs)[0] )

  dNde_f = convolution(e_array)
  norm_f = quad( lambda this_ef: np.interp(this_ef, e_array, dNde_f)  , min_ei , 1, epsrel=epsrel, epsabs=epsabs)[0]


  # Using En21 tidal track to return rmx, Vmx
  LL0 = norm_it / norm_i
  rr0 = rr0_track_from_M(MmxMmx0)
  VV0 = VV0_track(rr0)


  # Output ASCII file with three columns: 
  # (a) Energy e = 1-E/Phi0 
  # (b) initial dN/de : energies refer to the initial (NFW) potential minimum
  # (c) final   dN/de : energies refer to the final (exponentially truncated cusp) minimum, related to the final mass through E+22 Eq. 11 
  np.savetxt("final_dNde.dat", np.column_stack((e_array,dNde_i/norm_i, dNde_f/norm_f )))



  print ("Tidally stripped stellar tracer in initial NFW profile")
  print (" output 'final_dNde.dat' format  : |  e = 1-E/Phi0   |   dNde_i   |   dNde_f   | ")
  print (" initial dN/de_i with e referred to Phi0 = - 4.63 * Vmx0^2 of the initial NFW")
  print (" final   dN/de_f with e referred to Phi0 = - 3.35 * Vmx^2  of the final truncated cusp")
  print (" see 'warnings.log' for (integration) warnings")
  print (" relative change in Luminosity   : L/L0     = %.4f"%LL0 )
  print (" relative change in size (DM)    : rmx/rmx0 = %.4f"%rr0 )
  print (" relative change in velocity (DM): Vmx/Vmx0 = %.4f"%VV0 )
