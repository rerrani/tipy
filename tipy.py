#!/usr/bin/env python3
#
# Example implementation of the tidal stripping model introduced in Errani & Navarro 2020 
# 
# Two ways to use this file:  
#
#   i)  load it as a module: from tipy import *
#       you then have access to the functions
#       fecc_fit(rarp)         eccentricity delay factor fecc     = fecc_fit(rapo/rperi)     (paper Eq. 4)
#       VV0_track(rr0)         tidal track               Vmx/Vmx0 = VV0_track(rmx/rmx0)      (paper Eq. 5)
#       rr0_from_track(TT0)    tidal track               rmx/rmx0 = rr0_from_track(Tmx/Tmx0)
#
#   ii) run it as a program: python3 tipy.py
#       this will present a worked out example for the evolution of a subhalo


from math import pi, exp
import numpy as np
from scipy.optimize import fsolve



# tidal track Vmx/Vmx0 = VV0_track(rmx/rmx0)
def VV0_track(rr0): return (2.)**0.4 * rr0**0.65 / (1. + rr0**2.)**0.4

# numerically invert tidal track , rmx/rmx0 = rr0_from_track(Tmx/Tmx0)
def rr0_from_track(TT0): return fsolve(lambda rr0:  rr0/VV0_track(rr0) - TT0  , 1e-3)[0]


# eccentricity delay factor fecc = (2x/(x+1))^c with x = rapo/rperi
def fecc_fit(rarp): return ( 2*rarp / (rarp+1.) )**3.2



# Evolution of Tmx as a function of initial Tmx0 and time t 
# Note: here Tmx is in units of Tperi ,  and t is in units of Torb
def get_Tmx(Tmx0, t):  

  # Regime I ( "Heavy mass loss" ) 
  if Tmx0 > 0.66:  
    Tasy    = 0.22    # in units of Tperi    
    tau_asy = 0.65    # in units of Torb
    Y0      = Tmx0 - Tasy
    tau     = tau_asy / Y0                                   # Paper Eq. 13
    eta   = 1 - exp(-2.5 * Y0)                               # Paper Eq. 14
    Y       = Y0 * (1. + (t/tau)**eta )**(-1./eta)           # Paper Eq. 12
    Tmx     = Tasy + Y

  # Regime II ( "Moderate mass loss" ) 
  else:            
    Tasy  = Tmx0 / (1.+ Tmx0)**2.17   # in units of Tperi    # Paper Eq. 16
    Y0    = Tmx0 - Tasy
    tau   = 1.2 * Tmx0**-0.5          # in units of Torb     # Paper Eq. 17
    eta = 0.67
    Y     = Y0 * (1. + (t/tau)**eta )**(-1./eta)             # Paper Eq. 15
    Tmx   = Tasy + Y
    
  return Tmx     # in units of Tperi




# worked out example, if tipy is called as a program
if __name__ == "__main__":
  ### simulation units
  G = 1.
  mU= 1.0e10   # Msol
  rU= 1.0      # kpc
  tU= 4.714e-3 # Gyrs
  vU= 207.4    # km/s
  
  ########## PARAMETERS ##########
  ### HOST + ORBIT
  rperi     = 20. / rU                       # kpc -> sim
  raporperi = 5.                             # rapo / rperi  
  Torb      = 1.23 / tU                      # Gyrs -> sim
  Tperi     = 2 *pi * rperi / (220./vU)      # vc=220 km/s at rperi
  fecc      = fecc_fit(raporperi)
  
  ### SUBHALO
  rmx0      = 2.2/rU                         # kpc  -> sim
  Vmx0      = 13.95/vU                       # km/s -> sim
  Tmx0      = 2 *pi * rmx0  / Vmx0
  
  
  ##########   time loop    ##########
  for t in np.linspace(0,15,num=16) / tU :   # 15 Gyrs
    
    Tmx = Tperi * get_Tmx(Tmx0/Tperi, t/Torb/fecc)
    rmx = rmx0  * rr0_from_track(Tmx/Tmx0)
    Vmx = Vmx0  * VV0_track(rmx/rmx0)
    Mmx = rmx   * Vmx**2. / G
    
    print (t*tU, Tmx*tU, rmx*rU, Vmx*vU, Mmx*mU)

  
  
