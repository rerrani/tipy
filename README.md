# tipy.py 

Empirical model for the tidal evolution of subhalos as described in Errani & Navarro 2020

Two ways to use this file:


   *  load it as a module: ` from tipy import * `
   
       you then have access to the functions
       
       ` fecc_fit(rarp) `      eccentricity delay factor  ` fecc     = fecc_fit(rapo/rperi) `    (analytical, paper Eq. 4)
       
       ` VV0_track(rr0) `       tidal track              ` Vmx/Vmx0 = VV0_track(rmx/rmx0) `     (analytical, paper Eq. 5)
       
       ` rr0_from_track(TT0) `   tidal track              ` rmx/rmx0 = rr0_from_track(Tmx/Tmx0) ` (numerically solved)
       
       ` get_Tmx(Tmx0, t) ` time evolution of crossing time ` Tmx/Tperi = get_Tmx(Tmx0/Tperi, t/Torb)  ` (analytical, paper Eq. 12 + 15 )

   * run it as a program: ` python3 tipy.py `
   
       this will run a worked out example for the evolution of a subhalo





# dNde.py 

Empirical model to construct energy distribution of a tidally stripped stellar tracer as in E+22 Appendix G

Also here, two ways to use this file:

   *  load it as a module: ` from dNdE import * `
   
       this gives access to the functions
       
       ` dNde(e,es,alpha,beta) `    Empirical energy distribution in the initial conditions      (paper Eq. 13)
       
       ` f(e,emxt) `                Filter function, i.e., shape of the tidal energy truncation  (paper Eq. 9)
       
       ` ef(ei,emxt) `              Energy mapping initial ei to final ef                         (paper Eq. 12)
       
       ` emxt_fit(MmxMmx0) `        Tidal truncation energy vs. remnant mass                      (paper Eq. 10)
       
       ` VV0_track(rr0) `  tidal track              ` Vmx/Vmx0 = VV0_track(rmx/rmx0) `          (EN21 tital track)
       
       ` rr0_track_from_M(MmxMmx0) ` tidal track     ` rmx/rmx0 = rr0_track_from_M(MmxMmx0) `   (numerically solved)

   * run it as a program: ` python3 dNdE.py `
   
       this will run a worked out example for a stellar tracer embedded in an NFW that has been stripped to one per cent of its initial mass

