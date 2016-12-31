# ------------------------------------------------------------------------------------------ #
# Description : Additional functions for implementation of ray tracing equations by Rice (1997)    
#
# Author      : Miroslav Mocak (Slovak Organization For Space Activities)
# Date        : 31/August/2016
# References  : Rice W.K.M, 1997, "A ray tracing study of VLF phenomena", PhD thesis, 
#             : Space Physics Research Institute, Department of Physics, University of Natal
#             : Yabroff (1961), Kimura (1966)  
#             : INTERNATIONAL REFERENCE IONOSPHERE IRI-2012 http://omniweb.gsfc.nasa.gov/vitmo/iri2012_vitmo.html  
# ------------------------------------------------------------------------------------------ #

import numpy as np
import matplotlib.pyplot as plt
import ray_cmks
from scipy.interpolate import interp1d
	
# ------------------------------------ #
# DEFINE FUNCTIONS TO BE USED IN MAIN 
# ------------------------------------ #
	
def phase_refractive_index(rr,th,chi,freq,ion,ionosphere):
    height,ne,O_ions_per,H_ions_per,He_ions_per,O2_ions_per,NO_ions_per,N_ions_per = ionosphere    # unpack 
    ray_cmks.pconstants()	
	
    dip = np.arctan(2.*np.tan(np.pi/2.-th))
    phi = (3./2.)*np.pi - dip
    psi = chi - phi
	
    # angular frequency
    omega = 2.*np.pi*freq	
	
    rrE = rr/ray_cmks.Re

	# determines the magnetic field strength at r,theta using dipole field model	
    B0 = 3.0696381e-5 # in Tesla 	
    Bmag = B0*((ray_cmks.Re**3)/(rr**3))*np.sqrt(4.-3.*np.cos((np.pi/2.)-th)*np.cos((np.pi/2.)-th))
	
    # calculate free electron/proton density stratification Yabroff (1961)	
    if ion[0] == "iono_np":
	
        ne = 1.e6*(1.8e5*np.exp(-4.183119*(rrE-1.0471)))  
        nprot = 1.e6*(1.8e5*np.exp(-4.183119*(rrE-1.0471))) 
	
		# electron plasma angular frequency squared
        pie2 = (ne*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.me)
		# proton plasma angular frequency squared
        pip2 = (nprot*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mp)
	
		# electron angular gyrofrequency
        omegae = (ray_cmks.e*Bmag)/(ray_cmks.me)
		# proton angular gyrofrequency
        omegap = (ray_cmks.e*Bmag)/(ray_cmks.mp)
  
        piO2 = 0.
        piHe2 = 0.
        piO22 = 0.
        piNO2 = 0. 
        piN2 = 0.
		
        omegaO = 0.
        omegaHe = 0.
        omegaO2 = 0.
        omegaNO = 0.
        omegaN = 0. 		
  
    # calculate free electron density stratification Yabroff (1961)  
    if ion[0] == "iono_el":

        ne = 1.e6*(1.8e5*np.exp(-4.183119*(rrE-1.0471)))   
	
		#   print('rr,ne,Bmag',rr,ne,Bmag)
	
		# electron plasma angular frequency squared
        pie2 = (ne*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.me)
		# proton plasma angular frequency squared
        pip2 = 0.
	
		# electron angular gyrofrequency
        omegae = (ray_cmks.e*Bmag)/(ray_cmks.me)
		# proton angular gyrofrequency
        omegap = 0.	
		
        pip2 = 0.
        piO2 = 0.
        piHe2 = 0.
        piO22 = 0.
        piNO2 = 0. 
        piN2 = 0.
		
        omegap = 0.
        omegaO = 0.
        omegaHe = 0.
        omegaO2 = 0.
        omegaNO = 0.
        omegaN = 0.

    # read/intepolate IRI model  
    if ion[0] == "iono_ir":	
	
        height = ray_cmks.Re + np.asarray(height)*1000.
        ne = np.asarray(ne)		
        # convert percentage to density 
        O_ions  = np.asarray(O_ions_per)*np.asarray(ne)/100.
        H_ions  = np.asarray(H_ions_per)*np.asarray(ne)/100.
        He_ions = np.asarray(He_ions_per)*np.asarray(ne)/100.
        O2_ions = np.asarray(O2_ions_per)*np.asarray(ne)/100.
        NO_ions = np.asarray(NO_ions_per)*np.asarray(ne)/100.
        N_ions  = np.asarray(N_ions_per)*np.asarray(ne)/100.
		

        if (height[0] < rr < height[-1]):
          ne_int = np.interp(rr,height,ne)
          nH_ions_int = np.interp(rr,height,H_ions)
          nO_ions_int = np.interp(rr,height,O_ions)
          nHe_ions_int = np.interp(rr,height,He_ions)
          nO2_ions_int = np.interp(rr,height,O2_ions)
          nNO_ions_int = np.interp(rr,height,NO_ions)
          nN_ions_int = np.interp(rr,height,N_ions)
        else:
          ne_int = 0.
          nH_ions_int = 0.
          nO_ions_int = 0.
          nHe_ions_int = 0.
          nO2_ions_int = 0.
          nNO_ions_int = 0.
          nN_ions_int = 0.		
        
        if (rr > height[-1]):
          ne_int = ne[-1]
          nH_ions_int = H_ions[-1]
          nO_ions_int = O_ions[-1]
          nHe_ions_int = He_ions[-1]
          nO2_ions_int = O2_ions[-1]
          nNO_ions_int = NO_ions[-1]
          nN_ions_int = N_ions[-1]

        if (rr < height[0]):
          ne_int = 0.
          nH_ions_int = 0.
          nO_ions_int = 0.
          nHe_ions_int = 0.
          nO2_ions_int = 0.
          nNO_ions_int = 0.
          nN_ions_int = 0.		  
		  
#        print(rr,height[0],height[-1])
		  
		# electron plasma angular frequency squared
        pie2 = (ne_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.me)
		# proton plasma angular frequency squared
        pip2 = (nH_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mp)
		# atomic oxygen plasma angular frequency squared
        piO2 = (nO_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mO)
		# helium plasma angular frequency squared
        piHe2 = (nHe_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mHe)
		# molecular oxygen plasma angular frequency squared
        piO22 = (nO2_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mO2)
		# nitric oxide plasma angular frequency squared
        piNO2 = (nNO_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mNO)	
		# nitrogen plasma angular frequency squared
        piN2 = (nN_ions_int*(ray_cmks.e**2))/(ray_cmks.eps*ray_cmks.mN)				
		
		# electron angular gyrofrequency
        omegae = (ray_cmks.e*Bmag)/(ray_cmks.me)
		# proton angular gyrofrequency
        omegap = (ray_cmks.e*Bmag)/(ray_cmks.mp)
		# oxygen angular gyrofrequency
        omegaO = (ray_cmks.e*Bmag)/(ray_cmks.mO)
		# helium angular gyrofrequency
        omegaHe = (ray_cmks.e*Bmag)/(ray_cmks.mHe)
		# molecular oxygen angular gyrofrequency
        omegaO2 = (ray_cmks.e*Bmag)/(ray_cmks.mO2)
		# proton angular gyrofrequency
        omegaNO = (ray_cmks.e*Bmag)/(ray_cmks.mNO)
		# atomic nitrogen angular gyrofrequency
        omegaN = (ray_cmks.e*Bmag)/(ray_cmks.mN)
  
    R = 1.- (pie2/omega)*(1./(omega-omegae)) \
    - (pip2/omega)*(1./(omega+omegap)) \
    - (piO2/omega)*(1./(omega+omegaO)) \
    - (piHe2/omega)*(1./(omega+omegaHe)) \
    - (piO22/omega)*(1./(omega+omegaO2)) \
    - (piNO2/omega)*(1./(omega+omegaNO)) \
    - (piN2/omega)*(1./(omega+omegaN))
		  
    L = 1.- (pie2/omega)*(1./(omega+omegae)) \
    - (pip2/omega)*(1./(omega-omegap)) \
    - (piO2/omega)*(1./(omega-omegaO)) \
    - (piHe2/omega)*(1./(omega-omegaHe)) \
    - (piO22/omega)*(1./(omega-omegaO2)) \
    - (piNO2/omega)*(1./(omega-omegaNO)) \
    - (piN2/omega)*(1./(omega-omegaN))    		  
	
    P = 1.- (pie2/(omega**2)) \
    - (pip2/(omega**2)) \
    - (piO2/(omega**2)) \
    - (piHe2/(omega**2)) \
    - (piO22/(omega**2)) \
    - (piNO2/(omega**2)) \
    - (piN2/(omega**2))

#    R = 1.-(pie2/omega)*(1./(omega-omegae)) - (pip2/omega)*(1./(omega+omegap)) 
#    L = 1.-(pie2/omega)*(1./(omega+omegae)) - (pip2/omega)*(1./(omega-omegap)) 
#    P = 1.-(pie2/(omega**2))-(pip2/(omega**2))
	
    D = (1./2.)*(R-L)
    S = (1./2.)*(R+L)
    
    A  = S*np.sin(psi)*np.sin(psi) + P*np.cos(psi)*np.cos(psi)
    B  = R*L*np.sin(psi)*np.sin(psi) + P*S*(1.+np.cos(psi)*np.cos(psi))
    C  = P*R*L
    F2 = ((R*L-P*S)*(R*L-P*S))*(np.sin(psi)*np.sin(psi)*np.sin(psi)*np.sin(psi)) + 4.*(P**2)*(D**2)*(np.cos(psi)*np.cos(psi)) 
    F = np.sqrt(F2)  
				
    n2m = (B - F)/(2.*A)
    n_minus  = np.sqrt(n2m)
    n2p =(B + F)/(2.*A)
    n_plus  = np.sqrt(n2p)

#    print("n_plus","n_minus",n_plus,n_minus)
	
    n = n_minus
#    n = n_plus
	
    # dndpsi
    dAdpsi = 2.*(S-P)*np.sin(psi)*np.cos(psi)
    dBdpsi = 2.*(R*L-P*S)*np.sin(psi)*np.cos(psi)
    dCdpsi = 0.
    dndpsi = ((n**4)*dAdpsi-(n**2)*dBdpsi+dCdpsi)/(4.*A*(n**3)-2.*B*n)
	
#    print(n,dndpsi)
#    stop
	
    return n,dndpsi

def deriv_rr(rr,th,chi,freq,ion,ionosphere):
    drr = 1.e-11

    n_l = phase_refractive_index(rr-drr/2.,th,chi,freq,ion,ionosphere)
    n_r = phase_refractive_index(rr+drr/2.,th,chi,freq,ion,ionosphere)
	
    derivndr = (n_r[0] - n_l[0])/drr

    return derivndr
	
def deriv_th(rr,th,chi,freq,ion,ionosphere):
    dth = 1.e-11

    n_l = phase_refractive_index(rr,th-dth/2.,chi,freq,ion,ionosphere)
    n_r = phase_refractive_index(rr,th+dth/2.,chi,freq,ion,ionosphere)
	
    derivndth = (n_r[0] - n_l[0])/dth
	
    return derivndth

def deriv_ch(rr,th,chi,freq,ion,ionosphere):	
    dch = 1.e-11

    n_l = phase_refractive_index(rr,th,chi-dch/2.,freq,ion,ionosphere)
    n_r = phase_refractive_index(rr,th,chi+dch/2.,freq,ion,ionosphere)
	
    derivndch = (n_r[0] - n_l[0])/dch
	
    return derivndch
	
def deriv_fr(rr,th,chi,freq,ion,ionosphere):	
    df = 1.e-5

    n_l = phase_refractive_index(rr,th,chi,freq-df/2.,ion,ionosphere)
    n_r = phase_refractive_index(rr,th,chi,freq+df/2.,ion,ionosphere)
	
    derivndf = (n_r[0] - n_l[0])/df
	
    return derivndf




