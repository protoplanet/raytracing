# Implementation of ray tracing equations by Haselgrove (1957)  
# Ray tracing equations reformulated using delta by Shawhan (1966)

# Author: Miroslav Mocak
# Date: 21/July/2016
# Slovak Organization For Space Activities
# ODE solver syntax based on http://www.physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html

import numpy as np
import matplotlib.pyplot as plt
import ccgs
#from scipy.integrate import odeint
from scipy.integrate import ode
	
# ------------------------------------ #
# DEFINE FUNCTIONS TO BE USED IN MAIN 
# ------------------------------------ #
	
def phase_refractive_index(params,x,z,psi):
    Xm,Y,Bmag,freq,phi = params  # unpack parameters 
		
    ccgs.pconstants()	
		
    # angular frequency
    omega = 2.*np.pi*freq	
	
    freq = 3.e6
    Xm = 0.5
	
    # calculate free electron density stratification
    Zt = 1.e7 # half width in cm
    ne = ((np.pi*ccgs.me*(freq**2)*Xm)/(ccgs.e**2))*(1.-((1.-(z/Zt))**2))  

    if z > 200.e5:
        ne = ((np.pi*ccgs.me*(freq**2)*Xm)/(ccgs.e**2))*(1.-((1.-(2*Zt/Zt))**2)) 
	
    # electron plasma angular frequency squared
    pie2 = (4.*np.pi*ne*(ccgs.e**2))/(ccgs.me)
	
    # electron angular gyrofrequency
    omegae = (ccgs.e*Bmag)/(ccgs.me*ccgs.c)
  
    R = 1.-(pie2/omega)*(1./(omega-omegae))
    L = 1.-(pie2/omega)*(1./(omega+omegae))
    P = 1.-(pie2/(omega*omega))
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
		
    n = n_minus
    	
    return n
	
def f(t, xzd, freq):
    x, z, delta = xzd    # unpack current values 
    Xm,Y,Bmag,freq,phi = params  # unpack parameters  
	
    psi = (np.pi/2.) - delta - phi
	
    n = phase_refractive_index(params,x,z,psi)
 
#    zzs.append(z)   
#    nns.append(n)
	
#    print('z,n',z,n)
	
    dndx   = deriv_x(n,x,z,psi)
    dndz   = deriv_z(n,x,z,psi)
    dndpsi = deriv_p(n,x,z,psi)
	
    derivs = [(1./n**2)*(n*np.sin(delta) - dndpsi*np.cos(delta)),      
              (1./n**2)*(n*np.cos(delta) - dndpsi*np.sin(delta)),
	          (1./n**2)*(-dndz*np.sin(delta) + dndx*np.cos(delta))]
    return derivs

def deriv_x(n,x,z,psi):
    dx = 1.e-5

    n_l = phase_refractive_index(params,x-dx/2.,z,psi)
    n_r = phase_refractive_index(params,x+dx/2.,z,psi)
	
    derivndx = (n_r - n_l)/dx

    return derivndx
	
def deriv_z(n,x,z,psi):
    dz = 1.e-5

    n_l = phase_refractive_index(params,x,z-dz/2.,psi)
    n_r = phase_refractive_index(params,x,z+dz/2.,psi)
	
    derivndz = (n_r - n_l)/dz
	
    return derivndz

def deriv_p(n,x,z,psi):	
    dp = 1.e-5

    n_l = phase_refractive_index(params,x,z,psi-dp/2.)
    n_r = phase_refractive_index(params,x,z,psi+dp/2.)
	
    derivndp = (n_r - n_l)/dp
	
    return derivndp
	
def SetMatplotlibParams():
#   this routine sets some standard values for matplotlib to obtain publication-quality figures
	
#	plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('font',**{'family':'serif','serif':['Times New Roman']})
    plt.rc('font',size=22.)
    plt.rc('lines',linewidth=2,markeredgewidth=2.,markersize=10)
    plt.rc('axes',linewidth=1.5)
    plt.rcParams['xtick.major.size']=8.
    plt.rcParams['xtick.minor.size']=4.
    plt.rcParams['figure.subplot.bottom']=0.13
    plt.rcParams['figure.subplot.left']=0.17	
	
def ne(params):
    Xm,Y,Bmag,freq,phi = params  # unpack parameters 
		
    ccgs.pconstants()	
	
    aStop = 200.e5
    aInc = 1.e5
    alt = np.arange(0., aStop, aInc)
    freq = 3.e6
    Xm = 0.5	
    # calculate free electron density stratification
    Zt = 1.e7 # half width in cm
    ne = ((np.pi*ccgs.me*(freq**2)*Xm)/(ccgs.e**2))*(1.-((1.-(alt/Zt))**2)) 
	
    return ne,alt
	
# ------------------------------------ #
# MAIN 
# ------------------------------------ #
	
# load constants
ccgs.pconstants()
	
# input parameters 
Xm = 0.5
Bmag  = 0.5 # in Gauss	
freq  = 3.e6 # in Hz
dgr_phi   = 0.   # magnetic field inclination
phi = dgr_phi*np.pi/180.

fH = ((Bmag*ccgs.e)/(ccgs.me*ccgs.c))/(2.*np.pi)
Y = fH/freq

# initial conditions
x0 = 0.
z0 = 0.
dgr_delta0 = 0.
delta0 = dgr_delta0*np.pi/180.
	
# bundle initial conditions 
#xzd0 = [x0,z0,delta0] 

# bundle parameters for ODE solver
params = [Xm,Y,Bmag,freq,phi]

# make time array for the solver
t0 = 0. 
tStop = 2.e9 # maximum time for integration
nstep  = 1.e7 # number of integration steps
dt = tStop/nstep


# ------------------------------------ #
# PLOT RESULTS
# ------------------------------------ #

# set parameters for plotting
SetMatplotlibParams()

#xmin = 0.   # min x boundary for plotting
#xmax = 200. # max x boundary for plotting
#ymin = 0.   # min y boundary for plotting 
#ymax = 200.  # max y boundary for plotting

#plt.figure(1,figsize=(14,7))
#plt.axis([xmin,xmax,ymin,ymax])
	
xmin2 = 0.   # min x boundary for plotting
xmax2 = 200. # max x boundary for plotting
ymin2 = 0.   # min y boundary for plotting 
ymax2 = 200.  # max y boundary for plotting	
	
plt.figure(3,figsize=(14,7))
plt.axis([xmin2,xmax2,ymin2,ymax2])

angle = 5.*np.pi/180.
#angle = 11.25*np.pi/180.
for ii in range(1,15):
    delta0 = ii*angle # vary incident angle 
    print(str(delta0*180./np.pi))
	
    # bundle initial conditions 
    xzd0 = [x0,z0,delta0] 

    # call the ODE solver
    psoln = ode(f).set_integrator('lsoda',method='bdf')
    #psoln = ode(f).set_integrator('vode',method='bdf')
    psoln.set_initial_value(xzd0,t0).set_f_params(freq)
	
    xx = []
    zz = []
    nn = []
	
    while psoln.successful() and psoln.t < tStop and psoln.y[1] >= 0. and psoln.y[1] < 199.e5:
	    psoln.integrate(psoln.t+dt)
	    xx.append(psoln.y[0])
	    zz.append(psoln.y[1])
	    psi = (np.pi/2.) - psoln.y[2] - phi
	    ntmp = phase_refractive_index(params,psoln.y[0],psoln.y[1],psi)   
	    nn.append(ntmp)
#	    print(psoln.y[0],psoln.y[1],psoln.y[2])

    xx = np.asarray(xx)
    zz = np.asarray(zz)
    nn = np.asarray(nn)    

    dgr_delta0 = delta0*(180./np.pi)
	
    plt.text(120.,160.,r"Xm = "+str(Xm))
    plt.text(120.,140.,r"Y  = "+str(round(Y,1)))
    plt.text(120.,120.,r"f = "+str(round(freq,1))+" Hz")
    plt.text(120.,100.,r"B = "+str(Bmag)+" Gauss")
    plt.text(120.,80.,r"$\Phi$ = "+str(dgr_phi))	
		
    plt.plot(xx/1.e5,zz/1.e5,label=r"$\Delta$ = i = "+str(np.round(180.*delta0/np.pi,1)))
	
    plt.xlabel('x (km)')
    plt.ylabel('z (km)')
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'Haselgrove(1957) Shawhan(1966)')


#    plt.text(0.05,160.,r"Xm = "+str(Xm)+r" = $f_c^2/f^2 (f_c$ is critical freq)")
#    plt.text(0.05,140.,r"Y  = "+str(round(Y,1))+r" = $f_H/f (f_H$ is gyro freq)")
#    plt.text(0.05,120.,r"f = "+str(round(freq,1))+" Hz")
#    plt.text(0.05,100.,r"B = "+str(Bmag)+" Gauss")
#    plt.text(0.05,80.,r"$\Phi$ = "+str(dgr_phi))

#    plt.plot(nn,zz/1.e5,label=r"$\Delta$ = i = "+str(np.round(180.*delta0/np.pi,1)))
#    plt.xlabel(r"n")
#    plt.ylabel('z (km)')
#    plt.legend(loc=1,prop={'size':14})
#    plt.title(r'refractive index Haselgrove (1957)')	

name='rayhasel2e_phi'+str(dgr_phi)+'_freq'+str(round(freq,1))

#name='rayhasel2e_phi'+str(dgr_phi)+'_refractive_index_freq'+str(round(freq,1))

plt.savefig(name+'.png')
	
eled = ne(params)
a = eled[0]
b = eled[1]
	
xmin1 = 0.   # min x boundary for plotting
xmax1 = 120000. # max x boundary for plotting
ymin1 = 0.   # min y boundary for plotting 
ymax1 = 200.  # max y boundary for plotting

plt.figure(2,figsize=(8,7))
plt.axis([xmin1,xmax1,ymin1,ymax1])
plt.plot(eled[0],eled[1]/1.e5)
plt.xlabel(r"n$_e$ (cm$^{-3}$)")
plt.ylabel('z (km)')
plt.legend(loc=1,prop={'size':14})
plt.title(r'ne Haselgrove (1957)')
	
name='rayhasel2e_ne_'+str(round(freq,1))

plt.savefig(name+'.png')
	
plt.show()	

# ------------------------------------ #
# END
# ------------------------------------ #
