# Implementation of ray tracing equations of Haselgrove (1957) 

# Author: Miroslav Mocak
# Date: 18/July/2016
# ODE solver syntax based on http://www.physics.nyu.edu/pine/pymanual/html/chap9/chap9_scipy.html

import numpy as np
import matplotlib.pyplot as plt
import ccgs
from scipy.integrate import odeint

def deriv_x(n,x,y,psi):
    dx = 1.e-5

    n_l = phase_refractive_index(params,x-dx/2.,y,psi)
    n_r = phase_refractive_index(params,x+dx/2.,y,psi)
	
    derivndx = (n_r - n_l)/dx

    return derivndx
	
def deriv_y(n,x,y,psi):
    dy = 1.e-5

    n_l = phase_refractive_index(params,x,y-dy/2.,psi)
    n_r = phase_refractive_index(params,x,y+dy/2.,psi)
	
    derivndy = (n_r - n_l)/dy
	
    return derivndy

def deriv_c(n,x,y,psi):	
    dp = 1.e-5

    n_l = phase_refractive_index(params,x,y,psi-dp/2.)
    n_r = phase_refractive_index(params,x,y,psi+dp/2.)
	
    derivndc = (n_r - n_l)/dp
	
    return derivndc
	
def phase_refractive_index(params,x,y,psi):
    Y, Xm, phi = params  # unpack parameters 

    ht = 1.e7 # 100 km in cm 
	
    X = Xm*(1.-((1.-y/ht)**2))
	
    YT = -Y*np.sin(psi)
    YL = -Y*np.cos(psi)
	
    A = (YT**2)/(2.*(1.-X))
    B = (YT**4)/(4.*((1.-X)**2))+YL**2
	
    n2p = 1.-(X/((1.-A)+np.sqrt(B)))
    n2m = 1.-(X/((1.-A)-np.sqrt(B)))
	
    n = np.sqrt(n2p)
	
#    target.write(str(n))
#    target.write(" ")
#    target.write(str(y))
#    target.write("\n")
	
    return n

def f(xyi, t, params):
    x, y, chi = xyi     # unpack current values 
    Y, Xm, phi = params  # unpack parameters 
	
    psi = chi - phi
	
    n = phase_refractive_index(params,x,y,psi)

    dndc = deriv_c(n,x,y,psi)
    dndx = deriv_x(n,x,y,psi)
    dndy = deriv_y(n,x,y,psi)

#    print('n,dndx,dndy,dndc',n,dndx,dndy,dndc)
	
    derivs = [(1./n**2)*(n*np.cos(chi)    - dndc*np.sin(chi)),      
              (1./n**2)*(n*np.sin(chi)    - dndc*np.cos(chi)),
	          (1./n**2)*(dndy*np.cos(chi) - dndx*np.sin(chi))]
    return derivs

#filename='hasel1.txt'
#target = open(filename,'w')
	
# parameters 
Xm = 2.
Y  = 2.

#dgr_phi  = 45.  # dip angle in degrees. 
                 # dip angle is magnetic inclination or angle between horizontal and magnetic field line
dgr_incl = 60.   # inclination angle in degrees 

#phi  = dgr_phi*np.pi/180. # dip angle in radians
incl = dgr_incl*np.pi/180. # inclination angle in radians

ccgs.pconstants()

Bmag = 0.5 # in Gauss
fH = 2.*np.pi*Bmag*ccgs.e/(ccgs.me*ccgs.c) 
freq = fH/Y

# initial conditions
chi0 = (np.pi/2.)-incl  # angle between horizontal axis and wave normal
x0 = 0.
y0 = 0.
# x0 = 1.e-3 # initial x position
# y0 = 1.e-3 # initial y position
	
# Do some nice plotting
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font',**{'family':'serif','serif':['Times New Roman']})
plt.rc('font',size=22.)
plt.rc('lines',linewidth=2,markeredgewidth=2.,markersize=10)
plt.rc('axes',linewidth=1.5)
plt.rcParams['xtick.major.size']=8.
plt.rcParams['xtick.minor.size']=4.
plt.rcParams['figure.subplot.bottom']=0.13
plt.rcParams['figure.subplot.left']=0.17
	
xmin = 0.   # min x boundary for plotting
xmax = 130. # max x boundary for plotting
ymin = 0.   # min y boundary for plotting 
ymax = 35.  # max y boundary for plotting

#plt.close(1)
plt.figure(1,figsize=(14,7))
plt.axis([xmin,xmax,ymin,ymax])
	
plt.text(6.,31.,r"Xm = "+str(Xm))
plt.text(6.,28.,r"Y  = "+str(Y))
plt.text(6.,25.,r"i  = "+str(dgr_incl))
plt.text(6.,22.,r"$\chi$ = " +str(90.-dgr_incl))
plt.text(6.,19.,r"f = "+str(round(freq,1))+" Hz")
plt.text(6.,16.,r"B = "+str(Bmag)+" Gauss")
		
# bundle initial conditions 
xyi0 = [x0,y0,chi0] 

# make time array for the solver 
tStop = 2.e7 # maximum time for integration
tInc  = 1.e1 # dt increment for integration
t = np.arange(0., tStop, tInc) # time array for time integration 

angle = 22.5*np.pi/180.
#angle = 11.25*np.pi/180.
for ii in range(0,5):
    phi = ii*angle # vary dip angle 
    print('Phi: ',phi)	

	# bundle parameters for ODE solver
    params = [Y, Xm, phi]
	
	# call the ODE solver
	#psoln = odeint(f, xzd0, t, args=(), full_output=True)
    psoln = odeint(f, xyi0, t, args=(params,))

	# plot results
    plt.plot(psoln[:,0]/1.e5,psoln[:,1]/1.e5,label=r"$\Phi =$"+str(np.round(180.*phi/np.pi,1)))
	
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'Haselgrove(1957)')

name='rayhasel1'

#name='rayhasel2e_phi'+str(dgr_phi)+'_refractive_index_freq'+str(round(freq,1))

plt.savefig(name+'.png')	
	
#target.close()	
plt.show()	
