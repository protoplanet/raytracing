# ------------------------------------------------------------------------------------------ #
# Description : Plotting of results from ray tracing equations by Rice (1997)  
#
# Author      : Miroslav Mocak (Slovak Organization For Space Activities)
# Date        : 31/August/2016
# Usage       : run ray2d_polar.py (this code is OPEN-SOURCE and free to be used/modified by anyone)
# References  : Yabroff (1961), Kimura (1966), Rice (1997)  
# ------------------------------------------------------------------------------------------ #

import numpy as np
import matplotlib.pyplot as plt
import ray_cmks
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import json

# ------------------------------------------------------------------------------------------------ #
#  THIS ROUTINE SETS SOME STANDDARD VALUES FOR MATLPLOTLIB TO OBBTAIN PUBLICATION-QUALITY FIGURES
# ------------------------------------------------------------------------------------------------ #
def SetMatplotlibParams():
#	plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('font',**{'family':'serif','serif':['Times New Roman']})
    plt.rc('font',size=22.)
    plt.rc('lines',linewidth=2,markeredgewidth=2.,markersize=10)
    plt.rc('axes',linewidth=1.5)
    plt.rcParams['xtick.major.size']=8.
    plt.rcParams['xtick.minor.size']=4.
    plt.rcParams['figure.subplot.bottom']=0.13
    plt.rcParams['figure.subplot.left']=0.17	

# --------------------------------- #	
# DEFINE PARAMETERS FOR FINAL PLOT	
# --------------------------------- #

def finPlot(radius,latitude,gdt,freq,dth0,dchi0,ion,ii):
    xmin = -26.e3 # min x boundary for plotting
    xmax = +26.e3 # max x boundary for plotting
    ymin = -26.e3 # min y boundary for plotting 
    ymax = +26.e3 # max y boundary for plotting

    fig=plt.figure(1,figsize=(9,8))
    plt.axis([xmin,xmax,ymin,ymax])
	
    plt.hlines(y=0.,xmin=xmin,xmax=xmax,color='k',linewidth=1)
	
	
	
    xx = radius[:]*np.cos(latitude[:])
    yy = radius[:]*np.sin(latitude[:])

    xxr = np.asarray(radius[:])*np.asarray(np.cos(45.*np.pi/180.))
    yyr = np.asarray(radius[:])*np.asarray(np.sin(45.*np.pi/180.))	
	
    ax=fig.add_subplot(1,1,1)
    circ=plt.Circle((0,0), radius=ray_cmks.Re/1.e3, color='b', fill=True)
    ax.add_patch(circ)
	
    #thtgrid = np.arange(-np.pi/2.,np.pi/2.,0.1)
    thtgrid = np.arange(0.,np.pi,0.1)
    thtgridm = np.pi/2.-thtgrid
    rmag = ray_cmks.Re*np.cos(thtgridm)*np.cos(thtgridm)/1.e3
       
    #rw = 2.*ray_cmks.Re
    xmag1 = 1.*rmag*np.cos(thtgrid)
    ymag1 = 1.*rmag*np.sin(thtgrid)

    xmag2 = 2.*rmag*np.cos(thtgrid)
    ymag2 = 2.*rmag*np.sin(thtgrid)

    xmag3 = 3.*rmag*np.cos(thtgrid)
    ymag3 = 3.*rmag*np.sin(thtgrid)

    xmag4 = 4.*rmag*np.cos(thtgrid)
    ymag4 = 4.*rmag*np.sin(thtgrid)

    plt.plot(xmag4,ymag4,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag3,ymag3,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag2,ymag2,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag1,ymag1,color='r',linewidth=0.5,linestyle='--')


    plt.plot(xmag4,-ymag4,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag3,-ymag3,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag2,-ymag2,color='r',linewidth=0.5,linestyle='--')
    plt.plot(xmag1,-ymag1,color='r',linewidth=0.5,linestyle='--')

    plt.text(-24.e3,20.e3,r"$\theta_0$ = "+str(dth0)+"$^{o}$")	
    plt.text(-24.e3,17.e3,r"$\chi_0$ = " +str(dchi0)+"$^{o}$")
#    plt.text(-24.e3,14.e3,r"f = "+str(round(freq,1))+" Hz")
	
    plt.text(-24.e3,-3.e3,r"NORTH")
    plt.text(+16.e3,+1.e3,r"SOUTH")
	
    plt.plot(np.asarray(xx)/1.e3,np.asarray(yy)/1.e3)
#    plt.plot(np.asarray(xxr)/1.e3,np.asarray(yyr)/1.e3,linestyle=':')
	
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'')

    name='ray2d_polar_theta'+str(dth0)+'_chi'+str(dchi0)+'_freq'+str(freq)+'_'+str(ion)
    dir='results/'
    plt.savefig(dir+name+'.png')

def ligvalue(ii,jj,deltat,deltaf,time,freq):
  umax = len(time)
  vmax = len(freq)
  value = 0.3e-24
  for uu in range(0,umax-1):
#    for vv in range(0,vmax-1):
#	    print(ii,jj,uu)
	    if (ii*deltat <= time[uu] <= (ii+1)*deltat) and (jj*deltaf <= freq[uu] <= (jj+1)*deltaf):
		    value = Lig1(freq[uu])
#		    print(ii,jj,ii*deltat,(ii+1)*deltat,jj*deltat,(jj+1)*deltat,time[uu],freq[uu],uu)
  return value
	
def finGdt(radius,latitude,gdtb,freqb,dth0,dchi0,ion):

    fig=plt.figure(2,figsize=(9,8))

    xmin = 0.   # min x boundary for plotting
    xmax = 0.2 # max x boundary for plotting
#    xmax = np.asarray(gdtb[0])
    ymin = 0.   # min y boundary for plotting 
    ymax = 35.  # max y boundary for plotting

    plt.axis([xmin,xmax,ymin,ymax])	
	
#    plt.semilogx(np.asarray(gdtb-gdtb[-1]),np.asarray(freqb)/1.e3)
    plt.plot(np.asarray(gdtb),np.asarray(freqb)/1.e3)
    plt.ylabel(r'frequency (10$^3$ Hz)')
    plt.xlabel(r'group delay time (s)')
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'')

    name='ray2d_polar_theta'+str(dth0)+'_chi'+str(dchi0)+'_'+str(ion)
    dir='results/'
    plt.savefig(dir+name+'.png')
	
    f = open('frequency.dat', 'w')
    json.dump(freqb, f)
    f.close()

    f = open('gdt.dat_'+str(ion[0])+str(ion[1]), 'w')
    json.dump(gdtb, f)
    f.close()

def z_func(x,y):
 return (1-(x**2+y**3))*np.exp(-(x**2+y**2)/2)
	
def finGdtSpectrogram(radius,latitude,gdtb,freqb,dth0,dchi0,ion):

    fig=plt.figure(3,figsize=(9,8))

    xmin = gdtb[-1]
    xmax = gdtb[0]
    ymin = freqb[0]
    ymax = freqb[-1]
	
    npoints = 17
    deltat = (xmax - xmin)/npoints
    deltaf = (ymax - ymin)/npoints
	
    x = np.arange(xmin, xmax, deltat)
    y = np.arange(ymin, ymax, deltaf)
    X, Y = np.meshgrid(x, y)

    print(xmin,xmax,ymin,ymax,deltat,deltaf)
	 
#    Z = z_func(X, Y) # evaluation of the function on the grid

    Z = np.empty((npoints,npoints))
    Z.fill(0.)

    for ii in range(0,npoints-1):
      for jj in range(0,npoints-1):
          Z[npoints-ii-1,jj] = ligvalue(ii,jj,deltat,deltaf,gdtb,freqb)

    print(Z)	  
    im = imshow(Z,cmap=cm.gist_stern,interpolation="none") # drawing the function
    colorbar(im) # adding the colobar on the right
    # latex fashion title
    title('$z=(1-x^2+y^3) e^{-(x^2+y^2)/2}$')
		
def Lig():

    fig=plt.figure(4,figsize=(9,8))

    # mks units 

    z0 = 377. # intrinsic impedance in Ohms
    mu0 = 8.854e-12 # permeability of free space 
    he = 5.e3 # height of the cloud above the ground (set to 5 km)
    i0 = -10.53e3 # magnitude of the downward moving current 10.53 kA
    a = 5.e3 # model parameter
    b = 1.e5 # model parameter
    dgr_kappa = 10.   # angle of observer with respect to zenith
    kappa = dgr_kappa*np.pi/180.
    R = 500.e3 # distance to observer 100 km

    fStop = 3.e4
    fInc = 1.e2
    freq = np.arange(0., fStop, fInc)
    omg = 2.*np.pi*freq

    tmp = ((omg**2)*(a-b)**2)/((omg**2+a**2)*(omg**2+b**2))
    s = (1./z0)*(((mu0*he*i0)/(2.*np.pi))**2)*((np.sin(kappa)/R)**2)*tmp

    # set parameters for plotting
#    SetMatplotlibParams()

    xmin = 0.   # min x boundary for plotting
    xmax = 30. # max x boundary for plotting
    ymin = 0.   # min y boundary for plotting 
    ymax = 2.  # max y boundary for plotting

    plt.axis([xmin,xmax,ymin,ymax])	

    plt.plot(freq/1.e3,s/1.e-24)

    plt.xlabel(r"f (10$^3$ Hz)")
    plt.ylabel(r"S (10$^{-24}$ W m$^{-2}$ Hz$^{-1}$)")
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'Bortnik (2004) and ref therein')

    name='lightning_kappa'+str(dgr_kappa)
    dir='results/'
    plt.savefig(dir+name+'.png')	

	
def Lig1(ff):

    fig=plt.figure(4,figsize=(9,8))

    # mks units 

    z0 = 377. # intrinsic impedance in Ohms
    mu0 = 8.854e-12 # permeability of free space 
    he = 5.e3 # height of the cloud above the ground (set to 5 km)
    i0 = -10.53e3 # magnitude of the downward moving current 10.53 kA
    a = 5.e3 # model parameter
    b = 1.e5 # model parameter
    dgr_kappa = 10.   # angle of observer with respect to zenith
    kappa = dgr_kappa*np.pi/180.
    R = 500.e3 # distance to observer 100 km

    omgff = 2.*np.pi*ff
    tmpff = ((omgff**2)*(a-b)**2)/((omgff**2+a**2)*(omgff**2+b**2))
    sff = (1./z0)*(((mu0*he*i0)/(2.*np.pi))**2)*((np.sin(kappa)/R)**2)*tmpff
    
    return sff   

def plotne():
    ray_cmks.pconstants()	
	
    rrmin = ray_cmks.Re
    rrmax = 3.*ray_cmks.Re
    npoints = 100
    rrdelta = (rrmax-rrmin)/npoints
    rr = np.arange(rrmin, rrmax, rrdelta)
    rrE = rr/ray_cmks.Re
    
    ne = 1.e6*(1.8e5*np.exp(-4.183119*(rrE-1.0471)))	
	
    plt.plot(rr/1.e3,ne)
	
    plt.xlabel(r"r (km)")
    plt.ylabel(r"ne (10$^{-24}$ W m$^{-2}$ Hz$^{-1}$)")
    plt.legend(loc=1,prop={'size':14})
    plt.title(r'Yabroff')
	