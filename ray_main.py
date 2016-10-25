# ------------------------------------------------------------------------------------------ #
# Description : Python implementation of ray tracing equations stated in PhD thesis of Rice (1997)  
#               Electron/proton stratification according to Yabroff (1961) + IRI model available
#               Geometry is 2D polar   
#
# Author      : Miroslav Mocak (Slovak Organization For Space Activities)
# Date        : 14/October/2016
# Usage       : run ray_main.py (this code is OPEN-SOURCE and free to be used/modified by anyone)
# References  : Yabroff (1961), Kimura (1966), Rice (1997)  
# ------------------------------------------------------------------------------------------ #

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import warnings
import re  # python regular expressions
import ray_cmks
import ray_fncts 
import ray_plot
import sys

warnings.filterwarnings('ignore')

# ---------------------------------------------- #		
# READ INPUT PARAMETERS AND INITIAL CONDITIONS	
# ---------------------------------------------- #	

file=open('ray_param.dat','r')
next(file) # skip header line
next(file) # skip header line

input=[]
for line in file:
	prsvalue = re.search(r'\[(.*)\]', line).group(1) # parse out values from square brackets
	input.append(prsvalue)
file.close()

freq_in   = float(input[0])
freq_en   = float(input[1])
freq_de   = float(input[2]) 
orbit     = float(input[3]) 
frr0      = float(input[4])
dth0      = float(input[5])
dchi0     = float(input[6])
tt0       = float(input[7])
tstop     = float(input[8])
nsteps    = float(input[9])
pvode     = float(input[10])
pzvode    = float(input[11])
plsoda    = float(input[12])
pdopri5   = float(input[13])
pdop853   = float(input[14])
iono_el   = float(input[15])
iono_np	  = float(input[16])
iono_ir   = float(input[17])
iri_fna   = input[18]
pwhist    = float(input[19])
	
if (iono_el == 1. and iono_np == 1.) or (iono_el == 1. and iono_ir == 1.) or (iono_np == 1. and iono_ir == 1.):
    print('STOP (in ray_main.py): choose only one ionospheric model')
    sys.exit()
	
# load constants
ray_cmks.pconstants()

rr0   = frr0*ray_cmks.Re  # initial altitude :: approx 300 km = 1.0471*Re
th0   = dth0*np.pi/180.
chi0  = dchi0*np.pi/180.
G0 = tt0
t0 = 0.
dt = tstop/nsteps # calculate integration step
	
	
	
# bundle initial conditions 
rtcG0 = [rr0,th0,chi0,G0] 
	
# introduce some error handling !!
	
# select chosen ionospheric model
ion = ["0","0"]

# initialize ionosphere
height = []
ne = []
H_ions_per = []  
O_ions_per = [] 
He_ions_per = []
O2_ions_per = []
NO_ions_per = []
N_ions_per = []
ionosphere = [height,ne,O_ions_per,H_ions_per,He_ions_per,O2_ions_per,NO_ions_per,N_ions_per] 
 
if iono_el == 1:
    fname = ""
    ion = ["iono_el",fname]
	
if iono_np == 1:
    fname = ""
    ion = ["iono_np",fname]
	
if iono_ir == 1:
#    fname = "iri_2012_24537_night.lst"
#    fname = "iri_2012_25962_day.lst"
    fname = iri_fna
#    print(fname)
    ion = ["iono_ir",fname]
    height = []
    ne = []
    O_ions_per = []
    H_ions_per = []
    He_ions_per = []
    O2_ions_per = []
    NO_ions_per = []
    N_ions_per = []

    # Open file
    # fname = 'iri_2012_22651_night.lst'
    fname = ion[1]

    f = open('IRI/'+fname, 'r')
	
    # Loop over lines and extract variables of interest
    for line in f:
        line = line.strip()
        columns = line.split()
    #	if float(columns[5]) = -1.:
    #	 columns[5] = 0.
        if float(columns[8]) < 0.:
          columns[8] = 0.
        if float(columns[9]) < 0.:
         columns[9] = 0.
        if float(columns[10]) < 0.:
          columns[10] = 0.
        if float(columns[11]) < 0.:
          columns[11] = 0.
        if float(columns[12]) < 0.:
          columns[12] = 0. 
        if float(columns[13]) < 0.:
          columns[13] = 0. 
        if float(columns[14]) < 0.:
          columns[14] = 0.	 
        height.append(float(columns[5])) # height in km
        ne.append(float(columns[8])) # electron density in m-3
        O_ions_per.append(float(columns[9])) # atomic oxygen O+ ions percentage
        H_ions_per.append(float(columns[10])) # atomic hydrogen H+ ions percentage
        He_ions_per.append(float(columns[11])) # atomic helium He+ ions percentage
        O2_ions_per.append(float(columns[12])) # molecular oxygen O2+ ions percentage
        NO_ions_per.append(float(columns[13])) # nitric oxide ions NO+ percentage
        N_ions_per.append(float(columns[14])) # atomic nitrogen N+ ions percentage
    f.close()
	
    if np.asarray(height[-1]) < orbit:
         print('STOP (in ray_main.py): limiting orbit exceeds max altitude of IRI model')
         sys.exit()
	
    ionosphere = [height,ne,O_ions_per,H_ions_per,He_ions_per,O2_ions_per,NO_ions_per,N_ions_per]

#print(height[0],rr0-6371200.0,height[-1])	
	
if ion == ["0","0"]:
    print("Error in ionospheric model")
	
# ---------------------------------------------- #	
# CALCULATE RHS OF RAY TRACING EQUATIONS	
# ---------------------------------------------- #
	
def f(t, rtcG, freq):
    rr, th, chi, G = rtcG     # unpack current values 
 	
    c = 2.99792458e8
    n = ray_fncts.phase_refractive_index(rr,th,chi,freq,ion,ionosphere)
	
    dndrr = ray_fncts.deriv_rr(rr,th,chi,freq,ion,ionosphere)
    dndth = ray_fncts.deriv_th(rr,th,chi,freq,ion,ionosphere)
    dndfr = ray_fncts.deriv_fr(rr,th,chi,freq,ion,ionosphere)
    dndch = -n[1]

#   ngroup = n[0]+freq*dndfr
	
    derivs = [(1./n[0]**2)*(n[0]*np.cos(chi) + dndch*np.sin(chi)),	
              (1./(rr*n[0]**2))*(n[0]*np.sin(chi) - dndch*np.cos(chi)),
			  (1./(rr*n[0]**2))*(dndth*np.cos(chi) - (rr*dndrr+n[0])*np.sin(chi)),
			  (1./c)*(1.+(freq/n[0])*dndfr)]
    return derivs
	
# ---------------------------------------------- #
# MAIN CALLS ODE SOLVER AND STORES RESULTS
# ---------------------------------------------- #

if pvode == 1:
	intype="vode"
	
if pzvode == 1:
	intype="zvode"
	
if plsoda == 1:
	intype="lsoda"
	
if pdopri5 == 1:
	intype="dopri5"
	
if pdop853 == 1:
	intype="dop853"
	
print('Using ODE integrator: '+str(intype))
print('Limiting height: '+str(orbit)+str(' km'))

# set parameters for plotting	
ray_plot.SetMatplotlibParams()
	
fd = freq_de
fend = int((freq_en-freq_in)/freq_de)


# initial array for frequency and group delay time at chosen orbit
freqb = []	
gdtb = []

for ii in range(1,fend+1):
    freq = ii*fd # vary frequency 
    print('Calculating ray path: '+str("%.2g" % freq)+' Hz')
	
    psoln = ode(f).set_integrator(intype,method='bdf')
    psoln.set_initial_value(rtcG0,t0).set_f_params(freq)

    radius = []
    latitude = []
    gdt = []

    while psoln.successful() and psoln.t < tstop and psoln.y[0] > ray_cmks.Re and psoln.y[0] < (ray_cmks.Re+orbit*1.e3):
        psoln.integrate(psoln.t+dt)
        radius.append(psoln.y[0])
        latitude.append(psoln.y[1])
        gdt.append(psoln.y[3])
#        print(psoln.y[2],(180./np.pi)*psoln.y[2])
		
    xx = radius[:]*np.cos(latitude[:])
    yy = radius[:]*np.sin(latitude[:])
    freqb.append(freq)
    gdtb.append(gdt[-1])
	
    ray_plot.finPlot(radius,latitude,gdt,freq,dth0,dchi0,ion,ii)		
	
# ---------------------------------------------- #
# ray_plot RESULTS
# ---------------------------------------------- #

if pwhist == 1:
    ray_plot.finGdt(radius,latitude,gdtb,freqb,dth0,dchi0,ion)
	
plt.show()
plt.clf()	

# ---------------------------------------------- #
# END
# ---------------------------------------------- #
