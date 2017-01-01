import numpy as np
import matplotlib.pyplot as plt
import re
import cmks
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

height = []
ne = []
O_ions_per = []
H_ions_per = []
He_ions_per = []
O2_ions_per = []
NO_ions_per = []
N_ions_per = []

# Open file
#fname = 'iri_2012_25962_day.lst'
fname = 'iri_2012_24537_night.lst'

f = open('../IRI/'+fname, 'r')
#f = open(fname, 'r')

# Loop over lines and extract variables of interest
for line in f:
    line = line.strip()
    columns = line.split()
#	if float(columns[5]) < 0.:
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

#print(height)

# convert percentage to density 
O_ions  = np.asarray(O_ions_per)*np.asarray(ne)/100.
H_ions  = np.asarray(H_ions_per)*np.asarray(ne)/100.
He_ions = np.asarray(He_ions_per)*np.asarray(ne)/100.
O2_ions = np.asarray(O2_ions_per)*np.asarray(ne)/100.
NO_ions = np.asarray(NO_ions_per)*np.asarray(ne)/100.
N_ions  = np.asarray(N_ions_per)*np.asarray(ne)/100.

# Yabroff
cmks.pconstants()
rrE = (cmks.Re+np.asarray(height)*1000.)/cmks.Re
neY = 1.e6*(1.8e5*np.exp(-4.183119*(rrE-1.0471)))	
#print(neY)

SetMatplotlibParams()

xmin = 60. # min x boundary for plotting
xmax = 4000. # max x boundary for plotting
ymin = 1.e7 # min y boundary for plotting 
ymax = 1.e13 # max y boundary for plotting

fig=plt.figure(1,figsize=(9,8))
plt.axis([xmin,xmax,ymin,ymax])

plt.loglog(np.asarray(height),np.asarray(ne),label=r'n$_e$')
plt.loglog(np.asarray(height),H_ions,label='H_ions')
plt.loglog(np.asarray(height),He_ions,label='He_ions')
plt.loglog(np.asarray(height),O2_ions,label='O2_ions')
plt.loglog(np.asarray(height),O_ions,label='O_ions')
plt.loglog(np.asarray(height),NO_ions,label='NO_ions')
plt.loglog(np.asarray(height),N_ions,label='N_ions')
plt.plot(np.asarray(height),neY,label=r'n$_e$ (Yabroff, 1961)',linestyle='--')
#print(height,ne)

perigeum = 450. 
apogeum = 720.

plt.axvspan(perigeum, apogeum, color='y', alpha=0.5, lw=0)

plt.xlabel(r"r (km)")
plt.ylabel(r"plasma density (m$^{-3}$)")
plt.legend(loc=1,prop={'size':16})
#plt.title(fname)

name='ionosphere_'+fname
dir='../IRI/'
plt.savefig(dir+name+'.png')

plt.show()

#f = open('height.dat', 'w')
#json.dump(height, f)
#f.close()

#f = open('ne.dat', 'w')
#json.dump(ne, f)
#f.close()

#f = open('height.dat', 'r')
#heightread = json.load(f)
#f.close()

#f = open('ne.dat', 'r')
#neread = json.load(f)
#f.close()

#SetMatplotlibParams()

#xmin = 60. # min x boundary for plotting
#xmax = 2000. # max x boundary for plotting
#ymin = 1.e7 # min y boundary for plotting 
#ymax = 1.e13 # max y boundary for plotting

#fig=plt.figure(1,figsize=(9,8))
#plt.axis([xmin,xmax,ymin,ymax])

#plt.loglog(np.asarray(heightread),np.asarray(neread),label=r'n$_e$')

#plt.xlabel(r"r (km)")
#plt.ylabel(r"plasma density (m$^{-3}$)")
#plt.legend(loc=1,prop={'size':12})

#plt.show()

