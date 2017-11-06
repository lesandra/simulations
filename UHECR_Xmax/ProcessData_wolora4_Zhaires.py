import sys
from sys import argv

import numpy as np
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
import ownprocess_wolora4_Zhaires as ow
import process_func as prf
import ska_stijn_wolora4_Zhaires as ska
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


# Parabola model(x, a=1.0, b=700., c=1.65e9), b=sx, c=sy
def model(x, a, b, c):
    return a*(x-b)**2 +c

# Needed long-file

script = sys.argv[0]  ### name of script
datadir = sys.argv[1]  ### path to data read in
fileno= sys.argv[2] ### run nr
outputfile= sys.argv[3] ### Outputfile, Pickle

# Frequencies
lowco= int(sys.argv[5])#50
hico=int(sys.argv[6])#350 #80

# Output file, pickle
outputf=sys.argv[7]

# Antenna file GRAND
Map = sys.argv[8] #   skafile=open("/gluster/data/zilles/SKA2/scripts//SKA_FlowerAntenna.list")

print "\nOutput file: " + outputfile

# Handing over the simulation which should be used
filename = sys.argv[9]

l_sim=open(filename).read().splitlines()
txt = np.array(l_sim)
files = len(txt)

r=np.zeros([files,7])
r2=np.zeros([files,7])
r3=np.zeros([files,7])

## Filtering 50-200 MHz etc and writing pickle file;
## BUT no antenna model (used to project x,y,z component to reality, needed e.g. for triggering, polarisation study)
## TimeBin 1.25ns

# Script for filtering etc.
#ow.Process_Data(datadir, fileno, outputfile, files, lowco, hico, filename)


for i in np.arange(int(files)-1):
    print "\nSimulated event: %i" %i
    results = ska.reverseAnalysis(outputfile, eventno=fileno, eventno2=str(txt[i]), simevent=i, outputfolder=outputf, SKAmap=Map)
    print "Final results: reconstr new, reconstr old, real Xmax "
    print results [0:2]

    r2[i]=np.array([results[0][0], results[1], results[2], np.abs(results[0][0]- results[1]), np.abs(results[0][0]- results[2]), results[3], results[7]])
    r[i]=np.array([results[0][1], results[1], results[2], np.abs(results[0][1]- results[1]), np.abs(results[0][1]- results[2]), results[3], results[7]])
    r3[i]=np.array([results[0][2], results[1], results[2], np.abs(results[0][2]- results[1]), np.abs(results[0][2]- results[2]), results[3], results[7]])

#print results[0][0]
#print results[0][1]

#print r[:,3] # = xmax uncertainty for single data-simulation comparison

a = r[:,3]
a.sort() 

a2 = r2[:,3]
a2.sort()

a3 = r3[:,3]
a3.sort()

#print a
#print a2

total=files

proz=0
position=0
for m in np.arange(int(files)):  # sort them
  if(proz < 0.68*total): # get the 68% reconstruction uncertainty
    #proz+=a[m]
    proz+=1
  if(proz >= 0.68*total):
      #print a[m], m
    position = a[m]
    break

proz2=0
position2=0
for m in np.arange(int(files)):  # sort them
    if(proz2 < 0.68*total): # get the 68% reconstruction uncertainty
      proz2+=1
    if(proz2 >= 0.68*total):
      position2 = a2[m]
      break

proz3=0
position3=0
for m in np.arange(int(files)):  # sort them
    if(proz3 < 0.68*total): # get the 68% reconstruction uncertainty
        proz3+=1
    if(proz3 >= 0.68*total):
        position3 = a3[m]
        break


#print proz, 0.68*total
 
 
#### Plot the reconstruction uncertainty as a histogram: new fit
fig=plt.figure(1,figsize=(8,6)) 
 
plt.hist(r2.T[3], bins=20, histtype='step')
plt.axvline(position2,color='k', linestyle='--')

plt.ylabel("Nr. of Sim.")
plt.xlabel("New Abs(Xreco-Xreal) (g/cm$^2$)")

plt.show()
name = outputf+'Histtest2_{0}_unc{1}.pdf'.format(sys.argv[2], position2)
plt.savefig(name)

#### Save all of it as a text file for later analysis
#namefile = outputf+'Histtest_{0}_unc{1}.dat'.format(sys.argv[2], position)
#file= open(namefile, 'w')  
#file.write('reconstr, real Xmax, best sim xmax, rec. -real, rec.-best,  xmaxreco_stijn, primary\n')
#for s in range(0,int(files)): 
#	  file.write('{0}	{1}	{2}	{3}	{4}	{5}	{6}\n'.format(r[s,0],r[s,1],r[s,2],r[s,3],r[s,4],r[s,5],r[s,6]))
#file.close()


#### Plot the reconstruction uncertainty as a histogram: old fit
fig=plt.figure(2,figsize=(8,6))

plt.hist(r.T[3], bins=20, histtype='step')
plt.axvline(position,color='k', linestyle='--')

plt.ylabel("Nr. of Sim.")
plt.xlabel("Abs(Xreco-Xreal) (g/cm$^2$)")

plt.show()
name = outputf+'Histtest_{0}_unc{1}.pdf'.format(sys.argv[2], position)
plt.savefig(name)


#### Plot the reconstruction uncertainty as a histogram: new fit
#fig=plt.figure(33,figsize=(8,6))
#
#plt.hist(r3.T[3], bins=20, histtype='step')
#plt.axvline(position3,color='k', linestyle='--')
#
#plt.ylabel("Nr. of Sim.")
#plt.xlabel("New Abs(Xreco-Xreal) (g/cm$^2$), other method")
#
#plt.show()
#name = outputf+'Histtest3_{0}_unc{1}.pdf'.format(sys.argv[2], position3)
#plt.savefig(name)
#
