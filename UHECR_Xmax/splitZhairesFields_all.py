
# this python script should mimic the functionaliyt of processZhairesShowers.m
# this means it reads in the original ZHAires_output timefresnel_root.dat and splits it into single antenna files a_#.dat. It creates also antpos.dat with all antenna positions
# for filtering and hilbert envelope, please use the next script


# As input one needs the path to the folder containing the simulations like AzimuthScan4/AzimuthScan4_8 and the simulation name again like AzimuthScan4_8


import sys
from sys import argv

import numpy as np
import pylab as pl

import os
import re

pl.ion()
DISPLAY = 0


# hand over the path of the simulations
wkdir = sys.argv[1] #'/media/sf_work/Paris/scripts_GRAND/' # path where the simulation file is
fname = wkdir+'timefresnel-root.dat' #ZHAires_output

print 'file ' +fname+ ' gets split up into single antenna file'

# First load file
a = np.loadtxt(fname, dtype='float', comments='#')
shId = a[:,0]
antId =  a[:,1]
x =   a[:,2]  #X = S->N
y =   a[:,3]  #Y = E->W
z =   a[:,4]  # Z = Up
t =   a[:,5]  #ns 
Ex =  a[:,11]  # V/m 
Ey =  a[:,12]
#Ez =  a[:,13]



### needed for CR showers as well
print 'Ez and zpos flipped to correct for coordinate system used for simulations'
Ez = -1.* a[:,13]
for i in range(0, len(z)):
    z[i]= 100000 -2650 -z[i] # upper bound Aires atmosphere- ground altitude (Ulastai) - actual position in time frensel = input antenna position
    

    
#print 'Flipping z axis...'

#print len(z)

# file for antenna positions
file_antpos= wkdir+"antpos.dat"
FILE2 = open(file_antpos, "w" )

print "time in ns, Efield in muV/m" 
#print len(antId)/600
# Now split antenna data  
for i in range(0, len(antId)/120):
  sel = np.where(antId == i+1)[0]
  #print 'Antenna',int(antId[sel[0]]), ':',x[sel[0]],y[sel[0]],z[sel[0]],'m'
  try: 
   print >>FILE2,"%.2f	%.2f	%.2f" % (x[sel[0]], y[sel[0]],z[sel[0]] )

  except IndexError: # catch the error
    continue 

  ti = t[sel] #*1e-3  # [ns] 
  Exi = Ex[sel]*1e6   # [muV/m]
  Eyi = Ey[sel]*1e6   # [muV/m]
  Ezi = Ez[sel]*1e6   # [muV/m]
  #print antId[i]
  
  
  #filename = "a"+str(int(antId[i]))+".trace"
  filename = wkdir+"a"+str(i)+".trace"
  alld = np.transpose([ti,Exi,Eyi,Ezi])
  np.savetxt(filename,alld,fmt='%.6e')  
  ti0 = ti-ti[0]
  

  
  if DISPLAY:
    pl.figure(1)
    pl.plot(ti0,Exi,label='Ex (North-South)')
    pl.plot(ti0,Eyi,label='Ey (East-West)')
    pl.plot(ti0,Ezi,label='Ez (Up-Down)')
    pl.xlabel('Time [mus]')
    pl.ylabel('Amplitude [muV/m]') 
    pl.grid(True)
    pl.legend(loc='lower right');
    pl.show()
    raw_input()


FILE2.close()
print 'The end'
