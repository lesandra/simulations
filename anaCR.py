#Analyse ZHaireS + ComputeVoltage outputs 
# Produced by Nicolas (see email on 24/10/2017)

import sys
import numpy as np
import pylab as pl

pi = 3.1415927
nantmin = 6
ndays = 1  # duration of observation (days)

# Flux
J1 = 1.06e29  # /s/m2/sr/eV
gamma1 = 3.26 
J2 = 4.5e17
gamma2 = 2.65
eth=[pow(10,float(x)/10) for x in range(160,200)]
f1 = [J1*pow(x,-gamma1)*pow(x,3)/1e24 for x in eth]  # below ankle
f2 = [J2*pow(x,-gamma2)*pow(x,3)/1e24 for x in eth]  # above ankle

# Load results
#datadir = '/home/martineau/GRAND/GRAND/data/zhaires/CRs/'
resfile = 'CR_detection_count.txt'
res = np.loadtxt(resfile)
Elog = res[:,0]
eu = np.unique(Elog)
E= pow(10,Elog)
theta = res[:,1]
th = np.unique(theta)
sn = np.cos(theta*pi/180)*0.5*0.5  # Try area 
phi = res[:,2]
xc = res[:,3]
yc = res[:,4]
zc = res[:,5]
newc = res[:,6]
nnsc = res[:,7]
ntc = res[:,8]
newa = res[:,9]
nnsa = res[:,10]
nta = res[:,11]
  
## Plots
# Core positions
#pl.figure(0)
#pl.plot(xc,yc,'o')
#pl.show()
# Flux
#pl.figure(1)
#pl.loglog(eth,f1,label='flux below ankle')
#pl.loglog(eth,f2,label='flux above ankle')
#pl.xlabel('Energy (eV)')
#pl.ylabel('Flux*E$^3$ (eV$^2$/s/sr/m$^2$)')
#pl.legend(loc='best')
#pl.title("Fly's Eye flux - astro-ph 1511.07710")
#pl.grid(True)
#pl.show()

# Compute apperture
sefft = np.zeros(np.size(eu))
sefftopt = np.zeros(np.size(eu))
for i in range(np.size(eu)):
  e = eu[i]
  seffp = np.zeros(np.size(th))
  seffpopt = np.zeros(np.size(th))
  for t in range(np.size(th)):
    sel = np.where( (Elog == e) & (theta == th[t]) )
    sphi = sn[sel[0][0]]*200000*4  # GRAND physical area projected along traj  [km2]. Factor 4 comes from try area = 0.25km^2.
    ok = np.logical_or(newa[sel]>=nantmin,nnsa[sel]>=nantmin) 
    r = float(np.sum(ok))/np.size(ok)
    seffp[t] = 2*pi*sphi*r   # Effective area integrated over phi
    seffpopt[t] = 2*pi*sphi
    #print '\n##E = 10^',e,',eV, theta = ',th[t],'deg:'
    #print 'Sphi =',sphi,'km^2'
    #print 'phi=',phi[sel]
    #print 'Nants EW=',newa[sel]
    #print 'Nants NS=',nnsa[sel]
    #print ok
    #print 'Effective area',seffp[t],'km^2'
  print '\n##E = 10^',e,'eV'
  print seffp  
  sefft[i] = np.trapz(seffp*np.sin(th*pi/180),th*pi/180)
  sefftopt[i] = np.trapz(seffpopt*np.sin(th*pi/180),th*pi/180)
  
## Compute event rate
print ''
ee = pow(10,eu)
print 'Energy:',ee,'eV'
print 'Apperture:',sefft,'km2.sr'
print 'Optimal apperture:',sefftopt[0],'km2.sr'
sefft = sefft*1e6  #m2
sefftopt = sefftopt*1e6  #m2
dt = 3600*24*ndays  #1 day
evt1d = np.trapz(sefft*J1*pow(ee,-gamma1)*dt,ee)
evt1dopt = np.trapz(sefftopt*J1*pow(ee,-gamma1)*dt,ee)
print 'Expected daily event rate in 10^17-10^18eV:',evt1d
print 'Optimal daily event rate in 10^17-10^18eV:',evt1dopt
#pl.figure(2)
#pl.plot(eu,sefft)
#pl.plot(eu,sefftopt)
#pl.grid(True)
#pl.show() 

# HE events
e2 = [pow(10,float(x)/10) for x in range(190,200)]
f2 =[J2*pow(x,-gamma2) for x in e2]
evt1doptHE = sefftopt[0]*dt*np.trapz(f2,e2)
print 'Optimal daily event rate in 10^19-10^20eV:',evt1doptHE


# Apperture vs zenith
sefft = np.zeros(np.size(th))
for t in range(np.size(th)):
  sel = np.where( (theta == th[t]) )
  sphi = sn[sel[0][0]]*200000*4  # GRAND physical area projected along traj  [km2]. Factor 4 comes from try area = 0.25km^2.
  ok = np.logical_or(newa[sel]>=nantmin,nnsa[sel]>=nantmin) 
  # Warning: here all E count the same <=> assuming flat spectrum
  r = float(np.sum(ok))/np.size(ok)
  sefft[t] = 2*pi*sphi*r   # Effective area integrated over phi
  #print '\n##E = 10^',e,',eV, theta = ',th[t],'deg:'
  #print 'Sphi =',sphi,'km^2'
  #print 'phi=',phi[sel]
  #print 'Nants EW=',newa[sel]
  #print 'Nants NS=',nnsa[sel]
  #print ok
  #print 'Effective area',seffp[t],'km^2'

print '\n##theta =', th,'deg'
print 'Apperture (km^2.sr):',sefft
pl.figure(3)
pl.plot(th,sefft)
pl.xlabel('Zenith angle (deg)')
pl.ylabel('Apperture (km$^2$.sr)')
pl.title("GRAND apperture (for a flat spectrum)")
pl.grid(True)
pl.show()
