#Analyse ZHaireS + ComputeVoltage outputs 
# Produced by Nicolas (see email on 24/10/2017)

import sys
import numpy as np
import pylab as pl

pi = 3.1415927
nantmin = 6
ndays = 1  # duration of observation (days)

# Flux (TA - astro-ph 1511.07510)
J1 = 1.06e29  # /s/m2/sr/eV
gamma1 = 3.26 
J2 = 4.5e17
gamma2 = 2.65
eth=[pow(10,float(x)/10) for x in range(160,200)]
f1 = [J1*pow(x,-gamma1)*pow(x,3)/1e24 for x in eth]  # below ankle
f2 = [J2*pow(x,-gamma2)*pow(x,3)/1e24 for x in eth]  # above ankle

# Load results
#datadir = '/home/martineau/GRAND/GRAND/data/zhaires/CRs/'
resfile = 'detection_count.txt'
res = np.loadtxt(resfile)
Elog = res[:,0]
eu = np.unique(Elog)
sefftc = np.zeros([1,1])
eu[0] = 19
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
pl.figure(0)
pl.plot(xc,yc,'o')
#pl.show()
# Flux
#pl.figure(1)
#pl.loglog(eth,f1,label='flux below ankle')
#pl.loglog(eth,f2,label='flux above ankle')
#pl.xlabel('Energy (eV)')
#pl.ylabel('Flux*E$^3$ (eV$^2$/s/sr/m$^2$)')
#pl.legend(loc='best')
#pl.title("TA flux - astro-ph 1511.07510")
#pl.grid(True)
#pl.show()

# Compute apperture
sefftc = np.zeros(np.size(eu))
seffta = np.zeros(np.size(eu))
sefftopt = np.zeros(np.size(eu))
for i in range(np.size(eu)):
  e = eu[i]
  seffpa = np.zeros(np.size(th))
  seffpc = np.zeros(np.size(th))
  seffpopt = np.zeros(np.size(th))
  for t in range(np.size(th)):
    sel = np.where( (Elog == e) & (theta == th[t]) )
    sphi = sn[sel[0][0]]*200000*4  # GRAND physical area projected along traj  [km2]. Factor 4 comes from try area = 0.25km^2.
    oka = np.logical_or(newa[sel]>=nantmin,nnsa[sel]>=nantmin) 
    okc = np.logical_or(newc[sel]>=nantmin,nnsc[sel]>=nantmin) 
    ra = float(np.sum(oka))/np.size(oka)
    rc = float(np.sum(okc))/np.size(okc)
    seffpa[t] = 2*pi*sphi*ra   # Effective area integrated over phi
    seffpc[t] = 2*pi*sphi*rc   # Effective area integrated over phi
    seffpopt[t] = 2*pi*sphi
    #print '\n##E = 10^',e,',eV, theta = ',th[t],'deg:'
    #print 'Sphi =',sphi,'km^2'
    #print 'phi=',phi[sel]
    #print 'Nants EW=',newa[sel]
    #print 'Nants NS=',nnsa[sel]
    #print ok
    #print 'Effective area',seffp[t],'km^2'
  print '\n##E = 10^',e,'eV'
  print seffpa
  print seffpc
    
  seffta[i] = np.trapz(seffpa*np.sin(th*pi/180),th*pi/180)
  sefftc[i] = np.trapz(seffpc*np.sin(th*pi/180),th*pi/180)
  sefftopt[i] = np.trapz(seffpopt*np.sin(th*pi/180),th*pi/180)
  
## Compute event rate
print ''
ee = pow(10,eu)
print 'Energy:',ee,'eV'
print 'Apperture (agressive):',seffta,'km2.sr'
print 'Apperture (conservative):',sefftc,'km2.sr'
print 'Optimal apperture:',sefftopt[0],'km2.sr'
pl.figure(2)
pl.plot(eu,seffta,label='Agressive')
pl.plot(eu,sefftc,label='Conservative')
pl.plot(eu,sefftopt,label='Optimal')
pl.xlabel('Energy (eV)')
pl.ylabel('Apperture (km$^2$.sr)')
pl.grid(True)
pl.legend(loc='best')
pl.title("GRAND CR apperture")
#pl.show() 

seffta = seffta*1e6  #m2
sefftc = sefftc*1e6  #m2
sefftopt = sefftopt*1e6  #m2
dt = 3600*24*ndays  #1 day
evt1da = np.trapz(seffta*J1*pow(ee,-gamma1)*dt,ee)
evt1dc = np.trapz(sefftc*J1*pow(ee,-gamma1)*dt,ee)
evt1dopt = np.trapz(sefftopt*J1*pow(ee,-gamma1)*dt,ee)
print 'Expected daily event rate in 10^17-10^18eV (agressive):',evt1da
print 'Expected daily event rate in 10^17-10^18eV (conservative):',evt1dc
print 'Optimal daily event rate in 10^17-10^18eV:',evt1dopt

# HE events
e2 = [pow(10,float(x)/10) for x in range(190,200)]
f2 =[J2*pow(x,-gamma2) for x in e2]
evt1doptHE = sefftopt[0]*dt*np.trapz(f2,e2)
print 'Optimal daily event rate in 10^19-10^20eV:',evt1doptHE


# Apperture vs zenith
seffta = np.zeros(np.size(th))
sefftc = np.zeros(np.size(th))
for t in range(np.size(th)):
  sel = np.where( (theta == th[t]) )
  sphi = sn[sel[0][0]]*200000*4  # GRAND physical area projected along traj  [km2]. Factor 4 comes from try area = 0.25km^2.
  oka = np.logical_or(newa[sel]>=nantmin,nnsa[sel]>=nantmin) 
  okc = np.logical_or(newc[sel]>=nantmin,nnsc[sel]>=nantmin) 
  # Warning: here all E count the same <=> assuming flat spectrum
  ra = float(np.sum(oka))/np.size(oka)
  rc = float(np.sum(okc))/np.size(okc)
  seffta[t] = 2*pi*sphi*ra/(2*pi)   # Effective area integrated over phi and then averaged
  sefftc[t] = 2*pi*sphi*rc/(2*pi) 
  #print '\n##E = 10^',e,',eV, theta = ',th[t],'deg:'
  #print 'Sphi =',sphi,'km^2'
  #print 'phi=',phi[sel]
  #print 'Nants EW=',newa[sel]
  #print 'Nants NS=',nnsa[sel]
  #print ok
  #print 'Effective area',seffp[t],'km^2'

print '\n##theta =', th,'deg'
print 'Apperture (agressive):',seffta,'km^2.sr'
print 'Apperture (conservative):',sefftc,'km^2.sr'

pl.figure(3)

pl.plot(th,seffta,label='Agressive')
pl.plot(th,sefftc,label='Conservative')
pl.xlabel('Zenith angle (deg)')
pl.ylabel('Eff. area$_{<\phi>}$ (km$^2$)')
pl.title("GRAND CR apperture (for E=10$^{19}$eV)")
pl.grid(True)
pl.legend(loc='best')
pl.show()
