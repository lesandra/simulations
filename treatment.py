import numpy as np
from scipy.fftpack import rfft, irfft, rfftfreq
from random import *
import matplotlib.pyplot as plt
import sys


def Filtering(v,tstep,FREQMIN,FREQMAX):
    F=rfftfreq(len(v))/tstep #len(v) points between 0 and 1/2tstep
    V=rfft(v)
    V[F<FREQMIN]=0
    V[F>FREQMAX]=0
    return irfft(V)

def Digitization(v,t,tstep,TSAMPLING,SAMPLESIZE):
    vf=np.zeros(SAMPLESIZE)
    tf=np.zeros(SAMPLESIZE)  
    ratio=int(round(TSAMPLING/tstep))
    ind=np.arange(0,int(np.floor(len(v)/ratio)))*ratio
    if len(ind)>SAMPLESIZE:
        ind=ind[0:SAMPLING]
    vf[0:len(ind)]=v[ind]
    tf[0:len(ind)]=t[ind]
    for k in range(len(ind),SAMPLESIZE):
        tf[k]=tf[k-1]+TSAMPLING
    return vf,tf

def Addnoise(vrms,v):
    for i in range(len(v)):
        v[i]=v[i]+gauss(0,vrms)
    return v


#units:
TSAMPLING=2e-9 #sec, eq to 500 MHz, sampling of system 
FREQMIN=50e6 #Hz  # frequencies which will be used in the later analysis: 50-200MHz
FREQMAX=200e6 #Hz, 250MHz
tstep=1e-9 #=t[1]-t[0], sec, time bins in simulations
SAMPLESIZE=int(3e-6/TSAMPLING) #=1500, 3e-6sec length # traces lenth of system
vrms=15 #uvolts

#NOTE: hand over path to file and antenna idea
voltage_trace=str(sys.argv[1])+"/out_"+str(sys.argv[2])+".txt"
text=np.loadtxt(voltage_trace)#'out_128.txt')
t=text[:,0] # in s
vx=text[:,1] #EW axis antenna
vy=text[:,2] #NS axis antenna

DISPLAY=1
if DISPLAY==1:
    plt.plot(t*1e9,vx) # s*1e9 = ns
    plt.plot(t*1e9,vy)
    plt.xlabel('Time [ns]')
    plt.ylabel('Voltage [uV]')
    plt.show()

vx=Filtering(vx,tstep,FREQMIN,FREQMAX)
vy=Filtering(vy,tstep,FREQMIN,FREQMAX)
if DISPLAY==1:
    plt.plot(t*1e9,vx)
    plt.plot(t*1e9,vy)
    plt.xlabel('Time [ns]')
    plt.ylabel('Voltage [uV]')
    plt.show()

vx,tx=Digitization(vx,t,tstep,TSAMPLING,SAMPLESIZE)
vy,ty=Digitization(vy,t,tstep,TSAMPLING,SAMPLESIZE)
if DISPLAY==1:    
    plt.plot(vx)
    plt.plot(vy)
    plt.xlabel('Time [2 ns bins]')
    plt.ylabel('Voltage [uV]')
    plt.show()

vx=Addnoise(vrms,vx)
vy=Addnoise(vrms,vy)
if DISPLAY==1:
    plt.plot(vx)
    plt.plot(vy)
    plt.xlabel('Time [2 ns bins]')
    plt.ylabel('Voltage [uV]')
    plt.show()


outfile=str(sys.argv[1])+"/fake_"+str(sys.argv[2])+".txt"
f = file(outfile,"w")
for i in np.arange(len(tx)):
          print >>f,"%1.5e	%1.2e	%1.2e" % (tx[i], vx[i], vy[i] )  # .2e as in the output of the sim          
f.close()