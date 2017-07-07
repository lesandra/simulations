# This script shall mimic all the dunctionality of the script processSimTraces.m
# this means it reads in the files a*.traces and the corresponding antpos.dat (the loop over all antennas will be performed by a shell script)
# it should do a bandpass filter for a specific frequency band, do the Hilbert envelope and save the filtered traces (in t, x,y,z and vxB-vxvxB) and the max. amplitude (total and  for each component)



import sys
from sys import argv
import filters
from optparse import OptionParser
import numpy
from numpy import *
from scipy.signal import butter, lfilter
from scipy.signal import blackman,bartlett, hamming, hanning, kaiser
from scipy import signal
from scipy.signal import hilbert
import matplotlib.pyplot as pyplot
#import plotdata 
import pylab
import os
from matplotlib.pyplot import cm 
import re
#################################################################

def rfftfreq(n, d=1.0, nyquist_domain=1):
	'''calcs frequencies for rfft, exactly as numpy.fft.rfftfreq, lacking that function in my old numpy version.
	Arguments:
	---------
		n: int 
			Number of points.
		d: float
			Sample spacing, default is set to 1.0 to return in units of sampling freq. 
		
	Returns:
	-------
		f: array of floats
			frequencies of rfft, length is n/2 + 1
	'''
	if n % 2 == 0:
		f = array([n/2 - i for i in range(n/2,-1,-1)]) / (d*n)
	else:
		f = array([(n-1)/2 + 1 - i for i in range(n/2,-1,-1)]) / (d*n)
	# if nyquist_domain is 1 you're done and return directly
	if nyquist_domain != 1:
		# if nyquist_domain even, mirror frequencies
		if (nyquist_domain % 2) == 0: f = f[::-1]
		sampling_freq = 1./d
		fmax = 0.5*sampling_freq 
		f += (nyquist_domain-1)*fmax
	return f

####################################






script = sys.argv[0]  ### name of script



## define the frequency band    
x1= float(sys.argv[3])*1.e6 #30e6 #to get MHz
x2= float(sys.argv[4])*1.e6 #80e6 #

#print 'antenna trace nr'+ sys.argv[2] + ' gets filtered to ' +sys.argv[3]+ '-' +sys.argv[4]+ 'MHz and rotated to vxb-vxvxB frame'

path =sys.argv[1] #path to first simulation which shall later be rescaled or so

txt=numpy.loadtxt(path+ 'a'+sys.argv[2]+'.trace')
# t in ns, Ex in muV/m, Ey, Ez
# NOTE: Time binning always 1ns

totalline =len(txt)

steerfile_sim = path+ '/inp/{0}.inp'.format(str(sys.argv[5]))

zen=float(numpy.genfromtxt(re.findall("PrimaryZenAngle.*",open(steerfile_sim,'r').read()))[1])*numpy.pi/180.
az=float(numpy.genfromtxt(re.findall("PrimaryAzimAngle.*",open(steerfile_sim,'r').read()))[1])*numpy.pi/180.
zen= numpy.pi-zen
az=numpy.pi+az  
#print az*180./numpy.pi, zen*180./numpy.pi


time = (txt.T[0,1]-txt.T[0,0])*1.e-9
fsample=1./(time)# *1.e-9 to get time in s
#print fsample


############ Ex

freq = rfftfreq(totalline, 1./fsample)
Ex=filters.butter_bandpass_filter(txt.T[1],x1, x2, fsample, 6)  
FFT_Ex=numpy.fft.rfft(Ex)
############ Ey

Ey=filters.butter_bandpass_filter(txt.T[2],x1, x2, fsample, 6)  
FFT_Ey=numpy.fft.rfft(Ey)

############ Ez
Ez=filters.butter_bandpass_filter(txt.T[3],x1, x2, fsample, 6)  
FFT_Ez=numpy.fft.rfft(Ez)

hexf = abs(hilbert(Ex))
heyf = abs(hilbert(Ey))
hezf = abs(hilbert(Ez))
exm = max(hexf)
eym = max(heyf)
ezm = max(hezf)

amp = sqrt(exm*exm+ eym*eym+ezm*ezm)




##### rotation to vxvxB   
#taken from oliviers scripts: orientation/direction of magnetic field
phigeo =0*numpy.pi/180.  # 182.66#; (ie pointing 2.66 degrees East from full North) # phigeo= 0 from simulations inumpyutfile % In both EVA & Zhaires, North = magnetic North
thetageo =(180.-27.05)*numpy.pi/180.  

inc=thetageo

#same script as for the antenna psoition rotation
B = numpy.array([numpy.cos(phigeo)*numpy.sin(inc), numpy.sin(phigeo)*numpy.sin(inc),numpy.cos(inc)]) #from oliviers script including phigeo
B=B/numpy.linalg.norm(B)
v = numpy.array([numpy.cos(az)*numpy.sin(zen),numpy.sin(az)*numpy.sin(zen),numpy.cos(zen)]) # or *-1: change the direction
v=v/numpy.linalg.norm(v)
   #print v
vxB = numpy.cross(v,B) #numpy.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
vxB = vxB/numpy.linalg.norm(vxB)
vxvxB = numpy.cross(v,vxB) #numpy.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
vxvxB = vxvxB/numpy.linalg.norm(vxvxB)

Ev_full= txt.T[1]* v[0] +txt.T[2]*v[1]+ txt.T[3]*v[2]
EvxB_full= txt.T[1]* vxB[0] +txt.T[2]*vxB[1]+ txt.T[3]*vxB[2]
EvxvxB_full= txt.T[1]* vxvxB[0] +txt.T[2]*vxvxB[1]+ txt.T[3]*vxvxB[2]

# bandpass filtering
Ev=filters.butter_bandpass_filter(Ev_full,x1, x2, fsample, 6)  
EvxB=filters.butter_bandpass_filter(EvxB_full,x1, x2, fsample, 6)  
EvxvxB=filters.butter_bandpass_filter(EvxvxB_full,x1, x2, fsample, 6)  


hexf2 = abs(hilbert(Ev))
heyf2 = abs(hilbert(EvxB))
hezf2 = abs(hilbert(EvxvxB))
exm2 = max(hexf2)
eym2 = max(heyf2)
ezm2 = max(hezf2)

amp2 = sqrt(exm2*exm2+ eym2*eym2+ezm2*ezm2)








######Writing to file 
name2= sys.argv[1]+"a"+sys.argv[2]+'_'+str(int(x1*1e-6)) + '-' + str(int(x2*1e-6)) + 'MHz'
name3= name2+'.dat' 
FILE2 = open(name3, "w" )


for i in range( 0, len(Ez) ):
    print >>FILE2,"%3.2f	%1.5e	%1.5e	%1.5e	%1.5e	%1.5e	%1.5e" % (txt.T[0][i], Ex[i], Ey[i], Ez[i], Ev[i], EvxB[i], EvxvxB[i] )
    
    
# in the last line of the file we wanna write the max. ampl of the hilbert envelope    
# something like: amp exm eym ezm
print >>FILE2, ""
print >>FILE2,"%1.5f	%1.5e	%1.5e	%1.5e	%1.5f	%1.5e	%1.5e	%1.5e" % (amp, exm, eym, ezm, amp2,exm2, eym2, ezm2)

FILE2.close()

########################################## PLOTTING
fig1 = pylab.figure(1)

############ EField filtered
#pylab.set(left=0.02)
#pylab.tight_layout(0.4, 0.5,1.0)
pylab.subplot(2,2, 3)
name = 'filtered EField to' + str(x1*1e-6) + '-' + str(x2*1e-6) + 'MHz' 
pylab.title(name)
pylab.xlabel('Time [ns]')
pylab.ylabel('Electric Field [muV/m]')
#pylab.xlim(40, 100)
#pylab.plot(txt.T[0],txt.T[1], '-b', label='sim.')

pylab.plot(txt.T[0],Ez, '-g', label='Ez')
pylab.plot(txt.T[0],Ey, '-b', label='Ey')
pylab.plot(txt.T[0],Ex, '-r', label='Ex')
pylab.plot(txt.T[0], heyf, '--b',label='envelope')

pylab.legend(loc='best',  prop={'size':8}) ##'upper right'
pylab.tight_layout(0.4, 0.5,1.0)

#pylab.show() # show the plot

########## Freqency
pylab.subplot(2,2, 2)
pylab.title('Frequency spectra')
pylab.xlabel('Frequency [Hz]')
pylab.ylabel('Electric Field [muV/m/Hz]')
pylab.xlim(x1-100e6, x2+100e6)

pylab.plot(freq,abs(FFT_Ez), '-g',label='Ez')
pylab.plot(freq,abs(FFT_Ey), '-b',label='Ey')
pylab.plot(freq,abs(FFT_Ex), '-r',label='Ex')
#pylab.plot(freq,(fft.fft(txt.T[1])),'-r', label='filtered')
pylab.legend(loc='best',  prop={'size':8}) ##'upper right
pylab.tight_layout(0.4, 0.5,1.0)

########## EField sim
pylab.subplot(2,2, 1)
pylab.title('sim. EField')
pylab.xlabel('Time [ns]')
pylab.ylabel('Electric Field [muV/m]')
#pylab.xlim(40, 100)
#pylab.plot(txt.T[0],txt.T[1], '-b', label='sim.')

pylab.plot(txt.T[0],txt.T[3], '-g', label='Ez')
pylab.plot(txt.T[0],txt.T[2], '-b', label='Ey')
pylab.plot(txt.T[0],txt.T[1], '-r', label='Ex')
#pylab.legend(loc='best') ##'upper right'
pylab.tight_layout(0.4, 0.5,1.0)


#pylab.show() # show the plot

fig1.savefig(name2+'.png' )