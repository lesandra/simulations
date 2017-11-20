

# This script bases on the diploma thesis of Ewa Holt (KIT, 2013) in the context of AERA/AUGER 



#called via PulseShape.py
#needs as input antenna position 1 and 2, their traces (filtered or not) in one component, their time , and the desired antenna position
# return the trace ( in x,y,z coordinate system) and the time from the desired antenna position 
#upsamples the signal 
    
    ## IMPORTANT NOTE:
    ## The interpolation of the phases includes the
    ## interpolation of the signal arrival time. A linear interpolation implies a plane radio
    ## emission wave front, which is a simplification as it is hyperbolic in shape. However, the wave front can be estimated as a plane between two simu-
    ## lated observer positions for a sufficiently dense grid of observers, as then parts of
    ## the wave front are linear on small scales.
    
    
# ATTENTION:     sometimes the simualted ring positions starts first to form a ring, sometimes first to start a ray

import sys
from sys import argv

import numpy
from numpy import *
from numpy import linalg
#from scipy.signal import butter, lfilter
#from scipy.signal import blackman,bartlett, hamming, hanning, kaiser
from scipy import signal
#from scipy.signal import hilbert
import matplotlib.pyplot as plt
#import plotdata 
import pylab
import os





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


def mag(x):
    #return numpy.sqrt(x.dot(x))
    return numpy.linalg.norm(x) #numpy.abs(x)


####################################
def getn(h):
    n = 1.+325.e-6*numpy.exp(-0.1218*h*1e-3)
    return n


def Interpolate_PulseShape(t1, trace1, x1, t2, trace2, x2, xdes, upsampling=True, ontrue=None, flow=60.e6, fhigh=200.e6): #, nrdes):
    DISPLAY=0
    
    
    #hand over time traces of one efield component -t1=time, trace1=efield- and the position x1 of the first antenna, the same for the second antenna t2,trace2, x2.
    #xdes is the desired antenna position (m) where you would like to have the efield trace in time 
    # if necessary you have to do an upsampling of the trace: upsampling=On
    # onTrue=On would give you printings to the terminal to check for correctness
    # flow= lower freq in Hz, fhigh=higher freq in Hz, not necessarily needed
    
    factor_upsampling=1
    if upsampling is not None:
        factor_upsampling=8
    
    c= 299792458.e-9 #m/ns 
    
    
    #### calculating weights: should be done with the xyz coordinates
    #weight1= mag(x2-xdes)/mag(x2-x1)#0.5 #
    #weight2= mag(xdes-x1)/mag(x2-x1)#0.5 #
    #since in star shape pattern it is mor a radial function connection the poistion of same signal as linear go for that solution. 
    #if lines ar on a line, it will give the same result as before
    weight1= mag(x2-xdes)/(mag(x2-xdes) + mag(xdes-x1))#0.5 #
    weight2= mag(xdes-x1)/(mag(x2-xdes) + mag(xdes-x1))#0.5 #

    if numpy.isinf(weight1):
        print "weight = inf"
        print x1, x2, xdes
        weight1 = 1.
        weight2 = 0.
    if weight1 != weight1:
        print 'Attention: projected positions equivalent'
        weight1 = 1.
        weight2 = 0.
    if weight1 > 1. or weight2>1:
        #print 'x1, x2, xdes, mag(x2-x1), mag(x2-xdes), mag(xdes-x1)'
        print "weight larger 1: ", weight1, weight2, x1, x2, xdes, mag(x2-x1), mag(x2-xdes), mag(xdes-x1)
        #stop
    if weight1 + weight2 >1:
        print "PulseShape_Interpolation.py: order in simulated positions. Check whether ring or ray structure formed first"
        print weight1, weight2, weight1 + weight2
    
    #print "weights " ,weight1,weight2 #, ", positions: ", x1, x2, xdes, " distances: ", mag(x2-xdes),  mag(xdes-x1), mag(x2-x1)
    
    # get refractive indey at the antenna positions
    n1 =getn(x1[2])
    n2 =getn(x2[2])
    #print "refractive indizes" +str(n1) +" " + str(n2)

    #################################################################################
    #### linearly interpolation of the phases
    
    # first antenna
    # upsampling if necessary
    f = signal.resample(trace1, len(trace1)*factor_upsampling)
    xnew = numpy.linspace(t1[0], t1[-1], len(trace1)*factor_upsampling, endpoint=False)
    xnew=xnew*1.e-9# *1.e-9 to get time in s

    fsample=1./((xnew[1]-xnew[0])) #Hz

    freq = rfftfreq(len(xnew), 1./fsample)
    FFT_Ey=numpy.fft.rfft(f)

   
    Amp = numpy.abs(FFT_Ey)
    phi =numpy.angle(FFT_Ey)

    # unwrap the phase
    l=0
    phi_unwrapped= numpy.zeros([len(phi)])
    phi_unwrapped[0]=phi[0]
    for i in range(1, len(phi)):
        if (phi[i]< phi[i-1])  and l>=0:         
            if abs(phi[i]- phi[i-1])> abs(phi[i]- phi[i-1]) + numpy.pi:
                l=l-1
        else:
            if abs(phi[i]- phi[i-1])> abs(phi[i]- phi[i-1] )- numpy.pi:
                l=l+1
        phi_unwrapped[i]= phi[i]-l*2*numpy.pi
        if ontrue is not None:
            print i, phi[i],phi[i-1],l, phi_unwrapped[i], abs(phi[i]- phi[i-1]), abs(phi[i]- phi[i-1] + numpy.pi), abs(phi[i]- phi[i-1] - numpy.pi) , l
        
       
    phi =numpy.angle(FFT_Ey)

    #############################

    ### second antenna               
    #txt2=numpy.loadtxt(path+ 'a'+str(nr2)+'.trace')
    ## t in ns, Ex in muV/m, Ey, Ez
    ## NOTE: Time binning always 1ns

    # upsampling if needed
    f2 = signal.resample(trace2, len(trace2)*factor_upsampling)
    xnew2 = numpy.linspace(t2[0], t2[-1], len(trace2)*factor_upsampling, endpoint=False)

    fsample2=1./((xnew2[1]-xnew2[0])*1.e-9)# *1.e-9 to get time in s

    freq2 = rfftfreq(len(xnew2), 1./fsample2)
    FFT_Ey=numpy.fft.rfft(f2)
    
    
    
    Amp2 = numpy.abs(FFT_Ey)
    phi2 =numpy.angle(FFT_Ey)

    #unwrap phase
    l=0
    phi2_unwrapped=  numpy.zeros([len(phi2)])
    phi2_unwrapped[0]=phi2[0]
    for i in range(1, len(phi2)):
        if (phi2[i]< phi2[i-1])  and l>=0:             
            if abs(phi2[i]- phi2[i-1])> abs(phi2[i]- phi2[i-1]) + numpy.pi:
                l=l-1
        else:
            if abs(phi2[i]- phi2[i-1])> abs(phi2[i]- phi2[i-1]) - numpy.pi:
                l=l+1
        phi2_unwrapped[i]= phi2[i]-l*2*numpy.pi

        if ontrue is not None:
            print i, phi2[i],phi2[i-1],l, phi2_unwrapped[i], abs(phi2[i]- phi2[i-1]), abs(phi2[i]- phi2[i-1] + numpy.pi), abs(phi2[i]- phi2[i-1] - numpy.pi) , l

    
    phi2 =numpy.angle(FFT_Ey)
    
    
#################    # Get the pulsh sahpe at the desired antenna position

    ### get the phase
    
    # getnumpy.zeros([len(phi2)]) the angle for the desired position    
    phides= numpy.zeros([len(phi2)])
    for i in range(0, len(phi2)-1):  
        
        phides[i]= weight1*phi_unwrapped[i] + weight2* phi2_unwrapped[i] 
 

        if ontrue is not None:
            print phides[i]
        
            
    ### re-unwrap: get -pi to +pi range back and check whether phidesis inbetwwen
    k=0        
    n=0
    for i in range(0, len(phi2)):  
        if phides[i] >= numpy.pi:
            k=k-1
            for m in range(0, len(phides)-i):
                phides[m+i]=phides[m+i]-2*numpy.pi
        else:        
            if phides[i] < (-1*numpy.pi):
                k=k+1
                for m in range(0, len(phides)-i):
                    phides[m+i]=phides[m+i]+2*numpy.pi
            else:
                continue
        if (phides[i] >= numpy.pi) or  (phides[i] < (-1*numpy.pi)):
            i=i-1
        if ontrue is not None:
            print phides[i], i
            
#### just for later plotting            
    phi2 =numpy.angle(FFT_Ey)      
    phides2=numpy.zeros([len(phi2)])
    for i in range(0, len(phi2)):  
        phides2[i]= weight1*phi_unwrapped[i] + weight2* phi2_unwrapped[i]
                
                
    #################################################################################
    #### linearly interpolation of the amplitude

    #Amp, Amp2 
    #Since the amplitude shows a continuous unipolar shape, a linear interpolation is sufficient

    Ampdes=numpy.zeros([len(phi2)])
    Ampdes2=numpy.zeros([len(phi2)]) ## for later plotting
    for i in range(0, len(phi2)):  
        Ampdes[i]= weight1*Amp[i] + weight2* Amp2[i]
        Ampdes2[i]= weight1*Amp[i] + weight2* Amp2[i]        

####### inverse FFT for the signal at the desired position
    Ampdes= Ampdes.astype(complex64)
    phides= phides.astype(complex64)
    phides2= phides2.astype(complex64)
    for i in range(0, len(Amp)): 
        Ampdes[i]= Ampdes[i]*numpy.exp(1j *phides[i] )
        
    tracedes=(numpy.fft.irfft(Ampdes))
    
    # dirty hack: still some problems with the phase, need to include an flip by adding pi once:
    #print "dirty hack for phase flip"
    import operator
    index_1, value_1 = max(enumerate(abs(trace1)), key=operator.itemgetter(1))
    index_2, value_2 = max(enumerate(abs(trace2)), key=operator.itemgetter(1))
    index_des, value_des = max(enumerate(abs(tracedes)), key=operator.itemgetter(1))
  
    tracedes=tracedes.astype(float)
    ## NOTE : Phase flip removed at 19/10/2017
    ### do I really need that flip
    #if trace1[index_1] >0. and trace2[index_2]>0:
        #if tracedes[index_des] <0.:
            #tracedes=tracedes*-1.
            #print "dirty hack for phase flip"
 
    #if trace1[index_1] <0. and trace2[index_2]<0:
        #if tracedes[index_des] >0.:
            #tracedes=tracedes*-1.
            #print "dirty hack for phase flip"
    ##print value_1, value_2, value_des, trace1[index_1], trace2[index_2], tracedes[index_des]


    
    
    # TODO: correct for phase shift from phase gradient using index_1 and index_des: xnew= xnew- abs(index_des-index_1)* binning
    if mag(x1-xdes) >200.:
        #print "another dirty hack for time shifting instead of using phase gradient (doesnt work for this large distances in between antennas)"
        timeshift= mag(x1-xdes) * n1 / c # in ns
        xnew= xnew +timeshift*1e-9 + (xnew[index_des] - xnew[index_1])
        
        #print 'timeshift', timeshift*1e-9 , xnew[index_des], xnew[index_1], (xnew[index_des] - xnew[index_1])
        

        
    ########################### PLOTTING
    
    if DISPLAY==1:

        fig1=plt.figure(1, dpi=120, facecolor='w', edgecolor='k')
        plt.subplot(311)
        plt.plot(freq, phi, 'ro-', label= "first")
        plt.plot(freq2, phi2, 'bo-', label= "second")
        plt.plot(freq2, phides, 'go--', label= "interpolated")
        #plt.plot(freq2, phi_test, 'co--', label= "real")
        plt.xlabel(r"Frequency (Hz)", fontsize=16)
        plt.ylabel(r"phase (rad)", fontsize=16)
        plt.xlim(flow, fhigh)

        #pylab.legend(loc='upper left')

        plt.subplot(312)
        ax = fig1.add_subplot(3,1,2)
        plt.plot(freq, phi_unwrapped, 'r+')
        plt.plot(freq2, phi2_unwrapped, 'bx')
        plt.plot(freq2, phides2, 'g^')
        #plt.plot(freq2, phi_test_unwrapped, 'c^')
        plt.xlabel(r"Frequency (Hz)", fontsize=16)
        plt.ylabel(r"phase (rad)", fontsize=16)
        #plt.show()
        #plt.xlim([0,0.1e8])
        #plt.xlim([1e8,2e8])
        #plt.ylim([-10,10])
        #ax.set_xscale('log')
        plt.xlim(flow, fhigh)

        plt.subplot(313)
        ax = fig1.add_subplot(3,1,3)
        plt.plot(freq, Amp, 'r+')
        plt.plot(freq2, Amp2, 'bx')
        plt.plot(freq2, Ampdes2, 'g^')
        #plt.plot(freq2, Amp_test, 'c^')
        plt.xlabel(r"Frequency (Hz)", fontsize=16)
        plt.ylabel(r"Amplitude muV/m/Hz ", fontsize=16)
        #ax.set_xscale('log')
        #ax.set_yscale('log')
        plt.ylim([1e1,10e3])
        plt.xlim(flow, fhigh)

        #print len(txt.T[2]), len(txt2.T[2]), len(numpy.real(tracedes)), len(xnew)

        #fig2=plt.figure(2)
        #plt.plot(txt.T[0], txt.T[2], 'b-', label= "first")
        #plt.plot(txt2.T[0], txt2.T[2], 'r-', label= "second")
        #plt.plot(xnew, numpy.real(tracedes), 'g--', label= "interpolated")
        #plt.plot(txt2.T[0], txt_test.T[2], 'c:', label= "real")
        ##plt.plot(txt2.T[0], Amplitude, 'b--')
        #plt.xlabel(r"time (s)", fontsize=16)
        #plt.ylabel(r"Amplitude muV/m ", fontsize=16)
        ##plt.show()
        plt.show()

    if upsampling is not None:
        return xnew[0:-1:8]*1.e9, tracedes[0:-1:8]
    else:
        return xnew*1.e9, tracedes # back to ns

# at some point save the interpolated time traces for laer analysis