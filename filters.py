from numpy import *
from scipy.signal import butter, lfilter, bessel
from scipy.signal import blackman,bartlett, hamming, hanning, kaiser
from scipy import signal
import matplotlib.pyplot as pyp
#######################################################################
def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    #applies butterworth filter on timeseies and returns timeseries
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def bessel_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 *fs
    low = lowcut / fs
    high = highcut /fs
    b, a = bessel(order, [low,high], btype='bandpass')
    return b,a

def bessel_bandpass_filter(data, lowcut, highcut, fs, order=3):
    #applies bessel filter on timeseies and returns timeseries
    b, a = bessel_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butterworth_filter(x,GAIN,a,b,mode='td'):
    NSAMPLES = len(x)
    NCOEFF = len(a)-1
    xv=zeros(NCOEFF+1)
    yv=zeros(NCOEFF+1)
    
    y = zeros(len(x))
    for n in range(0,NSAMPLES):
        for i in range(NCOEFF):
            xv[i] = xv[i+1]
        xv[NCOEFF]=x[n] / GAIN
        for i in range(NCOEFF):
            yv[i] = yv[i+1]
        yv[NCOEFF] = sum(a*xv)
        yv[NCOEFF] += sum(b*yv)
        y[n] = yv[NCOEFF]
    if( mode == 'fd'):
        return fft.rfft(y)
    else:
        return y

def rectWindow(x, y, x1, x2):
    # x1 = start value 
    # x2 = stop value
    ind1 = where( x >= x1)[0][0]
    ind2 = where( x[ind1:] <  x2)[0][-1] + ind1
    
    ynew = zeros( len(y) )
    ynew[ind1:ind2] = y[ind1:ind2]
    return x, ynew

def rectWindow2(x, y, x1, x2):

    x = array(x, dtype=float64)
    y = array(y, dtype=complex64)
     
    if( abs( min(y) ) > abs( max(y) ) ):
        indm = where( y == min(y) )[0][0]
    else:
        indm = where( y == max(y) )[0][0]

    
    dx = x[1]-x[0]
    ind1 = int(indm - x1/dx)
    ind2 = int(indm + x2/dx)
    
    ynew = zeros( len(y), dtype=complex64 )
    ynew[ind1:ind2] = y[ind1:ind2]
    return x, ynew

# useful for windowing

def gaussian_filter(length, t0, width, dt=1.):
    t = arange(0, length)*dt
    window = exp(-0.5*((t-t0)/width)**2)
    return window
    
    
    
    
    
    