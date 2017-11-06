import numpy as np

def FreqFilter(data,lowfreq,hifreq,tstep):
   dlength=data.shape[0]
   spec=np.fft.rfft(data, axis=-2)
   freqhi = 0.5/tstep/1e6 # MHz
   freqstep = freqhi/(dlength/2+1) # MHz
   fb = int(np.floor(lowfreq/freqstep))
   lb = int(np.floor(hifreq/freqstep)+1)
   window = np.zeros([dlength/2+1,1])
   window[fb:lb+1,:]=1
   maxfreqbin= int(np.floor(tstep/5e-9 * dlength/2.)+1)
   shortspec=spec[0:maxfreqbin,:]*window[0:maxfreqbin,:]
   filt=np.fft.irfft(shortspec, axis=-2)
   return filt

def GetUVW(pos, cx, cy, cz, zen, az, Binc):
   relpos = pos-np.array([cx,cy,cz])
   B = np.array([0,np.cos(Binc),-np.sin(Binc)]) 
   v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
   vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
   return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T

def GetAlpha(zen,az,Binc):
   B = np.array([0,np.cos(Binc),-np.sin(Binc)])
   v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
   vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
   return np.arccos(np.inner(np.asarray(B), np.asarray(v)) / (np.linalg.norm(B) * np.linalg.norm(v)))

def stokes_parameters(x, y, hx, hy):
    """Stokes parameters given timeseries *x*, *y* in two orthogonal
polarisations and their Hilbert transforms *hx* and *hy*. The *x* and
*y* axis are along vxB and vxvxB respectively.
    """
    n = x.shape[0]

    I = (1./n) * np.sum(x*x + hx*hx + y*y + hy*hy)
    Q = (1./n) * np.sum(x*x + hx*hx - y*y - hy*hy)
    U = (2./n) * np.sum(x*y + hx*hy)
    V = (2./n) * np.sum(hx*y - x*hy)

    return np.array([I, Q, U, V])

def polarization_angle(S):

    Q = S[1]
    U = S[2]

    psi = (1./2.)*np.arctan2(U, Q)

    return psi
