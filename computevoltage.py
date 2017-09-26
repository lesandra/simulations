import sys
import numpy as np
pi = 3.1415927
import pylab as pl
import matplotlib.pyplot as plt
from scipy.fftpack import rfft, irfft, rfftfreq
from scipy.interpolate import interp1d

# Read zen and az value from input
zenith=float(sys.argv[1])
azim=float(sys.argv[2])
zen = zenith*pi/180;
az = azim*pi/180;

szen = np.sin(zen);
czen = np.cos(zen);
saz = np.sin(az);
caz = np.cos(az);

print 'Wave direction: zenith = ', zenith, ' deg, azimuth = ', azim, 'deg.'

efieldtxt='efield.txt'
fileleff='1ant10mavgnd_leff.npy'
freqscale=1 #freq*2 if h/2 and sizeant/2

# imp RLC
RL =np.array([2.946189699625225,   2.693640403293285,   2.368014134942750,   2.048063136874684,   1.763350109975618,   1.520403822949526,   1.316647804763135, 1.146675541086378,   1.004790940858134,   0.885905649139099,   0.785760827398088,   0.700893633058194,   0.628523271578350,   0.566425575932399,  0.512820527310258, 0.466279040988616])
RL=RL*100
XL =np.array([  -0.398164981757440,  -0.908417628420895,  -1.223336201352510,  -1.396290370230734,  -1.476701296666265,  -1.499861221583199,  -1.488751816960041,  -1.457793478766247,  -1.415898297104643,  -1.368535029964790,  -1.319038439288611,  -1.269420739673441,  -1.220871947347262,  -1.174069331311762,  -1.129370040642129, -1.086931910885215])
XL=XL*100

#===========================================================================================================
def get_voltage(amplitudet=None, amplitudep=None, Fs=1.0):
#===========================================================================================================


    freq,realimp,reactance,theta,phi,lefftheta,leffphi,phasetheta,phasephi=np.load(fileleff)
    nfreq=len(freq[:,0])
    f=np.zeros(nfreq)
    RA=np.zeros(nfreq)
    XA=np.zeros(nfreq)
    ltr=np.zeros(nfreq)
    lta=np.zeros(nfreq)
    lpr=np.zeros(nfreq)
    lpa=np.zeros(nfreq)
    roundazimuth=round(azim/10)*10+round((azim-10*round(azim/10))/5)*5
    if roundazimuth>=95 and roundazimuth<=180:
        roundazimuth=180-roundazimuth
    if roundazimuth>=185 and roundazimuth<=270:
        roundazimuth=roundazimuth-180
    if roundazimuth>=275 and roundazimuth<=360:
        roundazimuth=360-roundazimuth
    for i in range(nfreq):
        f[i]=freq[i,0]*freqscale
        RA[i]=realimp[i,0]
        XA[i]=reactance[i,0]
        indtheta=np.nonzero(theta[i,:]==round(zenith))[0]
        indphi=np.nonzero(phi[i,:]==roundazimuth)[0]
        indcom=np.intersect1d(indtheta,indphi)
        ltr[i]=lefftheta[i,indcom]
       	lta[i]=phasetheta[i,indcom]*np.pi/180
       	lpr[i]=leffphi[i,indcom]
       	lpa[i]=phasephi[i,indcom]*np.pi/180
        Rlefft=ltr[i]*np.cos(lta[i])
        Xlefft=ltr[i]*np.sin(lta[i])
        Rleffp=lpr[i]*np.cos(lpa[i])
        Xleffp=lpr[i]*np.sin(lpa[i])
        Rleqt=((Rlefft*RL[i]-Xlefft*XL[i])*(RA[i]+RL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
        Xleqt=((Rlefft*RL[i]+Xlefft*XL[i])*(XA[i]+XL[i]) + (Rlefft*XL[i]+Xlefft*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
        ltr[i]=np.sqrt(Rleqt**2+Xleqt**2)
        print(Rleqt,Xleqt,ltr[i])
        lta[i]=np.arccos(Rleqt/ltr[i])
        Rleqp=((Rleffp*RL[i]-Xleffp*XL[i])*(RA[i]+RL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(XA[i]+XL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
        Xleqp=((Rleffp*RL[i]+Xleffp*XL[i])*(XA[i]+XL[i]) + (Rleffp*XL[i]+Xleffp*RL[i])*(RA[i]+RL[i])) / ((RA[i]+RL[i])**2+(XA[i]+XL[i])**2)
        lpr[i]=np.sqrt(Rleqp**2+Xleqp**2)
        print(Rleqp,lpr[i])
        lpa[i]=np.arccos(Rleqp/lpr[i])

    #phases are not unwrap! so:
    for i in range(1,nfreq):
        while lpa[i]-lpa[i-1]<-180*np.pi/180:
            lpa[i]=lpa[i]+360*np.pi/180
        while lpa[i]-lpa[i-1]>180*np.pi/180:
            lpa[i]=lpa[i]-360*np.pi/180
        while lta[i]-lta[i-1]<-180*np.pi/180:
            lta[i]=lta[i]+360*np.pi/180
        while lta[i]-lta[i-1]>180*np.pi/180:
            lta[i]=lta[i]-360*np.pi/180

    print(round(zenith),roundazimuth,f,ltr,lta,lpr,lpa,Fs)

    fmin=f[0]
    fmax=f[-1]
    f=f*1e6

    nf  = int(2**np.floor(np.log(len(amplitudet))/np.log(2)))
    while Fs/nf > fmin*1e6:   # <== Make sure that the DFT resolution is at least fmin.
        nf *= 2
    F = rfftfreq(nf)*Fs


    modulust = interp1d(f, ltr, bounds_error=False, fill_value=0.0)(F)
    phaset   = interp1d(f, lta, bounds_error=False, fill_value=0.0)(F)
    modulusp = interp1d(f, lpr, bounds_error=False, fill_value=0.0)(F)
    phasep   = interp1d(f, lpa, bounds_error=False, fill_value=0.0)(F)

    #print(F,modulust)

    phaset -= phaset[0] # Switch the phase origin to be consistent with a real signal.
    phasep -= phasep[0] # Switch the phase origin to be consistent with a real signal.

    A = rfft(amplitudet, nf)
    ct = np.cos(phaset)
    st = np.sin(phaset)
    B = np.zeros(A.shape)
    B[1:-1:2] = modulust[1:-1:2]*(A[1:-1:2]*ct[1:-1:2]-A[2:-1:2]*st[2:-1:2])
    B[2:-1:2] = modulust[2:-1:2]*(A[1:-1:2]*st[1:-1:2]+A[2:-1:2]*ct[2:-1:2])
    B[0]  = A[0]*modulust[0]
    B[-1] = A[-1]*modulust[-1]

    C = rfft(amplitudep, nf)
    cp = np.cos(phasep)
    sp = np.sin(phasep)
    D = np.zeros(C.shape)
    D[1:-1:2] = modulusp[1:-1:2]*(C[1:-1:2]*cp[1:-1:2]-C[2:-1:2]*sp[2:-1:2])
    D[2:-1:2] = modulusp[2:-1:2]*(C[1:-1:2]*sp[1:-1:2]+C[2:-1:2]*cp[2:-1:2])
    D[0]  = C[0]*modulusp[0]
    D[-1] = C[-1]*modulusp[-1]


    vt=irfft(B)
    vp=irfft(D)

    return(vt,vp)






#===========================================================================================================
# Compute the time dependent voltage
#===========================================================================================================

# Model the input signal.
time1, Ex, Ey,Ez = np.loadtxt(efieldtxt,delimiter='  ',usecols=(1,2,3,4),unpack=True)
time1 = (time1-time1[0])*1e-6
print 'Efield signal length: ', time1[-1]*1e9 , 'ns.'

delt = time1[1]-time1[0];
Fs = 1/delt


Et = czen*(caz*Ex+saz*Ey)-szen*Ez
Ep = -saz*Ex+caz*Ey


print 'Now computing EW arm response...'
# Compute the output voltage.
vt,vp   = get_voltage(Et, Ep, Fs)
timet     = np.arange(0, len(vt))/Fs
timep     = np.arange(0, len(vp))/Fs

# Compute the total output voltage.
voltage = vt+vp

pl.savetxt('ew_out.txt', (timet, voltage,vt,vp))



'''print 'Now computing NS arm response...'
# Compute the output voltage.
vt,vp   = get_voltage(Et, Ep, Fs)
timet     = np.arange(0, len(vt))/Fs
timep     = np.arange(0, len(vp))/Fs

# Compute the total output voltage.
voltage = vt+vp

pl.savetxt('ns_out.txt', (timet, voltage,vt,vp))
'''

#get_voltage(1,1,1)


#plt.plot(time1,Et)
#plt.xlabel('Time [sec]')
#plt.ylabel('Etheta [muV]')
#plt.show()
#plt.plot(time1,Ep)
#plt.xlabel('Time [sec]')
#plt.ylabel('Ephi [muV/m]')
#plt.show()
'''plt.plot(timet,vt)
plt.xlabel('Time [sec]')
plt.ylabel('Vtheta [muV]')
plt.show()
plt.plot(timep,vp)
plt.xlabel('Time [sec]')
plt.ylabel('Vphi [muV]')
plt.show()'''
