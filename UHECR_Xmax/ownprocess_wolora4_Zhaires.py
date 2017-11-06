
#==============================================================================

# This script is friendly provided by Stijn Buitink!
# It was already adapted to perform Xmax reco for SKA-low.
# Now adapetd for GRAND Xmax reco.

#==============================================================================

import numpy as np
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
import process_func as prf


## Needs ProcessData.py to be called

def Process_Data(datadir, fileno1, outputfile, files, lowco, hico, filename1):#, nantennas):

    pickfile = open(outputfile,'w') # pickle file for saving all informations
     
    print "Number of sim: ", files
    
    Info = {}
    Info['zenith']=[]
    Info['azimuth']=[]
    Info['energy']=[]
    Info['hillas']=[]  # at the moment: look for xmax by scanning table of long. distribution, then it works also for parallel sim
    Info['longprofile']=[]
    Info['Xground']=[]
    Info['antenna_position']=[]
    Info['onskypower']=[]
    Info['filteredpower']=[]
    Info['power']=[]
    Info['power11']=[]
    Info['power21']=[]
    Info['power41']=[]
    Info['peak_time']=[]
    Info['peak_amplitude']=[]
    Info['particle_radius']=[]
    Info['energy_deposit']=[]
    Info['pol_angle']=[]
    Info['pol_angle_filt']=[]
    Info['primary']=[]

# Loop over all simulations
    
    filename = filename1
    l_sim=open(filename).read().splitlines()
    txt = np.array(l_sim)

    
    #print txt
    
    for s in np.arange(int(files)-1):
        #print s
      
      #print "WARNING: no loop on different files"
      fileno = str(txt.T[s])#int(fileno1) +s #comment out if just one file there, to plot sim
      #print lowco, hico, fileno
      print fileno
      
      
      
      longfile = '{0}/{1}/{2}.sry'.format(datadir,str(fileno),str(fileno)) #{0}/DAT{1}-999999999.long when from parallel
      #print longfile

      fname = '{0}/{1}/inp/{2}.inp'.format(datadir,str(fileno),str(fileno))
      if os.path.isfile(fname):
        steerfile = '{0}/{1}/inp/{2}.inp'.format(datadir,str(fileno),str(fileno))
      else:
        steerfile = '{0}/{1}/{2}.inp'.format(datadir,str(fileno),str(fileno))


      name1='{0}/{1}/antpos.dat'.format(datadir,str(fileno))
      #listfile = open(name1).read().splitlines()  #_2 for 160Ant
      #print name1
      
      #lines = listfile #.#readlines()
      #print lines, len(lines)
      lines=np.genfromtxt(name1)
      #print lines
      
      nantennas= len(lines)
      #print nantennas
    
      onskypower=np.zeros([nantennas,2])
      antenna_position=np.zeros([nantennas,3])
      filteredpower=np.zeros([nantennas,2])
      power=np.zeros([nantennas,2])
      power11=np.zeros([nantennas,2])
      power21=np.zeros([nantennas,2])
      power41=np.zeros([nantennas,2])
      peak_time=np.zeros([nantennas,2])
      peak_bin=np.zeros([nantennas,2])
      peak_amplitude=np.zeros([nantennas,2])
      pol_angle=np.zeros([nantennas])
      pol_angle_filt=np.zeros([nantennas])

     
      ### LORA FILE
      #lorafile = '{0}/{1}/DAT{2}.lora'.format(datadir,str(fileno).zfill(6),str(fileno).zfill(6)) # after running the Geant4 script
      
    
  #### GETTING XMAX     
  ########################       

      #longdata=np.genfromtxt(longfile, skip_header=3, usecols=(0,1,2))#skip_footer=5) ## 0:slantedheight; 1: photons, 2: positrons, 3: electrons
      #xlength=np.argmax(np.isnan(longdata[:,0]))
      #Xground=xlength*10.0 ## 
      #longprofile=np.zeros([len(longdata.T[0]),3]) #The profiles have different lengths, unlikely to exceed 400... (maybe inclined events??)

      #longprofile=np.array(  [longdata[:,0],longdata[:,1], longdata[:,2]] )
      #print "longfile ", longfile

      #m= max(longprofile[2,:]) # just get maximum
	
      #index= [i for i, j in enumerate(longprofile[2,:]) if j == m]
      #hillas=longprofile[1,index] # g/cm2# 0= slanted depth  
      
      longprofile=0.
      Xground=0.
  

	
      with open(longfile) as f:
        for line in f:
            keyword = ' Sl. depth of max. (g/cm2):' 
            if line.startswith(keyword):
                    try:
                        hillas =line.split('  ')[1]
                        print "xmax =", hillas, "g/cm2"
                    except ValueError:
                    # treat the error
                        print "error"
		
		
################	
## from Zhaires to GRAND coordinates
#zen= numpy.pi-zen
#az=numpy.pi+az  
### not for CR showers since no Raspass is used
      
      zenith=float(np.genfromtxt(re.findall("PrimaryZenAngle.*",open(steerfile,'r').read()))[1])*np.pi/180. #rad; CORSIKA coordinates
      azimuth=np.mod(np.genfromtxt(re.findall("PrimaryAzimAngle.*",open(steerfile,'r').read()))[1],360)*np.pi/180.  #rad; CORSIKA coordinates
      
      #print zenith*180./np.pi, azimuth*180./np.pi
      
      energy=np.genfromtxt(re.findall("PrimaryEnergy.*",open(steerfile,'r').read()))[1] #GeV
      
      with open(steerfile) as f:
        for line in f:
            keyword = 'PrimaryParticle' 
            if line.startswith(keyword):
                primary1=line.split(' ')[1]
                
      #primary= np.genfromtxt(re.findall("PrimaryParticle.*",open(steerfile,'r').read()))[1] # Primary RASPASS
      
      primary=0
      if primary1[0:-1] == "RASPASSelectron":
         primary =11
      if primary1[0:-1]=='Proton':
          primary =14
      if primary1[0:-1]== "Iron":
          primary =5626

      print "primary: %i" %primary

      
      
      
      #print "hardcoded primary"
      #primary=14
      #AddSpecialParticle      RASPASSProton    ./inp/RASPASSprimary Proton
      #AddSpecialParticle      RASPASSIron      ./inp/RASPASSprimary Iron
      #AddSpecialParticle      RASPASSelectron  ./inp/RASPASSprimary Electron


      
      ### read tables antenna model and initialize values
      #vt=np.loadtxt(os.environ["LOFARSOFT"] + "/data/lofar/antenna_response_model/LBA_Vout_theta.txt", skiprows=1)
      #vp=np.loadtxt(os.environ["LOFARSOFT"] + "/data/lofar/antenna_response_model/LBA_Vout_phi.txt", skiprows=1)
      #cvt = cr.hArray(vt[:, 3] + 1j * vt[:, 4])
      #cvp = cr.hArray(vp[:, 3] + 1j * vp[:, 4])
      #fstart = 10.0 * 1.e6
      #fstep = 1.0 * 1.e6
      #fn = 101
      #ttstart = 0.0
      #ttstep = 5.0
      #ttn = 19
      #pstart = 0.0
      #pstep = 10.0
      #pn = 37

      #lines = listfile.readlines()

	  
      for j in range(0,nantennas):
      	  antenna_position[j]= lines[j] #read antenna position...
	  #print antenna_position[j]
	  #print antenna_position
	  #antenna_file = lines[j].split()[3]   #... and output filename from the antenna list file
	  
	  
	  coreasfile = '{0}/{1}/a{2}.trace'.format(datadir,str(fileno), str(j)) #drop the \n from the string!
	  #print coreasfile 
	  
	  data=np.genfromtxt(coreasfile)
	  #print data
	  #data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter NOTE not needed after running split... script for Zhaires files
	  dlength=data.shape[0]
	  data[:,0]=data[:,0]*1.e-9 # ns to s
	  poldata=np.ndarray([dlength,2])
	  az_rot=azimuth #3*np.pi/2+   #conversion from CORSIKA coordinates to 0=east, pi/2=north to AUGER
	  zen_rot=zenith
	  if zen_rot ==0:
	    zen_rot+=0.001
	  XYZ=np.zeros([dlength,3])
	  XYZ[:,0]=data[:,1] # Ex
	  XYZ[:,1]=data[:,2] #Ey
	  XYZ[:,2]=data[:,3] #Ez
	  # Convert to, v, vxB, vxvxB coordinates to compute Stokes parameters and polarization angle
	  UVW=prf.GetUVW(XYZ,0,0,0,zen_rot,az_rot,(180.-27.05)*np.pi/180. )# ska  #1.1837)#LOFAR
	  Stokes=prf.stokes_parameters(UVW[:,0],UVW[:,1],fftp.hilbert(UVW[:,0]),fftp.hilbert(UVW[:,1]))
	  pol_angle[j]=prf.polarization_angle(Stokes)
	  UVWfilt=prf.FreqFilter(UVW,int(lowco),int(hico),data[1,0]-data[0,0])
	  Stokesfilt=prf.stokes_parameters(UVWfilt[:,0],UVWfilt[:,1],fftp.hilbert(UVWfilt[:,0]),fftp.hilbert(UVWfilt[:,1]))
	  pol_angle_filt[j]=prf.polarization_angle(Stokesfilt)
	  # Convert to on-sky coordinates (n, theta, phi) to prepare for application of antenna model
	  poldata[:,0] = -1.0/np.sin(zen_rot)*data[:,3] # -1/sin(theta) *z
	  poldata[:,1] = np.sin(az_rot)*data[:,2] + np.cos(az_rot)*data[:,1] # -sin(phi) *x + cos(phi)*y in coREAS 0=positive y, 1=negative x
	  spec=np.fft.rfft(poldata, axis=-2)
	  
	  
	  #Substitute Antenna model by simple Dipol antenna with flat bandpass (50-350MHz) NOTE: add antenna response model here
	  # Apply antenna model
	  tstep = data[1,0]-data[0,0] # simulation sampling
	  onskypower[j]=np.array([np.sum(poldata[:,0]*poldata[:,0]),np.sum(poldata[:,1]*poldata[:,1])])*tstep
	  freqhi = 0.5/tstep/1e6 # MHz
	  freqstep = freqhi/(dlength/2+1) # MHz
	  frequencies = np.arange(0,freqhi,freqstep)*1e6 # Hz
	  frequencies = np.arange(0,dlength/2+1)*freqstep*1e6

	  
	  ## no antenna model loaded -> cant calculate jones matrix
	  #jones_matrix = cr.hArray(complex, dimensions=(len(frequencies), 2, 2))
	  #for k,f in enumerate(frequencies):
	      #if (f>1e7 and f<1e8):
		  #cr.hGetJonesMatrix(jones_matrix[k], f, 180-(azimuth/np.pi*180),90-(zenith/np.pi*180), cvt, cvp, fstart, fstep, fn, ttstart, ttstep, ttn, pstart, pstep, pn)
		  	  
	  instr_spec=np.ndarray([dlength/2+1,2],dtype=complex)
	  
	  
	  #jm=jones_matrix.toNumpy()
	  #instr_spec[:,0] = jm[:,0,0] * spec[:,0] + jm[:,0,1] * spec[:,1]
	  #instr_spec[:,1] = jm[:,1,0] * spec[:,0] + jm[:,1,1] * spec[:,1]
	  instr_spec[:,0] = spec[:,0]
	  instr_spec[:,1] = spec[:,1]
	  
	  
	  #Apply window (filtering to desired frequency band) and reduce maximum frequency to acquire downsampled signal
	  fb = int(np.floor(lowco/freqstep))
	  lb = int(np.floor(hico/freqstep)+1)
	  window = np.zeros([1,dlength/2+1,1])
	  window[0,fb:lb+1,0]=1
	  pow0=np.abs(instr_spec[:,0])*np.abs(instr_spec[:,0]) # power = E*E
	  pow1=np.abs(instr_spec[:,1])*np.abs(instr_spec[:,1])
	  ospow0=np.abs(spec[:,0])*np.abs(spec[:,0])  # filtered time series pol0
	  ospow1=np.abs(spec[:,1])*np.abs(spec[:,1])  # filtered time series pol1
	  
	  ### NOTE: If I am right that is where you have to add the realistic noise time traces to pow0 and pow1
	  
	  
	  power[j]=np.array([np.sum(pow0[fb:lb+1]),np.sum(pow1[fb:lb+1])])/(dlength/2.)*tstep
	  filteredpower[j]=np.array([np.sum(ospow0[fb:lb+1]),np.sum(ospow1[fb:lb+1])])/(dlength/2.)*tstep
	  
	  
	  
	  ### DOWNSAMPLING NEEDED?
	  # assume that simulated time resolution is higher than LOFAR time resolution (t_step=5 ns)
	  tstep_data= 1.25e-9 #s   2*bandwitdh + extra ? NOTE: GRAND will trigger with 3ns ....
	  maxfreqbin= int(np.floor(tstep/tstep_data * dlength/2.)+1)
	  shortspec=np.array([instr_spec[0:maxfreqbin,0]*window[0,0:maxfreqbin,0],instr_spec[0:maxfreqbin,1]*window[0,0:maxfreqbin,0]])
	  filt=np.fft.irfft(shortspec, axis=-1)
	  # after downsampling, renormalize the signal!
	  dlength_new=filt.shape[1]
	  filt=filt*1.0*dlength_new/dlength
	  # to calculate the time of arrival upsample with a factor 5
	  filt_upsampled=resample(filt,5*dlength_new,axis=-1)
	  # compute hilbert enevelope
	  hilbenv=np.abs(hilbert(filt,axis=-1))
	  hilbenv_upsampled=np.abs(hilbert(filt_upsampled,axis=-1))
	  # peak_time is the bin where the maximum is located; NOT the actual time of the peak!
	  peak_bin[j]=np.argmax(hilbenv,axis=-1)
	  peak_time[j]=np.argmax(hilbenv_upsampled,axis=-1)*1e-9 #in seconds
	  peak_amplitude[j]=np.max(hilbenv_upsampled,axis=-1)
	  if (peak_amplitude[j,0]>peak_amplitude[j,1]):
	      pt=peak_bin[j,0]
	  else:
	      pt=peak_bin[j,1]
	  # for 3 different window size, the total power is calculated. The window is allowed to `wrap around', so some voodoo is needed to determine the range:
	  d=int(filt.shape[1])
	  rng=5
	  a=int(np.max([0,pt-rng]))
	  b=int(pt+rng+1)
	  c=int(np.min([d,pt+d-rng]))
	  power11[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))* 1.25e-9 #5e-9
	  rng=10
	  a=int(np.max([0,pt-rng]))
	  b=int(pt+rng+1)
	  c=int(np.min([d,pt+d-rng]))
	  power21[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))* 1.25e-9#*5e-9
	  rng=20
	  a=int(np.max([0,pt-rng]))
	  b=int(pt+rng+1)
	  c=int(np.min([d,pt+d-rng]))
	  power41[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))* 1.25e-9#*5e-9
	  
	  #HACK : Handmade Peak amplitude scanning
	  #peak_amplitude[j,0] = np.sqrt(np.amax(ospow0 + ospow1)) # squared of maximum in power, components added up in power
	  #peak_amplitude[j,1] = np.sqrt(np.amax(ospow0 + ospow1))
	  #peak_amplitude[j,0] = np.sqrt(np.amax(ospow0)) # squared of maximum in power, components added up in power
	  #peak_amplitude[j,1] = np.sqrt(np.amax(ospow1))
	  #peak_amplitude[j,0] += peak_amplitude[j,1]
	  #peak_amplitude[j,1] += np.sqrt(peak_amplitude[j,1]*peak_amplitude[j,1] + peak_amplitude[j,0]*peak_amplitude[j,0])
	  #print peak_amplitude[j,0]
	  
	  
      #### LORA FILE: has to be change from 0 to old setting if PD is included
      #pdata = np.genfromtxt(lorafile)
      particle_radius=0#5*pdata[:,0]+2.5
      energy_deposit=0#pdata[:,1]
    
      #### convert CORSIKA to AUGER coordinates (AUGER y = CORSIKA x, AUGER x = - CORSIKA y
      #temp=np.copy(antenna_position)
      #antenna_position[:,0], antenna_position[:,1], antenna_position[:,2] = -temp[:,1]/100., temp[:,0]/100., temp[:,2]/100.
    
      #temp=np.copy(peak_amplitude)
      #peak_amplitude[:,:,0], peak_amplitude[:,:,1] = -temp[:,:,1], temp[:,:,0]
    
      #temp=np.copy(peak_time)
      #peak_time[:,:,0], peak_time[:,:,1] = temp[:,:,1], temp[:,:,0]
    
      azimuth=azimuth #3*np.pi/2+  +x = east (phi=0), +y = north (phi=90)
      ###

      Info['zenith'].append(zenith) #0
      Info['azimuth'].append(azimuth) #1
      Info['energy'].append(energy)#2
      Info['hillas'].append(float(hillas)) # parallel -> tolist #3
      Info['longprofile'].append(longprofile)#4
      Info['Xground'].append(Xground)#5
      Info['antenna_position'].append(antenna_position)#6
      Info['onskypower'].append(onskypower)#7
      Info['filteredpower'].append(filteredpower)#8
      Info['power'].append(power)#9
      Info['power11'].append(power11)#10
      Info['power21'].append(power21)#11
      Info['power41'].append(power41)#12
      Info['peak_time'].append(peak_time)#13
      Info['peak_amplitude'].append(peak_amplitude)#14 hilbert pol0, pol1
      Info['particle_radius'].append(particle_radius)#15
      Info['energy_deposit'].append(energy_deposit)#16
      Info['pol_angle'].append(pol_angle)#17
      Info['pol_angle_filt'].append(pol_angle_filt)#18
      Info['primary'].append(primary)#19
    
    
    cPickle.dump(Info,pickfile)

