
#==============================================================================

# Creates Zhaires input file for antenna positions in star shape pattern on a
# mountain slope, for zenith/deg, azimuth/deg, energy/EeV, name
# Usage: python GetAntennas_mountainSlope.py 83 40 1.0 CR_0

#==============================================================================


import numpy as np
from optparse import OptionParser
from Tkinter import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.interpolate as intp
import os

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.mplot3d import Axes3D


def mag(x):
    y=0
    for i in range(0,len(x)):
        y=y+float(x[i])*float(x[i])
        #print i , float(x[i])*float(x[i]), y
    return float(np.sqrt(float(y)))


# NOTE: this functions returns other order of coordinates since the 3rd component should be zero in shower coordinates
def GetUVW(pos, cx, cy, cz, zen, az, phigeo, bfieldangle):
   
   relpos = pos-np.array([cx,cy,cz])
   inc=bfieldangle

   B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) #from oliviers script including phigeo
   B=B/np.linalg.norm(B)
   v = np.array([np.cos(az)*np.sin(zen),np.sin(az)*np.sin(zen),np.cos(zen)]) # or *-1: change the direction
   v=v/np.linalg.norm(v)
   #print v
   vxB = np.cross(v,B) #np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.cross(v,vxB) #np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
   vxvxB = vxvxB/np.linalg.norm(vxvxB)

   return np.array([np.dot(v,relpos),np.dot(vxB,relpos),np.dot(vxvxB,relpos)]).T # vector dot


def computeAntennaPos(zenith, azimuth, alpha, energy, dist_fromxmax, name):
  # Here (x,y,z) = (Northing,Westing, Up). Alt ref at sea level.
  # Use Zhaires (theta, phi) convention: defined wrt direction of origin
  # inside Zhaires a conversion to zaxis down is made
    degtorad= np.pi/180.
    radtodeg= 180./np.pi
  
    DISPLAY = 0
    PRINT = 0
    MAP = 1

    #name="CR_0"
    simdir="./CR-Sim/"
    fname = simdir+ name+'.inp'#"AntennaPositions.inp"
    primary= "Proton" #"Iron" #
    print "primary", primary
    #energy =0.1 #EeV
    
    ##diatnce antenna array to xmax, should be a parameter handed over
    #dist_fromxmax=100000. #m at the moment for 10^18eV proton  #NOTE will change with zenith, get this with density function and a rough value for xmax at this energy
    ##print "dist_fromxmax fixed to ", dist_fromxmax, "m" 
    
    n_ref=1.0003 # usually depends on density at Xmax, here neglected
    theta_ch=np.arccos(1./n_ref)
    #print "refractive index: " +str(n_ref) + " Cherenkovangle: " +str(theta_ch*radtodeg)


    # mountain slope 
    #alpha =0.*degtorad #10

    #Ulastai height
    h_ulas= 0. #2650.#m above sealevel. does Zhaires refer to observerlevel?
    #print 'height of Ulastai: ' + str(h_ulas)






    ## GRAND angles
    #zen2=83.*degtorad
    #az2=0.*degtorad # 0=coming from north

    ## ATTENTION use ZHAireS angle convention
    ### from GRAND to ZHAireS
    #zen_rad=np.pi-zen2 
    #az_rad=az2+np.pi
    #zen2=(180.-83.)*degtorad
    #az2=(180.+0.)*degtorad 
    zen2=float(zenith)*degtorad
    az2=float(azimuth)*degtorad 

    ## from ZHAireS to GRAND #### nt needed: since CR showers and no raspass needed
    ##taken from oliviers scripts: orientation/direction of magnetic field
    #phigeo =0*degtorad  # 182.66#; (ie pointing 2.66 degrees East from full North) # phigeo= 0 from simulations inputfile % In both EVA & Zhaires, North = magnetic North
    #thetageo =(180.-27.05)*degtorad # 152.95*degtorad #27.05*degtorad #; (pointing down)-62.95

    zen_rad=zen2
    az_rad=az2
    #zen_rad=np.pi-zen2
    #az_rad=np.pi+az2
    #print "Zhaires convention theta, phi " + str(zen2*radtodeg) + ', ' +str(az2*radtodeg)
    #print "GRAND convention theta, phi " + str(zen_rad*radtodeg) + ', ' +str(az_rad*radtodeg)




    ### Crate star shape in GRAND coordinates
    # Direction where Earth mag field points to @ Ulastai
    #az_B = 2.66*deg2rad
    az_B = 0.*degtorad  # North = Magnetic North
    zen_B = 152.95*np.pi/180. #Direction where we point to
    B = np.array([np.cos(az_B)*np.sin(zen_B),np.sin(az_B)*np.sin(zen_B),np.cos(zen_B)]) #in LOFAR coordinates

    v = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])
    v = v/np.linalg.norm(v)
    #print "shower diection " , v
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    vxvxB = vxvxB/np.linalg.norm(vxvxB)
    #print "vxB:",vxB
    #print "vxvxB:",vxvxB


    ## in principle, antenna array center should be at a fixed height on a slanted surface, just the star shape pattern should be orientated with the shower direction  
    #firstInter= (dist_toxmax+dist_fromxmax)

    ## this should be a fixes xyz on the slanted surface
    #plane0 = firstInter* np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])+decay
    #print 'Center of antenna plane:',plane0

    max_ang = 2.*degtorad  # Most distant antenans are 2degs from axis 
    d = dist_fromxmax*np.tan(max_ang)
    step = d/20.   # 20 antennas per arm (160 in total)
    #print 'distance (shower coord) between antennas in m:', step


    ####### projection on mountain
    ####plane0 = (dist_toxmax+dist_fromxmax)*np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])+decay

    # define mountan slope as plae which is always facing the shower
    # r=r0 + l*u + m*v, with u*u=v*v=1, u perp to v
    u=np.array([np.cos(az_rad+0.5*np.pi), np.sin(az_rad+0.5*np.pi),0.]) # always z=0, vector should be perp to shower axis = az_rad +0.5*pi
    #u=u/mag(u)
    #print "u:", u
    v=np.array([np.sin(0.5*np.pi-alpha)*np.cos(az_rad+np.pi), np.sin(0.5*np.pi-alpha)*np.sin(az_rad+np.pi), np.cos(0.5*np.pi-alpha) ]) # describes the mountain slope and should be perp to u,0.5*np.pi-alpha to account for mountain slope, az_rad+np.pi because slope pointing towards inverse dir of shower 
    
    n=np.cross(u,v)

    #v=v/mag(v)
    #print "v:", v
    #print "n:", n

    # translation , c includes already (0,0,h)
    y= -dist_fromxmax * np.sin(theta_ch)*np.cos(zen_rad-alpha+theta_ch)/( (np.sin(zen_rad-alpha))**2. -(np.cos(theta_ch))**2.  ) # major axis od ellipse
    h= np.sin(alpha) *y # height of antenna array center on the mountain, if complete shower should be on the mountain
    #print "y:", y , "h:", h
    c= y*np.cos(zen_rad) # along mountain position of shower core
    #NOTE: z component of alpha flipped 
    a =np.array([np.sin(zen_rad)*np.cos(az_rad), np.sin(zen_rad)*np.sin(az_rad), np.cos(zen_rad)]) # shower direction
    a=a/mag(a)
    a_op = np.array([a[0], a[1],0]) #ortogonal projection on ground to shift for shower axis
    a_op=a_op/mag(a_op)
    d= h* np.tan(zen_rad)# shift length to respect shower axis
    
    r0= np.array([0,0,h_ulas+h]) + d*a_op# +c* v #works like plane0, gives you the positions vector of the projection
    #print "a:",a, "shower direction", np.cross(vxB,vxvxB)
    #print "c:", c, "h:",h, "sin(pi-zen):",np.sin(0.5*np.pi-zen_rad), "y:",y, "sin(alpha):",np.sin(alpha)
    #print "r0:",r0

    
    
    #### star shape pattern in xyz
    xyz1=np.zeros([160,3]) 
    xyz=np.zeros([160,3]) 
    xyz2=np.zeros([160,3]) # back trafo in vxvxB
    #xyz3=np.zeros([160,3]) # back trafo in xyz
    ## Projection on mountain (loop over xi)
    # P(x) = (xi*u)*u + (xi*v)v
    
    #for i in np.arange(1,16):
    for i in np.arange(1,21):   #AZ setup
        for j in np.arange(8):
            #xyz0 = i*step*(np.cos(j/4.0*np.pi)*np.array([1,0,0])+np.sin(j/4.0*np.pi)*np.array([0,1,0])) # pattern in shower coordinates
            
    
            xyz0 = i*step*(np.cos(float(j/4.0)*np.pi)*vxB+np.sin(float(j/4.0)*np.pi)*vxvxB) # pattern in xyz     xyz0 # z*v=0, since z=0
            xyz1[(i-1)*8+j]=xyz0 # original starshape produced
            

            
            # intersection of shower and mountain plane
            b=-np.dot(n,xyz1[(i-1)*8+j])/ np.dot(n, a)
            xyz[(i-1)*8+j]= xyz1[(i-1)*8+j] +b*a +r0 # projected
            

  
    #file.close()
    for i in np.arange(160):
        xyz2[i]=GetUVW(xyz[i], r0[0], r0[1], r0[2], zen_rad, az_rad, az_B, zen_B)# as used later to fo in vxB
    

    
    shower=np.zeros([200,3])
    mount_u=np.zeros([200,3])
    mount_v=np.zeros([200,3])
    for i in np.arange(0,200):    
        shower[i]= (i-100)*100 *a +r0
        mount_u[i]= (i-100)*100 *u +r0
        mount_v[i]= (i-100)*100 *v +r0
    
    
    
    # ATTENTION back rotation has to work!    
    ###check: back to shower coordinates 
    #for i in np.arange(0,160):   #AZ setup
        ##r0[2]=xyz[i,2]
        #xyz[i,:]= GetUVW(xyz[i,:],r0[0],r0[1],r0[2], zen_rad, az_rad, az_B, zen_B)
        ##xyz[i]=xyz[i] -r0
        ##xyz[i]= np.dot(np.dot(xyz[i],u),u) +np.dot(np.dot(xyz[i],v),v) 

        
        ##xyz[i]=xyz[i]-r0
        ##xyz[i]=xyz[i]
        
    if PRINT:
        print 'produce input file ....', fname
        
        
        file= open(fname, 'w')

        #name=CR_1
        #energy =0.1 #EeV
        #zenith =83.0 # in deg
        import random
        seed=random.uniform(0, 1)

        task='TaskName '+str(name)+ '\n'
        file.write(task)
        prim='PrimaryParticle '+str(primary) + '\n'
        file.write(prim)
        file.write('PrimaryEnergy {0} EeV\n'.format(float(energy)))
        file.write('PrimaryZenAngle {0:2.1f} deg\n'.format(zenith))
        file.write('PrimaryAzimAngle {0:2.1f} deg Magnetic\n'.format(azimuth))
        file.write('RandomSeed {0:1.3f}\n'.format(seed)) # seems to be  not needed
        file.write('# mountain slope {0:2.1f} deg\n'.format(alpha*radtodeg))
        file.write('\n')
        file.write('PropagatePrimary On\n')
        file.write('\n')
        file.write('GeomagneticField On\n')
        file.write('\n')
        file.write('SetGlobal RASPASSTimeShift 0.00\n')
        file.write('SetGlobal RASPASSDistance 0.00\n')
        file.write('\n')
        file.write('Thinning   1e-4  Relative #Thinning is adjusted as per ZHAireS recommendations\n')
        file.write('ThinningWFactor 0.06\n')
        file.write('\n')
        #file.write('# Threshold energies. Particles are not followed below these energies.
        #file.write('#GammaCutEnergy     200 keV
        #file.write('#ElectronCutEnergy  200 keV
        #file.write('#MuonCutEnergy        1 MeV
        #file.write('#MesonCutEnergy     1.5 MeV
        #file.write('#NuclCutEnergy      150 MeV
        file.write('ElectronCutEnergy 10 MeV\n') # should be discussed
        file.write('ElectronRoughCut 10 MeV\n')
        file.write('GammaCutEnergy 10 MeV\n')
        file.write('GammaRoughCut 10 MeV\n')
        file.write('\n')
        file.write('ForceLowEDecay Never\n')
        file.write('ForceLowEAnnihilation Never\n')
        file.write('\n')
        file.write('#Location is Ulastai (TREND site 42circ55N 86.68E). where Bgeo is the same as the one used in your GRAND sims:\n')
        file.write('#phigeo = 182.66; (ie pointing 2.66 degrees East from full North)\n')
        file.write('#thetageo = 27.05; (pointing down)\n')
        file.write('#bgeo = 56.5; %muT\n')
        file.write('#The ground altitude is 2650m.\n')
        file.write('#Linsley\n')
        file.write('Atmosphere 1\n')
        file.write('AddSite Ulastai 42.55 deg 86.68 deg 2650 m\n')
        file.write('Site Ulastai\n')
        file.write('Date 1985 10 26\n')
        file.write('GeomagneticField On\n')
        file.write('GeomagneticField 56.5 uT 63.18 deg 2.72 deg\n')
        file.write('GroundAltitude 2650 m\n') # default =0m
        file.write('\n')
        file.write('# Save PerShowerData\n')
        file.write('# Parameters of output files.\n')
        file.write('# you wont be needing ground or longitudinal particles\n')
        file.write('\n')
        file.write('PerShowerData Full\n')
        file.write('SaveNotInFile lgtpcles All\n')
        file.write('SaveNotInFile grdpcles All\n')
        file.write('RLimsFile grdpcles 0.000001 m 10 km\n')
        file.write('ResamplingRatio 100\n')
        file.write('RLimsTables 10 m 10 km\n')
        file.write('ELimsTables 10 MeV 1 TeV\n')
        file.write('ObservingLevels 500  0 m 25000 m\n')
        file.write('\n')
        file.write('# All other parameters will be assigned a default value if not set.\n')
        file.write('TotalShowers 1\n')
        file.write('RunsPerProcess Infinite\n')
        file.write('ShowersPerRun 1\n')
        file.write('\n')
        file.write('# In RASPASS all tables are along the shower axis, using planes perpendicular to the axis\n')
        file.write('# The table X axis is in km, or in g/cm2 if Opt a is given.\n')
        file.write('# Opt s supresses header, for gnuplot paste command to work propperly.\n')
        file.write('# Current limitation is that it only works in run time. If tables are exported using AiresExport,\n')
        file.write('# the table content will be ok but the X axis will be wrong\n')
        file.write('\n')
        file.write('#Tables\n')
        file.write('All particles\n')
        file.write('ExportTable 1293 Opt s\n')
        file.write('ExportTable 1293 Opt as\n')
        file.write('ExportTable 1793 Opt s\n')
        file.write('ExportTable 1793 Opt as\n')
        file.write('ExportTable 7293 Opt s\n')
        file.write('ExportTable 7793 Opt s\n')
        file.write('ExportTable 7993 Opt s\n')
        file.write('\n')
        file.write('#e+/e-\n')
        file.write('ExportTable 1205 Opt s\n')
        file.write('ExportTable 1205 Opt as\n')
        file.write('ExportTable 1705 Opt s\n')
        file.write('ExportTable 1705 Opt as\n')
        file.write('ExportTable 7205 Opt s\n')
        file.write('ExportTable 7705 Opt s\n')
        file.write('ExportTable 7905 Opt s\n')
        file.write('\n')
        file.write('#gammas\n')
        file.write('ExportTable 1001 Opt s\n')
        file.write('ExportTable 1001 Opt as\n')
        file.write('ExportTable 1501 Opt s\n')
        file.write('ExportTable 1501 Opt as\n')
        file.write('ExportTable 7001 Opt s\n')
        file.write('ExportTable 7501 Opt s\n')
        file.write('ExportTable 7801 Opt s\n')
        file.write('\n')
        file.write('#mu+/mu-\n')
        file.write('ExportTable 1207 Opt s\n')
        file.write('ExportTable 1207 Opt as\n')
        file.write('ExportTable 1707 Opt s\n')
        file.write('ExportTable 1707 Opt as\n')
        file.write('ExportTable 7207 Opt s\n')
        file.write('ExportTable 7707 Opt s\n')
        file.write('ExportTable 7907 Opt s\n')
        file.write('\n')
        file.write('#pi+/pi-\n')
        file.write('ExportTable 1211 Opt s\n')
        file.write('ExportTable 1211 Opt as\n')
        file.write('ExportTable 1711 Opt s\n')
        file.write('ExportTable 1711 Opt as\n')
        file.write('\n')
        file.write('#k+/k-\n')
        file.write('ExportTable 1213 Opt s\n')
        file.write('ExportTable 1213 Opt as\n')
        file.write('ExportTable 1713 Opt s\n')
        file.write('ExportTable 1713 Opt as\n')
        file.write('\n')
        file.write('#neutrons\n')
        file.write('ExportTable 1021 Opt s\n')
        file.write('ExportTable 1021 Opt as\n')
        file.write('ExportTable 1521 Opt s\n')
        file.write('ExportTable 1521 Opt as\n')
        file.write('\n')
        file.write('#protons\n')
        file.write('ExportTable 1022 Opt s\n')
        file.write('ExportTable 1022 Opt as\n')
        file.write('ExportTable 1522 Opt s\n')
        file.write('ExportTable 1522 Opt as\n')
        file.write('\n')
        file.write('#antiprotons\n')
        file.write('ExportTable 1023 Opt s\n')
        file.write('ExportTable 1023 Opt as\n')
        file.write('ExportTable 1523 Opt s\n')
        file.write('ExportTable 1523 Opt as\n')
        file.write('\n')
        file.write('#nuclei\n')
        file.write('ExportTable 1041 Opt s\n')
        file.write('ExportTable 1041 Opt as\n')
        file.write('ExportTable 1541 Opt s\n')
        file.write('ExportTable 1541 Opt as\n')
        file.write('\n')
        file.write('#och\n')
        file.write('ExportTable 1591 Opt s\n')
        file.write('ExportTable 1591 Opt as\n')
        file.write('ExportTable 7091 Opt s\n')
        file.write('ExportTable 7591 Opt s\n')
        file.write('ExportTable 7891 Opt s\n')
        file.write('\n')
        file.write('#on\n')
        file.write('ExportTable 1592 Opt s\n')
        file.write('ExportTable 1592 Opt as\n')
        file.write('ExportTable 7092 Opt s\n')
        file.write('ExportTable 7592 Opt s\n')
        file.write('ExportTable 7892 Opt s\n')
        file.write('#\n')
        file.write('# ZHAireS v0.28r22\n')
        file.write('#\n')
        file.write('\n')
        file.write('#ZHAiresS Input\n')
        file.write('#ZHAireS On/Off (Default:Off)\n')
        file.write('#Enable the calculation of the electric field. as per the ZHS formalism\n')
        file.write('ZHAireS On\n')
        file.write('\n')
        file.write('#FresnelFreq On/Off (Default:Off)\n')
        file.write('#Perform the calculation in the frequency domain. using the fresnel aproximation\n')
        file.write('FresnelFreq Off\n')
        file.write('\n')
        file.write('#FresnelTime On (Default:Off)\n')
        file.write('#Perform the calculation in the time domain. using the fresnel aproximation\n')
        file.write('FresnelTime On\n')
        file.write('\n')
        file.write('#RefractionIndex n (Default:1.000325. Variable) (note that we are not using the default on purpose)\n')
        file.write('# Refraction Index in EVA?\n')
        file.write('#RefractionIndex 1.00035\n')
        file.write('#ConstRefrIndex\n')
        file.write('\n')
        file.write('#TimeDomainBin t (Default:0.5 ns)\n')
        file.write('#The widht of the bin to be used in time domain calculations\n')
        file.write('TimeDomainBin 1 ns\n')
        file.write('AntennaTimeMin -100 ns\n')
        file.write('AntennaTimeMax 500 ns\n')
        file.write('\n')
        file.write('\n')
        file.write('#\n')
        file.write('# End of the fixed part\n')
        file.write('# \n')
        file.write('# This input file has been generated using ProduceInputFile of the AJM suite.\n')
        file.write('# \n')
    
	
        for i in np.arange(160):
            file.write("AddAntenna  {0:11.3f} {1:11.3f} {2:11.3f}\n".format(xyz[i,0],xyz[i,1],xyz[i,2]))
        
        file.close()

    if DISPLAY: 
        fig1=plt.figure(1, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k') 
        #title="zen_G="+str(zen_rad*radtodeg) + " az_G="+str(az_rad*radtodeg) + " slope=" +str(alpha*radtodeg)
        #fig1.suptitle(title, fontsize=16)

        #plt.subplot(221)

        #plt.scatter(xyz[:,0],xyz[:,1], c='red') 
        #plt.scatter(shower[:,0],shower[:,1], c='blue') 

        #plt.xlabel(r"x")
        #plt.ylabel(r"y") #, fontsize=16
        #plt.axis('equal')



        ######
        #plt.subplot(222)
        #ax=plt.plot()

        #plt.scatter(xyz[:,0],xyz[:,2], c='red', label="antennas")
        #plt.scatter(shower[:,0],shower[:,2], c='blue', label="shower axis") 
        #plt.xlabel(r"x")
        #plt.ylabel(r"z") #, fontsize=16
        #plt.axis('equal')
        #plt.plot([min(shower[:,0]), max(shower[:,0])], [h_ulas, h_ulas], 'k-', label="Ulastai height")
        #plt.legend(loc='best')
        
        ######
        plt.subplot(121)
        ax=plt.plot()
        plt.scatter(xyz1[:,1],xyz1[:,2], c='green') # produce star shape pattern
        plt.scatter(xyz2[:,1],xyz2[:,2], c='yellow') # back transformed star shape pattern
        #plt.scatter(shower[:,1],shower[:,2], c='blue') 
        plt.xlabel(r"vxB - back trafo ")
        plt.ylabel(r"vxvxB - back trafo") #, fontsize=16
        plt.axis('equal')

        ######
        plt.subplot(122)
        #ax = fig1.add_subplot(224, projection='3d')
        ax = fig1.add_subplot(122, projection='3d')

        ax.scatter(xyz1[:,0],xyz1[:,1],xyz1[:,2], c='green') # produced star shape
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2], c='red')  # projected
        ax.scatter(shower[:,0],shower[:,1],shower[:,2], c='blue')  # showre
        ax.scatter(mount_u[:,0],mount_u[:,1],mount_u[:,2], c='black')  # mountain
        ax.scatter(mount_v[:,0],mount_v[:,1],mount_v[:,2], c='black')  # mountain
        #ax.scatter(xyz3[:,0],xyz3[:,1],xyz3[:,2], c='yellow') # backtransformed

        plt.axis('equal')
        plt.xlabel(r"x ")
        plt.ylabel(r"y") #, fontsize=16
        #plt.zlabel(r"z") #, fontsize=16

        plt.show()  
        
        
    if MAP: # build the antenna map for the script, should be centered around (000) , so dont shift by r0        
        u=u/mag(u)
        v=v/mag(v)

        # r0 center of array
        distance = 1000 # m
        mu=distance/mag(u) # 20times
        mv=distance/mag(v)  # 20times
        print u, v
        
        positions_grid=np.zeros([441,3]) 
        for i in np.arange(0, 21):
            for j in np.arange(0,21):
                positions_grid[(i)*21+(j)]=(i-10)*(mu*u) + (j-10)*(mv*v) +r0
             
        fig2=plt.figure(2, figsize=(12, 10), dpi=120, facecolor='w', edgecolor='k') 
        plt.subplot(111)
        ax = fig2.add_subplot(111, projection='3d')

        ax.scatter(positions_grid[:,0],positions_grid[:,1],positions_grid[:,2], c='green') # produced star shape
        ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2], c='red')  # projected

        plt.axis('equal')
        plt.xlabel(r"x ")
        plt.ylabel(r"y")
        
        fname='GRAND_antenna.list'
        print 'produce antenna file ....', fname
                
        file=open(fname, 'w')
        
        for i in np.arange(441):
            file.write("{0:11.3f} {1:11.3f} {2:11.3f}\n".format(positions_grid[i,0],positions_grid[i,1],positions_grid[i,2]))
        
        file.close()
        
  
                

if __name__ == '__main__':
    
  if np.size(sys.argv) <=2:
    print "Arguments = zen (deg) az (deg) slope (deg) energy (EeV) dist_fromxmax (m)."

  else:
    theta = float(sys.argv[1]) #in deg
    phi = float(sys.argv[2]) #in deg
    alpha = float(sys.argv[3]) #in deg
    
    #decay = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
    energy = float(sys.argv[4]) #in EeV
    
    ## not needed
    dist_fromxmax = 100000 #float(sys.argv[6])# in m
    
    name=str(sys.argv[5])
    
    # In meters, but not taken from simulations! Has to be constant, taken from density profile + mean Xmax for proton with zenith 83deg (GRAND) and E10^18.5eV, just to get a rough idea about size of footprint. ATTENTION: has to be adjusted for all E and theta
    
    print "****Shower direction (zen, az) = ("+str(theta)+','+str(phi) +") deg, Mountain slope = "+str(alpha)+" deg, energy= " ,energy, " EeV, dist_fromxmax="+ str(dist_fromxmax)+"m, sim file", str(name)

    computeAntennaPos(theta, phi, alpha, energy, dist_fromxmax, name)
