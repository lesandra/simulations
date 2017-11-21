"""
    
                    Script danton_to_zhaires
                        Version 1.0
    Written by N. Renault-Tinacci from a script provided by C. Medina
                Using danton.py developped by V. Niess

"""

import json
import os, glob
import random
import shutil
import struct
import subprocess
import sys
import time
import numpy as np
import pylab as pl
import danton
import StringIO
import random
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

##################################################################
# Define a few systematic others
part_dic={'221.0':'eta','211.0': 'pi+', '-211.0': 'pi-','111.0': 'pi0', '22.0':'gamma', '13.0':'muon', '11.0': 'electron', '15.0':'tau', '16.0':'nu(t)', '321.0': 'K+', '-321.0': 'K-','130.0':'K0L', '310.0':'K0S','-323.0':'K*+'}
GEOMAGNET = (56.5, 63.18, 2.72) # Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg]).
ZREF=1500. #Altitude of the array
ground_alt = 0. #Ground Altitude. If changed, all showers with injection height (provided by DANTON) below this altitude won't run with ZHAires
DISPLAY = False #To plot the 3D map of the radio array (in particular the selected antennas)

##########################################################################################################
def main():
    """ Main script allowing to produce the ZHAires input files for each one of the showers from a DANTON library """

    ##################################################################
    # Test arguments
    if (len(sys.argv)<7 or len(sys.argv)>8):
        print """\
    This script will allow to produce the ZHAires input files from the DANTON output libraries (in txt format). 
    It is dedicated to earth-skimming neutrino simulation preparation.
    It creates a regular rectangular array more elongated along the shower axis.
    The seed of each shower is uniformly randomized between 0 and 1.

    Inputs : 
        work_dir = directory where all DANTON libraries, ZHAires input files and the simulation results will be stored.
        danton_lib = path to danton library
        distance = distance between the decay point and the beginning of the radio array [in m]
        slope = angle between horizontal and the slope along which the radio array is distributed [in degrees]
        height = maximum height that antennas distributed along the slope can reach [in m]
        step = separation between antennas [in m]
        azimuth = you can:   (az=0deg <=> northward, az=90deg <=> westward)
                    _ leave at the default value = 0 deg
                    _ set it to a random value (randomly drawn between 0 and 360 deg) (rotation of the array correspondingly to be implemented)
                    _ set ot to the wanted azimuth defined in GRAND coordinates [in degrees] (rotation of the array correspondingly to be implemented)
        
    Ouput:
        The script will produce as many ZHAires input files as there are showers in the DANTON output library. 
        They will located in the work_dir+"/inp/" directory.

    Usage:  python danton_to_zhaires.py work_dir danton_lib distance slope height step azimuth(=option)

    Notes:
        The global variable DISPLAY allows you to turn on/off the display of the 3D map of the radio array and the print-out of the processed showers
        The "compute_antenna_pos" function can easily be modified/replaced to build a different antenna array configuration (star-shape pattern, circularly distributed, ...)

    """
        sys.exit(1)

    ##################################################################
    # Retrieve the input parameters
    work_dir=str(sys.argv[1])
    if work_dir=='.':
        work_dir = os.getcwd()
    danton_lib=str(sys.argv[2]) #Primary neutrino energy
    Dd=float(sys.argv[3]) #distance from decay point to beginning of radio array
    slope=float(sys.argv[4]) #slope 
    hz=int(sys.argv[5]) #Array maximum height
    sep=float(sys.argv[6]) #separation between antennas
    try:
        AZIMUTH = str(sys.argv[7]) #azimuth
    except:
        AZIMUTH = str(0.)
    Ny = int(np.round(25e3/sep)) #number of lines in Y direction

    ##################################################################
    #Output directory for ZhAires results.
    DATAFS = work_dir+'showerdata/'
    showerdata_file = DATAFS+os.path.splitext(os.path.basename(danton_lib+'-showerdata.txt'))[0]
    inp_dir = work_dir+'/inp_D'+str(int(Dd))+'m_Z'+str(int(slope))+'deg_h'+str(int(hz))+'m_'+str(int(sep))+'m/'
    if not(os.path.isdir(DATAFS)):
        os.mkdir(DATAFS)
    if not(os.path.isdir(inp_dir)):
        os.mkdir(inp_dir)

    ##################################################################
    #### Parse the events and build up some statistics.
    data=[]
    for event in danton.iter_event(danton_lib):
            lastid = event.id
            for decay in event.decay:
                ### Fill data arrays
                [dataprod,depth,height,theta,azim,delta,et] = parse_build(event,decay,AZIMUTH)
                data.append((event.id, decay.tau_f.energy, event.primary.energy, depth, height, theta, azim, delta, decay.generation,
                    et, et/decay.tau_f.energy))  

    ### Save data arrays
                neu_file = DATAFS + str(event.id)+'.part'    
                np.savetxt(neu_file,dataprod,fmt='%d %d %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f' ,header="Neu-ID Prod-ID ux uy uz EProd ThetaProd Zenith-ZHAires Azimuth_ZHAires Height-above-sea-level Depth NuEner TauEner")  
                
    data = np.array(data)
    np.savetxt(showerdata_file,data,fmt='%d %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f %1.1f', header="Id TauEnergy NeuEnergy  Depth Height-above-sea-level Zenith-ZHAires Azimuth_ZHAires Delta Generation  ProdEnergy  ")

    print "+ {:} tau decays for {:} incoming neutrinos".format(len(data), lastid+1)
    print('******************************')
    print os.getcwd()

    ##################################################################
    ### Initialize a too large radio antenna array
    ANTENNAS=[]
    ANTENNAS=compute_antenna_pos(Dd,slope,sep,Ny,hz)

    ### Compute parameters for ZHAires input files
    showers=glob.glob(DATAFS+'*.part')
    if DISPLAY:
        shower_list=showers[15:21]
    else:
        shower_list = showers
    for fname in shower_list:
        [showerID,etot,azim,theta,multip,alt] = compute_shower_parameters(fname)
        if DISPLAY:
            print 'showerID = ',showerID,' Eshower = ',etot,' azimuth_zhaires = ',azim,' zenith_zhaires = ',theta,' injection height =',alt

        fileZha = inp_dir+showerID+ '.inp'
        dir=os.path.dirname(fileZha)
        if not os.path.exists(dir):
            os.makedirs(dir)

        ### Randomize the position of the core of the array in X,Y and Z
        CORE = random_array_pos(slope,sep)
        
        ### Reduce the radio array to the shower geometrical footprint (we account for a footprint twice larger than the Cherenkov angle)
        ANTENNAS2 = reduce_antenna_array(alt,theta,ANTENNAS,CORE,DISPLAY)
        if azim!=180.:
            #ANTENNAS3 = rotate_antenna_array(ANTENNAS2,azim) #to be implemented
            ANTENNAS3 = np.copy(ANTENNAS2)
        else:
            ANTENNAS3 = np.copy(ANTENNAS2)

        ### Write the ZHAires input file
        inpfile = open(fileZha,"w+")
        totito  = generate_input(showerID, etot, azim, theta, multip, alt,ANTENNAS3)
        inpfile.write(totito)
        inpfile.close()





##########################################################################################################
##########################################################################################################
###                                  Let's define useful functions                                     ###
##########################################################################################################
##########################################################################################################
def array_display(ANTENNAS=None,datamap=None,title=None):
    fig1 = pl.figure(1,figsize=(5*3.13,3.8*3.13))
    binmap = ListedColormap(['white', 'black'], 'indexed')
    dar=(np.max(ANTENNAS[:,0])-np.min(ANTENNAS[:,0]))/(np.max(ANTENNAS[:,1])-np.min(ANTENNAS[:,1]))
    if dar==0:
        dar=1
    xlbl='X [m]'
    ylbl='Y [m]'
    zlbl='Z [m]'
    ax = pl.gca(projection='3d')
    ax.scatter(ANTENNAS[:,0],ANTENNAS[:,1],ANTENNAS[:,2],c=datamap)
    ax.set_title(title)
    ax.view_init(25,-130)
    pl.show()
    return

##########################################################################################################
def GRANDtoZHAires(zen_GRAND=None, azim_GRAND=0):
    """ Convert coordinates from GRAND convention to ZHAires convention """

    zen = (180.-zen_GRAND)
    azim = azim_GRAND - 180.
    if azim>360:
      azim = azim-360.
    elif azim<0.:
        azim = azim+360.
    return [zen,azim]

##########################################################################################################
def compute_antenna_pos(distance=None, inclin=0., step=1e3, nsidey=None,hz=None):
    """ Generate antenna positions in a flat or inclined plane @ a given distance from decay"""
    """ Return N positions (x,y,z) in Zhaires coordinates """

    if inclin!=0. and inclin!=90.:
        disty = step*nsidey
        distz = hz/(np.sin(np.radians(inclin)))
        nsidez = int(distz/step)
        nsidex = nsidez  
        distx = distz*np.cos(np.radians(inclin))
        xi,yi = (distance,-0.5*disty) 
        xf,yf = (distance+distx,0.5*disty)
        xx, yy= np.meshgrid(np.arange(xi,xf,step*np.cos(np.radians(inclin))),np.arange(yi,yf,step))
        zz=(xx-distance)*np.tan(np.radians(inclin))
        xxr = np.ravel(xx)
        yyr = np.ravel(yy+0.5*step)
        zzr = np.ravel(zz)
        xyz = np.array([xxr, yyr, zzr]).T
    elif inclin==0.:
        nsidex = int(np.round(250e3/step))
        xi,yi = (distance,-0.5*nsidey*step)
        xf,yf = (distance+(nsidex*step),0.5*nsidey*step)
        xx, yy= np.meshgrid(np.arange(xi,xf,step), np.arange(yi,yf,step))
        zz=xx*0.
        xxr = np.ravel(xx)
        yyr = np.ravel(yy+0.5*step)
        zzr = np.ravel(zz)
        xyz = np.array([xxr, yyr, zzr]).T

    elif inclin==90.:
        disty = step*nsidey
        distz = hz
        zi,yi = (0.,-0.5*disty)
        zf,yf = (hz,0.5*disty)
        zz, yy= np.meshgrid(np.arange(zi,zf,step),np.arange(yi,yf,step))
        xx= (yy*0.+distance)
        xxr = np.ravel(xx)
        yyr = np.ravel(yy+0.5*step)
        zzr = np.ravel(zz)
        xyz = np.array([xxr, yyr, zzr]).T

    return xyz

##########################################################################################################
def random_array_pos(slope=0.,sep=1e3):
    """ Compute a random offset for 2 or 3 of the space dimensions depending of the slope """

    CORE = np.array([0.,0.,0.])
    CORE[1] = random.uniform(-1.,1.)*sep/2 # always need an offset on Y (perpendicular to trajectory)
    if slope!=90. and slope!=0.: # random position along slope => x and z are random and their offset is related
        CORE[0] = random.uniform(0.,1.)*sep*np.cos(np.radians(slope))
        CORE[2] = CORE[0]*np.sin(np.radians(slope))
    elif slope==0.: # z = 0 => no offset in Z required
        CORE[0] = random.uniform(0.,1.)*sep*np.cos(np.radians(slope))
    elif slope==90.: # x = distance => no offset in X required
        CORE[2] = random.uniform(0,1.)*sep

    return CORE

##########################################################################################################
def getCerenkovAngle(h=100e3):
    """ Compute the Cherenkov angle of the shower at the altitude of injection """

    # h in meters
    n = 1.+325.e-6*np.exp(-0.1218*h*1e-3)      # Refractive index Zhaires (see email M. Tueros 25/11/2016)
    alphac = np.arccos(1./n)  
    return alphac

##########################################################################################################
def reduce_antenna_array(injh=None,theta=None,ANTENNAS=None,core=[0.,0.,0.],DISPLAY=False): 
    """ Reduce the size of the initialized radio array to the shower geometrical footprint by computing the angle between shower and decay-point-to-antenna axes """
    """ theta = zenith in ZHAires convention [in deg], injh = injection height [in m] """
    
    # Convert back to GRAND coordinates
    zenr = np.pi - np.radians(theta)

    # Shift antenna array with the randomized core position
    ANTENNAS[:,0] = ANTENNAS[:,0]+core[0]
    ANTENNAS[:,1] = ANTENNAS[:,1]+core[1]
    ANTENNAS[:,2] = ANTENNAS[:,2]+core[2]
    
    # Compute angle between shower and decay-point-to-antenna axes
    u_ant = ANTENNAS+[0.,0.,-injh]
    u_ant = (u_ant.T/np.linalg.norm(u_ant,axis=1)).T
    u_sh = [np.sin(zenr),0.,np.cos(zenr)]
    ant_angle = np.arccos(np.matmul(u_ant, u_sh))

    # Remove antennas of the initial array that are located outside the "footprint"
    omegar = getCerenkovAngle(injh)*2. #[in rad] # Accounting for a footprint twice larger than the Cherenkov angle
    angle_test = ant_angle<=omegar
    sel = np.where(angle_test)[0]
    ANTENNAS2 = ANTENNAS[sel,:]

    # Remove the farthest antennas to reduce the number of antenna positions to simulate so that this number falls below 1000
    while np.shape(sel)[0]>999:
        x_ant_max = np.max(ANTENNAS2[:,0])
        antisel = np.where(ANTENNAS2[:,0]==x_ant_max)[0]
        ANTENNAS2= np.delete(ANTENNAS2,antisel,0)
        sel= np.delete(sel,antisel,0)

    # 3D Display of the radio array
    if DISPLAY:
        ant_map_i = np.zeros(np.shape(ANTENNAS)[0])
        ant_map_i[sel]=1.
        cc = np.zeros((np.size(ant_map_i),3))
        cc[np.where(ant_map_i==0),:]=[1,1,1]
        array_display(ANTENNAS,ant_angle,'Shower axis to decay point-antenna axis angle map')
        array_display(ANTENNAS,cc,'Selected antenna map')

    return ANTENNAS2

##########################################################################################################
def rotate_antenna_array(ANTENNAS=None,azim=0.):
    """ For a azimuth different of 0, one need to rotate the radio array with azimuth """
    """ so that the longest side of the array is aligned with shower axis """

    #To implement

    return ANTENNAS2

##########################################################################################################
def parse_build(event=None,decay=None,AZIMUTH=0.):
    """ From a DANTON event, retrieve the shower characteristics and build some statistics """
    """ Produce a table for each shower with all its parameters """

    if AZIMUTH=='random':
        azim_i = random.uniform(0.,1.)*360.
    else:
        azim_i = float(AZIMUTH)

    r1 = decay.tau_i.position
    r2 = decay.tau_f.position
    R2 = np.linalg.norm(r2)
    u2 = decay.tau_f.direction
    c = np.dot(u2, event.primary.direction)
    if c > 1: delta = 0.
    else: delta = np.arccos(c)   
    depth = danton.EARTH_RADIUS - np.linalg.norm(r1)
    height = R2 - danton.EARTH_RADIUS 
    theta_danton = np.degrees(np.arccos(np.dot(u2, r2) / R2))
    theta, azim = GRANDtoZHAires(theta_danton,azim_i)
              
    dataprod = []
    et=0.0 #in GeV
    for i in decay.product:
        idp = i[0]
        pp=i[1]
        ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2) #in GeV ignoring mass energy
        et=et+ep #in GeV
        up=pp/ep
        thetap_danton=np.degrees(np.arccos(np.dot(up, r2) / R2))
        thetap, azim = GRANDtoZHAires(thetap_danton,azim_i)
        dataprod.append((event.id,idp,up[0],up[1],up[2],ep,thetap,theta,azim,height,depth,event.primary.energy,
            decay.tau_f.energy))
    dataprod = np.array(dataprod)

    return dataprod,depth,height,theta,azim,delta,et

##########################################################################################################
def compute_shower_parameters(fname=None):
    """ From the .part file produce for each shower, compute the parameters that will allow the production of the ZHAires input file """

    num_lines = sum(1 for line in open(fname))
    TASK0=fname.replace('/','.')
    showerID=TASK0.split('.')[-2]
    datap = np.loadtxt(fname)
    idpart =[]
    multip = []
    if (num_lines<3):
        part=datap[1]
        #print part
        etot = 1.0e+09*datap[5]
        theta = datap[6]
        azim = datap[8]
        alt = datap[9]
        #print alt
        prop = 1.0e+09*datap[5]/etot
        multip.append([part_dic[str(part)],prop])
    else:
        part = datap[:,1]
        etot = 1.0e+09*sum(datap[:,5])
        theta = datap[0,6]
        azim = datap[0,8]
        alt = datap[0,9]
        #print alt
        prop = 1.0e+09*datap[:,5]/etot
        for i in range(0,len(datap)):
            multip.append((part_dic[str(part[i])],prop[i]))

    return showerID,etot,azim,theta,multip,alt

##########################################################################################################
def generate_input(task=0,energy=None, azimuth=None, zenith=None, products=None, height=None, antennas=None):
    """Generate the input stream for ZHAIRES."""

    a=" ".join(map(str, products))
    b="".join( c for c in a if  c not in "(),[]''")

    seed = random.uniform(0.,1.)

    # Format the stream.
    stream = [
        "AddSpecialParticle      RASPASSProton    /home/renault/zhaires/RASPASSprimary/RASPASSprimary Proton",
        "AddSpecialParticle      RASPASSIron      /home/renault/zhaires/RASPASSprimary/RASPASSprimary Iron",
        "AddSpecialParticle      RASPASSelectron  /home/renault/zhaires/RASPASSprimary/RASPASSprimary Electron",
        "AddSpecialParticle      RASPASSpi+  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi+",
        "AddSpecialParticle      RASPASSpi-  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi-",
        "AddSpecialParticle      RASPASSpi0  /home/renault/zhaires/RASPASSprimary/RASPASSprimary pi0",
        "AddSpecialParticle      RASPASSMulti /home/renault/zhaires/RASPASSprimary/RASPASSprimary {:s}".format(b),
        "#########################",
        "TaskName {:s}".format(task),
        "PrimaryParticle RASPASSMulti",
        "PrimaryEnergy {:.5E} eV".format(energy),
        "PrimaryZenAngle {:.5f} deg".format(zenith),
        "PrimaryAzimAngle {:.5f} deg Magnetic".format(azimuth),
        "ForceModelName SIBYLL",
        "SetGlobal RASPASSHeight {:.5f} m".format(height),
        "RandomSeed {:.5f}".format(seed),
        "########################",
        "PropagatePrimary On",
        "SetGlobal RASPASSTimeShift 0.0",
        "SetGlobal RASPASSDistance 0.00"
    ]

    for a in antennas:    
        stream.append("AddAntenna {:1.2f} {:1.2f} {:1.2f}".format(a[0],a[1],a[2]))

    stream += [
        "##########################",
        "TotalShowers 1",
        "RunsPerProcess Infinite",
        "ShowersPerRun 1",
        "Atmosphere 1",
        "AddSite Ulastai 42.55 deg 86.68 deg {:.3f} m".format(ZREF),
        "Site Ulastai",
        "Date 1985 10 26",
        "GeomagneticField On",
        "GeomagneticField {:.4f} uT {:.2f} deg {:.2f} deg".format(*GEOMAGNET),
        "GroundAltitude {:.3f} m".format(ground_alt),
        "ObservingLevels 510 200 g/cm2   900 g/cm2",
        "PerShowerData Full",
        "SaveNotInFile lgtpcles All",
        "SaveNotInFile grdpcles All",
        "RLimsFile grdpcles 0.000001 m 10 km",
        "ResamplingRatio 100",
        "#########################",
        "RLimsTables 10 m 10 km",
        "ELimsTables 2 MeV 1 TeV",
        "ExportTables 5501 Opt a",
        "ExportTables 1293 Opt as",
        "ExportTable 1793 Opt s",
        "ExportTable 1793 Opt as",
        "ExportTable 7293 Opt s", 
        "ExportTable 7793 Opt s",
        "ExportTable 7993 Opt s",
        "########################",
        "ForceLowEDecay Never",
        "ForceLowEAnnihilation Never",
        "########################",
        "ZHAireS On",
        "FresnelTime On",
        "FresnelFreq Off",
        "TimeDomainBin 1 ns",
        "AntennaTimeMin -100 ns",
        "AntennaTimeMax 500 ns", #can be extended until 3e-6s but if then there is still nothing then there must be a problem somewhere
        "######################", 
        "ElectronCutEnergy 1 MeV",
        "ElectronRoughCut 1 MeV",
        "GammaCutEnergy 1 MeV",
        "GammaRoughCut 1 MeV",
        "ThinningEnergy 1.e-4 Relative", #It can be 1e-5, 1e-6 or below. But running time inversely proportional to it.
        "ThinningWFactor 0.06"
    ]

    return "\n".join(stream)

##########################################################################################################
if __name__ == '__main__':
    main()
