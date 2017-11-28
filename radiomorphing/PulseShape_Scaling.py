#python PulseShape_Scaling.py /media/sf_work/Paris/scripts_GRAND/Olivier_scripts/EnergyScan/ EnergyScan_2 /media/sf_work/Paris/scripts_GRAND/Olivier_scripts/EnergyScan/ EnergyScan_8 30 80

# NOTE: script for testing impact of scaling (and interpolation (along ray) for height an zenith scaling)


#'this script shall scale a whole traces like done for the peak amplitudes.
# Energy works
# azimuth works when including extra factor for flipping the traces - has to be tested, ATTENTION and included back again

# at the moment this script does the scaling in shower coordinates, meaning a scaling of Ev, EvxB, and EvxvxB which you obtain after filtering as columns 4,5,6
# if you wanna use Ex, Ey, Ez chosse columns 1,2,3
# ="= just filtered traces are used

# scaling with height and zenith: now included since position of the antennas would be slightly changed and then a pulse shape interpolation is needed -> has to be tested

# scaled traces will be then saved to a data files looking like the others ones, but just with entries i the v, vxB, and vxvxB columns (0, 4,5,6)

# 29Aug2017: corrections- kAz jus on EvxB, hxref=href +8000*tan(-90deg-zen)
# restructured script so that multiplying factors and strecting poition os now a function where one just have to hand over path and antenna number + all need primary informations+ all star shape positions
#30Aug 2017 included back transformation to Exyz, innclude hilbert etc,
#TODO:  make the outputfile containing the scaled traces aquivalent to the originals just in an folder caled scaled+...
# remove frequency dependency, read in simulated raw files for a later applying of the antenna response (another module of script chain)
# remove complete comparison

import os
import re
import sys
import numpy as np
from PulseShape.frame import UVWGetter, XYZGetter, get_rotation
from PulseShape.utils import getCerenkovAngle, load_trace

def getAirDensity(h):
  #% Using isothermal Model
  rho_0 = 1.225#; % kg/m3
  M = 0.028966#;  %kg/mol
  g = 9.81#; %ms-2
  T = 288.#; % K
  R = 8.32#;
  rho = rho_0*np.exp(-g*M*h/(R*T))
  return rho

###################################

def getXmax(primarytype, energy, zen2): # type of primary (electron or pion, energy in EeV, zen2 (GRAND) in rad
    if primarytype=='electron': # aprroximated by gamma shower
        a=82.5 # g/cm2
        c=342.5 #g/cm2
    if primarytype=='pion': # aprroximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2
    Xmax= a*np.log10(energy*10**6.)+c # E/EeV* 10**6. to be in TeV

    return Xmax#/abs(np.cos(np.pi-zen2)) # TODO: how to correct for slanted shower

def dist_decay_Xmax(zen2, injh2, Xmax_primary): #zen2: zenith of target shower
      #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    #zen2=89.5*np.pi/180

    hD=injh2
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    while X< Xmax_primary:
        i=i+1
        ai=i*100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        X=X+ rho_0*np.exp(-g*M*hi/(R*T)) * (deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values

    #print X, h, ai
    return h, ai # Xmax_height in m, Xmax_distance in m

def scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, thetageo):
    # 2: target shower, 1: generic shower
    #################### Energy scaling
    #% Energy
    #kE = Esh*1e18./Esh_ref
    kE = E2/E1 # both in 1e18eV

    ############## Azimuth scaling
    #% Azimuth
    Bgeo = [np.sin(thetageo), 0, np.cos(thetageo)]#  % With North = magnetic North
    vref = [np.cos(az1)*np.sin(zen1)*-1., np.sin(az1)*np.sin(zen1)*-1., np.cos(zen1)*-1.]# *-1: account for angle conventions
    vxB_ref = np.cross(vref,Bgeo)
    vxB_ref = np.linalg.norm(vxB_ref)/(np.linalg.norm(vref)*np.linalg.norm(Bgeo))
    v = [np.cos(az2)*np.sin(zen2)*-1., np.sin(az2)*np.sin(zen2)*-1., np.cos(zen2)*-1.]# *-1: account for angle conventions
    vxB = np.cross(v,Bgeo)
    vxB = np.linalg.norm(vxB)/(np.linalg.norm(v)*np.linalg.norm(Bgeo))
    kAz = vxB/vxB_ref   # =cos(az2)/cos(az1)

    #print 'included extra scaling factor cos(Delta azimuth) to respect flipping of components'
    kAz = kAz#* np.cos(az2-az1)

    h_ref=injh1
    h=injh2
    #%############## Height+Zenith, distance injection point to xmax rougly 8000m
    #hx_ref = h_ref+8000*np.cos(zen1) #   % Height at reference shower Xmax
    hx_ref = h_ref+8000*np.tan(0.5*np.pi-zen1) #   % Height at reference shower Xmax
    ac_ref = getCerenkovAngle(hx_ref)
    rho_ref = getAirDensity(hx_ref)
    #hx = h+8000*np.cos(zen2)#   % Height at  shower Xmax
    hx = h+8000*np.tan(0.5*np.pi-zen2)#   % Height at  shower Xmax
    ac = getCerenkovAngle(hx)
    rho = getAirDensity(hx)
    kStretch = float(ac)/float(ac_ref)#  % Stretch factor for the antenna grid
    kRho = np.sqrt(rho_ref/rho)
    kHeight = kRho/kStretch
    kAmp=kE*kAz*kHeight

    #print 'kStretch ', kStretch, ' kAmp ', kAmp,  ' kE ', kE, ' KAz ', kAz, ' kHeight ', kHeight

    return kStretch, kE, kAz, kHeight




def scalingpulse(dist1, E1, az1, zen1, injh1, E2, az2, zen2, injh2, primary, phigeo, thetageo, l,  positions, path): # hand over parameters from reference shower 1 and target shower 2, the number of the antenna in the star shape you would like to have  and all position of complete starshape (for strechting needed), the path to the folder containing the sim, and for now the frequencies (should be removed if one included the antenna response

#SCALING
    kStretch, kE, kAz, kHeight= scalingfactors(E1, az1, zen1, injh1, E2, az2, zen2, injh2, thetageo)
    kAmp=kE*kAz*kHeight
    if l==0:
        print 'kStretch ', kStretch, ' kAmp ', kAmp,  ' kE ', kE, ' KAz ', kAz, ' kHeight ', kHeight

    #print 'Amplitude changed by kAmp= ' + str(kAmp), ' Position changed by kStretch= '+ str(kStretch)



 ###############################################################################
 #### scaling electric fields amplitude
 ################################################

    try:
    ##read in full traces of antenna l: 0:time in ns, 1,2,3: vector field , 4,5,6: efield
        txt1 = load_trace(path, l)[:-2,:]
    except IOError:
        print("antenna ID ",str(int(l)), " file doesn't exist")
        sys.exit()

    # TODO: just 4 columns. use a temporary array for EvxB => Eshower A = np.array

    # Convert efield to shower coordinates to apply the scaling
    R = get_rotation(zen1, az1, phigeo, thetageo)
    EshowerA = np.dot(txt1[:,1:], R)

    ### Sciling, kHeight includes 1/kStretch
    EshowerA.T[0] *= kE * kHeight
    EshowerA.T[1] *= kE * kAz * kHeight
    EshowerA.T[2] *= kE * kHeight

    #Backtrafo of efield from shower coord (1,2,3) in xyz (4,5,6) after scaling and/or stretching using the target angles
    Rt = get_rotation(zen2, az2, phigeo, thetageo).T
    v2 = Rt[0,:]
    txt1[:,1:] = np.dot(EshowerA, Rt)


###################################################
#### stretching of positions
###############################

### Note: if kStrecht !=1 then one has to perform an interpolation of antenna positions to be able to compare the signals target to scaled

    # default parametes of star shape simulation
    #rings = 15
    angles= 8
    rings = len(positions[:,1])/angles
    beta= (360./angles)/180.*np.pi


################################
 # Calculating the new stretched antenna positions in the star shape
    offinz= np.mean(positions[:,2])
    offiny= np.mean(positions[:,1])
    offinx= np.mean(positions[:,0])
    pos= np.zeros([len(positions[:,1]),3])

    # rotate into shower coordinates for preparation to get the strechted antenna position to compare to - maybe there is a nicer way...
    GetUVW = UVWGetter(offinx,offiny,offinz, zen1, az1, phigeo, thetageo)
    for i in np.arange(0,len(positions[:,1])):
        pos[i,:] = GetUVW(positions[i,:])


    ###################################   substitue pos by pos and add if condition here
    if kStretch!=1.:
        r= np.linalg.norm(pos[6]- pos[5]) # first rays are formed
        step=1
        #print r
        ## get the distance between the antenna
        if np.linalg.norm(pos[6]-pos[5]) != 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
            #print ' not equal', np.linalg.norm(pos[6]-pos[5]), 0.5*np.linalg.norm(pos[7]-pos[5])
            #print "ATTENTION: check whether first rings formed or rays if streching is needed"
            r= np.linalg.norm(pos[16]- pos[8])
            step=8
            #print r

        r=r*kStretch # just STRETCHING the radius of the rings

        # default parametes of star shape simulation
        rings = 15
        angles= 8
        beta= (360./angles)/180.*np.pi

        ## get the stretched positions in teh star star shape pattern
        #pos = np.zeros([rings*angles,3])

        # first rays are formed
        if np.linalg.norm(pos[6]-pos[5]) == 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
            for n in range(0, angles): # second rings
                for m in range(0, rings): #first rays
                    #if n==0 and m==0:
                        #print " first rays than rings"

                    pos[n*rings+m,1 ]= (m+1)*r* np.cos(n*beta)
                    pos[n*rings+m,2 ]= (m+1)*r* np.sin(n*beta)
                    pos[n*rings+m,0 ]= 0.
                    #if l*angles+ n <8:
                    #print l, n*rings+l,  pos[n*rings+l], pos[n*rings+l]

        if np.linalg.norm(pos[6]-pos[5]) != 0.5*np.linalg.norm(pos[7]-pos[5]): # first rings are formed
                for m in range(0, rings): #sencond rays
                    for n in range(0, angles): # first rings
                        #if n==0 and m==0:
                            #print " first rings than rays"

                        pos[m*angles+n,1 ]= (m+1)*r* np.cos(n*beta) # vxB
                        pos[m*angles+n,2 ]= (m+1)*r* np.sin(n*beta) # vxvxB
                        pos[m*angles+n,0 ]= 0. # along v


################# CALCULATION OF NEW POSITION VECTOR

    ### the new/target position vector of the star shape plane
    ### pos_prime= p_prime + v_prime *(d_prime + toX)

    #print('distance to max to plane: ',dist1) ## ==toX: fixed distance in m from Xmax to the plane along the shower axis
    #print(v2) # == v_prime: vector of where target shower goes to
    Xmax_primary= getXmax(primary, E2, zen2)# approximation based on values from plots for gamma (=e) and protons (=pi) # g/cm2
    #print("xmax value " , Xmax_primary)
    Xmax_height, Xmax_distance = dist_decay_Xmax(zen2, injh2, Xmax_primary)# 8000.# d_prime: distance from decay point to Xmax

    #print("distance to decay to Xmax in m ", Xmax_distance)
    decay=np.array([0.,0.,injh2]) # ==: decay position as defined in zhaires sim, from DANTOn files

    # new position vector:
    x2= decay - v2 * (Xmax_distance+ dist1) # to account for going from Zhaires to GRAND conv


 ##############################   Backtrafo to XYZ
    ### Now the new 'stretched' positions are calculated in the xyz components, backrotation
    stretch2 = np.zeros([len(pos[:,1]),3])

    ### TODO: generell backtrafo of positions
    GetXYZ = XYZGetter(x2[0],x2[1],x2[2],zen2, az2,phigeo,thetageo)
    for m in range(0, len(pos[:,1])):
        stretch2[m,:] = GetXYZ(pos[m])
    #stretch2[l]  # the new wanted antenna position after stretching

    if l==0:
        print 'position ref to Xmax: ', x2,' position decay: ', decay, ' shower direction: ', v2, ' distance to Xmax: ', dist1, ' distance between Xmax and plane: ', Xmax_distance
        print len(x2), l, len(stretch2),


    #print ' scaling done , positions ', stretch2[l]
    return txt1, stretch2[l]



####################################





#taken from oliviers scripts: orientation/direction of magnetic field
phigeo =0*np.pi/180.  # 182.66#; (ie pointing 2.66 degrees East from full North) # phigeo= 0 from simulations inputfile % In both EVA & Zhaires, North = magnetic North
thetageo =(180.-27.05)*np.pi/180. # 152.95*np.pi/180. #27.05*np.pi/180. #; (pointing down)-62.95


### NOTE not needed anymore
start=0 #int(sys.argv[7]) # number of antenna in star shape which shall be read in
#end=start+1#int(sys.argv[8])

if start==0.:
    print '\n'
    print '########################################  new plane gets scales ... ' + sys.argv[2]

# Simulated shower# the one which gets scaled

run_sim = sys.argv[2]
path =sys.argv[1]  + str(run_sim) #path to first simulation which shall later be rescaled or so
steerfile_sim = path + "/../MasterIndex"




## reference shower, parameters are read-in from input file
Run='{0} .*'.format(str(run_sim))
E1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[2]) # energy in EeV
injh1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[5]) # injection height
dist1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[6]) # distance xmax to plane

zen1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[3])*np.pi/180. # from deg in rad
az1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[4])*np.pi/180. # from deg in rad

## conversion from Aires to GRAND convention
zen1=np.pi-zen1
az1=np.pi+az1

if start==0:
    print path

    print steerfile_sim

    print "Reference shower " + str(run_sim)
    print "Energy in EeV: " + str(E1) + ', antenna distance from Xmax in m: '+ str(dist1) + ", zenith: " +str("%.5f"%(zen1*180./np.pi)) +  ", azimuth: " + str("%.5f" %(az1*180./np.pi)) + ", height/m: " +str(injh1) #+ ", with dist. of Xmax to injection point in m: " +l1




##### target shower parameters, parameters handed over in GRAND conventions


E2= float(sys.argv[3])# energy in EeV
injh2= float(sys.argv[4]) # injection height in m above obs level
#dist2= # distance xmax to plane in m, horizontal distance --> has to be the same as the reference shower

zen2=float(sys.argv[5])*np.pi/180. # from deg in rad, Grand con
az2=float(sys.argv[6])*np.pi/180. # from deg in rad, Grand con
primary=str(sys.argv[7])

## conversion from Aires to GRAND convention, not needed any more
## hand over angles in GRAND not in aires
#zen2=np.pi-zen2
#az2=np.pi+az2

if start==0:

    print "Target shower parameters"
    print "Energy in EeV: " + str(E2) + ', antenna distance from Xmax in m: '+ str(dist1) + ", zenith: " +str("%.5f"%(zen2*180./np.pi)) +  ", azimuth: " + str("%.5f" %(az2*180./np.pi)) + ", height/m: " +str(injh2) #+ ", with dist. of Xmax to injection point in m: " +str(l2)




######################################### Stretching the position
################### Assume kStretch=1, no change of zen and height
##r=r*kStretch # just stretching the radius of the rings

posfile = path  +'/antpos.dat'
positions=np.genfromtxt(posfile)

pos_new= np.full_like(positions, 0.) #np.array([len(positions), 3])

 ################################################################################


end=len(positions[:,0])
#### loop over all antenna positions, outer positions should be skipped since then the interpolation doesnt work anymore
for l in np.arange(start,end):#0,120-1):
    if l==0:
        print path+'/a'+str(int(l)) + '.trace'

    txt3, pos_new[l,:] = scalingpulse(dist1, E1, az1, zen1, injh1, E2, az2, zen2, injh2, primary, phigeo, thetageo, l,  positions, path) # always hand over all need parameters,1 3D pulse, and all antenne positions


#########################

# comment in if hilbert envoleope needed one day
    #hexf2 = abs(hilbert(txt3.T[1])) # Ex
    #heyf2 = abs(hilbert(txt3.T[2]))# Ey
    #hezf2 = abs(hilbert(txt3.T[3])) # Ez
    #exm2 = max(hexf2)
    #eym2 = max(heyf2)
    #ezm2 = max(hezf2)

    #amp2 = sqrt(exm2*exm2+ eym2*eym2+ezm2*ezm2)





    ######Writing to file for later use
    directory= path+'/../scaled_'+str(sys.argv[2])+'/'
    if not os.path.exists(directory):
        os.makedirs(directory)

    #TODO: give it an equivalent name as original one jut in separate folder
    #name2= 'a'+str(sys.argv[8])+'_'+str(l) #sys.argv[1]+"a"+sys.argv[2]+'_'+str(int(x1*1e-6)) + '-' + str(int(x2*1e-6)) + 'MHz'
    #name3= name2+'.dat'
    name3= directory+'a'+str(l)+'.trace'
    FILE = open(name3, "w" )

    # 1,2,3 Ev,vxb,vxvxb then, 4,5,6 Exyz: since in original raw traces Exyz are also in columns 456
    for i in range( 0, len(txt3.T[0]) ):
        print >>FILE,"%3.2f	%1.5e	%1.5e	%1.5e" % (txt3.T[0,i], txt3.T[1,i], txt3.T[2,i], txt3.T[3,i] )


    ## in the last line of the file we wanna write the max. ampl of the hilbert envelope
    ## something like: amp exm eym ezm
    #print >>FILE, ""
    #print >>FILE,"%1.5f	%1.5e	%1.5e	%1.5e" % (amp2,exm2, eym2, ezm2)

    FILE.close()
    if l==0:
        print 'scaled traces saved like this: /../scaled_'+str(sys.argv[2])+'/'+'a'+str(l)+'.trace'






# save as well the posiion file somewhere if you scale the complete star shape pattern
posfile_new = directory  +'/antpos.dat'

#numpy.savetxt(posfile_new, (positions), fmt='%1.3e')

file_ant=open(posfile_new, 'w')
for i in range( 0, len(positions) ):
        print >>file_ant,"%.3f	%.3f	%.3f" % (pos_new[i,0], pos_new[i,1],pos_new[i,2] )
file_ant.close()

print len(positions), ' antennas scaled'
print 'antenna positions save in:', posfile_new
