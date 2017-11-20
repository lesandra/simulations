
#python PulseShape.py /media/sf_work/Paris/scripts_GRAND/Olivier_scripts/EnergyScan/antpos_desired1.dat /media/sf_work/Paris/scripts_GRAND/test/Simulations.dat /media/sf_work/Paris/scripts_GRAND/Olivier_scripts/EvetADetailed2/


# this script does a interpolation of a complete Pulse at any antenna position you desire.
#therfore you have to hand over a list of antenna position you would like to have, a file containing the simulations which should be use (names of the planes) and a path whre to find these simlations 
# the script calculates all the prjections on the exiting planes which are needed, and hand the traces and positions over to PulseShape_Interpolation.py which performs the interpoaltion alwys in between two positions
# whether you wanna use filtered traces is set in this script by hand at the beginning
# it returns files (t, Ex,Ey,Ez + peak amps) in a folder InterpoaltedSignals if it exists. It would make sense also to save the list of the positions in that folder too.


# at the moment the scripts get time, Ex, Ey, Ez components (filtered) from the files, not yet respecting projection to shower coord.:
#.T[0] =time
#.T[1] =Ex
#.T[2] =Ey
#.T[3] =Ez
#.T[4] =Ev
#.T[5] =EvxB
#.T[6] =EvxvxB


import sys
from sys import argv

import numpy as np
from numpy  import *
from numpy import linalg

#from scipy import signal
#from scipy.signal import hilbert #comment in if needed
import matplotlib.pyplot as plt
#import plotdata 
import pylab
import os
#from matplotlib.pyplot import cm 

import re

#from mpl_toolkits.mplot3d import Axes3D  ## just if 3D plotting is needed


import PulseShape_Interpolation as pulse

def ProjectPointOnLine(a, b, p):
    ap = p-a
    ab = b-a
    point = a + dot(ap,ab)/dot(ab,ab) * ab
    return point

def ProjectPointOnPlane(a,b,d, p):
    n= np.cross(a, b)
    n= n/np.linalg.norm(n)
    t= (np.dot(d,n)-np.dot(p,n))/np.dot(n,n)
    point= p+ t*n
    return point

def getCerenkovAngle(h):
   #% h in meters
   n = 1.+325.e-6*np.exp(-0.1218*h*1e-3)#;  % Refractive index Zhaires (see email M. Tueros 25/11/2016)
   alphac = np.arccos(1./n);
   return alphac


def GetUVW(pos, cx, cy, cz, zen, az, phigeo, bfieldangle):
   relpos = pos-np.array([cx,cy,cz])
   inc=bfieldangle#/180.*np.pi #magnetic field direction on SKA site
   #inc=-1.0554456843876574 #SKA
   #B = np.array([0,np.cos(inc),-np.sin(inc)]) 
   B = np.array([np.cos(phigeo)*np.sin(inc), np.sin(phigeo)*np.sin(inc),np.cos(inc)]) #from oliviers script including phigeo
   v = np.array([np.cos(az)*np.sin(zen)*-1.,np.sin(az)*np.sin(zen)*-1.,np.cos(zen)*-1.]) # or *-1: change the direction
   #print v
   vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]]) # crossproduct
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])# crossproduct
   return np.array([np.inner(v,relpos), np.inner(vxB,relpos),np.inner(vxvxB,relpos)]).T # vector dot

def GetXYZ(pos, zen, az, bfieldangle):
   inc=bfieldangle #magnetic field direction on SKA site
   #inc=-1.0554456843876574 #SKA
   B = np.array([0,np.cos(inc),-np.sin(inc)])
   v = np.array([np.cos(az)*np.sin(zen)*-1.,np.sin(az)*np.sin(zen)*-1.,np.cos(zen)*-1.]) # or *-1
   #print v
   vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
   vxB = vxB/np.linalg.norm(vxB)
   vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
   return pos[0]*vxB+pos[1]*vxvxB+pos[2]*v

DISPLAY=0
SCALED=1 # scaled traces shall be read in

#####
# add request whether f1 and f2 as freqeuncies are handed over
# if not: full=1, oherwise full =0
full=1 # fullband =1 == no filtering, raw pulses
f1= 60.e6 #MHz
f2=200.e6 # MHz
######



    
path0= sys.argv[1]# path to file with desired antenna positions
 #sys.argv[1] # path to file containing all the simulations you wanna use
path1= sys.argv[3] # path to simulations
path2= sys.argv[4] #path to folder for final traces




#taken from oliviers scripts: orientation/direction of magnetic field
phigeo =0*np.pi/180.  # 182.66#; (ie pointing 2.66 degrees East from full North) # phigeo= 0 from simulations inputfile % In both EVA & Zhaires, North = magnetic North
thetageo =(180.-27.05)*np.pi/180.  # 152.95*np.pi/180. #27.05*np.pi/180. #; (pointing down)-62.95



# Hand over a list file including the antenna positions you would like to have. 
posfile = path0 
positions=np.genfromtxt(posfile)# positions[:,0]:height, positions[:,1]:x,positions[:,2]:y  
if DISPLAY==1:
    print 'desired positions: '
    print positions, len(positions)
if len(positions) <=1:
    print "Files of desired positions has to consist of at least two positions, Bug to be fixed"

simfile = str(sys.argv[2])
#print simfile
sims = open(simfile).read().split('\n') 
#print len(sims), sims
print sims

try:
    ## should be the same for all the planes
    steerfile_sim = path1+"/MasterIndex" # '/{0}/inp/{1}.inp'.format(str(sims[0]),str(sims[0]))
    Run='{0} .*'.format(str(sims[0]))

    #E1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[2]) # energy in EeV
    #injh1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[5]) # injection height
    dist1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[6]) # distance xmax to plane
except: 
    
    ## should be the same for all the planes
    steerfile_sim = path1+"/../MasterIndex" # '/{0}/inp/{1}.inp'.format(str(sims[0]),str(sims[0]))
    Run='{0} .*'.format(str(sims[0]))

    #E1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[2]) # energy in EeV
    #injh1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[5]) # injection height
    dist1=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[6]) # distance xmax to plane

#print 'length arguments' , len(sys.argv)
if len(sys.argv) <= 5.: # in case of not using it for scaled traces
    print 'angles read out of masterindex file'
    zen=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[3])*np.pi/180. # from deg in rad
    az=float(np.genfromtxt(re.findall(Run,open(steerfile_sim,'r').read()))[4])*np.pi/180. # from deg in rad
    # since input file containes Aires convention --- convert to GRAND conventions
    zen=np.pi-zen
    az=np.pi+az

else: # angles in GRAND coordinates handed over,
    #print 'angles handed over in GRAND conv.'
    zen= float(sys.argv[5])*np.pi/180
    az=float(sys.argv[6])*np.pi/180
    
### scaled traces shall be read in
if SCALED==1:
    path1=path1+'/scaled_'


#zen=np.pi-zen
#az=np.pi+az
#print 'primary zen, az: ' + str(zen*180./np.pi) + ', ' + str(az*180./np.pi) , 'in GRAND coordninates'





########################## GET THE NEIGHBOURS
# Here the closets Neighbours should be found....

#### Finding the Neighbours:  In principal one would like to check which star shape pattern are the closest etc.
# it reads in all star shape pattern positions from the simulations you hand over via simulations.dat
positions_sims=np.zeros([len(sims),120,3])
for i in np.arange(0,len(sims)): # loop over simulated antenna positions
    print sims[i]
    #if i==0: ## Fixed 14.11.2017
        #print "WARNING: non-scaled pulses loaded... change if needed, eg by handing over another path where scaled are saved"
    #posfile = path1 + '/scaled_'+str(sims[i]) +'/antpos.dat' 
    posfile = path1 +str(sims[i]) +'/antpos.dat'
    print posfile
    if DISPLAY==1:
        print posfile
    positions_sims[i,:,:]=np.genfromtxt(posfile)
    
    
if DISPLAY==1:
    print "Antenna files loaded"    
#print positions_sims[0, 10]
    

#ATTENTION herethere should go a LOOP over b over desired antenna position in the list and find always the closest neighbours and get the interpolated pulse shape

for b in range(len(positions)):
    print "##############  begin of interpolation at desired position ", b, ' at ',  positions[b]


    # desired positions has to get projected on all the planes to get the orthogonal distances to the planes 
    # then one can check in between which planes the points is by choosing the two smallest distances

    # first find the two planes which are the closest
    dist_plane=np.zeros([len(sims)])  # 1dim distance array
    dist_value=np.zeros([len(sims)])  # 1dim distance array
    for i in np.arange(0,len(sims)): 

        PointPrime=ProjectPointOnPlane(positions_sims[i,10]-positions_sims[i,11], positions_sims[i,40]-positions_sims[i,41], positions_sims[i,1], positions[b])
        #print positions_sims[i,10]-positions_sims[i,11], positions_sims[i,12]-positions_sims[i,20]
        
        dist_plane[i] = np.linalg.norm(positions[b]-PointPrime)
        dist_value[i] = np.linalg.norm(positions[b]-PointPrime)

    dist_plane= np.argsort(dist_plane)# sort distances from min to max value and save the id 
    # The closest planes are planes:


    if DISPLAY==1:
        print 'nearest neighbour planes found'
        print dist_plane[0], dist_plane[1]


    #### reconstruct a lien between desired position and Xmax
    # since xmax positions is not given in coordinates you reconstruct its positions by defining aline from a plane and this line has a length of the given distance of the simulated plane to Xmax: dist1 which is the distance of the first plane in the simulations file
    # dist1 belongs to positions_sims[0,:], normal should always be the same


    #n= np.cross(positions_sims[0,10]-positions_sims[0,11], positions_sims[0,40]-positions_sims[0,41])
    #n= n/np.linalg.norm(n)
    # NOTE: in pricipal v and ne should be equal, but isnt at the moment
    v = np.array([np.cos(az)*np.sin(zen)*-1,np.sin(az)*np.sin(zen)*-1.,np.cos(zen)*-1.])#
    #v = np.array([np.cos(az)*np.sin(zen),np.sin(az)*np.sin(zen),np.cos(zen)])

    v= v/np.linalg.norm(v)

    s= dist1/np.linalg.norm(v) # in principal value should be one, but in any case calculate the absolute value
    p= [np.mean(positions_sims[0,:,0]), np.mean(positions_sims[0,:,1]), np.mean(positions_sims[0,:,2])] ## center of plane 1
    Xmax_pos= p- s* (v) # assuming that Xmax is "before" the planes
    
    if DISPLAY==1:
        print 'shower direction'
        print  v

        print dist1, np.linalg.norm(v)
        print 'postion Xmax, position desired'
        print Xmax_pos, positions[b]

        
        
    ## now you can construct a line given by Xmax_pos and your disired antenna positions. the Intersection points of this line with the planes gives you then the position for the pulseshape interpolation
    # plane is given by (point-p_0)*n=0, line given by point=s*(Xmax_pos-positions)+ Xmax_pos
    s0= np.dot( ( positions_sims[dist_plane[0],10]- Xmax_pos ), v)/ np.dot( (Xmax_pos-positions[b]), v) # intersection Point on plane dist_plane[0]
    Inter_plane0 = s0*  (Xmax_pos-positions[b]) + Xmax_pos
    s1= np.dot( ( positions_sims[dist_plane[1],10]- Xmax_pos ), v)/ np.dot( (Xmax_pos-positions[b]), v) # intersection Point on plane dist_plane[1]
    Inter_plane1 = s1*  (Xmax_pos-positions[b]) + Xmax_pos
    #print Inter_plane0, Inter_plane1
    ##############

    #if DISPLAY==1:


        ### Plot to check whether its working correctly
        #fig = plt.figure(1, facecolor='w', edgecolor='k')
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(positions_sims[dist_plane[0],:,0], positions_sims[dist_plane[0],:,1], positions_sims[dist_plane[0],:,2], c='red', marker='o', label="surrounding planes")
        #ax.scatter(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2], c='red', marker='o')
        ###ax.scatter(positions_sims[dist_plane[2],:,0], positions_sims[dist_plane[2],:,1], positions_sims[dist_plane[2],:,2], c='red', marker='o')
        ###ax.scatter(positions_sims[dist_plane[3],:,0], positions_sims[dist_plane[3],:,1], positions_sims[dist_plane[3],:,2], c='red', marker='o')
        ##ax.scatter(line[:,0],line[:,1],line[:,2],c='green', marker='o', lw = 0)# c='green', marker='+', s=80)
        ##ax.scatter(line_ortho[:,0],line_ortho[:,1],line_ortho[:,2], c='black', marker='o', lw = 0 )# c='green', marker='+', s=80)
        #ax.scatter(positions[b,0], positions[b,1], positions[b,2], c='blue', marker='x', label='desired position', s=80 )
        ##ax.scatter(test[0], test[1], test[2], c='green', marker='x', s=80) # orthogonal projection of point as test
        #ax.scatter(Xmax_pos[0], Xmax_pos[1], Xmax_pos[2], c='green', marker='x', label='Xmax positions' , s=80)
        #ax.scatter(Inter_plane0[0], Inter_plane0[1], Inter_plane0[2], c='green', marker='o', label='projection on planes' , s=80)
        #ax.scatter(Inter_plane1[0], Inter_plane1[1], Inter_plane1[2], c='green', marker='o', s=80 )
        #plt.legend(loc='upper right')
        #pylab.tight_layout(0.4, 0.5,1.0)
    
        #plt.show()

    #plt.close()
    


    ################## PulseShape Interpolation part
    
    ##rotate into shower coordninates
    
    #print positions_sims[dist_plane[0],:, 2], len(positions_sims[dist_plane[0],:])
    
    # find the 4 closest neighbours
    
    offinz= np.mean(positions_sims[dist_plane[0],:,2])
    offiny= np.mean(positions_sims[dist_plane[0],:,1])
    offinx= np.mean(positions_sims[dist_plane[0],:,0])
    pos_0= np.zeros([len(positions_sims[dist_plane[0],:]),3])
    
    Inter_0= GetUVW(Inter_plane0,offinx,offiny,offinz, zen, az, phigeo, thetageo) 
    Inter_0[0] = 0. # v-comp. correct for rounding errors or so   
    radius=np.linalg.norm(Inter_0)
 
    d0=np.zeros([4], dtype=int) 
    for i in np.arange(0,len(positions_sims[dist_plane[0],:])):
            pos_0[i,:] = GetUVW(positions_sims[dist_plane[0],i],offinx,offiny,offinz, zen, az, phigeo, thetageo) 
            pos_0[i,0] = 0. # v-comp. correct for rounding errors or so 
            if i==0:
                angle= np.arccos(np.dot(Inter_0, pos_0[0,:] ) / (np.linalg.norm(Inter_0) * np.linalg.norm(pos_0[0,:] )))
                if Inter_0[2] <0.: # vxvxB <0 -- plus pi in angle since function always return the smaller one
                    angle=angle+ np.pi
                #print angle
            if radius > np.linalg.norm(pos_0[i,:]) : # r > r1
               if i==0:
                   angle1 =0.
               #print radius, np.linalg.norm(pos_0[i,:])
               else:
                   angle1=np.arccos(np.dot(pos_0[i,:], pos_0[0,:] ) / (np.linalg.norm(pos_0[i,:]) * np.linalg.norm(pos_0[0,:] )))
                   if i%8>3:
                        angle1=2* np.pi- angle1
               if angle > ( angle1 ):   # look for clostest alpha, pos[0,:] refernce antenna
                   #print angle, np.arccos(np.dot(pos_0[i,:], pos_0[0,:] ) / (np.linalg.norm(pos_0[i,:]) * np.linalg.norm(pos_0[0,:] )))
                   d0[0]=(i) # antenna which radius is smaller and angle smaller
                   d0[1] =(i+1)# # antenna which radius is smaller and angle larger
                   d0[2]=(i+8)# antenna which radius is larger and angle smaller, based on 8 angles
                   d0[3] =(i+8+1)# # antenna which radius is larger and angle larger
                   
 
# plane 2 
    offinz= np.mean(positions_sims[dist_plane[1],:,2])
    offiny= np.mean(positions_sims[dist_plane[1],:,1])
    offinx= np.mean(positions_sims[dist_plane[1],:,0])
    pos_1= np.zeros([len(positions_sims[dist_plane[1],:]),3])
    
    Inter_1= GetUVW(Inter_plane1,offinx,offiny,offinz, zen, az, phigeo, thetageo) 
    Inter_1[0] = 0. # v-comp. correct for rounding errors or so   
    radius=np.linalg.norm(Inter_1)
 
    d1=np.zeros([4], dtype=int) 
    for i in np.arange(0,len(positions_sims[dist_plane[1],:])):
            pos_1[i,:] = GetUVW(positions_sims[dist_plane[1],i],offinx,offiny,offinz, zen, az, phigeo, thetageo) 
            pos_1[i,0] = 0. # v-comp. correct for rounding errors or so 
            if i==0:
                angle= np.arccos(np.dot(Inter_1, pos_1[0,:] ) / (np.linalg.norm(Inter_1) * np.linalg.norm(pos_1[0,:] )))
                if Inter_0[2] <0.: # vxvxB <0 -- plus pi in angle since function always return the smaller one
                    angle=angle+ np.pi
                #print angle
            if radius > np.linalg.norm(pos_1[i,:]) : # r > r1
               if i==0:
                   angle1 =0.
               #print radius, np.linalg.norm(pos_0[i,:])
               else:
                   angle1=np.arccos(np.dot(pos_1[i,:], pos_1[0,:] ) / (np.linalg.norm(pos_1[i,:]) * np.linalg.norm(pos_1[0,:] )))
                   if i%8>3:
                        angle1=2* np.pi- angle1
               if angle > ( angle1 ):   # look for clostest alpha, pos[0,:] refernce antenna
                   #print angle, np.arccos(np.dot(pos_0[i,:], pos_0[0,:] ) / (np.linalg.norm(pos_0[i,:]) * np.linalg.norm(pos_0[0,:] )))
                   d1[0]=(i) # antenna which radius is smaller and angle smaller
                   d1[1] =(i+1)# # antenna which radius is smaller and angle larger
                   d1[2]=(i+8)# antenna which radius is larger and angle smaller, based on 8 angles
                   d1[3] =(i+8+1)# # antenna which radius is larger and angle larger
                   
            #print angle, radius, '    ', i,  "angles, radius", angle1, np.linalg.norm(pos_0[i,:])
                   
    if d0.any()>120. or d1.any() >120:
        print "########  desired antenna position outside region in which interpolation works, no 4 neighbours.... antenna skipped"
        continue
    try:
        print pos_0[d0[3]]    
    except IndexError:
        print "########  desired antenna position outside region in which interpolation works, no 4 neighbours.... antenna skipped"
        continue
    try:
        print pos_1[d1[3]]    
    except IndexError:
        print "######## desired antenna position outside region in which interpolation works, no 4 neighbours.... antenna skipped"
        continue
    

 

 

    if DISPLAY==1:
        
        
            print d0, d1
                ### Plot to check whether its working correctly
            fig2 = plt.figure(2, facecolor='w', edgecolor='k')
            ax2 = fig2.add_subplot(111)#, projection='3d')
            ax2.scatter( pos_0[:,1], pos_0[:,2], c='red', marker='o', label="surrounding planes")#pos_0[:,0],
            ##ax.scatter(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2], c='red', marker='o')
            ax2.plot(pos_0[:,1], pos_0[:,2])#pos_0[:,0], 
            

            ax2.scatter( pos_0[d0[0],1], pos_0[d0[0],2], c='blue', marker='x', s=80)#pos_0[d0[0],0],
            ax2.scatter( pos_0[d0[1],1], pos_0[d0[1],2], c='blue', marker='x', s=80)#pos_0[d0[1],0],        
            ax2.scatter( pos_0[d0[2],1], pos_0[d0[2],2], c='blue', marker='x', s=80)#pos_0[d0[2],0],       
            ax2.scatter( pos_0[d0[3],1], pos_0[d0[3],2], c='blue', marker='x', s=80) #  pos_0[d0[3],0],    
            
            ax2.scatter( Inter_0[1], Inter_0[2], c='green', marker='+', label='projection on planes' , s=80)#Inter_0[0],
            plt.legend(loc='upper right')
            plt.tight_layout(0.4, 0.5,1.0)
            plt.axis('equal')
            
            plt.xlabel(r"vxB", fontsize=16)
            plt.ylabel(r"vxvxB", fontsize=16)
            
            
            
            #ax3 = fig2.add_subplot(132)#, projection='3d')
            #ax3.scatter( pos_0[:,0], pos_0[:,1], c='red', marker='o', label="surrounding planes")#pos_0[:,0],
            ##ax.scatter(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2], c='red', marker='o')
            #ax3.plot(pos_0[:,0], pos_0[:,1])#pos_0[:,0], 
            

            #ax3.scatter( pos_0[d0[0],0], pos_0[d0[0],1], c='blue', marker='x', s=80)#pos_0[d0[0],0],
            #ax3.scatter( pos_0[d0[1],0], pos_0[d0[1],1], c='blue', marker='x', s=80)#pos_0[d0[1],0],        
            #ax3.scatter( pos_0[d0[2],0], pos_0[d0[2],1], c='blue', marker='x', s=80)#pos_0[d0[2],0],       
            #ax3.scatter( pos_0[d0[3],0], pos_0[d0[3],1], c='blue', marker='x', s=80) #  pos_0[d0[3],0],    
            #ax3.scatter( Inter_0[0], Inter_0[1], c='green', marker='o', label='projection on planes' , s=80)#Inter_0[0],

            #plt.xlabel(r"v", fontsize=16)
            #plt.ylabel(r"vxv", fontsize=16)
            
            #ax3 = fig2.add_subplot(133)#, projection='3d')
            #ax3.scatter( pos_0[:,0], pos_0[:,2], c='red', marker='o', label="surrounding planes")#pos_0[:,0],
            ##ax.scatter(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2], c='red', marker='o')
            #ax3.plot(pos_0[:,0], pos_0[:,2])#pos_0[:,0], 
            

            #ax3.scatter( pos_0[d0[0],0], pos_0[d0[0],2], c='blue', marker='x', s=80)#pos_0[d0[0],0],
            #ax3.scatter( pos_0[d0[1],0], pos_0[d0[1],2], c='blue', marker='x', s=80)#pos_0[d0[1],0],        
            #ax3.scatter( pos_0[d0[2],0], pos_0[d0[2],2], c='blue', marker='x', s=80)#pos_0[d0[2],0],       
            #ax3.scatter( pos_0[d0[3],0], pos_0[d0[3],2], c='blue', marker='x', s=80) #  pos_0[d0[3],0],    
            #ax3.scatter( Inter_0[0], Inter_0[2], c='green', marker='o', label='projection on planes' , s=80)#Inter_0[0],

            #plt.xlabel(r"v", fontsize=16)
            #plt.ylabel(r"vxvxB", fontsize=16)
        
            plt.show()


                
                
     
    #if DISPLAY==1:
        #print '\n cloest antennas on ecach plane, Plane 1 and Plane 2'
        #print d0[0], d0[1], d0[2], d0[3]
        #print d1[0], d1[1], d1[2], d1[3]


        ### Plot to check whether its working correctly
        #fig = plt.figure(1, facecolor='w', edgecolor='k')
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(positions_sims[dist_plane[0],:,0], positions_sims[dist_plane[0],:,1], positions_sims[dist_plane[0],:,2], c='red', marker='o', label="surrounding planes")
        #ax.scatter(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2], c='red', marker='o')
        ##ax.plot(positions_sims[dist_plane[1],:,0], positions_sims[dist_plane[1],:,1], positions_sims[dist_plane[1],:,2])
        
        
        
        ###ax.scatter(positions_sims[dist_plane[2],:,0], positions_sims[dist_plane[2],:,1], positions_sims[dist_plane[2],:,2], c='red', marker='o')
        ###ax.scatter(positions_sims[dist_plane[3],:,0], positions_sims[dist_plane[3],:,1], positions_sims[dist_plane[3],:,2], c='red', marker='o')
        ##ax.scatter(line[:,0],line[:,1],line[:,2],c='green', marker='o', lw = 0)# c='green', marker='+', s=80)
        ##ax.scatter(line_ortho[:,0],line_ortho[:,1],line_ortho[:,2], c='black', marker='o', lw = 0 )# c='green', marker='+', s=80)
        #ax.scatter(positions[b,0], positions[b,1], positions[b,2], c='blue', marker='o', label='desired position', s=80 )
        #ax.scatter(positions_sims[dist_plane[0],d0[0]][0], positions_sims[dist_plane[0],d0[0]][1], positions_sims[dist_plane[0],d0[0]][2], c='blue', marker='x', s=80)
        #ax.scatter(positions_sims[dist_plane[0],d0[1]][0], positions_sims[dist_plane[0],d0[1]][1], positions_sims[dist_plane[0],d0[1]][2], c='blue', marker='x', s=80)        
        #ax.scatter(positions_sims[dist_plane[0],d0[2]][0], positions_sims[dist_plane[0],d0[2]][1], positions_sims[dist_plane[0],d0[2]][2], c='blue', marker='x', s=80)       
        #ax.scatter(positions_sims[dist_plane[0],d0[3]][0], positions_sims[dist_plane[0],d0[3]][1], positions_sims[dist_plane[0],d0[3]][2], c='blue', marker='x', s=80)     
        #ax.scatter(positions_sims[dist_plane[1],d1[0]][0], positions_sims[dist_plane[1],d1[0]][1], positions_sims[dist_plane[1],d1[0]][2], c='blue', marker='x', s=80)
        #ax.scatter(positions_sims[dist_plane[1],d1[1]][0], positions_sims[dist_plane[1],d1[1]][1], positions_sims[dist_plane[1],d1[1]][2], c='blue', marker='x', s=80)        
        #ax.scatter(positions_sims[dist_plane[1],d1[2]][0], positions_sims[dist_plane[1],d1[2]][1], positions_sims[dist_plane[1],d1[2]][2], c='blue', marker='x', s=80)       
        #ax.scatter(positions_sims[dist_plane[1],d1[3]][0], positions_sims[dist_plane[1],d1[3]][1], positions_sims[dist_plane[1],d1[3]][2], c='blue', marker='x', s=80)   
        
        #ax.scatter(Xmax_pos[0], Xmax_pos[1], Xmax_pos[2], c='green', marker='x', label='Xmax positions' , s=80)
        #ax.scatter(Inter_plane0[0], Inter_plane0[1], Inter_plane0[2], c='green', marker='o', label='projection on planes' , s=80)
        #ax.scatter(Inter_plane1[0], Inter_plane1[1], Inter_plane1[2], c='green', marker='o', s=80 )
        #plt.legend(loc='upper right')
        #plt.tight_layout(0.4, 0.5,1.0)
    
        #plt.show()
           
                
                

    if DISPLAY==1:

        print '\n\n PLANE1'

    ## PLANE 1

    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '

        print '\n Interpolate x'
    
    


    point_online1=ProjectPointOnLine(positions_sims[dist_plane[0],d0[0]], positions_sims[dist_plane[0],d0[1]], Inter_plane0)# Project Point on line 1
    if DISPLAY==1:
        print positions_sims[dist_plane[0],d0[0]], positions_sims[dist_plane[0],d0[1]], point_online1
    if full==1:
        #txt0=np.loadtxt(path1+'/scaled_'+str(sims[dist_plane[0]]) +'/a'+str(d0[0])+'.trace')
        #txt1=np.loadtxt(path1 +'/scaled_'+str(sims[dist_plane[0]]) +'/a'+str(d0[1])+'.trace')
        txt0=np.loadtxt(path1+str(sims[dist_plane[0]]) +'/a'+str((d0[0]))+'.trace')
        txt1=np.loadtxt(path1 +str(sims[dist_plane[0]]) +'/a'+str((d0[1]))+'.trace')
    else:
        txt0=np.genfromtxt(path1+str(sims[dist_plane[0]]) +"/a"+str((d0[0]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
        txt1=np.genfromtxt(path1+str(sims[dist_plane[0]]) +"/a"+str((d0[1]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[1], positions_sims[dist_plane[0],d0[0]] , txt1.T[0], txt1.T[1], positions_sims[dist_plane[0],d0[1]], point_online1 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2

    point_online2=ProjectPointOnLine(positions_sims[dist_plane[0],d0[2]], positions_sims[dist_plane[0],d0[3]], Inter_plane0)# Project Point on line 2
    if DISPLAY==1:
        print '\n\n Projection 2 '
        print positions_sims[dist_plane[0],d0[2]], positions_sims[dist_plane[0],d0[3]], point_online2

    if full==1:
        #txt2=np.loadtxt(path1 +'/scaled_'+str(sims[dist_plane[0]]) +'/a'+str(d0[2])+'.trace')
        #txt3=np.loadtxt(path1 +'/scaled_'+str(sims[dist_plane[0]]) +'/a'+str(d0[3])+'.trace')
        txt2=np.loadtxt(path1 +str(sims[dist_plane[0]]) +'/a'+str((d0[2]))+'.trace')
        txt3=np.loadtxt(path1 +str(sims[dist_plane[0]]) +'/a'+str((d0[3]))+'.trace')
    else:
        txt2=np.genfromtxt(path1+str(sims[dist_plane[0]]) +"/a"+str((d0[2]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
        txt3=np.genfromtxt(path1+str(sims[dist_plane[0]]) +"/a"+str((d0[3]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)

    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[1], positions_sims[dist_plane[0],d0[2]] , txt3.T[0], txt3.T[1], positions_sims[dist_plane[0],d0[3]], point_online2  ,upsampling=True) #switch on upsamling by factor 8

    if DISPLAY==1:
        print '\n interpolation plane 1'
    ##### Get the pulse shape of the desired position (projection on plane0) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planex0, tracedes_planex0 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online1, xnew2, tracedes2, point_online2, Inter_plane0) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)



    if DISPLAY==1:
        print '\n Interpolate y'


    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[2], positions_sims[dist_plane[0],d0[0]] , txt1.T[0], txt1.T[2], positions_sims[dist_plane[0],d0[1]], point_online1 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2 ---- some wrong
    if DISPLAY==1:
        print '\n\n Projection 2 '
    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[2], positions_sims[dist_plane[0],d0[2]] , txt3.T[0], txt3.T[2], positions_sims[dist_plane[0],d0[3]], point_online2  ,upsampling=True) #switch on upsamling by factor 8

    if DISPLAY==1:
        print '\n interpolation plane 1'
    ##### Get the pulse shape of the desired position (projection on plane0) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planey0, tracedes_planey0 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online1, xnew2, tracedes2, point_online2, Inter_plane0) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)


    if DISPLAY==1:
        print '\n Interpolate z'


    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[3], positions_sims[dist_plane[0],d0[0]] , txt1.T[0], txt1.T[3], positions_sims[dist_plane[0],d0[1]], point_online1 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2
    if DISPLAY==1:
        print '\n\n Projection 2 '
    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[3], positions_sims[dist_plane[0],d0[2]] , txt3.T[0], txt3.T[3], positions_sims[dist_plane[0],d0[3]], point_online2  ,upsampling=True) #switch on upsamling by factor 8

    
    if DISPLAY==1:
     print '\n interpolation plane 1'
    ##### Get the pulse shape of the desired position (projection on plane0) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planez0, tracedes_planez0 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online1, xnew2, tracedes2, point_online2, Inter_plane0) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)








    if DISPLAY==1:
        print '\n\n PLANE2'

    ## PLANE 2


        print '\n Interpolate x'

    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '


    point_online12=ProjectPointOnLine(positions_sims[dist_plane[1],d1[0]],positions_sims[dist_plane[1], d1[1]], Inter_plane1)# Project Point on line 1
    if DISPLAY==1:
        print positions_sims[dist_plane[1],d1[0]], positions_sims[dist_plane[1],d1[1]], point_online12
    if full==1:
        #txt0=np.loadtxt(path1+ '/scaled_'+str(sims[dist_plane[1]]) +'/a'+str(d1[0])+'.trace')
        #txt1=np.loadtxt(path1+ '/scaled_'+str(sims[dist_plane[1]]) +'/a'+str(d1[1])+'.trace')
        txt0=np.loadtxt(path1+ str(sims[dist_plane[1]]) +'/a'+str((d1[0]))+'.trace')
        txt1=np.loadtxt(path1+ str(sims[dist_plane[1]]) +'/a'+str((d1[1]))+'.trace')
    else:
        txt0=np.genfromtxt(path1+str(sims[dist_plane[1]]) +"/a"+str((d1[0]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
        txt1=np.genfromtxt(path1+str(sims[dist_plane[1]]) +"/a"+str((d1[1]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[1], positions_sims[dist_plane[1],d1[0]] , txt1.T[0], txt1.T[1], positions_sims[dist_plane[1],d1[1]], point_online12 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2
    if DISPLAY==1:
        print '\n\n Projection 2 '
    point_online22=ProjectPointOnLine(positions_sims[dist_plane[1],d1[2]], positions_sims[dist_plane[1],d1[3]], Inter_plane1)# Project Point on line 2
    if DISPLAY==1:
        print positions_sims[dist_plane[1],d1[2]], positions_sims[dist_plane[1],d1[3]], point_online22

    if full==1:
        #txt2=np.loadtxt(path1+ '/scaled_'+str(sims[dist_plane[1]]) +'/a'+str(d1[2])+'.trace')
        #txt3=np.loadtxt(path1+ '/scaled_'+str(sims[dist_plane[1]]) +'/a'+str(d1[3])+'.trace')
        txt2=np.loadtxt(path1+ str(sims[dist_plane[1]]) +'/a'+str((d1[2]))+'.trace')
        txt3=np.loadtxt(path1+ str(sims[dist_plane[1]]) +'/a'+str((d1[3]))+'.trace')
    else:
        txt2=np.genfromtxt(path1+str(sims[dist_plane[1]]) +"/a"+str((d1[4]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)
        txt3=np.genfromtxt(path1+str(sims[dist_plane[1]]) +"/a"+str((d1[5]))+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat', skip_footer=2)

    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[1], positions_sims[dist_plane[1],d1[2]] , txt3.T[0], txt3.T[1], positions_sims[dist_plane[1],d1[3]], point_online22 ,upsampling=True) #switch on upsamling by factor 8

    if DISPLAY==1:
        print '\n interpolation plane 2'
    ##### Get the pulse shape of the desired position (projection on plane1) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planex1, tracedes_planex1 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online12, xnew2, tracedes2, point_online22, Inter_plane1 ) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)


    if DISPLAY==1:
        print '\n Interpolate y'

    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[2], positions_sims[dist_plane[1],d1[0]] , txt1.T[0], txt1.T[2], positions_sims[dist_plane[1],d1[1]], point_online12 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2
    if DISPLAY==1:
        print '\n\n Projection 2 '
    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[2], positions_sims[dist_plane[1],d1[2]] , txt3.T[0], txt3.T[2], positions_sims[dist_plane[1],d1[3]], point_online22 ,upsampling=True) #switch on upsamling by factor 8

    if DISPLAY==1:
        print '\n interpolation plane 2'
    ##### Get the pulse shape of the desired position (projection on plane1) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planey1, tracedes_planey1 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online12, xnew2, tracedes2, point_online22, Inter_plane1 ) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)



    if DISPLAY==1:
        print '\n Interpolate z'

    ## Get the pulseshape for the projection on line 1
        print ' Projection 1 '
    ## the interpolation of the pulse shape is performed  
    xnew1, tracedes1 =pulse.Interpolate_PulseShape(txt0.T[0], txt0.T[3], positions_sims[dist_plane[1],d1[0]] , txt1.T[0], txt1.T[3], positions_sims[dist_plane[1],d1[1]], point_online12 ,upsampling=True) #switch on upsamling by factor 8


    ### Get the pulseshape for the projection on line 2
    if DISPLAY==1:
        print '\n\n Projection 2 '
    ## the interpolation of the pulse shape is performed  
    xnew2, tracedes2 =pulse.Interpolate_PulseShape(txt2.T[0], txt2.T[3], positions_sims[dist_plane[1],d1[2]] , txt3.T[0], txt3.T[3], positions_sims[dist_plane[1],d1[3]], point_online22 ,upsampling=True) #switch on upsamling by factor 8

    if DISPLAY==1:
        print '\n interpolation plane 2'
    ##### Get the pulse shape of the desired position (projection on plane1) from projection on line1 and 2
    #print ' Pulse Shape '
    xnew_planez1, tracedes_planez1 =pulse.Interpolate_PulseShape(xnew1, tracedes1, point_online12, xnew2, tracedes2, point_online22, Inter_plane1 ) #(t1, trace1, x1, t2, trace2, x2, xdes, path, nrdes)





    if DISPLAY==1:
        print '\n\n final interpolation'
    xnew_desiredx, tracedes_desiredx =pulse.Interpolate_PulseShape(xnew_planex0, tracedes_planex0, Inter_plane0 ,xnew_planex1, tracedes_planex1, Inter_plane1, positions[b])

    xnew_desiredy, tracedes_desiredy =pulse.Interpolate_PulseShape(xnew_planey0, tracedes_planey0, Inter_plane0 ,xnew_planey1, tracedes_planey1, Inter_plane1, positions[b])

    xnew_desiredz, tracedes_desiredz =pulse.Interpolate_PulseShape(xnew_planez0, tracedes_planez0, Inter_plane0 ,xnew_planez1, tracedes_planez1, Inter_plane1, positions[b])




    if DISPLAY==1:
        print ' length of time traces: ', len(txt2.T[0]), len(xnew_desiredx)

    #hexf = abs(hilbert(tracedes_desiredx))
    #heyf = abs(hilbert(tracedes_desiredy))
    #hezf = abs(hilbert(tracedes_desiredz))
    #exm = max(hexf)
    #eym = max(heyf)
    #ezm = max(hezf)
    #amp = sqrt(exm*exm+ eym*eym+ezm*ezm)
    
    


    ## NOTE: This traces should be saved in some way similar to a#if.trace for and easier inclusion into the filtering process
    print ' interpolated signal belonging to positions in ' +str(path0) +' saved as '

    #### lop over b as number of desired positions
    if full==1:
        name=path2+ '/a'+str(b)+'.trace'
        print name
    else:
        name=path2+ "/a"+str(b)+'_'+str((f1*1e-6)) + '-' + str((f2*1e-6)) + 'MHz.dat'
        # add the Hilbert envelope later.... but just possible if all three components interpolated ### add thsi to name
        print name

    FILE = open(name, "w+" )
    for i in range( 0, len(xnew_desiredx) ):
            
        #print >>FILE,"%3.2f	%1.5e	%1.5e	%1.5e	%1.5e	%1.5e	%1.5e" % (txt.T[0][i], Ex[i], Ey[i], Ez[i], Ev[i], EvxB[i], EvxvxB[i] )
            print >>FILE,"%3.2f %1.5e %1.5e %1.5e" % (xnew_desiredx[i], tracedes_desiredx[i], tracedes_desiredy[i], tracedes_desiredz[i])

        
        
    ## in the last line of the file we wanna write the max. ampl of the hilbert envelope    
    ## something like: amp exm eym ezm
    #print >>FILE, ""
    ##print >>FILE,"%1.5f	%1.5e	%1.5e	%1.5e	%1.5f	%1.5e	%1.5e	%1.5e" % (amp, exm, eym, ezm, amp2,exm2, eym2, ezm2)
    #print >>FILE,"%1.5f %1.5e %1.5e %1.5e" % (amp, exm, eym, ezm)


    FILE.close()

    if DISPLAY==1:
        #### PLOTTING
        #fig2=plt.figure(2)
        #plt.plot(txt0.T[0], txt0.T[2], 'b--', label= "first")
        #plt.plot(txt1.T[0], txt1.T[2], 'r--', label= "second")
        #plt.plot(xnew1, np.real(tracedes1), 'g--', label= "interpolated1")

        #plt.plot(txt2.T[0], txt2.T[2], 'b:', label= "third")
        #plt.plot(txt3.T[0], txt3.T[2], 'r:', label= "fourth")
        fig2 = plt.figure(2, facecolor='w', edgecolor='k')
        plt.plot(xnew_planey0, np.real(tracedes_planey0), 'g:', label= "plane 0")

        plt.plot(xnew_planey1, np.real(tracedes_planey1), 'b:', label= "plane 1")

        plt.plot(xnew_desiredy, np.real(tracedes_desiredy), 'r-', label= "desired")

        #plt.plot(txt2.T[0], txt_test.T[2], 'c:', label= "real")
            #plt.plot(txt2.T[0], Amplitude, 'b--')
        plt.xlabel(r"time (s)", fontsize=16)
        plt.ylabel(r"Amplitude muV/m ", fontsize=16)
        plt.legend(loc='best')

        plt.show()
