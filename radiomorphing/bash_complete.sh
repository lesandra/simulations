#!/bin/bash

#### TODO to be substituted by script reading in DANTON files
#
# parameters of target shower: AzimuthScan4_26
E=0.96 #in EeV 
injh=2000 # in m above ground
zen=89.5 # =0 90.5 in zhaires# in deg (aires convention) ----> grand
az=0 #==-180  in zhaires  # in deg (aires convention) ----> grand
alpha=3 # in deg, mountain slope for antenna repsonse
primary=electron


fullpath=/media/sf_work/Paris/scripts_GRAND/ 
#### information about the reference shower
file_des_pos=${fullpath}/Olivier_scripts/InterpolatedSignals/antpos_desired2.dat  # file containing desired antenna position in m, at least 2
file_sim=${fullpath}/test/Simulations.dat # file containing runnames/numbers which should be used = star shape pattern planes
# path_sim=/media/sf_work/Paris/scripts_GRAND/Olivier_scripts/EvetADetailed2/ # path to needed sim (see file_sim)
RUN=GrandEventADetailed2/
SIMDIR=${fullpath}/Olivier_scripts/${RUN}/ # path to needed sim (see file_sim)

### folder where the scaled results can be found
path_inter=${fullpath}/Olivier_scripts/InterpolatedSignals/ # path to folder where integrated traces at desired positions should be saved



##### START ######

### TODO: include mountain shadowing for the desired antenna psoitions

##################### loop over the single stare shape plane belonging to one REFERENCE shower

    ##### 0st step: PRE-PROCESSING: this has to be done just once and can then be comment out

    #####Splits up ZHAireS output in single antenna files containing electic field traces + an antenna file
    # TODO: reduce numbers of digits to .1nuV/m, remove the first zeroes from all traces, study the effect of reducing digits
    
filename=${file_sim}
while read line || [[ -n $line ]];
# while read -r line
do
    name="$line"
    echo "Name read from file - $name"
    RUNSIM=$name
#     echo $RUNSIM
    
    
    FILE=${SIMDIR}/${RUNSIM}/a0.trace

    if [ ! -f "$FILE" ]
    then
        echo "File $FILE does not exists. Create...."
        python splitZhairesFields_all.py ${SIMDIR}/${RUNSIM}/ ${RUNSIM}
    fi
    


    
    
    
##### 1st step: Scale the reference shower to the desired parameters
                                     
    #remove the comparison in pulse shape as well freq dependency
    # hand over parameters to which you would like to scale to
    ### scale the traces of all star shapes which will be used for the interpoaltion to the parameters of the desired shower A--->B via scaling and stretching in shower coordinates and backtransformation, gives you instead of ExyzA at xi---> ExyzB at xi
    #arguments: path to reference shower, name of ref. sim, energy target in EeV, injection height target refereeing to obs level in m, zen target in deg(Aires), azimuth target in deg(Aires), identification number of antenna in star shape
    # saves scaled taces in scaled_${RUNSIM}

    # TODO: Are there any preferred shower direction, several shower parameters to run
    # TODO: add the correct position vector for scaled plane, so far angles of target and reference has to be the same
     python  PulseShape_Scaling.py ${SIMDIR}/ ${RUNSIM} ${E} ${injh} ${zen} ${az} ${primary} #${l} 

done < "$filename"    
### loop end    


    
#### 2nd step: interpolation of the traces : #scaled traces have to be read in

# # NOTE: the scaled traces of the planes are saved in another folder: scaled_Planename, it is not SIMDIR/scaled_...
###### Interpolation: read in the scaled reference shower planes and a list of the desired antenna positions, save the interpolated traces

## interpolated traces of desired antennas positions via reading in a bunch of star shape planes (change to scaled), gives you Exyz at desired antenna position x
# TODO: correct for time shift due to phase shift and then do the geometric time shift, so that Valentin can use it for wavefront timing
# --- just the timing between the antennas has to somehow
# TODO: check for the best distances between antenna positions along the rays and distances between single star-shape planes
# NOTE: conversion from aires to GRAND convention if less thean 5 arguments, then angles read in from simu, otherwise hand over angles in GRAND conv.
 python PulseShape.py ${file_des_pos} ${file_sim} ${SIMDIR} ${path_inter} ${zen} ${az}




###### NOTE: In principle here Nicolas' script can be used.... antenna response, trigger rate etc.

# # #### 3rd sstep: Apply antenna response : read interpolated traces in Exyz, returns Etheta,phi or the power or voltage respectively, antenna sensitive to 20-300MHz, EW+NS
# # # cd in antenna response folder
# # cd ./antennaresponse_2/antennaresponse_work/
# # 
# # 
# # #TODO: read in file of desired antenna positions and handover the position and antenna ID automatically
# # # loop over s= the number of desired antenna positions = len(antpos_desired2.dat)
# #     # python computevoltage zen(referring to aires zen) az((referring to aires az) path to traces(file containing Exyz)   ID number of antenna (a#.trace)
# #     python computevoltage.py ${zen} ${az} ${alpha} ${path_inter} 0 1 43499.34 692.75  2379.61
# #     # returns timet in s, voltage in muV,vt,vp into ew_out_a#.
# #     # internal converion from aires to antenna response angles included
# #     # Done: check for timing multiply by 10^-3 or 10^3, hand over ns
# #     # TODO: antenna response for vertical needed
# #     # DONE: hand over GRAND coordinates, internally converted to NEC conv
# #     # TODO: read in desired antenna positions and read them in
# #     python computevoltage.py ${zen} ${az} ${alpha} ${path_inter} 1 1 43499.34 1187.58 2379.61
# # 
# #     
