
#==============================================================================

# Script created to analyze Xmax GRAND simulations

#==============================================================================

SCRIPTDIR=/Users/guepin/Documents/GRAND/Xmax_new2/ # Folder where your scripts are saved


# GEANTDIR=/home/zilles/Geant4-SKA/LORA

# How to proceed with numbering the simulations is not consistent so far. If there is one number missing the script just skips it.
# For running the analysis script you have nevertheless to hand over a list file with the simuations you like to use

NR=CR190_77deg_flat
CODE=${NR}_0 #370350 #number of your first simulation
CODE2=${NR}_0 #370350 #the same number again


FILES=10 #in generell 50 protons + 20 iron
filename=SimlistTest190_flat.txt

################ STARTING THE ANALYSIS ################

  cd ${SCRIPTDIR}

  DATADIR=./CR-Sim/${NR} # the folder of all the simulations
  OTFD=./E${CODE2}_1000m/ # Outputfile folder
  


  
#### This has to be done ONE time after the simulations are finished: renames files and put inputfiles in an special folder
 for (( i=0; i < ${FILES}; i++ )); do
      SIMDIR=${DATADIR}/${NR}_${i}

    FILE=${SIMDIR}/a0.trace

    if [ ! -f "$FILE" ]
    then
        echo "File $FILE does not exists"
        python splitZhairesFields_all.py ${SIMDIR}/ ${NR}_${i}
    fi
done

  
  if ! [ -d "$OTFD" ]; then
 	# Control will enter here if $DIRECTORY exists.
 	  mkdir ${OTFD}
 
  fi
  

  MAPDIR=GRAND_antenna.list # List of GRAND antenna positions
  python ProcessData_wolora4_Zhaires.py ${DATADIR} ${CODE2} ${OTFD}/DAT${CODE2}-pickle.dat ${FILES} 50 200 ${OTFD} ${MAPDIR} ${filename}

### particle file not included
### lower freq:50 - high freq.:350  in MHz
