

RUN=AzimuthScan4
RUNSIM=${RUN}_20
SIMDIR=PATH/${RUN}/${RUNSIM}/ #/media/sf_work/Paris/scripts_GRAND/ # folder where your sims are saved


FREQ1=60
FREQ2=200


#####Splits up ZHAireS output in single antenna files containing electic field traces + an antenna file

FILE=${SIMDIR}/a0.trace

if [ ! -f "$FILE" ]
then
    echo "File $FILE does not exists. Create...."
    python splitZhairesFields_all.py ${SIMDIR} ${RUNSIM}
fi




# # # this script has to loop over all antenna files
# # # does filtering, procudes Hilbert envelope and gets the max, creates file containing filtered efield + hilbert max, does plotting
# # # you can hand over the range of the frequency band as parameter 3 and 4 in MHz,
# # # it uses a butterworth band pass filter of 6th order
FILE=${SIMDIR}/a0_${FREQ1}-${FREQ2}MHz.dat

if [ ! -f "$FILE" ]
then
    echo "File $FILE does not exists"
    CODE=0
    #at the moment #ant=120 is hardcoded, has to be substituted at some point
    for (( i=0; i < 120; i++ )); do
        python processSim.py ${SIMDIR} ${CODE} ${FREQ1} ${FREQ2} ${RUNSIM}
        let CODE+=1
    done
fi