# simulations
all the scripts needed to set up the end-to-end simulation
Please, write as well a short description of the script at the top in each script and how to start it (whoch input is needed)

splitZhairesFields_all.py: reads in the timefresnel-root.dat file and splits it up into singles antenna files. It also produces a file with all antenna psoitions.
filters.py: contains some basic filters
processSim.py: Does the filtering and produce a file containing the filtered traces in xyz and in shower coordinates. It calculates as well the hilbert envelope and gets the maximum amplitudes from that for later analysis.
SimGRAND_Filtering.sh: Just a simple bash file to start the procedure
