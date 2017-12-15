#!/usr/bin/env python

"""
                    Script filter_danton
                        Version 1.0
    Written by N. Renault-Tinacci from a script provided by C. Medina
                Using danton.py developped by V. Niess
"""

import os, glob
import sys
import numpy as np
import pylab as pl
import danton
import StringIO

##########################################################################################################
def main():
    """ Main script allowing to filter a DANTON library according to some criterions """

    ##################################################################
    # Test arguments
    if (len(sys.argv)<2 or len(sys.argv)>2):
        print """\
    This script will allow to filter a DANTON library accordingly to some criterions. Four of them are implemented and listed here below.
    It is designed to filter the libraries once. Run more times will overwrite the previous filtered shower list.
    Currently implemented criterions : Minimum shower energy, minimum and maximum injection height, 
                                    minimum and maximum elevation and earth-skimming vs atmospheric neutrino

    Usage:  python filter_danton.py
    """
        sys.exit(1)
    ##################################################################
    # Retrieve the input parameters and set the criterions
    file=str(sys.argv[1]) #Initial DANTON library
    Emin=5e16 #minimum total energy of a shower in eV
    injh_min=0. #minimum injection height in m
    injh_max=13000. #maximum injection height in m
    theta_min=65. #zenith in GRAND convention (0 for horizontal showers and <90deg for up-going showers) in degrees
    theta_max=90. #zenith in GRAND convention (0 for horizontal showers and >90deg for down-going showers) in degrees
    depth_type=1 #1 for earth-skimming neutrino, -1 for atmospheric neutrino
    depth_min=0.

    # Parse the events, build up some statistics and perform the selection of the events.
    selevt_ID=[]
    selevt=[]
    Nevt=0
    for event in danton.iter_event(file):
            lastid = event.id
            for decay in event.decay:
                [depth,height,theta,et] = parse_build(event,decay)

            # Selecting the events meeting the criterions
            if ((theta<=theta_max) & (theta>=theta_min)) & (et>=Emin) & ((height>=injh_min) & (height<=injh_max)) & (depth>0. and np.sign(depth)==depth_type):
                selevt_ID.append(event.id)
                selevt.append(event)
                Nevt=Nevt+1
            elif ((theta<=theta_max) & (theta>=theta_min)) & (et>=Emin) & ((height>=injh_min) & (height<=injh_max)) & (depth<=0. and np.sign(depth-0.01)==depth_type):
                selevt_ID.append(event.id)
                selevt.append(event)
                Nevt=Nevt+1
                    
    selevt_ID=np.array(selevt_ID)
    print "+ {:} tau decays for {:} incoming neutrinos".format(len(selevt_ID), lastid+1)

    #Re-writing the DANTON library with only the filtered showers
    selevt_file=os.path.dirname(file)+'/'+os.path.splitext(os.path.basename(file))[0]+'_filtered.txt'
    print selevt_file
    evtfile = open(selevt_file,"w+")
    totito  = generate_input(selevt,selevt_ID)
    evtfile.write(totito)
    evtfile.close()
    
##########################################################################################################
##########################################################################################################
###                                  Let's define useful functions                                     ###
##########################################################################################################
##########################################################################################################
def generate_input(selevt,selevt_ID):
    """Generate the stream of filtered danton events"""
    stream = [
"    Event   PID    Energy             Direction or Momentum                       Position                    Weight",
"                    (GeV)                 (1 or GeV/c)                               (m)",
"                                ux or Px     uy or Py    uz or Pz        X            Y            Z",
    ]

    for event in selevt:  
        stream.append("      {:d} {:4d} {:12.5E} {:12.5E} {:12.5E} {:12.5E} {:12.3f} {:12.3f} {:12.3f} {:12.5E}".format(event.id,event.primary.pid,event.primary.energy,event.primary.direction[0],event.primary.direction[1],event.primary.direction[2],event.primary.position[0],event.primary.position[1],event.primary.position[2],event.weight))
        for decay in event.decay:  
            stream.append("         {:d} {:4d} {:12.5E} {:12.5E} {:12.5E} {:12.5E} {:12.3f} {:12.3f} {:12.3f}".format(decay.generation,decay.tau_i.pid,decay.tau_i.energy,decay.tau_i.direction[0],decay.tau_i.direction[1],decay.tau_i.direction[2],decay.tau_i.position[0],decay.tau_i.position[1],decay.tau_i.position[2]))
            stream.append("                {:12.5E} {:12.5E} {:12.5E} {:12.5E} {:12.3f} {:12.3f} {:12.3f}".format(decay.tau_f.energy,decay.tau_f.direction[0],decay.tau_f.direction[1],decay.tau_f.direction[2],decay.tau_f.position[0],decay.tau_f.position[1],decay.tau_f.position[2]))

            for product in decay.product:
                stream.append("           {:4d}              {:12.5E} {:12.5E} {:12.5E}".format(product.pid,product.momentum[0],product.momentum[1],product.momentum[2]))
    return "\n".join(stream)

##########################################################################################################
def parse_build(event=None,decay=None):
    """ From a DANTON event, retrieve the shower characteristics and build some statistics """

    r1 = decay.tau_i.position
    r2 = decay.tau_f.position
    R2 = np.linalg.norm(r2)
    u2 = decay.tau_f.direction
    depth = danton.EARTH_RADIUS - np.linalg.norm(r1)
    height = R2 - danton.EARTH_RADIUS 
    theta_danton = np.degrees(np.arccos(np.dot(u2, r2) / R2))
              
    et=0.0 #in GeV
    for i in decay.product:
        pp=i[1]
        ep=np.sqrt(pp[0]**2+pp[1]**2+pp[2]**2) #in GeV ignoring mass energy
        et=et+ep #in GeV
    et=et*1e9 #eV

    return depth,height,theta_danton,et    

##########################################################################################################
if __name__ == '__main__':
    main()
