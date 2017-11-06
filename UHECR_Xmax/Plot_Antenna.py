#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 11:46:58 2017

@author: guepin
"""

import numpy as np
import matplotlib.pyplot as plt


filename = 'GRAND_antenna.list'
f_nu = open(filename, 'r') 
table_nu = f_nu.read()
c_nu = StringIO(unicode(table_nu))
info_nu = np.loadtxt(c_nu)

tab_X = info_nu[:,0]
tab_Y = info_nu[:,1]
tab_Z = info_nu[:,2]



origin = 'lower'
#
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.figure()
ax = plt.gca()

#ax.set_xscale('log')
#ax.set_yscale('log')

#ax.set_xlim([1e12,1e19])
#ax.set_ylim([1e-11,1e-6]) # GeV

plt.plot(tab_X, tab_Y, linestyle='', marker='o')

ax.tick_params(labelsize=16)

#plt.xlabel(r"${\rm Energy}$",fontsize=16)
#plt.ylabel(r"${\rm Cross section (mb)}$",fontsize=16)

plt.show()