import numpy as np
from optparse import OptionParser
from Tkinter import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import cPickle
import scipy.interpolate as intp
import scipy.optimize as opt
from scipy.special import gamma
import random
import os
import re
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.optimize import curve_fit
  
def model(x, a, b, c):# model(x, a=1.0, b=700., c=1.65e9)
  return a*(x-b)**2 +c #parabola, b=sx, c=sy
 
label_size = 16
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
 
# ATTENTION Check for all the hardcoded paths before using the scripts! eg search for CR185


nfiles= 68 #len(txt1)  
#energy=3.1e17
#zenith=36.87 # 36.87  45.57#0.0 #
#nrp=50
#nrFe=20
#nrp=41
#nrFe=17
antnr=1 # =secondline#4 #=fourthline
path='./ECR185_1/'# HERE: set the folder of the simulations you wanna study
datadir='./'






 
 ### TODO: need list with all real Xmax from simulations 
  
#### loop ueber l  , Xreco[l]  
Xreal=np.zeros([nfiles])
Xreco=np.zeros([nfiles])
Diff=np.zeros([nfiles])
for l in np.arange(int(nfiles)):



      
   longfile = datadir+'CR-Sim/CR185_{0}/CR185_{0}.sry'.format(l) #{0}/DAT{1}-999999999.long when from parallel
   #print longfile

   longprofile=0.
   Xground=0.
   hillas=0.  


   with open(longfile) as f:
        for line in f:
            keyword = ' Sl. depth of max. (g/cm2):' 
            if line.startswith(keyword):
                    #try:
                        hillas =line.split('  ')[1]
   #print "xmax=", hillas, "g/cm2"
                    #except ValueError:
                    ## treat the error
                        #print "error"
 
 
   Xreal[l]= hillas
#print hillas
 

  ##################### 
#GRAND_xmaxcurve_CR185_1_ev32_realXmax822.051_CtCR185_32.txt 
   

 #read in: primary, hillas, combchi2 
   #name = path+ 'GRAND_xmaxcurve_{0:06d}_ev{1}_realXmax{2:4.2f}.txt'.format(int(fileno1), l, Xreal[l])
   #txt=np.loadtxt(name)
   try:
    name = path+ 'GRAND_xmaxcurve_CR185_1_ev{0}_realXmax{1:4.3f}_CtCR185_{0}.txt'.format( l, Xreal[l])
    #name = path+ 'GRAND_xmaxcurve_{0:06d}_ev{1}_realXmax{2:4.2f}_Ct370141.txt'.format(int(fileno1), l, Xreal[l])
    txt=np.loadtxt(name)
   except IOError:
     try: 
      name = path+ 'GRAND_xmaxcurve_CR185_1_ev{0}_realXmax{1:3.3f}_CtCR185_{0}.txt'.format( l, Xreal[l])
      #name = path+ 'GRAND_xmaxcurve_{0:06d}_ev{1}_realXmax{2:4.1f}_Ct370141.txt'.format(int(fileno1), l, Xreal[l])

      txt=np.loadtxt(name)
     except:
        name = path+ 'GRAND_xmaxcurve_CR185_1_ev{0}_realXmax{1:3.2f}_CtCR185_{0}.txt'.format( l, Xreal[l])
        txt=np.loadtxt(name)
        #print "error"

#print name

    
   #print txt
   primary  = txt.T[0]
   hillas = txt.T[1]
   combchi2  = txt.T[2]
   nsim = len(txt)
   #print nsim
   
   combchi2[l]=100000 # punish the simulated event so that is does not go into fitting

   #print "primary", primary, "xmax sim ", hillas, "Chi2", combchi2

### for plotting

#proton   
   h_p=[]
   c_p=[]
   for i in np.arange(nsim):
     if (primary[i] <20):
      h_p.append(hillas[i])# for plotting
      c_p.append(combchi2[i])# for plotting

   #plt.scatter(h_p,c_p[:],50,color='r',  label="Proton")
   #print hillas
   #print combchi2[:]
   #plt.scatter(hillas[fit_selection>0],combchi2[fit_selection>0],50,color='m')
   hfit_p=[]
   cfit_p=[]
   for i in np.arange(nsim):
    if (primary[i] <20):
	hfit_p.append(hillas[i])#for fitting
	cfit_p.append(combchi2[i])#for fitting
   #plt.scatter(hfit_p,cfit_p,50,color='r')
   
#Iron
   h_Fe=[]
   c_Fe=[]
   for i in np.arange(nsim):
     if (primary[i] >20):
      h_Fe.append(hillas[i]) # for plotting
      c_Fe.append(combchi2[i])# for plotting

   #plt.scatter(h_Fe,c_Fe[:],50,color='b', marker='s', label="Iron")
   #print hillas
   #print combchi2[:]
   #plt.scatter(hillas[fit_selection>0],combchi2[fit_selection>0],50,color='m')
   hfit_Fe=[]
   cfit_Fe=[]
   for i in np.arange(nsim):
    if (primary[i] >20):
	hfit_Fe.append(hillas[i]) #for fitting
	cfit_Fe.append(combchi2[i])#for fitting
   #plt.scatter(hfit_p,cfit_p,50,color='b')
   

   
##all for fitting: here put some condituons to chose for ehich events should be included besides chi2 smaller than 5
   #h=[]
   #c=[]
   #for i in np.arange(nsim):
      #h.append(hillas[i])
      #c.append(combchi2[i])
   ##plt.scatter(h,c[:],50,color='r')
   ##print hillas
   ##print combchi2[:]
   ##plt.scatter(hillas[fit_selection>0],combchi2[fit_selection>0],50,color='m')
   
   bestsim= np.argmin(combchi2)
   urange=hillas[bestsim]+100 # parallel
   drange=hillas[bestsim]-100 # parallel
   chirange=combchi2[bestsim]+20
   
   #print "ranges", urange, drange, chirange, "bestim", hillas[bestsim], combchi2[bestsim], bestsim 
   
   fit_selection=np.zeros(nsim) # all points that have lower chi2 values on on side only
   
   for i in np.arange(nsim):
      if (np.sum(combchi2[(hillas[:]>hillas[i])]<combchi2[i])==0): fit_selection[i]=fit_selection[i]+1 # parallel
      if (np.sum(combchi2[(hillas[:]<hillas[i])]<combchi2[i])==0): fit_selection[i]=fit_selection[i]+1 # parallel
     
   #print 'combchi2...'
   #print combchi2#, hillas, np.sum(combchi2), fit_selection
   #print 'hillas...'
   #print hillas[:]>drange, hillas[:]<urange, combchi2<chirange

   
   #fit_selection=fit_selection*(hillas[2,:]>drange)*(hillas[2,:]<urange)*(combchi2<chirange) # non parallel
   fit_selection=fit_selection*(hillas[:]>drange)*(hillas[:]<urange)*(combchi2<chirange) # parallel
   #print 'fit_selection ', fit_selection

   #if (simmode): fit_selection[simevent]=0 # for fake data
   
   
   #hillas[(fit_selection>0)],combchi2[(fit_selection>0) # fit_selection=liste
   
#print "selected", hillas[(fit_selection>0)],combchi2[(fit_selection>0)], len(hillas[(fit_selection>0)])
   
   hfit=hillas[(fit_selection>0)]
   cfit= combchi2[(fit_selection>0)]
   
   #hfit=[]
   #cfit=[]
   #for i in np.arange(nsim):
    ##if(fit_selection[i]>0):

    #if(combchi2[i]<5): #20):## NOTE: Change this to a more reasonable number like 4
    ##if( abs(Xreal[i]-hillas[i])<=200):  
	#hfit.append(hillas[i])
	#cfit.append(combchi2[i])
   ##plt.scatter(hfit,cfit,50,color='r')
   
   #####################
 
 
   try:
      popt, pcov = curve_fit(model, hfit, cfit, maxfev=20000)
  
      Xreco[l]=popt[1]
   except: # RuntimeError as exc: #except
      print "exc" 
      ind = np.where(combchi2 == combchi2.min()) 
      print ind, len(ind)


      try: 
	Xreco[l]=hillas[ind]
      except ValueError as exc:
	Xreco[l]=10000
 
  #####################   
   #if (l==l):#29):  # choose your favorite simulation
   fig1=plt.figure(4,figsize=(8,8))
   ##axins = fig.add_axes([0.57, 0.15, 0.3, 0.3])
   	
	#plt.xlim(Xreco[l]-150,Xreco[l]+150)
   	
   plt.ylim(0.8,10.0)
   t = np.linspace(min(hfit), max(hfit), 50)
   try:
       plt.plot(t, model(t, *popt), label="Fitted Curve")
   except:
       continue
	
   plt.ylabel("Chi2/ndf", fontsize=20)
   plt.xlabel("Xsim ( g/cm$^2$)", fontsize=20)

   plt.scatter(h_Fe,c_Fe,50,color='b', marker='s')### written in File as: c_Fe[:]/(ndf_comb+1e-25)
   plt.scatter(h_p,c_p,50,color='r', marker='o')

   name1= path+ 'Fit_xmaxcurve_ev{0}_realXmax{1:4.2f}.svg'.format( l, Xreal[l])
   fig1.savefig(name1)  
   plt.xlim(Xreco[l]-150,Xreco[l]+150)
   plt.ylim(0.,5.0)
   name1= path+ 'Fit_xmaxcurve_ev{0}_realXmax{1:4.2f}_zoom.svg'.format(l, Xreal[l])
   fig1.savefig(name1)
   
   plt.clf()
   del fig1
   del hfit, cfit, h_Fe, h_p, c_Fe, c_p

   print "eventnr", l,  "Xreco", Xreco[l], "Xreal", Xreal[l], "bestim", hillas[bestsim]
   Diff[l]=abs(Xreco[l]- Xreal[l])
   print "diff", abs(Xreco[l]- Xreal[l])
   
   
a = Diff[:]
a.sort() 

total=nfiles

proz=0
position=0
for m in np.arange(int(nfiles)):  
  if(proz < 0.68*total):
    #proz+=a[m]
    proz+=1
  if(proz >= 0.68*total):
      #print a[m], m
    position = a[m]
    break

fig=plt.figure(1,figsize=(8,6))
plt.hist(Diff, bins=int(max(Diff)/2), histtype='step')
plt.axvline(position,color='k', linestyle='--')

plt.xlim(0,60)
plt.ylabel("Nr. of Sim.", fontsize=20)
plt.xlabel("Abs(Xreco-Xreal) ( g/cm$^2$)", fontsize=20)
plt.show()

name =path+'Fittest_CR185_unc{0}.svg'.format(position)
fig.savefig(name, format='svg', dpi=1000)
