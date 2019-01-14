import numpy as np
import math
import matplotlib.pyplot as plt
import pyfits as fits
import scipy.integrate as integrate
from scipy.special import kn
import time
from simpson import * # in ~/mypythonlib

image_file="gyoto_spectrum.fits"

image_cube = fits.getdata(image_file)

# Lower and upper values of the exponent of your frequency
# as chosen in the XML Spectrometer field.
low_freq_exp=11.
high_freq_exp=18.

convertSItoCGS=1e3 # SI to cgs conversion for I_nu (gyoto speaks in SI)

fov=500 # observer's field of view in µas, as specified
# in the XML FieldOfView field

# Sgr A* distance and mass
DD_kpc=8.1
MM_sunmass=4.1e6

# Choose whether to plot or not Sgr A* observed data
plotobsradio=1
plotobsIR=1
plotobsX=1

###################################################
############ End of input parameters ##############
############  Nothing to edit below  ##############
###################################################

kiloparsec_cgs  = 1e3*3.08e18
sunmass_cgs = 1.99e33
DD=DD_kpc*kiloparsec_cgs
MM=MM_sunmass*sunmass_cgs # M used only if unit_length=0

GG=6.67e-8;c2=2.99792458e10*2.99792458e10
GMoc2 = GG*MM/c2

angle=0.5*fov*1e-6/3600.*math.pi/180.

nbnuobs = image_cube.shape[0] # number of frequencies computed
npix = image_cube.shape[1] # number of pixels along one dimension of Screen

dfovreal=2*angle/npix # number of rad per pxl

spectrum = np.zeros(nbnuobs, dtype=float)

# Put in spectrum, for each frequency, the sum of all intensities
# over screen pixel, multiplied by the solid angle subtended by one pixel
for ii in range(0,nbnuobs):
    imnu = image_cube[ii,:,:]
    spectrum[ii] = mysimps2D(imnu,dfovreal,dfovreal)
    
spectrum*=convertSItoCGS

# Create the array containing the values of frequencies computed, in Hz
# Remember the frequencies are evenly separated in log space
nuobs = np.zeros(nbnuobs, dtype=float)
freq_exp_array = np.zeros(nbnuobs+1, dtype=float)
for ii in range(0,nbnuobs+1):
    freq_exp_array[ii] = low_freq_exp + ii*(high_freq_exp-low_freq_exp)/nbnuobs
for ii in range(1,nbnuobs+1):
    nuobs[ii-1] = 0.5*(pow(10,freq_exp_array[ii-1]) + pow(10,freq_exp_array[ii]))

plt.figure(0)
plt.clf()
plt.ion()

plt.loglog(nuobs,4*math.pi*DD*DD*nuobs*spectrum,color="red")
plt.xlabel('$\\nu$ (Hz)',fontsize=14)
plt.ylabel('$\\nu L_\\nu$ (erg/s)',fontsize=14)
plt.show()

# SgrA* observation

if plotobsradio==1:
    # RADIO
    # Data from Bower+15 Table 7
    # To these points, are added (last 4 points of the vector):
    # 2 points at 100GHz from Brinkerink+15 Fig. 2
    # 1 point at 492 GHz from Liu+16 (A&A 593 A44) see abstract
    # 1 point at 690GHz from Marrone+06
    nuobsdata=np.asarray([1.6,3.1,5.4,9.,14.,21.1,32.,40.9,218.,220.,231.9,233.8,341.6,343.6,351.7,353.6,216.8,223.9,238.2,266.8,274.,331.1,338.3,352.6,98.,108.,492.,690.])*1e9
    fluxobs=np.asarray([0.592,0.702,0.87,0.932,1.075,1.164,1.382,1.485,3.667,3.661,3.676,3.704,3.602,3.609,3.595,3.553,3.677,3.391,3.310,3.369,3.526,3.205,3.436,4.89,2.41,2.6,3.6,3.8])*1e-23 # Jy
    errobs=np.asarray([0.028,0.032,0.118,0.129,0.135,0.052,0.087,0.073,0.65,0.652,0.664,0.68,0.866,0.87,0.884,0.86,0.762,0.489,0.424,0.096,0.697,1.074,0.863,0.721,0.18,0.2,0.72,2.2])*1e-23
    plt.scatter(nuobsdata,4*math.pi*DD*DD*nuobsdata*fluxobs,color="black")
    plt.errorbar(nuobsdata,4*math.pi*DD*DD*nuobsdata*fluxobs,yerr=4*math.pi*DD*DD*nuobsdata*errobs,linestyle='None',ecolor="black",elinewidth=3)

if plotobsIR==1:
    # FIR point from Fellenberg+18 (upper lims)
    ll=np.asarray([160.,100.])*1e-6 # lambda in m
    c_SI=299792458.
    nuobsdata=c_SI/ll # nu in Hz
    fluxFIR=np.asarray([1.06,0.64])*1e-23 # Jy
    errFIR=np.asarray([0.24,0.4])*1e-23
    plt.scatter(nuobsdata,4*math.pi*DD*DD*nuobsdata*fluxFIR,color="black")
    plt.errorbar(nuobsdata,4*math.pi*DD*DD*nuobsdata*fluxFIR,yerr=4*math.pi*DD*DD*nuobsdata*errFIR,uplims=4*math.pi*DD*DD*nuobsdata*errFIR,linestyle='None',ecolor="black",elinewidth=3)

    # NIR points taken from Witzel+18
    ll=np.asarray([4.5,2.18])*1e-6 # lambda in m
    nuobsdata_IR=c_SI/ll # nu in Hz
    nnuLnu_IR_8d3kpc=np.asarray([3.2,2.6])*1e34
    nnuLnu_IR_DD=nnuLnu_IR_8d3kpc/8.3**2 * DD_kpc**2
    fluxobs_IR=nnuLnu_IR_DD/(nuobsdata_IR*4.*math.pi*DD**2) #cgs flux
    errobs_IR_8d3kpc=np.asarray([1.4,1.2])*1e34 # in nu*Lnu for 8.3kpc
    errobs_IR_DD=errobs_IR_8d3kpc/8.3**2 * DD_kpc**2
    errobs_IR=errobs_IR_DD/(nuobsdata_IR*4.*math.pi*DD**2) # cgs flux error
    plt.scatter(nuobsdata_IR,4*math.pi*DD*DD*nuobsdata_IR*fluxobs_IR,color="black")
    plt.errorbar(nuobsdata_IR,4*math.pi*DD*DD*nuobsdata_IR*fluxobs_IR,yerr=4*math.pi*DD*DD*nuobsdata_IR*errobs_IR,linestyle='None',ecolor="black",elinewidth=3)
    
if plotobsX==1:
    # From Baganoff+01: 2-10 keV integrated luminosity is 2.2[+0.4 -0.3]e33 erg/s
    numin=2. * 1e3*1.6e-19/6.62e-34 # in Hz
    numax=10. * 1e3*1.6e-19/6.62e-34
    luminosityX=2.2e33
    errlum_low=0.3e33
    errlum_high=0.4e33
    Gamma=2.2
    index=2.-Gamma
    integral=(numax**index/index - numin**index/index)
    Aavg=luminosityX/integral
    Ahigh=(luminosityX+errlum_high)/integral
    Alow=(luminosityX-errlum_low)/integral
    plt.plot([numin,numax],[Ahigh*numin**index,Alow*numax**index],color="black",linewidth=3)
    plt.plot([numin,numax],[Alow*numin**index,Ahigh*numax**index],color="black",linewidth=3)
    plt.plot([numin,numin],[Alow*numin**index,Ahigh*numin**index],color="black",linewidth=3)
    plt.plot([numax,numax],[Alow*numax**index,Ahigh*numax**index],color="black",linewidth=3)

