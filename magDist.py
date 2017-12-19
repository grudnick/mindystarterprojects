# Program to read in and plot magnitude distribution

from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pylab as plt 
import matplotlib.lines as mlines
# from matplotlib.legend import Legend

# Read in data
hdulist = fits.open('/home/mtownsend/anaconda3/survey-dr3-DR12Q.fits') # survey
hdulist2 = fits.open('/home/mtownsend/anaconda3/specObj-dr13.fits') # sdss
tbdata = hdulist[1].data
tbdata2 = hdulist2[1].data

# Put data in arrays
# Object ID from survey file; value -1 for non-matches
objid = []
objid = tbdata.field('OBJID') 

# Number of observations of source from legacy file
obs = []
obsmatch = []
obs = tbdata.field('DECAM_NOBS') 
obsmatch = obs[np.where(objid > -1)]

# Put number of observations per filter into arrays that match the filter
uobs = []
gobs = []
robs = []
iobs = []
zobs = []
yobs = []

b = np.array(obsmatch)
uobs = b[:,0]
gobs = b[:,1]
robs = b[:,2]
iobs = b[:,3]
zobs = b[:,4]
yobs = b[:,5]

# Put flux data in an array from legacy file
# Flux has ugrizY, so needs to be divided into 6 arrays
flux =[]
fluxmatch = []
flux = tbdata.field('DECAM_FLUX')
fluxmatch = flux[np.where(objid > -1)]

# Divide flux arrays into 6 arrays
uflux = []
gflux = []
rflux = []
iflux = []
zflux = []
yflux = []

a = np.array(fluxmatch)
uflux = a[:,0]
gflux = a[:,1]
rflux = a[:,2]
iflux = a[:,3]
zflux = a[:,4]
yflux = a[:,5]

# Put inverse flux variance data in an array from legacy file
flux_ivar =[]
flux_ivar_match = []
flux_ivar = tbdata.field('DECAM_FLUX_IVAR')
flux_ivar_match = flux[np.where(objid > -1)]

# Divide flux arrays into 6 arrays
uflux_ivar = []
gflux_ivar = []
rflux_ivar = []
iflux_ivar = []
zflux_ivar = []
yflux_ivar = []

c = np.array(flux_ivar_match)
uflux_ivar = c[:,0]
gflux_ivar = c[:,1]
rflux_ivar = c[:,2]
iflux_ivar = c[:,3]
zflux_ivar = c[:,4]
yflux_ivar = c[:,5]

# Match flux array index to nobs = 3
gflux_obs = []
rflux_obs = []
zflux_obs = []
gflux_obs = gflux[np.where(gobs == 3)]
rflux_obs = rflux[np.where(robs == 3)]
zflux_obs = zflux[np.where(zobs == 3)]

# Calculate magnitudes
gmag = 22.5 - 2.5 * np.log10(gflux_obs[np.where(gflux_obs > 0.)])
rmag = 22.5 - 2.5 * np.log10(rflux_obs[np.where(rflux_obs > 0.)])
zmag = 22.5 - 2.5 * np.log10(zflux_obs[np.where(zflux_obs > 0.)])

# Make sure there are good exposures
# gexp = []
# rexp = []
# zexp = []
# 
# gexp = gflux[np.where(gflux > 0.)] / (1./(gflux_ivar[np.where(gflux > 0.)])**0.5)
# rexp = rflux[np.where(rflux > 0.)] / (1./(rflux_ivar[np.where(rflux > 0.)])**0.5)
# zexp = zflux[np.where(zflux > 0.)] / (1./(zflux_ivar[np.where(zflux > 0.)])**0.5)

# plt.xlim((10.0, 30.0)) 

# plt.hist(zmag, bins = 50, color = 'blue')
# plt.hist(rmag[np.where(rexp > 5.)], bins = 150, color = 'red')
# plt.hist(gmag[np.where(gexp > 5.)], bins = 50, color = 'green')

# plt.hist(gmag, bins = 50, color = 'green')
plt.hist(rmag, bins = 50, color = 'red')
# plt.hist(zmag, bins = 50, color = 'blue')

plt.grid(True)
plt.title('Magnitude Distribution (g)')
plt.xlabel(r'$magnitude$')
plt.ylabel(r'$counts$')

plt.show()

