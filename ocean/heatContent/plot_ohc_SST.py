'''
Compute and plot moc, as a post-processing computation
August 2016, Mark Petersen, LANL
'''
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

def movingaverage (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma

########################################################################
#
#  Open files, load variables
#
########################################################################
print('** load variables: \n')

data = np.genfromtxt('output', delimiter=',')

timeInYears = data[:,0]-1. + data[:,1]/12. - 1/24. 
timeInYears2 = data[:,0]-1. + data[:,1]/12. - 1/24. + 1948-4*60
timeInYears3 = data[:,0]-1. + data[:,1]/12. - 1/24. + 1948-60
sst = data[:,2]
temp = data[:,3]

tempSmooth = movingaverage (temp,12)
sstSmooth = movingaverage (sst,12)

# Specific heat capacity, Cp=3985J/kg/K
# density of 1027 kg/m3
# volume of ocean: 1,347,000,000 cu km
ohc = (temp+273.15)*1027*3985 * 1.347e9*1.e9/1.e23
print temp
print temp+273.15
print ohc
# units: K * kg/m3 * J/kg/K *m3 = J

########################################################################
#
#  Plots
#
########################################################################
print('** plot moc: \n')

plt.clf()
plt.figure(figsize=(8,4.5))
plt.plot(timeInYears,sst)
plt.xlabel('time,years')
plt.ylabel('averge SST, C')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'sst'
plt.savefig('f/'+imageFileName+'.png')

plt.clf()
plt.figure(figsize=(8,4.5))
plt.plot(timeInYears,temp)
plt.xlabel('time,years')
plt.ylabel('Average global temperature, C')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'temperature'
plt.savefig('f/'+imageFileName+'.png')


plt.clf()
plt.figure(figsize=(8,4.5))
plt.plot(timeInYears,ohc)
plt.xlabel('time,years')
plt.ylabel('ocean heat content, J x1e23')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'ohc'
plt.savefig('f/'+imageFileName+'.png')

## convert J/time into W/m^2:
#2e23J/15yrs/5.1e8km2

2e23/15/1e6/5.1e8/86400/365

2.2e23/20/1e6/5.1e8/86400/365
1e23/20/1e6/5.1e8/86400/365
1e23/80/1e6/5.1e8/86400/365

plt.clf()
plt.figure(figsize=(8,4.5))
plt.plot(timeInYears[6:np.size(tempSmooth)+6],tempSmooth)
#plt.plot(timeInYears,temp)
plt.ylim([3.3,4.1])
plt.xlabel('time,years')
plt.ylabel('temperature, C, one year smoothing')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'temp2'
plt.savefig('f/'+imageFileName+'.png')

plt.clf()
plt.figure(figsize=(8,4.5))
#plt.plot(timeInYears,sst)
plt.plot(timeInYears3[6:np.size(sstSmooth)+6],sstSmooth)
plt.ylim([17.78,18.6])
plt.xlim([1948,2008])
plt.xlabel('time,years')
plt.ylabel('global average SST, C, one year smoothing')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'sst2'
plt.savefig('f/'+imageFileName+'.png')


plt.clf()
plt.figure(figsize=(8,4.5))
plt.plot(timeInYears2[6:np.size(tempSmooth)+6],tempSmooth)
#plt.plot(timeInYears,temp)
plt.ylim([3.3,4.1])
plt.xlabel('time,years')
plt.ylabel('temperature, C, one year smoothing')
plt.title('ACME EC60to30v2, CORE-II G case')
plt.grid()
imageFileName = 'temp3'
plt.savefig('f/'+imageFileName+'.png')
