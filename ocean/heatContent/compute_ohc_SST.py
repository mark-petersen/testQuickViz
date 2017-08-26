'''
Compute and plot moc, as a post-processing computation
August 2016, Mark Petersen, LANL
'''
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

########################################################################
#
#  Set directory, file, and variable names
#
########################################################################

fileNameMesh = 'restarts/mpaso.rst.0002-01-01_00000.nc'
print 'Reading data from: ' + fileNameMesh
ncfileMesh = Dataset(fileNameMesh,'r') 

nCells = ncfileMesh.dimensions['nCells'].size
nVertLevels = ncfileMesh.dimensions['nVertLevels'].size
areaCell = ncfileMesh.variables['areaCell'][:]
maxLevelCell = ncfileMesh.variables['maxLevelCell'][:]

sumArea = 0.0
for iCell in range(nCells):
   sumArea += areaCell[iCell]

########################################################################
#
#  Open files, load variables
#
########################################################################
print('** load variables: \n')
avgSST = np.zeros(12*200)
avgTemp = np.zeros(12*200)
i=0
for iYear in range(70,200):
   for iMonth in range(1,13):

	fileName = 'ocnMonthlyAverage/mpaso.hist.am.timeSeriesStats.'+'%0*d'%(4, iYear)+'-'+'%0*d'%(2,iMonth)+'-01.nc'
	#print 'Reading data from: ' + fileName
	ncfile = Dataset(fileName,'r') 
	temperature = ncfile.variables['time_avg_activeTracers_temperature'][:,:,:]
	layerThickness = ncfile.variables['time_avg_layerThickness'][:,:,:]

	sumSST = 0.0
	sumTemp = 0.0
	sumVol = 0.0
	for iCell in range(nCells):
	   sumSST += temperature[0,iCell,0] * areaCell[iCell]
	   for k in range(maxLevelCell[iCell]):
	      sumTemp += temperature[0,iCell,k] * areaCell[iCell] * layerThickness[0,iCell,k]
	      sumVol += areaCell[iCell] * layerThickness[0,iCell,k]
	avgSST[i] = sumSST/sumArea
	avgTemp[i] = sumTemp/sumVol
        #print "year=",iYear,"month=",iMonth,"SST=",avgSST[i],"T=",avgTemp[i]
        print iYear,",",iMonth,",",avgSST[i],",",avgTemp[i]
	i += 1
	ncfile.close()


print 'avgSST',avgSST
print 'avgTemp',avgTemp
exit()

########################################################################
#
#  Plots
#
########################################################################
print('** plot moc: \n')

plt.clf()
plt.figure(figsize=(8,4.5))
X, Y = np.meshgrid(mocLat,refTopDepth) # change to zMid
contourSet = plt.contourf(X,Y,mocTop.T,levels=contourLines)
contourSet2 = plt.contour(X,Y,mocTop.T,levels=contourLines[::2],colors='k',hold='on')
plt.xlabel('latitude')
plt.ylabel('depth, m')
plt.title('MPAS-Ocean Atlantic MOC, Sv')
plt.gca().invert_yaxis()
#plt.clabel(contourSet, fontsize=10)
# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(contourSet)
cbar.add_lines(contourSet2)
imageFileName = 'moc'
plt.savefig('f/'+imageFileName+'.png')

ncfileMask.close()
