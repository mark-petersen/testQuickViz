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

fileName = '/lustre/scratch3/turquoise/maltrud/ACME/cases/GMPAS_oEC60to30_test03/run/am.mpas-o.timeSeriesStats.0002-02-01.nc'
fileNameMask = '/lustre/scratch1/turquoise/mpeterse/ACME/cases/a09b/oEC60to30.151031.SingleRegionAtlantic.nc'
contourLines = np.arange(-8,24.1,2)

timeIndex = 0
horVelName = 'time_avg_normalVelocity'
vertVelName = 'time_avg_vertVelocityTop'

# This is specific to Atlantic MOC:
mocLat = np.arange(-34.5,70.1,0.5)
#contourLines = np.arange(-6,18.1,2)

########################################################################
#
#  Open files, load variables
#
########################################################################
print('** load variables: \n')

print 'Reading data from: ' + fileName
ncfile = Dataset(fileName,'r') 
nCells = ncfile.dimensions['nCells'].size
nVertLevels = ncfile.dimensions['nVertLevels'].size

horVel = ncfile.variables[horVelName][:,:,:]
vertVel = ncfile.variables[vertVelName][:,:,:]

# delete:
#print np.shape(vertVel)
#print np.size(vertVel)

print 'Reading data from: ' + fileNameMesh
ncfileMesh = Dataset(fileNameMesh,'r') 
dvEdge = ncfileMesh.variables['dvEdge'][:]
areaCell = ncfileMesh.variables['areaCell'][:]
latCell = ncfileMesh.variables['latCell'][:]*180./np.pi
refBottomDepth = ncfileMesh.variables['refBottomDepth'][:]

refLayerThickness = np.zeros(nVertLevels)
refLayerThickness[0] = refBottomDepth[0]
refLayerThickness[1:nVertLevels] = refBottomDepth[1:nVertLevels]-refBottomDepth[0:nVertLevels-1]

refTopDepth = np.zeros(nVertLevels+1)
refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]

print 'Reading transects from: ' + fileNameMask
ncfileMask = Dataset(fileNameMask,'r') 
nTransects = ncfileMask.dimensions['nTransects'].size
maxEdgesInTransect = ncfileMask.dimensions['maxEdgesInTransect'].size

transectEdgeMaskSigns = ncfileMask.variables['transectEdgeMaskSigns'][:,:]
transectEdgeGlobalIDs = ncfileMask.variables['transectEdgeGlobalIDs'][:,:]

print 'Reading regionss from: ' + fileNameMask
ncfileMask = Dataset(fileNameMask,'r') 
nRegions = ncfileMask.dimensions['nRegions'].size
regionCellMasks = ncfileMask.variables['regionCellMasks']

########################################################################
#
#  Compute transport through transects
#
########################################################################
print('** compute transport: \n')

m3ps_to_Sv = 1e-6; # m^3/sec flux to Sverdrups

# the volume transport
transport = np.zeros(nTransects);
transportZ = np.zeros([nVertLevels,nTransects]);
transportZEdge = np.zeros([nVertLevels,maxEdgesInTransect,nTransects]);

for iTransect in range(nTransects):
   for i in range(maxEdgesInTransect):
       if transectEdgeGlobalIDs[iTransect, i]==0: 
           break
       iEdge = transectEdgeGlobalIDs[iTransect, i] - 1 # subtract 1 because of python 0-indexing
       for k in range(nVertLevels):
           transportZEdge[k,i,iTransect] = ( transectEdgeMaskSigns[iEdge,iTransect] 
              * horVel[timeIndex,iEdge,k] * dvEdge[iEdge] * refLayerThickness[k]*m3ps_to_Sv )
           transportZ[k,iTransect] = transportZ[k,iTransect] + transportZEdge[k,i,iTransect]
           transport[iTransect] = transport[iTransect] + transportZEdge[k,i,iTransect]

########################################################################
#
#  Compute MOC
#
########################################################################
print('** compute moc: \n')

nLat = np.size(mocLat)
mocTop = np.zeros((nLat,nVertLevels+1))

for iRegion in range(nRegions):
   # assume transects and regions have the same index ordering:
   iTransect = iRegion

   for k in range(1,nVertLevels+1):
      mocTop[0,k] = mocTop[0,k-1] + transportZ[k-1,iTransect]/m3ps_to_Sv

   for iLat in range(1,nLat):
      ind =np.logical_and(np.logical_and(regionCellMasks[:,iRegion]==1,latCell >= mocLat[iLat-1] ), latCell <  mocLat[iLat])

      for k in range(1,nVertLevels+1):
         mocTop[iLat,k] = mocTop[iLat-1,k] + np.sum(vertVel[timeIndex,ind,k]*areaCell[ind])

   # convert m^3/s to Sverdrup
   mocTop = mocTop * m3ps_to_Sv

########################################################################
#
#  Plot MOC
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

ncfile.close()
ncfileMask.close()
