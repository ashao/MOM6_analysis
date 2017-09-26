#!/usr/bin/env python
# Estimate vertical velocity (w) by the divergence of the horizontal mass transports. The vertical mass transport across an
# interface is the cumulative integral starting from the bottom of a water column. The sign convention is w>0 is an upward velocity,
# (i.e., towards the surface of the ocean). By this convention, then div(u,v)<0 implies a negative (downward) velocity and vice vers.

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

try: import argparse
except: raise Exception('This version of python is not new enough. python 2.7 or newer is required.')

parser = argparse.ArgumentParser(description=
         '''Script for calculating vertical velocity from divergence of horizontal mass transports.''')
parser.add_argument('infile', type=str, help='''Daily file containing ssu,ssv.''')
parser.add_argument('--uname', type=str, default='umo', help='''Name of the u-component of mass transport''')
parser.add_argument('--vname', type=str, default='vmo', help='''Name of the v-component of mass transport''')
parser.add_argument('--wrapx', type=bool, default=True, help='''True if the x-component is reentrant''')
parser.add_argument('--wrapy', type=bool, default=False, help='''True if the x-component is reentrant''')
parser.add_argument('--gridspec', type=str,
  help='''File containing variables wet,areacello. Usually the ocean_static.nc from diag_table.''')
parser.add_argument('--plot', type=bool, default=False, help='''Plot w at the bottom of the first layer''')
cmdLineArgs = parser.parse_args()

# Check for presence of variables
rootGroup = netCDF4.Dataset( cmdLineArgs.infile )
if cmdLineArgs.uname not in rootGroup.variables:
  raise Exception('Could not find "%s" in file "%s"' % (cmdLineArgs.uname, cmdLineArgs.infile))
if cmdLineArgs.vname not in rootGroup.variables:
  raise Exception('Could not find "%s" in file "%s"' % (cmdLineArgs.vname, cmdLineArgs.infile))

# Get grid-related metrics
msk = netCDF4.Dataset(cmdLineArgs.gridspec).variables['wet'][:,:]
area = msk*netCDF4.Dataset(cmdLineArgs.gridspec).variables['areacello'][:,:]
# Get mass transports
# Add a slice for the halo in case there's re-entrance
ntime, nk, nlat, nlon = netCDF4.Dataset(cmdLineArgs.infile).variables[cmdLineArgs.uname].shape
u_var = np.zeros( (ntime, nk, nlat, nlon+1) )
u_var[:,:,:,1:] = msk*netCDF4.Dataset(cmdLineArgs.infile).variables[cmdLineArgs.uname][:,:,:,:]
if cmdLineArgs.wrapx:
  u_var[:,:,:,0] = u_var[:,:,:,-1]
ntime, nk, nlat, nlon = netCDF4.Dataset(cmdLineArgs.infile).variables[cmdLineArgs.vname].shape
v_var = np.zeros( (ntime, nk, nlat+1, nlon) )
v_var[:,:,1:,:] = msk*netCDF4.Dataset(cmdLineArgs.infile).variables[cmdLineArgs.vname][:,:,:,:]
if cmdLineArgs.wrapy:
  v_var[:,:,0,:] = v_var[:,:,-1,:]
# Dimensions on h-points
nlat, nlon = area.shape
w = np.zeros( (ntime, nk+1, nlat, nlon) )
w[:,1:,:,:] += u_var[:,:,:,1:] - u_var[:,:,:,:-1]
w[:,1:,:,:] += v_var[:,:,1:,:] - v_var[:,:,:-1,:]
w = w[:,::-1,:,:].cumsum(axis=1)[:,::-1,:,:]
# Convert from volume flux to m/s
w = w/area

if cmdLineArgs:
  plt.pcolormesh(w[0,1,:,:])
  plt.colorbar()
  plt.show()
