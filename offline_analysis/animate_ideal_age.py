#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.basemap import Basemap

agefile = '/archive/aes/projects/offline_online_comparison/offline_vgrid/5day/age.0075.0100.nc'
gridfile = '/archive/aes/projects/offline_online_comparison/ocean_annual.static.nc'

geolat = Dataset(gridfile).variables['geolat'][:,:]
geolon = Dataset(gridfile).variables['geolon'][:,:]
ntime = np.size(Dataset(agefile).variables['time'][:])

def plot_spatial_map(data,fig):
  m.pcolormesh(x,y,data,vmin=0,vmax=60,cmap=plt.get_cmap('nipy_spectral'))
  fig.canvas.draw_idle()
  plt.pause(0.1)
  fig.show()

def main(arguments):
  t = 0
  fig = plt.figure(figsize=(12,9))
  
  age = Dataset(agefile).variables['age'][:,9,:,:].squeeze()
  while t < ntime:
    plt.clf()
    m = Basemap(projection='eck4',lon_0=(-300+60)/2,resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    x,y = m(geolon,geolat)
    m.pcolormesh(x,y,age[t,:,:],vmin=0,vmax=60,cmap=plt.get_cmap('nipy_spectral'))
    print("Making timelevel %d/%d" % (t,ntime-1))
    fig.savefig('age_%04d.png' % t)
    t = t + 1

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))

