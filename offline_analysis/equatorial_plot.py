from netCDF4 import Dataset
import matplotlib.pyplot as plt
import argparse
import sys
import numpy as np

def main(arguments):

  parser = argparse.ArgumentParser(description=
  '''
  Plot a variable along the equator from a given run. Written because of offline ideal age biases
  ''',
  epilog="Written by A.E.Shao 2017")

  parser.add_argument('varfile', help="Path to file with the variable to plot", type=str)
  parser.add_argument('varname', help="Variable to plot",                   type = str)
  parser.add_argument('tidx',    help="Time level of input file",           type = float, default = -1.)
  parser.add_argument('plotlat', help="Latitude to make the zonal transect",type = float, default = 0.)
  parser.add_argument('maxdepth',help="Maximum depth of transect",          type = float, default = 2000.)
  args = parser.parse_args()

  ncfile = Dataset(args.varfile)
  depth = ncfile.variables['z_l'][:]
  lon = ncfile.variables['xh'][:]
  lat = ncfile.variables['yh'][:]
  latidx = np.argmin( np.abs(lat-args.plotlat) )
  plotvar = ncfile.variables[args.varname][args.tidx,:,latidx,:]

  plt.figure()
  plt.pcolormesh(lon, depth, plotvar.squeeze())
  plt.xlim( (-240, -80) )
  plt.ylim( (0, args.maxdepth) )
  plt.gca().invert_yaxis()
  plt.colorbar()
  plt.title(args.varfile[-40:])
  plt.show()

if __name__ == "__main__":
  sys.exit(main(sys.argv[1:]))
