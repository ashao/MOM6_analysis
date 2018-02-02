#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse
import numpy as np
from netCDF4 import Dataset
#import matplotlib.pyplot as plt

def main(arguments):
  args = parse_input_arguments(arguments)

  if (args.tmax==-1):
    args.tmax = np.size(Dataset(args.infile).variables['time'][:])

#  plt.figure()
  for tidx in range(args.tmin,args.tmax):
    lhs = Dataset(args.infile).variables['dhdt'][tidx,:,:,:]
    temp_rhs = calculate_rhs(args.infile,tidx)
    print("---Cell-by-cell statistics---")
    max_idx = np.unravel_index(np.argmax(np.abs(lhs-temp_rhs)),lhs.shape)
    print("Maximum difference: %e of %e at %s" % (np.max(np.abs(lhs-temp_rhs)), lhs[max_idx], max_idx))
    print("Mean difference: %f" % np.mean(np.abs(lhs-temp_rhs)))
    print("Sum difference: %f" % np.sum(np.abs(lhs-temp_rhs)))
    print("Global sum of total tendency: %f Error: %f" % (np.abs(lhs).sum(), np.sum(np.abs(lhs-temp_rhs))))
    print("---Column integral statistics---")
    col_err = np.sum(np.abs( lhs - temp_rhs),axis=0)
    max_idx = np.unravel_index(np.argmax(col_err),lhs.shape[1:])
    print("Maximum difference: %e at %s"  % (col_err.max(),max_idx))
    print("Mean difference: %e" % np.mean( np.abs(lhs.sum(axis=0)-temp_rhs.sum(axis=0)) ) )
    print("Sum difference: %e" % np.sum( np.abs(lhs.sum(axis=0)-temp_rhs.sum(axis=0)) ) )

#    plt.pcolormesh(np.sum(np.abs(lhs-temp_rhs),axis=0),cmap=plt.cm.coolwarm)
#    plt.colorbar()
#    plt.show()



def parse_input_arguments(arguments):
  parser = argparse.ArgumentParser(description=
  '''
  Calculates the difference between the total tendency of temp content cell by cell and returns the globally integrated error summed
  over all time levels in the file unless otherwise specified.
  Assumes the following variables exist in the input file:
  boundary_forcing_temp_tendency
  S_advection_xy
  Sh_tendency_vert_remap
  otempdiff
  otemppmdiff
  Sh_tendency
  ''',
  epilog="Written by A.E.Shao 2017")

  parser.add_argument('infile', help="Input file", type=str)
  parser.add_argument('--tmin', help="First time level in calculation", type=int, default=0)
  parser.add_argument('--tmax', help="Last time level in calculation", type=int, default=-1)
  args = parser.parse_args(arguments)
  return args

# Right hand side of equation is the sum of all individual tendencies
def calculate_rhs(infile,tidx):
  rhs =  Dataset(infile).variables['vert_remap_h_tendency'][tidx,:,:,:]
  rhs += Dataset(infile).variables['dynamics_h_tendency'][tidx,:,:,:]
  rhs += Dataset(infile).variables['boundary_forcing_h_tendency'][tidx,:,:,:]
  return rhs

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
