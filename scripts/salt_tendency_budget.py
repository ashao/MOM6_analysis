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
    lhs = Dataset(args.infile).variables['osalttend'][tidx,:,:,:]
    salt_rhs = calculate_rhs(args.infile,tidx)
    print("---Cell-by-cell statistics---")
    max_idx = np.unravel_index(np.argmax(np.abs(lhs-salt_rhs)),lhs.shape)
    print("Maximum difference: %e of %e at %s" % (np.max(np.abs(lhs-salt_rhs)), lhs[max_idx], max_idx))
    print("Mean difference: %f" % np.mean(np.abs(lhs-salt_rhs)))
    print("Sum difference: %f" % np.sum(np.abs(lhs-salt_rhs)))
    print("Global sum of total tendency: %f Error: %f" % (np.abs(lhs).sum(), np.sum(np.abs(lhs-salt_rhs))))
    print("---Column integral statistics---")
    print("Maximum difference: %e" % np.max( np.abs( lhs.sum(axis=0)-salt_rhs.sum(axis=0) ) ))
    print("Mean difference: %e" % np.mean( np.abs(lhs.sum(axis=0)-salt_rhs.sum(axis=0)) ) )
    print("Sum difference: %e" % np.sum( np.abs(lhs.sum(axis=0)-salt_rhs.sum(axis=0)) ) )

#    plt.pcolormesh(np.sum(np.abs(lhs-salt_rhs),axis=0),cmap=plt.cm.coolwarm)
#    plt.colorbar()
#    plt.show()



def parse_input_arguments(arguments):
  parser = argparse.ArgumentParser(description=
  '''
  Calculates the difference between the total tendency of salt content cell by cell and returns the globally integrated error summed
  over all time levels in the file unless otherwise specified.
  Assumes the following variables exist in the input file:
  boundary_forcing_salt_tendency
  S_advection_xy
  Sh_tendency_vert_remap
  osaltdiff
  osaltpmdiff
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
  rhs = Dataset(infile).variables['Sh_tendency_vert_remap'][tidx,:,:,:]
  rhs += Dataset(infile).variables['osaltdiff'][tidx,:,:,:]
#  rhs += Dataset(infile).variables['osaltpmdiff'][tidx,:,:,:]
  rhs += Dataset(infile).variables['S_advection_xy'][tidx,:,:,:]
  rhs +=  Dataset(infile).variables['boundary_forcing_salt_tendency'][tidx,:,:,:]
  return rhs

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
