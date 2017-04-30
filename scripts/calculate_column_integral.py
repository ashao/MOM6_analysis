#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
from time import time as
from numpy import sum, save
from netCDF4 import Dataset

def main(arguments):

  parser = argparse.ArgumentParser(description=
  '''
  Calculates the column integral of a field from MOM6 output
  ''',
  epilog="Written by A.E.Shao 2017")

  parser.add_argument('infile', help="Input file", type=str)
  parser.add_argument('invar', help="Variable to calculate the column integral", type=str)
  parser.add_argument('hfile', help="File containing layer thicknesses", type=str)
  parser.add_argument('hvarname', help="Name of thickness variable", type=str)
  args = parser.parse_args(arguments)

  print(args)

  time_in = Dataset(arguments.infile).variables['time'][:]
  ntime = time_in.size

  h = Dataset(arguments.hfile).variables[hvarname][:,:,:,:].squeeze()

  nk, nlat, nlon = h.shape
  # Allocate the output array
  column_integral = np.zeros((ntime,nlat,nlon))

  print('Starting column integration')
  time();
  tidx = 0
  while (tidx < ntime):
    beg_time = time()
    conc = Dataset(infile).variables[invar][tidx,:,:,:].squeeze()
    column_integral[tidx,:,:] = np.squeeze(np.sum( conc*h, axis=0 ))
    end_time = time()
    print('Time elapsed processing timelevel %03d: %f' % (tidx,beg_time-end_time))

  print('Complete. Total time elapsed: %f min' % time()/60.)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
