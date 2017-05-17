#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
from time import time
from numpy import sum, save, squeeze, zeros
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
  parser.add_argument('outfile', help="Path to file where output will be stored", type=str)
  args = parser.parse_args(arguments)

  print(args)

  time_in = Dataset(args.infile).variables['time'][:]
  ntime = time_in.size

  h = Dataset(args.hfile).variables[args.hvarname][:,:,:,:].squeeze()

  nk, nlat, nlon = h.shape
  # Allocate the output array
  column_integral = zeros((ntime,nlat,nlon))

  print('Starting column integration')
  time();
  tidx = 0
  while (tidx < ntime):
    beg_time = time()
    conc = Dataset(args.infile).variables[args.invar][tidx,:,:,:].squeeze()
    column_integral[tidx,:,:] = squeeze(sum( conc*h, axis=0 ))
    end_time = time()
    print('Time elapsed processing timelevel %03d: %f' % (tidx,end_time-beg_time))
    tidx = tidx + 1

  print('Saving output field as %s' % args.outfile)
  save(args.outfile, column_integral)
  print('Complete. Total time elapsed: %f min' % time())

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
