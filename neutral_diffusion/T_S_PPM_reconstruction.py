#!/usr/bin/env python

import os
import sys
import argparse
from time import time
import numpy as np
from netCDF4 import Dataset
from gsw import CT_from_pt, rho, p_from_z
import matplotlib.pyplot as plt
from ppm_routines import ppm_routines

def reconstruct_cell_edges(nk, h, pressure, s):
  s_int = np.zeros(nk+1,dtype=np.float32)
  s_L = np.zeros(nk,dtype=np.float32)
  s_R = np.zeros(nk,dtype=np.float32)
  ppm_routines.interface_scalar(h,s,s_int,2)
  ppm_routines.ppm_left_right_edge_values(s_L, s_R, s.astype(np.float32),s_int)
  return s_L, s_R

def main(arguments):

  parser = argparse.ArgumentParser(description=
  '''
  Calculate and plot the PPM reconstructions of T, S and subsequent density based on layer averages
  ''',
  epilog="Written by A.E.Shao 2017")

  parser.add_argument('infile',   help="Input file",                    type=str)
  parser.add_argument('tempvar',  help="Name of temperature variable",  type=str)
  parser.add_argument('saltvar',  help="Name of salinity variable",     type=str)
  parser.add_argument('hvar',     help="Name of thickness variable",    type=str)
  parser.add_argument('latidx',   help="Latitude index",                type=int)
  parser.add_argument('lonidx',   help="Longitude index",               type=int)
  parser.add_argument('--save',   help="Save the resulting figure", action='store_true')
  args = parser.parse_args(arguments)

  print(args)

  h =     Dataset(args.infile).variables[args.hvar][0,:,args.latidx,args.lonidx].squeeze()
  salt =  Dataset(args.infile).variables[args.saltvar][0,:,args.latidx,args.lonidx].squeeze()
  temp =  Dataset(args.infile).variables[args.tempvar][0,:,args.latidx,args.lonidx].squeeze()
  lat =   Dataset(args.infile).variables['lath'][args.latidx]
  depth  = p_from_z(-(h.cumsum()-h/2.),lat)

  nk = salt.size
  pressure = np.zeros(nk+1,dtype=np.float32)
  pressure[1:] = p_from_z(-h.cumsum(),np.ones(nk)*lat)

  temp_L, temp_R = reconstruct_cell_edges(nk,h,pressure,temp)
  salt_L, salt_R = reconstruct_cell_edges(nk,h,pressure,salt)
  rho_avg = rho(salt,temp,depth)
  rho_L = rho(salt_L, temp_L, pressure[0:-1])
  rho_R = rho(salt_R, temp_R, pressure[1:])

  plt.figure()
  plt.plot(rho_avg,depth,marker='o')
  plt.scatter(rho_L,pressure[0:-1],marker='+')
  plt.scatter(rho_R,pressure[1:], marker='x')
  plt.xlim(1030,1037)
  plt.ylim(500,2000)
  plt.gca().invert_yaxis()
  plt.grid()

#  plt.figure()
#  plt.subplot(1,3,1)
#  plt.plot(temp,depth,marker='o')
#  plt.scatter(temp_L,pressure[0:-1],marker='+')
#  plt.scatter(temp_R,pressure[1:], marker='x')
#  plt.subplot(1,3,2)
#  plt.plot(salt,depth,marker='o')
#  plt.scatter(salt_L,pressure[0:-1],marker='+')
#  plt.scatter(salt_R,pressure[1:], marker='x')
#  plt.subplot(1,3,3)
#  plt.plot(rho_avg,depth,marker='o')
#  plt.scatter(rho_L,pressure[0:-1],marker='+')
#  plt.scatter(rho_R,pressure[1:], marker='x')
  if (args.save):
    plt.savefig("%03d.%03d.png" % (args.latidx,args.lonidx))
  plt.show(block=True)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


