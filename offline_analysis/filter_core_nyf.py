from netCDF4 import Dataset
from scipy import signal
from os import sep
from glob import glob
import sys
import argparse

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data, axis=0)
    return y

def main(arguments):
  parser = argparse.ArgumentParser(description=
            '''
            Filter CORE forcing for MOM6 using a lowpass Butterworth filter using a specified critical value.
            WARNING: for simplicity, this overwrite the fields in the input files
            ''',
            epilog = "Written by A.E. Shao 2017")
  parser.add_argument('inpath' , help='Path to the directory where the CORE fields are stored', type = str)
  parser.add_argument('cutfreq', help='Any frequencies below this will be filtered. Should be 1/day', type = float)
  args = parser.parse_args()

  core_files = glob(args.inpath + sep + '*.nc')

  for infile in core_files:
    print("Processing variables in: " + infile)
    nc_file = Dataset(infile,'r+')
    print(len(nc_file.variables['TIME']))
    if len(nc_file.variables['TIME']) >= 73:
      fs = 1./(nc_file.variables['TIME'][1] - nc_file.variables['TIME'][0])
      b, a = butter_lowpass(args.cutfreq, fs, order = 2)
      for var in nc_file.variables:
        # Only work on variables which have time, lat, lon dimensions
        if len(nc_file.variables[var].shape) == 3:
          print("Filtering: " + var)
          nc_file.variables[var][:,:,:] = signal.filtfilt(b, a, nc_file.variables[var][:,:,:], axis=0)
    nc_file.close()

if __name__ == '__main__':
  sys.exit(main(sys.argv[1:]))
