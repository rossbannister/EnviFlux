# ===================================================================
def CalcFluxes (flux, xs, ys):
  #############################################################
  ### Calculate and report the total flux                   ###
  #############################################################
  RE            = 6371000.0 # in m
  deg2rad       = np.pi / 180.0
  nlon          = len(xs)
  nlat          = len(ys)
  dx_rad        = deg2rad * 360.0 / float(nlon)
  nd            = 30.0 # Number of days in our month

  flux_tot         = 0.0
  flux_tot_source  = 0.0
  flux_tot_sink    = 0.0
  flux_peak_source = 0.0
  flux_peak_sink   = 0.0
  lon_peak_source  = 0.0
  lat_peak_source  = 0.0
  lon_peak_sink    = 0.0
  lat_peak_sink    = 0.0

  # Compute the total fluxes
  # Note: longitudes are ascending, latitudes are descending
  for lat in range(nlat):
    if (lat == 0):
      uppery = 90.0
    else:
      uppery = (ys[lat-1] + ys[lat]) / 2.0
    if (lat == nlat-1):
      lowery = -90.0
    else:
      lowery = (ys[lat] + ys[lat+1]) / 2.0
    #print ('Latitude', lat, ' =', ys[lat], ' lower, upper =', lowery, uppery)
    dy = uppery - lowery   # Size of this grid box in y (degrees)
    # Area of this grid box in square m
    dA = RE * RE * np.cos(deg2rad * ys[lat]) * dx_rad * dy * deg2rad
    
    for lon in range(nlon):
      flux_cont  = flux[lat,lon] * dA  # micrograms per second
      flux_tot  += flux_cont
      # Divide into pos and neg
      if flux_cont > 0.0:
        flux_tot_source += flux_cont
      else:
        flux_tot_sink   += flux_cont
  
  # The flux variables are in micrograms of methane per second
  # We want mass (in Tg) of methane per month
  flux_tot        *= 24.0 * 3600.0 * nd / 1.0E18
  flux_tot_source *= 24.0 * 3600.0 * nd / 1.0E18
  flux_tot_sink   *= 24.0 * 3600.0 * nd / 1.0E18


  # Find the minimum and maximum flux values, and the positions
  for lat in range(nlat):
    for lon in range(nlon):
      flux_cont = flux[lat,lon]
      if flux_cont > 0.0:
        # This is a source point
        if flux_cont > flux_peak_source:
          flux_peak_source = flux_cont
          lon_peak_source  = xs[lon]
          lat_peak_source  = ys[lat]
      if flux_cont < 0.0:
        # This is a sink point
        if np.abs(flux_cont) > flux_peak_sink:
          flux_peak_sink   = np.abs(flux_cont)
          lon_peak_sink    = xs[lon]
          lat_peak_sink    = ys[lat]

  if lon_peak_sink > 180.0:
    lon_peak_sink -= 360.0
  if lon_peak_source > 180.0:
    lon_peak_source -= 360.0


  #print (flux_peak_source, flux_peak_sink, flux_tot_source, flux_tot_sink, lon_peak_source, lat_peak_source, lon_peak_sink, lat_peak_sink)
  print ('%14.8f %14.8f %14.8f %14.8f %10.1f %10.1f %10.1f %10.1f' %(flux_peak_source, flux_peak_sink, flux_tot_source, flux_tot_sink, lon_peak_source, lat_peak_source, lon_peak_sink, lat_peak_sink))


  return



################################################################################
################################################################################
################################################################################
################################################################################

import numpy as np
from netCDF4 import Dataset
import os
import datetime
import sys


# How to run
# python3 ExtractTotalsPeaksLocs.py


# Set directories
base_dir   = '/media/ross/banny/EnviFlux/RecoverBlobs/Sept2025_experiments'
truth_dir  = 'Truth'
bg_dir     = 'Background'
#assim_dirs = ['TracerFac1.0_FluxFac0.01_wFac0.6_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac0.7_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac0.8_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac0.95_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac0.9_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.05_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.1_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.2_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.3_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.4_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac1.0_Bias0.0', 'TracerFac1.0_FluxFac0.01_wFac0.6_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac0.7_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac0.8_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac0.95_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac0.9_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac1.05_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac1.0_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac1.1_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac1.2_Bias0.0_tob_z0', 'TracerFac1.0_FluxFac0.01_wFac1.0_Bias0.0_ReverseSignFlux', 'TracerFac1.0_FluxFac0.01_wFac0.8_Bias0.0_ReverseSignFlux', 'TracerFac1.0_FluxFac0.01_wFac1.3_Bias0.0_tob_z0']
assim_dirs = ['TracerFac1.0_FluxFac0.01_wFac1.4_Bias0.0_tob_z0']


# The data come from netcdf data file (including the longitudes, latitudes, and heights)


# Truth (and get longitudes, latiudes, levels, and times)
print (base_dir + '/' + truth_dir + '/Truth.nc')
nc_file          = Dataset(base_dir + '/' + truth_dir + '/Truth.nc')
Field_flux_truth = nc_file.variables['source'][:][:][:]
longs            = nc_file.variables['longitude'][:]
lats             = nc_file.variables['latitude'][:]
#print ('Longitudes : ', longs[:])
#print ('Latitudes : ', lats[:])
nc_file.close()
print ('Data for the truth')
CalcFluxes (Field_flux_truth[0][:][:], longs[:], lats[:])


# Background
print (base_dir + '/' + bg_dir + '/Background.nc')
nc_file       = Dataset(base_dir + '/' + bg_dir + '/Background.nc')
Field_flux_bg = nc_file.variables['source'][:][:][:]
nc_file.close()
print ('Data for the background')
CalcFluxes (Field_flux_bg[0][:][:], longs[:], lats[:])

# Analyses
for anal in assim_dirs:
  print (base_dir + '/Assim/' + anal + '/Anal.nc')
  nc_file         = Dataset(base_dir + '/Assim/' + anal + '/Anal.nc')
  Field_flux_anal = nc_file.variables['source'][:][:][:]
  nc_file.close()
  print ('Data for ' + anal)
  CalcFluxes (Field_flux_anal[0][:][:], longs[:], lats[:])
