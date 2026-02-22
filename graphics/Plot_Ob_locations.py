# ===================================================================
# ===================================================================
def ReadObs (ObsFile):
  print ('### Reading the observational information, including model observations')
  Obs_file = open (ObsFile, 'r')
  line  = Obs_file.readline()
  line  = Obs_file.readline()
  print (line)
  split = line.split()
  nlevs = int(split[1])   # Needed to read-in the right number of lines of a profile at the start of the file
  for lev in range(nlevs+1):
    line = Obs_file.readline()

  # Set-up the arrays of observation information
  ob_time_s    = []
  ob_lon       = []
  ob_lat       = []
  ob_lev       = []
  ob_val       = []
  ob_stddev    = []
  ob_model     = []
  line         = Obs_file.readline()
  keep_reading = line != ''

  while keep_reading:
    line_ob_of   = Obs_file.readline()
    split        = line_ob_of.split()
    ob_of        = split[1]
    line_ob_type = Obs_file.readline()
    line_time    = Obs_file.readline()
    line_lon     = Obs_file.readline()
    line_lat     = Obs_file.readline()
    if ob_of == 't':
      line_lev   = Obs_file.readline()
    line_ob      = Obs_file.readline()
    line_modelob = Obs_file.readline()
    line_innov   = Obs_file.readline()
    line_grad    = Obs_file.readline()

    split = line_time.split()
    ob_time_s.append(float(split[5]))
    split = line_lon.split()
    ob_lon.append(float(split[2]))
    split = line_lat.split()
    ob_lat.append(float(split[2]))
    if ob_of == 't':
      split = line_lev.split()
      ob_lev.append(float(split[2]))
    else:
      ob_lev.append(0.0)
    split = line_ob.split()
    ob_val.append(float(split[1]))
    ob_stddev.append(float(split[2]))
    split = line_modelob.split()
    ob_model.append(float(split[1]))

    # In preparation for the next observation
    line         = Obs_file.readline()
    keep_reading = line != ''

  Obs_file.close()
  return ob_time_s, ob_lon, ob_lat, ob_lev, ob_val, ob_stddev, ob_model



# ===================================================================
def plot_ob_hts (obs_longs, obs_lats, obs_alts, quantity, output_type, plot_dir, html_file):
  # Subroutine to plot observation heights

  obs_longs_station = [] 
  obs_lats_station  = []
  obs_alts_station  = []
  n_stations        = 0
  nobs              = len(obs_longs)
  crit              = 0.00001
  max_alt           = 0.0

  if (output_type != 'png'):
    print ('Eliminating repeated observation locations')
    for ob in range(nobs):
      this_long = obs_longs[ob]
      this_lat  = obs_lats[ob]
      this_alt  = obs_alts[ob]
      found     = False
      for check in range(n_stations):
        # Is there a match?
        found = (abs(this_long - obs_longs_station[check]) < crit) and (abs(this_lat - obs_lats_station[check]) < crit) and (abs(this_alt - obs_alts_station[check]) < crit)
        if (found):
          break
      if not(found):
        # Add a new station
        obs_longs_station.append(this_long)
        obs_lats_station.append(this_lat)
        obs_alts_station.append(this_alt)
        if (this_alt > max_alt):
          max_alt = this_alt
        n_stations += 1
  else:
    # Keep all observations
    obs_longs_station = obs_longs
    obs_lats_station  = obs_lats
    obs_alts_station  = obs_alts

  print ('Maximum station alitude ', max_alt)

  # Find size of dot
  obs_dot_station = []
  for this_alt in obs_alts_station:
    if max_alt != 0.0:
      dot_size = 100.0 * this_alt / max_alt
    else:
      dot_size = 100.0 * this_alt
    if (dot_size < 0.01):
      dot_size = 0.01
    obs_dot_station.append(dot_size)
  

  if (output_type == 'png'):
    filesuffix = '.png'
    html_file.write ('<h4>Observation locations for ' + quantity + '</h4>\n')
  else:
    filesuffix = '.eps'


  m = Basemap(projection='hammer',lon_0=0)

  # Convert the longitudes and latitudes to co-ordinates that basemap can use
  x, y = m(obs_longs_station, obs_lats_station)

  # Draw stuff
  m.drawmapboundary(fill_color='#99ffff')
  m.fillcontinents(color='#cc9966',lake_color='#99ffff')

  # Plot the locations
  m.scatter(x, y, obs_dot_station, marker='o', color='black')
  plt.title('Locations and alt (max alt ' + str(max_alt) + 'm)', fontsize=16)

  figfilename = quantity + '_locations' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')

  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300></td><br>\n')
  return



################################################################################
################################################################################
################################################################################
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
import os
from mpl_toolkits.basemap import Basemap
import datetime
import sys


# How to run
# python3 Plot_Ob_locations.py <Observations dir> <png or eps>

# Observations dir: Directory continaing the actial observations

# Plots will be placed inside the <Observations dir>/plots
# If "png" is selected, the plots will be png and will be viewable in a web page
# If "eps" is selected, the plots will be eps

num_args = len(sys.argv) - 1
if num_args != 2:
  print ("Incorrect number of arguments")
  print ("Expecting <Observations dir> <png or eps>")
  print ("Exiting")
  exit (0)

# Set directories
obs_dir     = sys.argv[1]
plot_dir    = obs_dir + '/plots'
output_type = sys.argv[2]

# Make the required directories
os.system('mkdir -p ' + plot_dir + '/graphics')

if (output_type != "png") and (output_type != "eps"):
  print ("The output type specified must be either png or eps")
  print ("Exiting")
  exit (0)



# Find the time now
now          = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")

# Set-up html file for output
if output_type == 'png':
  filesuffix = '.png'
  html_file = open (plot_dir + '/plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write (current_time + '\n<br>\n')
  html_file.write ("Observations dir : " + obs_dir + '\n<br>\n')
else:
  filesuffix = '.eps'
  html_file  = 0






################################################################################
## CODE TO DO WITH OBSERVATIONS ################################################
################################################################################

ob_time_s, ob_lon, ob_lat, ob_lev, ob_val, ob_stddev, ob_model = ReadObs (obs_dir + '/' + 'Observations.dat')
nobs  = len(ob_val)
print ('Number of observations read in: ', nobs)

plot_ob_hts (ob_lon, ob_lat, ob_lev, 'tracer', output_type, plot_dir, html_file)


if (output_type == 'png'):
  html_file.write ('</html>')
  html_file.close()
  print ('An html file has been created to view the figures')
  print ('Please view the following html file with your browser')
  print (plot_dir + '/plots.html')
  print ('')
