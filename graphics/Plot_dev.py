################################################################################
def PlotCostFn (output_type, html_file, plot_dir, filename, xvals, y1vals, y2vals, y1label, y2label, col1, col2):
  # Plot one or two cost functions
  fig, ax = plt.subplots()

  ax.set_xlabel("iteration", fontsize=16)
  ax.set_ylabel("cost function", fontsize=16)
  ax.set_yscale("log")
  plt.title("Cost function with iteration", fontsize=16)

  # Set the tick label font sizes
  matplotlib.rc("xtick", labelsize=16)
  matplotlib.rc("ytick", labelsize=16)

  ax.plot(xvals[:], y1vals[:], linewidth=2, ls="solid",  color=col1, label=y1label)
  if (len(y2vals) > 0):
    ax.plot(xvals[:], y2vals[:], linewidth=2, ls="solid",  color=col2, label=y2label)

  # Show the legend
  ax.legend()

  if (output_type == "web"):
    graphics_filename = filename + ".png"
  else:
    graphics_filename = filename + ".eps"

  plt.savefig(plot_dir + "/" + graphics_file_name, bbox_inches="tight")
  plt.close("all")

  if (output_type == "web"):
    html_file.write ("<img src=" + graphics_file_name + " width=350>\n")
    html_file.write ("<br>")
  return


################################################################################
def ReadObs (filename):
  # Read-in data in observation space
  print (filename)
  ob_file = open (filename, "r")
  
  obs_tracer         = []
  model_obs_tracer   = []
  obs_tracer_lon     = []
  obs_tracer_lat     = []
  obs_tracer_lev     = []

  obs_flux           = []
  model_obs_flux     = []
  obs_flux_lon       = []
  obs_flux_lat       = []

  obs_totalcol       = []
  model_obs_totalcol = []
  obs_totalcol_lon   = []
  obs_totalcol_lat   = []

  # Read in header
  line    = ob_file.readline()
  line    = ob_file.readline()
  nlevs   = int(line[6:])
  line    = ob_file.readline()
  for l in range(nlevs):
    line = ob_file.readline()

  # Read in the observations
  line = ob_file.readline()
  while (len(line) > 0):
    line   = ob_file.readline()
    obof   = line[7]            # "t" tracer ob, "f" flux ob, "x" total column tracer ob
    #print ("obof = ", obof)
    for l in range(3):
      line = ob_file.readline()
    # Extract the longitude
    lon  = float(line[4:].split()[1])
    #print ("lon = ", lon)
    line = ob_file.readline()
    # Extract the latitude
    lat  = float(line[4:].split()[1])
    #print ("lat = ", lat)
    if (obof == "t"):
      line = ob_file.readline()
      lev  = float(line[4:].split()[1])
      #print ("lev = ", lev)
    line = ob_file.readline()
    obs  = float(line[3:].split()[0])   # Element 0 is actual ob, element 1 is err
    #print ("obs = ", obs)
    line = ob_file.readline()
    mobs = float(line[9:])              # Model obs
    #print ("mobs = ", mobs)
    for l in range(3):
      line = ob_file.readline()

    # Store the observation
    if (obof == "t"):
      # Tracer observation
      obs_tracer_lon.append(lon)
      obs_tracer_lat.append(lat)
      obs_tracer_lev.append(lev)
      obs_tracer.append(obs)
      model_obs_tracer.append(mobs)
    if (obof == "f"):
      # Flux observation
      obs_flux_lon.append(lon)
      obs_flux_lat.append(lat)
      obs_flux.append(obs)
      model_obs_flux.append(mobs)
    if (obof == "x"):
      # Flux observation
      obs_totalcol_lon.append(lon)
      obs_totalcol_lat.append(lat)
      obs_totalcol.append(obs)
      model_obs_totalcol.append(mobs)

  ob_file.close()
  return obs_tracer, model_obs_tracer, obs_tracer_lon, obs_tracer_lat, obs_tracer_lev, obs_flux, model_obs_flux, obs_flux_lon, obs_flux_lat, obs_totalcol, model_obs_totalcol, obs_totalcol_lon, obs_totalcol_lat  


################################################################################
def make_histogram_ob (output_type, html_file, plot_dir, nbins, obs, other_obs, title_part, quantity):
  # Subroutine to generate histogram data of obs - other_obs
  # other_obs could be truth or model obs

  # Number of observations
  Nobs = len(obs)

  # Convert the inputs to numpy arrays
  obs       = np.asarray(obs)
  other_obs = np.asarray(other_obs)

  # Compute the differences
  diff = obs - other_obs

  # Compute the means
  mean = np.mean(diff)

  # Compute the standard deviations
  std = np.std(diff)
  
  # Find the maximum absolute value of these differences
  maxval = 0.0
  maxval = max([maxval, max(abs(diff))])
  maxval = max([maxval, max(abs(diff))])
  maxval = flexi_round(maxval)

  # Generate the bins
  bins   = np.linspace(-1.0*maxval, maxval, nbins+1)
  print ("The following bins have been generated")
  print (bins)

  # Plot histograms
  fig, ax = plt.subplots()
  ax.hist(diff, bins=bins, histtype="bar", facecolor="b")
  ax.set_title(title_part + " for " + quantity + "\nmean=" + str.format("%.5f" % mean) + ", stddev=" + str.format("%.5f" % std) + ", Nobs=" + str(Nobs), fontsize=16)
  ax.set_xlabel(title_part, fontsize=16)
  ax.set_ylabel("frequency", fontsize=16)
  if (output_type == "web"):
    graphics_filename = quantity + "_" + title_part + ".png"
    html_file.write ("<img src=" + graphics_filename  + " width=300></td>\n")
  else:
    graphics_filename = quantity + "_" + title_part + ".eps"
  plt.savefig(plot_dir + "/" + graphics_filename, bbox_inches="tight")
  plt.close("all")

  return

################################################################################
def flexi_round (value):
  s      = np.log10(value)
  si     = int(s-1.0)
  sm     = 10.0 ** float(si)
  modval = int(value/sm + 0.999999999999) * sm
  return modval


################################################################################
def make_scatterplot_ob (output_type, html_file, plot_dir, obs, other_obs, title_part, quantity):
  # Subroutine to do scatter plots of obs - other_obs
  # other_obs could be truth or model obs

  # Find the maximum and minimum values over all data
  maxval = max((obs))
  maxval = max([maxval, max(other_obs)])

  minval = min((obs))
  minval = min([minval, min(other_obs)])

  print ("Minimum value : " + str(minval))
  print ("Maximum value : " + str(maxval))

  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(obs[:], other_obs[:])

  # Do scatter plot
  fig, ax = plt.subplots()
  ax.scatter(obs[:], other_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color="cyan", label="x=y")
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color="red", label="best fit")
  ax.set_title(quantity + " observations vs " + title_part, fontsize=16)
  ax.set_xlabel("observations", fontsize=16)
  ax.set_ylabel(quantity, fontsize=16)
  ax.legend()
  if (output_type == "web"):
    graphics_filename = quantity + "_ob_" + title_part + ".png"
    html_file.write ("<img src=" + graphics_filename  + " width=300><br>\n")
  else:
    graphics_filename = quantity + "_ob_" + title_part + ".eps"

  plt.savefig(plot_dir + "/" + graphics_filename, bbox_inches="tight")
  plt.close("all")
  if (output_type == "web"):
    html_file.write ("m = " + str(slope) + "<br>\n")
    html_file.write ("c = " + str(intercept) + "<br>\n")
    html_file.write ("r = " + str(r_value) + "<br>\n")
    html_file.write ("p = " + str(p_value) + "<br>\n")

  return



################################################################################
def PlotMap (output_type, html_file, plot_dir, field, xs, ys, minplot, maxplot, title_part, level_num, level, quantity, symm=False):
# Plot a single 2d field using basemap

  # Compute some stats for this field
  minfield = min_field2d(field[:,:])
  maxfield = max_field2d(field[:,:])
  # Calculate the rms of term 1
  mean     = np.mean(field[:,:])
  rms      = np.sqrt(np.mean(field[:,:]*field[:,:]))

  print ("\n", field.shape)
  print (minfield, maxfield, mean, rms)
  print (level)

  if (minfield == maxfield):
    if (output_type == "web"):
      html_file.write (name + " field has constant value " + str(minfield) + "<br>\n")
  else:
    # Set-up the basemap projection
    m          = Basemap(projection="hammer",lon_0=0)
    # Make arrays (the same shape as field) containing the longitudes and latitudes of each point
    # Have one more longitude to avoid zero contour along Greenwich meridian
    nlons = len(xs);  nlats = len(ys)
    lons  = np.zeros([nlons+1,nlats])
    lats  = np.zeros([nlons+1,nlats])
    fld   = np.zeros([nlons+1,nlats])
    
    
    # Shift arrays so that -180 longitude is first
    # Find location of first instance at or above 180
    #found = False
    #for lo in range(nlons):
    #  if (not(found) and xs[lo] >= 180.0):
    #    firstindex = lo
    #    found      = True
    #for lo in range(nlons+1):
    #  lo_old_index = lo + firstindex
    #  if (lo_old_index >= nlons):
    #    lo_old_index -= nlons
    #    mod           = 0.0
    #  else:
    #    mod           = 360.0
    #  for la in range(nlats):
    #    lons[lo,la] = xs[lo_old_index] - mod
    #    lats[lo,la] = ys[la]
    #    fld[lo,la]  = field[la,lo_old_index]

    for lo in range(nlons):
      for la in range(nlats):
        lons[lo,la] = xs[lo]
        lats[lo,la] = ys[la]
        fld[lo,la]  = field[la,lo]
    for la in range(nlats):
      lons[nlons,la] = 2.0 * xs[nlons-1] - xs[nlons-2]
      lats[nlons,la] = ys[la]
      fld[nlons,la]  = field[0,la]

    #print ("Longitudes")
    #print (lons[:,10])
    #print ("Latitudes")
    #print (lats[10,:])

    # Is the scale required to be symmetric about zero?
    if symm:
      if (minplot != 0) and (maxplot != 0):
        # The min and max have been specified in the argument list
        # Now make sure they have the same absolute value
        if abs(minplot) > abs(maxplot):
          maxplot = abs(minplot)
        else:
          minplot = -1.0 * abs(minplot)
      else:
        # The min and max are found from the data itself
        if abs(minfield) > abs(maxfield):
          minplot = -1.0 * abs(minfield)
          maxplot = abs(minfield)
        else:
          minplot = -1.0 * abs(maxfield)
          maxplot = abs(maxfield)

    # Adjust the field so it has required min and max (dirty but works)
    minset = False
    maxset = False
    if (minplot != 0) and (maxplot != 0):
      for lo in range(nlons+1):
        for la in range(nlats):
          if fld[lo,la] > maxplot:
            fld[lo,la] = maxplot
            maxset     = True
          else:
            if fld[lo,la] < minplot:
              fld[lo,la] = minplot
              minset     = True

    # Choose the colour scheme
    if symm:
      # Must be a radiant colour scheme
      colsch = "coolwarm"
    else:
      # Not a radiant colour scheme
      #colsch = "Greys"
      #colsch = "gist_rainbow"
      colsch = "YlGnBu"

    # Convert the longitudes and latitudes into basemap co-ordinates
    lons_bm, lats_bm = m(lons, lats)
    # Set-up the colour map
    cmap = cm.get_cmap(colsch, 11)
    # Do the plot
    m.drawcoastlines(linewidth=1.25)
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    cs = m.contourf(lons_bm, lats_bm, fld, cmap=cmap)
    cbar = m.colorbar(cs)
    cs = m.contour(lons_bm, lats_bm, fld)
    title = quantity + " " + title_part + \
            " level ", level #+ "\n" #+ \
    print (title)
    #        "min=" + str.format("%e" % minfield) + \
    #        ", max=" + str.format("%e" % maxfield) + ",\n" + \
    #        "mean=" + str.format("%e" % mean) + \
    #        ", rms=" + str.format("%e" % rms) + "\n"
    plt.title(title, fontsize=14)

    if output_type == "web":
      graphics_filename = quantity + "_" + title_part + "lev" + str(level_num) + ".png"
    else:
      graphics_filename = quantity + "_" + title_part + "lev" + str(level_num) + ".eps"
    plt.savefig(plot_dir + "/" + graphics_filename, bbox_inches="tight")
    plt.close("all")
    if (output_type == "web"):
      html_file.write ("<br><img src=" + graphics_filename + " width=350>\n")
      html_file.write ("<br><hr>")
  return

################################################################################
def max_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

################################################################################
def min_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin



################################################################################
################################################################################
################################################################################
################################################################################

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
from scipy import stats
import os
from mpl_toolkits.basemap import Basemap
import datetime
import sys
#from os.path import exists


# How to run
# python3 Plot_PriorPostInc.py <Assim dir> <Backg dir> <Truth dir> <Backg obs dir> <True obs dir> <web or eps>

# E.g.
# BASE_DIR=/media/ross/banny/Enviflux/RecoverBlobs
# python3 Plot_dev.py $BASE_DIR/Assim/Assim_02 \
#                     $BASE_DIR/Background/Background_02 \
#                     $BASE_DIR/Truth/Truth_02 \
#                     $BASE_DIR/Observations/Bg_02_Observations_02 \
#                     $BASE_DIR/Observations/Observations_02 \
#                     web

# Plots will be placed inside the <Assim dir>/plots
# If "web" is selected, the plots will be png and will be viewable in a web page
# If "eps" is selected, the plots will be eps

num_args = len(sys.argv) - 1
if (num_args != 6):
  print ("Incorrect number of arguments")
  print ("Expecting <Assim dir> <Backg dir> <Truth dir> <Backg obs dir> <True obs dir> <web or eps>")
  print ("Exiting")
  exit (0)

# Set directories
assim_dir   = sys.argv[1]  # Contains anal, anal incs, obs, anal obs
bg_dir      = sys.argv[2]  # Contains bg, bg err
truth_dir   = sys.argv[3]  # Contains truth
bg_obs_dir  = sys.argv[4]  # Contains obs, bg obs
obs_dir     = sys.argv[5]  # Contains obs, true obs
plot_dir    = assim_dir + "/plots"
os.system("mkdir -p " + plot_dir)
output_type = sys.argv[6]

if (output_type != "web") and (output_type != "eps"):
  print ("The output type specified must be either web or eps")
  print ("Exiting")
  exit (0)


# ==============================================================
# \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/
# USER OPTIONS

# Plot cost function?
plot_cost = False

# Plot truth?
plot_truth = False

# Plot background?
plot_bg = False

# Plot analysis?
plot_anal = False

# Plot analysis increments?
plot_anal_inc = False

# Plot errors?
plot_error_field = False

# Plot observation histograms?
plot_ob_hist = True
# Number of bins for histograms
nbins        = 20

# Plot observation scatter plots?
plot_ob_scatter = True

# Plot observation locations?
plot_ob_locations = False

# Plot tracer at vertical levels?
plot_vert_levs = False
# Choose vertical level indices to plot (start at 1)
vert_levs2plot = [1] #, 5] #, 10, 40, 55]

# Plot flux at times?
plot_flux_times = False
# Choose time indices to plot (start at 1)
times2plot     = [1, 2 , 3, 4, 5, 6, 7, 8]

# Plot flux at longitudes?
plot_flux_longs = False
# Choose longitude indices to plot (start at 1)
longs_2plot = [1, 15, 30, 45]

# /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
# ==============================================================

# Find the time now
now          = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")

# Set-up html file for output
if (output_type == "web"):
  filesuffix = ".png"
  html_file = open (plot_dir + "/plots.html", "w")
  html_file.write ("<html>\n")
  html_file.write (current_time + "\n<br>")
  html_file.write ("Assim dir : " + assim_dir)
  html_file.write ("<br>\nBg dir : " + bg_dir)
  html_file.write ("<br>\nTruth dir : " + truth_dir)
  html_file.write ("<br>\n")
else:
  filesuffix = ".eps"
  html_file  = 0




# ==============================================================
if (plot_cost):
  # Read and plot the cost function with iteration
  print ("### Reading and plotting the cost function")
  J_file = open (assim_dir + "/AnalysisDiags.dat", "r")

  # Read in the cost function values for each iteration
  Jb    = []
  Jo    = []
  J     = []
  Jgrad = []
  its   = []
  it    = 0
  line  = J_file.readline()
  line  = J_file.readline()
  line  = J_file.readline()
  while (line != ""):
    split = line.split()
    its.append (float(split[0]))
    Jb.append (float(split[1]))
    Jo.append (float(split[2]))
    J.append (float(split[3]))
    Jgrad.append (float(split[5]))
    line = J_file.readline()

  J_file.close()

  # Plot the total and observational part of the cost function
  PlotCostFn (output_type, html_file, plot_dir, "JJo", its[:], J[:], Jo[:], "J", "Jo", "black", "red")
  # Plot the background part of the cost function
  PlotCostFn (output_type, html_file, plot_dir, "Jb", its[:], Jb[:], [], "J", "", "blue", "")



# ==============================================================
if (plot_ob_hist or plot_ob_scatter):
  # Read observation histograms

  # Read-in and make plots for the actual observations and the true observations
  obs_tracer, true_obs_tracer, obs_tracer_lon, obs_tracer_lat, obs_tracer_lev, \
    obs_flux, true_obs_flux, obs_flux_lon, obs_flux_lat, \
    obs_totalcol, true_obs_totalcol, obs_totalcol_lon, obs_totalcol_lat = ReadObs (obs_dir+"/Observations.dat")

  print ("There are", len(obs_tracer), " actual observations of the tracer")
  print ("There are", len(true_obs_tracer), " truth observations of the tracer")
  print ("There are", len(obs_flux), " actual observations of the flux")
  print ("There are", len(true_obs_flux), " truth observations of the flux")
  print ("There are", len(obs_totalcol), " actual observations of the total column tracer")
  print ("There are", len(true_obs_totalcol), " truth observations of the total column tracer")

  if (plot_ob_hist):
    if len(obs_tracer) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_tracer, true_obs_tracer, "o-t", "tracer")
    if len(obs_flux) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_flux, true_obs_flux, "o-t", "flux")
    if len(obs_totalcol) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_totalcol, true_obs_totalcol, "o-t", "x")
  if (plot_ob_scatter):
    if len(obs_tracer) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_tracer, true_obs_tracer, "truth", "tracer")
    if len(obs_flux) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_flux, true_obs_flux, "truth", "flux")
    if len(obs_totalcol) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_totalcol, true_obs_totalcol, "truth", "x")


  # Read-in and make plots for the background observations
  nil, bg_obs_tracer, nil, nil, nil, \
    nil, bg_obs_flux, nil, nil, \
    nil, bg_obs_totalcol, nil, nil = ReadObs (bg_obs_dir+"/Observations.dat")
  print ("There are", len(bg_obs_tracer), " background observations of the tracer")
  print ("There are", len(bg_obs_flux), " background observations of the flux")
  print ("There are", len(bg_obs_totalcol), " background observations of the total column tracer")
  if (plot_ob_hist):
    if len(obs_tracer) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_tracer, bg_obs_tracer, "o-b", "tracer")
    if len(obs_flux) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_flux, bg_obs_flux, "o-b", "flux")
    if len(obs_totalcol) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_totalcol, bg_obs_totalcol, "o-b", "x")
  if (plot_ob_scatter):
    if len(obs_tracer) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_tracer, bg_obs_tracer, "background", "tracer")
    if len(obs_flux) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_flux, bg_obs_flux, "background", "flux")
    if len(obs_totalcol) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_totalcol, bg_obs_totalcol, "background", "x")

  # Read-in and make plots for the analysis observations
  nil, anal_obs_tracer, nil, nil, nil, \
    nil, anal_obs_flux, nil, nil, \
    nil, anal_obs_totalcol, nil, nil = ReadObs (assim_dir+"/ObsPostAnal.dat")
  print ("There are", len(anal_obs_tracer), " analysis observations of the tracer")
  print ("There are", len(anal_obs_flux), " analysis observations of the flux")
  print ("There are", len(anal_obs_totalcol), " analysis observations of the total column tracer")
  if (plot_ob_hist):
    if len(obs_tracer) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_tracer, anal_obs_tracer, "o-a", "tracer")
    if len(obs_flux) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_flux, anal_obs_flux, "o-a", "flux")
    if len(obs_totalcol) > 0:
      make_histogram_ob (output_type, html_file, plot_dir,
                         nbins, obs_totalcol, anal_obs_totalcol, "o-a", "x")
  if (plot_ob_scatter):
    if len(obs_tracer) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_tracer, anal_obs_tracer, "analysis", "tracer")
    if len(obs_flux) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_flux, anal_obs_flux, "analysis", "flux")
    if len(obs_totalcol) > 0:
      make_scatterplot_ob (output_type, html_file, plot_dir,
                           obs_totalcol, anal_obs_totalcol, "analysis", "x")

  
# ==============================================================

# Read in tracer fields
if plot_truth or plot_bg or plot_anal or plot_anal_inc or plot_error_field:
  # Read in longitudes, latitudes, levels, and times
  nc_file = Dataset(truth_dir + "/Truth.nc")
  lons    = nc_file.variables["longitude"][:]
  lats    = nc_file.variables["latitude"][:]
  levs    = nc_file.variables["level"][:]
  times   = nc_file.variables["time"][:]
  nc_file.close()

if plot_truth:
  # Read in the truth
  nc_file = Dataset(truth_dir + "/Truth.nc")
  truth_tracer = nc_file.variables["tracer0"][:][:][:]
  truth_source = nc_file.variables["source"][:][:][:]

  

# Plot tracer on horizontal levels

if plot_vert_levs:
  for lev in vert_levs2plot:
    if plot_truth:
      PlotMap (output_type, html_file, plot_dir,
               truth_tracer[lev-1][:][:], lons[:], lats[:],
               0.0, 0.0,
               "truth", lev, str.format("%f4.1", levs[lev-1]), "tracer", False)







if (output_type == "web"):
  html_file.write ("</html>")
  html_file.close()
  print ("Go to a web browser and paste the following as a URL")
  print ("file://" + plot_dir + "/plots.html")
else:
  print ("eps files are in the following directory")
  print (plot_dir)
