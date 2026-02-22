# ===================================================================
def max_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

# ===================================================================
def min_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin

# ===================================================================
def max_field3d (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max_field2d(lev))
  globalmax = max(maxlev)
  return globalmax

# ===================================================================
def min_field3d (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min_field2d(lev))
  globalmin = min(minlev)
  return globalmin


# ===================================================================
def flexi_round (value):
  s      = np.log10(value)
  si     = int(s-1.0)
  sm     = 10.0 ** float(si)
  modval = int(value/sm + 0.999999999999) * sm
  return modval


# ===================================================================
def PlotJ (J, Jb, Jo, its, html_file, plotdir, output_type):
  fig, ax = plt.subplots()
  ax.set_xlabel('iteration', fontsize=16)
  ax.set_ylabel('cost function', fontsize=16)
  ax.set_yscale('log')
  plt.title('Cost function with iteration', fontsize=16)

  # Set the tick label font sizes
  matplotlib.rc('xtick', labelsize=16)
  matplotlib.rc('ytick', labelsize=16)

  ax.plot(its[:], J[:],       linewidth=2, ls='solid',  color='black', label='J')
  ax.plot(its[:], Jb[:],      linewidth=2, ls='solid',  color='blue', label='Jb')
  ax.plot(its[:], Jo[:],      linewidth=2, ls='solid',  color='red', label='Jo')

  # Show the legend
  ax.legend()

  graphics_file_name = 'graphics/CostFn.' + output_type
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<img src=' + graphics_file_name + ' width=350>\n')
    html_file.write ('<br>')



# ===================================================================
def Plot2d (field, xs, ys, xname, yname, title, name, reverse_axes, html_file, plotdir, output_type, minplot=0, maxplot=0, symm=False):
  #############################################################
  ### Plot a single 2d field                                ###
  #############################################################

  #print ('Field shape')
  #print (field.shape)

  # Compute some stats for this field
  minfield = min_field2d(field[:,:])
  maxfield = max_field2d(field[:,:])
  # Calculate the rms of term 1
  mean     = np.mean(field[:,:])
  rms      = np.sqrt(np.mean(field[:,:]*field[:,:]))
  diff     = field[:,:] - mean
  rmsdev   = np.sqrt(np.mean(diff[:,:]*diff[:,:]))

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
  else:
    if (minplot == 0) and (maxplot == 0):
      minplot = minfield
      maxplot = maxfield

  # Choose the colour scheme
  if symm:
    # Must be a radiant colour scheme
    colsch = 'coolwarm'
  else:
    # Not a radiant colour scheme
    #colsch = 'Greys'
    #colsch = 'gist_rainbow'
    colsch = 'YlGnBu'





  if (minfield == maxfield):
    if (output_type == 'png'):
      html_file.write (name + ' field has constant value ' + str(minfield) + '<br>\n')
  else:
    matplotlib.rc('xtick', labelsize=16)
    matplotlib.rc('ytick', labelsize=16)
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    levels = np.linspace(minplot, maxplot, 10)
    cmap = cm.get_cmap(colsch, 11)
    cax  = ax.contourf(xs[:], ys[:], field[:,:], cmap=cmap, levels=levels, vmin=minplot, vmax=maxplot)
    cbar = fig.colorbar(cax)
    cax  = ax.contour (xs[:], ys[:], field[:,:], colors='k')
    ax.set_title(title, fontsize=14)
    #ax.set_title(title + '\n' +
    #             'min=' + str.format('%e' % minfield) +
    #             ', max=' + str.format('%e' % maxfield) + ',\n' +
    #             'mean=' + str.format('%e' % mean) +
    #             ', rms=' + str.format('%e' % rms) + '\n' +
    #             ', rmsdev=' + str.format('%e' % rmsdev) + '\n'
    #             , fontsize=14)

    if (reverse_axes):
      ax.invert_xaxis()
      ax.invert_yaxis()

    ax.set_xlabel(xname, fontsize=16)
    ax.set_ylabel(yname, fontsize=16)
    #plt.show()
    graphics_file_name = 'graphics/' + name + '.' + output_type
    plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
    plt.close('all')
    if (output_type == 'png'):
      html_file.write ('<img src=' + graphics_file_name + ' width=600>\n')


# ===================================================================
def PlotMap (field, xs, ys, xname, yname, title, name, html_file, plotdir, output_type, minplot=0, maxplot=0, symm=False):
  #############################################################
  ### Plot a single 2d field using basemap                  ###
  #############################################################

  # Compute some stats for this field
  minfield = min_field2d(field[:,:])
  maxfield = max_field2d(field[:,:])
  # Calculate the rms of term 1
  mean     = np.mean(field[:,:])
  rms      = np.sqrt(np.mean(field[:,:]*field[:,:]))
  diff     = field[:,:] - mean
  rmsdev   = np.sqrt(np.mean(diff[:,:]*diff[:,:]))

  #print ('\n*** ', minfield, maxfield)

  if minfield == maxfield:
    if output_type == 'png':
      html_file.write (name + ' field has constant value ' + str(minfield) + '<br>\n')
  else:
    # Set-up the basemap projection
    #m = Basemap(projection='hammer', lon_0=0)
    m = Basemap(projection='moll', lon_0=0)

    # Make arrays (the same shape as field) containing the longitudes and latitudes of each point
    # Have one more longitude to avoid zero contour along Greenwich meridian
    nlons = len(xs);  nlats = len(ys)
    lons  = np.zeros([nlons,nlats])
    lats  = np.zeros([nlons,nlats])
    fld   = np.zeros([nlons,nlats])
    #for lo in range(nlons):
    #  for la in range(nlats):
    #    lons[lo,la] = xs[lo]
    #    lats[lo,la] = ys[la]
    #    fld[lo,la]  = field[la,lo]
    ## Add extra longitude
    #for la in range(nlats):
    #  lons[nlons,la]   = xs[nlons-1] + xs[nlons-1] - xs[nlons-2]
    #  lats[nlons,la]   = ys[la]
    #  field[nlons,la]  = field[0,la]


    # Shift arrays so that -180 longitude is first
    # Find location of first instance at or above 180
    found = False
    for lo in range(nlons):
      if (not(found) and xs[lo] >= 180.0):
        firstindex = lo
        found      = True
    for lo in range(nlons):
      lo_old_index = lo + firstindex
      if (lo_old_index >= nlons):
        lo_old_index -= nlons
        mod           = 0.0
      else:
        mod           = 360.0
      for la in range(nlats):
        lons[lo,la] = xs[lo_old_index] - mod
        lats[lo,la] = ys[la]
        fld[lo,la]  = field[la,lo_old_index]



    print ('\n')
    print ('Min of basemap field =', min_field2d(fld[:,:]))
    print ('Max of basemap field =', max_field2d(fld[:,:]))


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
    else:
      if (minplot == 0) and (maxplot == 0):
        minplot = minfield
        maxplot = maxfield

    # Choose the colour scheme
    if symm:
      # Must be a radiant colour scheme
      colsch = 'coolwarm'
    else:
      # Not a radiant colour scheme
      #colsch = 'Greys'
      #colsch = 'gist_rainbow'
      colsch = 'YlGnBu'

    # Convert the longitudes and latitudes into basemap co-ordinates
    lons_bm, lats_bm = m(lons, lats)
    # Set-up the colour map
    cmap = cm.get_cmap(colsch, 11)
    # Do the plot

    m.drawcoastlines(linewidth=1.25)
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    #print ('\n*** Using min and max: ', minplot, maxplot,)
    levels = np.linspace(minplot, maxplot, 10)
    cs = m.contourf(lons_bm, lats_bm, fld, cmap=cmap, levels=levels, vmin=minplot, vmax=maxplot)
    cbar = m.colorbar(cs)
    #cs = m.contour(lons_bm, lats_bm, fld)

    # Conventional plot
    #fig, ax = plt.subplots()
    #cax  = ax.contourf(xs, ys, field, cmap=cmap, levels=levels, vmin=minplot, vmax=maxplot)
    #cbar = fig.colorbar(cax, orientation='vertical', cmap=cmap)


    # Labels etc
    plt.title(title, fontsize=14)
    #plt.title(title + '\n' +
    #          'min=' + str.format('%e' % minfield) +
    #          ', max=' + str.format('%e' % maxfield) + ',\n' +
    #          'mean=' + str.format('%e' % mean) +
    #          ', rms=' + str.format('%e' % rms) + '\n' +
    #          ', rmsdev=' + str.format('%e' % rmsdev) + '\n'
    #          , fontsize=14)

    graphics_file_name = 'graphics/' + name + '.' + output_type
    plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
    plt.close('all')
    if output_type == 'png':
      html_file.write ('<br><img src=' + graphics_file_name + ' width=600>\n')

  return

# ===================================================================
def CalcTotalFlux (flux, xs, ys, landsea, name, html_file, output_type):
  #############################################################
  ### Calculate and report the total flux                   ###
  #############################################################
  RE           = 6371000.0 # in m
  deg2rad      = np.pi / 180.0
  nlon         = len(xs)
  nlat         = len(ys)
  dx_rad       = deg2rad * 360.0 / float(nlon)
  flux_tot     = 0.0
  flux_land    = 0.0
  flux_sea     = 0.0
  flux_tot_pos = 0.0
  flux_tot_neg = 0.0
  nd           = 30.0 # Number of days in our month

  #print ('nlon, nlat', nlon, nlat)
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
      # Divide into land and sea
      if landsea[lon,lat] == 1:
        # Contribution to land
        flux_land += flux_cont
      else:
        # Contribution to land
        flux_sea  += flux_cont
      # Divide into source or sink
      if flux_cont > 0.0:
        flux_tot_pos += flux_cont
      elif flux_cont < 0.0:
        flux_tot_neg += flux_cont

  # The flux variables are in micrograms of methane per second
  # We want mass (in Tg) of methane per month
  flux_tot     *= 24.0 * 3600.0 * nd / 1.0E18
  flux_land    *= 24.0 * 3600.0 * nd / 1.0E18
  flux_sea     *= 24.0 * 3600.0 * nd / 1.0E18
  flux_tot_pos *= 24.0 * 3600.0 * nd / 1.0E18
  flux_tot_neg *= 24.0 * 3600.0 * nd / 1.0E18

  # Find the minimum and maximum flux values
  min_flux = min_field2d (flux)
  max_flux = max_field2d (flux)

  #print ('Global flux = ', flux_tot)
  if (output_type == 'png'):
    html_file.write ('<br>Total flux calculations for ' + name + ' (Tg methane per month)\n')
    html_file.write ('<br>Total flux   = ' + str(flux_tot) + '\n')
    html_file.write ('<br>Land flux    = ' + str(flux_land) + '\n')
    html_file.write ('<br>Sea flux     = ' + str(flux_sea) + '\n')
    html_file.write ('<br>Total source = ' + str(flux_tot_pos) + '\n')
    html_file.write ('<br>Total sink   = ' + str(flux_tot_neg) + '\n')
    html_file.write ('<br>Min and max fluxes for ' + name + ' (ug methane per sq m per s)')
    html_file.write ('<br>Min flux     = ' + str(min_flux))
    html_file.write ('<br>Max flux     = ' + str(max_flux))


  return flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg


# ===================================================================
def ReadObsBg (ObsFile):
  print ('### Reading the observational information, including background model observations')
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
  ob_bg        = []
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
    ob_bg.append(float(split[1]))

    # In preparation for the next observation
    line         = Obs_file.readline()
    keep_reading = line != ''

  Obs_file.close()
  return ob_time_s, ob_lon, ob_lat, ob_lev, ob_val, ob_stddev, ob_bg
  
# ===================================================================
def ReadObsModelObs (ObsFile, ob_val):
  print ('### Reading the model observations')
  # ob_val is array of observation values read-in beforehand (when reading background obs)
  # This is input to check that the observations are same (i.e. that the analysis obs file is the correct one)
  Obs_file = open (ObsFile, 'r')
  line  = Obs_file.readline()
  line  = Obs_file.readline()
  split = line.split()
  nlevs = int(split[1])   # Needed to read-in the right number of lines of a profile at the start of the file
  for lev in range(nlevs+1):
    line  = Obs_file.readline()

  # Set-up the arrays of observation information
  ob_model     = []
  line         = Obs_file.readline()
  keep_reading = line != ''
  item         = -1

  while keep_reading:
    item        += 1
    line_ob_of   = Obs_file.readline()
    split        = line_ob_of.split()
    ob_of        = split[1]
    line_ob_type = Obs_file.readline()
    line_time    = Obs_file.readline()
    line_lon     = Obs_file.readline()
    line_lat     = Obs_file.readline()
    if ob_of == 't':
      line_lev   = Obs_file.readline()
    else:
      line_lev   = 0.0
    line_ob      = Obs_file.readline()
    line_modelob = Obs_file.readline()
    line_innov   = Obs_file.readline()
    line_grad    = Obs_file.readline()

    split = line_ob.split()
    # Check that the ob is the same as was read in before
    if np.abs(float(split[1]) - ob_val[item]) > 0.00000001:
      print ('Error: the model observation file is not compatible with the background observation file')
      print ('Observation', item)
      print (ObsFile)
      exit(0)
    split = line_modelob.split()
    ob_model.append(float(split[1]))

    # In preparation for the next observation
    line         = Obs_file.readline()
    keep_reading = line != ''

  Obs_file.close()
  return ob_model


# ===================================================================
def make_histogram_ob (nbins, obs, back_obs, anal_obs, true_obs, quantity, append_this, append_file, output_type, plot_dir, html_file):
  # Subroutine to generate histogram data of observations - background observations,
  #                                          observations - analysis observations,
  #                                          observations - true observations,
  #                                          background - true observations,
  #                                          analysis - true observations

  if (output_type == 'png'):
    filesuffix = '.png'
    html_file.write ('<h4>Observation differences for ' + quantity + '</h4>\n')
    html_file.write ('<table>\n')
    html_file.write ('<tr>\n')
  else:
    filesuffix = '.eps'

  # Number of observations
  Nobs = len(obs)

  # Convert the inputs to numpy arrays
  obs      = np.asarray(obs)
  back_obs = np.asarray(back_obs)
  anal_obs = np.asarray(anal_obs)
  true_obs = np.asarray(true_obs)

  # Compute the differences
  omb      = obs - back_obs
  oma      = obs - anal_obs
  omt      = obs - true_obs
  bmt      = back_obs - true_obs
  amt      = anal_obs - true_obs

  # Compute the means
  omb_mean = np.mean(omb)
  oma_mean = np.mean(oma)
  omt_mean = np.mean(omt)
  bmt_mean = np.mean(bmt)
  amt_mean = np.mean(amt)

  # Compute the standard deviations
  omb_std = np.std(omb)
  oma_std = np.std(oma)
  omt_std = np.std(omt)
  bmt_std = np.std(bmt)
  amt_std = np.std(amt)

  # Append these files if requested
  if (append_this):
    append_file.write (', ' + str(omb_mean) + ', ' + str(omb_std) + ', ' + str(oma_mean) + ', ' + str(oma_std))

  # Find the maximum absolute value of these differences
  maxval = 0.0
  maxval = max([maxval, max(abs(omb))])
  maxval = max([maxval, max(abs(oma))])
  maxval = max([maxval, max(abs(omt))])
  maxval = max([maxval, max(abs(bmt))])
  maxval = max([maxval, max(abs(amt))])
  maxval = flexi_round(maxval)

  # Generate the bins
  bins   = np.linspace(-1.0*maxval, maxval, nbins+1)
  print ('The following bins have been generated')
  print (bins)


  # Plot histograms

  # Observations - background
  fig, ax = plt.subplots()
  ax.hist(omb, bins=bins, histtype='bar', facecolor='b')
  ax.set_title('O-B for ' + quantity + '\nmean=' + str.format('%.5f' % omb_mean) + ', stddev=' + str.format('%.5f' % omb_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-b', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = quantity + '_o-b' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Observations - background<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')


  # Observations - analysis
  fig, ax = plt.subplots()
  ax.hist(oma, bins=bins, histtype='bar', facecolor='r')
  ax.set_title('O-A for ' + quantity + '\nmean=' + str.format('%.5f' % oma_mean) + ', stddev=' + str.format('%.5f' % oma_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-a', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = quantity + '_o-a' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Observations - analysis<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Observations - truth
  fig, ax = plt.subplots()
  ax.hist(omt, bins=bins, histtype='bar', facecolor='b')
  ax.set_title('O-T for ' + quantity + '\nmean=' + str.format('%.5f' % omt_mean) + ', stddev=' + str.format('%.5f' % omt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = quantity + '_o-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Observations - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Background - truth
  fig, ax = plt.subplots()
  ax.hist(bmt, bins=bins, histtype='bar', facecolor='b')
  ax.set_title('B-T for ' + quantity + '\nmean=' + str.format('%.5f' % bmt_mean) + ', stddev=' + str.format('%.5f' % bmt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('b-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = quantity + '_b-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Background - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Analysis - truth
  fig, ax = plt.subplots()
  ax.hist(amt, bins=bins, histtype='bar', facecolor='b')
  ax.set_title('A-T for ' + quantity + '\nmean=' + str.format('%.5f' % amt_mean) + ', stddev=' + str.format('%.5f' % amt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('a-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = quantity + '_a-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Analysis - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')


  if (output_type == 'png'):
    html_file.write ('</tr>\n')
    html_file.write ('</table>\n')


  return





# ===================================================================
def make_scatterplot_ob (obs, back_obs, anal_obs, true_obs, quantity, append_this, append_file, output_type, plot_dir, html_file):
  # Subroutine to do scatter plots of observations vs background observations,
  #                                   observations vs analysis observations,
  #                                   observations vs true observations,
  #                                   background vs true observations
  #                                   analysis vs true observations

  if (output_type == 'png'):
    filesuffix = '.png'
    html_file.write ('<h4>Observation scatter plots for ' + quantity + '</h4>\n')
    html_file.write ('<table>\n')
    html_file.write ('<tr>\n')
  else:
    filesuffix = '.eps'

  # Convert the inputs to numpy arrays
  obs      = np.asarray(obs)
  back_obs = np.asarray(back_obs)
  anal_obs = np.asarray(anal_obs)
  true_obs = np.asarray(true_obs)

  # Find the maximum and minimum values over all data
  maxval = max((obs))
  maxval = max([maxval, max(back_obs)])
  maxval = max([maxval, max(anal_obs)])
  maxval = max([maxval, max(true_obs)])

  minval = min((obs))
  minval = min([minval, min(back_obs)])
  minval = min([minval, min(anal_obs)])
  minval = min([minval, min(true_obs)])

  print ('Minimum value : ' + str(minval))
  print ('Maximum value : ' + str(maxval))


  # Do scatter plots

  # ===== Observations vs background =====
  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(obs[:], back_obs[:])
  fig, ax = plt.subplots()
  ax.scatter(obs[:], back_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color='cyan', label='x=y')
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color='red', label='best fit')
  ax.set_title(quantity + ' observations vs background', fontsize=16)
  ax.set_xlabel('observations', fontsize=16)
  ax.set_ylabel('background', fontsize=16)
  ax.legend()
  figfilename = quantity + '_o_vs_b' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Background vs observations<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300><br>\n')
    html_file.write ('m = ' + str(slope) + '<br>\n')
    html_file.write ('c = ' + str(intercept) + '<br>\n')
    html_file.write ('r = ' + str(r_value) + '<br>\n')
    html_file.write ('p = ' + str(p_value) + '\n')
  if (append_this):
    append_file.write (", " + str(slope) + ', ' + str(intercept))

  # Plot data with 1:1 line removed
  fig, ax = plt.subplots()
  ax.scatter(obs[:], back_obs[:]-obs[:])
  ax.set_xlabel('observations', fontsize=16)
  ax.set_ylabel('background - observations', fontsize=16)
  figfilename = quantity + '_o_vs_b-o' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300><br></td>\n')



  # ===== Observations vs analysis =====
  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(obs[:], anal_obs[:])
  fig, ax = plt.subplots()
  ax.scatter(obs[:], anal_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color='cyan', label='x=y')
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color='red', label='best fit')
  ax.set_title(quantity + ' observations vs analysis', fontsize=16)
  ax.set_xlabel('observations', fontsize=16)
  ax.set_ylabel('analysis', fontsize=16)
  ax.legend()
  figfilename = quantity + '_o_vs_a' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Analysis vs observations<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300><br>\n')
    html_file.write ('m = ' + str(slope) + '<br>\n')
    html_file.write ('c = ' + str(intercept) + '<br>\n')
    html_file.write ('r = ' + str(r_value) + '<br>\n')
    html_file.write ('p = ' + str(p_value) + '\n')
  if (append_this):
    append_file.write (", " + str(slope) + ', ' + str(intercept))

  # Plot data with 1:1 line subtracted
  fig, ax = plt.subplots()
  ax.scatter(obs[:], anal_obs[:]-obs[:])
  ax.set_xlabel('observations', fontsize=16)
  ax.set_ylabel('analysis - observations', fontsize=16)
  figfilename = quantity + '_o_vs_a-o' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300></td><br>\n')



  # ===== Observations vs truth =====
  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(obs[:], true_obs[:])
  fig, ax = plt.subplots()
  ax.scatter(obs[:], true_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color='cyan', label='x=y')
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color='red', label='best fit')
  ax.set_title(quantity + ' observations vs truth', fontsize=16)
  ax.set_xlabel('observations', fontsize=16)
  ax.set_ylabel('truth', fontsize=16)
  ax.legend()
  figfilename = quantity + '_o_vs_t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Observations vs truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300><br>\n')
    html_file.write ('m = ' + str(slope) + '<br>\n')
    html_file.write ('c = ' + str(intercept) + '<br>\n')
    html_file.write ('r = ' + str(r_value) + '<br>\n')
    html_file.write ('p = ' + str(p_value) + '\n')
  if (append_this):
    append_file.write (", " + str(slope) + ', ' + str(intercept))

  # Plot data with 1:1 line subtracted
  fig, ax = plt.subplots()
  ax.scatter(true_obs[:], anal_obs[:]-true_obs[:])
  ax.set_xlabel('truth', fontsize=16)
  ax.set_ylabel('analysis - truth', fontsize=16)
  figfilename = quantity + '_o_vs_a-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300></td><br>\n')





  # ===== Background vs truth =====
  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(back_obs[:], true_obs[:])
  fig, ax = plt.subplots()
  ax.scatter(back_obs[:], true_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color='cyan', label='x=y')
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color='red', label='best fit')
  ax.set_title(quantity + ' background vs truth', fontsize=16)
  ax.set_xlabel('background', fontsize=16)
  ax.set_ylabel('truth', fontsize=16)
  ax.legend()
  figfilename = quantity + '_b_vs_t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Background vs truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300><br>\n')
    html_file.write ('m = ' + str(slope) + '<br>\n')
    html_file.write ('c = ' + str(intercept) + '<br>\n')
    html_file.write ('r = ' + str(r_value) + '<br>\n')
    html_file.write ('p = ' + str(p_value) + '\n')
  if (append_this):
    append_file.write (", " + str(slope) + ', ' + str(intercept))

  # Plot data with 1:1 line subtracted
  fig, ax = plt.subplots()
  ax.scatter(true_obs[:], back_obs[:]-true_obs[:])
  ax.set_xlabel('truth', fontsize=16)
  ax.set_ylabel('background - truth', fontsize=16)
  figfilename = quantity + '_t_vs_b-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300></td><br>\n')




  # ===== Analysis vs truth =====
  # Find line of best fit
  slope, intercept, r_value, p_value, std_err = stats.linregress(anal_obs[:], true_obs[:])
  fig, ax = plt.subplots()
  ax.scatter(anal_obs[:], true_obs[:])
  x = [minval, maxval]
  y = [minval, maxval]
  ax.plot(x[:], y[:], linewidth=1, color='cyan', label='x=y')
  x = [minval, maxval]
  y = [slope*minval+intercept, slope*maxval+intercept]
  ax.plot(x[:], y[:], linewidth=1, color='red', label='best fit')
  ax.set_title(quantity + ' analysis vs truth', fontsize=16)
  ax.set_xlabel('analysis', fontsize=16)
  ax.set_ylabel('truth', fontsize=16)
  ax.legend()
  figfilename = quantity + '_a_vs_t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<td>Analysis vs truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300><br>\n')
    html_file.write ('m = ' + str(slope) + '<br>\n')
    html_file.write ('c = ' + str(intercept) + '<br>\n')
    html_file.write ('r = ' + str(r_value) + '<br>\n')
    html_file.write ('p = ' + str(p_value) + '\n')
  if (append_this):
    append_file.write (", " + str(slope) + ', ' + str(intercept))

  # Plot data with 1:1 line subtracted
  fig, ax = plt.subplots()
  ax.scatter(true_obs[:], anal_obs[:]-true_obs[:])
  ax.set_xlabel('truth', fontsize=16)
  ax.set_ylabel('analysis - truth', fontsize=16)
  figfilename = quantity + '_t_vs_a-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'png'):
    html_file.write ('<br><img src=' + figfilename  + ' width=300></td><br>\n')



  if (output_type == 'png'):
    html_file.write ('</tr>\n')
    html_file.write ('</table>\n')

  return




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



# ===================================================================
def hexstring2 (dec):
  # Take the decimal number (0 to 255) and convert into 2-digit hex string
  h = str(hex(int(dec)))
  if (len(h) == 4):
    hs = h[-2:]
  else:
    hs = '0' + h[-1:]
  return hs


# ===================================================================
def Err_in_Tracer (Estimate_tracer, Truth_tracer):
  # Compute the RMSE in the Estimate wrt the truth
  # Don't worry about grid weighting for now
  deviation = Estimate_tracer[:][:][:] - Truth_tracer[:][:][:]
  rmse      = np.sqrt(np.sum(deviation[:][:][:] * deviation[:][:][:]))
  return rmse


# ===================================================================
def Err_in_Flux (Estimate_tracer, Truth_tracer):
  # Compute the RMSE in the Estimate wrt the truth
  # Don't worry about grid weighting for now
  deviation = Estimate_tracer[:][:][:] - Truth_tracer[:][:][:]
  rmse      = np.sqrt(np.sum(deviation[:][:][:] * deviation[:][:][:]))
  return rmse


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
from os.path import exists


# How to run
# python3 Plot_PriorPostInc.py <Assim dir> <Backg dir> <Observations dir> <Truth dir> <Bg obs dir> <CVT file> <min tracer> <max tracer> <min flux> <max flux> <png or eps>

# Assim dir       : Directory containing the analysis state
# Backg dir       : Directory containing the background state
# Observations dir: Directory continaing the actial observations (used to get the true model observations only)
# Truth dir       : Directory containing the truth
# Bg obs dir      : Directory containing the background observations
# CVT file        : CVT file used by the DA
# min tracer      : Minimum value for tracer plot scale (0 if auto)
# min tracer      : Maximum value for tracer plot scale (0 if auto)
# min flux        : Minimum value for flux plot scale (0 if auto)
# min flux        : Maximum value for flux plot scale (0 if auto)
# png or eps      : Output graphics type

# Plots will be placed inside the <Assim dir>/plots
# If "png" is selected, the plots will be png and will be viewable in a web page
# If "eps" is selected, the plots will be eps

# E.g. python3 Plot_PriorPostInc.py /media/ross/banny/Enviflux/Bias/Bias0Assim/TracerFac_1.0__FluxFac_1.0_NEW /media/ross/banny/Enviflux/Bias/Background /media/ross/banny/Enviflux/Bias/Bias0Assim/UnbiasedObs /media/ross/banny/Enviflux/Bias/Truth /media/ross/banny/Enviflux/Bias/Bias0Assim/TracerFac_1.0__FluxFac_1.0_NEW /media/ross/banny/Enviflux/Bias/Covs/CVT__TracerFac_0.2__FluxFac_1.2.nc png
# E.g. python3 Plot_PriorPostInc.py /storage/research/nceo/ross/EnviFlux/Bias/Bias0Assim/TracerFac_1.0__FluxFac_1.0_NEW /storage/research/nceo/ross/EnviFlux/Bias/Background /storage/research/nceo/ross/EnviFlux/Bias/UnbiasedObs /storage/research/nceo/ross/EnviFlux/BiasTruth /storage/research/nceo/ross/EnviFlux/Bias/Bias0Assim/TracerFac_1.0__FluxFac_1.0_NEW /storage/research/nceo/ross/EnviFlux/Bias/Covs/CVT__TracerFac_0.2__FluxFac_1.2.nc png

num_args = len(sys.argv) - 1
if num_args != 11:
  print ("Incorrect number of arguments")
  print ("Expecting <Assim dir> <Backg dir> <Observations dir> <Truth dir> <Bg obs dir> <CVT file> <min tracer> <max tracer> <min flux> <max flux> <png or eps>")
  print ("Exiting")
  exit (0)

# Set directories
assim_dir        = sys.argv[1]
bg_dir           = sys.argv[2]
obs_dir          = sys.argv[3]
truth_dir        = sys.argv[4]
bg_obs_dir       = sys.argv[5]
plot_dir         = assim_dir + '/plots'
cvt_file         = sys.argv[6]
min_tracer_scale = float(sys.argv[7])
max_tracer_scale = float(sys.argv[8])
min_flux_scale   = float(sys.argv[9])
max_flux_scale   = float(sys.argv[10])
output_type      = sys.argv[11]

# Make the required directories
os.system('mkdir -p ' + plot_dir + '/graphics')

if (output_type != "png") and (output_type != "eps"):
  print ("The output type specified must be either png or eps")
  print ("Exiting")
  exit (0)


# ==============================================================
# \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/  \/
# USER OPTIONS

# Plot cost function?
plot_costfn = True

# Copy cost function info to html file (if png is selected)
copy_costfn_info = False


# Copy observations?
copy_ob_info = False

# Plot observation histograms?
plot_ob_hist = True
# Number of bins for histograms
nbins        = 40

# Plot observation scatter plots?
plot_ob_scatter = True

# Plot observation locations?
plot_ob_locations = True


# Plot flux fields?
plot_flux = True

# Plot background?
plot_tracer = True


# Plot tracer at vertical levels?
plot_vert_levs = True
# Choose vertical level indices to plot (start at 1)
vert_levs2plot = [1, 10, 55]

# Plot flux at times?
plot_flux_times = True
# Choose time indices to plot (start at 1)
times2plot     = [1]

# Plot tracer at longitudes?
plot_tracer_longs = True
# Choose longitude indices to plot (start at 1)
longs_2plot = [1, 15, 30, 45]

# /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\  /\
# ==============================================================

# Find the time now
now          = datetime.datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")

# Set-up html file for output
if output_type == 'png':
  filesuffix = '.png'
  html_file = open (plot_dir + '/plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write (current_time + '\n<br>\n')
  html_file.write ("Assim dir : " + assim_dir + '\n<br>\n')
  html_file.write ("Bg dir : " + bg_dir + '\n<br>\n')
  html_file.write ("Truth dir : " + truth_dir + '\n<br>\n')
else:
  filesuffix = '.eps'
  html_file  = 0



################################################################################
## CODE TO DO WITH COST FUNCTION ###############################################
################################################################################

if plot_costfn or copy_costfn_info:
  # Read the cost function with iteration
  print ('### Reading the cost function')
  J_file = open (assim_dir + '/AnalDiags.dat', 'r')

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
  while (line != ''):
    split = line.split()
    its.append (float(split[0]))
    Jb.append (float(split[1]))
    Jo.append (float(split[2]))
    J.append (float(split[3]))
    Jgrad.append (float(split[5]))
    line  = J_file.readline()

  J_file.close()

  #print (J)

if plot_costfn:
  # Plot the cost function
  PlotJ (J[:], Jb[:], Jo[:], its[:], html_file, plot_dir, output_type)


if copy_costfn_info and (output_type == 'png'):
  # Copy the cost function information to the web page
  html_file.write ('<h2>Cost function data</h2>\n')
  html_file.write ('<table cols=5 cellspacing=5 cellpadding=3 border=1>\n')
  html_file.write ('<tr><td>iter</td><td>Jb</td><td>Jo</td>' +
                   '<td>J</td><td>J grad</td></tr>\n')
  for l in range(len(J)):
    html_file.write ('<tr><td>' + str(l) + '</td><td>' + str(Jb[l]) + '</td><td>' + str(Jo[l]) + '</td>' +
                     '<td>' + str(J[l]) + '</td><td>' + str(Jgrad[l]) + '</td></tr>\n')
  html_file.write ('</table>\n')





################################################################################
## CODE TO DO WITH OBSERVATIONS ################################################
################################################################################

if plot_ob_hist or plot_ob_scatter or plot_ob_locations or copy_ob_info:
  # Read the observation information
  # Read observations, background observations, and meta data
  ob_time_s, ob_lon, ob_lat, ob_lev, ob_val, ob_stddev, ob_bg = ReadObsBg (bg_obs_dir + '/' + 'ObsBg.dat')

  # Read analysis observations
  print ('Reading analysis model observations')
  ob_anal   = ReadObsModelObs (assim_dir + '/' + 'ObsAnal.dat', ob_val)


  # Read true values at observation locations (these are the model obs in the actual observation files)
  print ('Reading true observation values')
  ob_true   = ReadObsModelObs (obs_dir + '/' + 'Observations.dat', ob_val)

  #for item in range(len(ob_bg)):
  #  print (item, ob_time_s[item], ob_lon[item], ob_lat[item], ob_lev[item], ob_val[item], ob_stddev[item], ob_bg[item], ob_anal[item])

  nobs  = len(ob_val)

  print ('Number of observations read in: ', nobs)


  if (copy_ob_info and (output_type == 'png')):

    html_file.write ('<h2>Observation list</h2>\n')
    html_file.write ('Number of observations = ' + str(nobs) + '<br>\n')
    html_file.write ('<table cols=8 cellspacing=5 cellpadding=3 border=1>\n')
    html_file.write ('<tr><td>lon (deg)</td><td>lat (deg)</td><td>lev (m)</td>' +
                     '<td>time (s)</td><td>obs</td><td>sigma</td><td>prior</td>' + 
                     '<td>post</td></tr>\n')
    for l in range(nobs):
      html_file.write ('<tr>')
      html_file.write ('<td>' + str(ob_lon[l]) + '</td>' +
                       '<td>' + str(ob_lat[l]) + '</td>' +
                       '<td>' + str(ob_lev[l]) + '</td>' +
                       '<td>' + str(ob_time_s[l]) + '</td>' +
                       '<td>' + str(ob_val[l]) + '</td>' +
                       '<td>' + str(ob_stddev[l]) + '</td>' +
                       '<td>' + str(ob_bg[l]) + '</td>' +
                       '<td>' + str(ob_anal[l]) + '</td>')
      html_file.write ('</tr>\n')
    html_file.write ('</table>\n')



  if (plot_ob_hist):
    # Generate histograms
    make_histogram_ob (nbins,
                       ob_val[0:nobs], ob_bg[0:nobs], ob_anal[0:nobs], ob_true[0:nobs],
                       'Tracer_obs', False, 0, output_type, plot_dir, html_file)


    # Generate histograms for normalised observation data
    #make_histogram_ob (nbins, ob_val[0:nobs]/ob_stddev[0:nobs],
    #                          ob_bg[0:nobs]/ob_stddev[0:nobs],
    #                          ob_anal[0:nobs]/ob_stddev[0:nobs],
    #                   'Tracer_obs_norm', output_type, plot_dir, html_file)


  if (plot_ob_scatter):
    # Generate scatter plots for the raw observation data
    make_scatterplot_ob (ob_val[0:nobs],
                         ob_bg[0:nobs],
                         ob_anal[0:nobs],
                         ob_true[0:nobs],
                         'Tracer_obs', False, 0, output_type, plot_dir, html_file)


  if (plot_ob_locations):
    plot_ob_hts (ob_lon, ob_lat, ob_lev, 'tracer', output_type, plot_dir, html_file)






################################################################################
## CODE TO DO WITH FIELDS ######################################################
################################################################################


# The data come from netcdf data file (including the longitudes, latitudes, and heights)

if plot_flux or plot_tracer:

  # Truth (and get longitudes, latiudes, levels, and times)
  nc_file = Dataset(truth_dir + '/Truth.nc')
  longs = nc_file.variables['longitude'][:]
  #print ('Longitudes : ', longs[:])
  lats  = nc_file.variables['latitude'][:]
  #print ('Latitudes : ', lats[:])
  levs  = nc_file.variables['level'][:] / 1000.0
  #print ('Levels : ', levs[:])
  times = nc_file.variables['time'][:]
  # Change from seconds to 30-day periods
  times = times[:] / 2592000.0
  #print ('Times : ', times[:])
  if plot_tracer:
    Field_tracer_truth = nc_file.variables['tracer0'][:][:][:]
    print ('Shape of tracer')
    print (Field_tracer_truth.shape)
  if plot_flux:
    Field_flux_truth   = nc_file.variables['source'][:][:][:]
    print ('Shape of flux')
    print (Field_flux_truth.shape)
  nc_file.close()

  # Background
  nc_file = Dataset(bg_dir + '/Background.nc')
  if plot_tracer:
    Field_tracer_bg = nc_file.variables['tracer0'][:][:][:]
  if plot_flux:
    Field_flux_bg   = nc_file.variables['source'][:][:][:]
  nc_file.close()

  # Analysis and analysis increment
  nc_file = Dataset(assim_dir + '/Anal.nc')
  if plot_tracer:
    Field_tracer_anal    = nc_file.variables['tracer0'][:][:][:]
  if plot_flux:
    Field_flux_anal      = nc_file.variables['source'][:][:][:]
  nc_file.close()

  nc_file = Dataset(assim_dir + '/AnalInc.nc')
  if plot_tracer:
    Field_tracer_analinc = nc_file.variables['tracer0'][:][:][:]
  if plot_flux:
    Field_flux_analinc   = nc_file.variables['source'][:][:][:]
  nc_file.close()

  # Background error standard deviations
  nc_file = Dataset(cvt_file)
  x_tomcat         = nc_file.variables['x_tomcat'][:]
  y_tomcat         = nc_file.variables['y_tomcat'][:]
  if plot_tracer:
    Field_tracer_std = nc_file.variables['std_tracer'][:][:][:]
    print ('Shape of tracer standard deviation')
    print (Field_tracer_std.shape)
  if plot_flux:
    Field_flux_std   = nc_file.variables['std_flux'][:][:][:][:]
    print ('Shape of flux standard deviation')
    print (Field_flux_std.shape)
  nc_file.close()


  # ====== Compute errors in the background and analysis ======
  if plot_tracer:
    RMSE_tracer_bg   = Err_in_Tracer (Field_tracer_bg[:][:][:], Field_tracer_truth[:][:][:])
    RMSE_tracer_anal = Err_in_Tracer (Field_tracer_anal[:][:][:], Field_tracer_truth[:][:][:])
    if output_type == 'png':
      html_file.write ('RMSE in the background tracer : ' + str(RMSE_tracer_bg) + '\n<br>\n')
      html_file.write ('RMSE in the analysis tracer : ' + str(RMSE_tracer_anal) + '\n<br>\n')
      html_file.write ('Reduction in RMSE : ' + str(RMSE_tracer_bg - RMSE_tracer_anal) + '\n<br>\n')
      html_file.write ('As a percentage : ' + str(100.0 * (RMSE_tracer_bg - RMSE_tracer_anal) / RMSE_tracer_bg) + '\n<br>\n')
    else:
      print ('RMSE in the background tracer : ', RMSE_tracer_bg)
      print ('RMSE in the analysis tracer   : ', RMSE_tracer_anal)
      print ('Reduction in RMSE             : ', RMSE_tracer_bg - RMSE_tracer_anal)
      print ('As a percentage               : ', 100.0 * (RMSE_tracer_bg - RMSE_tracer_anal) / RMSE_tracer_bg)

  if plot_flux:
    RMSE_flux_bg   = Err_in_Flux (Field_flux_bg[:][:][:], Field_flux_truth[:][:][:])
    RMSE_flux_anal = Err_in_Flux (Field_flux_anal[:][:][:], Field_flux_truth[:][:][:])
    if output_type == 'png':
      html_file.write ('RMSE in the background flux : ' + str(RMSE_tracer_bg) + '\n<br>\n')
      html_file.write ('RMSE in the analysis flux : ' + str(RMSE_tracer_anal) + '\n<br>\n')
      html_file.write ('Reduction in RMSE : ' + str(RMSE_flux_bg - RMSE_flux_anal) + '\n<br>\n')
      html_file.write ('As a percentage : ' + str(100.0 * (RMSE_flux_bg - RMSE_flux_anal) / RMSE_flux_bg) + '\n<br>\n')
    else:
      print ('RMSE in the background flux   : ', RMSE_flux_bg)
      print ('RMSE in the analysis flux     : ', RMSE_flux_anal)
      print ('Reduction in RMSE             : ', RMSE_flux_bg - RMSE_flux_anal)
      print ('As a percentage               : ', 100.0 * (RMSE_flux_bg - RMSE_flux_anal) / RMSE_flux_bg)


  if plot_flux:
    # Infer the land-sea mask from the background error data
    nlon    = len(longs)
    nlat    = len(lats)
    landsea = np.zeros([nlon,nlat], dtype=int)
    av      = np.mean(Field_flux_std[0,0,:,:])
    for yy in range(nlat):
      if yy == nlat-1:
        y = nlat-2
      else:
        y = yy
      for x in range(nlon):
        if (x == 0) or (x == nlon-1):
          x1 = nlon-2
          x2 = 0
        else:
          x1 = x - 1
          x2 = x
        if ((Field_flux_std[0][0][x1][y] + Field_flux_std[0][0][x2][y])/2.0 > av) or (y == nlat-2) or (y == nlat-1):
          landsea[x][y] = 1
        #print (landsea[x][y], end='')
      #print ()



  # Plot the flux data for each field type and month
  if plot_flux and plot_flux_times:
    if output_type == 'png':
      html_file.write ('<table>\n')

    for m in times2plot:

      # ===== Truth =====
      if output_type == 'png':
        html_file.write ('<tr><td>\n')

      # Calculate the total, land, and sea fluxes for the truth
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Field_flux_truth[m-1,:,:], longs[:], lats[:], landsea[:,:],
                                                     'truth_m'+str(m), html_file, output_type)
      print ('Truth total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the truth
      PlotMap (Field_flux_truth[m-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux truth month ' + str(m),
              'flux_truth_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Background =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Calculate the total, land, and sea fluxes for the background
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Field_flux_bg[m-1,:,:], longs[:], lats[:], landsea[:,:],
                                                     'prior_m'+str(m), html_file, output_type)
      print ('Background total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the background
      PlotMap (Field_flux_bg[m-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux prior month ' + str(m),
              'flux_prior_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Analysis =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Calculate the total, land, and sea fluxes for the analysis
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Field_flux_anal[m-1,:,:], longs[:], lats[:], landsea[:,:],
                                                     'analysis_m'+str(m), html_file, output_type)
      print ('Analysis total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the analysis
      PlotMap (Field_flux_anal[m-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux anal month ' + str(m),
              'flux_anal_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Analysis increment =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Calculate the total, land, and sea fluxes for the analysis increment
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Field_flux_analinc[m-1,:,:], longs[:], lats[:], landsea[:,:],
                                                     'analysis_inc_m'+str(m), html_file, output_type)
      print ('Analysis increment total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the analysis increment
      PlotMap (Field_flux_analinc[m-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux anal inc month ' + str(m),
              'flux_anal_inc_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Background error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Calculate the total, land, and sea fluxes for the background error
      Err = Field_flux_bg[m-1,:,:] - Field_flux_truth[m-1,:,:]
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Err[:,:], longs[:], lats[:], landsea[:,:],
                                                     'bgerr_m'+str(m), html_file, output_type)
      print ('Background error total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the background error
      PlotMap (Err[:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux prior error month ' + str(m),
              'flux_prior_err_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Analysis error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Calculate the total, land, and sea fluxes for the analysis error
      Err = Field_flux_anal[m-1,:,:] - Field_flux_truth[m-1,:,:]
      flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg = CalcTotalFlux (Err[:,:], longs[:], lats[:], landsea[:,:],
                                                     'analerr_m'+str(m), html_file, output_type)
      print ('Analysis error total, land, sea, total source, total sink fluxes for month', m, ':')
      print (flux_tot, flux_land, flux_sea, flux_tot_pos, flux_tot_neg)

      # Plot the flux map for the analysis error
      PlotMap (Err[:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Flux anal error month ' + str(m),
              'flux_anal_err_m' + str(m),
              html_file, plot_dir, output_type,
              minplot=min_flux_scale, maxplot=max_flux_scale, symm=True)


      # ===== Background error standard deviation =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')

      # Prior flux error stddev
      PlotMap (np.transpose(Field_flux_std[0,m-1,:,:]), x_tomcat[:], y_tomcat[:], 'longitude', 'latitude',
               'flux prior err std month ' + str(m),
               'flux_err_std_month' + str(m),
               html_file, plot_dir,
               output_type, minplot=min_flux_scale, maxplot=max_flux_scale)




      if (output_type == 'png'):
        html_file.write ('</td></tr>\n<tr>\n')

    if output_type == 'png':
      html_file.write ('</table>\n')





  # Plot the tracer data for each chosen level
  if plot_tracer and plot_vert_levs:
    if output_type == 'png':
      html_file.write ('<table>\n')

    for l in vert_levs2plot:

      # ===== Truth =====
      if output_type == 'png':
        html_file.write ('<tr><td>\n')
      PlotMap (Field_tracer_truth[l-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer truth level ' + str(l),
              'tracer_truth_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Background =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      PlotMap (Field_tracer_bg[l-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer prior level ' + str(l),
              'tracer_prior_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Analysis =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      PlotMap (Field_tracer_anal[l-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer anal level ' + str(l),
              'tracer_anal_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Analysis increment =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      PlotMap (Field_tracer_analinc[l-1,:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer anal inc level ' + str(l),
              'tracer_anal_inc_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale, symm=True)
              
      # =====  Background error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Err = Field_tracer_bg[l-1,:,:] - Field_tracer_truth[l-1,:,:]
      PlotMap (Err[:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer prior error level ' + str(l),
              'tracer_prior_err_inc_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale, symm=True)

      # =====  Analysis error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Err = Field_tracer_anal[l-1,:,:] - Field_tracer_truth[l-1,:,:]
      PlotMap (Err[:,:], longs[:], lats[:], 'longitude', 'latitude',
              'Tracer anal error level ' + str(l),
              'tracer_anal_err_inc_l' + str(l),
              html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale, symm=True)

      # ===== Background error standard deviation =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      PlotMap (np.transpose(Field_tracer_std[l-1,:,:]), x_tomcat[:], y_tomcat[:], 'longitude', 'latitude',
               'tracer prior err std level ' + str(l),
               'tracer_err_std_l' + str(l),
               html_file, plot_dir, output_type,
               minplot=min_tracer_scale, maxplot=max_tracer_scale)



      if output_type == 'png':
        html_file.write ('</td></tr>\n<tr>\n')

    if output_type == 'png':
      html_file.write ('</table>\n')



  # Plot the tracer data for each chosen longitude
  if plot_tracer and plot_tracer_longs:

    if output_type == 'png':
      html_file.write ('<table>\n')

    for l in longs_2plot:

      # ===== Truth =====
      if output_type == 'png':
        html_file.write ('<tr><td>\n')
      Plot2d (Field_tracer_truth[:,:,l], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer truth lon ' + str(longs[l]),
              'tracer_truth_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Background =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Plot2d (Field_tracer_bg[:,:,l], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer prior lon ' + str(longs[l]),
              'tracer_prior_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Analysis =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Plot2d (Field_tracer_anal[:,:,l], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer anal lon ' + str(longs[l]),
              'tracer_anal_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # ===== Analysis increment =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Plot2d (Field_tracer_analinc[:,:,l], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer anal inc lon ' + str(longs[l]),
              'tracer_anal_inc_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)

      # =====  Background error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Err = Field_tracer_bg[:,:,l] - Field_tracer_truth[:,:,l]
      Plot2d (Err[:,:], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer prior error lon ' + str(longs[l]),
              'tracer_prior_err_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale, symm=True)

      # =====  Analysis error =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Err = Field_tracer_anal[:,:,l] - Field_tracer_truth[:,:,l]
      Plot2d (Err[:,:], lats[:], levs[:], 'latitude', 'level (km)',
              'Tracer anal error lon ' + str(longs[l]),
              'tracer_anal_err_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale, symm=True)

      # ===== Background error standard deviation =====
      if output_type == 'png':
        html_file.write ('</td><td>\n')
      Plot2d (Field_tracer_std[:,l,:], y_tomcat[:], levs[:], 'latitude', 'level (km)',
              'Tracer prior err std lon ' + str(longs[l]),
              'tracer_err_std_lon' + str(l),
              False, html_file, plot_dir, output_type,
              minplot=min_tracer_scale, maxplot=max_tracer_scale)


      if output_type == 'png':
        html_file.write ('</td></tr>\n<tr>\n')

    if output_type == 'png':
      html_file.write ('</table>\n')


if (output_type == 'png'):
  html_file.write ('</html>')
  html_file.close()
  print ('An html file has been created to view the figures')
  print ('Please view the following html file with your browser')
  print (plot_dir + '/plots.html')
  print ('')
