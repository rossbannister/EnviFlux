def inlist (sequence, value):
  # Function returns True if the list 'sequence' contains 'value'
  # False otherwise
  present = False
  for item in sequence:
    if (item == value):
      present = True
  return present


def posinlist (sequence, value):
  # Function to return the position of 'vale' in list 'sequence'
  counter = 0
  pos     = -1
  for item in sequence:
    if (item == value):
      pos      = counter
    counter += 1
  return pos


################################################################
### File to study stability of advection schemes
### This code visualises the output of Master_RunForecast.cpp
################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
import os


# Set input file containing maximum maximum values of tracer
max_file = '../TestRuns/StabilityAnalysis3a'
plotdir  = '../TestRuns/Plots3a'




# Values of independent variables
dt      = []
kappa_h = []
kappa_v = []


################################################################
# First read through to extract dt, kappa_h, and kappa_v values
input_file = open (max_file, 'r')
proceed = True
counter = 0

while (proceed):
  lod = input_file.readline()
  if ((lod == '') or (lod == '\n')):
    proceed = False
  else:
    counter += 1
    print (counter)
    data     = lod[13:-1]
    dt_val   = float(data)
    #print (dt_val, ':', dt)
    if (not(inlist(dt, dt_val))):
      dt.append(dt_val)

    lod        = input_file.readline()
    #print ('lod (kappa_h): ', lod)
    data       = lod[13:-1]
    kappah_val = float(data)
    if (not(inlist (kappa_h, kappah_val))):
      kappa_h.append(kappah_val)

    lod        = input_file.readline()
    #print ('lod (kappa_v): ', lod)
    data       = lod[13:-1]
    kappav_val = float(data)
    if (not(inlist (kappa_v, kappav_val))):
      kappa_v.append(kappav_val)

    lod = input_file.readline()
    lod = input_file.readline()
    lod = input_file.readline()
    lod = input_file.readline()

input_file.close()

print ('There are ', counter, ' entries')

print ('Here are the dt values extracted')
print (dt)

print ('Here are the kappa_h values extracted')
print (kappa_h)

print ('Here are the kappa_v values extracted')
print (kappa_v)

ndt      = len(dt)
nkappa_h = len(kappa_h)
nkappa_v = len(kappa_v)


################################################################
# Now go round and collect the maximum tracer value data
# (large values will be due to instability of the advection scheme)

maxval = np.zeros((ndt, nkappa_h, nkappa_v))
input_file = open (max_file, 'r')
for c in range(counter):
  print (c)
  lod        = input_file.readline()
  data       = lod[13:-1]
  dt_val     = float(data)
  #print (dt_val, ':', dt)
  pos_dt     = posinlist (dt, dt_val)
  lod        = input_file.readline()
  #print ('lod (kappa_h): ', lod)
  data       = lod[13:-1]
  kappah_val = float(data)
  pos_kappah = posinlist (kappa_h, kappah_val)
  lod        = input_file.readline()
  #print ('lod (kappa_v): ', lod)
  data       = lod[13:-1]
  kappav_val = float(data)
  pos_kappav = posinlist (kappa_v, kappav_val)
  lod        = input_file.readline()
  lod        = input_file.readline()
  lod        = input_file.readline()
  data       = lod[13:-1]
  if (data[0:2] == 'inf'):
    maxval[pos_dt, pos_kappah, pos_kappav] = np.nan
  else:
    val = float(data)
    if (val > 1000000.0):
      maxval[pos_dt, pos_kappah, pos_kappav] = np.nan
    else:
      maxval[pos_dt, pos_kappah, pos_kappav] = val

  lod = input_file.readline()
  print (dt_val, kappah_val, kappav_val, ':', pos_dt, pos_kappah, pos_kappav, ':', maxval[pos_dt, pos_kappah, pos_kappav])

input_file.close()


################################################################
# Systematically plot the results

os.system('mkdir -p ' + plotdir)
html_file = open (plotdir + '/Stability.html', 'w')
html_file.write ('<html>\n')


# Output the values of dt, kappa_h, and kappa_v
maxnoparams = ndt
if (nkappa_h > maxnoparams):
  maxnoparams = nkappa_h
if (nkappa_v > maxnoparams):
  maxnoparams = nkappa_v
html_file.write ('<table cols=4 border=1>\n')
html_file.write ('<tr><td>parameter</td><td>dt</td><td>kappa_h</td><td>kappa_v</td></tr>')
for p in range(maxnoparams):
  html_file.write ('<tr><td>' + str(p) + '</td>')
  if (ndt >= p+1):
    html_file.write ('<td>' + str(dt[p]) + '</td>')
  else:
    html_file.write ('<td>.</td>')
  if (nkappa_h >= p+1):
    html_file.write ('<td>' + str(kappa_h[p]) + '</td>')
  else:
    html_file.write ('<td>.</td>')
  if (nkappa_v >= p+1):
    html_file.write ('<td>' + str(kappa_v[p]) + '</td>')
  else:
    html_file.write ('<td>.</td>')
  html_file.write ('</tr>\n')
html_file.write ('</table>\n')


# 1. Slices through dt
html_file.write ('Different values of dt')
html_file.write ('<table cols=' + str(ndt) + '>')
html_file.write ('<tr>\n')

for number in range(ndt):
  fig, ax = plt.subplots()
  cax     = ax.imshow(maxval[number,:,:], interpolation='None', cmap='rainbow')
  cbar    = fig.colorbar(cax, orientation='vertical')
  ax.set_title('Maxvals for dt=' + str(dt[number]))
  ax.set_xlabel('kappa v', fontsize=16)
  ax.set_ylabel('kappa h', fontsize=16)
  #plt.show()
  filename = 'dt_%04i.png'% (number)
  plt.savefig(plotdir + '/' + filename, bbox_inches='tight')
  plt.close('all')
  html_file.write ('<td><img src=' + filename + ' width=350></td>')
html_file.write ('</tr>\n')
html_file.write ('</table>\n')



# 2. Slices through kappa_h
html_file.write ('Different values of kappa_h')
html_file.write ('<table cols=' + str(nkappa_h) + '>')
html_file.write ('<tr>\n')

for number in range(nkappa_h):
  fig, ax = plt.subplots()
  cax     = ax.imshow(maxval[:,number,:], interpolation='None', cmap='rainbow')
  cbar    = fig.colorbar(cax, orientation='vertical')
  ax.set_title('Maxvals for kappah=' + str(kappa_h[number]))
  ax.set_xlabel('kappa_v', fontsize=16)
  ax.set_ylabel('dt', fontsize=16)
  #plt.show()
  filename = 'kappah_%04i.png'% (number)
  plt.savefig(plotdir + '/' + filename, bbox_inches='tight')
  plt.close('all')
  html_file.write ('<td><img src=' + filename + ' width=350></td>')
html_file.write ('</tr>\n')
html_file.write ('</table>\n')


# 3. Slices through kappa_v
html_file.write ('Different values of kappa_v')
html_file.write ('<table cols=' + str(nkappa_v) + '>')
html_file.write ('<tr>\n')

for number in range(nkappa_v):
  fig, ax = plt.subplots()
  cax     = ax.imshow(maxval[:,:,number], interpolation='None', cmap='rainbow')
  cbar    = fig.colorbar(cax, orientation='vertical')
  ax.set_title('Maxvals for kappav=' + str(kappa_v[number]))
  ax.set_xlabel('kappa_h', fontsize=16)
  ax.set_ylabel('dt', fontsize=16)
  #plt.show()
  filename = 'kappav_%04i.png'% (number)
  plt.savefig(plotdir + '/' + filename, bbox_inches='tight')
  plt.close('all')
  html_file.write ('<td><img src=' + filename + ' width=350></td>')
html_file.write ('</tr>\n')
html_file.write ('</table>\n')

html_file.write ('</html>\n')
html_file.close ()
