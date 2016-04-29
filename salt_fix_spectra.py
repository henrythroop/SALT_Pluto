# This routine will read in the fully processed and merged / combined SALT spectra.
# It will then do science on them: measure bands, do smoothing, etc.
#
#
# HBT 29-Oct-2014

import pdb
import sys # for sys.stdout.write

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
# import astropy
# from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 

import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *
from   astropy.convolution import convolve, Box1DKernel

#########
# Function for wheremin()
#########

def wheremin( arr ):
   "Determines the index at which an array has its minimum value"
   index = where(arr == amin(arr))
   return index[0]
   
##########
# Function to align_wavelengths()
##########
   
def align_wavelengths( wavelengths, fluxes, wavelength_range, wavelength_center, smoothing=1):
    "Take a bunch of spectra, and slide them right/left until they align"
    
# All spectra share the same identical wavelength axis
# Procedure: 
#   Loop over all wavelengths
#   Extract a portion between min and max
#   If requested, smooth it by 'smoothing' with boxcar    
#   In that spectrum, calc the minimum value between min and max set in wavelength_range
#   Shift it so that position of minimum is set to the wavelength desired
#   ** During the shift, no data are lost. They roll from one end to the other.
#    
# NB: Python parameters are all passed by reference. So the original array gets changed, as intended.   

    bin_start = wheremin(abs(wavelengths - wavelength_range[0]))[0] # Index for starting wavelength in region 
    bin_end   = wheremin(abs(wavelengths - wavelength_range[1]))[0] # Index for ending wavelenngths
    wavel_extract = wavelengths[bin_start : bin_end]      
    index_right = wheremin(abs(wavel_extract - wavelength_center))[0] # Get the proper position

    for i in range(len(fluxes)):   
        fluxes_extract = fluxes[i][bin_start : bin_end].copy()               # Extract the small region to examine
        fluxes_extract_smooth = convolve(fluxes_extract, Box1DKernel(smoothing))
#        fluxes_extract_smooth = roll(fluxes_extract_smooth, -(smoothing-1)/2) # Weird. Sometimes the output array is dfft size from input.
        hw = (smoothing-1)/2 + 1 # halfwidth of the smoothing
        index_wrong = wheremin(fluxes_extract_smooth[hw:-hw])[0] + hw # Get the central position 
        droll = index_right - index_wrong
#        print "index_right = {0}, index_wrong = {1}".format(index_right, index_wrong)
#        print "droll = " + repr(droll)
#        quit       
#        print 'Rolling spectrum #' + repr(i) + ' by ' + repr(droll) + " to match " + repr(wavelength_center)
        fluxes[i] = roll(fluxes[i], droll)
   
    return # No return arguments -- everything is by reference
        
##########
# Start of main program
##########
   
# Set the directory for the data and all analysis. 

dir_data = "/Users/throop/Data/SALT_Pluto_Jun14/product/" 

files_hd = glob.glob(dir_data + '/spect_hd*txt')
files_pluto = glob.glob(dir_data + '/spect_pluto*txt')

# Each  

flux_hd     = []
flux_pluto  = []
wavelength    = [] # I think that all of the files should have an identical wavelength now.
date_pluto    = []
date_hd       = []

# *** Should look up 'genfromtext' or 'getfromcsv' in numpy
# Actually, SciPy FAQ says that numpy.loadtxt() is probably the way to do it.
# It will deal with header lines, CSV, etc.
# Astropy also has a 'unified table reader'

# Read all of the data into a list of arrays.
# It might be possible to do this as an array of arrays -- not sure.
# See 'advanced topics' section in NumPy book.
# File is a simple three-column list of numbers, separated by spaces

for file in files_hd: 
    print "Loading file " + file
    d = loadtxt(file, delimiter=',')
    flux_hd.append(d[:,1])
    l = len(file)
    date_hd.append(file[l-14:l-4])
    wavelength.append(d[:,0])

wavelength = wavelength[0]
    
for file in files_pluto: 
    print "Loading file " + file
    d = loadtxt(file, delimiter=',')
    flux_pluto.append(d[:,1])
    l = len(file)
    date_pluto.append(file[l-14:l-4])

# Set up plotting parameters

rcParams['figure.figsize'] = 20, 10
colors = ['red', 'blue', 'green', 'purple', 'orange', 'black', 'darkblue', 'salmon', 'pink', 'brown', 'olive', 'violet', 'brown', 'tan']

# Make a plot of all HD obs, so we can pull out the best

xlim = (3500, 9500)
ylim = ((0, 1.2))
index_good_spectra_hd = {0, 1, 2, 3, 4} # 3 looks the best (26-Aug)
for i in range(len(files_hd)):
  if i in index_good_spectra_hd:
      plot(wavelength, flux_hd[i] + i*0.1, color=colors[i])
      plt.text(amax(wavelength), 
               np.median(flux_hd[i][-100:-1]) + i*0.1, 
               '  ' + date_hd[i])

plt.title('HD, SALT June-August 2014', fontsize=24)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

show()

##########
# Make a plot of region around 0.73 um band
##########

# ** This is the important one. Check RSS Simulator to see where it puts the gap.
# PG0900, actual tilt in FITS files = 20.375
# A: RSS simulator says gap should be at 7150 .. 7200 (for tilt 20.375)
# Data put the gap at                    7175 .. 7225 (for tilt 20.375)
# So it's off by about a half a gap width. That is quite unfortunate. But it's more my fault than SALT's.

# If I would have requested a tilt angle of 40 deg (20 deg), then the gap would be clearly off of the CH4 band, at 7000 .. 7050.
# That is clearly what I intended from my proposal. Fuck.
 
xlim = (6700, 7500)

binning = 10

for i in range(len(files_pluto)):
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.1, Box1DKernel(width=binning)), color=colors[i])
      plt.text(amax(wavelength), 
               np.median(flux_pluto[i][-100:-1]) + i*0.1, 
               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT. Band = CH$_4$, 0.71 $\mu$m - 0.74 $\mu$m', fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim((0.7,2.5))
show()

##########
# Make plot around solar Ha (H alpha) line at 6563 A
##########

xlim = (6520, 6620)
ylim = (0.4, 3.0)
for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
plt.title("Pluto (bottom) and HD (top). Band H$\\alpha$ = 6563 $\AA$", fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make plot around Cruikshank's 0.62 um CH4 band
##########

xlim = (6100, 6300)
ylim = (0.4, 3.0)
for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
plt.title("Pluto (bottom) and HD (top). Band = 6200 $\AA$ CH$_4$", fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Now center on Ha, and plot again
##########

#wavelengths = wavelength
#fluxes = flux_pluto
#wavelength_range = (6540, 6580)
#wavelength_center = 6563
#smoothing = 10
#
#stop

align_wavelengths(wavelength, flux_pluto, (6540, 6580), 6563, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (6540, 6580), 6563, smoothing=10)

xlim = (6520, 6620)
ylim = (0.4, 3.0)
for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
plt.title("Pluto (bottom) and HD (top). Band H$\\alpha$ = 6563 $\AA$, aligned @ 6563", fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

#stop

##########
# Make plot around Telluric O2 line at 7620 A. See Grundy 1996 @ 332. Used for wavelength calibration.
##########

align_wavelengths(wavelength, flux_pluto, (7580, 7620), 7600, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (7580, 7620), 7600, smoothing=10)

xlim = (7570, 7670)
ylim = (0.2, 2.6)
binning = 1

for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
plt.title('Pluto (bottom) and HD (top). Band O$_2$ = 7620 $\AA$, aligned @ 7600', fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim(ylim)
show()
     
#stop
      
      
##########
# Make plot around 0.89 um. HD and Pluto separately. Used for wavelength calibration.
##########

align_wavelengths(wavelength, flux_pluto, (8530, 8550), 8540, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (8530, 8550), 8540, smoothing=10)

xlim = (8300, 9300)

ylim = (0.2, 2.6)
binning = 1

for i in range(len(files_pluto)):
      plot(wavelength, flux_pluto[i] + i*0.1, color=colors[i])
      
for i in range(len(files_hd)):
      plot(wavelength, flux_hd[i] + i*0.1 + 1.5, color=colors[i])      
      
plt.title('Pluto (bottom) and HD (top). Band O$_2$ = 7620 $\AA$, aligned @ 8540', fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim(ylim)
show()


##########
# Make a plot of region around 0.89 um band. Ratio.
##########
              
xlim = (8300, 9300)
ylim = (1.0, 4)
binning = 5

index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
#      plt.text(amax(wavelength), 
#               np.median(flux_pluto[i][-200:-50]) + i*0.2, 
#               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT, 8500 - 8700 and 8700 - 8900', fontsize=24)
plt.xlim(xlim)
plt.ylim((ylim))

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make a plot of region around 0.89 um band, BINNED
##########
align_wavelengths(wavelength, flux_pluto, (8530, 8550), 8540, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (8530, 8550), 8540, smoothing=10)        

xlim = (8300, 9300)
ylim = (0.9, 4)
binning = 10

for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT, 8500-8700 and 8700-8900, binned x ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)

file_out = 'spect_ratio_8900_bin{0}.png'.format(binning)
plt.savefig(file_out)
print "Wrote: " + file_out

show()

##########
# Make a plot to try to match that of Grundy & Fink 1996
##########
              
xlim = (5000, 10000)
ylim = (0.7, 4)
binning = 30

for i in range(len(files_pluto)):
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i], color=colors[i], fontweight='bold')

plt.title('Pluto/HD, SALT, Compared to Grundy & Fink, binning = ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
show()

##########
# Make a plot to try to match that of Grundy & Fink 1996 -- but over the whole wavelength range
##########
              
# Align on some feature near the middle of wavelength range. Looks like a good one at ~ 6870

align_wavelengths(wavelength, flux_pluto, (6850, 6890), 6865, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (6850, 6890), 6865, smoothing=10)
  
xlim = (3200, 10000)

ylim = (0.3, 4)
binning = 1

for i in range(len(files_pluto)):
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])
              
plt.title('Pluto/HD, SALT, binning = ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
file_out = 'spect_ratio_full_{0}-{1}_bin{2}.png'.format(xlim[0], xlim[1], binning)
plt.savefig(file_out)
print "Wrote: " + file_out

show()


##########
# Make a plot of 5500 - 6500 to try to see O2 bands of Calvin & Spencer
# Should see bands at 5600-5700 and 6200 - 6300.
# SHould see a bit peak up at 5750. I don't see it...
##########
              
# Align on some feature near the middle of wavelength range. Looks like a good one at ~ 6870

align_wavelengths(wavelength, flux_pluto, (6850, 6890), 6865, smoothing=10)
align_wavelengths(wavelength, flux_hd,    (6850, 6890), 6865, smoothing=10)
  
xlim = (5500, 6500)

ylim = (0.3, 4)
binning = 1

for i in range(len(files_pluto)):
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
      plt.text(xlim[1] - 0.1 * (xlim[1]-xlim[0]), 
               flux_pluto[i][wheremin(abs(wavelength - 9000))[0]] + i*0.2 + 0.75, 
              '  ' + date_pluto[i])
              
plt.title('Pluto/HD, SALT, binning = ' + repr(binning), fontsize=24)
plt.xlim(xlim)
plt.ylim(ylim)

plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
file_out = 'spect_ratio_full_{0}-{1}_bin{2}.png'.format(xlim[0], xlim[1], binning)
plt.savefig(file_out)
print "Wrote: " + file_out

show()

# Make plots of region around 0.89 um band
# Loop over this, trying each HD observation, to find the best one for each Pluto spectrum
# Make 12 plots, each with 5 lines

#rcParams['figure.figsize'] = 20, 6
#
#xlim = (8300, 9300)
#index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy. Does it really go to 14?
#for i in range(len(files_pluto)):
#  if i in index_good_spectra_pluto:
#      for j in range(len(files_hd)):
#        plot(wavelength, flux_pluto[i] / flux_hd[j] + j*0.3, color=colors[j], 
#             label = 'Pluto #' + repr(i) + ' / HD #' + repr(j))
#
#      plt.title('Pluto #' + repr(i) + ' / HD #n', fontsize=24)
#      plt.xlim(xlim)
#      plt.ylabel('Intensity [arbitrary]', fontsize=24)
#      plt.xlabel('Wavelength [$\AA$]', fontsize=24)
#      plt.ylim((1.0,3))
#      legend()
#      show()

# Make plot with all Pluto / HD, scaled to match that in Grundy paper

xlim = (5000, 9300)
binning = 5
# rcParams['figure.figsize'] = 12, 6

index_good_spectra_pluto = {0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14} # 5 is messed up; 11 is noisy
for i in range(len(files_pluto)):
  if i in index_good_spectra_pluto:
      plot(wavelength, convolve(flux_pluto[i] / flux_hd[3] + i*0.2, Box1DKernel(width=binning)), color=colors[i])
      plt.text(amax(wavelength), 
               np.median(flux_pluto[i][-200:-50]) + i*0.2, 
               '  ' + date_pluto[i])

plt.title('Pluto/HD, SALT June-August 2014', fontsize=24)
plt.xlim(xlim)
plt.ylabel('Intensity [arbitrary]', fontsize=24)
plt.xlabel('Wavelength [$\AA$]', fontsize=24)
plt.ylim((0.7, 4))

      
      
# Now make plot of the ratios of the HD stars night-to-night
# This shows us the atmospheric variability.
# Concl: it looks like a few percent

##########
# HD stars, mutual ratios, narrow wavelength range
##########

xlim = (8300, 9300)
for i in range(len(files_hd)):
      for j in range(len(files_hd)):
        plot(wavelength, flux_hd[i] / flux_hd[j] + j*1, color=colors[j], 
             label = 'HD #' + repr(i) + ' / HD #' + repr(j))

      plt.title('HD #' + repr(i) + ' / HD #n', fontsize=24)
      plt.xlim(xlim)
      plt.ylabel('Intensity [arbitrary]', fontsize=24)
      plt.xlabel('Wavelength [$\AA$]', fontsize=24)
      plt.ylim((0.5,5.5))
      legend()
      show()

# Now make plot of the ratios of the HD stars night-to-night
# This shows us the atmospheric variability.
# ** Same as above, but different wavelength range

##########
# HD stars, mutual ratios, full wavelength range
##########

xlim = (3300, 9300)
for i in range(len(files_hd)):
      for j in range(len(files_hd)):
        plot(wavelength, flux_hd[i] / flux_hd[j] + j*1, color=colors[j], 
             label = 'HD #' + repr(i) + ' / HD #' + repr(j))

      plt.title('HD #' + repr(i) + ' / HD #n', fontsize=24)
      plt.xlim(xlim)
      plt.ylabel('Intensity [arbitrary]', fontsize=24)
      plt.xlabel('Wavelength [$\AA$]', fontsize=24)
      plt.ylim((0.5,5.5))
      legend()
      show()

# Make a plot of airmasses vs. date, for both HD and Pluto

plot(jd[is_pluto], airmass[is_pluto], color= 'blue', label = 'Pluto', linestyle = 'none', marker = 'o')
plot(jd[is_hd], airmass[is_hd], color='red', label = 'HD', linestyle = 'none', marker = 'o')

plt.ylabel('Airmass')
plt.xlabel('JD')
plt.legend()
show()

# Measure the slope of the blue spectrum
# Do this by fitting a line to the 

# Now experiment with binning


plot(wavelength, convolve(flux_pluto[2] / flux_hd[4], Box1DKernel(width=2)))
plt.ylim((0.4, 1.4))
plt.xlim((8300,9300))
show()

        
        
