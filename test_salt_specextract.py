# These are batch commands to calibrate SALT data using PyRAF / PySALT

# Program takes the raw FITS files, and applies the SALT pipeline to them.
# Output is .txt file for each spectrum, one per file.
# It does not combine these or make (non-trivial) plots.
#
# HBT 15-Aug-2014

from pyraf import iraf
from iraf import pysalt
from iraf import saltspec
from iraf import specextract
from iraf import specprepare
from iraf import specidentify
from iraf import specrectify
from iraf import specsky
from iraf import saltred

import sys # for sys.sdout.write()

import pdb

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
import astropy
from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 

import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *

from matplotlib.path import Path
import matplotlib.patches as patches


dir_data      = "/Users/throop/Data/SALT_Pluto_Jun14/product"
dir           = '/Users/throop/python'
dir_linelists = dir + "/linelists"

## Select one of these 'file' settings below, to choose betwteen HD, Pluto, Red, and Blue

# file = dir_data + "/mbxgpP201407010021_cr_rect_sky.fits" # A 'typical' 300-sec blue-tilt Pluto spectrum [13.25  deg]. Concl: 20 rows is fine.
# file = dir_data + "/mbxgpP201407010034_cr_rect_sky.fits" # A 'typical' 300-sec red- tilt Pluto spectrum [20.375 deg]. Concl: 20 rows is fine.
# file = dir_data + "/mbxgpP201408120034_cr_rect_sky.fits" # A 'typical' 300-sec blue-tilt HD spectrum [13.375 deg]. 20 rows is fine.
file = dir_data + "/mbxgpP201408120047_cr_rect_sky.fits" # A 'typical' 300-sec red- tilt Pluto spectrum [20.375 deg]. 20 rows is fine.

file = dir_data + "/mbxgpP201407010034_cr_rect.fits"

os.chdir(dir_data)

# Open the FITS file

hdulist = fits.open(file)
image = hdulist['SCI'].data

# Look up the brightest row

rowsum = image.sum(1)
row_center = where(rowsum == amax(rowsum))[0][0]  # 1022 for this dataset

list_height = []
files_txt = []
files_fits = []

num_lines_target_max = 20  # How many lines to sum over. This is the maximum.
step_num_lines_target = 20 # Stepsize -- not number of steps

num_lines_sky_max = 50 # How many lines of the sky to include
step_num_lines_sky = 5
y0_sky = 50 # Offset (up/down) from the center of the signal line


# Concl: ideal  num_lines is 20-40 or so? SNR keeps increasing even as I go into 100+, but I'm guessing I'm picking up 
# sky lines there and not actual target spectra, so it's hard to say. Clearly SNR is very high at 40 rows, so let's do that.
# For safety, let's go to 30 rows.... weaker signal but better odds of avoiding the sky. (And we have enough SNR as it is!)

DO_VARY_HEIGHT_SKY    = True
DO_VARY_DIST_SKY      = False
DO_VARY_HEIGHT_TARGET = False

# Vary the sky extraction.
# In looking at the FITS files, it looks like if I go anywhere within +- 200 pixels of the central image row, that is fine.
# Beyond that the spectrum starts to get curved. So keep it < 200 pixels away. But fine until then.

##########
# Loop over the height of the sky to extract. Create a FITS file for each of several heights.
##########

# Concl: height = 20 rows of pixels seems to be pretty good. And perhaps no coincidence, the actual signal (in FITS file) looks about this high, too.

if (DO_VARY_HEIGHT_SKY):
   for height in range(1,num_lines_sky_max, step_num_lines_sky):
#      w = np.core.defchararray.startswith(name_files_positions, files_short[i].replace('.fits', '')) 

        outfile = file.replace(".fits", "_hskyextract" + repr(int(height)) + '_y0skyextract' + repr(int(y0_sky)) + '.fits')
          
        section = '[' + repr(row_center + y0_sky) + ':' + repr(row_center + y0_sky + height-1) + ']'

        print "\n ** Calling SPECSKY(" + file + ', section = ' + section + ", outfile = " + outfile + "\n\n"
            
        s = iraf.specsky(images=file, outimages=outfile, section=section, clobber=1,
                    logfile='salt_specsky.log', outpref='', verbose=1)
        files_fits.append(outfile) # Append to a list of files we've written
   
##########           
# Loop over the height (i.e., number of rows) and call SPECEXTRACT for each one
# SPECEXTRACT takes a FITS file in, and creates a TXT file. 
# It only extracts the spectrum -- it does *not* subtract the sky. That has already been done in SPECSKY
##########
# Concl: 50 pixels up from center looks good. And height of the sky box doesn't really matter.
# Well, height does matter some... spikes definitely become better defined with a wider sky box. But I don't know if these are real, or an artifact!            
# But, 10 pixel height is certainly sufficient.
   
if (DO_VARY_HEIGHT_TARGET):
  for height in range(1,num_lines_target_max,step):
    
    section = repr(int(row_center - height/2)) + ':' + repr(int(row_center + height/2))
    outfile = file.replace(".fits", "_hextract" + repr(int(height)) + '.txt')
      
    print "\n ** Calling SPECEXTRACT(" + file + "), section = " + section + "\n\n"      
    iraf.specextract(images=file, outfile=outfile, method='normal', 
      section=section, 
      thresh=3.0, minsize=3.0, 
      outformat='ascii', convert=1, clobber=True, logfile='salt_spectextract.log', verbose=True)

    list_height.append(height)
    files_txt.append(outfile)

if not(DO_VARY_HEIGHT_TARGET):
  height_target = 20
  for file in files_fits:
    
    section = repr(int(row_center - height_target/2)) + ':' + repr(int(row_center + height_target/2))
    outfile = file.replace(".fits", "_hextract" + repr(int(height_target)) + '.txt')
      
    print "\n ** Calling SPECEXTRACT(" + file + "), section = " + section + "\n\n"      
    iraf.specextract(images=file, outfile=outfile, method='normal', 
      section=section, 
      thresh=3.0, minsize=3.0, 
      outformat='ascii', convert=1, clobber=True, logfile='salt_spectextract.log', verbose=True)

    list_height.append(height)
    files_txt.append(outfile)
    

        
##########
# Now plot the spectra
##########

fluxes        = []
dfluxes       = []
wavelengths   = []

# *** Should look up 'genfromtext' or 'getfromcsv' in numpy
# Actually, SciPy FAQ says that numpy.loadtxt() is probably the way to do it.
# It will deal with header lines, CSV, etc.
# Astropy also has a 'unified table reader'

# Read all of the data into a list of arrays.
# It might be possible to do this as an array of arrays -- not sure.
# See 'advanced topics' section in NumPy book.
# File is a simple three-column list of numbers, separated by spaces

for i in range(len(files_txt)):  # Files_txt is created when writing output. It is a list of all the newly-created files from SPECEXTRACT
    infile = files_txt[i]    
    d = loadtxt(infile)
    print "Reading file " + repr(i) + "/" + repr(len(files_txt)) + ": " + files_txt[i]
    wavelengths.append(d[:,0])
    fluxes.append(d[:,1])
    dfluxes.append(d[:,2])

rcParams['figure.figsize'] = 20, 15 # Make plot  size

colors = np.array(['grey', 'brown', 'orange', 'yellow', 'blue', 'red', 'pink', 'green'])

for i in range(len(fluxes)):
    f = fluxes[i].astype(float)
    f = f / np.amax(f)  # normalize it
    plot(wavelengths[i],f + i * 0.3, color=colors[i % len(colors)], ls='-', label = files_txt[i].split('/')[-1])
#    plt.text(amax(wavelengths[i]), f[-10] + i*0.3, "     " + files_txt[i].split("/")[-1])
legend()
plt.xlim( (amin(wavelengths[0]), amax(wavelengths[0])+2000) )
plt.ylim((0,len(fluxes) * 0.3 + 1))
show()
       
for i in range(len(fluxes)):
    f = fluxes[i].astype(float)
    f = f / np.amax(f)  # normalize it
    plot(wavelengths[i],f + i * 0.1, color='black', ls='-')
#    plt.xlim((6800,8000))
    if (amin(wavelengths[0]) < 4000): 
        plt.xlim((4000,4500))
    else:
        plt.xlim((7500,8000))
#    plt.ylim((0.5,len(fluxes) * 0.1 + 1.0))
#    plt.text(7200, f[-10] + i*0.1, "     " + repr(list_height[i]))

show()
    