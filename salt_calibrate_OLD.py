# These are batch commands to calibrate SALT data using PyRAF / PySALT

# Program takes the raw FITS files, and applies the SALT pipeline to them.
# Output is .txt file for each spectrum, one per file.
# It does not combine these or make (non-trivial) plots.
#
# HBT 15-Aug-2014

# Set the prompt color. Can't figure out how to put magic command into a config file in .ipython...

# %config PromptManager.in_template = r'{color.LightGreen}In [{color.LightGreen}{count}{color.LightGreen}]: '

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
import subprocess

from   subprocess import call
import string
import glob
import os       # for chdir()
import os.path  # for isfile()
import astropy
from   astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt # pyplot 
from astropy.convolution import convolve, Box1DKernel

import numpy as np
from   pandas import DataFrame, Series
import pandas as pd
from   pylab import *

from matplotlib.path import Path
import matplotlib.patches as patches

##########
# Function to determine if an entered string is a number. 
# From http://stackoverflow.com/questions/354038/how-do-i-check-if-a-string-is-a-number-in-python
##########

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
dir_data_2014 = "/Users/throop/Data/SALT_Pluto_Jun14/product"
dir_data_2015 = "/Users/throop/Data/SALT_Pluto_May15/product"

dir           = '/Users/throop/python'
dir_pysalt    = '/Users/throop/iraf/pysalt'
dir_linelists = dir_pysalt + "/data/linelists"

##########
# Set flags for which of the tasks we want to do on this run. Can set any/all of these.
##########

DO_SPECIDENTIFY = True
DO_CRCLEAN      = True
DO_SPECRECTIFY  = False
DO_SPECSKY      = False
DO_SPECEXTRACT  = False
DO_PLOT_IMAGE_RAW= True
DO_PLOT_IMAGE_SKYSUB = False
DO_PLOT_IMAGE = True
DO_PLOT_SPECT   = True # Make a plot showing the output image
DO_FORCE        = True  # Normally, if output file already exists, don't do the step. This flag forces it to be redone.
DO_QUERY_USER   = True   # Ask the user after every image to go on?
DO_FORCE_CAL_SINGLE = False  # Force calibration for this file (including SPECTRECTIFY, SPECEXTRACT, etc.)?

DO_ARCS_ONLY    = True  # Downselect to only ARCs
DO_ALL_BUT_ARCS = False   # Downselect to only science frames

##########
# Set Parameters for all of the steps
##########

year_obs        = 2015  # Set this to 2014 or 2015, depending on which dataset

mode_cr         = 'edge'  # Mode for CR rejection in SALTCRCLEAN. Can be 'edge' or 'median'.
mode_cr         = 'median'
mode_cr         = 'fast' # Fast still takes 60-90 sec.
iter_cr         = 3
y0_sky_default  = 50      # Position for sky spectral extraction start. 0 = center of Pluto / HD spectrum. Can be overridden.

dy_sky          = 20      # Number of rows to use for the sky spectrum
dy_spect        = 20      # Number of rows to use for the spectrum itself (full width)


# CD into python. Just dump everything there for now.
# I don't know how PyRAF handles directories, so better to cd into the right one from the start.

os.chdir(dir)

# Create a list for each of the keywords. After we do the loop, 
# each list will be populated -- e.g., FITS_EXPTIME = [1, 1, 10, 1], etc

fits_exptime = []	# new list (same idea as array) 
fits_object  = [] 
fits_date_obs= []
fits_utc_obs = [] # starting time of exposure
fits_jd      = [] 
fits_obsmode = [] 
fits_moonang = []  
fits_instrume= [] 
fits_gainset = [] 
fits_rospeed = [] 
fits_ccdsum  = [] # '2 2' -> 2x2 bininng 
fits_pixscale= [] # before binning 
fits_proposer= [] 
fits_propid  = [] 
fits_airmass = [] 
fits_masktyp = [] 
fits_maskid  = []
 
fits_ra      = [] 
fits_dec     = [] 
fits_telra   = [] 
fits_teldec  = [] 

fits_grating = [] 
fits_gr_sta  = [] 
fits_gr_angle= []
fits_lampid=   []

# Get a list of raw FITS files. Do not include *_rect_sky.fits, etc -- just the 'original' fits files.
# Also, it's fine to exclude FLATs, since we can basically ignore them entirely in the pipeline as per Steve Crawford 26-Aug-14.

if (year_obs == 2014):
  file_list = glob.glob(dir_data_2014 + '/mbxgpP2014[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits') # glob won't take regexp

if (year_obs == 2015):
  file_list = glob.glob(dir_data_2015 + '/mbxgpP2015[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].fits') # glob won't take regexp

files = np.array(file_list, dtype='S150') # Need to make string fields long enough for expansion!
 
# Read the JD from each file. Then sort the files based on JD.

# files      = np.array(file_list,dtype='S100')
 
print "Reading " + repr(size(file_list)) + " files"
jd = []
object = []
for file in files:
    sys.stdout.write('.') # print w/o newline
    hdulist = fits.open(file)
    jd.append(hdulist[0].header['JD']) # For now, only read JD and OBJECT
    object.append(hdulist[0].header['OBJECT'])

    hdulist.close()
    
print
print "Sorting by JD, removing flats, removing bad files"

object = np.array(object)


if (DO_ARCS_ONLY):
    arcs = [object == 'ARC']
    files = files[arcs]
    object = object[arcs]
    jd = np.array(jd)
    jd = jd[arcs]

    print
    print 'Doing ARCs only...'

if (DO_ALL_BUT_ARCS):
    arcs = [object != 'ARC']
    files = files[arcs]
    object = object[arcs]
    jd = np.array(jd)
    jd = jd[arcs]

    print
    print 'Doing all but ARCs...'

indices = np.argsort(jd)
files = files[indices]          # Resort by JD    
files = files[object != 'FLAT'].copy() # Remove all of the flats

# Remove some bad files on a one-off basis.
# NB: My philosophy here is to remove as little as possible. Process what I can. Then in the next routine, I can 
#     mix & match.

# files = files[files != dir_data + '/mbxgpP201409050032.fits'] # Pluto spectra. Blank. 0029, 0030, 0031 are all OK.

# 20140905 0{189, 190, 191} : All three of the these are on-target but very faint.
# 20140905 0{198, 199, 200, 201}: All four of these are usable. Not sure why we have four. 198, 199 = best of the bunch.

# files = files[files != dir_data + '/mbxgpP201409050189.fits'] # Pluto spectra. Took four. First two are clearly off-target
# files = files[files != dir_data + '/mbxgpP201409050190.fits'] # Pluto spectra. Took four. First two are clearly off-target
# files = files[files != dir_data + '/mbxgpP201409050198.fits'] # Very faint. Off target? 200 / 201 are OK. Took four.
# files = files[files != dir_data + '/mbxgpP201409050199.fits'] # Very faint. Off target? 200 / 201 are OK. Took four.

# *** 191 looks off-target also??

print

print "Reading " + repr(size(files)) + " files"

for file in files:
    sys.stdout.write('.')
        
    hdulist = fits.open(file)
    
    header = hdulist[0].header
    keys = header.keys()

    fits_object.append(header['OBJECT'])
    fits_exptime.append(header['EXPTIME'])
    fits_date_obs.append(header['DATE-OBS'])
    fits_utc_obs.append(header['UTC-OBS'])
    fits_jd.append(header['JD'])
    fits_obsmode.append(header['OBSMODE'])
    fits_moonang.append(header['MOONANG'])
    fits_instrume.append(header['INSTRUME'])
    fits_gainset.append(header['GAINSET'])
    fits_rospeed.append(header['ROSPEED'])
    fits_ccdsum.append(header['CCDSUM'])
    fits_pixscale.append(header['PIXSCALE'])
    fits_proposer.append(header['PROPOSER'])
    fits_propid.append(header['PROPID'])
    fits_airmass.append(header['AIRMASS'])
    
    fits_dec.append(header['RA'])
    fits_ra.append(header['DEC'])
    fits_telra.append(header['TELRA'])
    fits_teldec.append(header['TELDEC'])
    
    if ('GRATING' in keys):    
        fits_grating.append(header['GRATING'])
        fits_gr_sta.append(header['GR-STA'])
        fits_gr_angle.append(header['GR-ANGLE'])
        fits_masktyp.append(header['MASKTYP']) # not always there
        fits_maskid.append(header['MASKID'])
    else:
       fits_grating.append('')
       fits_gr_sta.append('')
       fits_gr_angle.append('')
       fits_masktyp.append('')
       fits_maskid.append('')

    if ('LAMPID' in keys):
        fits_lampid.append(header['LAMPID'])
    else:
        fits_lampid.append('')

# Convert from lists to NumPy arrays as needed

jd          = np.array(fits_jd)
utc         = np.array(fits_utc_obs)
date        = np.array(fits_date_obs)
airmass     = np.array(fits_airmass)
object      = np.array(fits_object) # np.array can use string arrays as easily as float arrays
instrument  = np.array(fits_instrume)
proposer    = np.array(fits_proposer)
exptime     = np.array(fits_exptime)
lampid	  = np.array(fits_lampid)
grating     = np.array(fits_grating)
tilt        = np.array(fits_gr_angle)
obsmode     = np.array(fits_obsmode)  # Remove extra spaces

# Extract just the good ones
# ** When initializing the np string array, need to set the maximum string length.
# ** I choose a long default here, so we won't exceed it.

# files      = np.array(file_list,dtype='S100')

files_short = np.array(files)

for i in range(files.size):
    files_short[i] = files[i].split('/')[-1]  # Get just the filename itself
   
is_bad = ((proposer != 'Throop') | 
          (instrument != 'RSS') | 
          (object == 'BIAS') | 
          (obsmode == 'FABRY-PEROT') | 
          (files_short == "mbxgpP201408260016.fits"))  # This file is an ARC but it is bad... all noise. 
                                                       # Don't know why. It was retaken, so skip it.
          
is_good = (is_bad == False)

# pdb.set_trace()

files      = files[is_good]  # Is this legal in python? Does it make sense, or is it a view into itself??
files_short=files_short[is_good]
jd         = jd[is_good]
utc        = utc[is_good]
date       = date[is_good]
airmass    = airmass[is_good]
object     = object[is_good]
instrument = instrument[is_good]
proposer   = proposer[is_good]
exptime    = exptime[is_good]
lampid     = lampid[is_good]
grating    = grating[is_good]
obsmode    = obsmode[is_good]
tilt       = tilt[is_good].astype(float) # For some reason this reads as a string from the FITS file...  

# Make arrays for the various output files

print "\ndone reading files\n"

# Get a list of Arcs, Flats, etc. Every image is exactly one of these.

is_arc   = (object == 'ARC')
is_pluto = (object == 'Pluto')
is_hd    = (object == 'HD 146233')
is_flat  = (object == 'FLAT')
                 
## Call specidentify to identify individual lines in the cal spectrum

# "SPECIDENTIFY identies Arc lines in a spectrum and calculates the wavlength
# identification of the lines or it can be run automatically where the task will
# identify the lines without any assistance from the user."
#
# Q: Where do I get the line lists from A:  SAAO-SALT webpage. 

# Display a file
# iraf.imexamine(files[2] + '[1]')  # display an arc spectrum in a ds9 window.
# for file in files do:
# iraf.imhead(filesarr[is_arc][2]+'[0]',long=1)

calfile         = "salt_pluto_cal.txt" # File which gets written by SPECIDENTIFY with the parameters of the wavelength solution.
                                       # This is a running file, which gets appended every time SPECRECTIFY is run. It has many entries.
file_positions  = "salt_positions_extract.csv" # File which I made from Excel and ds9, of which rows to extract.

# Read it into an array of lists
# Problem is that numpy arrays have a single datatype

# f = loadtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, 
#            dtype = {'names' : ('djd', 'utc', 'object', 'exptime', 'tilt', 'file', 'sky', 'spectrum', 'notes'),
#                     'formats': ('f4', 'S10', 'S10',     'f4',     'f4',    'S10', 'S10', 'S10',      'S10')} )

# Read it in as a string array. Then extract the proper columns, and use those. Not very slick.

# f = loadtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, 
#            dtype = 'S300')
# genfromtxt() is more general than loadtxt()

f = genfromtxt(dir + '/' + file_positions, delimiter=',', skiprows=1, dtype = 'str', invalid_raise=False, usecols=range(9))

name_files_positions  = f[:,5]
positions_extract_sky     = f[:,6] # Which rows of each image do I use for the sky subtraction?
positions_extract_spect   = f[:,7] # Which rows of each image do I use for the spectrum extraction?

#########
# Create the array for the sky y-positions 
#########

y0_sky          = y0_sky_default + np.zeros(size(files))

# Change these one-off values for the Y position of the sky. This is because there are stars in these positions.

y0_sky[ np.char.find(files, '06280032') > 0 ] = 40
y0_sky[ np.char.find(files, '06280033') > 0 ] = 40
y0_sky[ np.char.find(files, '06280034') > 0 ] = 40
y0_sky[ np.char.find(files, '06280041') > 0 ] = 40
y0_sky[ np.char.find(files, '06280042') > 0 ] = 40
y0_sky[ np.char.find(files, '06280054') > 0 ] = 30
y0_sky[ np.char.find(files, '06280055') > 0 ] = 30

y0_sky[ np.char.find(files, '07020053') > 0 ] = 25
y0_sky[ np.char.find(files, '07020054') > 0 ] = 25
y0_sky[ np.char.find(files, '07020061') > 0 ] = 25
y0_sky[ np.char.find(files, '07020062') > 0 ] = 25
y0_sky[ np.char.find(files, '09050191') > 0 ] = 40


##########
# Create filenames for all the other files we might load
##########

files_cr              = files.copy()
files_cr_rect         = files.copy()
files_cr_rect_sky     = files.copy()
files_cr_rect_sky_ext = files.copy()

section_sky           = files.copy()
section_spect         = files.copy()

for i in range(files.size):
    files_cr[i]              = string.replace(files[i], '.fits', '_crmode=' + mode_cr + '_criter=' + repr(iter_cr) + '.fits')
    files_cr_rect[i]         = string.replace(files_cr[i], '.fits', '_rect.fits')
    files_cr_rect_sky[i]     = string.replace(files_cr_rect[i], '.fits', '_y0sky=' + repr(y0_sky[i]) + '_dysky=' + repr(dy_sky) + 
                                               '.fits') 
    files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky[i], '.fits', '_dyspect=' + repr(dy_spect) + 
                                               '_y0spect=XXX' + '.txt')
    section_sky[i]           = ''
    section_spect[i]         = ''
    
# Order of operation:
#  1. Spec identify (interactive). No parameters.
#  2. CR rejection / cleaning. Parameters: crmode=edge | fast | median; maxiter=3
#  3. Spec rectify. No parameters. (Based on fit in spec identify)
#  4. Sky subtraction
#  5. Spectral extraction
#  6. Make plot
#
# Example filename:
# image_cr.edge.iter5_rect_sky.50.10_extract

##########    
# Now start the main loop
##########

IS_FIRST_CALL = True

i = 0        

while (i < files.size): # Can't do a for loop

    IS_DONE_PROMPT = False
    
    if (DO_QUERY_USER == False):
        i += 1;
    if IS_FIRST_CALL:    
      print "#    File                   Type     Exptime SPCIDNTFY _cr  _cr_rect #"
      for j in range(size(files)):
          print "{0:3d}. {1} {2:10} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d}.".format(
	          j, 
		  files_short[j], 
		  object[j], 
		  int(exptime[j]), 
		  (subprocess.call(['grep', files_short[j], calfile]) == 0),	# Has file been run thru SPECIDENTIFY yet?
		  os.path.isfile(files_cr[j]), 					# Does _cr.fits file exist?
		  os.path.isfile(files_cr_rect[j]), j)				# Does _cr_rect.fits file exist?
      IS_FIRST_CALL = False;
        
    while (IS_DONE_PROMPT== False):
        if (DO_QUERY_USER):

           inp = raw_input("File " + repr(i) + ": (p)rocess, (#)number, (c)ontinuous w/o questions,\n" + 
              "(q)uit, (l)ist, (f)orce full calibration of file, (n)ext, (" + repr(i) + ") ?: ")



              
        if (inp == 'q'):
              quit()
              
        if (inp == 'c'):                  # 'C'ontinue w/o asking more questions
            DO_QUERY_USER = False        # Actually, I want to ask questions... otherwise, can get in loop that ctrl-c cannot kill!
            IS_DONE_PROMPT = True
              
        if (inp == 'p') | (inp == ''):
            IS_DONE_PROMPT = True
        
        if (inp == 'n'):
            i+=1
            IS_DONE_PROMPT = True;
            
        if ((inp == 'l') | (inp == '?')):
            print "#    File                   Type     Exptime SPCIDNTFY _cr  _cr_rect #"
            for j in range(size(files)):
             print "{0:3d}. {1} {2:10} {3:5d} {4:5d} {5:5d} {6:5d} {7:5d}.".format(
	          j, 
		  files_short[j], 
		  object[j], 
		  int(exptime[j]), 
		  (subprocess.call(['grep', files_short[j], calfile]) == 0),	# Has file been run thru SPECIDENTIFY yet?
		  os.path.isfile(files_cr[j]), 					# Does _cr.fits file exist?
		  os.path.isfile(files_cr_rect[j]), j)				# Does _cr_rect.fits file exist?
          
        if (inp == 'f'):
            DO_FORCE_CAL_SINGLE = True   # Force calibration thru all of the steps for this individual file
            DO_QUERY_USER = True         # If we force calibration, we don't want to run forever.
            IS_DONE_PROMPT = True
    
        if (is_number(inp)):
            i = int(inp)
            IS_DONE_PROMPT = True
	
	print "??\n"
      
    
       
##########
# Step 1: SPECIDENTIFY. Do this for ARCs *only*.
#
# This is the one where we identify spectral lines, and do the actual wavelength calibration
# 
# *** Works fine for i=2 (Ar). Doesn't look too good for i=10 (Ne, dfft tilt).
# *** Looks like I match more lines if I set mdiff = 10, not mdiff=5. This defines the 
#     max difference between catalog and observed wavelength.
#
# How to use SPECIDENTIFY: Hit 'Arc', then 'a', 'b', 'f', 'X'. Residual. 'q'. Residuals are typically < 1.
##########
    print "starting specidentify: i = " + repr(i) + "\n"
    
    DO_FORCE_SPECIDENTIFY = True;			# Redo SPECIDENTIFY, even if already done?
    
    file = files[i]  # Pathname to original FITS file
    
    if (object[i] == "ARC") & (DO_SPECIDENTIFY):

# Search the output file to see if we've already process this

        found = (call(["grep", files[i], calfile]) == 0) # Has this file already been processed?
        print "SPECIDENTIFY: ARC file " + files[i] + " been processed? " + repr(found)

    
        if (found == 0) | (DO_FORCE_SPECIDENTIFY == True):
            
            print "\n ** Calling SPECIDENTIFY, " + files[i] + ", i = " + repr(i) + ", object = " + object[i] + \
              ", lampid = " + \
              lampid[i].replace(" ", "") + ", grating = " + grating[i] + ", tilt = " + repr(tilt[i]) + ' deg' + "\n\n"
         
	    a = tilt[i]
	    if (tilt[i] < 15):
	      guessfile = 'salt_pluto_cal_guess_angle_13.txt'
	    else:
	      guessfile = 'salt_pluto_cal_guess_angle_20.txt'

            mdiff = 10 # differnce between observed line position and catalog line position. In tutorial, this was 5, 
                       # and I missed a lot of lines as a result.

# Figure out which ARC light file to use. As per Steve Crawford, I should use the .salt ones, not the .txt ones.

            if (year_obs == 2015):
	      rmiddle = 342			# 2015 data seemed to be binned more aggressively than 2014, with fewer rows in the output!
	      rstep = 50
	      rstart = rmiddle - rstep

            out = iraf.specidentify(images=files[i], linelist=dir_linelists + '/' + lampid[i].replace(" ", "") + '.salt', 
                outfile=calfile, automethod='MatchZero', 
#                   guesstype='file', guessfile=guessfile, # Works for 2014 obs, but crashes 2015 ones
                guesstype='rss', guessfile=None,
		function='polynomial', order=3, rstep=rstep, rstart=rstart,	
                mdiff=mdiff, thresh=3, startext=0, niter=5, inter=1, clobber=0, logfile='salt.log', verbose=1)

#                guesstype='rss', guessfile=' ', function='polynomial', order=3, rstep=100, rstart='middlerow',
##########
# Step 2: Process cosmic rays. 
# Not really necessary on ARCs, but we do it just for consistency to keep filenames the same.
##########

    if (object[i] in ["Pluto", "HD 146233", "ARC"]) & (DO_CRCLEAN):
      e = os.path.exists(files_cr[i])
      print "SALTCRCLEAN: File " + repr(i) + '/' + repr(files.size) + ' ' + files_cr[i] + " exists: " + repr(e)
      if (not e):
          print "\n ** Calling CRCLEAN(" + files[i] + "), i = " + repr(i) + ", object = " + object[i] + "\n\n"
          s = iraf.saltcrclean(files[i], files_cr[i], '', crtype=mode_cr, maxiter = iter_cr, 
                               logfile = 'salt_crclean.log', verbose=1, clobber=1)

##########
# Step 3: Spec Rectify. Call this for all SCIENCE data (HD and Pluto).
# Can also do on ARCs. It is not necessary, but it shows us how well the SPECIDENTIFY worked.
# This uses the polynomial transformations above to rewrite new FITS files which are wavelength calibrated.
# Seems like a backwards way to do it (aren't there artifacts introduced here?), but that's the process, so OK.
# ** Steve's walkthru used function=polynomial above and function=legendre below. That is an error.
#
# This routine *automatically* finds the proper solution, from the 
#    "The task will find the calibration using the same instrument setup that is closest in time 
#     and use that calibration to correct the data."
#
# This just runs. No input required. Takes about 30 sec per image frame.
# I did get an error on one file: "Improper input: N=4 must not exceed M=0"
# ** Problem image: mbxgpP201406240026.fits: "TypeError: Improper input: N=4 must not exceed M=0" in minpack.pyc
# ** Solved: problem was that this was a Fabrey-Perot image. That was an error to take them.
##########

    if (object[i] in ["Pluto", "HD 146233", "ARC"]) & (DO_SPECRECTIFY) :
      e = os.path.exists(files_cr_rect[i])
      print "SPECRECTIFY: File " + files_cr_rect[i] + " exists: " + repr(e)
      if (not e) | (DO_FORCE_CAL_SINGLE):
        print "\n ** Calling SPECRECTIFY(" + files[i] + "), i = " + repr(i) + ", object = " + object[i] + "\n\n"
        s = iraf.specrectify(images=files_cr[i], solfile=calfile, outimages=files_cr_rect[i], caltype='line', 
          function='polynomial', order=3, inttype='interp', outpref='', w1='None', w2='None', dw='None', nw='None', 
          blank=0.0, clobber=1, logfile='salt_specrectify.log', verbose=1)

      rcParams['figure.figsize'] = 20, 30 
      hdulist = fits.open(files_cr_rect[i])
      image = hdulist['SCI'].data
      dy_spect_plot = 400		# How many rows vertically is the spectrum?
      ny = image.shape[0]	# Look up the central row of the spectrum, and take +- dy/2 rows from that
      plt.imshow(log(image))
      plt.ylim((ny/2 - dy_spect_plot/2, ny/2 + dy_spect_plot/2))

      plt.title(files_short[i] + ", " + object[i] + ", " + date[i] + ' RECT')
      plt.show()
      hdulist.close()


##########       
# Step 4: SPECSKY: Remove the sky from a portion of the image.
# Apply this to science data (ARC and Pluto). 
# Looks like this *cannot* take an input file which is a list of 'section's to use for the sky subtraction.
# I need to manage that myself.
# *** To do: Read the 'section' values from a file. They should be computed for every image.
##########

# Open the FITS file. We do this here because we need to get the position of the spectrum.
# We could do this earlier, but we'd prefer to do this *after* the CR's have been processed.

    hdulist = fits.open(file)
    image = hdulist['SCI'].data

# Look up the brightest row. Or something close to that, using a more intelligent algorithm.
# ** For observations on 20140816, there are two objects in the field. For these, we must limit the search
#    width to +- 90 pixels; otherwise, the star will be chosen, not Pluto.

    row_center = 1022               # Row that the center of the spectrum is often on for RSS
    row_min    = row_center - 200   # Search only near the center.
    row_max    = row_center + 200

    if (file.find('20140816') > 1): # Special case here since bright star in nearby field
        row_min = row_center - 100
        row_max = row_center + 100
        
    rowsum = image.sum(1)
    
# Set it to zero outside the center, so to limit the search to bright lines near the center

    rowsum[0:row_min-1] = 0
    rowsum[row_max:] = 0
    
# Smooth it

    height_smooth = 11 # Must be odd
    rowsum = convolve(rowsum, Box1DKernel(height_smooth)) # Convolving it shifts it because output array is longer than input
    rowsum = roll(rowsum, -(height_smooth-1)/2)           # Shift it back...
#    rowsum = rowsum - np.median(rowsum) # subtract off a background
  
    row_max = where(rowsum == amax(rowsum))[0][0]

    hdulist.close()
    
# Construct a string for PySALT about the rows to use for extraction
  
    section_sky[i] = "[{0:d}:{1:d}]".format(int(row_max + y0_sky[i]), int(row_max + y0_sky[i] + dy_sky))
    
    print "Image " + repr(i) + '/' + repr(size(files)) + ': ' + file + \
      '; object = ' + object[i] + ', exptime = ' + repr(trunc(exptime[i])) + ', brightest row = ' + repr(row_max) 
    
    if (object[i] in ["Pluto", "HD 146233"]) & (DO_SPECSKY):
      e = os.path.exists(files_cr_rect_sky[i])
      print "SPECSKY: File " + files_cr_rect_sky[i] + " exists: " + repr(e) 
      if (not e) | (DO_FORCE_CAL_SINGLE):
          
          print "\n ** Calling SPECSKY(" + files[i] + "), i = " + repr(i) + ", object = " + object[i] + \
                ', section = ' + section_sky[i] + "\n\n"
        
          s = iraf.specsky(images=files_cr_rect[i], outimages=files_cr_rect_sky[i], section=section_sky[i], clobber=1,
                    logfile='salt_specsky.log', outpref='', verbose=1)

     
##########
# SPECEXTRACT: Now that everything is finished, extract the spectrum and write to txt file
# *** To do: read the 'section' values from a file. They should be computed for every image.
##########

    # Now that we have calculated the position of the spectrum, put it into the output filename
    # Spectral position is used for both sky and spectrum extraction, but we only list it here.

    files_cr_rect_sky_ext[i] = string.replace(files_cr_rect_sky_ext[i], 'XXX', repr(row_max))
    section_spect[i] = '[' + repr(row_max + -dy_spect/2) + ':' + repr(row_max + dy_spect/2) + ']'
    
    if (object[i] in ["Pluto", "HD 146233"]) & (DO_SPECEXTRACT):
      e = os.path.exists(files_cr_rect_sky_ext[i])
      print "SPECEXTRACT: File " + files_cr_rect_sky_ext[i] + " exists: " + repr(e)
      if (not e) | (DO_FORCE_CAL_SINGLE):
# Look up from text file which rows to use ** NO
        
#        w = np.core.defchararray.startswith(name_files_positions, files_short[i].replace('.fits', '')) 
#        position_extract_sky = positions_extract_spect[w][0]
          
        print "\n ** Calling SPECEXTRACT(" + files_cr_rect_sky[i] + "), i = " + repr(i) + ", object = " + object[i] + \
          ", section = " + section_spect[i] + "\n\n"      
        iraf.specextract(images=files_cr_rect_sky[i], outfile=files_cr_rect_sky_ext[i], method='normal', 
          section=section_spect[i],
          thresh=3.0, minsize=3.0, 
          outformat='ascii', convert=1, clobber=True, logfile='salt_spectextract.log', verbose=True)

##########
# Step 6: Make a plot of the image, showing where the sky and spectral extraction regions are
##########

    if (object[i] in ["Pluto", "HD 146233"]) & (DO_PLOT_IMAGE):

# Load the FITS file
 
#        hdulist = fits.open(files_cr_rect_sky[i])
        hdulist = fits.open(files_cr_rect_sky[i])
        image = hdulist['SCI'].data

# Define the size of the plot to make, in pixels

        dx = 300
        dy = 300 # 300 -> 150 pixels above and below the spectrum itself

# Define the size of the plot, in cm

        rcParams['figure.figsize'] = 10, 10 # Make plot normal size

# Calc the LHS and RHS of the image
        
        x0 = 5
        x1 = image.shape[1] - x0

# Calc the top/bottom of the regions to extract for spectrum and sky

# Extract two halves of the image, and paste them together
        
        image_extract_1 = image[:, 0:dx/2] # Left portion of image
        image_extract_2 = image[:, -dx/2:] # Right portion of image
        image_extract = np.hstack((image_extract_1, image_extract_2))

        fig = plt.imshow(log(image_extract), origin='lower', interpolation='none') # 'lower' means that row order matchs that in DS9
         
        x1 = dx - x0

# Parse the section_sky variable and extract its values
         
        y_sky_plot   = np.array(string.split(string.replace(string.replace(section_sky[i], '[', ''), ']', ''), ':'),dtype='float')
        y_spect_plot = np.array(string.split(string.replace(string.replace(section_spect[i], '[', ''), ']', ''), ':'), dtype='float')
        
        plot([x0, x1, x1, x0, x0], [y_spect_plot[0], y_spect_plot[0], y_spect_plot[1], y_spect_plot[1], y_spect_plot[0]], 
             color='black')
        plot([x0, x1, x1, x0, x0], [y_sky_plot[0],   y_sky_plot[0],   y_sky_plot[1],  y_sky_plot[1],    y_sky_plot[0]],   
             color='blue')
        plt.ylim((y_spect_plot[0]-dy/2, y_spect_plot[0] + dy/2))
        plt.xlim((0, dx))
        plt.title(files_short[i] + ", " + object[i] + ", " + date[i] + ', SKY = ' + section_sky[i] + 
                  ", SPECTRUM = " + section_spect[i], fontsize=25)
        plt.text(dx/2, y_sky_plot[0],   "SKY: "   + section_sky[i], color='white', fontsize=30)
        plt.text(dx/2, y_spect_plot[0], "SPECT: " + section_spect[i], color='white', fontsize=30)

        plt.show()
        
        hdulist.close()

##########
# End of loop: increase counter and loop again
##########

#    i += 1  # Except we really only want to increase the counter if we are looping continuously...
    DO_FORCE_CAL_SINGLE = False
        
##########
# Now plot the spectra
##########

files_spec = files_cr_rect_sky_ext # We've already created the filenames

fluxes        = []
dfluxes       = []
wavelengths   = []

for i in range(size(files_spec)):
    file = files_spec[i]
    if (object[i] == 'ARC'):
        print "Skipping ARC: " + file
        wavelengths.append([])
        fluxes.append([])
        dfluxes.append([])
        
    else: 
        print "Reading spectrum " + file
        d = loadtxt(file)
        wavelengths.append(d[:,0])
        fluxes.append(d[:,1])
        dfluxes.append(d[:,2])

rcParams['figure.figsize'] = 20, 30 # Make plot normal size

colors = np.array(['grey', 'brown', 'orange', 'yellow', 'blue', 'red', 'pink', 'green'])
for i in range(len(fluxes)):
    if (object[i] == 'Pluto'):
        f = np.array(fluxes[i]).astype(float)
        f = f / np.amax(f)  # normalize it
        plot(wavelengths[i],f+i, color=colors[i % size(colors)], ls='-')
    plt.axis([3000,9000,0,100])
  
plt.xlim((3000,9000))
        
        
##########
# Now extract the HD spectra
##########
        
w_spect_hd_blue  = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 13) & (tilt < 14))[0]
w_spect_hd_red   = where((instrument == 'RSS') & (object == 'HD 146233') & (tilt > 20) & (tilt < 21))[0]

w_spect_pluto_blue  = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 13) & (tilt < 14))[0]
w_spect_pluto_red   = where((instrument == 'RSS') & (object == 'Pluto') & (tilt > 20) & (tilt < 21))[0]



        
        
