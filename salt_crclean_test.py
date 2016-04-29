file = 'mbxgpP201506300017.fits'
dir = '/Users/throop/data/SALT_Pluto_2015/test/'

# Program to test what parameters for SALT_CRCLEAN are best. 
# Concl: 4 iterations of median.
#
# HBT 22-Jul-2015

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
from   subprocess import call   # Shell output
from   subprocess import check_output   # Shell output

import numpy as np
from   pandas import DataFrame, Series
from   pylab import *
from   astropy.convolution import convolve, Box1DKernel
import cspice

from scipy.optimize import curve_fit

from pyraf import iraf
from iraf import pysalt
from iraf import saltspec
from iraf import specextract
from iraf import specprepare
from iraf import specidentify
from iraf import specrectify
from iraf import specsky
from iraf import saltred

mode_cr         = 'edge' # Fast still takes 60-90 sec.
#   mode_cr         = 'fast' 
#   mode_cr         = 'median'
iter_cr         = 4

file_in = dir + file
file_out = file_in.replace('.fits', '_cr_' + mode_cr + '_' + repr(iter_cr) + '.fits')

verbose = 2

s = iraf.saltcrclean(dir + file, file_out, '', crtype=mode_cr, maxiter = iter_cr, 
                                 logfile = 'salt_crclean.log', verbose=(verbose > 1), clobber=1)


quit

####

# Now try stacking: median, and then edge
# Concl: This is not useful. 'edge' introduces more artifacts.
# Best deal is just to use 4 iterations of 'median'

mode_cr = 'median'

file_in = dir + file
file_out = file_in.replace('.fits', '_cr_' + mode_cr + '_' + repr(iter_cr) + '.fits')

s = iraf.saltcrclean(dir + file, file_out, '', crtype=mode_cr, maxiter = iter_cr, 
                                 logfile = 'salt_crclean.log', verbose=(verbose > 1), clobber=1)

mode_cr = 'edge'

file_in = file_out

file_out = file_in.replace('.fits', '_cr_' + mode_cr + '_' + repr(iter_cr) + '.fits')

s = iraf.saltcrclean(file_in, file_out, '', crtype=mode_cr, maxiter = iter_cr, 
                                 logfile = 'salt_crclean.log', verbose=(verbose > 1), clobber=1)
