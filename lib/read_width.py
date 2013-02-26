#!/usr/local/bin/python

from sys import argv, exit
import numpy as np
from psrread import read_widths


def read_width_rad(data_file, rfile=None):

    # Open data file

#    try:
#        f_data = open(data_file, 'r')
#    except IOError as (errno, strerror):
#        if (errno == 2):  # file not found
#            print "IOError ({0}): File".format(errno), proffile, "not found"
#        else:
#            print "IOError ({0}): {1}".format(errno, strerror)
#            exit

#    mjd, width, width_err = np.loadtxt(f_data, dtype='double', comments='#', \
#                                           usecols=(0,1,2), unpack=True)
 
# Use my standard width reader -- returns widths in units of phase
    w = read_widths(data_file, rfile=rfile)

# If we want to adjust widths by comparing means of a reference width file and 
# subtracting by difference of means, then user must provide a reference width 
# file (same format), and this adjustment will be made.
#    if(rfile != None):
#        rw = read_widths(rfile)
#        w['width'] = w['width'] - \
#            (np.average(w['width'], weights=1./(w['werr']**2)) - \
#                np.average(rw['width'], weights=1./(rw['werr']**2)))
#        w['werr'] = np.sqrt(w['werr']**2 + \
#                                (np.std(w['width']))**2 + \
#                                (np.std(rw['width']))**2 )
             
   
    # get widths in radians from units of pulse phase..
    y = 2.*np.pi*w['width']
#    y = np.cos(y_rad)
    
    # Find corresponding y_err. First convert to half-width error, in radians...
    y_err = 2.*np.pi*w['werr']
    # Now get error of y=cos(phi0)...
#    y_err = (y_err_rad*y_err_rad*np.cos(y_rad)/2.) + np.sin(y_rad)*y_err_rad

    return w['mjd'], y, y_err


def read_width_cos_phi(data_file, rfile=None):

    # Open data file
#    try:
#        f_data = open(data_file, 'r')
#    except IOError as (errno, strerror):
#        if (errno == 2):  # file not found
#            print "IOError ({0}): File".format(errno), proffile, "not found"
#        else:
#            print "IOError ({0}): {1}".format(errno, strerror)
#            exit

#    mjd, width, width_err = np.loadtxt(f_data, dtype='double', comments='#', \
#                                           usecols=(0,1,2), unpack=True)

# Use my standard width reader -- returns widths in units of phase
    w = read_widths(data_file, rfile=rfile)

# If we want to adjust widths by comparing means of a reference width file and 
# subtracting by difference of means, then user must provide a reference width 
# file (same format), and this adjustment will be made.
#    if(rfile != None):
#        rw = read_widths(rfile)
#        w['width'] = w['width'] - \
#            (np.average(w['width'], weights=1./(w['werr']**2)) - \
#                np.average(rw['width'], weights=1./(rw['werr']**2)))
#        w['werr'] = np.sqrt(w['werr']**2 + \
#                                (np.std(w['width']))**2 + \
#                                (np.std(rw['width']))**2 )



    # Convert widths to cos(phi0) = cos(half width), where phi0 is in radians
    # First convert widths from units of pulse phase to radians, then divide by 2 to
    # get half-widths.  Then take cosine of these to get y array:
    phi = 2.*np.pi*w['width']/2.  # Multiply by 2 then divide by 2 for clarity :)
    cos_phi = np.cos(phi)

    # Find corresponding y_err. First convert to half-width error, in radians...
    phi_err = 2.*np.pi*w['werr']/2.
    # Now get error of y=cos(phi0)...
    # cos_phi_err = np.abs(np.sin(phi)*phi_err - (phi_err*phi_err*np.cos(phi)/2.))

# Simple error method: average the diff between max/min of cos and cos itself.
    cos_max = np.cos(phi + phi_err)
    cos_min = np.cos(phi - phi_err)

    # cos_phi_err = np.mean([np.abs(cos_max - cos_phi), np.abs(cos_min - cos_phi)])
    
    cos_phi_err = np.zeros_like(cos_phi)
    for i_cos in np.arange(len(cos_phi_err)):
        cos_phi_err[i_cos] = np.mean([np.abs(cos_max[i_cos] - cos_phi[i_cos]), np.abs(cos_min[i_cos] - cos_phi[i_cos])])

#    print "phi = ", phi
#    print "phi_err = ", phi_err
    
#    print ""
#    print "cos_phi = ", cos_phi
#    print "cos_phi_err = ", cos_phi_err


    return w['mjd'], cos_phi, cos_phi_err
### Some global parameters
#prec_period = 297.52*365. # precession period in days - 1913
# prec_period = 25915.0 # precession period in days

# incl = (47.2)*np.pi/180. # orbital inclination in radians -- 1913
#incl = (180. - 47.2)*np.pi/180. # orbital inclination in radians -- 1913
# incl = 1.55 # orbital inclination in radians

#data_file = 'widths.dat'

#mjd, y, y_err = read_data(data_file)

#mjd_mid = np.mean([np.min(mjd), np.max(mjd)])
    
    

