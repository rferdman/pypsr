#!/usr/local/bin/python

# Calculate correction to pbdot given errors on input variables
# Output various components to correction, final correction, on both 
# pbdot/pb and pbdot.

# Will only take into account errors on distance and proper motion, and 
# solar galactocentric velocity and distance, since
# errors on position and pb are really tiny and would not contribute.

# final error on pb will be quadrature sum of error on pbdot correction and 
# error on measured pbdot

import sys
from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
import psrcalc as pc
import argparse
import utils
# Do argument parsing.  Errors are not mandatory on input variables.  If errors 
# on all values are not given, then will just do a straight calculation without 
# any fancy business.  If errors are given on all variables required, will do 
# standard error propagation to get final uncertainties.

def get_opt(progname):
    parser = argparse.ArgumentParser( 
            prog=progname,
            description='Calculate correction to pbdot given errors on input variables. Output various components to correction, final correction, on both pbdot/pb and pbdot.')
     
    parser.add_argument('--ra', 
                        nargs=3,
                        type=float,
                        required=True,
                        help='Right ascension, given as tuple (hr, min, sec.sec)')
    parser.add_argument('--dec',
                        nargs=3,
                        type=float,
                        required=True,
                        help='Declination, given as tuple (deg, min, sec.sec)')
    parser.add_argument('--dist',
                        type=float,
                        required=True,
                        help='Pulsar distance in kpc')
    parser.add_argument('--disterr',
                        type=float,
                        help='Uncertainty in pulsar distance in kpc')
    parser.add_argument('--pm',
                        type=float,
                        required=True,
                        help='Total pulsar proper motion in mas/yr')
    parser.add_argument('--pmerr',
                        type=float,
                        help='Uncertainty in total pulsar proper motion in mas/yr')
    parser.add_argument('--pb',
                        type=float,
                        required=True,
                        help='Pulsar orbital period in days')
    parser.add_argument('--pbdot',
                        type=float,
                        required=True,
                        help='Pulsar orbital period derivative in units of 10^12 s/s')
    parser.add_argument('--pbdoterr',
                        type=float,
                        help='Uncertainty on pulsar orbital period derivative in units of 10^12 s/s')
    parser.add_argument('--v0',
                        type=float,
                        default=220.0,
                        help='Solar velocity around Galactic centre, in km/s')
    parser.add_argument('--v0err',
                        type=float,
                        default=20.0,
                        help='Uncertainty in solar velocity around Galactic centre, in km/s')
    parser.add_argument('--R0',
                        type=float,
                        default=7.7,
                        help='Galactocentric distance of Sun, in kpc')
    parser.add_argument('--R0err',
                        type=float,
                        default=0.7,
                        help='Uncertainty in Galactocentric distance of Sun, in kpc')
        
    args=parser.parse_args()

    # Create ra, dec tuples
    args.ra = (args.ra[0], args.ra[1], args.ra[2])
    args.dec = (args.dec[0], args.dec[1], args.dec[2])
    
    # Create pbdot by multiplying by 10^12
    args.pbdot *= 10**(-12.)
    if(args.pbdoterr != None):
        args.pbdoterr *= 10**(-12.)

    return args
  
    
def main():

    progname = 'find_pbdot_corr.py'
    args = get_opt(progname)
    
    ra = 360.0 * (args.ra[0] + args.ra[1]/60. + args.ra[2]/3600.)/24.0
    if(args.dec[0] < 0):
        sgn = -1.
    else:
        sgn = 1.
    dec =  sgn * (np.abs(args.dec[0]) + args.dec[1]/60. + args.dec[2]/3600.)
    print 'RA   = ', ra,  ' deg  =  ', np.radians(ra),  ' rad'
    print 'Dec  = ', dec, ' deg  =  ', np.radians(dec), ' rad'
#    print 'l    = ', l,   ' deg'
#    print 'b    = ', b,   ' deg'
    if(args.disterr != None):
        derr = ' +/- '+str(args.disterr)
    else:
        derr = ''    
    print 'dist = ', args.dist,   derr, ' kpc'
    if(args.pmerr != None):
        pmerr = ' +/- '+str(args.pmerr)
    else:
        pmerr = ''        
    print 'PM   = ', args.pm, pmerr, ' mas/yr'

    print 'Pbdot = ', args.pbdot
    if(args.pbdoterr != None):
        print 'Pbdot error = ', args.pbdoterr

        
        
    # Run calcuation of pbdot correction through each iteration and build 
    # up array of values for pbdot correction
    # for i_iter in np.arange(args.niter):
    pbdot_corr, pdbot_corr_err = pc.get_pbdot_corr(args.ra, args.dec, 
                                            args.pm, args.dist,
                                            args.pb,
                                            pmerr=args.pmerr, 
                                            derr=args.disterr,
                                            v0=args.v0, R0=args.R0,
                                            v0err=args.v0err,
                                            R0err=args.R0err, 
                                            verbose=True)
                    
    # Create histogram of pbdot correction, and fit a gsaussian to it.
    # Extract median (mean?) and sigma as our final value and error.

    print ''
    print 'Original Pbdot = ', args.pbdot

    # Check whether all error bars required are given for distance and proper
    # motion.  If they are, go ahead with the random number generation thing.
    # If not, just do a quick straight calculation of the pbdot correction.
    
    if((args.disterr != None) & (args.pmerr != None) & (args.pbdoterr != None)):
        pbdot_corr_err = ''
        pbdot_final_err = ''

    print 'Correction to Pbdot = ', pbdot_corr, ' +/- ', pbdot_corr_err

    # Now calculate final Pbdot error


    ###print 'FINAL corrected Pbdot = ', pbdot_final, ' +/- ', pbdot_final_err

        
    return
        
main()

