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
from pdfplot import *
# Do argument parsing.  Errors are not mandatory on input variables.  If errors 
# on all values are not given, then will just do a straight calculation without 
# any fancy business.  If errors are given on all variables required, will do 
# a Monte-Carlo-style calculation over some number of iterations, and get 
# values and errors out that way.

prob_intervals = np.array([0.683, 0.954, 0.9973])


def get_opt(progname):
    parser = argparse.ArgumentParser( 
            prog=progname,
            description='Calculate correction to pbdot given errors on input variables. Output various components to correction, final correction, on both pbdot/pb and pbdot.')
     
    parser.add_argument('--niter',
                        type=int,
                        default=32000,
                        help='Number of iterations over which to process calculation [3200]')
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
    parser.add_argument('--corrlim',
                        nargs=2,
                        type=float,
                        help='Limits over which to determine orbital decay correction through histogram, in units of 10^-12 s/s')
                                
        
    args=parser.parse_args()

    # Create ra, dec tuples
    args.ra = (args.ra[0], args.ra[1], args.ra[2])
    args.dec = (args.dec[0], args.dec[1], args.dec[2])
    
    # Create pbdot by multiplying by 10^12
    args.pbdot *= 10**(-12.)
    if(args.pbdoterr != None):
        args.pbdoterr *= 10**(-12.)

    return args
  
  
def display_status(iteration, n_iter,):
    utils.restart_line()
#    sys.stdout.write('{0:<10d}[{1:>3d}%]  {2:8.4f}  {3:8.4f}  {4:10.4f}  {5:8.4f} +/- {6:8.4f}  {7:8.4f} +/- {8:8.4f}  {9:14.2g}  {10:14.2g}'.format( \
    sys.stdout.write('{0:<10d}[{1:>3d}%]'.format( \
             iteration, int(100*float(iteration)/float(n_iter))))
    sys.stdout.flush()
  
    
def main():

    progname = 'find_pbdot_corr_mc.py'
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

    # Check whether all error bars required are given for distance and proper
    # motion.  If they are, go ahead with the random number generation thing.
    # If not, just do a quick straight calculation of the pbdot correction.
    
    if((args.disterr != None) & (args.pmerr != None) & (args.pbdoterr != None)):
        
        # For each of the varying parameters (dist, pm, v0, R0), generate
        # args.niter gaussian deviates, with mean = value and sigma = error
        dist = np.random.normal(args.dist, args.disterr, args.niter)
        pm   = np.random.normal(args.pm, args.pmerr, args.niter)
        v0   = np.random.normal(args.v0, args.v0err, args.niter)
        R0   = np.random.normal(args.R0, args.R0err, args.niter)

        # Plot distributions for each parameter:
        fig_param = plt.figure()
        # distance
        ax_param = fig_param.add_subplot(221)
        hist, bin_val = np.histogram(dist, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        A, mu, sig = utils.fitgauss(bin_val, hist)
        #ax_param.plot(bin_val, utils.gaussian(bin_val, A, mu, sig), color='red')
        ax_param.set_xlabel('Distance (kpc)')
        # proper motion
        ax_param = fig_param.add_subplot(222)
        hist, bin_val = np.histogram(pm, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Proper motion (mas/yr)')
        # v0
        ax_param = fig_param.add_subplot(223)
        hist, bin_val = np.histogram(v0, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Solar velocity (km/s)')
        # R0
        ax_param = fig_param.add_subplot(224)
        hist, bin_val = np.histogram(R0, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Galactocentric distance (kpc)')
        
        plt.show()
        
        
        # Run calcuation of pbdot correction through each iteration and build 
        # up array of values for pbdot correction
        #for i_iter in np.arange(args.niter):
        sys.stdout.write('\n {0:<7s}  {1:<7s} \n'.format('Iter', '% done'))
        pbdot_corr_array = []
        for i_iter in np.arange(args.niter):
            pbdot_corr_array.append(pc.get_pbdot_corr(args.ra, args.dec, 
                                    pm[i_iter], dist[i_iter],
                                    args.pb,
                                    v0=v0[i_iter],
                                    R0=R0[i_iter]))
            if(np.fmod(i_iter, 16)==0):
                display_status(i_iter, args.niter)
                
        pbdot_corr_array = np.array(pbdot_corr_array)
        
        # Create histogram of pbdot correction, and fit a gsaussian to it.
        # Extract median (mean?) and sigma as our final value and error.
        if(args.corrlim == None):
            corrlim = (np.amin(pbdot_corr_array)-0.1*np.abs(np.amin(pbdot_corr_array)), 
                       np.amax(pbdot_corr_array)+0.1*np.abs(np.amax(pbdot_corr_array)))
        else:
            corrlim = (args.corrlim[0] * 10.0**(-12), args.corrlim[1] * 10**(-12))
           
        pbdot_pdf, bin_val = np.histogram(pbdot_corr_array, 
                      range=corrlim, bins=192, density=True)
        bin_size = bin_val[1]-bin_val[0]
        pbdot_x = bin_val[0:len(bin_val)-1] + 0.5*bin_size
        #fig_pbdot = plt.figure()
        #ax_pbdot = fig_pbdot.add_axes([0.12, 0.1, 0.8, 0.85])
        #ax_pbdot.plot(pbdot_x, pbdot_pdf, linestyle='steps-mid', linewidth=1.4, 
        #            color='blue')
        
        pbdot_corr_med, pbdot_corr_min, pbdot_corr_max = \
            get_pdf_prob(pbdot_x, pbdot_pdf, prob_intervals, norm=True) #, \
        plot_pdf(pbdot_x, pbdot_pdf, \
###                     weights=alpha_weights, \
             xlabel='$\\dot{P}_{b, corr}}$', \
             ylabel='Probability density',\
             prob_lines=np.append(pbdot_corr_min, pbdot_corr_max),\
             prob_linestyle=['dashed','dashdot','dotted', \
                                 'dashed','dashdot','dotted'], \
             )
        plt.axvline(pbdot_corr_med)
        pbdot_corr_err_high = np.abs(pbdot_corr_max-pbdot_corr_med)
        pbdot_corr_err_low  = np.abs(pbdot_corr_min-pbdot_corr_med)
        
        #A, pbdot_corr, pbdot_corr_err = utils.fitgauss(pbdot_x, pbdot_pdf)
        #ax_pbdot.plot(pbdot_x, utils.gaussian(pbdot_x, A, pbdot_corr, pbdot_corr_err), 
                    #linestyle='solid', linewidth=1.4, color='red')
        plt.savefig('pbdot_new.png')
        
        print ''
        print 'Original Pbdot = ', args.pbdot
        print ''
        #print 'Corrected Pbdot = ', pbdot_corr, ' +/- ', pbdot_corr_err
        print 'Pbdot correction = ', pbdot_corr_med, ' +', pbdot_corr_err_high[0], \
        ' -', pbdot_corr_err_low[0]
        print 'Pbdot correction peak = ', pbdot_x[np.argmax(pbdot_pdf)]
        
        pbdot_new = args.pbdot + pbdot_corr_med
        pbdot_corr_err_mean = np.mean(np.array([pbdot_corr_err_low[0],pbdot_corr_err_high[1]]))
        pbdot_new_err = np.sqrt(args.pbdoterr**2.0 + pbdot_corr_err_mean**2.0)
        print 'Pbdot new = ', pbdot_new, ' +/- ', pbdot_new_err

        
    else:
        
        pbdot_corr = pc.get_pbdot_corr(args.ra, args.dec, args.pm, args.dist, 
                                    args.pb, v0=args.v0, R0=args.R0, 
                                    verbose=True)
            
        pbdot_new = args.pbdot - pbdot_corr
        
        print '\nWill not calculate final uncertainties...\n'
        print 'Original Pbdot = ', args.pbdot
        print 'Corrected Pbdot = ', pbdot_new
        
        return
        
main()

