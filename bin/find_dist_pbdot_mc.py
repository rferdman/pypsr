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
from astropy import units as u
from astropy import constants as cst
from astropy import coordinates as coord
from scipy.optimize import newton, ridder
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
    parser.add_argument('--raerr', 
                        type=float,
                        help='Uncertainty in right ascension, in seconds')
    parser.add_argument('--dec',
                        nargs=3,
                        type=float,
                        required=True,
                        help='Declination, given as tuple (deg, min, sec.sec)')
    parser.add_argument('--decerr', 
                        type=float,
                        help='Uncertainty in declination, in arcseconds')
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
    parser.add_argument('--xpbdot',
                        type=float,
                        required=True,
                        help='Observed minus GR-predicted orbital period derivative in units of 10^-14 s/s')
    parser.add_argument('--xpbdoterr',
                        type=float,
                        help='Uncertainty on observed minus GR-predicted orbital period derivative in units of 10^-14 s/s')
    parser.add_argument('--v0',
                        type=float,
                        default=240.0,
                        help='Solar velocity around Galactic centre, in km/s')
    parser.add_argument('--v0err',
                        type=float,
                        default=8.0,
                        help='Uncertainty in solar velocity around Galactic centre, in km/s')
    parser.add_argument('--R0',
                        type=float,
                        default=8.34, # From Reid et al. 2009
                        help='Galactocentric distance of Sun, in kpc')
    parser.add_argument('--R0err',
                        type=float,
                        default=0.16,  # From Reid et al. 2009
                        help='Uncertainty in Galactocentric distance of Sun, in kpc')
    parser.add_argument('--distguess',
                        type=float,
                        default=2.4,
                        help='Initial guess for distance, in kpc')
    parser.add_argument('--stepguess',
                        type=float,
                        default=0.05,
                        help='Step size for moving distance guess if unwanted result (e.g. distance < 0), in kpc')
    parser.add_argument('--distlim',
                        nargs=2,
                        type=float,
                        help='Limits over which to determine distance through histogram, in kpc')
                                
        
    args=parser.parse_args()

    # Create ra, dec tuples
    args.ra = (args.ra[0], args.ra[1], args.ra[2])
    args.dec = (args.dec[0], args.dec[1], args.dec[2])
    
    # Create pbdot by multiplying by 10^12
    args.xpbdot *= 10**(-14.)
    if(args.xpbdoterr != None):
        args.xpbdoterr *= 10**(-14.)

    return args
  
  
def display_status(iteration, n_iter, n_bad, n_negative):
    utils.restart_line()
#    sys.stdout.write('{0:<10d}[{1:>3d}%]  {2:8.4f}  {3:8.4f}  {4:10.4f}  {5:8.4f} +/- {6:8.4f}  {7:8.4f} +/- {8:8.4f}  {9:14.2g}  {10:14.2g}'.format( \
    sys.stdout.write('{0:<10d}[{1:>3d}%]      {2}          {3}'.format( \
             iteration, int(100*float(iteration)/float(n_iter)),
             n_bad, n_negative))
    sys.stdout.flush()
  
    
def main():

    progname = 'find_dist_pbdot_mc.py'
    args = get_opt(progname)
    
    ra_deg = 360.0 * (args.ra[0] + args.ra[1]/60. + args.ra[2]/3600.)/24.0
    if(args.raerr != None):
        raerr = (args.raerr/3600.)*360./24.0  # Convert to degrees
        raerr_str = ' +/- '+str(raerr)
    else:
        raerr_str = ''    

    if(args.dec[0] < 0):
        sgn = -1.
    else:
        sgn = 1.
    dec_deg =  sgn * (np.abs(args.dec[0]) + args.dec[1]/60. + args.dec[2]/3600.)
    if(args.decerr != None):
        decerr = args.decerr/3600.   # Convert to degrees
        decerr_str = ' +/- '+ str(decerr)
    else:
        decerr_str = ''
        
    print 'RA   = ', ra_deg,  raerr_str, ' deg  '
    print 'Dec  = ', dec_deg, decerr_str, ' deg  '
    # Convert nominal values of ra and dec to l and b here just for show
    c = coord.ICRSCoordinates(ra=ra_deg, dec=dec_deg, unit=(u.degree, u.degree))
    l = c.galactic.l.degrees
    b = c.galactic.b.degrees
    print 'l    = ', c.galactic.l.degrees, ' deg   =  ', c.galactic.l.radians, ' rad'
    print 'b    = ', c.galactic.b.degrees, ' deg   =  ', c.galactic.b.radians, ' rad'

    if(args.xpbdoterr != None):
        xpbdoterr = ' +/- '+str(args.xpbdoterr)
    else:
        xpbdoterr = ''    
    print 'Xpbdot = ', args.xpbdot,  xpbdoterr
    if(args.pmerr != None):
        pmerr = ' +/- '+str(args.pmerr)
    else:
        pmerr = ''        
    print 'PM   = ', args.pm, pmerr, ' mas/yr'
    
    # convert pb to seconds
    args.pb = u.day.to(u.s, args.pb)
    
    # convert velocity/ velocity error to kpc/s
    # No need to convert R0, which is already in R0
    args.v0 = u.km.to(u.kpc, args.v0)
    args.v0err = u.km.to(u.kpc, args.v0err)

    # convert pm to rad/s
    args.pm = (u.mas.to(u.radian, args.pm)) / u.year.to(u.s, 1.0)
    if(args.pmerr != None):
        args.pmerr = (u.mas.to(u.radian, args.pmerr)) / u.year.to(u.s, 1.0)

    

    # Check whether all error bars required are given for distance and proper
    # motion.  If they are, go ahead with the random number generation thing.
    # If not, just do a quick straight calculation of the pbdot correction.
    
    if((args.pmerr != None) & (args.xpbdoterr != None) & (args.raerr != None) & (args.decerr != None)):
        
        # For each of the varying parameters (dist, pm, v0, R0), generate
        # args.niter gaussian deviates, with mean = value and sigma = error
        xpbdot = np.random.normal(args.xpbdot, args.xpbdoterr, args.niter)
        pm   = np.random.normal(args.pm, args.pmerr, args.niter)
        v0   = np.random.normal(args.v0, args.v0err, args.niter)
        R0   = np.random.normal(args.R0, args.R0err, args.niter)
        ra   = np.random.normal(ra_deg, raerr, args.niter)
        dec  = np.random.normal(dec_deg, decerr, args.niter)

        # Plot distributions for each parameter:
        fig_param = plt.figure()
        # xpbdot
        ax_param = fig_param.add_subplot(321)
        hist, bin_val = np.histogram(xpbdot, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        A, mu, sig = utils.fitgauss(bin_val, hist)
        #ax_param.plot(bin_val, utils.gaussian(bin_val, A, mu, sig), color='red')
        ax_param.set_xlabel('Difference in orbital decay ($10^{14}$ s/s)')
        
        # proper motion
        ax_param = fig_param.add_subplot(322)
        hist, bin_val = np.histogram(pm, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Proper motion (mas/yr)')
        
        # v0
        ax_param = fig_param.add_subplot(323)
        hist, bin_val = np.histogram(v0, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Solar velocity (km/s)')
        
        # R0
        ax_param = fig_param.add_subplot(324)
        hist, bin_val = np.histogram(R0, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Galactocentric distance (kpc)')
        
        # RA
        ax_param = fig_param.add_subplot(325)
        hist, bin_val = np.histogram(ra, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Right ascension (deg)')

        # Dec
        ax_param = fig_param.add_subplot(326)
        hist, bin_val = np.histogram(dec, bins=32)
        bin_size = (bin_val[1]-bin_val[0])
        bin_val = bin_val[0:len(bin_val)-1] + bin_size/2.
        ax_param.plot(bin_val, hist, linestyle='steps-mid', linewidth=1.4, 
                    color='blue')
        ax_param.set_xlabel('Declination (deg)')
 
        plt.show()
        
        
        # Run calcuation of pbdot correction through each iteration and build 
        # up array of values for pbdot correction
        #for i_iter in np.arange(args.niter):
        sys.stdout.write('\n {0:<7s}  {1:<7s}    {2:5s}   {3:<10s}\n'.format('Iter', '% done', 'n_bad', 'n_negative'))
        dist_array = []
        n_bad = 0
        n_negative = 0
        n_step_worked = 0
        for i_iter in np.arange(args.niter):
            c = coord.ICRSCoordinates(ra=ra[i_iter], dec=dec[i_iter], unit=(u.degree, u.degree))
            l = c.galactic.l.radians
            b = c.galactic.b.radians
            n_tries = 0
            sgn_guess = 1.0
            dist_guess = args.distguess
            dist_trial = -1.0
            try:
               # while(dist_trial < 0.001 and n_tries <= 25):             
               dist_trial = newton(pc.dist_func, dist_guess, 
                                    args=(l, b, 
                                    v0[i_iter], R0[i_iter], 
                                    pm[i_iter], args.pb, 
                                    xpbdot[i_iter]))
            #         n_tries += 1
            #        sgn_guess *= -1.0
            #        dist_guess = args.distguess + sgn_guess*np.ceil(float(n_tries)/2.0)*args.stepguess
                    
            except RuntimeError:
                n_bad+=1
            else:  # If there are no errors, then append to array
#                if(dist_trial > 0.000001 and n_tries <= 25):
                if(dist_trial > 0.001):
                    dist_array.append(dist_trial)
                #    if(n_tries > 1):
                #        n_step_worked += 1
                else:
                    n_negative += 1
                
            
#            dist_array.append(pc.get_dist_pbdot(ra[i_iter], dec[i_iter], 
#                                    pm[i_iter], xpbdot[i_iter],
#                                    args.pb,
#                                    v0=v0[i_iter],
#                                    R0=R0[i_iter],
#                                    dist_guess=args.distguess))
            if(np.fmod(i_iter, 16)==0):
                display_status(i_iter, args.niter, n_bad, n_negative)
                
        dist_array = np.array(dist_array)
        print '\nNumber of Runtime Errors = ', n_bad
        print 'Number of times stepping guess value worked = ', n_step_worked
        # Create histogram of pbdot correction, and fit a gsaussian to it.
        # Extract median (mean?) and sigma as our final value and error.
        if(args.distlim == None):
            distlim = (np.amin(dist_array)-0.1*np.abs(np.amin(dist_array)), 
                       np.amax(dist_array)+0.1*np.abs(np.amax(dist_array)))
        else:
            distlim = (args.distlim[0], args.distlim[1])
           
        dist_pdf, bin_val = np.histogram(dist_array, 
                      range=distlim, bins=96, density=True)
        bin_size = bin_val[1]-bin_val[0]
        dist_x = bin_val[0:len(bin_val)-1] + 0.5*bin_size
        #fig_pbdot = plt.figure()
        #ax_pbdot = fig_pbdot.add_axes([0.12, 0.1, 0.8, 0.85])
        #ax_pbdot.plot(pbdot_x, pbdot_pdf, linestyle='steps-mid', linewidth=1.4, 
        #            color='blue')
        
        # Get upper limit from distance distribution
#        dist_med, dist_min, dist_max = \
        dist_upper = get_pdf_prob(dist_x, dist_pdf, prob_intervals, 
                        norm=True, upper=True) #, \
        plot_pdf(dist_x, dist_pdf, \
###                     weights=alpha_weights, \
             xlabel='Distance (kpc)', \
             ylabel='Probability density',\
             prob_lines=dist_upper,\
#             prob_lines=np.append(dist_min, dist_max),\
             prob_linestyle=['dashed','dashdot','dotted'] #, \
                #                 'dashed','dashdot','dotted'], \
             )
             
        print " "
        print "DISTANCE (kpc): "
        print "  68%: < ", dist_upper[0]
        print "  95%: < ", dist_upper[1]
        print "  99%: < ", dist_upper[2]
        print " "
        
             
#        plt.axvline(dist_med)
#        dist_err_high = np.abs(dist_max-dist_med)
#        dist_err_low  = np.abs(dist_min-dist_med)
        
        #A, pbdot_corr, pbdot_corr_err = utils.fitgauss(pbdot_x, pbdot_pdf)
        #ax_pbdot.plot(pbdot_x, utils.gaussian(pbdot_x, A, pbdot_corr, pbdot_corr_err), 
                    #linestyle='solid', linewidth=1.4, color='red')
        plt.savefig('dist_pbdot.pdf')
        
        #print 'Corrected Pbdot = ', pbdot_corr, ' +/- ', pbdot_corr_err
#        print 'Distance (kpc) = ', dist_med, ' -', dist_err_low[0], \
#                                            ' +', dist_err_high[0]
#        print 'Distance peak = ', dist_x[np.argmax(dist_pdf)]
        
#        pbdot_new = args.pbdot + pbdot_corr_med
        #pbdot_corr_err_mean = np.mean(np.array([pbdot_corr_err_low[0],pbdot_corr_err_high[0]]))
        #pbdot_new_err = np.sqrt(args.pbdoterr**2.0 + pbdot_corr_err_mean**2.0)
#        pbdot_new_err_low = np.sqrt(args.pbdoterr**2.0 + pbdot_corr_err_low**2.0)
#        pbdot_new_err_high = np.sqrt(args.pbdoterr**2.0 + pbdot_corr_err_high**2.0)
        # print 'Pbdot new (median and mean error) = ', pbdot_new, ' +/- ', pbdot_new_err
        # print 'Pbdot new range (68%) = ', pbdot_corr_min[0], ' , ', pbdot_corr_max[0]
#        pbdot_new_mid = np.mean(np.array([pbdot_new - pbdot_new_err_low[0], pbdot_new + pbdot_new_err_high[0]]))
#        print 'Pbdot new (using midpoint between ranges) = ', pbdot_new_mid, ' - ', pbdot_new_mid - (pbdot_new-pbdot_new_err_low[0]), '   + ', (pbdot_new+pbdot_new_err_high[0])-pbdot_new_mid 
#        print 'Pbdot new (median corr and asymmetric error)   = ', pbdot_new, ' - ', pbdot_new_err_low[0], '   + ', pbdot_new_err_high[0]
        
    else:

        c = coord.ICRSCoordinates(ra=ra_deg, dec=dec_deg, 
                                unit=(u.degree, u.degree))
        l = c.galactic.l.radians
        b = c.galactic.b.radians
        
        dist_pbdot = newton(pc.dist_func, args.distguess, 
                                        args=(l, b, args.v0, args.R0, 
                                              args.pm, args.pb, args.xpbdot))
            
        #pbdot_new = args.pbdot - pbdot_corr
        
        print '\nWill not calculate final uncertainties...\n'
        print 'Distance (kpc) = ', dist_pbdot
        
        return
        
main()

