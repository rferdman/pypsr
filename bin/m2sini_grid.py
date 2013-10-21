#!/usr/local/bin/python

# Running tempo over a grid of cosi and m2 values to eventually plot probability
# contours.

import sys
import subprocess
import numpy as np
from mjd import *
from scipy.optimize import curve_fit
# from lmfit import minimize, Parameters
from pdfplot import *
from utils import *
from mass_contours import *
import argparse


prob_intervals = np.array([0.683, 0.954, 0.9973])

def get_opt(progname):
    parser = argparse.ArgumentParser( 
            prog=progname,
            description='Plots residuals against time, orbital phase, or serial number.')
 
    parser.add_argument('timfile', 
                        nargs=1,
                        help='input TOA data file')
    parser.add_argument('-f', '--parfile', dest='parfile',
                        nargs=1,
                        help='input base ephemeris file')
    parser.add_argument('--psr', dest='psrname',
                        default=None,
                        help='Pulsar name.  Only used to create file names')
    parser.add_argument('--m2lim',
                        nargs=2,
                        type=float,
                        default=None,
                        help='limits of companion mass over which to create grid')
    parser.add_argument('--cosilim',
                        nargs=2,
                        type=float,
                        default=None,
                        help='limits of total mass over which to create grid')
    parser.add_argument('--nm2',
                        type=int,
                        default=32,
                        help='Number of m2 grid points')
    parser.add_argument('--ncosi',
                        type=int,
                        default=32,
                        help='Number of cosi grid points')
    parser.add_argument('--tempo2',
                        action='store_true',
                        default=False,
                        help='Parameter/TOA files are in tempo2 format')
    parser.add_argument('--m2gr',
                        type=float,
                        help='GR-derived comapnion mass')
    parser.add_argument('--binsize',
                        type=float,
                        help='Size in x axis units to bin data')


    args=parser.parse_args()

    # Error checking not done by argparse:
    if(args.parfile):
        args.parfile = args.parfile[0]  # otherwise stored as a list...
    else:
        print 'Must provide par file with --parfile option. Exiting...'
        exit()

    if(args.timfile):
        args.timfile = args.timfile[0]  # otherwise stored as a list...
    else:
        print 'Must provide TOA file with --timfile option. Exiting...'
        exit()

    # m2lim in correct order
    if(args.m2lim != None):
        try:
            m2lim_backwards = (args.m2lim[0] > args.m2lim[1])
            if(m2lim_backwards):
                raise ArgError(progname, '--m2lim', 0)
        except ArgError, err:
            print err.output
        if(args.m2lim!=None):
            args.m2lim=(args.m2lim[0], args.m2lim[1])
    else:
        args.m2lim = (0.81, 4.6) # Based on 1756-2251 for now as default

    # cosilim in correct order
    if(args.cosilim != None):
        try:
            cosilim_backwards = (args.cosilim[0] > args.cosilim[1])
            if(cosilim_backwards):
                raise ArgError(progname, '--cosilim', 0)
        except ArgError, err:
            print err.output
        if(args.cosilim!=None):
            args.cosilim=(args.cosilim[0], args.cosilim[1])
    else:
        args.cosilim = (0.19, 0.64)

    return args


def main():
 
    progname = 'm2mtot_grid.py'
    args = get_opt(progname)

    if(args.psrname == None):
        outfile_base = ''
    else:
        outfile_base = args.psrname
    
# Set par and tim files
#    par_base = '/Users/ferdman/Work/pulsar/1756-2251/timing/tempo/1756.dd.par.BASE'
#    tim_file = '/Users/ferdman/Work/pulsar/1756-2251/timing/tempo/1756.tempo.tim'

# Prepare par file for grid fitting
    # First, read in par file:
    par_base_contents = []
    f_par = open(args.parfile, 'r')
    for par_line in f_par.readlines():
        if ((par_line.split()[0] != 'SINI') & (par_line.split()[0]!='M2')):
            par_base_contents.append(par_line.split())
    #par_base_contents = [par_line.split() for par_line in f_par.readlines()]
    f_par.close()

    parfile_base = 'par_base.par'
    f_par_base = open(parfile_base, 'w')
    for par_line in par_base_contents:                
        f_par_base.write(' '.join(par_line)+'\n')
    f_par_base.close()    
    

    p_out =  grid_fit_shapiro_tempo(parfile_base, args.timfile, 
                                    m2_range=args.m2lim, 
                                    cosi_range=args.cosilim,
                                    n_m2=args.nm2, n_cosi=args.ncosi,
                                    tempo2=args.tempo2)


# Now make contour plots
    plot_contour_pdf(p_out['m2'], p_out['cosi'], np.transpose(p_out['norm_like']),
                     xlabel='Companion mass ($M_\\odot$)', 
                     ylabel='Cosine of inclination angle')
    plt.savefig('1756_m2cosi_contours.png')
    

# Now plot 1D pdfs for m2 and cosi
    m2_pdf = np.sum(p_out['norm_like'], axis=1)
    m2_med, m2_prob_min, m2_prob_max = \
        get_pdf_prob(p_out['m2'], m2_pdf, prob_intervals)
    plot_pdf(p_out['m2'], m2_pdf, 
             xlabel='Companion mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m2_prob_min, m2_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig('1756_m2_m2sini_pdf.png')
    print 'M2 = ', m2_med
    print '   68%: ', m2_prob_min[0], m2_prob_max[0]
    print '   95%: ', m2_prob_min[1], m2_prob_max[1]
    print '   99%: ', m2_prob_min[2], m2_prob_max[2]
    print ' '

    sini_pdf = np.sum(p_out['norm_like'], axis=0)
    sini_med, sini_prob_min, sini_prob_max = \
        get_pdf_prob(p_out['sini'], sini_pdf, prob_intervals)
    plot_pdf(p_out['sini'], sini_pdf, 
             xlabel='Sine of inclination angle', ylabel='Probability density',
             prob_lines=np.append(sini_prob_min, sini_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig('1756_sini_m2sini_pdf.png')
    print 'SINI = ', sini_med
    print '   68%: ', sini_prob_min[0], sini_prob_max[0]
    print '   95%: ', sini_prob_min[1], sini_prob_max[1]
    print '   99%: ', sini_prob_min[2], sini_prob_max[2]
    print ' '


    cosi_pdf = np.sum(p_out['norm_like'], axis=0)
    cosi_med, cosi_prob_min, cosi_prob_max = \
        get_pdf_prob(p_out['cosi'], cosi_pdf, prob_intervals)
    plot_pdf(p_out['cosi'], cosi_pdf, 
             xlabel='Total system mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(cosi_prob_min, cosi_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig('1756_cosi_m2sini_pdf.png')
    print 'COSI = ', cosi_med
    print '   68%: ', cosi_prob_min[0], cosi_prob_max[0]
    print '   95%: ', cosi_prob_min[1], cosi_prob_max[1]
    print '   99%: ', cosi_prob_min[2], cosi_prob_max[2]
    print ' '



main()
