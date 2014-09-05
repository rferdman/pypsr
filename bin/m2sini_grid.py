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
from psrread import read_par
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
    parser.add_argument('-s', '--save', dest='savefile',
                        action='store_true',
                        default=False,
                        help='Save results of grid into a python save file')
    parser.add_argument('-l', '--load', dest='loadfile',
                        help='Load results from file given as argument')
    parser.add_argument('--m1curve',
                        nargs='+',
                        type=float,
                        default=None,
                        help='Vales of pulsar mass in solar units to draw over contour plot')
    parser.add_argument('--m1bins',
                        type=float,
                        default=16,
                        help='Number of histogram bins for m1')
    parser.add_argument('--plotformat',
                        default='png',
                        help='File format for output plots')


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
        
        
    if(args.savefile & (args.loadfile!=None)):
        print 'Cannot load a file and save it at the same time!  Please choose one or the other.'
        print 'Exiting...'
        sys.exit()
        
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
    
    # read in pb and asini from par file, adn calculation mass function:
    if(args.tempo2):
        tempo_ver = 'tempo2'
    else:
        tempo_ver = 'tempo1'
    params=read_par(args.parfile, file_format=tempo_ver)
    # asini is in (light) seconds.  Convert pb to secs as well.
    # tsun is in us to give f_mass in solar units.
    tsun = 4.925490947
    tsun_secs = tsun*10**(-6)
    pb_secs = float(params['pb'])*86400.0
    f_mass = 4*np.pi*np.pi*(float(params['a1'])**3.)/(tsun_secs*(pb_secs**2.))
#    f_mass = 4*np.pi*np.pi*(float(params['a1'])**3.)/(tsun*(float(params['pb'])**2.))
    print 'pb = ', float(params['pb']), ' = ', pb_secs, ' sec'
    print 'a1 = ', float(params['a1'])
    print 'fmass = ', f_mass
    
    
    
    if(args.loadfile==None):        
        p_out =  grid_fit_shapiro_tempo(parfile_base, args.timfile, 
                                        m2_range=args.m2lim, 
                                        cosi_range=args.cosilim,
                                        n_m2=args.nm2, n_cosi=args.ncosi,
                                        fmass=f_mass,
                                        tempo2=args.tempo2)

        if(args.savefile):
            save_file = 'm2cosi_grid_{0}_{1}'.format(args.nm2, args.ncosi)
            save_array = np.array([p_out['m2'], p_out['cosi'], p_out['sini'], 
                                    p_out['m1'], p_out['m1_prob']])
            np.save(save_file+'_params', save_array)
            np.save(save_file+'_prob', p_out['norm_like'])
    else:
        load_array = np.load(args.loadfile+'_params.npy')
        # for the purposes of this routine, only need the following 
        # things in p_out
        p_out = {'m2':load_array[0],
                 'cosi':load_array[1],
                 'sini':load_array[2],
                 'm1':load_array[3],
                 'm1_prob':load_array[4]}
        p_out['norm_like'] = np.load(args.loadfile+'_prob.npy')

# Now make contour plots
    plot_contour_pdf(p_out['cosi'], p_out['m2'], p_out['norm_like'],
                         xlabel='|cos $i$|',
                         ylabel='Companion mass ($M_\\odot$)')
    # Add in m1 curves
    if(args.m1curve != None):
        for i_m1 in np.arange(len(args.m1curve)):
            sini_plot = ( (f_mass*(args.m1curve[i_m1]+p_out['m2'])**2.)**(1./3.) )/p_out['m2']
            cosi_plot = np.sqrt(1.0 - sini_plot**2.)
            plt.plot(cosi_plot, p_out['m2'], linestyle='dashed', color='black')
    
    plt.savefig('1756_m2cosi_contours.'+args.plotformat)
    

# Now plot 1D pdfs for m2 and cosi
    m2_pdf = np.sum(p_out['norm_like'], axis=1)
    m2_med, m2_prob_min, m2_prob_max = \
        get_pdf_prob(p_out['m2'], m2_pdf, prob_intervals)
    plot_pdf(p_out['m2'], m2_pdf, 
             xlabel='Companion mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m2_prob_min, m2_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig('1756_m2_m2sini_pdf.'+args.plotformat)
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
    plt.savefig('1756_sini_m2sini_pdf.'+args.plotformat)
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
    plt.savefig('1756_cosi_m2sini_pdf.'+args.plotformat)
    print 'COSI = ', cosi_med
    print '   68%: ', cosi_prob_min[0], cosi_prob_max[0]
    print '   95%: ', cosi_prob_min[1], cosi_prob_max[1]
    print '   99%: ', cosi_prob_min[2], cosi_prob_max[2]
    print ' '


    # Now deal with m1's:  create histogram weighted by likelihood
    m1_pdf, bin_edges = np.histogram(p_out['m1'], args.m1bins, 
                                    density=True, weights=p_out['m1_prob'])
    # We can define the bin centres as follows since our call to np/histogram gives 
    # back evenly spaced bins
    bin_size = bin_edges[1] - bin_edges[0]
    m1_val = bin_edges[0:len(bin_edges)-1] + 0.5*bin_size
    # Get PDF intervals and values:
    # pdf_rho = rho_hist/np.sum(rho_hist)
    m1_med, m1_prob_min, m1_prob_max = \
                        get_pdf_prob(m1_val, m1_pdf,prob_intervals, norm=True)
    plot_pdf(m1_val, m1_pdf, 
             xlabel='Pulsar mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m1_prob_min, m1_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig('1756_m1_m2sini_pdf.'+args.plotformat)

    print ' '
    print 'M1 = ', m1_med
    print '  68%: ', m1_prob_min[0], '  ', m1_prob_max[0]
    print '  95%: ', m1_prob_min[1], '  ', m1_prob_max[1]
    print '  99%: ', m1_prob_min[2], '  ', m1_prob_max[2]
    print ' '




main()
