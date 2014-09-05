#!/usr/local/bin/python

# Running tempo over a grid of cosi and m2 values to eventually plot probability
# contours.

import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
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
    parser.add_argument('--mtotlim',
                        nargs=2,
                        type=float,
                        default=None,
                        help='limits of total mass over which to create grid')
    parser.add_argument('--m2gr',
                        type=float,
                        help='Measured companion mass as found through other means, e.g. timing with DDGR binary model')
    parser.add_argument('--mtotgr',
                        type=float,
                        help='Measured total mass as found through other means, e.g. timing with DDGR binary model')
    parser.add_argument('--m2gr_err',
                        type=float,
                        help='Measured uncertainty on companion mass as found through other means, e.g. timing with DDGR binary model')
    parser.add_argument('--mtotgr_err',
                        type=float,
                        help='Measured uncertainty on total mass as found through other means, e.g. timing with DDGR binary model')
    parser.add_argument('--nm2',
                        type=int,
                        default=32,
                        help='Number of m2 grid points')
    parser.add_argument('--nmtot',
                        type=int,
                        default=32,
                        help='Number of mtot grid points')
    parser.add_argument('--tempo2',
                        action='store_true',
                        default=False,
                        help='Parameter/TOA files are in tempo2 format')
    parser.add_argument('--binsize',
                        nargs=1,
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
        sys.exit()

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
        args.m2lim = (1.203, 1.256) # Based on 1756-2251 for now as default

    # mtotlim in correct order
    if(args.mtotlim != None):
        try:
            mtotlim_backwards = (args.mtotlim[0] > args.mtotlim[1])
            if(mtotlim_backwards):
                raise ArgError(progname, '--mtotlim', 0)
        except ArgError, err:
            print err.output
        if(args.mtotlim!=None):
            args.mtotlim=(args.mtotlim[0], args.mtotlim[1])
    else:
        args.mtotlim = (2.569525, 2.57028)

    if(args.savefile & (args.loadfile!=None)):
        print 'Cannot load a file and save it at the same time!  Please choose one or the other.'
        print 'Exiting...'
        exit()

    return args


def main():
   
    progname = 'm2mtot_grid.py'
    args = get_opt(progname)

    rc('font',**{'family':'sans-serif','sans-serif':['sans-serif']})
    #rc('text', usetex=True)
    
# Set par and tim files
    #par_base = '/Users/ferdman/Work/pulsar/1756-2251/timing/tempo/1756.ddgr.par.BASE'
    #tim_file = '/Users/ferdman/Work/pulsar/1756-2251/timing/tempo/1756.tempo.tim'
    
    if(args.psrname == None):
        outfile_base = ''
    else:
        outfile_base = args.psrname
        
# Prepare par file for grid fitting
    # First, read in par file:
    par_base_contents = []
    f_par = open(args.parfile, 'r')
    for par_line in f_par.readlines():
        if ((par_line.split()[0] != 'MTOT') & (par_line.split()[0]!='M2')):
            par_base_contents.append(par_line.split())
    #par_base_contents = [par_line.split() for par_line in f_par.readlines()]
    f_par.close()

    parfile_base = 'par_base.par'
    f_par_base = open(parfile_base, 'w')
    for par_line in par_base_contents:                
        f_par_base.write(' '.join(par_line)+'\n')
    f_par_base.close()    
    
    

    if(args.loadfile==None):        
        p_out =  grid_fit_m2mtot_tempo(parfile_base, args.timfile, 
                                       m2_range=args.m2lim, mtot_range=args.mtotlim,
                                       n_m2=args.nm2, n_mtot=args.nmtot, 
                                       tempo2=args.tempo2)
 
        if(args.savefile):
            save_file = 'm2mtot_grid_{0}_{1}'.format(args.nm2, args.nmtot)
            save_array = np.array([p_out['m2'], p_out['mtot'], 
                                p_out['m1'], p_out['m1_prob']]) 
            np.save(save_file+'_params', save_array)
            np.save(save_file+'_prob', p_out['norm_like'])

    else:
        load_array = np.load(args.loadfile+'_params.npy')
        # for the purposes of this routine, only need the following 
        # things in p_out
        p_out = {'m2':load_array[0],
                 'mtot':load_array[1],
                 'm1':load_array[2],
                 'm1_prob':load_array[3]}
        p_out['norm_like'] = np.load(args.loadfile+'_prob.npy')

    
# Now make contour plots
    plot_contour_pdf(p_out['m2'], p_out['mtot'], np.transpose(p_out['norm_like']),
                     xlim=args.m2lim, ylim=args.mtotlim,
                     xlabel='Companion mass (M$_\\odot$)', 
                     ylabel='Total system mass (M$_\\odot$)')
    # Add in m1 curves
    if(args.m1curve != None):
        for i_m1 in np.arange(len(args.m1curve)):
            m2_plot = p_out['mtot'] - args.m1curve[i_m1]
            plt.plot(m2_plot, p_out['mtot'], linestyle='dashed', color='black')
            # annotate
            y_text = args.mtotlim[1] - 0.015*abs(args.mtotlim[1] - args.mtotlim[0])
            x_text = y_text - args.m1curve[i_m1] - 0.005*abs(args.m2lim[1] - args.m2lim[0])
            text_str = '$\\mathsf{m_p} = $ '+str(args.m1curve[i_m1])+' M$_\\odot$'
            plt.text(x_text, y_text, text_str, 
                            fontsize=16, \
                            horizontalalignment='right', \
                            verticalalignment='top',
                            rotation=88)
                            
    if(args.m2gr != None and args.mtotgr != None):
        plt.plot(args.m2gr, args.mtotgr, 'o', markersize=2.8, color='black')
        if(args.m2gr_err != None and args.mtotgr_err != None):
            plt.errorbar(args.m2gr, args.mtotgr, 
                        xerr=args.m2gr_err, yerr=args.mtotgr_err, 
                        ecolor='black')


    plt.savefig(outfile_base+'_m2mtot_contours.'+args.plotformat)


# Try making a m1-m2 contour plot...
    m1_contours = p_out['mtot']-p_out['m2']
    m1lim = (args.mtotlim[0] - args.m2lim[0], args.mtotlim[1] - args.m2lim[1])
    plot_contour_pdf(m1_contours, p_out['m2'], p_out['norm_like'],
                     xlim=m1lim, ylim=args.m2lim, 
                     ylabel='Companion mass (M$_\\odot$)', 
                     xlabel='Pulsar mass (M$_\\odot$)')
                            
#    if(args.m2gr != None and args.mtotgr != None):
#        plt.plot(args.m2gr, args.mtotgr, 'o', markersize=2.8, color='black')
#        if(args.m2gr_err != None and args.mtotgr_err != None):
#            plt.errorbar(args.m2gr, args.mtotgr, 
#                        xerr=args.m2gr_err, yerr=args.mtotgr_err, 
#                        ecolor='black')


    plt.savefig(outfile_base+'_m2m1_contours.'+args.plotformat)
    

# Now plot 1D pdfs for m2 and mtot
    m2_pdf = np.sum(p_out['norm_like'], axis=1)
    m2_med, m2_prob_min, m2_prob_max = \
        get_pdf_prob(p_out['m2'], m2_pdf, prob_intervals)
    plot_pdf(p_out['m2'], m2_pdf, 
             xlabel='Companion mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m2_prob_min, m2_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_m2_m2mtot_pdf.'+args.plotformat)
    print 'M2 = ', m2_med
    print '   68%: ', m2_prob_min[0], m2_prob_max[0]
    print '   95%: ', m2_prob_min[1], m2_prob_max[1]
    print '   99%: ', m2_prob_min[2], m2_prob_max[2]
    print ' '

    mtot_pdf = np.sum(p_out['norm_like'], axis=0)
    mtot_med, mtot_prob_min, mtot_prob_max = \
        get_pdf_prob(p_out['mtot'], mtot_pdf, prob_intervals)
    plot_pdf(p_out['mtot'], mtot_pdf, 
             xlabel='Total system mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(mtot_prob_min, mtot_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_mtot_m2mtot_pdf.'+args.plotformat)
    print 'MTOT = ', mtot_med
    print '   68%: ', mtot_prob_min[0], mtot_prob_max[0]
    print '   95%: ', mtot_prob_min[1], mtot_prob_max[1]
    print '   99%: ', mtot_prob_min[2], mtot_prob_max[2]
    print ' '

# And now make weighted histogram of m1 values
    m1_pdf, bin_edges = np.histogram(p_out['m1'], args.m1bins, 
                                     range=(np.amin(p_out['m1']), np.amax(p_out['m1'])),
                                     density=True, weights=p_out['m1_prob'])
    bin_size = bin_edges[1] - bin_edges[0]
    m1_x = bin_edges[0:len(bin_edges)-1] + 0.5*bin_size
    
    m1_med, m1_prob_min, m1_prob_max = \
        get_pdf_prob(m1_x, m1_pdf, prob_intervals, norm=True)
    plot_pdf(m1_x, m1_pdf, 
             xlabel='Pulsar mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m1_prob_min, m1_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_m1_m2mtot_pdf.'+args.plotformat)
    print 'M1 = ', m1_med
    print '   68%: ', m1_prob_min[0], m1_prob_max[0]
    print '   95%: ', m1_prob_min[1], m1_prob_max[1]
    print '   99%: ', m1_prob_min[2], m1_prob_max[2]
    print ' '

main()
