#!/usr/local/bin/python2

# Running tempo over a grid of m1 and m2 values to eventually plot probability
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
    parser.add_argument('--m1lim',
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
    parser.add_argument('--nm1',
                        type=int,
                        default=32,
                        help='Number of m1 grid points')
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
    parser.add_argument('--mtotbins',
                        type=int,
                        default=16,
                        help='Number of histogram bins for mtot')
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

    # m1lim in correct order
    if(args.m1lim != None):
        try:
            m1lim_backwards = (args.m1lim[0] > args.m1lim[1])
            if(m1lim_backwards):
                raise ArgError(progname, '--m1lim', 0)
        except ArgError, err:
            print err.output
        if(args.m1lim!=None):
            args.m1lim=(args.m1lim[0], args.m1lim[1])
    else:
        args.m1lim = (2.569525, 2.57028)

    if(args.savefile & (args.loadfile!=None)):
        print 'Cannot load a file and save it at the same time!  Please choose one or the other.'
        print 'Exiting...'
        exit()

    return args


def main():
   
    progname = 'm1m2mtot_grid.py'
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

    parfile_base = 'par_base_m2mtot.par'
    f_par_base = open(parfile_base, 'w')
    for par_line in par_base_contents:                
        f_par_base.write(' '.join(par_line)+'\n')
    f_par_base.close()    
    
    

    if(args.loadfile==None):        
        p_out =  grid_fit_m1m2_tempo(parfile_base, args.timfile, 
                                       m1_range=args.m1lim, m2_range=args.m2lim,
                                       n_m1=args.nm1, n_m2=args.nm2, 
                                       tempo2=args.tempo2)
 
        if(args.savefile):
            save_file = 'm1m2mtot_grid_{0}_{1}'.format(args.nm1, args.nm2)
            save_array = np.array([p_out['m1'], p_out['m2'], 
                                p_out['mtot'], p_out['mtot_prob']]) 
            np.save(save_file+'_params', save_array)
            np.save(save_file+'_prob', p_out['norm_like'])

    else:
        load_array = np.load(args.loadfile+'_params.npy')
        # for the purposes of this routine, only need the following 
        # things in p_out
        p_out = {'m1':load_array[0],
                 'm2':load_array[1],
                 'mtot':load_array[2],
                 'mtot_prob':load_array[3]}
        p_out['norm_like'] = np.load(args.loadfile+'_prob.npy')

    
# Now make contour plots
    plot_contour_pdf(p_out['m1'], p_out['m2'], np.transpose(p_out['norm_like']),
                     xlim=args.m1lim, ylim=args.m2lim,
                     xlabel='Pulsar mass (M$_\\odot$)', 
                     ylabel='Companion mass (M$_\\odot$)')
                     
    ########
    if(False):
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
                                fontsize=12, \
                                horizontalalignment='right', \
                                verticalalignment='top',
                                rotation=88)
                                
        if(args.m2gr != None and args.mtotgr != None):
            plt.plot(args.m2gr, args.mtotgr, 'o', markersize=2.8, color='black')
            if(args.m2gr_err != None and args.mtotgr_err != None):
                plt.errorbar(args.m2gr, args.mtotgr, 
                            xerr=args.m2gr_err, yerr=args.mtotgr_err, 
                            ecolor='black')


    plt.savefig(outfile_base+'_m1m2_contours.'+args.plotformat)

    if(False):
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
    

    m1_pdf = np.sum(p_out['norm_like'], axis=0)
    m1_med, m1_prob_min, m1_prob_max = \
        get_pdf_prob(p_out['m1'], m1_pdf, prob_intervals)
    plot_pdf(p_out['m1'], m1_pdf, 
             xlabel='Pulsar mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m1_prob_min, m1_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_m1_m1m2_pdf.'+args.plotformat)
    print 'M1= ', m1_med
    print '   68%: ', m1_prob_min[0], m1_prob_max[0]
    print '   95%: ', m1_prob_min[1], m1_prob_max[1]
    print '   99%: ', m1_prob_min[2], m1_prob_max[2]
    print ' '

# Now plot 1D pdfs for m1 and m2:
    m2_pdf = np.sum(p_out['norm_like'], axis=1)
    m2_med, m2_prob_min, m2_prob_max = \
        get_pdf_prob(p_out['m2'], m2_pdf, prob_intervals)
    plot_pdf(p_out['m2'], m2_pdf, 
             xlabel='Companion mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m2_prob_min, m2_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_m2_m1m2_pdf.'+args.plotformat)
    print 'M2 = ', m2_med
    print '   68%: ', m2_prob_min[0], m2_prob_max[0]
    print '   95%: ', m2_prob_min[1], m2_prob_max[1]
    print '   99%: ', m2_prob_min[2], m2_prob_max[2]
    print ' '


# And now make weighted histogram of mtot values
    mtot_pdf, bin_edges = np.histogram(p_out['mtot'], args.mtotbins, 
                                     range=(np.amin(p_out['mtot']), np.amax(p_out['mtot'])),
                                     density=True, weights=p_out['mtot_prob'])
    bin_size = bin_edges[1] - bin_edges[0]
    mtot_x = bin_edges[0:len(bin_edges)-1] + 0.5*bin_size
    
    mtot_med, mtot_prob_min, mtot_prob_max = \
        get_pdf_prob(mtot_x, mtot_pdf, prob_intervals, norm=True)
    plot_pdf(mtot_x, mtot_pdf, 
             xlabel='Total system mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(mtot_prob_min, mtot_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'])
    plt.savefig(outfile_base+'_mtot_m1m2_pdf.'+args.plotformat)
    print 'MTOT = ', mtot_med
    print '   68%: ', mtot_prob_min[0], mtot_prob_max[0]
    print '   95%: ', mtot_prob_min[1], mtot_prob_max[1]
    print '   99%: ', mtot_prob_min[2], mtot_prob_max[2]
    print ' '


# Now make a plot of m1 and m2 PDFs on same plot
    plot_pdf(p_out['m1'], m1_pdf, 
             xlabel='Pulsar mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m1_prob_min, m1_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'], xlim=(1.10, 1.75), ylim=(0.00001,0.0065))

    plot_pdf(p_out['m2'], m2_pdf, 
             xlabel='Companion mass ($M_\\odot$)', ylabel='Probability density',
             prob_lines=np.append(m2_prob_min, m2_prob_max),
             prob_linestyle=['dashed','dashdot','dotted', 
                             'dashed','dashdot','dotted'], overplot=True)

    plt.savefig(outfile_base+'_m1m2_both_pdf_old.'+args.plotformat)


    fig = plt.figure(figsize=(16,5))
    ax = fig.add_axes([0.04, 0.15, 0.92, 0.81])
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    # Actually set xlabel manually so that companion and pulsar mass fall under respective PDFs
    # ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel("Probability density", fontsize=18)
    # Have no ticks or tick labels for the y-axis
    for tick in ax.yaxis.get_major_ticks():
        tick.label1On = False
        tick.label2On = False
        
    # Now plot m1 and m2 data
    plt.plot(p_out['m1'], m1_pdf, color='black', linestyle='steps-mid')
    plt.plot(p_out['m2'], m2_pdf, color='black', linestyle='steps-mid')
    pdf_append = np.append(m1_pdf,m2_pdf)
    mass_append = np.append(p_out['m1'], p_out['m2'])
    xmin = np.min(mass_append)
    xmax = np.max(mass_append)
    ymin = np.min(pdf_append)
    ymax = np.max(pdf_append)
    xspan = abs(xmax - xmin)
    yspan = abs(ymax - ymin)

    # Plot confidence intervals
    prob_linestyle=['dashed','dashdot','dotted', 
                    'dashed','dashdot','dotted']
    prob_lines_m1 = np.append(m1_prob_min, m1_prob_max)
    ax.vlines(prob_lines_m1[[0,2,3,5]], ymin-0.5, ymax+0.5,  # Include only 1 and 3-sig
               color="black", linestyle=prob_linestyle)
    prob_lines_m2 = np.append(m2_prob_min, m2_prob_max)
    ax.vlines(prob_lines_m2[[0,2,3,5]], ymin-0.5, ymax+0.5,  # Include only 1 and 3-sig
               color="black", linestyle=prob_linestyle)

    ax.fill_between(p_out['m1'], 0., m1_pdf, where=(p_out['m1']>=prob_lines_m1[0]) & (p_out['m1']<=prob_lines_m1[3]), facecolor='red', alpha=0.3) 
    #ax.fill_between(p_out['m1'], 0., m1_pdf, where=(p_out['m1']>=prob_lines_m1[1]) & (p_out['m1']<=prob_lines_m1[0]), facecolor='red', alpha=0.3) 
    ax.fill_between(p_out['m1'], 0., m1_pdf, where=(p_out['m1']>=prob_lines_m1[2]) & (p_out['m1']<=prob_lines_m1[0]), facecolor='yellow', alpha=0.3) 
    #ax.fill_between(p_out['m1'], 0., m1_pdf, where=(p_out['m1']>=prob_lines_m1[3]) & (p_out['m1']<=prob_lines_m1[4]), facecolor='red', alpha=0.3) 
    ax.fill_between(p_out['m1'], 0., m1_pdf, where=(p_out['m1']>=prob_lines_m1[3]) & (p_out['m1']<=prob_lines_m1[5]), facecolor='yellow', alpha=0.3) 

    ax.fill_between(p_out['m2'], 0., m2_pdf, where=(p_out['m2']>=prob_lines_m2[0]) & (p_out['m2']<=prob_lines_m2[3]), facecolor='red', alpha=0.3) 
    #ax.fill_between(p_out['m2'], 0., m2_pdf, where=(p_out['m2']>=prob_lines_m2[1]) & (p_out['m2']<=prob_lines_m2[0]), facecolor='red', alpha=0.3) 
    ax.fill_between(p_out['m2'], 0., m2_pdf, where=(p_out['m2']>=prob_lines_m2[2]) & (p_out['m2']<=prob_lines_m2[0]), facecolor='yellow', alpha=0.3) 
    #ax.fill_between(p_out['m2'], 0., m2_pdf, where=(p_out['m2']>=prob_lines_m2[3]) & (p_out['m2']<=prob_lines_m2[4]), facecolor='red', alpha=0.3) 
    ax.fill_between(p_out['m2'], 0., m2_pdf, where=(p_out['m2']>=prob_lines_m2[3]) & (p_out['m2']<=prob_lines_m2[5]), facecolor='yellow', alpha=0.3) 

    ax.set_xlim(xmin-0.02*xspan, xmax+0.02*xspan)
    ax.set_ylim(0., ymax + 0.02*yspan)
    
    
    

# Use this to manually do x-axis labels...
    ax.text(0.215, -0.12, 'Companion mass', fontsize=18, 
        horizontalalignment='center', verticalalignment='center', 
        transform=ax.transAxes)
# Use this to manually do x-axis labels...
    ax.text(0.82, -0.12, 'Pulsar mass', fontsize=18, 
        horizontalalignment='center', verticalalignment='center', 
        transform=ax.transAxes)

    plt.savefig(outfile_base+'_m1m2_both_pdf.'+args.plotformat)


main()
