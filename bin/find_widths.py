#!/usr/local/bin/python

from sys import argv, exit
import glob
import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import mjd
from utils import *
from psrread import *
from psrprof import *
from psrplot import *
import argparse


xtick_scale = 10000.
xtick_offset = 0.6475

def tick_scale_formatter(x, pos):
     val_str = '{:.1f}'.format((x-xtick_offset)*xtick_scale)
     return val_str

# Determine the number of rows and columns for a subplot array 
# depending on the number of subplots N. We have it so that 
# n_cols >= n_rows.
def get_rows_cols(N):
     lower_square_root = int(np.floor(np.sqrt(N)))
     if(N - lower_square_root**2 > 0):
          if(N - lower_square_root**2 <= lower_square_root):
               n_rows = lower_square_root
               n_cols = lower_square_root + 1
          else:
               n_rows = lower_square_root + 1
               n_cols = lower_square_root + 1
     else: # perfect square
          n_rows = lower_square_root      
          n_cols = lower_square_root

     return n_rows, n_cols



def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Calculates widths for a set of pulse profiles, at a given fractional pulse height.')
     
     parser.add_argument('profiles',
                         nargs='+',
                         help='Input ascii profile files')
     parser.add_argument('--percent',
                         nargs='?',
                         type=float,
                         default=50.0,
                         help='Percent pulse height')
     parser.add_argument('--phasecut',
                         nargs='*',
                         type=float,
                         default=None,
                         help='Phase at which to cut profile for calculating width on separate components')
     parser.add_argument('--npolyfit',
                         nargs='?',
                         type=int,
                         default=12,
                         help='Number of points on each side of phase corresponding to chosen fractional pulse height to include in polynomial fitting, so there are a total of 2*npolyfit + 1 points used')
     parser.add_argument('--npolyorder',
                         nargs='?',
                         type=int,
                         default=None, # be sure to set this if None
                         help='Force the order of polynomial to attempt to fit to data points. Default is to use n_poly_fit-1 points')
     parser.add_argument('--niter',
                         nargs='?',
                         type=int,
                         default=32768,
                         help='Number of bootstrap iterations to run in width determination')
     parser.add_argument('--nhistbins',
                         nargs='?',
                         type=int,
                         default=32,
                         help='Number of bins to use in histogram for calculating widths and width errors')
     parser.add_argument('--ntestphase',
                         nargs='?',
                         type=int,
                         default=4,
                         help='Number of points to include on each side of estimated fractional pulse height phase when performing root-finding to determine exact phase value')
     parser.add_argument('--peakphase',
                         nargs='?',
                         type=float,
                         default=None,
                         help='Force phase of profile peak, for the purposes of calculating the phase corresponding to chosen fractional pulse height')
     
     args=parser.parse_args()


     # Fix up stuff that can't be done in (read: I haven't bothered to 
     # find out if I can do it in) add_argument:
     if(args.npolyorder==None):
          args.npolyorder = args.npolyfit-1

     if(args.phasecut==None):
          args.phasecut=[]
#     print 'args = ', args

     return args


def main():

     progname = 'find_widths.py'
     args = get_opt(progname)

# Cast command-line arguments into separate variables:
     in_files = args.profiles
     percent_height = args.percent
     phase_cut = args.phasecut
     n_poly_fit = args.npolyfit
     n_poly_order = args.npolyorder
     n_pts_omit = n_poly_fit - 1
     n_bootstrap = args.niter
     n_hist_bins = args.nhistbins
     n_pts_test = args.ntestphase
     peak_phase = args.peakphase
     

##     in_files = ['/Users/ferdman/Work/pulsar/0737-3039A/profs/add_campaigns/0737-3039A.20??.??.add.asc']
###     in_files = ['/Users/ferdman/Work/pulsar/1756-2251/profs/add_epochs/*.add.asc']
##     in_files = ['/Users/ferdman/Work/pulsar/1756-2251/profs/singleprofs/gasp/TemplateTestRot.asc']
##     percent_height = float(argv[1])
#    percent_height = 25.
#     phase_cut = [0.5]
##     phase_cut = []
    # number of points on *each* side of the profile in fit, so there are 
     # 2*n_poly_fit + 1 points to choose from
##     n_poly_fit = 12
##     n_poly_fit = 12  
     # number of iterations of fit
##     n_bootstrap = 32768
##     n_hist_bins = 48
##     n_hist_bins = 32
     # number of points to omit for each iteration of fit 
     # (= half total # of points -1)
##     n_pts_omit = n_poly_fit-1
     # number of points to use around expected phase to do root finding
     # for getting exact phase
##     n_pts_test = 4
     # order of polynomial in fit
##     n_poly_order = 6
##     n_poly_order = 10

#####     peak_phase = 0.00321002 # known from finding peak in high-S/N template which we assume is aligned

#     n_poly_fit = 16
#     n_bootstrap = 32768
#     n_hist_bins = 48
#     n_pts_omit = n_poly_fit
#     n_pts_test = 6
#     n_poly_order = 8
# n_poly_order = 2*n_poly_fit - n_pts_omit - 

##     n_subplot_rows = 3  # 0737A
##     n_subplot_cols = 5  # 0737A
     n_subplot_rows = 6
     n_subplot_cols = 6


#     matplotlib.rc('font', size=22)

     input_files = []
     for file_name in in_files:
# For some reason this works and ".append()" doesn't:
          input_files[len(input_files):] = glob.glob(file_name)

     n_files = len(input_files)

# Determine number of rows and columns for subplotting 
# bootstrap results histograms.  Always have it so n_cols >= n_rows
     n_subplot_rows, n_subplot_cols = get_rows_cols(n_files)

     print "n_profs = ", n_files
     print "Input files = ", input_files
     print " "
     print "Percent height = ", percent_height
     print "Num rows = ", n_subplot_rows
     print "Num cols = ", n_subplot_cols


     if(len(phase_cut) > 0):
          print "Will separate profiles into separate phase ranges before processing."
          if(len(phase_cut) ==1):
               print "   Profile will be split along phase ", phase_cut[0]
          else:
               print "   Profile will be split along phases ", phase_cut
          print " "


     mjd_width = np.zeros(len(input_files))
     width = np.zeros((len(input_files), len(phase_cut)+1))
     width_err = np.zeros((len(input_files), len(phase_cut)+1))
#     A = np.zeros((len(input_files), len(phase_cut)+1))
#     pdf_width = np.zeros((len(input_files), len(phase_cut)+1), n_hist_bins)
#     x_width = np.zeros((len(input_files), len(phase_cut)+1), n_hist_bins)


     for i_prof in np.arange(n_files):
          prof_data = read_asc_prof(input_files[i_prof])
          prof_head = read_asc_header(input_files[i_prof])
          psr_name = prof_head['psrname']
          mjd_width[i_prof] = prof_head['imjd'] + prof_head['smjd']/86400.
          n_bins = len(prof_data['i'])
          phase_cut_bins = np.append(0., phase_cut)*n_bins
          phase_cut_bins = np.append(phase_cut_bins, n_bins)
          phase_cut_bins = phase_cut_bins.astype(int)
          #print "phase_cut_bins = ", phase_cut_bins
          phase_cuts = phase_cut_bins.astype(float)/n_bins
       
          #print "arange = ", np.arange(len(phase_cut_bins) - 1)
#          print "File ", input_files[i_prof]
          print " "
          mjdtext = 'MJD {0:8.2f}'.format(mjd_width[i_prof])
          datetext = '{0:10s}'.format(mjd.mjdtodate(mjd_width[i_prof], dateformat="%Y-%b-%d"))
          print mjdtext, '  ', datetext
          
          prof_in = prof_data.copy()
#          x_peak, y_peak = get_peak(prof_in, n_pts_fit=10, n_order=19, n_test_fit=4)
#          plt.savefig('peak_fit.'+'{0:5.0f}'.format(mjd_width[i_prof])+'.png')
    
          for i_phase in np.arange(len(phase_cut_bins) - 1):
               
               print "   Phase range {0:.2f} to {1:.2f}".format(prof_data['phase'][phase_cut_bins[i_phase]], \
                                                                prof_data['phase'][phase_cut_bins[i_phase+1]-1])
               print " "
               prof_in['i'] = prof_data['i'][phase_cut_bins[i_phase]:phase_cut_bins[i_phase+1]-1]
               prof_in['phase'] = prof_data['phase'][phase_cut_bins[i_phase]:phase_cut_bins[i_phase+1]-1]
               
               width[i_prof,i_phase], width_err[i_prof,i_phase], A, pdf_width, x_width = \
                   get_width(prof_in, psr_name, percent_height, x_peak=peak_phase, \
                                  n_pts_fit=n_poly_fit, n_order=n_poly_order, \
                                  n_omit=n_pts_omit, n_test_fit=n_pts_test, \
                                  n_boot=n_bootstrap, hist_bins=n_hist_bins, \
                                  return_more=True)
               
               # Plot gaussian fit to bootstrap histograms as we go along as a grid of subplots:
               fig = plt.figure(i_phase)
               print 'N_ROWS = ', n_subplot_rows
               print 'N_COLS = ', n_subplot_cols
               print 'I_PROF = ', i_prof
               ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, i_prof+1)
#              plt.subplot(3, 4, i_prof)
               ax.plot(x_width*360., pdf_width, linestyle='steps-mid', linewidth=1.4, color='blue')
               ax.plot(x_width*360., gaussian(x_width, A, width[i_prof, i_phase], width_err[i_prof, i_phase]), color='red')
               # Set y limits here since we can't retrieve them later:
               ax.set_ylim(0., 1.1*np.amax(np.append(pdf_width, A)))
               # set the number of ticks, putting them at nice locations, pruning off edge labels if they 
               # happen on the plot corners:
               ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4, prune='both'))
               ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4, prune='both'))
               # No y ticklabels:
               ax.yaxis.set_ticklabels([])

               # Add date as text to the top left corner of each subplot
               ax.text(0.04, 0.96, datetext, horizontalalignment='left', verticalalignment='top', \
                            fontsize='6', transform=ax.transAxes)

               print "   Width = ", width[i_prof, i_phase], "+/-", width_err[i_prof, i_phase]
               print ""
               print ""
            
#        plt.show()

     # Set up formatter for scaling xtick labels:
     # major_formatter = FuncFormatter(tick_scale_formatter)
     print "len(phase_cut_bins) = ", len(phase_cut_bins)

#      print "phase_cuts = ", phase_cuts
     for i_phase in np.arange(len(phase_cut_bins) - 1):
          width_data = {'mjd':mjd_width, 'width':width[:, i_phase], 'werr':width_err[:, i_phase]}
          outfile_root = '{0}.w{1:.1f}.{2:3.1f}_{3:3.1f}'.format(psr_name, percent_height, phase_cuts[i_phase], phase_cuts[i_phase+1])
          # Plot widths we just calculated
          plot_widths(width_data, yunits='deg')
          plt.savefig(outfile_root+'.png')
          print np.std(width_data['width'])*360. # std dev in degrees
          print np.mean(width_data['werr'])*360.  # Mean error
          print np.median(width_data['werr'])*360.  # Mean error

          # Output plot of bootstrap histogram
          fig = plt.figure(i_phase)
          # Find overall lowest and highest x values, and calculate differnce
          x_lo = 360.0*np.amin((width_data['width'] - 6.0*width_data['werr']))
          x_hi = 360.0*np.amin((width_data['width'] + 6.0*width_data['werr']))
          x_span = np.abs(x_hi - x_lo)
          # Run through each profile and adjust limits to be consistent and show variation, and fix up labels
          for i_prof in np.arange(n_files):
               ax = fig.add_subplot(n_subplot_rows, n_subplot_cols, i_prof+1)
               # Centre on histogram, but set limits so that all subplots are on same scale:
               ax.set_xlim(360.*(width_data['width'][i_prof]) - 0.5*x_span, \
                                360.0*(width_data['width'][i_prof]) + 0.5*x_span)
               # Set text size for tick labels, and get rid of all labels other than lower x axis:
               ax.tick_params(axis='x', labelsize=8, top='off')
               ax.tick_params(axis='y', left='off', right='off')
               # Rescale tick labels to some reasonable numbers by removing the lowest x value of all width hists
               major_formatter = FuncFormatter(lambda x, pos: ('%.1f')%(x - x_lo))
               ax.xaxis.set_major_formatter(major_formatter)

#               x_tick_labels = ax.get_xticklabels()
#                    tick_label.set_label(str( float(tick_label.get_text())/x_base_val) ) 
#               x_tick_vals = [str( float(x_tick_labels[i_tick])/x_base_val ) for i_tick in len(x_ticks)]
#               ax.xticks(x_tick_locs, x_tick_labels)

#               ax.ticklabel_format(useOffset=0.6475, axis='x')
#               for tick in ax.xaxis.get_major_ticks():
#                    tick.label1.set_fontsize(8)
#               for tick in ax.yaxis.get_major_ticks():
#                    tick.label1.set_fontsize(8) 

          # Set common labels by setting text on current plot (easiest way I found) -- x label only here:
          fig.text(0.5, 0.04, \
                        '{0:2d}% profile width (- {1:.1f} degrees)'.format(int(percent_height), x_lo), \
                        ha='center', va='center')


          # Output bootstrap histogram plots:
          plt.savefig(outfile_root+'.boot_fit.png')

    # Now output widths to file.
          outfile = '{0}.w{1:.1f}.{2:3.1f}_{3:3.1f}'.format(psr_name, percent_height, phase_cuts[i_phase], phase_cuts[i_phase+1])+'.dat'
          f_out = open(outfile, 'w')
          np.savetxt(outfile_root+'.dat', \
                          np.transpose((width_data['mjd'], width_data['width'], width_data['werr'])), \
                          fmt='%-23.15f %14.10f   %14.10f') # x,y,z equal sized 1D arrays

        


main()
