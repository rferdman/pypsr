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

def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Calculates peak height and phase for input pulse profiles.')

     parser.add_argument('profiles',
                         nargs='+',
                         help='Input ascii profile files')
     parser.add_argument('--phasecut',
                         nargs='*',
                         type=float,
                         default=None,
                         help='Phase at which to cut profile for calculating width on separate components')
     parser.add_argument('-v', '--verbose', dest='verbose',
                         action='store_true',
                         default=False,
                         help='Verbose output mode')

     args=parser.parse_args()

     if(args.phasecut==None):
          args.phasecut = []

     return args



def tick_scale_formatter(x, pos):
     val_str = '{:.1f}'.format((x-xtick_offset)*xtick_scale)
     return val_str

def main():
     progname = 'find_peak.py'
     args = get_opt(progname)
#     in_files = ['/Users/ferdman/Work/pulsar/0737-3039A/profs/add_campaigns/0737-3039A.20??.??.add.asc']
#     in_files = ['/Users/ferdman/Work/pulsar/1756-2251/profs/add_epochs/*.add.asc']
#     in_files = ['/Users/ferdman/Work/pulsar/1756-2251/profs/template/gasp/1756-2251.1400.2048_bins.std']
#     in_files = ['/Users/ferdman/Work/pulsar/1756-2251/profs/singleprofs/gasp/TemplateTestRot.asc']
     #phase_cut = []
    # number of points on *each* side of the profile in fit, so there are 
     # 2*n_poly_fit + 1 points to choose from

#     matplotlib.rc('font', size=22)

     input_files = []
     for file_name in args.profiles:
# For some reason this works and ".append()" doesn't:
          input_files[len(input_files):] = glob.glob(file_name)

     n_files = len(input_files)
     print "Number of input profiles: ", n_files, '\n'
     # print "Input files = ", input_files

     if(len(args.phasecut) > 0):
          print "Will separate profiles into separate phase ranges before processing."
          if(len(args.phasecut) ==1):
               print "   Profile will be split along phase ", args.phasecut[0]
          else:
               print "   Profile will be split along phases ", args.phasecut
          print " "


     mjd_obs = np.zeros(len(input_files))
#     A = np.zeros((len(input_files), len(phase_cut)+1))
#     pdf_width = np.zeros((len(input_files), len(phase_cut)+1), n_hist_bins)
#     x_width = np.zeros((len(input_files), len(phase_cut)+1), n_hist_bins)


     for i_prof in np.arange(n_files):
          prof_data = read_asc_prof(input_files[i_prof])
          prof_head = read_asc_header(input_files[i_prof])
          psr_name = prof_head['psrname']
          mjd_obs[i_prof] = prof_head['imjd'] + prof_head['smjd']/86400.
          n_bins = len(prof_data['i'])
          phase_cut_bins = np.append(0., args.phasecut)*n_bins
          phase_cut_bins = np.append(phase_cut_bins, n_bins)
          phase_cut_bins = phase_cut_bins.astype(int)
          #print "phase_cut_bins = ", phase_cut_bins
          phase_cuts = phase_cut_bins.astype(float)/n_bins
       
          #print "arange = ", np.arange(len(phase_cut_bins) - 1)
#          print "File ", input_files[i_prof]
          print " "
          mjdtext = 'MJD {0:8.2f}'.format(mjd_obs[i_prof])
          datetext = '{0:10s}'.format(mjd.mjdtodate(mjd_obs[i_prof], dateformat="%Y-%b-%d"))
          print mjdtext, '  ', datetext
          
          prof_in = prof_data.copy()

    
          for i_phase in np.arange(len(phase_cut_bins) - 1):
               
               plt.figure()
               print "   Phase range {0:.2f} to {1:.2f}".format(prof_data['phase'][phase_cut_bins[i_phase]], \
                                                                prof_data['phase'][phase_cut_bins[i_phase+1]-1])
               print " "
               prof_in['i'] = prof_data['i'][phase_cut_bins[i_phase]:phase_cut_bins[i_phase+1]-1]
               prof_in['phase'] = prof_data['phase'][phase_cut_bins[i_phase]:phase_cut_bins[i_phase+1]-1]
               
               x_peak, y_peak = get_peak(prof_in, n_pts_fit=12, n_order=5, n_test_fit=4)
               if(len(args.phasecut) > 0):
                    plt.savefig('peak_fit.'+'{0:5.0f}'.format(mjd_obs[i_prof])+
                                '.phase_{0:d}'.format(i_phase)+'.png')
               else:
                    plt.savefig('peak_fit.'+'{0:5.0f}'.format(mjd_obs[i_prof])+
                                '.png')
               if(x_peak < 0):
                    print 'Could not choose one.  Exiting: '
                    exit()


               print 'Peak phase: ', x_peak[0]
               print 'Peak height: ', y_peak[0]
    # Now output widths to file.
###          outfile = '{0}.w{1:.1f}.{2:3.1f}_{3:3.1f}'.format(psr_name, percent_height, phase_cuts[i_phase], phase_cuts[i_phase+1])+'.dat'
###          f_out = open(outfile, 'w')
###          np.savetxt(outfile_root+'.dat', \
###                          np.transpose((width_data['mjd'], width_data['width'], width_data['werr'])), \
###                          fmt='%-23.15f %14.10f   %14.10f') # x,y,z equal sized 1D arrays

        


main()
