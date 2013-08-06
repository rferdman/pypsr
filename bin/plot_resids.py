#!/usr/local/bin/python

# Reads in residuals from ascii file ("print_resids" format for now), and
# plots them 

from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
from psrread import read_resid
from psrplot import plot_resid
import argparse

resid_file = 'c.tmp'

def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Plots residuals against time, orbital phase, or serial number.')
     
     parser.add_argument('resfile', 
                         nargs='?',
                         default='c.tmp',
                         help='input residual data file')
     parser.add_argument('--infofile',
                         nargs=1,
                         help='input info data file to distinguish data sets from tempo1')
     parser.add_argument('--outfile',
                         nargs='?',
                         const='resid.png',
                         default=None,
                         help='output plot file name, with extension')
     parser.add_argument('--info',
                         nargs='*',
                         help='info numbers/values to plot based on info file (tempo1), or telescope IDs/flag values (tempo2)')
     parser.add_argument('--infoflag',
                         nargs=1,
                         help='use this flag\'s arguments to for colour-coding with tempo2 residuals')
     parser.add_argument('--xunits',
                         nargs=1,
                         default=['mjd'],
                         choices=['mjd', 'mjd0', 'year', 'orbphase', 'serial'],
                         help='units of x axis')
     parser.add_argument('--yunits',
                         nargs=1,
                         default=['us'],
                         choices=['s', 'ms', 'us', 'ns'],
                         help='units of y axis')
     parser.add_argument('--xlim',
                         nargs=2,
                         type=float,
                         default=None,
                         help='limits of x axis in units given by xunits')
     parser.add_argument('--ylim',
                          nargs=2,
                          type=float,
                          default=None,
                          help='limits of y axis in units given by yunits')
     parser.add_argument('--tempo2',
                          action='store_true',
                          default=False,
                          help='Residual file comes from tempo2 run')
     

     args=parser.parse_args()

     # Error checking not done by argparse:
     print 'infofile = ', args.infofile
     if(not args.infofile):
          args.infofile=None
          if(args.info):
               print 'No info file given.  Will ignore --info command line option.'
     else:
          args.infofile = args.infofile[0]  # otherwise stored as a list...

     if(args.infoflag):
          args.infoflag = '-'+args.infoflag[0]

     if(args.xunits): # true for now, for xunits being only one value
          args.xunits = args.xunits[0]

     if(args.yunits): # true for now, for yunits being only one value
          args.yunits = args.yunits[0]

     print 'INFO = ', args.info

     # xlim in correct order
     if(args.xlim):
          try:
               xlim_backwards = (args.xlim[0] > args.xlim[1])
               if(xlim_backwards):
                    raise ArgError(progname, '--xlim', 0)
          except ArgError, err:
               print err.output
          if(args.xlim!=None):
               args.xlim=(args.xlim[0], args.xlim[1])

     # ylim in correct order
     if(args.ylim):
          try:
               ylim_backwards = (args.ylim[0] > args.ylim[1])
               if(ylim_backwards):
                    raise ArgError(progname, '--ylim', 0)
          except ArgError, err:
               print err.output
          if(args.ylim!=None):
               args.ylim=(args.ylim[0], args.ylim[1])


     return args



def main():

     progname = 'plot_resids.py'
     args = get_opt(progname)

# Input print_resids-format data file will be first argument.
# (command-line options to come later)
#     resid_file = argv[1]

# First, read in residuals data file, and assign each column to a separate
# numpy array     
     resid_data = read_resid(args.resfile, 
                             tempo2=args.tempo2, info_file=args.infofile, info_flag=args.infoflag)
#     print res_data['mjd']

#     if (len(argv) > 2):   # meaning 2nd argument is the desired output plot file name
#          plot_file = argv[2]
#     else:

     print "OUTFILE = ", args.outfile
     if(args.outfile == None):
          fig_size = (16, 6)
     else:
          fig_size = (14, 5)
          plot_file = args.outfile


# If --info is not used and this is a tempo2 input file, then make info_plot==None
     if(args.tempo2):
          if(args.info==None):
               resid_data['info']=None
               resid_data['info_val']=None
               resid_data['info_instr']=None
          elif(args.info==[]): # gave flag but no arguments
               args.info = resid_data['info_val']
     else:
          if(args.info==[]):
               args.info=None

     print 'xlim = ', args.xlim
     plot_resid(resid_data, info_plot=args.info,
                canvassize=fig_size,
                xunits=args.xunits, yunits=args.yunits, 
                xlim=args.xlim,
                ylim=args.ylim)

     if(args.outfile):
          plt.savefig(plot_file)
     else:
          plt.show()
     

main()
