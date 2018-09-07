#!/usr/local/bin/python2

# Reads in residuals from ascii file ("print_resids" format for now), and
# plots them 

from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
from psrplot import plot_m1m2
import argparse

def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Plots post-Keplerian parameter curves in pulsar mass vs companion mass space, assuming GR')
     
     parser.add_argument('m1m2file', 
                         nargs='?',
                         default='m1m2.dat',
                         help='input m1m2 data file')
     parser.add_argument('--m1gr',
                         type=float,
                         nargs=1,
                         help='Measured pulsar mass as found through other means, e.g. timing with DDGR binary model')
     parser.add_argument('--m2gr',
                         type=float,
                         nargs=1,
                         help='Measured companion mass as found through other means, e.g. timing with DDGR binary model')
     parser.add_argument('--m1gr_err',
                         type=float,
                         nargs=1,
                         help='Measured uncertainty on pulsar mass as found through other means, e.g. timing with DDGR binary model')
     parser.add_argument('--m2gr_err',
                         type=float,
                         nargs=1,
                         help='Measured uncertainty on companion mass as found through other means, e.g. timing with DDGR binary model')
     parser.add_argument('--m1m2contour', dest='m1m2_contour_file',
                         default=None,
                         help='Input numpy file for plotting m1-m2 1-sigma contour instead of simply error bars')
     parser.add_argument('--outfile',
                         nargs='?',
                         const='m1m2.png',
                         default=None,
                         help='Output plot file name, with extension')
     parser.add_argument('--pkparams',
                         nargs='*',
                         default=['omdot', 'gamma'],
                         choices=['omdot', 'gamma', 'pbdot', 'r', 's'],
                         help='Post-Keplerian parameters to plot')
     parser.add_argument('--pkcoords',
                         nargs='*',
                         type=float,
                         default=None,
                         help='Label post-Keplerian parameters on plot at these coordinates (given is same order as --pkparams argument)')
     parser.add_argument('--plot_sin1',
                         action='store_true',
                         default=False,
                         help='Plot area restricted by sin(i)<=1 condition')
     parser.add_argument('--parfile', 
                         nargs='?',
                         default=None,
                         help='input pulsar ephemeris file, needed for plotting sin(i) restriction')
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
     parser.add_argument('--inset',
                          action='store_true',
                          default=False,
                          help='Residual file comes from tempo2 run')


     args=parser.parse_args()

     # Error checking not done by argparse:

     if(args.pkparams==None): # true for now, for yunits being only one value
          print 'ERROR: Must provide post-Keplerian parameters for plotting.'
          exit()

     if(args.pkcoords):
          print 'No. of pkparams: ', len(args.pkparams)
          print 'No. of pkcoords: ', len(args.pkcoords)
          if(len(args.pkcoords) != 2*len(args.pkparams)):
               print('ERRO: Must give an x and y value for each PK param plotted.')
               exit()
          else:
               pkcoords = []
               # Organize into tuples:
               for i_coord in range(0, len(args.pkcoords), 2):
                    pkcoords.append((args.pkcoords[i_coord], args.pkcoords[i_coord+1]))
               args.pkcoords=pkcoords
               print args.pkcoords

#     if(args.m1m2_contour_file!=None):
#         if(args.m1gr!=None | args.m2gr!=None):
#             print 'Cannot plot both m1/m2 GR values AND m1/m2 contours'
#             exit()

     if(args.m1gr):
          args.m1gr = args.m1gr[0]
          if(args.m1gr_err==None):
               print 'WARNING: No error to m1gr has been provided.'
          else:
               args.m1gr_err = args.m1gr_err[0]
          if(args.m2gr==None):
               print 'ERROR: Must provide an --m2gr to correspond to the --m1gr argument.'
               exit()
          else:
               args.m2gr = args.m2gr[0]
               if(args.m2gr_err==None):
                    print 'WARNING: No error to m2gr has been provided.'
               else: 
                    args.m2gr_err = args.m2gr_err[0]
                
     if(args.plot_sin1):
          if(args.parfile==None):
              print 'Plotting sin(i) restricted region requires a par file.'
              exit()
          

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

     progname = 'plot_m1m2.py'
     args = get_opt(progname)


     plot_m1m2(args.m1m2file, plot_pk=args.pkparams, 
               m1gr=args.m1gr, m1gr_err=args.m1gr_err, 
               m2gr=args.m2gr, m2gr_err=args.m2gr_err, 
               m1m2_contour=args.m1m2_contour_file,
               plot_inset=args.inset, xlim=args.xlim, ylim=args.ylim,
               plot_sin1=args.plot_sin1, parfile=args.parfile,
               pk_label_coord=args.pkcoords)

     if(args.outfile):
          print 'Plotted to file ', args.outfile
          plt.savefig(args.outfile)
     else:
          plt.show()
     

main()
