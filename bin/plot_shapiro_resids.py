#!/usr/local/bin/python

# Reads in residuals from ascii file ("print_resids" format for now), and
# plots them 

from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
from psrread import read_resid
from psrplot import plot_resid
from utils import exec_cmd
import argparse

tempo_resid_file = 'c.tmp'
tempo2_resid_file = 'tempo2_resid.dat'
tempo_info_file = 'info.tmp'

def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Plots residuals against time, orbital phase, or serial number.')
     
     parser.add_argument('timfile', 
                         nargs=1,
                         help='input TOA data file')
     parser.add_argument('-f', '--parfile', dest='parfile',
                         nargs=1,
                         help='input ephemeris file')
     parser.add_argument('-o', '--outfile', dest='outfile',
                         nargs='?',
                         const='shapiro_resid.png',
                         default=None,
                         help='output plot file name, with extension')
     parser.add_argument('--info',
                         nargs='*',
                         help='info numbers/values to plot based on info file (tempo1), or telescope IDs/flag values (tempo2)')
     parser.add_argument('--infoflag',
                         help='use this flag\'s arguments to for colour-coding with tempo2 residuals')
#     parser.add_argument('--xunits',
#                         nargs=1,
#                         default=['orbphase'],
#                         choices=['mjd', 'mjd0', 'year', 'orbphase', 'serial'],
#                         help='units of x axis')
     parser.add_argument('--yunits',
                         nargs=1,
                         default=['us'],
                         choices=['s', 'ms', 'us', 'ns'],
                         help='units of y axis')
     parser.add_argument('--xlim',
                         nargs=2,
                         type=float,
                         default=None,
                         help='limits of x axis in orbital phase')
     parser.add_argument('--ylim',
                          nargs=2,
                          type=float,
                          default=None,
                          help='limits of y axis in units given by yunits')
     parser.add_argument('--tempo2',
                          action='store_true',
                          default=False,
                          help='Residual file comes from tempo2 run')
     parser.add_argument('--binsize',
                         nargs=1,
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

     if(args.infoflag):
          args.infoflag = '-'+args.infoflag


     if(args.yunits): # true for now, for yunits being only one value
          args.yunits = args.yunits[0]

     print 'INFO = ', args.info

     # xlim in correct order
     if(args.xlim != None):
          try:
               xlim_backwards = (args.xlim[0] > args.xlim[1])
               if(xlim_backwards):
                    raise ArgError(progname, '--xlim', 0)
          except ArgError, err:
               print err.output
          if(args.xlim!=None):
               args.xlim=(args.xlim[0], args.xlim[1])

     # ylim in correct order
     if(args.ylim != None):
          try:
               ylim_backwards = (args.ylim[0] > args.ylim[1])
               if(ylim_backwards):
                    raise ArgError(progname, '--ylim', 0)
          except ArgError, err:
               print err.output
          if(args.ylim!=None):
               args.ylim=(args.ylim[0], args.ylim[1])

     if(args.binsize):
          print 'Binning data along x axis with bins of size ', args.binsize, 'in orbital phase.'
          args.binsize=args.binsize[0]

     return args



def main():

     progname = 'plot_shapiro_resids.py'
     args = get_opt(progname)


     # First, read in par file:
     f_par = open(args.parfile, 'r')
     par_data = [par_line.split() for par_line in f_par.readlines()]
     f_par.close()


     if(args.tempo2):
          resid_file = tempo2_resid_file
          infofile = None
     else:
          resid_file = tempo_resid_file
          infofile = tempo_info_file

     # Will run tempo/tempo2 3 times on best fit profile:
     
     # (1) Full solution with par file as is:
     if(args.tempo2):
          print 'par file = ', args.parfile
          tempo_command = 'tempo2 -output print_resids -f '+args.parfile+' '+\
                          args.timfile+' -file resid_output_format.dat '+ \
                          '-outfile '+tempo2_resid_file
          print 'tempo_command: ', tempo_command
          cmd_out = exec_cmd(tempo_command)
          print cmd_out
     else:
          tempo_command = 'tempo -f '+args.parfile+' '+args.timfile
          cmd_out = exec_cmd(tempo_command)
          exec_cmd('extract')
     # Read in residuals
     resid_data_allfit = read_resid(resid_file, 
                                    tempo2=args.tempo2, 
                                    info_file=infofile, 
                                    info_flag=args.infoflag)
     if(args.tempo2):
          if(args.info==None):
               resid_data_allfit['info']=None
               resid_data_allfit['info_val']=None
               resid_data_allfit['info_instr']=None
          elif(args.info==[]): # gave flag but no arguments
               args.info = resid_data_allfit['info_val']


     # (2) Get rid of SINI and M2 lines and fit everything else
     # Check as we go that SINI and M2 are in file
     temp_par_file = 'temp.par'
     f_temp_par = open(temp_par_file, 'w')
     for par_line in par_data:
          if((par_line[0] != 'SINI') & (par_line[0] != 'M2')):
               f_temp_par.write(' '.join(par_line)+'\n')
     f_temp_par.close()

     if(args.tempo2):
          tempo_command = 'tempo2 -output print_resids -f '+temp_par_file+' '+\
                          args.timfile+' -file resid_output_format.dat '+ \
                          '-outfile '+tempo2_resid_file
          cmd_out = exec_cmd(tempo_command)
     else:
          tempo_command = 'tempo -f '+temp_par_file+' '+args.timfile
          cmd_out = exec_cmd(tempo_command)
          exec_cmd('extract')
     # Read in residuals
     resid_data_orbfit = read_resid(resid_file, 
                                    tempo2=args.tempo2, 
                                    info_file=infofile, 
                                    info_flag=args.infoflag)
     if(args.tempo2):
          if(args.info==None):
               resid_data_orbfit['info']=None
               resid_data_orbfit['info_val']=None
               resid_data_orbfit['info_instr']=None
          elif(args.info==[]): # gave flag but no arguments
               args.info = resid_data_orbfit['info_val']
   

     # (3) Sam as (2), but turn off all fitting
     f_temp_par = open(temp_par_file, 'w')
     for par_line in par_data:
          if(len(par_line) > 2):
               if(par_line[2]=='1'):
                    par_line[2]='0'               
               if(par_line[0]=='JUMP'):
                    if(par_line[4]=='1'):
                         par_line[4]='0'
          if((par_line[0] != 'SINI') & (par_line[0] != 'M2')):
               f_temp_par.write(' '.join(par_line)+'\n')
     f_temp_par.close()
    
     if(args.tempo2):
          tempo_command = 'tempo2 -output print_resids -f '+temp_par_file+' '+ \
                          args.timfile+' -file resid_output_format.dat '+ \
                          '-outfile '+tempo2_resid_file
          exec_cmd(tempo_command)
     else:
          tempo_command = 'tempo -f '+temp_par_file+' '+args.timfile
          cmd_out = exec_cmd(tempo_command)
          exec_cmd('extract')
     # Read in residuals
     resid_data_nofit = read_resid(resid_file, 
                                   tempo2=args.tempo2, 
                                   info_file=infofile, 
                                   info_flag=args.infoflag)
     if(args.tempo2):
          if(args.info==None):
               resid_data_nofit['info']=None
               resid_data_nofit['info_val']=None
               resid_data_nofit['info_instr']=None
          elif(args.info==[]): # gave flag but no arguments
               args.info = resid_data_nofit['info_val']
     

     


# First, read in residuals data file, and assign each column to a separate
# numpy array     
#     print res_data['mjd']

#     if (len(argv) > 2):   # meaning 2nd argument is the desired output plot file name
#          plot_file = argv[2]
#     else:

     print "OUTFILE = ", args.outfile
     if(args.outfile == None):
          fig_size = (16, 18)
     else:
          fig_size = (14, 15)
          plot_file = args.outfile


               
     resid_data_list = [resid_data_allfit, resid_data_orbfit, resid_data_nofit]
     # resid_data_list = resid_data_allfit

     print 'xlim = ', args.xlim
     print 'ylim = ', args.ylim
     plot_resid(resid_data_list, info_plot=args.info, binsize=args.binsize,
                canvassize=fig_size, symsize=2.0,
                xunits='orbphase', yunits=args.yunits, 
                xticks=[True, False, False],
                xlabel=[True, False, False],
                ylabel=[False, True, False],
                xlim=args.xlim, ylim=args.ylim, gridlines=[0.],
                axislabelsize=32, ticklabelsize=32)

     if(args.outfile):
          plt.savefig(plot_file)
     else:
          plt.show()
     

main()
