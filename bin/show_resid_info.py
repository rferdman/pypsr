#!/usr/local/bin/python

from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
from psrread import read_resid, read_par
from psrinfo import get_resid_info
import argparse

def get_opt(progname):
     parser = argparse.ArgumentParser( 
          prog=progname,
          description='Plots residuals against time, orbital phase, or serial number.')
     
     parser.add_argument('resfile', 
                         nargs='?',
                         default='c.tmp',
                         help='Input residual data file')
     parser.add_argument('--parfile',
                         default=None,
                         help='Input parameter file')
     parser.add_argument('--infofile',
                         help='input info data file to distinguish data sets from tempo1 residuals')
     parser.add_argument('--infoflag',
                         help='use this flag\'s arguments for dividing up information from tempo2 residuals')
     parser.add_argument('--tempo2',
                          action='store_true',
                          default=False,
                          help='Residual file comes from tempo2 run')

     args=parser.parse_args()

     if(not args.infofile):
          args.infofile=None
     else:
          args.infofile = args.infofile  # otherwise stored as a list...

     if(args.infoflag):
          args.infoflag = '-'+args.infoflag



     return args


def main():

     progname = 'show_resid_info.py'
     args = get_opt(progname)

# First, read in residuals data file, and assign each column to a separate
# numpy array     
     resid_data = read_resid(args.resfile, 
                             tempo2=args.tempo2, info_file=args.infofile, 
                             info_flag=args.infoflag)


     # For now worry about ntoa in calculating chi2's but will implement 
     # par file reader to determine the number of DOF
     # ( = n_toa + n_free_param + 1 for fit for phase)
     if(args.parfile==None):
         n_param = 0
     else:
         if(args.tempo2):
             ffmt='tempo2'
         else:
             ffmt='tempo1'
         param_name, para_val, param_fit = read_par(args.parfile, file_format=ffmt)
         n_param = param_fit.count(True)
#         param_data = read_par(args.parfile, file_format=ffmt)



     # Now get information from residuals
     rinfo = get_resid_info(resid_data, nparam=n_param)
               
     # Now print out results, to stdout for now:
     print ''
     print 'Residual file:  ', args.resfile
     print 'Par file:       ', args.parfile
     print 'Number of TOAs: ', rinfo['ntoa']
     print 'Number of parameters:  ', rinfo['nparam'] 
     print 'Number of DOF:         ', rinfo['ndof']
     print '\n'
     
     print 'Info        Total   Avg weight  Number     Chi^2    Adjusted chi^2     rms         rms        MJD range     Years   Centre'
     print '            weight   per TOA    of TOAs   per TOA     (per DOF)     unweighted   weighted                            freq '
     print ''
     
     for i_info in range(len(resid_data['info_val'])):
         print '{0:8}    {1:7.5f}  {2:7.5f}  {3:7d}  {4:10.4f}  {5:10.4f}     {6:10.4f}  {7:10.4f}    {8:5d} - {9:5d}  {10:5.2f}    {11:6.1f}'.format(
             resid_data['info_val'][i_info], rinfo['normwgt'][i_info], 
             rinfo['avgwgt'][i_info], rinfo['npts'][i_info], 
             rinfo['rchi2'][i_info], rinfo['rchi2x'][i_info], 
             rinfo['resrms'][i_info], rinfo['resrmsw'][i_info], 
             int(rinfo['mjdstart'][i_info]), int(rinfo['mjdend'][i_info]), 
             (rinfo['mjdend'][i_info]-rinfo['mjdstart'][i_info])/365.25,
             rinfo['cfreq'][i_info] )

     print ''
     print '{0:8}    {1:7.5f}  {2:7.5f}  {3:7d}  {4:10.4f}  {5:10.4f}     {6:10.4f}  {7:10.4f}    {8:5d} - {9:5d}  {10:5.2f}    {11:6.1f}'.format(
         'Total', rinfo['sum_normwgt'], 
         rinfo['sum_avgwgt'], rinfo['sum_npts'], 
         rinfo['sum_rchi2'], rinfo['sum_rchi2x'], 
         rinfo['sum_resrms'], rinfo['sum_resrmsw'], 
         int(rinfo['sum_mjdstart']), int(rinfo['sum_mjdend']), 
         (rinfo['sum_mjdend']-rinfo['sum_mjdstart'])/365.25,
         rinfo['sum_cfreq'] )
     
     
main()

