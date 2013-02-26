#!/usr/local/bin/python

from sys import argv, exit
import glob
import numpy as np
import matplotlib.pyplot as plt
from psrread import *
from psrprof import *
from psrplot import plot_prof
import argparse
import mjd


#def get_opt(progname):
     #parse command line
#     parser = argparse.ArgumentParser(
#          prog=progname,
#          description='Makes a tiled or overlap plot of multiple lines.')

#     parser.add_argument('--offset',
#                         nargs='+',
#                         type=float,
#                         help='amount of offset to apply to each successive plot')
#     parser.add_argument('--outfile',
#                         nargs=1,
                       #  type=argparse.FileType('w'),
#                         default='plot_out.eps',
#                         help='Name of output file')

#     args=parser.parse_args()




def main():
# Input asp-style ascii profile will be first argument.
# (command-line options to come later)
     prog_name = argv[0].strip('.py')
     
#     arg = get_opt(prog_name)

#     prof_file = argv[1]

     input_files = []
     for file in argv[1:]:
# For some reason this works and ".appen()" doesn't:
          input_files[len(input_files):] = glob.glob(file)
#          input_files.append(glob.glob(file))

     input_files.reverse()
          
     prof_data = []
#     prof_date = []
     date_text = []
# First, read in profile data file, and assign each column to a separate
# numpy array     
#     prof_file = input_files[0]
#     prof_header = read_asc_header(prof_file)
#     prof_data[0] = read_asc_prof(prof_file)
      
     duty = get_duty('0737-3039A')

#     if(len(input_files) > 1):
#     for prof_file in input_files:
     for i_prof in np.arange(len(input_files)):
          prof_header = read_asc_header(input_files[i_prof])
          prof_data_temp = read_asc_prof(input_files[i_prof])
          prof_data_temp['i'] = norm(prof_data_temp['i'], duty)
          prof_data_temp['i'] = prof_data_temp['i'] + i_prof
          print "Index = ", i_prof, \
                   ", Min = ", np.min(prof_data_temp['i']), \
                   ", Max = ", np.max(prof_data_temp['i']) 
          prof_data.append(prof_data_temp)         
# Set up labelling for each profile:
          prof_date = mjd.mjdtodate(prof_header['imjd'], \
                                              dateformat='%Y %b %d')
          date_text.append((0.5, i_prof+0.16, prof_date))
          
     print "Date = ", date_text

          

# Do this just to make the ordering such that the first alphanumerically
# is at the top...
     # prof_data.reverse()
     # date_text.reverse()

     plot_prof(prof_data, yticks=False, canvassize=(8,10), \
                    hgrid=False, vgrid=False, \
                    ylim=(np.min(prof_data[0]['i'])-0.1, len(input_files)+0.1), \
                    figtext=date_text)



# meaningThe following means that the 2nd argument is the 
# desired output plot file name
     # plot_file = 'multi_profile.png'
     plot_file = 'multi_profile.eps'

     plt.savefig(plot_file)


main()
