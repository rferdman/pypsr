#!/usr/local/bin/python

from sys import argv, exit
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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

# Scale factor for plotting difference plots
     scale = 2.0
     
#     arg = get_opt(prog_name)

#     prof_file = argv[1]
#     ref_file = argv[1]
     ref_file = argv[len(argv)-1]

     input_files = []
     for file in argv[1:]:
# For some reason this works and ".append()" doesn't:
          input_files[len(input_files):] = glob.glob(file)
#          input_files.append(glob.glob(file))

     input_files.reverse()

# Remove reference file from list if it is in the plotting list
     if(input_files.count(ref_file) > 0):
          input_files.remove(ref_file)
   
     duty = get_duty('0737-3039A')
# Read in reference profile
     ref_header = read_asc_header(ref_file)
     ref_data = read_asc_prof(ref_file)
     # ref_data['i'] = norm(ref_data['i'], duty)
    
# Make first prof in list for plot to be the ref profile
     prof_data = [ref_data]
     ref_date = mjd.mjdtodate(ref_header['imjd'], \
                                   dateformat='%Y %b %d')
     date_text = [(0.5, 0.20, ref_date)]

     nobs = []      
# Now calculate difference profiles and append to plotting list
     for i_prof in np.arange(len(input_files)):
          prof_header = read_asc_header(input_files[i_prof])
          nobs.append(prof_header['obscode'])
          print 'NOBS = ', nobs
          prof_data_temp = read_asc_prof(input_files[i_prof])
          # prof_data_temp['i'] = norm(prof_data_temp['i'], duty)
          diff_prof = remove_base((prof_data_temp['i'] - ref_data['i']), duty)
          prof_data_temp['i'] = scale*(diff_prof) + i_prof + 1.2 
          print "Index = ", i_prof
          prof_data.append(prof_data_temp)         
# Set up labelling for each profile:
          prof_date = mjd.mjdtodate(prof_header['imjd'], \
                                              dateformat='%Y %b %d')
          date_text.append((0.5, i_prof+1.4+scale*0.02, prof_date))
          
     print "Date = ", date_text

     nobs_unique = list(set(nobs))
     clr=[]
     for i_nobs in range(len(nobs)):
         if(nobs_unique <= 1):
             clr.append('black')
         else:
             clr.append(cm.gist_heat(float(nobs_unique.index(nobs[i_nobs]))/float(len(nobs_unique))))
          

# Do this just to make the ordering such that the first alphanumerically
# is at the top...
#      prof_data.reverse()
#      date_text.reverse()

     plot_prof(prof_data, yticks=False, canvassize=(8,10), vgrid=False, \
                    ylim=(np.min(prof_data[0]['i'])-0.1, len(input_files)+1 +0.1), \
                    figtext=date_text, linecolour=clr)



# meaningThe following means that the 2nd argument is the 
# desired output plot file name
     # plot_file = 'diff_profile.png'
     plot_file = 'diff_profile.pdf'

     plt.savefig(plot_file)


main()
