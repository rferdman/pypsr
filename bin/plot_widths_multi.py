#!/usr/local/bin/python

# Reads in widths from ascii file and plots them 

from sys import argv, exit
from psrread import read_widths
from psrplot import plot_widths
import numpy as np
import mjd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, MultipleLocator, FormatStrFormatter
import glob

def main():
     input_files = argv[1:]

# Parse command line to get all width files we wish to plot
     width_file = []
     for file_name in input_files:
          # For some reason this works and ".append()" doesn't:
          width_file[len(width_file):] = glob.glob(file_name)

     n_subplots = len(width_file)
     print "N_SUBPLOTS = ", n_subplots
     
# Set up the plot:
     fig = plt.figure()#figsize=canvassize)

     fig_top = 0.95
     fig_bottom = 0.13
     fig_left = 0.13
     fig_right = 0.95

     for i_w in np.arange(n_subplots):
          print "I_W = ", i_w
#          width_data =[]
#          for wfile in width_files:
          width_data = read_widths(width_file[i_w])
#     print res_data['mjd']

# Now plot width in degrees vs time in years for each subplot:
#          for width in width_data:
                    # date_out = [mjd.mjdtodate(m) for m in width['mjd']]
          width_data['year'] = [mjd.mjdtoyear(m) for m in width_data['mjd']]
          width_data['width'] *= 360.
          width_data['werr']  *= 360.
                    # width['year'] = [d.year + d.day/365. + \
                        #                    d.hour/(365.*24.) + \
                        #                    d.minute/(365.*24.*60.) + \
                        #                    d.second/(365.*24.*60.*60.) \
               #                    for d in date_out]
# Set up plot limits now
          xmin = np.amin(width_data['year'])-0.25
          xmax = np.amax(width_data['year'])+0.25
          ymin = np.amin(width_data['width'] - width_data['werr'])
          ymax = np.amax(width_data['width'] + width_data['werr'])
          xspan = abs(xmax - xmin)
          yspan = abs(ymax - ymin)
          
#          ax = fig.add_subplot(n_subplots, 1, i_w+1)
          if(i_w == 0):
               max_yspan = yspan
          else:
# Keep track of max yspan to later set all plots to have same scale
               if(yspan > max_yspan):
                    max_yspan = yspan
          ax = fig.add_axes([fig_left, \
                    fig_bottom+(fig_top-fig_bottom)*(float(n_subplots-(i_w+1))/float(n_subplots)),\
                    fig_right-fig_left, \
                    (fig_top-fig_bottom)*(1./float(n_subplots))])
#          ax.set_ylabel('Pulse width (degrees)')
          ax.set_xlim(xmin, xmax)
# Set y limits so that all plots have same scale
          ax.set_ylim(0.5*(ymin+ymax)-0.5*max_yspan, \
                           0.5*(ymin+ymax)+0.5*max_yspan)
          if(i_w < n_subplots-1):
               ax.xaxis.set_ticklabels([])
          else:
               xMajorFormatter = FormatStrFormatter('%d')
               ax.xaxis.set_major_formatter(xMajorFormatter)
          yMajorFormatter = FormatStrFormatter('%.1f')
          ax.yaxis.set_major_formatter(yMajorFormatter)
          ax.yaxis.set_major_locator(MaxNLocator(5, prune='both'))


# Now plot the widths
          ax.plot(width_data['year'], width_data['width'], 'o')
          ax.errorbar(width_data['year'], width_data['width'], \
                           width_data['werr'], fmt=None)
#          plot_widths(width_data, yunits='deg')


#  Finally, put axis labels for entire figure, and ensure that there aren't corner tick labels... maybe adjust y limits to avoid this?  There was a way of setting it, but maybe too much waster of time to find it.
     fig.text(0.5*(fig_left+fig_right), 0.06, 'Year', fontsize=16,  ha='center', va='center')
     fig.text(0.04, 0.5*(fig_top+fig_bottom), 'Pulse width (deg)', fontsize=16,  ha='center', va='center', rotation='vertical')
     #fig.text(0.1, 0.5, )

     plot_file = 'widths_multi.png'
          
     plt.savefig(plot_file)
     
main()
