#!/usr/local/bin/python

from sys import argv, exit
import numpy as np
import matplotlib.pyplot as plt
from psrread import read_asc_prof
from psrplot import plot_prof

def main():
# Input asp-style ascii profile will be first argument.
# (command-line options to come later)
     prof_file = argv[1]

# First, read in profile data file, and assign each column to a separate
# numpy array     
     prof_data = read_asc_prof(prof_file)
#     print res_data['mjd']

# meaningThe following means that the 2nd argument is the 
# desired output plot file name
     if (len(argv) > 2):   
          plot_file = argv[2]
     else:
          plot_file = 'profile.png'

     plot_prof(prof_data, hgrid=True, vgrid=True)
     plt.show()
#     plt.savefig(plot_file)


main()
