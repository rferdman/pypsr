#!/usr/local/bin/python

# Reads in widths from ascii file and plots them 

from sys import argv, exit
from psrread import read_widths
from psrplot import plot_widths
import numpy as np
import matplotlib.pyplot as plt


def main():
# Input print_resids-format data file will be first argument.
# (command-line options to come later)
     width_file = argv[1:]

     width_data =[]
     for wfile in width_file:
          width_data.append(read_widths(wfile))
#     print res_data['mjd']

     plot_widths(width_data, yunits='deg')

     plot_file = 'widths.png'

     plt.savefig(plot_file)

main()
