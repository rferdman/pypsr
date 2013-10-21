#!/usr/local/bin/python

# Reads in widths from ascii file and plots them 

from sys import argv, exit
from psrread import read_dmx
from psrplot import plot_dmx
import numpy as np
import matplotlib.pyplot as plt


def main():
# Input print_resids-format data file will be first argument.
# (command-line options to come later)
     dmx_file = argv[1:]

     dmx_data =[]
     for dfile in dmx_file:
          dmx_data.append(read_dmx(dfile))
#     print res_data['mjd']

     plot_dmx(dmx_data)

     plot_file = 'dmx.png'

     plt.savefig(plot_file)

main()
