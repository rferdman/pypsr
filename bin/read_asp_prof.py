#!/usr/local/bin/python

# Script to read in asp fits-format profile data and plot a simple 
# Stokes-I pulse profile

from sys import argv, exit
import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt

def read_asp_info(asphdu):

# make primary header its own thing for later use
    asphead = asphdu[0].header
# make all keywords lowercase
    asp_info = dict((k.lower(), v) for k,v in dict(asphead).items())
  
# number of dumps
    n_scan = (len(asphdu)-3)/2  # asp/gasp (not nancay)
    asp_info['n_scan'] = n_scan

# Get channels
    chanhdu = 1
    
# Read this table into an array
# the 0 element is because it is stored as a 2D array (in this case with only one column)
    chandata = asphdu[chanhdu].data[0] 

# Number of channels is the first entry
    n_chan = chandata[0]
    asp_info['n_chan'] = n_chan

# Now read the next n_chan array elements into an array of channel frequencies
    chan_freq = chandata[1:n_chan+1]    
    asp_info['chan_freq'] = chan_freq

# Calculate centre frequency and store it
    centre_freq = (chan_freq[0] + chan_freq[n_chan-1])/2.
    asp_info['centre_freq'] = centre_freq

# Now get DM from second table after primary header
    dmhdu = 2
    dmdata = asphdu[dmhdu].data[0]
    dm = dmdata[0]
    asp_info['dm'] = dm

    return asp_info


def read_asp_scan(i_dump, asphdu):
    dumpref_head = asphdu['dumpref'+str(i_dump)].header
    dumpref_data = asphdu['dumpref'+str(i_dump)].data
    stokes_head  = asphdu['stokes'+str(i_dump)].header
    stokes_data  = asphdu['stokes'+str(i_dump)].data
    
    n_chan = stokes_head['tfields']/4

# Make list of dictionaries, each having Stokes profiles for each channel 
# Remember column 
    prof = [{'i':stokes_data.field(4*i_chan + 0), \
                 'q':stokes_data.field(4*i_chan + 1),  \
                 'u':stokes_data.field(4*i_chan + 2), \
                 'v':stokes_data.field(4*i_chan + 3)} for i_chan in range(n_chan)]
    mid_secs = dumpref_head['midsecs']
    ref = {'mid_secs':mid_secs, 'phase':dumpref_data.field(0), 'period':dumpref_data.field(1)}

    return ref, prof

# Read profile from ascii text file (in pspmtoa/asp format)
def read_asc_prof(prof_file):
    # Open ascii file
     try:
          f_prof = open(prof_file,'r')
     except IOError as (errno, strerror):
          if (errno == 2):  # file not found
               print "IOError ({0}): File".format(errno), prof_file, "not found"
          else:
               print "IOError ({0}): {1}".format(errno, strerror)
          return

     print "File", prof_file, "open."

# Read columns in file into numpy arrays
     prof_array = \
         np.loadtxt(f_prof, dtype='double', comments='#', usecols=(1,2,3,4))

# Get number of columns 
     n_cols = prof_array.shape[1]
# Assuming this is the order of Stokes parameters (and that, e.g., if there is only 
# one column, it is Stokes I
     stokes_keys_temp = ['i', 'q', 'u', 'v']
# Make sure that we are mapping the correct number of columns' worth of profile data
     stokes_keys = stokes_keys_temp[0:n_cols]
# Create a dictionary that maps each Stokes paramater to teh appropriate array
     prof_data = dict(zip(stokes_keys, prof_array.transpose()))
     
# Create dictionary entries for linear polarisation and position angle (and error)
     prof_data['lin'] = np.sqrt(prof_data['q']*prof_data['q'] + \
                                    prof_data['u']*prof_data['u'])
     prof_data['phi'] = np.arctan2(prof_data['u'], prof_data['q'])
     prof_data['phi_err'] = np.sqrt(( prof_data['q']*prof_data['q'] /   \
                                       (np.std(prof_data['u'])*np.std(prof_data['u'])) + \
                                      prof_data['u']*prof_data['u'] /   \
                                       (np.std(prof_data['q'])*np.std(prof_data['q'])))/ \
                                      (prof_data['lin']*prof_data['lin']*  \
                                           prof_data['lin']*prof_data['lin']))

# Now rewind file position and read header line if it exists (i.e. if line begins 
# with '#')
     f_prof.seek(0)
     header = f_prof.readline().split()
     if(header[0] == '#'):  # Header is in expected format
         del header[0]    # get rid of leading '#'
         del header[7:9]  # get rid of junk '1's
# Now convert appropriate values to int, float, etc., and leave rest (nobs, src_name) 
# as strings
         header[0] = float(header[0]) # imjd
         header[1] = float(header[1]) # smjd
         header[2] = float(header[2]) # p0
         header[4] = float(header[4]) # freq
         header[5] = float(header[5]) # dm
         header[6] =   int(header[6]) # nbins
         header[8] = float(header[8]) # phase
         header_keys = ['imjd', 'smjd', 'p0', 'n_obs', 'freq', 'dm', 'n_bins', 'src_name', 'phase']
         prof_info = dict(zip(header_keys[0:len(header)-1], header))
     else:
         print 'No header line found in file '+prof_file
         prof_info = None

# That's it! Close and return profile data:
     f_prof.close()

     return prof_data, prof_info


def plot_asc_prof(prof_file, plot_file='prof_out.eps', pol=False):
# First, read in ascii profile data file:
    prof_data, prof_info = read_asc_prof(prof_file)
    phase = np.arange(float(prof_info['n_bins']))/float(prof_info['n_bins'])

    if(pol):
        plt.plot(phase, prof_data['lin'], 'red')
        plt.plot(phase, prof_data['v'], 'blue')
    plt.plot(phase, prof_data['i'], 'black')
    min_i = min(prof_data['i'])
    max_i = max(prof_data['i'])
    max_diff = max_i - min_i
    plt.axis([min(phase), max(phase)+0.001, min_i - 0.01*max_diff, max_i + 0.01*max_diff])
    plt.grid(True)
    plt.xlabel('Pulse phase')
    plt.ylabel('Flux density')
    
    plt.savefig(plot_file)


#def main():

#    prof_file  = argv[1]
#    plot_file  = 'profile_plot.eps'

#    plot_asc_prof(prof_file, pol=True)


#main()



