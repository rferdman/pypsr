# Scripts for dealing with XTE fits files and data

import numpy as np
import scipy as sp
from astropy.io import fits
import matplotlib.pyplot as plt
import time
from datetime import date, datetime, timedelta
import mjd
import psr_utils as pu
from psrread import read_par
from psrprof import prof_align
import glob
from utils import *
from fftfit import *
from fluxtool import *


chan_table_file = '/Users/ferdman/Work/pulsar/data_tables/xte_energy_channels.dat'
# Min/max energies for XTE in keV
xte_emin = 0.
xte_emax = 130.0

def xte_chan2energy(ph_mjd, ph_chan, ph_pcuid, return_chans=False):

    # Set up start and stop dates for channel-energy table
    start_date = [datetime(1800,1,1,0,0,0),  # just pick a time far in the past here
                  datetime(1996, 3, 21, 18, 33, 0),
                  datetime(1996, 4, 15, 23, 5, 0),
                  datetime(1999, 3, 22, 17, 37, 0),
                  datetime(2000, 5, 13, 0, 0, 0)]
    
    stop_date  = [datetime(1996, 3, 21, 18, 33, 0),
                  datetime(1996, 4, 15, 23, 5, 0),
                  datetime(1999, 3, 22, 17, 37, 0),
                  datetime(2000, 5, 13, 0, 0, 0),
                  datetime(3000, 1, 1, 0, 0, 0)]  # just pick a time far in the future here
    n_date = len(start_date)
    
    # Set chan_list to lower limit of channel range instead of a string with dashes...
    chan_list = np.loadtxt(chan_table_file, dtype=str, skiprows=13, usecols=(0,))
    #chan_list = np.append(chan_list, 256)
    for i_chan in range(len(chan_list)): 
        n_dash = chan_list[i_chan].count('-')
        if(n_dash > 0):
            chan_list[i_chan] = chan_list[i_chan][0:chan_list[i_chan].index('-')]
    chan_list = np.array(map(int, chan_list))
    
    # Now fill in 'in between' indices and map them to energy index.
    energy_ind = []
    for i_chan in range(len(chan_list)): # Not counting the last element, in which case I will append, not insert
        if(i_chan < len(chan_list)-1):
            n_diff = chan_list[i_chan+1] - chan_list[i_chan]
        else:
            n_diff = 256 - chan_list[len(chan_list) - 1]            
        for i_ind in range(n_diff):
            energy_ind.append(i_chan)
    energy_ind = np.array(energy_ind)

#    print 'n_energy = ', len(energy_ind)
#    print 'energy_ind = ', energy_ind            
    
    
    
    # add a final number on the end here for ease of testing ranges
    #print 'n_chan = ', len(chan_list)
    #print 'n_chan = ', len(chan_list)
    #print 'CHAN_LIST = ', chan_list
    # Now read in table of energies, which depend on date the datra was taken
    chan_table = np.loadtxt(chan_table_file, skiprows=13, usecols=(2,3,4,5,6,7))

    # Channel table is now indexed as chan_table[pha_range][date/pcuid]
    # Need to differentiate between PCUOD = 0 OR > 0 for start/stop_date[4]
    
    n_event = len(ph_mjd)
##    ph_chan_new = np.zeros_like(ph_chan)
    # Map ph_chan input (0-255) to a row number from table (0-128) using chan_list extracted above
##    for i_event in range(n_event):
#        for i_chan in range(len(chan_list)):
##        i_chan = 0   
##        while((i_chan < len(chan_list)-1) and (not (chan_list[i_chan] <= ph_chan[i_event] < chan_list[i_chan+1]))):
           #print 'i_chan = ', i_chan
##           i_chan += 1
##        if(i_chan >= len(chan_list)):
##            print 'Warning: PHA value at MJD {0:10.4f} is {1:d}'.format(ph_mjd, ph_chan)
##            print ' Setting to -1.'
##            ph_chan_new[i_event] = -1
##        else:
##            ph_chan_new[i_event] = i_chan 
                
    # Make array of datetime-format dates from MJDs
    ph_datetime = mjd.mjdtodate(ph_mjd)
    
    ph_chan_col = np.zeros_like(ph_datetime) - 1 # initialize to -1
    for i_date in range(n_date):
        if(i_date < 4):
            match_ind = np.where((ph_datetime > start_date[i_date]) & (ph_datetime < stop_date[i_date]))    
            ph_chan_col[match_ind] = i_date
        else:  # to accoutn for final date split between PCUID 0 vs 1/2/3/4
            match_ind = np.where((ph_datetime > start_date[i_date]) & (ph_datetime < stop_date[i_date]) & (ph_pcuid == 0))
            ph_chan_col[match_ind] = 4
            match_ind = np.where((ph_datetime > start_date[i_date]) & (ph_datetime < stop_date[i_date]) & (ph_pcuid > 0))
            ph_chan_col[match_ind] = 5
                    
    # Read off energies from table based on PHA channel keyword and date index. 
    # Energies in table are EMAX, so EMIN for a given channel is given as the energy corresponding to channel-1
    # If channel = 0, then EMIN = 0.
    ph_emin = np.array([chan_table[energy_ind[ph_chan[i_event]]-1][ph_chan_col[i_event]] if(energy_ind[ph_chan[i_event]]>0) else 0. for i_event in range(n_event)])
    ph_emin[ph_chan_col < 0] = -1.
    ph_emax = np.array([chan_table[energy_ind[ph_chan[i_event]]][ph_chan_col[i_event]] for i_event in range(n_event)])
    ph_emax[ph_chan_col < 0] = -1.
    ph_erange = np.array([(ph_emin[i_event], ph_emax[i_event]) for i_event in range(n_event)])
    
    if(return_chans==True):
        return ph_erange, ph_chan, energy_ind[ph_chan]
    else:
        return ph_erange
    
    
def read_xte_fits(fits_file, efilter=None, etol=0.0):

    hdulist = fits.open(fits_file, 'readonly')
    data_header = hdulist[1].header
    data_table = hdulist[1].data
    n_event = len(data_table)
#    ph_time = np.array([data_table[i_event]['time'] for i_event in range(n_event)])
#    ph_chan = np.array([data_table[i_event]['pha'] for i_time in range(n_event)])
#    ph_pcuid = np.array([data_table[i_event]['pcuid'] for i_time in range(n_event)])
    ph_time = data_table['time']
    ph_chan = data_table['pha']
    ph_pcuid = data_table['pcuid']
    # Close fits file 
    hdulist.close()
    ### Convert times to MJDs ###
    ph_mjd = ph_time/86400.  # first convert to units of days
    # To cover header value possibilities
    try:
      ph_mjd += data_header['mjdrefi'] + data_header['mjdreff']
    except (Exception):
      ph_mjd += data_header['mjdref'] # where only one number is given; not split into integer/fractional
            
    # Filter data by energy
    # Default if not provided (= all XTE energies) -- will filter out badly marked chans which have been denoted "-1"
    if (efilter==None):  
        efilter=[xte_emin, xte_emax]
                
    ph_erange, ph_chan, energy_ind = xte_chan2energy(ph_mjd, ph_chan, ph_pcuid, return_chans=True)
#    ph_erange = xte_chan2energy(ph_mjd, ph_chan, ph_pcuid)

    in_range = np.array([(efilter[0]-etol < ph_erange[i_event][0] < efilter[1]+etol) | (efilter[0]-etol < ph_erange[i_event][1] < efilter[1]+etol) for i_event in range(n_event)])

    ph_mjd = ph_mjd[in_range]
    # ph_energy = ph_energy[in_range]
    ph_erange = ph_erange[in_range]
    ph_time = ph_time[in_range]
    ph_chan = ph_chan[in_range]
    ph_pcuid = ph_pcuid[in_range]
    
    # Now sort by MJD
    ind_sort = np.argsort(ph_mjd)
    ph_mjd = ph_mjd[ind_sort]
    ph_erange = ph_erange[ind_sort]
    ph_time = ph_time[ind_sort]
    ph_chan = ph_chan[ind_sort]
    ph_pcuid = ph_pcuid[ind_sort]
    
    
    ph_data = {'time':ph_time, 'mjd':ph_mjd, 'chan':ph_chan, 'energy_ind':energy_ind, 'pcuid':ph_pcuid, 'erange':ph_erange,
                'datafile':fits_file}
    
    return ph_data


# Pluck out pointing information from XTE fits data files
def get_xte_pointing(xte_fits_file):
    hdulist = fits.open(xte_fits_file, 'readonly')
    primary_header = hdulist[0].header
    data_header = hdulist[1].header
    data_table = hdulist[1].data


    mean_time = 0.5 * (np.min(data_table['time']) + np.max(data_table['time']) )/86400.  # in days
    try:
        mean_mjd = mean_time + data_header['mjdrefi'] + data_header['mjdreff']
    except (Exception):
        mean_mjd = mean_time + data_header['mjdref'] # where only one number is given; not split into integer/fractional
    
    # Start dictionary, reading pointing coordinates and get mean MJD of observation
    pnt_info = {'mjd':mean_mjd}
    pnt_info['ra_pnt'] = data_header['ra_pnt']
    pnt_info['dec_pnt'] = data_header['dec_pnt']
    if('ra_nom' in data_header and 'dec_nom' in data_header):
        pnt_info['ra_bary'] = data_header['ra_nom']
        pnt_info['dec_bary'] = data_header['dec_nom']
    else:
        print xte_fits_file + ' has no ra_nom/dec_nom keywords'
        pnt_info['ra_bary'] = -1.
        pnt_info['dec_bary'] = -1.
        

    # Close fits file 
    hdulist.close()
   
    return pnt_info
   
    

# Get start and end dates covered by given filter file
# Works for data files and filter files
def get_xte_mjd_range(xte_fits_files, strip_path=False):
    outfile = []
    mjd_start = []
    mjd_end = []
    for ffile in xte_fits_files:
        hdulist = fits.open(ffile, 'readonly')
        data_header = hdulist[1].header
        data_table = hdulist[1].data
        mjd0 = data_table[0]['time']/86400.
        mjd1 = data_table[-1]['time']/86400.
        try:
          mjd0 += data_header['mjdrefi'] + data_header['mjdreff']
          mjd1 += data_header['mjdrefi'] + data_header['mjdreff']
        except (Exception):
          mjd0 += data_header['mjdref'] # where only one number is given; not split into integer/fractional
          mjd1 += data_header['mjdref'] # where only one number is given; not split into integer/fractional
    
        mjd_start.append(mjd0)
        mjd_end.append(mjd1)
        if(strip_path):
            outfile.append(ffile[ffile.rindex('/') + 1:]) # remove path from filename string
        else:
            outfile.append(ffile)
    
        hdulist.close()
    xte_file_mjd = {'file':np.array(outfile), 'mjd_start':np.array(mjd_start), 'mjd_end':np.array(mjd_end)}
    return xte_file_mjd    
    
    
# Get number of PCUs working at a given epoch, via XTE filter file
# Returns corresponding MJD and n_pcu_on a
# Allow for multiple files
def get_xte_num_pcu_on(filter_file, filter_pcu_mjd=True):
    # Open filter file (which is a FITS file)
    mjd = []
    num_pcu_on = []
    for ffile in filter_file:
        hdulist = fits.open(ffile, 'readonly')
        data_header = hdulist[1].header
        data_table = hdulist[1].data
        n_epoch = len(data_table)
        #n_pcu = 5
        
        for i_epoch in range(n_epoch):
            mjd_cur = data_table[i_epoch]['time']/86400.
            
            #mjd_cur = np.array([data_table[i_epoch]['time']/86400. for i_epoch in range(n_epoch)])
            try:
                mjd_cur += data_header['mjdrefi'] + data_header['mjdreff']
            except (Exception):
                mjd_cur += data_header['mjdref'] # where only one number is given; not split into integer/fractional
                
            mjd.append(mjd_cur)
            #pcu_on = []
            #for i_pcu in range(n_pcu):
            #    pcu_on.append([data_table[i_epoch]['pcu'+str(i_pcu)+'_on'] for i_epoch in range(n_epoch)])
            #pcu_on = np.array(pcu_on)   
    
            num_pcu_on_cur = data_table[i_epoch]['num_pcu_on']
        #num_pcu_on_cur = np.array([data_table[i_epoch]['num_pcu_on'] for i_epoch in range(n_epoch)])
        ### Convert times to MJDs ###
        #  mjd_cur = epochs/86400.  # first convert to units of days
        # To cover header value possibilities
        
        # Hardcode exclusion of PCU0 data on or after 12 May 2000 (MJD 51676), due to loss of propane layer
        # and exclusion of PCU1 data on or after 25 Dec 2006 (MJD 54094), also due to loss of propane layer
        # BE FOREWARNED that you should ensure that profile data has also accounted for these exclusions!  Otherwise they won't match!
            if(filter_pcu_mjd):
                n_exclude_pcu = 0
                if( (mjd_cur >= 51676.0) & data_table[i_epoch]['pcu0_on']):
                    n_exclude_pcu += 1
                if( (mjd_cur >= 54094.0) & data_table[i_epoch]['pcu0_on']):
                    n_exclude_pcu += 1
                num_pcu_on_cur = num_pcu_on_cur - n_exclude_pcu
                
            num_pcu_on.append(num_pcu_on_cur)
    
        hdulist.close()
        mjd = np.array(mjd)
        num_pcu_on = np.array(num_pcu_on)
    
    n_pcu_on = {'mjd':mjd, 'num':num_pcu_on}
    
    return n_pcu_on
    #return mjd, num_pcu_on
    
    
    
      
# Given date range and get_xte_mjd_range() output (list of filter files and start/end dates for each),
# find mean # of PCUs that are on from within the specified date range
def get_mean_num_pcu(filter_file_mjd, mjd_start, mjd_end):

    mjd_condition =  (filter_file_mjd['mjd_end'] >= mjd_start) & (filter_file_mjd['mjd_start'] <= mjd_end)

    ind_file = np.where(mjd_condition)[0]
    filter_files_select = filter_file_mjd['file'][ind_file]
    
    n_filter_files = len(filter_files_select)

    if(n_filter_files > 0):
        # Now grab all num_PCUs on from each file, and order them by MJD.  
        n_pcu_on = get_xte_num_pcu_on(filter_files_select)
        # Include only dates of interest:  (epoch <= start_mjd & epoch+1 > start_mjd) & (epoch > end_mjd & epoch-1 < end_mjd)
        # Get rid of anything < 1 and > 5 (will usually be set to 255 in the case of a default bad value)
        pcu_condition = (n_pcu_on['mjd'] >= mjd_start) & (n_pcu_on['mjd'] <= mjd_end) & (n_pcu_on['num'] > 0) & (n_pcu_on['num'] <= 5)
        n_pcu_on_subset = n_pcu_on['num'][pcu_condition]
        # Calculate mean num_pcu_on from these values
        mean_num_pcu = np.mean(n_pcu_on_subset)
    else:
        mean_num_pcu = -1
    
    return mean_num_pcu


def get_filter_file_info(filter_file, n_epoch_per_bin=None, n_bin=1):
    
    hdulist = fits.open(filter_file, 'readonly')
    data_header = hdulist[1].header
    data_table = hdulist[1].data
    n_epoch = len(data_table)
    # Use single_file version if don't want to bin data (dangerous, since this will be a lot of saved data!)
    if(n_bin==n_epoch):
        filter_info = get_filter_file_info_nobin(filter_file)
    else:
        if(n_epoch_per_bin!=None):
            n_bin = n_epoch/n_epoch_per_bin
        else:
            n_epoch_per_bin = n_epoch/n_bin
   
    
        # collect data in each bin and take means for plotting
        mjd = []
        mean_num_pcu = []
        mean_num_pcu_err = []
        elv_angle = []
        elv_angle_err = []
        pnt_offset = []
        pnt_offset_err = []
        time_saa = []
        time_saa_err = []
        electron0 = []
        electron0_err = []
        electron1 = []
        electron1_err = []
        electron2 = []
        electron2_err = []
        electron3 = []
        electron3_err = []
        electron4 = []
        electron4_err = []
    
        for i_bin in np.arange(n_bin):
            # Now time in seconds
            mjd.append( np.mean(data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['time']/86400.) ) 
        
            # Number of PCUs on
            n_pcu_on = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['num_pcu_on'] 
            # get rid of anything that has weird values
            pcu_condition = (n_pcu_on >= 0) & (n_pcu_on <= 5)
            mean_num_pcu.append( np.mean(n_pcu_on[pcu_condition]) )
            # take std dev to be error on quantity, and for other quantities
            mean_num_pcu_err.append( np.std(n_pcu_on[pcu_condition]) )
        
            # Elevation angle, to see if there were any Earth occultations
            elv = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['elv']
            elv_angle.append( np.mean(elv) )
            elv_angle_err.append( np.std(elv) )
        
            # Offset, to determine pointing stability
            offset = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['offset']
            pnt_offset.append( np.mean(offset) )
            pnt_offset_err.append( np.std(offset) )
        
            # time_since_saa, to determine if the Southern Atlantic Anomaly significantly affected the background levels
            t_saa = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['time_since_saa']
            time_saa.append( np.mean(t_saa) )
            time_saa_err.append( np.std(t_saa) )
        
            e0 = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['electron0']
            electron0.append( np.mean(e0) )
            electron0_err.append( np.std(e0) )

            e1 = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['electron1']
            electron1.append( np.mean(e1) )
            electron1_err.append( np.std(e1) )

            e2 = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['electron2']
            electron2.append( np.mean(e2) )
            electron2_err.append( np.std(e2) )

            e3 = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['electron3']
            electron3.append( np.mean(e3) )
            electron3_err.append( np.std(e3) )

            e4 = data_table[i_bin*n_epoch_per_bin:(i_bin+1)*n_epoch_per_bin]['electron4']
            electron4.append( np.mean(e4) )
            electron4_err.append( np.std(e4) )
    
        # Now have a set of values for this file.  Make into numpy arrays, sorted by MJD:
        mjd = np.array(mjd)
        # convert times to MJDs
        try:
            mjd += data_header['mjdrefi'] + data_header['mjdreff']
        except (Exception):
            mjd += data_header['mjdref'] # where only one number is given; not split into integer/fractional
        # Sort by MJD
        ind_sort = np.argsort(mjd)
        mjd = mjd[ind_sort]
        mean_num_pcu = np.array(mean_num_pcu)[ind_sort]
        mean_num_pcu_err = np.array(mean_num_pcu_err)[ind_sort]
        elv_angle = np.array(elv_angle)[ind_sort]
        elv_angle_err = np.array(elv_angle_err)[ind_sort]
        pnt_offset = np.array(pnt_offset)[ind_sort]
        pnt_offset_err = np.array(pnt_offset_err)[ind_sort]
        time_saa = np.array(time_saa)[ind_sort]
        time_saa_err = np.array(time_saa_err)[ind_sort]
        electron0 = np.array(electron0)[ind_sort]
        electron0_err = np.array(electron0_err)[ind_sort]
        electron1 = np.array(electron1)[ind_sort]
        electron1_err = np.array(electron1_err)[ind_sort]
        electron2 = np.array(electron2)[ind_sort]
        electron2_err = np.array(electron2_err)[ind_sort]
        electron3 = np.array(electron3)[ind_sort]
        electron3_err = np.array(electron3_err)[ind_sort]
        electron4 = np.array(electron4)[ind_sort]
        electron4_err = np.array(electron4_err)[ind_sort]
    
        
     
        # Time since various breakdowns, for each PCU
    #    t_brk0 = data_table['time_since_brk0']
    #    t_brk1 = data_table['time_since_brk1']
    #    t_brk2 = data_table['time_since_brk2']
    #    t_brk3 = data_table['time_since_brk3']
    #    t_brk4 = data_table['time_since_brk4']
    #    t_brk = data_table['time_since_brk']  # This is the union of the times since breakdown for all PCUs
    
        filter_info = {'mjd':mjd, 'elv':elv_angle, 'elv_err':elv_angle_err, 'offset':pnt_offset, 'offset_err':pnt_offset_err, 
                       't_saa':time_saa, 't_saa_err':time_saa_err, 'mean_num_pcu':mean_num_pcu, 'mean_num_pcu_err':mean_num_pcu_err,
                       'electron0':electron0, 'electron1':electron1, 'electron2':electron2, 'electron3':electron3, 'electron4':electron4,
                       'electron0_err':electron0_err, 'electron1_err':electron1_err, 'electron2_err':electron2_err, 
                       'electron3_err':electron3_err, 'electron4_err':electron4_err}                  
                       #'t_brk0':t_brk0, 't_brk1':t_brk1, 't_brk2':t_brk2, 't_brk3':t_brk3, 't_brk4':t_brk4, 't_brk':t_brk}
                   
                   
                   
    
    return filter_info


# One file no binning version

def get_filter_file_info_nobin(filter_file):
    
    hdulist = fits.open(filter_file, 'readonly')
    data_header = hdulist[1].header
    data_table = hdulist[1].data
    n_epoch_file = len(data_table)
    
    # Now time in seconds
    mjd = data_table['time']/86400.
    # convert times to MJDs
    try:
        mjd += data_header['mjdrefi'] + data_header['mjdreff']
    except (Exception):
        mjd += data_header['mjdref'] # where only one number is given; not split into integer/fractional
    
    # Sort by MJD
    ind_sort = np.argsort(mjd)
    mjd = mjd[ind_sort]
    
        
    # Number of PCUs on.  Initialize output arrays as False and then only change if value == 1
    # Doing it this way since values can be 255 as well as 0 for False...
    pcu0_on = np.zeros_like(mjd, dtype=bool)
    pcu0_on[data_table['pcu0_on'] == 1] = True
    pcu0_on = pcu0_on[ind_sort]
    pcu1_on = np.zeros_like(mjd, dtype=bool)
    pcu1_on[data_table['pcu1_on'] == 1] = True
    pcu1_on = pcu0_on[ind_sort]
    pcu2_on = np.zeros_like(mjd, dtype=bool)
    pcu2_on[data_table['pcu2_on'] == 1] = True
    pcu2_on = pcu0_on[ind_sort]
    pcu3_on = np.zeros_like(mjd, dtype=bool)
    pcu3_on[data_table['pcu3_on'] == 1] = True
    pcu3_on = pcu0_on[ind_sort]
    pcu4_on = np.zeros_like(mjd, dtype=bool)
    pcu4_on[data_table['pcu4_on'] == 1] = True
    pcu4_on = pcu0_on[ind_sort]

    num_pcu_on = data_table['num_pcu_on'][ind_sort]
    # get rid of anything that has weird values
    #        pcu_condition = (n_pcu_on >= 0) & (n_pcu_on <= 5)
    #        mean_num_pcu.append( np.mean(n_pcu_on[pcu_condition]) )
    # take std dev to be error on quantity, and for other quantities
    #        mean_num_pcu_err.append( np.std(n_pcu_on[pcu_condition]) )
        
    # Elevation angle, to see if there were any Earth occultations
    elv = data_table['elv'][ind_sort]
    #        elv_angle.append( np.mean(elv) )
    #        elv_angle_err.append( np.std(elv) )
        
    # Offset, to determine pointing stability
    offset = data_table['offset'][ind_sort]
    #        pnt_offset.append( np.mean(offset) )
    #        pnt_offset_err.append( np.std(offset) )
        
    # time_since_saa, to determine if the Southern Atlantic Anomaly significantly affected the background levels
    t_saa = data_table['time_since_saa'][ind_sort]
#        time_saa.append( np.mean(t_saa) )
#        time_saa_err.append( np.std(t_saa) )
        
    electron0 = data_table['electron0'][ind_sort]
    electron1 = data_table['electron1'][ind_sort]
    electron2 = data_table['electron2'][ind_sort]
    electron3 = data_table['electron3'][ind_sort]
    electron4 = data_table['electron4'][ind_sort]
    
    # Time since various breakdowns, for each PCU
    if('t_brk0' in data_table):
        t_brk0 = data_table['time_since_brk0'][ind_sort]
    else:
        t_brk0 = -1
    if('t_brk1' in data_table):
        t_brk1 = data_table['time_since_brk1'][ind_sort]
    else:
        t_brk1 = -1
    if('t_brk2' in data_table):
        t_brk2 = data_table['time_since_brk2'][ind_sort]
    else:
        t_brk2 = -1
    if('t_brk3' in data_table):
        t_brk3 = data_table['time_since_brk3'][ind_sort]
    else:
        t_brk3 = -1
    if('t_brk4' in data_table):
        t_brk4 = data_table['time_since_brk4'][ind_sort]
    else:
        t_brk4 = -1
    if('t_brk' in data_table):
        t_brk = data_table['time_since_brk'][ind_sort]
    else:
        t_brk = -1
    
    # Time since various breakdowns, for each PCU
#    t_brk0 = data_table['time_since_brk0']
#    t_brk1 = data_table['time_since_brk1']
#    t_brk2 = data_table['time_since_brk2']
#    t_brk3 = data_table['time_since_brk3']
#    t_brk4 = data_table['time_since_brk4']
#    t_brk = data_table['time_since_brk']  # This is the union of the times since breakdown for all PCUs

    
    
    filter_info = {'mjd':mjd, 'elv':elv, 'offset':offset, 't_saa':t_saa, 
                   'pcu0_on':pcu0_on, 'pcu1_on':pcu1_on, 'pcu2_on':pcu2_on, 'pcu3_on':pcu3_on, 'pcu4_on':pcu4_on, 'num_pcu_on':num_pcu_on,
                   'electron0':electron0, 'electron1':electron1, 'electron2':electron2, 'electron3':electron3, 'electron4':electron4,
                   't_brk0':t_brk0, 't_brk1':t_brk1, 't_brk2':t_brk2, 't_brk3':t_brk3, 't_brk4':t_brk4, 't_brk':t_brk}
    
    return filter_info

# Find mask for each event based on exclusion criteria suggested by the RXTE cookbook:
    #   1)  Include only data with number of PCUs between 1 and 5, inclusive.  This will also take care of SAA passages, for which it is given a zero value
    #   2)  Include only data with elevation angle greater than 10 degrees due to Earth occultation and bright Earth effects
    #   3)  Include only data with pointing offset < 0.02 degrees
    #   4)  Include only data with ELECTRON2 <= 0.1, to filter out electron contamination
    #   4)  Exclude PCU0 data on or after 12 May 2000 (MJD 51676), due to loss of propane layer
    #   5)  Exclude PCU1 data on or after 25 Dec 2006 (MJD 54094), due to loss of propane layer
    # For criteria 4-5, check if PCU0_on or PCU1_on is True.  If so, then it was counted in NUM_PCU_ON, and that number must be reduced by 1.
# Input will be ph_data output from read_xte_fits function and output dictionary arrays from get_xte_mjd_range(filter_files)
# to know start and end times of each filter file.
# Output will be a numpy.ma compatible mask array
def get_xte_event_filter_mask(ph_data, filter_file_mjd):
# Now figure out which filter file(s) cover this particular data file
    mjd_range = np.array([np.min(ph_data['mjd']), np.max(ph_data['mjd'])])
    mjd_condition =  (filter_file_mjd['mjd_end'] >= mjd_range[0]) & (filter_file_mjd['mjd_start'] <= mjd_range[1])
    # print mjd_range
    ind_mjd = np.where(mjd_condition)[0]
    filter_files_select = filter_file_mjd['file'][ind_mjd]
    # filter_files_select_mjd_range = np.array([filter_file_mjd['mjd_start'][ind_mjd], filter_file_mjd['mjd_end'][ind_mjd]]).flatten()
    # print filter_files_select
    # print filter_files_select_mjd_range
    # Now generalize for more than one filter file, since it can and will happen
    # Initialize so that all events have equal weight of one
    # Weights will be binary for now; zero means they are to be included in analysis, whereas ones are to be excluded.
    # This will make the output array compatible as a mask
    ph_mask = np.zeros_like(ph_data['mjd'], dtype=int)
    for i_filter in range(len(filter_files_select)):
    
        # Get filter file info for this particular filter file
        ffinfo = get_filter_file_info_nobin(filter_files_select[i_filter])
        # Initialize filter_weights to be zero for current file
        filter_weight = np.zeros_like(ffinfo['mjd'], dtype=int)
        # Assign weight of 0 or 1 to each mjd in current filter file depending on exclusion criteria
        # Remember that the following are criteria to EXCLUDE data points.
        pcu_condition = (np.isnan(ffinfo['num_pcu_on'])) | (ffinfo['num_pcu_on'] == 0) | (ffinfo['num_pcu_on'] > 5)   # Number of PCUs on
        elv_condition = (np.isnan(ffinfo['elv'])) | (ffinfo['elv'] <= 10.)  # Elevation angle - assuming this is in degrees
        offset_condition = (np.isnan(ffinfo['offset'])) | (ffinfo['offset'] >= 0.02)  # Pointing offset -- assume in degrees
        elec2_condition = (np.isnan(ffinfo['electron2'])) | (ffinfo['electron2'] > 0.1) # Electron contamination
        all_condition = (pcu_condition | elv_condition | offset_condition | elec2_condition)# | pcu0_condition | pcu1_condition))
        filter_weight[all_condition] = 1  # 1 means that it is to be masked
            
        # copy weight to all MJDs in event data that correspond to given filter file MJD
        ph_mjd_condition = (ph_data['mjd'] > np.min(ffinfo['mjd'])) & (ph_data['mjd'] <= np.max(ffinfo['mjd']))
        ind_ph_mjd = np.where( ph_mjd_condition )[0]
    #    for i_mjd in ind_mjd_filter
        # Find next highest MJD in filter data to each MJD in event data, and use that weight -- assuming arrays are already sorted, which they are:
        # Assume that relevant info is for next MJD in filter file that is larger than data MJD
        # i.e., if data MJD is between MJD1 and MJD2 in filter file, then go with MJD2 values
        ind_filter = np.searchsorted(ffinfo['mjd'], ph_data['mjd'][ind_ph_mjd], side='right')
        ph_mask[ind_ph_mjd] = filter_weight[ind_filter]
    
#        print np.min(ffinfo['mjd'])
#        print np.max(ffinfo['mjd'])
#        print ph_data['mjd']
    #    print filter_weight[ind_filter]
    #    print len(filter_weight[ind_filter])
    #    print len(filter_weight[ind_filter][filter_weight[ind_filter] ==0])
#        print ffinfo['mjd']
    #    print ph_data['mjd'][ind_ph_mjd]
    #    print ind_filter
#        print ffinfo['mjd'][ind_filter]
        
    # Finally set weights to zero if they meet the PCU0 and PCU1 exclusion criteria due to propane layer loss:
    pcu0_condition = (np.isnan(ph_data['pcuid'])) | ((ph_data['mjd'] >= 51676.0) & (ph_data['pcuid'] == 0))
    pcu1_condition = (np.isnan(ph_data['pcuid'])) | ((ph_data['mjd'] >= 54094.0) & (ph_data['pcuid'] == 1))
    ph_mask[pcu0_condition | pcu1_condition] = 1  # 1 means it is to be masked

    return ph_mask


    
# Now read in a par file and then fold data into profiles
# Read par file
def fold_event_prof(ph_mjd, par_file, nbins=16):
    params = read_par(par_file)
    # Fold data with psr_utils by first calculating phases for each MJD based on par file's F0, F1, ... only good for isolated puslar
    ph_phase = pu.calc_phs(ph_mjd, params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])
    # Cast all negative phases to be between 0 and 1
    while (len(np.where(ph_phase <= 0.)[0]) > 0 ): # Need the [0] because output is a tuple...
        below_zero_ind = np.where(ph_phase <= 0.)
        ph_phase[below_zero_ind] += 1.0
    # Now do the same to ensure anything greater than 1.0 is between (0., 1.0)
    while (len(np.where(ph_phase > 1.0)[0]) > 0 ): # Need the [0] because output is a tuple...
        above_zero_ind = np.where(ph_phase > 1.0)
        ph_phase[below_zero_ind] -= 1.0
    
        
    # Now histogram phases to create profile for a given number of bins (and thus bin size), and give user option to split data into N equal parts
    prof, bin_edges = np.histogram(ph_phase, bins=nbins, range=(0.,1.))
    bin_size = (bin_edges[1] - bin_edges[0])
    bin_val = bin_edges[0: len(bin_edges)-1] + bin_size/2.
    prof_err = np.sqrt(prof)
    # Calculate central MJD and MJD span for profile from min/max event MJDs 
    if(len(ph_mjd) > 0):
        mjd_mean = 0.5*(np.min(ph_mjd) + np.max(ph_mjd))
        mjd_span = np.max(ph_mjd) - np.min(ph_mjd)
        # spin phase at representative MJD of profile
        ref_phase = pu.calc_phs(mjd_mean, params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])
        # Cast to be between 0 and 1
        while(ref_phase <= 0.):
            ref_phase += 1.0
        while(ref_phase > 1.0):
            ref_phase -= 1.0
        ref_freq  = pu.calc_freq(mjd_mean, params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])
    else:
        mjd_mean = 0.
        ref_phase = 0.
        ref_freq = 0.
    profile = {'i':prof, 'i_err':prof_err, 'phase':bin_val, 'mjd':mjd_mean, 'mjd_span': mjd_span, 'psrname':params['psr'], 'ref_phase':ref_phase, 'ref_freq':ref_freq}
    return profile
    

def read_xte_prof(proffile):
    # Read header
    f_prof = open(proffile, 'r')
    header = f_prof.readline()
    header = header.split()
    prof_mjd = float(header[1]) + float(header[2])/86400.0
    if(float(header[3])):
        ref_freq = 1.0/float(header[3])
    else:
        ref_freq = 0.0
    # Read time span if in header
    if(len(header) > 12): 
        mjd_span = float(header[12])
    psr_name = header[10]
    ref_phase = float(header[11])
    f_prof.close()

    # read profile data
    prof_i = np.loadtxt(proffile, skiprows=1, usecols=(1,), unpack=True)
    prof_i_err = np.sqrt(prof_i)
    n_bins = len(prof_i)
    phase, bin_size = np.linspace(0., 1., num=n_bins, endpoint=False, retstep=True)
    phase += bin_size/2.
    if(len(header) > 12):
        prof = {'i':prof_i, 'i_err':prof_i_err, 'phase':phase, 'mjd':prof_mjd, 'psrname':psr_name, 'ref_phase':ref_phase, 'ref_freq':ref_freq, 'mjd_span':mjd_span}
    else:
        prof = {'i':prof_i, 'i_err':prof_i_err, 'phase':phase, 'mjd':prof_mjd, 'psrname':psr_name, 'ref_phase':ref_phase, 'ref_freq':ref_freq}
    
    return prof
    
# Write profile to file
def write_xte_prof(prof, proffile=None):
    # Calculate integer and fractional components of MJD:
    imjd = np.int_(prof['mjd'])
    fmjd = prof['mjd'] - np.float_(imjd)
    # Seconds after UT 00:00 
    smjd = fmjd*86400.0
    if(proffile==None):
        proffile = prof['psrname']+'_'+str(imjd)+'_'+str(np.int_(smjd))+'_xte_prof.asc'
    # open file for writing
    f_prof = open(proffile, 'w')
    if(prof['ref_freq'] > 0.0):
        ref_period = 1.0/prof['ref_freq'] 
    else:
        ref_period = 0.0  
    if prof.has_key('mjd_span'):  
        header = '# {0:5d}.0 {1:13.7f} {2:.12f} 1 0.0 0.0 {3:d} @ 1 {4} {5:.12f} {6:.7f}\n'.format(imjd, smjd, ref_period, len(prof['i']), prof['psrname'], prof['ref_phase'], prof['mjd_span'])
    else:
        header = '# {0:5d}.0 {1:13.7f} {2:.12f} 1 0.0 0.0 {3:d} @ 1 {4} {5:.12f}\n'.format(imjd, smjd, ref_period, len(prof['i']), prof['psrname'], prof['ref_phase'])            
#    print header
    f_prof.write(header)
    for i_bin in range(len(prof['i'])):
        f_prof.write('{0:6d} {1:.1f}\n'.format(i_bin, prof['i'][i_bin]))
    f_prof.close()
    
    
# Simple profile plot
def plot_xte_prof(prof):
    # Plot twice for clarity :)
    plt.plot(np.append(prof['phase'], prof['phase']+1.0), np.append(prof['i'], prof['i']), linestyle='steps-mid', color='blue')
    plt.errorbar(np.append(prof['phase'], prof['phase']+1.0), np.append(prof['i'], prof['i']), yerr=np.append(prof['i_err'], prof['i_err']), ecolor='blue', fmt=None, capsize=0)
    plt.xlabel('Pulse phase')
    plt.ylabel('Counts')
    

# Add profs together by date.  Input is list of profile dictionaries
# Output is list of added profile dictionaries
# Also, need parameters as read from a par file or otherwise, to recalculate reference phases and frequencies
def add_xte_prof(profile_in, params, days_add=0.5, align=False, template=None, filter_files=None):
    
    # Assume we are not adding all profs together to start
    add_all = False
    
    n_prof_in = len(profile_in)
    prof = []
    prof_mjd = []
    for prof_in in profile_in:
        prof.append(prof_in['i'])
        prof_mjd.append(prof_in['mjd'])
        
    prof = np.array(prof)
    prof_mjd = np.array(prof_mjd)
    
    # Sort profs by MJD
    sort_ind = np.argsort(prof_mjd)
    prof_mjd = prof_mjd[sort_ind]
    prof = prof[sort_ind]
    
    # If days_add is negative, then add everything in list
    if(days_add < 0):
        days_add = np.max(prof_mjd) - np.min(prof_mjd) + 10.0 # add some extra days just to be safe
        add_all = True
    
    # Initialize starting index
    i_first = 0
    # Same number of phase bins for all profiles, pre- and post-added
    prof_phase_out, bin_size = np.linspace(0., 1., num=np.shape(prof)[1], endpoint=False, retstep=True)
    prof_phase_out += bin_size/2.

#    print 'prof_mjd = ', prof_mjd

#    print 'new_mjd = '
    # While we are still before the end of the MJD array
    profile_out = []
    while(i_first < len(prof_mjd)):
        bin_inds = (prof_mjd >= prof_mjd[i_first]) & ( prof_mjd < (prof_mjd[i_first] + days_add) )
        n_prof_out = len(np.where(bin_inds)[0]) # Number of profs going into added prof
        #print prof_mjd[bin_inds]
    
        current_profs = prof[bin_inds]
        #print np.shape(current_profs), np.shape(prof), n_prof_out
        # Now, if we want to align profiles before adding, do so, rotating profiles in place:
        if(align==True):
            # In this case, choose the profile with the best 'signal', which I will take to be the largest
            # difference between peak and minimum
            if(template==None):
                i_template = 0
                max_peak = -1.
                for i_prof in range(n_prof_out):
                    peak_height = np.max(current_profs[i_prof]) - np.min(current_profs[i_prof])
                    if(peak_height > max_peak):
                        max_peak = peak_height
                        i_template = i_prof
                # Now run through and align profiles with our found template.  Skip if profile is template.
                for i_prof in range(n_prof_out):
                    if(i_prof != i_template):
                        current_profs[i_prof] = prof_align(current_profs[i_prof], current_profs[i_template])
            # If we have a template, just use it:
            else:                
                for i_prof in range(n_prof_out):
                    current_profs[i_prof] = prof_align(current_profs[i_prof], template)
    
        # Now, just add!
        prof_out = np.sum(current_profs, 0)  # I think 0 dim is time, and 1 is phase bins.  This way, preserve phase bins.\
        prof_err_out = np.sqrt(prof_out)
        # Calculate central MJD based on min and max (central) MJD in sub-array
        prof_mjd_out = 0.5*(np.min(prof_mjd[bin_inds]) + np.max(prof_mjd[bin_inds])) 
        # Calculate ref phase and period for this MJD:
        prof_ref_phase = pu.calc_phs(prof_mjd_out, params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])
            # Cast to be between 0 and 1
        while(prof_ref_phase <= 0.):
            prof_ref_phase += 1.0
        while(prof_ref_phase > 1.0):
            prof_ref_phase -= 1.0
       # prof_phase_out = ref_phase
        prof_ref_freq  = pu.calc_freq(prof_mjd_out, params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])

        # Create file name
        profile_out.append({'i':prof_out, 'i_err':prof_err_out, 'phase':prof_phase_out, 'mjd':prof_mjd_out, 'psrname':params['psr'], 'ref_phase':prof_ref_phase, 'ref_freq':prof_ref_freq})
        # Now increase i_first to point to the beginning of next set
        i_first = np.where(bin_inds)[0][-1] + 1  # = last index of previous set of indices, plus 1
        
    # If there is only one output file, then just make profile_out a single profile (i.e. not an array of profiles)
    if(add_all):
        profile_out = profile_out[0]

    return profile_out
    

def xte_toa(profile, template, params, n_trials=512):
    
    print 'Number of TOA trials per profile: ', n_trials
    
    # If given only one profile, just make it into a one-element list:
    if type(profile) is not list:
        profile = [profile]
    
    # First get amplitudes and phases from FFTing the input template profile
    t_prof_fft, t_prof_amp, t_prof_phase = cprof(template['i'])
        
    n_prof = len(profile)
    toa_int = []
    toa_frac = []
    toa_err = []
    time_lapse = []
    i_prof = 0
    for prof in profile:
        # Run a number of trials, sampling profiles from a Poisson distribution each time, and get out an 
        # mean TOA and TOA err that is equal to the std dev of the resulting distribution of TOAs
        shift = []
        stt_time = time.time()
        for i_trial in range(n_trials):
            
            # Create trial profile by sampling from a Poisson distribution based on counts in each bin
            prof_trial = np.random.poisson(prof['i'])
            # Get shift between template and profile, in bins:            
            # Now run fftfit to match the profiles and out the shifts and scales
            shift_trial,eshift,snr,esnr,b,errb,ngood = fftfit(prof_trial,t_prof_amp,t_prof_phase)
        
            n_bins = len(prof['i'])
            # ensure that shift is not negative, and less than nbins-1
            while(shift_trial <= 0.0):
                shift_trial += float(n_bins)
            while(shift_trial > float(n_bins-1)):
                shift_trial -= float(n_bins)

            shift.append(shift_trial)


        end_time = time.time()
        time_lapse.append(end_time-stt_time)
        # Now calculate a mean shift and error, we can calculate final TOA.
        shift = np.array(shift)
        shift_mean = np.mean(shift)
        shift_err = np.std(shift)
        
        # Convert to a shift in phase, and then time in days based on current pulse period:
        shift_phase = shift_mean/float(n_bins)
        shift_err_phase = shift_err/float(n_bins)
        # Get MJD closest to reference MJD that is at zero spin phase



        mjd0 = pu.calc_t0(prof['mjd'], params['pepoch'], params['f'][0], params['f'][1], params['f'][2], params['f'][3])
        # Separate out into integer and fractional times
        mjd0_int = int(mjd0)
        mjd0_frac = mjd0 - float(mjd0_int)

#        toa_f = mjd0_frac - shift_phase/prof['ref_freq']/86400.0
        toa_f = mjd0_frac + shift_phase/prof['ref_freq']/86400.0

        # Identify extra days in getting TOA calculation
        extra_days = np.floor(toa_f)
        toa_int.append(mjd0_int + int(extra_days))
        toa_frac.append(toa_f - extra_days)
        # Get final TOA error
        toa_err.append(1.0e6*shift_err_phase/prof['ref_freq']) # convert from secs to microsecs
        
#        i_prof += 1
#        print 'prof ', i_prof

 
    time_lapse = np.array(time_lapse)
    mean_time_lapse = np.mean(time_lapse) # in minutes
    print 'Mean time per profile: ', mean_time_lapse*1000., 'ms'
    print 'Total time taken: ', np.sum(time_lapse), 'sec'
    
    # Finally, package results into dictionary format:
    toa_int = np.array(toa_int)
    toa_frac = np.array(toa_frac)
    toa_err = np.array(toa_err)
    toa = {'toa_int':toa_int, 'toa_frac':toa_frac, 'toa_err':toa_err}
        
    return toa

   
# Functiopn to apply Anne's flux code to my XTE data
def get_pulsed_counts(profile, template, smooth=True):
    
    n_harmonics = 5
    
    profile_err = np.sqrt(profile)
    template_err = np.sqrt(template)
    
    # Input profile
    total_flux = np.mean(profile)
    rms_value, rms_uncertainty = rms_estimator(n_harmonics)(profile, profile_err)
    print "RMS pulsed flux:            \t%#0.7g\t+/-\t%#0.7g" % (rms_value, rms_uncertainty)
    print "RMS pulsed fraction:        \t%#0.7g\t+/-\t%#0.7g" % (rms_value/total_flux, rms_uncertainty/total_flux)
    
    # Now, template:
    t_value, t_uncertainty = rms_estimator(n_harmonics)(template, np.zeros(len(template)))
    scal = np.mean(template-amin(template))/t_value
    print "RMS to area scaling factor: \t\t%#0.7g" % scal
    print "Scaled RMS pulsed flux:     \t%#0.7g\t+/-\t%#0.7g" % (rms_value*scal, rms_uncertainty*scal)
    print "Scaled RMS pulsed fraction: \t%#0.7g\t+/-\t%#0.7g" % (rms_value*scal/total_flux, rms_uncertainty*scal/total_flux)
    
    # Set up off_pulse bin calculation:
    off_pulse_bins = None
    off_pulse_threshold = True
#    compute_area = True
    
    # default is to estimate baseline via smoothing algortithm
    if smooth:
        E = smoothed_minimum_estimator(template,
                                       off_pulse_bins,
                                       off_pulse_threshold,
                                       n_harmonics)
    else:
        E = minimum_estimator(template,
                              off_pulse_bins,
                              off_pulse_threshold,
                              None)
    (v, u) = E(profile, profile_err)
    print "Area pulsed flux:           \t%#0.7g\t+/-\t%#0.7g" % (v,u)
    print "Area pulsed fraction:       \t%#0.7g\t+/-\t%#0.7g" % (v/total_flux,u/total_flux)
    
    return rms_value, rms_uncertainty, rms_value/total_flux, rms_uncertainty/total_flux
    
##########################################################################################################################
 

    
