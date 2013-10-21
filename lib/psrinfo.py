# Scripts for taking in various pulsar-related data products and returning
# useful information about those data.  

# can also have a prof_info routine, of a data_info routine, which looks at the raw data?

from sys import argv, exit
import os
import numpy as np

# Get info about the current fit based on residual information.  This is
# basically a re-write of David's obswgt script.
def get_resid_info(res_data, nparam=0):
    
    # If parameter data is given, get number of fitted parameters:
    if(nparam==0):
        print 'WARNING: no parameter data is given; chi^2 values will not'
        print '   reflect true number of degrees of freedom'
        n_param = 1
    else:
        n_param = nparam + 1  # +1 is to account for fit for phase
        
        
    # Set up info information for doing separate calculations for each value
    if (res_data['info'] != None):
        # Find common info numbers between the command-line 
        # requested info numbers, and those in the info file itself 
        info_common = res_data['info_val']
        n_info = len(info_common)
    else:
        n_info=1
        
    
    n_toa = len(res_data['res'])
    n_dof = n_toa - n_param    
        
    sum_weight = np.empty(n_info)
    chi2 = np.empty(n_info)
    res2 = np.empty(n_info)
    res2_weight = np.empty(n_info)
    mjd_start = np.empty(n_info)
    mjd_end = np.empty(n_info)
    centre_freq = np.empty(n_info)
    n_pts = np.empty(n_info, dtype=int)
    for i_info in np.arange(n_info):
        if(n_info > 1):
            info_condition = \
            res_data['info']==info_common[i_info]
            info_ind = np.where(info_condition)
            #res_x = res_data[xunits[i_plot]][info_ind]
            # np.where(res_data['info']==res_data['info_val'][i_info])
            res = res_data['res'][info_ind] 
            res_err = res_data['reserr'][info_ind]  
            mjd = res_data['mjd'][info_ind]   
            obsfreq = res_data['obsfreq'][info_ind]   
        else:       
            res = res_data['res']
            res_err = res_data['reserr']
            mjd = res_data['mjd']
            obsfreq = res_data['obsfreq']
        
        weight = 1.0/(res_err**2.0)
        sum_weight[i_info] = np.sum(weight)
        chi2[i_info] = np.sum((res**2.) * weight)
        res2[i_info] = np.sum(res**2.)
        res2_weight[i_info] = np.sum((res**2.) * weight)
        mjd_start[i_info] = np.amin(mjd)
        mjd_end[i_info] = np.amax(mjd)
        centre_freq[i_info] = np.mean(obsfreq)
        n_pts[i_info] = len(res)
     
    # Now calculate stats based on the above, for each info value.
    # Be sure to check what happens if n_info=1...should be fine.   
    tot_weight = np.sum(sum_weight)
    
    norm_weight = sum_weight/tot_weight
    avg_weight = norm_weight/n_pts
    red_chi2 = chi2/n_pts
    red_chi2x = red_chi2 * (n_toa/n_dof)
    res_rms = np.sqrt(res2/n_pts)
    res_rms_weight = np.sqrt(res2_weight/sum_weight)
    
    
    sum_norm_weight = np.sum(norm_weight)
    sum_avg_weight = np.sum(avg_weight*n_pts)/n_toa
    sum_red_chi2 = np.sum(chi2)/n_toa
    sum_red_chi2x = np.sum(chi2)/(n_toa - n_param)
    sum_res_rms = np.sqrt(np.sum(res2)/n_toa)
    sum_res_rms_weight = np.sqrt(np.sum(res2_weight)/np.sum(sum_weight))
    sum_n_pts = np.sum(n_pts)
    sum_mjd_start = np.amin(mjd_start)
    sum_mjd_end = np.amax(mjd_end)
    sum_centre_freq = np.mean(centre_freq)
    
    # OK, package results:
    resid_info = {'info':res_data['info_val'],
                  'npts':n_pts,
                  'ntoa':n_toa,
                  'nparam':n_param,
                  'ndof':n_dof,
                  'normwgt':norm_weight,
                  'avgwgt':avg_weight,
                  'rchi2':red_chi2,
                  'rchi2x':red_chi2x,
                  'resrms':res_rms,
                  'resrmsw':res_rms_weight,
                  'mjdstart':mjd_start,
                  'mjdend':mjd_end,
                  'cfreq':centre_freq,
                  'sum_normwgt':sum_norm_weight,
                  'sum_avgwgt':sum_avg_weight,
                  'sum_npts':sum_n_pts,
                  'sum_rchi2':sum_red_chi2,
                  'sum_rchi2x':sum_red_chi2x,
                  'sum_resrms':sum_res_rms,
                  'sum_resrmsw':sum_res_rms_weight,
                  'sum_mjdstart':sum_mjd_start,
                  'sum_mjdend':sum_mjd_end,
                  'sum_cfreq':sum_centre_freq}
    
    return resid_info
        

        

    
    
