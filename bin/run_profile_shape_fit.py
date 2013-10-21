#!/usr/local/bin/python
from sys import exit
import numpy as np
import glob
import argparse
import psrcalc as pc
from read_width import read_width_cos_phi
from psrread import read_widths
from psrplot import plot_widths
import matplotlib.cm as cm
from mjd import *
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
from pdfplot import *
from astropy import units as u


from profile_shape_fit import *

# from grid_setup_1913 import *
# from grid_setup_1756 import *

# gloabel constants -- precession period and orbital inclination:
#from width_cos_phi import prec_period, incl, mjd, y, y_err, mjd_mid

n_incl = 1 # Set to 2 if you want to combine probs from incl AND 180-incl

use_parkes_data = False
# parkes_data_file = parkes_data_file_all

# data_file = data_file_all

# Set up command line arguments
def get_opt(progname):
    parser = argparse.ArgumentParser( 
        prog=progname,
        description='Given input width files and pulsar parameters, this will fit the Rafikov and Lai model for the geometry of the system')
        
    parser.add_argument('wfiles',
                        nargs='+',
                        help='Input width data files.  Each file represents a given fractional pulse height, so that all data corresponding to a given pulse height should be pre-concatenated into a single file.')
    parser.add_argument('-f', '--filebase', dest='data_filebase',
                        default=None,
                        help='Base name for saved data files.  Will be the root of the probability, parameter, and weights file names.')
    parser.add_argument('-o', '--outdir', dest='outdir',
                        default='.',
                        help='Directory in which to place resulting files.')
    parser.add_argument('-s', '--save', dest='save_results',
                        action='store_true',
                        default=False,
                        help='Save current run\'s results in numpy save file format. Overwrites existing data files. File names taken from -f or --filebase option.  If not given, will use default file name.')
    parser.add_argument('-l', '--load', dest='load_results',
                        action='store_true',
                        default=False,
                        help='Load a previous run\'s result data saved in numpy save file format for the current run. File names taken from -f or --filebase option.  If not given, will use default file name.')
    parser.add_argument('--nalpha', dest='alpha_sampling',
                        nargs='+',
                        type=int,
                        default=64,
                        help='Set of number of samples to use in alpha, corresponding to number of separate alpha ranges.')
    parser.add_argument('--ndelta', dest='delta_sampling',
                        nargs='+',
                        type=int,
                        default=64,
                        help='Set of number of samples to use in delta, corresponding to number of separate delta ranges.')
    parser.add_argument('--nt1', dest='T1_sampling',
                        nargs='+',
                        type=int,
                        default=64,
                        help='Set of number of samples to use in T1, corresponding to number of separate T1 ranges.')
    # For alpha and delta, we will convert to radians below.
    parser.add_argument('--arange', dest='alpha_ranges',
                        nargs='+',
                        type=float,
                        default=[0.0, 180.0],
                        help='Define range boundaries for alpha samples in grid fit, in degrees.  Apart from the first and last argument, each argument will be high boundary to one range and low boundary to the following range.')   
    parser.add_argument('--drange', dest='delta_ranges',
                        nargs='+',
                        type=float,
                        default=[0.0, 90.0],
                        help='Define range boundaries for delta samples in grid fit, in degrees.  Apart from the first and last argument, each argument will be high boundary to one range and low boundary to the following range.')   
    parser.add_argument('--trange', dest='T1_ranges',
                        nargs='+',
                        type=float,
                        help='Define range boundaries for delta samples in grid fit, in years.  Apart from the first and last argument, each argument will be high boundary to one range and low boundary to the following range.')   
    parser.add_argument('--rholim', dest='rho_lim',
                        nargs=2,
                        type=float,
                        help='Set minimum and maximum values of rho to fit for, in degrees')
    parser.add_argument('--rhoguess', dest='rho_guess',
                        type=float,
                        help='Set of initial guesses for rho values. Must have same number of arguments as number of width files (or number of fractional pulse heights being fit simultaneously).')
    # Or, for T1, ask user for central year; will then calculate ranges based on 
    # precession period and overwrite T1_ranges argument value.  This means that
    # of this option is selected, we will only have one T1 range: 
    # T1_centre +/- prec_period/2.
    parser.add_argument('--tcentre', '--tcenter', dest='T1_centre', 
                        type=float,
                        help='Central time in years for T1 range. Overrides --trange option.')   
    # Now some pulsar parameters:
    parser.add_argument('-p', '--psr', dest='psr_name',
                        required=True,
                        help='Optionally input the pulsar name.  Will be used for file name construction, among other things')
    # Inclination in degrees.  Must know this value for now.  In the future, may
    # include possibility of Monte-Carlo-type functionalituy with this value if 
    # it is not constrained.
    parser.add_argument('-i', '--inclination', dest='incl',
                        type=float,
                        help='Orbital inclination in degrees.  For now, must know this value.')
    # For precession period, user can either input the value themself, or can 
    # input parameters to enable its calculation.
    parser.add_argument('--pprec', dest='prec_period', # will convert to days
                        type=float,
                        help='Precession period in years.') 
    parser.add_argument('--ma',
                        type=float,
                        help='Pulsar mass in solar masses.')
    parser.add_argument('--mb',
                        type=float,
                        help='Companion mass in solar masses.')
    parser.add_argument('--ecc',
                        type=float,
                        help='Binary eccentricity.')
    parser.add_argument('--pb',
                        type=float,
                        help='Binary orbital period in days.')
                        
                        
    args=parser.parse_args()
    
    # Now deal with argument input.
    
    # Start with precession period. Convert to days.
    if(args.prec_period is not None):
        args.prec_period = u.year.to(u.day, args.prec_period)
    else: # Make sure ALL pulsar parameters are given
        if((args.ma is not None) & (args.mb is not None) & \
          (args.ecc is not None) & (args.pb is not None) ):
            args.prec_period = pc.get_prec_period(ma=args.ma, mb=args.mb, 
                        ecc=args.ecc, pb=args.pb, units='days')
        else:
            print 'ERROR: If --pprec is not given, then all of --ma, --mb, --ecc, and --pb must be given.'
            exit()        
    
    # Ensure that inclination is provided. Convert to radians.
    if(args.incl is not None):
        args.incl = np.radians(args.incl)
    else:
        print 'ERROR: Must provided orbital inclination via -i or --inclination.'
        exit()
                     
    # Set data file name.  Include pulsar name is given   
    if(args.data_filebase is None):
        args.data_filebase='grid_data_rl'
        if(args.psr_name):
            args.data_filebase += '_'+args.psr_name
    # Add output directory to name
    args.data_filebase = args.outdir+'/'+args.data_filebase

    # Make sure user doesn't save and load data at the same time
    if( (args.save_results == True) & (args.load_results == True) ):
        print 'ERROR: Must choose either to save OR load result data, not both'
        exit()
        
    # Now set up alpha, delta, and T1 ranges
    if(len(args.alpha_ranges) >= 2):
        temp_range = []
        for i_range in np.arange(len(args.alpha_ranges)-1):
            temp_range.append([args.alpha_ranges[i_range], 
                            args.alpha_ranges[i_range+1]])
        args.alpha_ranges = np.radians(np.array(temp_range))
    else: # only 1 argemnt -- not allowed
        print 'ERROR: --arange accepts 2 or more arguments.'
        exit()
 
    if(len(args.delta_ranges) >= 2):
        temp_range = []
        for i_range in np.arange(len(args.delta_ranges)-1):
            temp_range.append([args.delta_ranges[i_range], 
                            args.delta_ranges[i_range+1]])
        args.delta_ranges = np.radians(np.array(temp_range))
    else: # only 1 argemnt -- not allowed
        print 'ERROR: --drange accepts 2 or more arguments.'
        exit()
    
    # For T1, first check if both --trange and --tcentre was used.  If so, 
    # show warning that we will use --tcentre
    if( (args.T1_centre is not None) & (args.T1_ranges is not None) ):
        print 'WARNING: Both --tcentre and --trange were used.'
        print 'Will use --tcentre for determining T1 ranges.'
        
    if(args.T1_centre is not None):
        args.T1_ranges = np.array([[yeartomjd(args.T1_centre) - args.prec_period/2.,
                            yeartomjd(args.T1_centre) + args.prec_period/2.]])
    else:
        if(len(args.T1_ranges) >= 2):
            temp_range = []
            for i_range in np.arange(len(args.T1_ranges)-1):
                temp_range.append([args.T1_ranges[i_range], 
                                args.T1_ranges[i_range+1]])
            args.T1_ranges = np.array(temp_range)
        else: # only 1 argemnt -- not allowed
            print 'ERROR: --drange accepts 2 or more arguments.'
            exit()
 
    limiting_year=1902.1
    if(np.amin(args.T1_ranges) < yeartomjd(limiting_year)):
        date_diff = yeartomjd(limiting_year) - np.amin(args.T1_ranges)
        args.T1_ranges += date_diff
        print '\nWARNING: chosen minimum T1 range value is less than the year '
        print '1902, which (for now) causes problems.  Will increase centre '
        print 'T1 year by the difference between the two. This gives us new '
        print 'range(s) for T1 sampling.\n'
 
 
        
    # Finally, check sampling arrays to ensure they have the same length as the
    # range arrays
    if(len(args.alpha_sampling) != len(args.alpha_ranges)):
        print 'ERROR: number of argument(s) to --nalpha and --arange must be equal.'
        exit()
    if(len(args.alpha_sampling) != len(args.alpha_ranges)):
        print 'ERROR: number of argument(s) to --ndelta and --drange must be equal.'
        exit()
    if(len(args.alpha_sampling) != len(args.alpha_ranges)):
        print 'ERROR: number of argument(s) to --nt1 and --trange must be equal.'
        exit()
    
    # Set default value for rho_lim:
    if(args.rho_lim is None):
        args.rho_lim = np.array([0.0, np.pi])
    else:
        args.rho_lim = np.radians(np.array(args.rho_lim))

    if(args.rho_guess is None):
        # Just choose a default rho value guess for all heights
        args.rho_guess = np.radians(np.array(len(args.wfiles)*[10.0]))
    else:
        if(len(args.rho_guess) == len(args.wfiles)):
            # convert to radians
            args.rho_guess=np.radians(np.radians(args.rho_guess))
        else:
            print 'ERROR: Number of arguments to --rhoguess MUST equal number of input data files.'
            exit()
    
    return args
        


def main():
    
    progname = 'run_profile_shape_fit.py'
    args = get_opt(progname)
    
    input_files = []
    for infile in args.wfiles:
        # For some reason this works and ".appen()" doesn't:
        input_files[len(input_files):] = glob.glob(infile)
    
    n_alpha = np.sum(args.alpha_sampling)
    n_delta = np.sum(args.delta_sampling)
    n_T1    = np.sum(args.T1_sampling)
    
    print 'Pulsar parameters:'
    print '   Pulsar:  '+args.psr_name
    print '   Precession period:  '+str(u.day.to(u.year, args.prec_period))+' years'
    print '   Inclination:  '+str(np.degrees(args.incl)) + ' deg'
    print 'Fitting parameters:'
    print '   Alpha: '+str(n_alpha)+' points spanning the ranges '+ \
        str(np.degrees(args.alpha_ranges)) + ' with corresponding sampling '+ \
        str(args.alpha_sampling)
    print '   Delta: '+str(n_delta)+' points spanning the ranges '+ \
        str(np.degrees(args.delta_ranges)) + ' with corresponding sampling '+ \
        str(args.delta_sampling)
    print '   T1:    '+str(n_T1)+' points spanning the ranges '+ \
        str(args.T1_ranges) + ' with corresponding sampling '+ \
        str(args.T1_sampling) + '\n'
    
    # We will have a list of dictionaries, each representing a given input file
    wdata=[]
    # Run through data files, and creates list of dictionaries of width data
    for i_file in np.arange(len(args.wfiles)):
        print 'Reading in width data file '+args.wfiles[i_file]
        # Read in data from each file and append cos_phi values
        wdata.append(read_widths(args.wfiles[i_file], 
                    units_in='phase', units_out='cos_phi'))
        
        
 
# Start doing the fit to the Rafikov and Lai model
   ####### mjd, y, y_err = read_width_cos_phi(data_file)

# Finally, plot cos(phi):
    #print "Y: ", y
    #print "YERR: ", y_err
##    wdata = {'mjd':mjd, 'width':y, 'werr':y_err}
##    plot_widths(wdata)
##    plt.savefig('cos_phi_rl.png')

    #######mjd_mid = np.mean([np.min(mjd), np.max(mjd)])
#    init_guess = [1.55, 0.14, .38, mjd_mid+1056.]

# If we want to include the Parkes data, we must read in the corresponding
# width file, then subtract the difference, adding the appropriate errors
# on the means to the Parkes widths in quadrature: 
    #####if(use_parkes_data == True):
        #####mjd_parkes, y_parkes, y_err_parkes = \
            ######read_width_cos_phi(parkes_data_file, rfile=data_file)
        
#        y_parkes = y_parkes - (np.average(y_parkes, weights=1./(y_err_parkes**2)) - np.average(y, weights=1./(y_err**2)))
#        y_err_parkes = np.sqrt(y_err_parkes**2 + \
#                                (np.std(y_parkes))**2 + (np.std(y))**2)
        
        ##### mjd = np.append(mjd_parkes, mjd)
        ##### y = np.append(y_parkes, y)
        ##### y_err = np.append(y_err_parkes, y_err)

    # Set up data file names -- leaving off '.npy' because when saving it is 
    # not included, whereas when loading, it IS included...
    grid_param_file = args.data_filebase + '_param'
    grid_prob_file = args.data_filebase + '_prob'

##### Set up data for fit #######
    mjd=[]
    y = []
    y_err = []
    for w in wdata:
        mjd.append(w['mjd'])
        y.append(w['width'])
        y_err.append(w['werr'])
    mjd = np.array(mjd)
    y = np.array(y)
    y_err = np.array(y_err)

# If not loading data in, then do grid fit:
    if (args.load_results != True):

# Make input variables into on nice easy-to-read dictionary:
        input_data = {'mjd':mjd, 'y':y, 'y_err':y_err, 
                      'incl':args.incl, 'prec_period':args.prec_period}
        input_params = {'rho_lim':args.rho_lim, 
                        'rho_guess':args.rho_guess,
                        'alpha_ranges':args.alpha_ranges, 
                        'delta_ranges':args.delta_ranges, 
                        'T1_ranges':args.T1_ranges, 
                        'alpha_sampling':args.alpha_sampling, 
                        'delta_sampling':args.delta_sampling, 
                        'T1_sampling':args.T1_sampling}

##### Do fit: #####   
        p_out = profile_shape_fit(input_params, input_data)

# Expand output dictionary:
        norm_like = p_out['norm_like']
        norm_vol = p_out['norm_vol']
#        alpha_weights = p_out['alpha_weights']
#        delta_weights = p_out['delta_weights']
#        T1_weights = p_out['T1_weights']
        alpha = p_out['alpha']
        delta = p_out['delta']
        T1 = p_out['T1']
        rho = p_out['rho']
        rho_prob = p_out['rho_prob']

# If saving results to file for later use, save pdf array separately from 
# parameter values
        if(args.save_results == True):
            prob_array = np.array([norm_like, norm_vol])
            np.save(grid_prob_file, prob_array)
#           weight_array = np.array([alpha_weights, delta_weights, T1_weights])
#            np.save(grid_weight_file_rl, weight_array)
            param_array = np.array([alpha, delta, T1, rho, rho_prob])
            np.save(grid_param_file, param_array)
            print "norm_like shape  = ", norm_like.shape
            print "norm_vol shape  = ", norm_vol.shape
            print "param_array shape = ", param_array.shape

# Otherwise load data:
    else:
        prob_array = np.load(grid_prob_file_rl+'.npy')
        norm_like = prob_array[0]
        norm_vol = prob_array[1]
        print 'Loading data from {0:s}.'.format(grid_prob_file_rl+'.npy')
#        weight_array = np.load(grid_weight_file_rl+'.npy')
#        print 'Loading data from {0:s}.'.format(grid_weight_file_rl+'.npy')
#        alpha_weights = weight_array[0]
#        delta_weights = weight_array[1]
#        T1_weights = weight_array[2]
        param_array = np.load(grid_param_file_rl+'.npy')
        print 'Loading data from {0:s}.'.format(grid_param_file_rl+'.npy')
        alpha = param_array[0]
        delta = param_array[1]
        T1 = param_array[2]
        rho = param_array[3]
        rho_prob = param_array[4]

# Also need to generally redefine n_alpha in, etc. since file may be different 
# than input parameter file
    n_alpha = len(alpha)
    n_delta = len(delta)
    n_T1 = len(T1)
    n_dims = (n_alpha, n_delta, n_T1)

    print "n_dims = ", n_dims
    # print "alpha = ", np.degrees(alpha)
    # print "delta = ", np.degrees(delta)
    # print "T1 = ", T1

# Collapse 3-d chi2 grid into 3 2-d grids (alpha/delta, alpha/T1, delta/T1)

# The likelihood returned by grid_fit is normalized.  So, in summing 
# along an axis, :
#	 np.sum(norm_like * weight[i_axis], axis = i_axis) = 1

# Sum along each direction in turn, giving a 2-d array for each pair of parameters:

# Convert T1 to years
    T1_year = np.zeros_like(T1)
    for i_T1 in np.arange(len(T1)):
        T1_year[i_T1] = mjdtoyear(T1[i_T1])

# Get delta/T1 array:
    prob_delta_T1 = np.sum(norm_vol, axis=0)
    plot_contour_pdf(np.degrees(delta), T1_year, \
                         np.transpose(prob_delta_T1), \
###                         weights=np.transpose(delta_T1_weights), \
                         xlabel='$\\delta$ (degrees)', ylabel='$T_1$ (year)')
    plt.savefig('contour_delta_T1_rl.png')

# Get alpha/T1 array:
    prob_alpha_T1 = np.sum(norm_vol, axis=1)
    plot_contour_pdf(T1_year, np.degrees(alpha), \
                         prob_alpha_T1, \
###                          weights=alpha_T1_weights, \
                         xlabel='$T_1$ (year)', ylabel='$\\alpha$ (degrees)')
    plt.savefig('contour_alpha_T1_rl.png')

# Get alpha/delta array:
    prob_alpha_delta = np.sum(norm_vol, axis=2)
    plot_contour_pdf(np.degrees(delta), np.degrees(alpha), \
                         prob_alpha_delta, \
###                         weights=alpha_delta_weights, \
                         xlabel='$\\delta$ (degrees)', ylabel='$\\alpha$ (degrees)')
    plt.savefig('contour_alpha_delta_rl.png')


# (3) Collapse 3-d grid into 3 1-d histograms to get posteriors for each of
# alpha, delta, and T1

# Finally, get individual PDFs by summing along one more axis in the above arrays
# (NB there will thus be two ways to get each PDF, but it doesn't matter which is 
# chosen for each)
    prob_intervals = np.array([0.683, 0.954, 0.9973])

# alpha:
    alpha_vol = np.sum(prob_alpha_delta, axis=1)
    # Sum normalized likelihood first over the T1 then delta axes to 
    # get alpha pdf:
    alpha_pdf = norm_like.sum(axis=2).sum(axis=1)
    alpha_med, alpha_prob_min, alpha_prob_max = \
        get_pdf_prob(np.degrees(alpha), alpha_vol, prob_intervals) ### , \
###                          weights=alpha_weights)
    plot_pdf(np.degrees(alpha), alpha_pdf, \
###                  weights=alpha_weights, \
                 xlabel='$\\alpha$ (degrees)', ylabel='Probability density',\
                 prob_lines=np.append(alpha_prob_min, alpha_prob_max),\
                 prob_linestyle=['dashed','dashdot','dotted', \
                                     'dashed','dashdot','dotted'])
    plt.savefig('pdf_alpha_rl.png')
    
    print " "
    print "ALPHA = ", alpha_med
    print "  68%: ", alpha_prob_min[0], "  ", alpha_prob_max[0]
    print "  95%: ", alpha_prob_min[1], "  ", alpha_prob_max[1]
    print "  99%: ", alpha_prob_min[2], "  ", alpha_prob_max[2]
    print " "

# delta:
###    pdf_delta = np.sum(prob_delta_T1*delta_T1_weights, axis=1)
    delta_vol = np.sum(prob_delta_T1, axis=1)
    # Sum normalized likelihood over T1, then alpha axis to get delta pdf:
    delta_pdf = norm_like.sum(axis=2).sum(axis=0)
    if (np.degrees(np.amax(args.delta_ranges)) > 100.):
        delta_upper = \
        get_pdf_prob(np.degrees(delta[0:n_delta/2]), \
                         delta_vol[0:n_delta/2], \
                         prob_intervals, \
###                         weights=delta_weights[0:n_delta/2], \
                         upper=True)
        delta_lower = \
            get_pdf_prob(np.degrees(delta[n_delta/2:n_delta]), \
                             delta_vol[n_delta/2:n_delta], \
                             prob_intervals, \
###                             weights=delta_weights[n_delta/2:n_delta], \
                             lower=True)

        plot_pdf(np.degrees(delta), delta_pdf, \
###                     weights=delta_weights, \
                     xlabel='$\\delta$ (degrees)', \
                     ylabel='Probability density', \
                     prob_lines=np.append(delta_upper, delta_lower),\
                     prob_linestyle=['dashed','dashdot','dotted', \
                                         'dashed','dashdot','dotted'])
        print " "
        print "DELTA: "
        print "  68%: < ", delta_upper[0], " / > ", delta_lower[0]
        print "  95%: < ", delta_upper[1], " / > ", delta_lower[1]
        print "  99%: < ", delta_upper[2], " / > ", delta_lower[2]
        print " "
    else:
        delta_upper = \
            get_pdf_prob(np.degrees(delta[0:n_delta]), \
                             delta_vol[0:n_delta], \
                             prob_intervals, \
###                             weights=delta_weights[0:n_delta], \
                             upper=True)

        plot_pdf(np.degrees(delta), delta_pdf, \
###                     weights=delta_weights, \
                     xlabel='$\\delta$ (degrees)', \
                     ylabel='Probability density', \
                     prob_lines=delta_upper,\
                     prob_linestyle=['dashed','dashdot','dotted'])
        print " "
        print "DELTA: "
        print "  68%: < ", delta_upper[0]
        print "  95%: < ", delta_upper[1]
        print "  99%: < ", delta_upper[2]
        print " "

    plt.savefig('pdf_delta_rl.png')
    


# T1:
###    pdf_T1 = np.sum(prob_alpha_T1*alpha_T1_weights, axis=0)
    T1_vol = np.sum(prob_alpha_T1, axis=0)
    # Sum normalized likelihood over delta, then alpha to get T1 pdf:
    T1_pdf = norm_like.sum(axis=1).sum(axis=0)
# Now figure out peaks in T1 distribution, and overplot our data span using 
# boundary lines.  Do this by taking first peak, then find second by removing 
# peak MJD +/- prec_period/4 (so that prec_period/2. of T1 remains)..
    T1_peak = np.zeros(2)
    peak_ind = T1_pdf.argmax()
    T1_peak[0] = T1_year[peak_ind]
# peak is too close to first index:
    if(peak_ind - n_T1/4 < 0):
        extra_ind = n_T1/4 - peak_ind
        temp_T1 = T1_year[peak_ind + n_T1/4 : n_T1-extra_ind]
        temp_pdf = T1_pdf[peak_ind + n_T1/4 : n_T1-extra_ind]
# peak is too close to last index:
    elif(peak_ind + n_T1/4 > n_T1):
        extra_ind = (peak_ind + n_T1/4) - n_T1
        temp_T1 = T1_year[extra_ind : peak_ind - n_T1/4]
        temp_pdf = T1_pdf[extra_ind : peak_ind - n_T1/4]
# peak falls within +/- prec_period of edges:
    else:
        temp_T1 = np.append(T1_year[peak_ind + n_T1/4 :  n_T1], \
                                     T1_year[0:peak_ind - n_T1/4])
        temp_pdf = np.append(T1_pdf[peak_ind + n_T1/4 :  n_T1], \
                                     T1_pdf[0:peak_ind - n_T1/4])

    peak_ind = temp_pdf.argmax()
    T1_peak[1] = temp_T1[peak_ind]
    mjd_mid = np.mean([np.amin(mjd), np.amax(mjd)])
    print " "
    print "T1: "
    print "  peak 1: ", T1_peak[0], " = MJD ", yeartomjd(T1_peak[0])
    print "  peak 2: ", T1_peak[1], " = MJD ", yeartomjd(T1_peak[1])
    print "  data range =  ", np.amin(mjd), np.amax(mjd)
    print "  mjd mid = ", mjd_mid
    
    
    plot_pdf(T1_year, T1_pdf, \
###                 weights=T1_weights, \
                 xlabel='$T_1$ (year)',ylabel='Probability density',\
                 prob_lines = [mjdtoyear(np.amin(mjd)), \
                                   mjdtoyear(np.amax(mjd))],\
                 prob_linestyle = 'dashed')
    plt.savefig('pdf_T1_rl.png')
    

###### Now deal with rho values #####
    n_rho = len(args.rho_guess)
    for i_rho in np.arange(n_rho):
        
        # Make a histogram of rho values to see where majority of probability 
        # lies, using the fit reduced chi squared as a weight
        pdf_rho, bin_edges = np.histogram(rho[i_rho], 256, 
                                        range=(args.rho_lim[0], args.rho_lim[1]), 
                                        density=True, weights=rho_prob)
        # We can define the bin centres as follows since our call to np/histogram gives 
        # back evenly spaced bins
        bin_size = bin_edges[1] - bin_edges[0]
        rho_val = bin_edges[0:len(bin_edges)-1] + 0.5*bin_size
        # Get PDF intervals and values:
        # pdf_rho = rho_hist/np.sum(rho_hist)
        rho_med, rho_prob_min, rho_prob_max = \
                                get_pdf_prob(np.degrees(rho_val), pdf_rho, 
                                                     prob_intervals, norm=True)

        print ' '
        print 'RHO '+str(i_rho)+' = ', rho_med
        print '  68%: ', rho_prob_min[0], '  ', rho_prob_max[0]
        print '  95%: ', rho_prob_min[1], '  ', rho_prob_max[1]
        print '  99%: ', rho_prob_min[2], '  ', rho_prob_max[2]
        print ' '

        # Plot rho PDFs:
        # Set colour
        if(n_rho==1):
            clr='black'
        else:
            clr = cm.jet(float(i_rho)/float(n_rho)-1)
        if(i_rho==0):
            min_pdf_rho = np.amin(pdf_rho)
            max_pdf_rho = np.amax(pdf_rho)
            plot_pdf(np.degrees(rho_val), pdf_rho, xlabel='$\\rho$ (degrees)', 
                        ylabel='Probability density', linecolour=clr)
        else:
            min_pdf_rho = np.amin(np.append(min_pdf_rho, np.amin(pdf_rho)))
            max_pdf_rho = np.amax(np.append(max_pdf_rho, np.amax(pdf_rho)))
            plot_pdf(np.degrees(rho_val), pdf_rho, linecolour=clr, overplot=True)
    plt.xlim(np.degrees(np.amin(rho_val)), np.degrees(np.amax(rho_val)))
    plt.ylim(min_pdf_rho, max_pdf_rho)
    plt.savefig('pdf_rho_rl.png')


# Finally, plot widths in degrees, making the Parkes widths (if used)
# a different symbol/colour than the GBT data.
   ### width_data = []
   ### width_data.append(read_widths(data_file))
   ### colour = ['red']
   ### mf_colour = ['red']

    ### if(use_parkes_data == True):
# Remember that we are adjusting parkes widths to match mean of GBT widths
       ### width_data.append(read_widths(parkes_data_file, rfile=data_file))
       ### colour.extend(['red'])
       ### mf_colour.extend(['white']) # open circles

    ### plot_widths(width_data, yunits='deg', colour=colour, \
        ###            mf_colour=mf_colour, me_colour=colour)

###    plt.savefig('widths_rl.png')



main()
