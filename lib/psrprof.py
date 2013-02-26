# Scripts for doing various pulsar-related profile manipulation.  Based on 
# passing data through dictionary objects

from sys import argv, exit
import math
import numpy as np
import scipy as sp
from utils import *
from psrread import get_duty
import matplotlib.pyplot as plt

import warnings


# Routine that creates a mask (of 1s and 0s) to denote which parts of the 
# profile are off and on-peak, respectively
def bin_mask(prof_data, duty):
 
#     thresh_factor = 1.8
     thresh_factor = 3.0

# Get various values and bins of inpute profile    
     n_bin = len(prof_data)
     prof_sort = prof_data[prof_data.argsort()]
     i_mid = int(math.floor(0.5*duty*n_bin+0.5))
     prof_mid = prof_sort[i_mid - 1]
     prof_peak = prof_sort[n_bin - 1]

# Get rms (std dev actually) of lowest i_mid bins of profile
# Remember in python the second index doesn't need the "-1" there
     rms = np.std(prof_sort[0:i_mid])
#     print "rms = ", rms
#     print "prof_mid = ", prof_mid

# Determine number of nearest neighbours to use
     n_check = n_bin/128
     if (n_check < 2):
          n_check = 2
     
     n_test = n_check/2 - 1

# Now cycle through each bin, compare it to its neighbours, and determine 
# whether it can be considered to be off-pulse
     
# First make a vector of 1s and 0s, depending on how each bin value compares 
# to the midpoint value
    # big = np.zeros(n_bin, dtype='int') # start with array of zeros
     mask = np.ones_like(prof_data, dtype=int)
#     print mask

     big_condition = ((prof_data - prof_mid) > (thresh_factor*rms))
#     print big_condition

     # First create array of profile array nearest neighbour elements, to see if
     # if each element will pass criteria for masking, and change 
     # corresponding mask array element if so.
     for i_bin in range(n_bin):
          if (i_bin < n_check):
               i_test = np.append(np.arange(n_bin-n_check+i_bin,n_bin), \
                                       np.arange(0,i_bin+n_check+1))
          elif (i_bin > n_bin-n_check-1):
               i_test = np.append(np.arange(i_bin-n_check,n_bin), \
                                  np.arange(0,n_check-(n_bin-i_bin)+1))
          else:          
               i_test = np.arange(i_bin-n_check,i_bin+n_check+1)



          # extract() elements of array where big_condition is true, get length
          # and compare to n_test.  If prof_data is not over threshole, then it
          # is part of off-pulse region.  Otherwise, mask it.

          # *Note* second argument to extract in this case can be anything with 
          # same dimensions as prof_data[i_test] or big_condition[i_test], 
          # since we are just after the number of elements in that sub-array
          # that pass the threshold condition.
          
          n_thresh = len(np.extract(big_condition[i_test], prof_data[i_test]))
          if(n_thresh > n_test):
               mask[i_bin] = 0
          
     return mask


def get_base(prof_data, duty):
     # get profile mask to determine off-pulse bins
     mask = bin_mask(prof_data, duty)
     # select those with mask==0, i.e. baseline bins
     baseline = prof_data[mask==1]
     # get mean and rms of baseline
     base_mean = np.mean(baseline)
     base_rms = np.std(baseline)

     # return tuple consisting of mean and rms of baseline
     return base_mean, base_rms


def remove_base(prof_data, duty):
     # get profile mask to determine off-pulse bins
#     mask = bin_mask(prof_data, duty)

     # select those with mask==0, i.e. baseline bins
#     baseline = prof_data[mask]
     # get mean and rms of baseline
#     base_mean = np.mean(baseline)

     baseline, base_rms = get_base(prof_data, duty)

     # remove baseline mean from profile in place
     prof_data = prof_data - baseline

     return prof_data


# Normalise a profile to have values between 0 and 1, in place. 
# One input profile at a time.
def norm(prof_data, duty):
     # get baseline mean and rms
     base_mean, base_rms = get_base(prof_data, duty)

     prof_peak = np.amax(prof_data)
     snr = (prof_peak - base_mean)/base_rms
     
     prof_norm = (prof_data - base_mean)/(prof_peak - base_mean)

     # return SNR     
     return prof_norm


# Will take in a profile and determine the width and width error of the 
# profile at the given PERCENTAGE height through bootstrapping a polynomial
# fit and interpolation using n_interp_pts points around the phase bin 
# corresponding to the given height.  The error on eahc side of the profile
# will be added in quadrature.
def get_width(prof_data, psr_name, percent_height, x_peak=None, \
                   n_pts_fit=12, n_order=6, n_omit=3, n_test_fit=4, \
                   n_boot=16384, hist_bins=64, return_more=False):

     PHASE_TOL = 0.005
     # Ignore warnings about polynomial rank:
     warnings.simplefilter('ignore', np.RankWarning)

# Remove baseline and normalize profile:
#     mask = bin_mask(prof_data['i'], duty)
#     prof_data['i'] = norm(prof_data['i'], duty)

# Find the peak value of the profile, and determine the height in absolute 
# terms. Assume we are given a baselined-out profile.
#   TO DO: interpolate for this as well
# First, remove baseline in order to get the height calculated correctly
     n_bin_prof = len(prof_data['i'])
     duty = get_duty(psr_name)
     prof_data['i'] = remove_base(prof_data['i'], duty)
#     x_peak, y_peak = get_peak(prof_data, n_pts_fit=10, n_order=19, n_test_fit=4)
     x_peak, y_peak = get_peak(prof_data, x_peak=x_peak, 
                               n_pts_fit=12, n_order=16, n_test_fit=5)
     print "x_peak, y_peak = ", x_peak, y_peak
     if(x_peak < 0):
          print 'Could not choose one.  Exiting: '
          exit()

     height = 0.01*percent_height*y_peak
#     peak_height = np.amax(prof_data['i'])
#     height = 0.01*percent_height*peak_height

############# do this another time ###########################################
# If there is one section to the baseline, then:
#   - if baseline wraps, profile doesn't
#   - if baseline is continuous, then profile is wrapping around
#
# If there is more than one section to the baseline, for now, just assume that 
# The profile is not wrapped and is rotated to the correct phase.  If not, do 
# it yourself!  
#   - Come up with an to deal with this case algorithm at a later time
##############################################################################

###################################################################################
     
# Find the nearest profile bin on each side of the profile, corresponding to 
# the given pulse height
###     prof_extra = np.append(prof_data['i'], prof_data['i'][0])
###     prof_bool = np.zeros_like(prof_data['i']).astype(bool)
###     for i in np.arange(len(prof_bool)):
###          if (prof_extra[i] < height  and prof_extra[i+1] > height \
###                   or prof_extra[i] > height  and prof_extra[i+1] < height):
###               prof_bool[i]= True

# These are the x phases of the profile where the width at given height passes
# nearesst:
###     x_val = prof_data['phase'][prof_bool]
# If there are several values very nearby to each other it is likely due to 
# slightly noisy profile, and so just take the first of the nearby values
# and shorten the array.
###     x_val_diff = np.array([x_val[i+1] - x_val[i] for i in np.arange(len(x_val)-1)])
###     x_val_diff = np.append(x_val_diff, np.abs(1.0+x_val[0]-x_val[len(x_val)-1]))
###     x_val = x_val[x_val_diff > PHASE_TOL]
# specify [0] -- np.where output is tuple, and match bins to phase values
# just found
###     x_bins = np.where(prof_bool)[0][x_val_diff > PHASE_TOL]
     x_bins = nearest_index(prof_data['i'], height, ind_tol=int(PHASE_TOL*n_bin_prof))
     n_width_bins = len(x_bins)

###################################################################################

# Now test quick slope of each point returned to determine is profile has wrapped.  
# If first slope is -ve and second is +ve, then we have a wrapped profile
     wrap = False
# Make test profile, adding wrapped points before and after, just in case.
     testprof = np.append(prof_data['i'][n_bin_prof-n_pts_fit:n_bin_prof], prof_data['i'])
     testprof = np.append(testprof, prof_data['i'][0:n_pts_fit])
     testphase = np.append(prof_data['phase'][n_bin_prof-n_pts_fit:n_bin_prof]- 1.0, prof_data['phase'])
     testphase = np.append(testphase, 1.0 + prof_data['phase'][0:n_pts_fit])
     if(n_width_bins == 2):
          if(testprof[x_bins[0]-n_pts_fit] > testprof[x_bins[0]+n_pts_fit] and \
                  testprof[x_bins[1]-n_pts_fit] < testprof[x_bins[1]+n_pts_fit]):
               wrap = True

     testbin = x_bins + n_pts_fit   
     print 'WRAP = ', wrap

# Now for each side of the profile, run a bootstrap of polynomial fits, 
# omitting n_omit points each time, and keeping track of the fit value for x at the given height
     
# May need an error check in here if either x_left or x_right gives a larger than one-value solution...
# ... perhaps use a while loop around fit/x-finding to keep doing until there is only one x-value...
# ... and set a limit to some number of tries, after which routine gives up.

# First set up set of all points from which to sample, on each side of profile:
     x_pts_left = testphase[testbin[0]-n_pts_fit : testbin[0]+n_pts_fit+1]
     y_pts_left = testprof[testbin[0]-n_pts_fit : testbin[0]+n_pts_fit+1]
     x_pts_right = testphase[testbin[len(testbin)-1]-n_pts_fit : testbin[len(testbin)-1]+1+n_pts_fit]
     y_pts_right = testprof[testbin[len(testbin)-1]-n_pts_fit : testbin[len(testbin)-1]+1+n_pts_fit]

# Initialize array of solutions for each side.  These are arrays since there may be 
# more than one solution, in which case we can flag an error:
     x_left = np.array([])
     x_right = np.array([])
# Main bootstrap loop:
     for i_boot in np.arange(n_boot):
# Left
          # choose indices at random
          inds = pick(x_pts_left, len(x_pts_left)-1-n_omit)
          x_pts_fit_left = x_pts_left[inds]
          y_pts_fit_left = y_pts_left[inds]
          poly_coeff_left = np.polyfit(x_pts_fit_left, y_pts_fit_left, n_order)
# In order to find the roots for y = height, subtract the height value from the final 
# coefficient, which is the zeroth-order coefficient
# As bounds, use point before and after in phase where y=height occurs, but now 
# re-sample y using the polynomial we just found, to get a smoother sampling:
#          x_left_bins = nearest_index(np.polyval(poly_coeff_left, x_pts_left), height)
          x_fit_left = polyxval(poly_coeff_left, height, \
          # using 0th element in x_left_bins since there is (hopefully) only one
#                                     x1=x_pts_left[x_left_bins[0]-n_test_fit], \
#                                     x2=x_pts_left[x_left_bins[0]+n_test_fit])
                                     x1=prof_data['phase'][x_bins[0]-n_test_fit], \
                                     x2=prof_data['phase'][x_bins[0]+n_test_fit])

# Right.  Same procedure:
          inds = pick(x_pts_right, len(x_pts_right)-1-n_omit)
          x_pts_fit_right = x_pts_right[inds]
          y_pts_fit_right = y_pts_right[inds]
          poly_coeff_right = np.polyfit(x_pts_fit_right, y_pts_fit_right, n_order)
#          x_right_bins = nearest_index(np.polyval(poly_coeff_right, x_pts_right), height)
          x_fit_right = polyxval(poly_coeff_right, height, \
#                            x1=x_pts_right[x_right_bins[0]-n_test_fit], \
#                            x2=x_pts_right[x_right_bins[0]+n_test_fit])
                            x1=prof_data['phase'][x_bins[len(x_bins)-1]-n_test_fit], \
                            x2=prof_data['phase'][x_bins[len(x_bins)-1]+n_test_fit])

# Choose solutions where there is only one x value per side:
          if(len(x_fit_left) == 1 and len(x_fit_right) == 1): # and \
# and where polynomial fit gives a sensible fit, i.e. only one value per side should
# be nearest to requested height for smooth fit:
#                  len(x_left_bins) == 1 and len(x_right_bins) == 1):
               x_left = np.append(x_left, x_fit_left)
               x_right = np.append(x_right, x_fit_right)
          else:
               if(len(x_fit_left) == 0):
                    print "Found zero solutions for x_left at iteration ", i_boot, \
                        ". Omitting from bootstrap."
          #               exit()
               elif(len(x_fit_left) > 1):
                    print "x_left has more than one value!  Found x_left = ", x_fit_left, \
                        " at iteration ", i_boot, ". Omitting from bootstrap."
#               exit()
                    
#          if(len(x_fit) == 1):
#               x_right = np.append(x_right, x_fit)
               if(len(x_fit_right) == 0):
                    print "Found zero solutions for x_right at iteration ", i_boot, \
                        ". Omitting from bootstrap."
#               exit()
               elif(len(x_fit_right) > 1):
                    print "x_right has more than one value!  Found x_right = ", x_fit_right, \
                        " at iteration ", i_boot, ". Omitting from bootstrap."
#               exit()
#               if(len(x_left_bins) != 1):
#                    print 'Left side has', len(x_left_bins), 'nearest bins to height:', x_left_bins, 'at iteration', i_boot, '. Omitting from bootstrap.'
#               if(len(x_right_bins) != 1):
#                    print 'Right side has', len(x_right_bins), 'nearest bins to height:', x_right_bins, 'at iteration', i_boot, '. Omitting from bootstrap.'




# Do histogram of widths:
     x_widths = x_right - x_left
     if (wrap):
          x_widths = 1.0 - (x_right - x_left)

# Perhaps weight by chisq of fit?
     hist_width, bins_width = np.histogram(x_widths, bins=72, normed=True, weights=None)

# Now fit a Gaussian to the width histogram, determine mean, median and width of distribution
     pdf_width = hist_width/np.sum(hist_width)
# Get bin size and add half to each bins_width, since the latter array is the edges of the hist bins
     bin_size = bins_width[1] - bins_width[0]
# Don't include the last (bin edge) + (half bin size)/2.
     bin_centres = bins_width[0:len(bins_width)-1] + (bin_size/2.)

     A, mu, sigma = fitgauss(bin_centres, pdf_width, p0=[np.max(pdf_width), np.average(bin_centres, weights=pdf_width), np.std(x_widths)])

# Possible additions:
     # add ability to plot pdf and gaussian fit of widths, or to return pdf/bin arrays and A/mu/sigma

     if (return_more):
          return mu, sigma, A, pdf_width, bin_centres
     else:
          return mu, sigma


# Simple routine to fit the peak of a profile using a polynomial fit, the using the derivative of
# that fit to determine the maximum around the roughly estimated maximum.  This is done by first
# finding the roots of the derivative, then evaluating the original fit polynomial at the found 
# x-value to get the y-value.  
# Default is an order-2 polynomial to approximate a quadratic around the peak.  This may break 
# down, so user is free to choose order.  No error estimation or bootstrap here.  May include that
# later on.
def get_peak(prof_data, x_peak=None, n_pts_fit=8, n_order=3, n_test_fit=2, return_more=False, warn=False):

# Supress warnings (default)
     if (warn==False):
          warnings.simplefilter('ignore', np.RankWarning)

     n_bin_prof = len(prof_data['i'])

     i_peak = np.argmax(prof_data['i'])
     rough_peak_x = prof_data['phase'][i_peak]
     rough_peak_y = prof_data['i'][i_peak]
     print 'ROUGH PEAK = ', [rough_peak_x, rough_peak_y]

# Make test profile, adding wrapped points before and after, just in case.
##     testprof = np.append(prof_data['i'][n_bin_prof-n_pts_fit:n_bin_prof], prof_data['i'])
##     testprof = np.append(testprof, prof_data['i'][0:n_pts_fit])
##     testphase = np.append(prof_data['phase'][n_bin_prof-n_pts_fit:n_bin_prof]- 1.0, prof_data['phase'])
##     testphase = np.append(testphase, 1.0 + prof_data['phase'][0:n_pts_fit])

# Each bin in testphase/testprof must be advanced by n_pts_test to match the actual points we are fitting
##     i_peak_test = i_peak + n_pts_fit

# Now choose points to perform polynomial fit; take into account possible
# phase wrap in used part of profile
     print 'i_peak = ', i_peak
     print 'n_pts_fits  = ', n_pts_fit
     print 'i_peak - n_pts_fit = ', i_peak - n_pts_fit 
     print 'max phase = ', np.max(prof_data['phase'])
     # In all cases, bring peak region we are using to ~0.5 in phase.
     # It seems that polyfit's results are affected when close to zero 
     # for some reason...     
     artificial_offset = (0.5 - rough_peak_x)
     if(i_peak-n_pts_fit < 0):
          print 'BELOW ZERO'
          n_overlap = np.abs(i_peak-n_pts_fit)
          x_peak_pts = np.append(
               prof_data['phase'][n_bin_prof - n_overlap : n_bin_prof] - 1.0,
               prof_data['phase'][0:i_peak+n_pts_fit+1]) + \
               artificial_offset
          y_peak_pts = np.append(
               prof_data['i'][n_bin_prof - n_overlap : n_bin_prof],
               prof_data['i'][0:i_peak+n_pts_fit+1])
     elif(i_peak + n_pts_fit >= n_bin_prof):
          print 'ABOVE N_BIN_PROF'
          n_overlap = np.abs(i_peak + n_pts_fit - n_bin_prof)
          x_peak_pts = np.append(
               prof_data['phase'][i_peak-n_pts_fit:n_bin_prof],
               prof_data['phase'][0:n_overlap])+ \
               artificial_offset
          y_peak_pts = np.append(
               prof_data['i'][i_peak-n_pts_fit:n_bin_prof],
               prof_data['i'][0:n_overlap])
     else:
          print 'NORMAL'
          x_peak_pts = prof_data['phase'][i_peak-n_pts_fit : i_peak+n_pts_fit+1]+ \
               artificial_offset
          y_peak_pts = prof_data['i'][i_peak-n_pts_fit : i_peak+n_pts_fit+1]
     
     print 'n_pts   = ', len(x_peak_pts)
     print 'n_bin   = ', n_bin_prof
     print 'n_order = ', n_order
# Perform polynomial fit
     poly_peak = np.polyfit(x_peak_pts, y_peak_pts, n_order)

# fit peak may have shifted from data peak (still based on data points):
     y_poly_test = np.polyval(poly_peak, x_peak_pts)
     i_new_peak = np.argmax(y_poly_test)
     x_new_peak = x_peak_pts[i_new_peak]
#     plt.plot(prof_data['phase'],prof_data['i'])
     print 'x_peak_pts  = ', x_peak_pts
#     print 'y_peak_pts  = ', y_peak_pts
#     print 'y_poly_test = ', y_poly_test
     plt.plot(x_peak_pts, y_peak_pts, 'o')
     plt.plot(x_peak_pts, y_poly_test)
#     plt.show()
     
# Find derivative of polynomial we just fit.  Default is 1st derivative
     poly_peak_deriv = np.polyder(poly_peak)

# n_test_fit in units of phase
     n_test_phase = (float(n_test_fit)/(float(n_bin_prof)))#/np.amax(prof_data['phase'])))
# Find roots of derivative around peak, in range given by user
     print 'x1 = ', x_new_peak - n_test_phase
     print 'x2 = ', x_new_peak + n_test_phase
     print 'xmin = ', np.min(x_peak_pts)
     print 'xmax = ', np.max(x_peak_pts)
     if(x_peak==None):
          x_peak = real_roots(
               poly_peak_deriv, 
               x1=(x_new_peak - n_test_phase), 
               x2=(x_new_peak + n_test_phase))
#                             x1=prof_data['phase'][i_peak-n_test_fit], \
#                             x2=prof_data['phase'][i_peak+n_test_fit])
     else:
          # Don't forget to artificially offset the x_peak provided by user:
          # (put into array format)
          x_peak = np.array([x_peak + artificial_offset])

     if(len(x_peak) == 1):
          # Evaluate original polynomial at above fit value     
          y_peak = np.polyval(poly_peak, x_peak)    
          # return x_peak to original phase by removing offset originally added
          # to make it around 0.5 phase:
          x_peak = x_peak - artificial_offset
          plt.plot(x_peak, y_peak, 'ro')
     # All done! Return values:
          if (return_more):
               return x_peak, y_peak, poly_peak
          else:
               return x_peak, y_peak

     elif(len(x_peak) == 0):
          print "Found zero solutions for x_peak. "
          return -1, -1
     else:
          print "x_peak has more than one value!  Found x_peak = ", x_peak
          return -1, -1
