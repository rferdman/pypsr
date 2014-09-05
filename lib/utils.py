# Scripts for reading in various pulsar data  products and packaging them into dictionary objects
import sys
import subprocess
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings

# Given a set of parameters p = (A, mu, sigma), return the corresponding
# Gaussian value.
def gaussian(x, A, mu, sigma):
#     A, mu, sigma = p
     return A*np.exp(-(x-mu)**2/(2*sigma**2))
     
def straightline(x, m, b):
    return m*x + b

# Perform a fit of an x-y data set to a gaussian profile.
# Add in a chi-sq calculation as well to return back to the parent routine.
def fitgauss(xdata, ydata, yerr=None, p0=None):
     sig = None
     if(yerr!=None):
         sig = 1.0/yerr**2.0
     popt, pcov = curve_fit(gaussian, xdata, ydata, sigma=sig, p0=p0)
     A, mu, sigma = popt
     return A, mu, sigma
     
def fitline(xdata, ydata, yerr=None, p0=None):
    sig = None
    if(yerr!=None):
        sig = 1.0/yerr**2.0
    popt, pcov = curve_fit(straightline, xdata, ydata, sigma=sig, p0=p0)
    m, b = popt
    merr = np.sqrt(pcov[0,0])
    berr = np.sqrt(pcov[1,1])
    return m, b, pcov

# Pick out and return n random elements (non-repeating) of an array a, 
# keeping the order in which they appear in the original input array.
def pick(a, n):
     arr_inds = np.random.permutation(len(a))[0:n-1]
     arr_inds.sort()
     return arr_inds


def restart_line():
     sys.stdout.write('\r')
     sys.stdout.flush()
     

# Given an array of polynomial coefficients poly_coeffs, find all real 
# x value(s) corresponding to a given y value by root finding, within the 
# range [x1, x2]. 
# Returns a 1-D numpy array, which may have one or more elements.
def polyxval(poly_coeffs, y, x1=None, x2=None, warn=False):
     if (x1==None or x2==None):
          print "Must assign range [x1, x2] in which to search for x values."
          return None
     
     if(x1==x2):
          print "x1 must not equal x2."
          return None

     if (x1 > x2):
          temp_x = x1
          x1 = x2
          x2 = temp_x
     
# Supress warnings (default)
     if (warn==False):
          warnings.simplefilter('ignore', np.RankWarning)

     poly_coeffs_y = poly_coeffs
     # subtract y-value from zeroth order coefficient (i.e. no x's)
     poly_coeffs_y[len(poly_coeffs_y)-1] -= y
     re_roots = \
         np.roots(poly_coeffs_y)[np.roots(poly_coeffs_y).imag == 0.].real

# restrict solution to range [x1, x2]
     x_val = re_roots[(re_roots >= x1) & (re_roots <= x2)]

     return x_val
             
           
# Find the real roots of a polynomial within a given x range
def real_roots(poly_coeffs, x1=None, x2=None, warn=False):

# Supress warnings (default)
     if (warn==False):
          warnings.simplefilter('ignore', np.RankWarning)

# Evaluate roots, keeping only the real parts or those without an imaginary component
     re_roots = \
         np.roots(poly_coeffs)[np.roots(poly_coeffs).imag == 0.].real

# Go through limit possibilities, returning the appropriate values
# If no limits were given then return all real roots
     if (x1==None and x2==None):
          return re_roots
# The following are cases where either or both limits are given
     elif (x2==None):  # If only lower limit was given
          return re_roots[(re_roots >= x1)]
     elif (x1==None):  # If only upper limit was given
          return re_roots[(re_roots <= x2)]
     else:             # If both limits are given
          # Check that x1 < x2 and fix if necessary
          if (x1 > x2):
               temp_x = x1
               x1 = x2
               x2 = temp_x
          return re_roots[(re_roots >= x1) & (re_roots <= x2)]
          


def nearest_index(y, y_val, ind_tol=None):

     if(y_val==None):
          print "Please choose y_val keyword"
          sys.exit()

# Find the nearest profile ind on each side of the profile, corresponding to 
# the given pulse height
    # y_extra = np.append(y, y[0])
     y_bool = np.zeros_like(y).astype(bool)
     for i in np.arange(len(y_bool)-1): # -1 to avoid out of range array element for i+1:
         if (y[i] < y_val  and y[i+1] > y_val or y[i] > y_val  and y[i+1] < y_val):
             y_bool[i]= True
               
# These are the inds where the width at given y_val passes nearest:
     near_ind = np.where(y_bool)[0]
     print 'near_ind = ', near_ind

     if(ind_tol!=None):
# If there are several values very nearby to each other it is likely due to 
# slightly noisy profile, and so just take the first of the nearby values
# and shorten the array.
          ind_diff = np.array([near_ind[i+1] - near_ind[i] for i in np.arange(len(near_ind)-1)])
          ind_diff = np.append(ind_diff, np.abs(len(y)+near_ind[0]-near_ind[len(near_ind)-1]))
          near_ind = near_ind[ind_diff > ind_tol]
# specify [0] -- np.where output is tuple, and match inds to phase values
# just found
    # x_inds = np.where(prof_bool)[0][x_val_diff > PHASE_TOL]

     return near_ind

def bounding_index(x, x_val):

     n = nearest_index(x, x_val)
     if(x[n] > x_val):
          n_low = n-1
          n_high = n
     else:
          n_low = n
          n_high = n+1

     near_inds = np.array([n_low, n_high])
     return near_inds
     
# Simple command execution, returning standard output
def exec_cmd(cmd_str):
     cmd_str = cmd_str.split()
     proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
     proc_out = proc.communicate()[0]
     return proc_out

def div_check(num, denom):
  try: 
     float(num) / float(denom)
  except ZeroDivisionError:
     return True
  else:
     return False
