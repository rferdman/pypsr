#!/usr/local/bin/python

# Testing two different implementations of Levenberg-Marquardt in Python.  The 
# first is curve_fit, which is a wrapper for the scipy.optimize.leastsq 
# routine, and the second is lmfit-py, which uses leastsq as an engine, but is 
# moreflexible in that it allows for fixing and setting boundaries on 
# parameters.  I have been able to weight the fit by the measurement errors, 
# in both cases.
import sys
import numpy as np
from read_width import read_width_cos_phi
from width_models import width_rl
from mjd import *
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
from pdfplot import *
from utils import *

# from grid_setup_1913 import *
#from grid_setup_0737A import *

# gloabel constants -- precession period and orbital inclination:
#from width_cos_phi import prec_period, incl, mjd, y, y_err, mjd_mid

n_incl = 1 # Set to 2 if you want to combine probs from incl AND 180-incl

#define model to be fit.  First up is Rafikov and Lai
"""Given certain geometrical parameters, determine the cosine of the
half-pulse widths at each epoch t.  Parameters are:

t = epoch, array of times
a - alpha, the angle between the pulsar spin and magnetic axes
d = delta, the misalignment angle, between the pulsar spin vector and the
total orbital angular momentum
r = rho, the pulsar beam half-opening angle
T - T1, the epoch of zero precession phase

t and T must be in same units.

Based on equation (17) in Rafikov and Lai (2006), ApJ, 641, 438."""

def residual(params, t, y, y_err, inclination, prec_period):
    alpha = params['alpha'].value
    delta = params['delta'].value
    rho = params['rho'].value
    T1 = params['T1'].value


    cos_phi0 = width_rl(t, alpha, delta, T1, rho, inclination, prec_period)
    
 #   print y
 #   print cos_phi0
 #   print y_err
    # Divide by data uncertainty, so that minimizer will minimize
    # sum((y-cos_phi0)**2/y_err**2)
    return np.abs(y-cos_phi0)/y_err
    
# Derivatives of above function
# def dfunc_rl(t, a, d, r, T):
    
    

# Main fitting routine.  Outside wrapper will pass a dictionary object
# containing all the necessary input parameters.
def grid_fit_rl(params_in, data_in):

# Expand input parmeter dictionary:
    n_alpha = params_in['n_alpha']
    n_delta = params_in['n_delta']
    n_T1 = params_in['n_T1']
    alpha_min = params_in['alpha_lim'][0]
    alpha_max = params_in['alpha_lim'][1]
    delta_min = params_in['delta_lim'][0]
    delta_max = params_in['delta_lim'][1]
    T1_min = params_in['T1_lim'][0]
    T1_max = params_in['T1_lim'][1]
    rho_min = params_in['rho_lim'][0]
    rho_max = params_in['rho_lim'][1]
    alpha_ranges = params_in['alpha_ranges']
    delta_ranges = params_in['delta_ranges']
    T1_ranges = params_in['T1_ranges']
    alpha_sampling = params_in['alpha_sampling']
    delta_sampling = params_in['delta_sampling']
    T1_sampling = params_in['T1_sampling']
    

    mjd = data_in['mjd']
    y = data_in['y']
    y_err = data_in['y_err']
    incl = data_in['incl']
    prec_period = data_in['prec_period']
    
    print 'Precession period:  {0:.3f} days'.format(prec_period)

# Now get arrays of grid steps for loops, and for later plotting.
# Shift each *_samples array by half the step size.  This will avoid having 
# zero values, which when run through the function, can give divisions by zero.
# Also: sample in parts, based on the setup file parameters, usually used to 
# sample the "important" parts of each parameter space:
    alpha_samples = np.array([])
    alpha_weights = np.array([])
    alpha_step    = np.array([])
#    sample_size = (alpha_ranges[:,1] - alpha_ranges[:,0]) / \
#            alpha_sampling.astype(np.float)
    for i_range in np.arange(len(alpha_sampling)):
        alpha_samp_temp, alpha_step_temp = \
            np.linspace(alpha_ranges[i_range, 0], alpha_ranges[i_range, 1], \
                            alpha_sampling[i_range], \
                            endpoint=False, retstep=True) 
        alpha_step = np.append(alpha_step, alpha_step_temp)
        alpha_samples = np.append(alpha_samples, \
                                      alpha_samp_temp + alpha_step_temp/2.)
        alpha_weights = np.append(alpha_weights, \
                          np.ones_like(alpha_samp_temp)*alpha_step_temp)
# To get ratio of each step size to the max step size, to get weight
    alpha_step_max = np.amax(alpha_step)
    alpha_weights = alpha_weights/alpha_step_max
    n_alpha = len(alpha_samples)
    print "ALPHA weights = ", alpha_weights

    delta_samples = np.array([])
    delta_weights = np.array([])
    delta_step    = np.array([])
    for i_range in np.arange(len(delta_sampling)):
        delta_samp_temp, delta_step_temp = \
            np.linspace(delta_ranges[i_range, 0], delta_ranges[i_range, 1], \
                            delta_sampling[i_range], \
                            endpoint=False, retstep=True) 
    
        delta_step = np.append(delta_step, delta_step_temp)
        delta_samples = np.append(delta_samples, \
                                      delta_samp_temp + delta_step_temp/2.)
        delta_weights = np.append(delta_weights, \
                          np.ones_like(delta_samp_temp)*delta_step_temp)
    delta_step_max = np.amax(delta_step)
    delta_weights = delta_weights/delta_step_max
    n_delta = len(delta_samples)
    print "DELTA weights = ", delta_weights
        
    T1_samples  = np.array([])
    T1_weights  = np.array([])
    T1_step     = np.array([])
    for i_range in np.arange(len(T1_sampling)):
        T1_samp_temp, T1_step_temp = \
            np.linspace(T1_ranges[i_range, 0], T1_ranges[i_range, 1],\
                            T1_sampling[i_range], \
                            endpoint=False, retstep=True) 
        T1_step = np.append(T1_step, T1_step_temp)
        T1_samples = np.append(T1_samples, \
                                   T1_samp_temp + T1_step_temp/2.)
        T1_weights = np.append(T1_weights, \
                          np.ones_like(T1_samp_temp)*T1_step_temp)
    T1_step_max = np.amax(T1_step)
    T1_weights = T1_weights/T1_step_max
    n_T1 = len(T1_samples)
    print "T1 weights = ", T1_weights

	
# Set up parameters with arbitrary values, and min/max boundary values.
    params=Parameters()
    params.add('alpha', value=1.5, min=alpha_min, max=alpha_max)
    params.add('delta', value=0.1, min=delta_min, max=delta_max)
    params.add('rho', value=np.radians(88.0), min=rho_min, max=rho_max)
    params.add('T1', value=yeartomjd(1990.), min=T1_min, max=T1_max)

# rho will vary in each fit: 
    params['rho'].vary = True
# The other parameters will be fixed at each grid point value:
    params['alpha'].vary = False
    params['delta'].vary = False
    params['T1'].vary = False

# Set up 3-d array of ones to use for inputting resulting chi2's for each fit
    n_dims = (n_alpha, n_delta, n_T1)
    print "n_alpha, n_delta, n_T1 = ", n_alpha, n_delta, n_T1

    n_iter = n_alpha*n_delta*n_T1

# Initialize array to keep fitted rho values for later histogram-ation
    rho_val = np.array([])
    rho_chi2 = np.array([])
    rho_wgt = np.array([])

# Now run grid.  Opting not to run a list comprehension-type thing, for clarity.
    
# Header for standard output:
    sys.stdout.write('\n {0:<15s}  {1:>7s}  {2:>8s}  {3:>10s}  {4:^21s}  {5:>14s}  {6:>14s}\n'.format("Iter", "alpha", "delta", "T1", "rho", "Chi2", "Red. Chi2"))

    for i_incl in np.arange(n_incl):
# Just to be paranoid, reset chi2 array
        chi2    = np.ones(n_dims, dtype=float)
        redchi2 = np.ones(n_dims, dtype=float)
        weight  = np.ones(n_dims, dtype=float)

# alpha and delta each will run from 0. to pi:
        for i_alpha in np.arange(n_alpha):
        # Set alpha value for current iteration of fit:
            params['alpha'].value = alpha_samples[i_alpha]

            for i_delta in np.arange(n_delta):
            # Set delta value for current iteration of fit:
                params['delta'].value = delta_samples[i_delta]

# T1 will run from mjd_mid-prec_period/2. to mjd_mid+prec_period/2.
                for i_T1 in np.arange(n_T1):
                # Set T1 value for current iteration of fit:
                    params['T1'].value = T1_samples[i_T1]

                    # Rest initial guess for rho1 and rho2 each fit:
                    params['rho'].value = np.radians(88.0)

                    # Do fit twice: one for inclination as is, and 
                    # once for 180-incl
                    if(i_incl==0):
                        fit = minimize(residual, params, \
                                           args=(mjd, y, y_err, \
                                                     incl, prec_period) )
                    else:
                        fit = minimize(residual, params, \
                                           args=(mjd, y, y_err, \
                                                     (np.pi-incl, prec_period)) )
                    

                    chi2[i_alpha, i_delta, i_T1] = fit.chisqr
                    redchi2[i_alpha, i_delta, i_T1] = fit.redchi
                    weight[i_alpha, i_delta, i_T1] = \
                    	alpha_weights[i_alpha] * \
                    	delta_weights[i_delta] * \
                    	T1_weights[i_T1]
# Append the fit value of rho, and its associated reduce chisq, 
# to existing array
                    rho_val = np.append(rho_val, params['rho'].value)
                    rho_chi2 = np.append(rho_chi2, fit.chisqr)
#####
                    # print to screen every 8 iterations of alpha (= 128*128*8)
  #                  if (np.fmod(i_alpha, 4)==0. and i_delta==0 and i_T1==0):

                    # Print a status every 32 iterations
#########                #    chi2 = chisq(func_rl, mjd, y, popt, sigma=y_err)
                    iteration = i_alpha*n_delta*n_T1 + i_delta*n_T1 + i_T1
                    if(np.fmod(iteration, 64)==0):
                        restart_line()
                        sys.stdout.write('{0:<10d}[{1:>3d}%]  {2:8.4f}  {3:8.4f}  {4:10.4f}  {5:8.4f} +/- {6:8.4f}  {7:14.2g}  {8:14.2g}'.format(iteration, int(100*float(iteration)/float(n_iter)), np.degrees(params['alpha'].value), np.degrees(params['delta'].value), \
                            params['T1'].value, np.degrees(params['rho'].value), \
                            np.degrees(params['rho'].stderr), fit.chisqr, fit.redchi))
                        sys.stdout.flush()
                        
#                        print "Iter ", i_alpha * n_delta * n_T1
#                        print "Chisq = ", fit.chisqr
#                        print "Reduced Chisq = ", fit.redchi
#                        print 'Best-Fit Values:'
#                        for name, par in params.items():
#                            if (name=='T1'):
#                                print '  %s = %.4f +/- %.4f ' % (name, par.value, par.stderr)
#                            else:
#                                print '  %s = %.4f +/- %.4f ' % (name, par.value*180./np.pi, par.stderr*180./np.pi)


# Now have grid done.  Will need to now:
#     (1) Convert chi2 grid to an exp(-(chi2 - chi2_min)) grid.  Normalize to make
#         sum(3-d grid values) = 1

# Find minimum chi2 value:
        chi2_min = np.amin(chi2)
        redchi2_min = np.amin(redchi2)

        print " "
        print "Min chisq = ", chi2_min
        print "Min reduced chisq = ", redchi2_min
        indices = np.where(chi2 == chi2_min)
        print "indices = ", indices
        print "min_alpha = ", alpha_samples[indices[0]]
        print "min_delta = ", delta_samples[indices[1]]
        print "min_T1    = ", T1_samples[indices[2]]

# create likelihood and normalize:
        if(i_incl==0):
            likelihood = np.exp(-(chi2 - chi2_min)/2.0)
        else:  # multiply together and combine probs from both inclination cases
            likelihood = likelihood*np.exp(-(chi2 - chi2_min)/2.0)

#    norm_like = likelihood/np.sum(likelihood)
    norm_vol = (likelihood * weight)/np.sum(likelihood * weight)
    bin_vol_max = alpha_step_max * delta_step_max * T1_step_max
    norm_like = likelihood/(bin_vol_max * np.sum(likelihood * weight))

    rho_prob = np.exp(-(rho_chi2 - chi2_min)/2.0)
    #rho_prob = rho_prob*rho_wgt/np.sum(rho_prob*rho_wgt)
#    rho_prob = rho_prob/np.sum(rho_prob)

    print "Number of elements in likelihood array = ", np.size(likelihood)
    print "Sum of normalized array = ", np.sum(norm_like)


# Make up output dictionary to pass back to wrapping routine:
    p_out = {'norm_like':norm_like, \
                 'norm_vol':norm_vol, \
                 'alpha':alpha_samples, \
                 'delta':delta_samples, \
                 'T1':T1_samples, \
#                 'alpha_weights':alpha_weights, \
#                 'delta_weights':delta_weights, \
#                 'T1_weights':T1_weights, \
                 'rho':rho_val, 'rho_prob':rho_prob}
    
    return p_out
