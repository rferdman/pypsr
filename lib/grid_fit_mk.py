#!/usr/local/bin/python

# Testing two different implementations of Levenberg-Marquardt in Python.  The 
# first is curve_fit, which is a wrapper for the scipy.optimize.leastsq 
# routine, and the second is lmfit-py, which uses leastsq as an engine, but is 
# moreflexible in that it allows for fixing and setting boundaries on 
# parameters.  I have been able to weight the fit by the measurement errors, 
# in both cases.

import numpy as np
from read_width import read_width_rad
from width_models import width_mk
from mjd import *
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters
from pdfplot import *

# global constants -- precession period and orbital inclination:
#from grid_setup_1913 import *
# from grid_setup_0737A import *

#prec_period, incl, mjd, y, y_err, mjd_mid

n_incl = 1 # Set to 2 if you want to combine probs from incl AND 180-incl
large_res = 1.0e8  # set residual to this if argument inside sqrt is <=0

#define model to be fit.  First up is Rafikov and Lai
"""Given certain geometrical parameters, determine the cosine of the
half-pulse widths at each epoch t.  Parameters are:

t = epoch, array of times
a - alpha, the angle between the pulsar spin and magnetic axes
l = lambda, the misalignment angle, between the pulsar spin vector and the
total orbital angular momentum
r = rho, the pulsar beam half-opening angle
T = T0, the epoch of zero precession phase

t and T0 must be in same units.

Based on equation (3) in Kramer (1998), ApJ, 509, 856."""

def residual(params, t, y, y_err, inclination, prec_period):
    alpha = params['alpha'].value
    lambd = params['lambda'].value
    rho = params['rho'].value
    T0 = params['T0'].value

    
#    print "phi = ", phi
#    print "cos_beta = ", cos_beta
#    print "sigma = ", sigma
#    print "sin_beta = ", np.sin(beta)
#    print "alpha = ",alpha
#    print "beta = ", beta
#    print "sin_alpha = ", np.sin(alpha)
#    print "sin_alpha * sin_beta = ", np.sin(alpha) * np.sin(beta)
#    print "sin(sigma/2) = ", np.sin(sigma/2.)
#    print "sin(rho/2) = ", np.sin(rho/2.)
#    print "incl = ", inclination
#    print "lambda = ", lamb
#    print "T0 = ", T0
#    print "rho = ", rho
#    print "numer = ", \
#        (np.sin(rho/2.))**2. - (np.sin(sigma/2.))**2.
#    print "denom = ", \
#        (np.sin(alpha) * np.sin(beta))
#    print "sqrt_arg = ", ( (np.sin(rho/2.))**2. - (np.sin(sigma/2.))**2. ) / \
#        (np.sin(alpha) * np.sin(beta))


#    print "W = ", W
#    print " "
#    print " "
#    print " "


#    cos_zeta = -np.cos(delta)*np.cos(inclination) + np.sin(delta)*np.sin(inclination)*np.cos(phi_SO)
#    sin_zeta = np.sqrt(1.0 - cos_zeta*cos_zeta)

#    cos_phi0 = (np.cos(rho) - cos_zeta*np.cos(alpha))/(sin_zeta*np.sin(alpha))
    

    # Divide by data uncertainty, so that minimizer will minimize
    # sum((y-cos_phi0)**2/y_err**2)

    W = width_mk(t, alpha, lambd, T0, rho, inclination, prec_period)
    
# First, if any of the values of the width array W are not finite (likely due 
# to negative argument inside the sqrt), then return a large array for the 
# residual so that the chi squared will be very large, making the pdf at that 
# point to be negligible
    finite_inds = np.where(np.isfinite(W))
#    neg_arg = np.where(sqrt_arg <= 0.)

    if(len(finite_inds[0]) == len(W) ):
# Divide by data uncertainty, so that minimizer will minimize
# sum((y-cos_phi0)**2/y_err**2)
        res =  (y-W)/y_err
    else:
        res = large_res * np.ones_like(W)

    return res
    
# Derivatives of above function
# def dfunc_rl(t, a, d, r, T):
    
    


# Main fitting routine.  Outside wrapper will pass a dictionary object
# containing all the necessary input parameters.
def grid_fit_mk(params_in, data_in):

# Expand input parmeter dictionary:
    n_alpha = params_in['n_alpha']
    n_lambda = params_in['n_lambda']
    n_T0 = params_in['n_T0']
    alpha_min = params_in['alpha_lim'][0]
    alpha_max = params_in['alpha_lim'][1]
    lambda_min = params_in['lambda_lim'][0]
    lambda_max = params_in['lambda_lim'][1]
    T0_min = params_in['T0_lim'][0]
    T0_max = params_in['T0_lim'][1]
    rho_min = params_in['rho_lim'][0]
    rho_max = params_in['rho_lim'][1]

    mjd = data_in['mjd']
    y = data_in['y']
    y_err = data_in['y_err']
    incl = data_in['incl']
    prec_period = data_in['prec_period']


# Now get arrays of grid steps for loops, and for later plotting.
# Shift each *_samples array by half the step size.  This will avoid having 
# zero values, which when run through the function, can give divisions by zero.
    alpha_samples, alpha_step = np.linspace(alpha_min, alpha_max, n_alpha, \
                                                endpoint=False, retstep=True) 
    alpha_samples = alpha_samples + alpha_step/2.
    lambda_samples, lambda_step = np.linspace(lambda_min, lambda_max, n_lambda, \
                                                endpoint=False, retstep=True) 
    lambda_samples = lambda_samples + lambda_step/2.
    T0_samples, T0_step = np.linspace(T0_min, T0_max, n_T0, \
                                                endpoint=False, retstep=True) 
    T0_samples = T0_samples + T0_step/2.

# Set up parameters with arbitrary values, and min/max boundary values.
    params=Parameters()
    params.add('alpha', value=np.radians(152.3), min=alpha_min, max=alpha_max)
    params.add('lambda', value=np.radians(12.4), min=lambda_min, max=lambda_max)
    params.add('rho', value=np.radians(9.), min=rho_min, max=rho_max)
    params.add('T0', value=yeartomjd(1990.), min=T0_min, max=T0_max)

# rho will vary in each fit: 
    params['rho'].vary = True
# The other parameters will be fixed at each grid point value:
    params['alpha'].vary = False
    params['lambda'].vary = False
    params['T0'].vary = False

# Set up 3-d array of ones to use for inputting resulting chi2's for each fit
    n_dims = (n_alpha, n_lambda, n_T0)

# Initialize array to keep fitted rho values for later histogram-ation
    rho_val = np.array([])
    rho_redchi2 = np.array([])

# Now run grid.  Opting not to run a list comprehension-type thing, for clarity.

    for i_incl in np.arange(n_incl):
# Just to be paranoid, reset chi2 array
        chi2 = np.ones(n_dims, dtype=float)
        redchi2 = np.ones(n_dims, dtype=float)

# alpha and lambda each will run from 0. to pi:
        for i_alpha in np.arange(n_alpha):
        # Set alpha value for current iteration of fit:
            params['alpha'].value = alpha_samples[i_alpha]

            for i_lambda in np.arange(n_lambda):
            # Set lambda value for current iteration of fit:
                params['lambda'].value = lambda_samples[i_lambda]

# T0 will run from mjd_mid-prec_period/2. to mjd_mid+prec_period/2.
                for i_T0 in np.arange(n_T0):
                # Set T0 value for current iteration of fit:
                    params['T0'].value = T0_samples[i_T0]

                # Do fit twice: one for inclination as is, and once for 180-incl
                    if(i_incl==0):
                        fit = minimize(residual, params, args=(mjd, y, y_err, \
                                                                   incl, prec_period) )
                    else:
                        fit = minimize(residual, params, args=(mjd, y, y_err, \
                                                               (np.pi-incl, prec_period)) )
                    
                    #chi2[i_alpha, i_lambda, i_T0] = \
                    #    np.sum(residual(params, mjd, y, y_err, incl)**2.)


                    chi2[i_alpha, i_lambda, i_T0] = fit.chisqr
                    redchi2[i_alpha, i_lambda, i_T0] = fit.redchi
 
# Append the fit value of rho to existing array
                    rho_val = np.append(rho_val, params['rho'].value)
                    rho_redchi2 = np.append(rho_redchi2, fit.redchi)
                   

                    # print to screen every 8 iterations of alpha (= 128*128*8)
                    if (np.fmod(i_alpha, 4)==0. and i_lambda==0 and i_T0==0):
                    #    chi2 = chisq(func_rl, mjd, y, popt, sigma=y_err)
                        print "Iter ", i_alpha * n_lambda * n_T0
                        print "Chisq = ", fit.chisqr
                        print "Reduced Chisq = ", fit.redchi
                        print 'Best-Fit Values:'
                        for name, par in params.items():
                            if (name=='T0'):
                                print '  %s = %.4f +/- %.4f ' % (name, par.value, par.stderr)
                            else:
                                print '  %s = %.4f +/- %.4f ' % (name, par.value*180./np.pi, par.stderr*180./np.pi)

# Append the fit value of rho to existing array


# Now have grid done.  Will need to now:
#     (1) Convert chi2 grid to an exp(-(chi2 - chi2_min)) grid.  Normalize to make
#         sum(3-d grid values) = 1

# Find minimum chi2 value:
        chi2_min = np.amin(chi2)
        redchi2_min = np.amin(redchi2)

        print "Min chisq = ", chi2_min
        print "Min reduced chisq = ", redchi2_min
        indices = np.where(chi2 == chi2_min)
        print "indices = ", indices
        print "min_alpha = ", alpha_samples[indices[0]]
        print "min_lambda = ", lambda_samples[indices[1]]
        print "min_T0    = ", T0_samples[indices[2]]

# create likelihood and normalize:
        if(i_incl==0):
            likelihood = np.exp(-(chi2 - chi2_min)/2.0)
        else:  # multiply together and combine probs from both inclination cases
            likelihood = likelihood*np.exp(-(chi2 - chi2_min)/2.0)

    norm_like = likelihood/np.sum(likelihood)

    print "Number of elements in likelihood array = ", np.size(likelihood)
    print "Sum of normalized array = ", np.sum(norm_like)

# Make up output dictionary to pass back to wrapping routine:
    p_out = {'norm_like':norm_like, \
                 'alpha':alpha_samples, 'lambda':lambda_samples, \
                 'T0':T0_samples, 'rho':rho_val, 'rho_redchi2':rho_redchi2}



    return p_out
