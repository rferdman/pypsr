#!/usr/local/bin/python

# Running tempo over a grid of cosi and m2 values to eventually plot probability
# contours.

import sys
import subprocess
import numpy as np
from mjd import *
# from scipy.optimize import curve_fit
# from lmfit import minimize, Parameters
from pdfplot import *
from utils import *


"""
Given certain parameters, perform a gridded fit over sin(i) and m2 using 
tempo.  Parameters are:

parfile_base = input parfile, assumes no M2 or SINI lines
timfile = file containing TOAs
m2_range = range of m2 values to grid
cosi_range = range of cosi values to grid
n_m2 = number of m2 values to grid
n_cosi = number of cosi (and thus sini) values to grid
"""
# tempo call.  Return chi squared of fit.
def tempo_chi2(parfile, timfile, tempo2=False):

    # Call tempo
    if(tempo2):
        tempo_cmd = 'tempo2 -f '+parfile+' '+timfile
        exec_out = exec_cmd(tempo_cmd).split()
        i_chi2 = exec_out.index('Chisqr/nfree')+2
        chi2_str = exec_out[i_chi2]
        chi2 = float(chi2_str[0:chi2_str.index('/')])
        redchi2 = float(exec_out[i_chi2+2])
        
    else:   
        tempo_cmd = 'tempo -f '+parfile+' '+timfile
        exec_out = exec_cmd(tempo_cmd)
        # Read last line of output to get chi_squared and reduced
        # chi_squared
        f_lis = open('tempo.lis', 'r')
        file_lines = f_lis.readlines()
        n_lines = len(file_lines)
        last_line = file_lines[n_lines-1].split()
        chi2 = float(last_line[1].replace('/', '')) # get rid of slash...
        redchi2 = float(last_line[4])
        f_lis.close()

    return chi2, redchi2
    
    
    

# Main fitting routine.  Outside wrapper will pass all the necessary input parameters.
def grid_fit_shapiro_tempo(parfile_base, timfile, 
                           m2_range=[1.1, 1.4], cosi_range=[0.2, 0.4],
                           n_m2=64, n_cosi=64, fmass=None, tempo2=False):

# Expand input parmeter dictionary:

# Now get arrays of grid steps for loops, and for later plotting.
# Shift each array of samples by half the step size.  This will avoid having 
# zero values, which when run through the function, can give divisions by zero.
# Also: sample in parts, based on the setup file parameters, usually used to 
# sample the "important" parts of each parameter space:
    m2, m2_step = \
        np.linspace(m2_range[0], m2_range[1], n_m2, endpoint=False, retstep=True) 
    m2 = m2 + m2_step/2.

    cosi, cosi_step = \
        np.linspace(cosi_range[0], cosi_range[1], n_cosi, endpoint=False, retstep=True) 
    cosi = cosi + cosi_step/2.

    # create array of corresponding sin(i) values to actually put into parfile:
    sini = np.sqrt(1.0 - cosi*cosi)

# Set up 2-d array of ones to use for inputting resulting chi2's for each fit
    n_dims = (n_m2, n_cosi)
    print 'n_m2, n_cosi = (', n_m2, n_cosi, ')'

    n_iter = n_m2*n_cosi

# Now run grid.  Opting not to run a list comprehension-type thing, for clarity.
    
# Header for standard output:
    sys.stdout.write('\n {0:<15s}  {1:>7s}  {2:>8s}  {3:>10s}   {4:>14s}  {5:>14s}\n'.format("Iter", "m2", "cosi(i)", "sin(i)", "Chi2", "Red. Chi2"))

# Just to be paranoid, reset chi2 array
    chi2    = np.ones(n_dims, dtype=float)
    redchi2 = np.ones(n_dims, dtype=float)
    #        weight  = np.ones(n_dims, dtype=float)
    
    # First, read in par file:
    # First, read in par file:
    #par_base_contents = []
    #f_par_base = open(parfile_base, 'r')
    #for par_line in f_par_base.readlines():
    #    if ((par_line.split()[0] != 'SINI') & (par_line.split()[0]!='M2')):
    #        par_base_contents.append(par_line.split())
    #par_base_contents = [par_line.split() for par_line in f_par.readlines()]
    #f_par_base.close()
# Read in base par file:
    f_par_base = open(parfile_base, 'r')
    par_base_contents = f_par_base.read()  # no arguments, so read entire file
    f_par_base.close()
    
# Set name for temporary par file on which to run tempo:
    parfile_temp = 'par_temp.par'

    if(fmass!=None):
        m1 = np.zeros(n_iter, dtype=float)
        m1_chi2 = np.zeros(n_iter, dtype=float)
# alpha and delta each will run from 0. to pi:
    for i_m2 in np.arange(n_m2):
        
        for i_cosi in np.arange(n_cosi):            
            # Add the M2 and SINI lines to the base par file to create a temporary
            # par file to be fed to tempo:
#            f_par = open(parfile_temp, 'w')
#            for par_line in par_base_contents:                
#                f_par.write(' '.join(par_line)+'\n')
            f_par = open(parfile_temp, 'w')
            f_par.write(par_base_contents)
            f_par.write('M2    {0}  0\n'.format(m2[i_m2]))
            f_par.write('SINI  {0}  0\n'.format(sini[i_cosi]))
            f_par.close()

            # Run tempo, get chi squared:            
            chi2_cur, redchi2_cur = tempo_chi2(parfile_temp, timfile, 
                tempo2=tempo2)
            chi2[i_m2, i_cosi] = chi2_cur
            redchi2[i_m2, i_cosi] = redchi2_cur


            i_iter = i_m2*n_cosi + i_cosi
            # array of corresponding values of m1 at each grid point, given mass function
            if(fmass!=None):
                m1[i_iter] = np.sqrt(((m2[i_m2]*sini[i_cosi])**3.)/fmass) - m2[i_m2]
                m1_chi2[i_iter] = chi2_cur

#####
# print to screen every 8 iterations of m2 (= 128*128*8)
#                  if (np.fmod(i_alpha, 4)==0. and i_delta==0 and i_T1==0):
            
            # Print a status every 32 iterations
#########                #    chi2 = chisq(func_rl, mjd, y, popt, sigma=y_err)
#            if(np.fmod(iteration, 64)==0):
            restart_line()
            sys.stdout.write('{0:<10d}[{1:>3d}%]  {2:8.4f}  {3:8.4f}  {4:8.4f}   {5:14.5g}  {6:14.5g}'.format(
                    i_iter, int(100*float(i_iter)/float(n_iter)), 
                    m2[i_m2], cosi[i_cosi], sini[i_cosi], 
                    chi2_cur, redchi2_cur))
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
    print "min_m2 = ", m2[indices[0]]
    print "min_cosi = ", cosi[indices[1]]
    print "min_sini    = ", sini[indices[1]]

    likelihood = np.exp(-(chi2 - chi2_min)/2.0)

    norm_like = likelihood/np.sum(likelihood)  # because not varying bin sizes
 
#    print "Number of elements in likelihood array = ", np.size(likelihood)
#    print "Sum of normalized array = ", np.sum(norm_like)

    #m1_prob is basically the same as likelihood, but 1D
    if(fmass != None):
        m1_prob = np.exp(-(m1_chi2 - chi2_min)/2.0)/np.sum(likelihood)

# Make up output dictionary to pass back to wrapping routine:
    p_out = {'norm_like':norm_like, 
             'm2':m2, 
             'cosi':cosi, 
             'sini':sini,
             'm1':m1,
             'm1_prob':m1_prob}

    return p_out



# Main fitting routine.  Outside wrapper will pass all the necessary input parameters.
def grid_fit_m2mtot_tempo(parfile_base, timfile, 
                          m2_range=[1.15, 1.45], mtot_range=[2.56, 2.58],
                          n_m2=64, n_mtot=64, tempo2=False):

# Expand input parmeter dictionary:

# Now get arrays of grid steps for loops, and for later plotting.
# Shift each array of samples by half the step size.  This will avoid having 
# zero values, which when run through the function, can give divisions by zero.
# Also: sample in parts, based on the setup file parameters, usually used to 
# sample the "important" parts of each parameter space:
    m2, m2_step = \
        np.linspace(m2_range[0], m2_range[1], n_m2, endpoint=False, retstep=True) 
    m2 = m2 + m2_step/2.

    mtot, mtot_step = \
        np.linspace(mtot_range[0], mtot_range[1], n_mtot, endpoint=False, retstep=True) 
    mtot = mtot + mtot_step/2.


# Set up 2-d array of ones to use for inputting resulting chi2's for each fit
    n_dims = (n_m2, n_mtot)
    print 'n_m2, n_mtot = (', n_m2, n_mtot, ')'

    n_iter = n_m2*n_mtot

# Now run grid.  Opting not to run a list comprehension-type thing, for clarity.
    
# Header for standard output:
    sys.stdout.write('\n {0:<15s}  {1:>7s}  {2:>8s}  {3:>10s}   {4:>14s}  {5:>14s}\n'.format("Iter", "m2", "mtot", "m1", "Chi2", "Red. Chi2"))

# Just to be paranoid, reset chi2 array
    chi2    = np.ones(n_dims, dtype=float)
    redchi2 = np.ones(n_dims, dtype=float)
#    m1      = np.ones(n_dims, dtype=float)
    #        weight  = np.ones(n_dims, dtype=float)
    
    # First, read in par file:
#    par_base_contents = []
#    f_par_base = open(parfile_base, 'r')
##    for par_line in f_par_base.readlines():
 #       if ((par_line.split()[0] != 'MTOT') & (par_line.split()[0]!='M2')):
 #           par_base_contents.append(par_line.split())
   #   f_par_base.close()

# Read in base par file:
    f_par_base = open(parfile_base, 'r')
    par_base_contents = f_par_base.read()  # no arguments, so read entire file
    f_par_base.close()
    
# Set name for temporary par file on which to run tempo:
    parfile_temp = 'par_temp.par'

# create array of corresponding m1 values:
    m1 = np.zeros(n_iter, dtype=float)
    m1_chi2 = np.zeros(n_iter, dtype=float)

# alpha and delta each will run from 0. to pi:
    for i_m2 in np.arange(n_m2):
        
        for i_mtot in np.arange(n_mtot):     
       
            # Add the M2 and SINI lines to the base par file to create a temporary
            # par file to be fed to tempo:
            f_par = open(parfile_temp, 'w')
            f_par.write(par_base_contents)
            f_par.write('M2    {0}  0\n'.format(m2[i_m2]))
            f_par.write('MTOT  {0}  0\n'.format(mtot[i_mtot]))
            f_par.close()

            # Run tempo, get chi squared:            
            chi2_cur, redchi2_cur = tempo_chi2(parfile_temp, timfile, 
                tempo2=tempo2)
            chi2[i_m2, i_mtot] = chi2_cur
            redchi2[i_m2, i_mtot] = redchi2_cur
#####
# print to screen every 8 iterations of m2 (= 128*128*8)
#                  if (np.fmod(i_alpha, 4)==0. and i_delta==0 and i_T1==0):
            
            # Print a status every 32 iterations
#########                #    chi2 = chisq(func_rl, mjd, y, popt, sigma=y_err)
            i_iter = i_m2*n_mtot + i_mtot
            # Calculate m1:
            m1[i_iter] = mtot[i_mtot] - m2[i_m2]
            m1_chi2[i_iter] = chi2_cur

#            if(np.fmod(iteration, 64)==0):
            restart_line()
            sys.stdout.write('{0:<10d}[{1:>3d}%]  {2:8.4f}  {3:8.4f}  {4:8.4f}   {5:14.5g}  {6:14.5g}'.format(
                    i_iter, int(100*float(i_iter)/float(n_iter)), 
                    m2[i_m2], mtot[i_mtot], m1[i_iter], 
                    chi2_cur, redchi2_cur))
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
    print "min_m2 = ", m2[indices[0]]
    print "min_mtot = ", mtot[indices[1]]
    min_m1_ind = np.where(m1_chi2 == chi2_min)
    print "min_m1    = ", m1[min_m1_ind]

    likelihood = np.exp(-(chi2 - chi2_min)/2.0)

    norm_like = likelihood/np.sum(likelihood)  # because not varying bin sizes
 
    m1_prob = np.exp(-(m1_chi2 - chi2_min)/2.0)/np.sum(likelihood)
 
#    print "Number of elements in likelihood array = ", np.size(likelihood)
#    print "Sum of normalized array = ", np.sum(norm_like)


# Make up output dictionary to pass back to wrapping routine:
    p_out = {'norm_like':norm_like, 
             'm2':m2, 
             'mtot':mtot, 
             'm1':m1,
             'm1_prob':m1_prob}

    return p_out
