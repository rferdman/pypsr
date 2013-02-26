# Scripts for doing various probability distribution-related plots

from sys import exit
import numpy as np
import matplotlib.pyplot as plt

# Will take pdf data and return values where certain fractions of the total
# probability are enclosed, given as an array by the user.
# Assuming that the x-axis values are in numerical order (and pdf values 
# correspond to these) so that no initial sorting needs to be done.

# In the case where user is after upper limit, this routine returns 
# an array of upper limits corresponding to each probability interval input 
# by the user, with same dimensions as that input vector.

# In the case where user is after lower limit, this routine returns 
# an array of lower limits corresponding to each probability interval input 
# by the user, with same dimensions as that input vector.

# In the case where the user request probability intervals (by NOT making True
# the upper and/or lower keywords in the function call), this routine will
# return the median x value (scalar), the minimum x value vector, 
# and the maximum x value vector, each corresponding to each probability 
# interval input by the user, with same dimensions as that input vector.

# NOTE that it is up to the user to provide the appropriate set of intervals
# and data arrays.

def get_pdf_prob(x, pdf, prob_intervals, norm=False, \
                     weights=None, upper=False, lower=False):

# If weights are None, assign them to ones, with the same shape as the 
# input z array:

    if (weights==None):
        weights = np.ones_like(x, dtype=float)
# If weights are given, ensure they are the same shape as z:
    else:
        if(weights.shape != x.shape):
            print 'Shape of weight array ', weights.shape, \
                ' does not match the input data array ', x.shape
            return None

    if(norm==True):
        pdf = (pdf*weights)/np.sum(pdf*weights)

# Create an array for the cumularive probability with same dimension as x
    cumul_prob = np.zeros_like(x)

# Normalize the pdf:

# Deal with three cases here:  upper limit, lower limit, 
# or probability intervals


# Start with upper limit:
    if(upper==True):
        
# Cannot ask for upper AND lower limit:
        if(lower==True):
            print "Cannot do upper AND lower limits.  Must call each case "+\
                "individually.  Will proceed with upper limit only."
            print " "
        
# Here we do intervals by adding probability starting from minimum x value, and 
# interpolating cumulative probability distribution at limit values we 
# are after.
        
        x_upper = np.zeros_like(prob_intervals)
        
        for i_x in np.arange(len(x)):
            cumul_prob[i_x] = np.sum(pdf[0:i_x]*weights[0:i_x])
  
        x_upper =  np.interp(prob_intervals, cumul_prob, x)

        return x_upper


# Next, we move to lower limit:    
    elif(lower==True):
        
# No need to worry about upper being True; have checked for both being 
# simulataneously True in the previous if statement.

        
# Here we do intervals by adding probability starting from minimum x value, and 
# interpolating cumulative probability distribution at limit values we 
# are after.
        
        x_lower = np.zeros_like(prob_intervals)
        
        for i_x in np.arange(len(x)):
            cumul_prob[i_x] = np.sum(pdf[0:i_x]*weights[0:i_x])
  
# Here we will look for points corresponding to 1 - limiting prob, since this
# is a lower limit and we are working from left to right:
# (NOTE that we are NOT dividing by 2 here!)
        exclude_intervals = 1. - prob_intervals    
        
        x_lower =  np.interp(exclude_intervals, cumul_prob, x)

        return x_lower


# If the user doesn't want upper or lower limits, we will do probability 
# intervals    
    else:

# Arrays of min/max of probability intervals
        x_min = np.zeros_like(prob_intervals)
        x_max = np.zeros_like(prob_intervals)


# array of probability to be excluded at each side for each of the intervals
# given by user
        exclude_intervals = (1. - prob_intervals)/2.
    
        for i_x in np.arange(len(x)):
       # interval = exclude_intervals[i_x]
            cumul_prob[i_x] = np.sum(pdf[0:i_x]*weights[0:i_x])

# Now interpolate for the minimum values of each interval:
        x_min = np.interp(exclude_intervals, cumul_prob, x)
        x_max = np.interp(1.-exclude_intervals, cumul_prob, x)

# Calculate the median x value of the distribution
    
# Median is halfway point in area under distribution 
        x_med = np.interp(0.5, cumul_prob, x)


# return the min's and the max's for each prob interval as two separate but corresponding arrays to the parent routine.
        return x_med, x_min, x_max


def plot_pdf(x_val, pdf_data, canvassize=None, xticks=True, yticks=True, \
             xlabel=None, ylabel=None, norm=False, weights=None, \
                 linecolour='black', \
                 prob_lines=None, prob_linecolour='black', \
                 prob_linestyle= ['dashed', 'dashdot', 'dotted'], \
                 xlim=None, ylim=None, figtext=None, \
                 hgrid=False, vgrid=False, \
                 overplot=False):

# Normalize the pdf:
    if(norm==True):
        # If weights are None, assign them to ones, with the same 
        # shape as the input z array:
        if (weights==None):
            weights = np.ones_like(pdf_data, dtype=float)
         # If weights are given, ensure they are the same shape as z:
        else:
            if(weights.shape != pdf_data.shape):
                print 'Shape of weight array ', weights.shape, \
                    ' does not match the input data array ', pdf.shape
                return None
        pdf_data = (pdf_data*weights)/np.sum(pdf_data*weights)



# Start by setting up lot limits.  Assuming 1-D array input:
    xmin = np.min(x_val)
    xmax = np.max(x_val)
    ymin = np.min(pdf_data)
    ymax = np.max(pdf_data)
    xspan = abs(xmax - xmin)
    yspan = abs(ymax - ymin)
    
# Set up the plot:
    if(overplot==False):
        fig = plt.figure(figsize=canvassize)
        ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        if(xlim==None):
            ax.set_xlim(xmin, xmax)
        else:
            ax.set_xlim(xlim)
        if(ylim==None):
            ax.set_ylim(ymin - 0.01*yspan, ymax + 0.02*yspan)
        else:
            ax.set_ylim(ylim)

        if (xlabel!=None):
            ax.set_xlabel(xlabel, fontsize=18)
        if (ylabel!=None):
            ax.set_ylabel(ylabel, fontsize=18)
       
        if(not xticks):
            for tick in ax.xaxis.get_major_ticks():
                tick.label1On = False
                tick.label2On = False

        if(not yticks):
            for tick in ax.yaxis.get_major_ticks():
                tick.label1On = False
                tick.label2On = False

        if(hgrid):
            ax.yaxis.grid(linestyle='--', color='black', \
                              linewidth=0.4)
        if(vgrid):
            ax.xaxis.grid(linestyle='--', color='black', \
                              linewidth=0.4)          

# Now plot the pdf data
    plt.plot(x_val, pdf_data, color=linecolour, linestyle='steps-mid')


# plot vertical lines if asked for by user.  Choose colour or linestyles:
    if(prob_lines != None):
        plt.vlines(prob_lines, np.min(pdf_data)-0.5, np.max(pdf_data)+0.5, \
                      color=prob_linecolour, linestyle=prob_linestyle)

# include text is provided by user
    if(figtext!=None):
        for txt in figtext:
            plt.text(txt[0], txt[1], txt[2], fontsize=10, \
                        horizontalalignment='center', \
                        verticalalignment='center',)


#    return ax
    


# Will take in array of desired contour probability intervals and return 
# corresponding z-levels.
def get_prob_2D_levels(z, prob_intervals, norm=False, \
                           n_steps=32, weights=None):
    
# If weights are None, assign them to ones, with the same shape as the 
# input z array:
    if (weights==None):
        weights = np.ones_like(z, dtype=float)
# If weights are given, ensure they are the same shape as z:
    else:
        if(weights.shape != z.shape):
            print 'Shape of weight array ', weights.shape, \
                ' does not match the input data array ', z.shape
            return None

# Do normalization
    if(norm==True):
        z = (z*weights)/np.sum(z*weights)
# Find maximum value of normalized 2D probability function
    z_max = np.amax(z)
# Set up contour_levels array
    contour_level = np.zeros_like(prob_intervals) # Ensures same dimenstions

# Initialize step size to half the maximum probability value, as well as 
# intial pdf value
    step_size = z_max/2.
    z_level = z_max - step_size

# Now, for each of the given contour levels, determine z level that corresponds to
# an enclosed probability of that contour level.  Do this by stepping around the 
# final value until we arrive at the step threshold.
    for i_prob in np.arange(len(prob_intervals)):
# Run through number of steps given by user, dividing in half each time        
        for i_step in np.arange(n_steps):
            test_ind = np.where(z >= z_level)
            test_prob = np.sum(z[test_ind]*weights[test_ind])
            step_size = step_size/2.
            if(test_prob > prob_intervals[i_prob]):
                z_level = z_level + step_size
            else:
                z_level = z_level - step_size
# Now reset step_size to half the current prob interval (e.g. 0.683, 0.954, 0.9973)
        step_size = z_level/2.
# Now that we have gone down to desired step threshold, set current z_level to
# the level corresponding to desired current interval in loop
        contour_level[i_prob] = z_level


    return contour_level


def plot_contour_pdf(x_val, y_val, contour_data, n_steps=32,\
                         norm=False, weights=None, \
                         canvassize=None, xticks=True, yticks=True, \
                         xlabel=True, ylabel=True, linecolour='black', \
                         xlim=None, ylim=None, figtext=None, \
                         hgrid=False, vgrid=False):

# If weights are None, assign them to ones, with the same shape as the 
# input z array:
    if (weights==None):
        weights = np.ones_like(contour_data, dtype=float)
# If weights are given, ensure they are the same shape as z:
    else:
        if(weights.shape != contour_data.shape):
            print 'Shape of weight array ', weights.shape, \
                ' does not match the input data array ', contour_data.shape
            return None

# Start by setting up lot limits.  Assuming 1-D array input:
    xmin = np.min(x_val)
    xmax = np.max(x_val)
    ymin = np.min(y_val)
    ymax = np.max(y_val)
    xspan = abs(xmax - xmin)
    yspan = abs(ymax - ymin)
    
# Set up the plot:
    fig = plt.figure(figsize=canvassize)
    ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    if(xlim==None):
        ax.set_xlim(xmin - 0.01*xspan, xmax + 0.02*xspan)
    else:
        ax.set_xlim(xlim)
    if(ylim==None):
        ax.set_ylim(ymin - 0.01*yspan, ymax + 0.02*yspan)
    else:
        ax.set_ylim(ylim)

    if (xlabel!=None):
        ax.set_xlabel(xlabel, fontsize=18)
    if (ylabel!=None):
        ax.set_ylabel(ylabel, fontsize=18)
       
    if(not xticks):
        for tick in ax.xaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False

    if(not yticks):
        for tick in ax.yaxis.get_major_ticks():
            tick.label1On = False
            tick.label2On = False

    if(hgrid):
        ax.yaxis.grid(linestyle='--', color='black', \
                          linewidth=0.4)
    if(vgrid):
        ax.xaxis.grid(linestyle='--', color='black', \
                          linewidth=0.4)          

    prob_intervals = [0.683, 0.954, 0.9973]

# Create levels at which to plot contours at each of the above intervals.  
# Will not assume going in that Z values are normalized to total volume of 1. 
    contour_level = get_prob_2D_levels(contour_data, prob_intervals, n_steps=n_steps)

    if (norm==True):
        z_val = (contour_data*weights)/np.sum(contour_data*weights)
    else:
        z_val = contour_data

# Now plot the pdf data
    ax.contour(x_val, y_val, z_val, levels=contour_level, \
                   colors=('red', 'blue', 'green'))

    if(figtext!=None):
        for txt in figtext:
            ax.text(txt[0], txt[1], txt[2], fontsize=10, \
                        horizontalalignment='center', \
                        verticalalignment='center',)


