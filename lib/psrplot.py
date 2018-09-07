# Scripts for doing various pulsar-related plotting.  Based on passing data through dictionary objects

from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from datetime import date, datetime, timedelta
import mjd
from psrcalc import bin_data, get_mass_func
from pdfplot import get_prob_2D_levels


# This routine will bin resids, and recalculate the error in binned residuals.
# Will default to 64 bins unless given a bin size.
# Here x is the x-axis, which can be in one on many units; y and yerr are the 
# unput residual values and errors, respectively.
def bin_resids(x, y, yerr, binsize=None, weigh_x=True, even_bins=False):
     
     # Default will be to start with first point, bin all points within binsize after that 
     # point, then find the next point after that, and do the same, and so on.
     if(binsize==None):
         n_bin = 64
         binsize = (np.amax(x) + np.amin(x))/float(n_bin)
     else:
         n_bin = int((np.amax(x) + np.amin(x))/float(binsize))


     # Best to sort x, y, and y_err to start
     ind_sort = np.argsort(x)
     x = x[ind_sort]
     y = y[ind_sort]
     yerr = yerr[ind_sort]
     # Have weight aarray ready
     weight = 1.0/(yerr**2.) 
  
     xbin = []
     ybin = []
     ybinerr = []
     
     # If we are doing evenly divided bins, then do it this way
     if(even_bins==True):
         # Set up bin edge values (adding the endpoint value so that there is 
         # n_bins+1 values in the array:
         bin_edges = np.linspace(np.amin(x), np.amax(x), num=n_bin+1, 
                                 endpoint=False)                     
         bin_min = bin_edges[0:n_bin]
         bin_max = bin_edges[1:n_bin+1]

         for i_bin in np.arange(n_bin):
             bin_inds = (x>=bin_min[i_bin]) & (x<bin_max[i_bin])
             count = len(x[bin_inds])
             if(len(x[bin_inds]) > 0):
                 bin_weight = np.sum(weight[bin_inds])
                 res = np.sum(y[bin_inds]*weight[bin_inds])/bin_weight
                 res2sum = np.sum((y[bin_inds]**2.0)*weight[bin_inds])
                 chi2 = (res2sum - bin_weight*(res**2.0))/(count - 1)
                 ybin.append(res)
                 ybinerr.append(np.sqrt(1.0/weight[bin_inds])*np.sqrt(chi2)/bin_weight)
#                 ybin.append(np.sum(y[bin_inds]*weight[bin_inds])/bin_weight)
#                 ybinerr.append(np.sqrt(1.0/bin_weight))
                 if(weigh_x==True):
                     xbin.append(np.sum(x[bin_inds]*weight[bin_inds])/bin_weight)
          
         if(weigh_x==False):
             xbin = bin_edges[0:n_bin]+binsize/2.
         else:
             xbin = np.array(xbin)

     # Otherwise do it as we go along, which is the default behaviour
     else:
         i_first = 0
         # Run binning until we use up all the data
         while(i_first < len(x)):
             bin_inds = (x>=x[i_first]) & ( x<(x[i_first]+binsize) )
             count = len(x[bin_inds])
             print 'X = ', x[i_first], ' --> ', x[i_first]+binsize, ':  COUNT = ', count
             bin_weight = np.sum(weight[bin_inds])
             res = np.sum(y[bin_inds]*weight[bin_inds])/bin_weight
             res2sum = np.sum((y[bin_inds]**2.0)*weight[bin_inds])
             if(count > 1):
                 chi2 = (res2sum - bin_weight*(res**2.0))/float(count - 1)
             else:
                 chi2 = 1.0
             ybin.append(res)
             ybinerr.append(np.sqrt(bin_weight)*np.sqrt(chi2)/bin_weight)
             print '   RES = ', res, '   ERR = ', np.sqrt(bin_weight)*np.sqrt(chi2)/bin_weight
#             ybinerr.append(np.sqrt(1.0/bin_weight))
             if(weigh_x==True):
                 xbin.append(np.sum(x[bin_inds]*weight[bin_inds])/bin_weight)
             else: # do straight mean
                 xbin.append(np.mean(x[bin_inds]))
             # Set next first index from which to count binsize. Go to next element after
             # last index of previous condition.  No need to sort since we already sorted x
             i_first = np.where(bin_inds)[0][-1] + 1
         
         xbin = np.array(xbin)

     ybin = np.array(ybin)
     ybinerr = np.array(ybinerr)

     return xbin, ybin, ybinerr



# Routine to plot the residuals against date. 
# Can handle multiple data sets, but for now will default output plot
# to be a stacked plot in y, i.e. vertical only
def plot_resid(resid_data, info_plot=None, canvassize=None, 
               axis_limits=None, binsize=None, resoffset=0.,
               preres=False, xunits='year', yunits='us', 
               xticks=True, yticks=True, xlabel=True, ylabel=True, 
               sym='o', symsize=1.8, colour=None, csize=1.3, 
               xlim=None, ylim=None, figtext=None, gridlines=None,
               ticklabelsize=18, axislabelsize=18):
#               mjdlim=None, yearlim=None, orbphaselim=None):

     
     
     
     # First determine if certain keywords are one-dimensional.  If so, 
     # will need to list-ify them
     if(type(resid_data) is not list):
          resid_data = [resid_data]
     else:
          print 'len(resid_data) = ', len(resid_data)

     if(type(xunits) is list):
          # Case where we want to plot multiple xunits for same 
          # residual data set
          if(len(xunits) > 1 & len(resid_data) == 1):
               resid_data = resid_data*len(xunits)
          # Case where we want multiple residual data sets with same 
          # xunits.
          elif(len(xunits) == 1 & len(resid_data) > 1):
               xunits = xunits*len(resid_data)
          # Otherwise, will need lengths of xunits and resid_data to
          # be equal
          elif(len(xunits) != len(resid_data)):
               print 'plot_resid ERROR: xunits and resid_data must have same dimensions if they are both > 1 in length'
     else:
          # Case where we want one set of xunits, no matter the amount
          # of residual data sets
          xunits = [xunits]*len(resid_data)

     # Deal with yunits the same as xunits:
     if(type(yunits) is list):
          # Case where we want to plot multiple yunits for same 
          # residual data set
          if(len(yunits) > 1 & len(resid_data) == 1):
               resid_data = resid_data*len(yunits)
          # Case where we want multiple residual data sets with same 
          # yunits.
          elif(len(yunits) == 1 & len(resid_data) > 1):
               yunits = yunits*len(resid_data)
          # Otherwise, will need lengths of yunits and resid_data to
          # be equal
          elif(len(yunits) != len(resid_data)):
               print 'plot_resid ERROR: yunits and resid_data must have same dimensions if they are both > 1 in length'
     else:
          # Case where we want one set of yunits, no matter the amount
          # of residual data sets
          yunits = [yunits]*len(resid_data)

     # For xlabels, xticks, ylabel, yticks:
     #     - if they are not passed as lists, then default will be to 
     #       make one label/tick label for x,  and one label for y axes.
     #     - if they are passed as lists, then it is up to the user to
     #       ensure that they correctly correspond to the residual 
     #       data list

     if(type(xlabel) is list):
          if(len(xlabel) != len(resid_data)):
               print 'plot_resid ERROR: xlabel list must have same dimensions as resid_data'
               exit()
     else:
          xlabel = [xlabel]*len(resid_data)

     if(type(ylabel) is list):
          if(len(ylabel) != len(resid_data)):
               print 'plot_resid ERROR: ylabel list must have same dimensions as resid_data'
               exit()
     else:
          ylabel = [ylabel]*len(resid_data)
          


     if (type(xticks) is list):
          if(len(xticks) != len(resid_data)):
               print 'plot_resid ERROR: xticks list must have same dimensions as resid_data'
               exit()
     else:
          xticks = [xticks]*len(resid_data)

     if (type(yticks) is list):
          if(len(yticks) != len(resid_data)):
               print 'plot_resid ERROR: yticks list must have same dimensions as resid_data'
               exit()
     else:
          yticks = [yticks]*len(resid_data)
    
     # Default is to plot same info IDs on each plot in the same way.
     # However, may have different data sets with different info IDs; in
     # such a case, ensure dimensions are the same as resid_data
     if(info_plot != None):
          if(type(info_plot[0]) is list):
               if(len(info_plot) == 1):   
                    # replicate this list to have the same dimensions 
                    # as resid_data (need square brackets to have 
                    # n *lists*; otherwise will continue on same list)
                    info_plot = [info_plot]*len(resid_data)
               elif(len(info_plot) != len(resid_data)):
                    print 'plot_resid ERROR: info_plot keyword much has same dimensions as resid_data.'
               exit()
          else:
               info_plot = [[info_plot]]*len(resid_data)

     
     # Set axis limits.  If they are not given, default behaviour is to 
     # divide evenly in vertical direction
     # We assume that for default plotting, first data set is plotted
     # on the bottome, and last is on the top opf the canvas
     # Based on a full plot being [0.12, 0.1, 0.8, 0.85]
     if(axis_limits is None):
          axis_limits = []
          for i_plot in np.arange(len(resid_data)):
               x1 = 0.12
               xwidth = 0.8
               ywidth=0.85/len(resid_data)
               y1 = 0.12 + i_plot*ywidth
               axis_limits.append([x1, y1, xwidth, ywidth])
     else:
          if(type(axis_limits) is list):
               if(len(axis_limits) == 1 & len(resid_data)>1):
                    axis_limits = [axis_limits]*len(resid_data)
               elif(len(axis_limits) != len(resid_data)):
                    print 'plot_resid ERROR: If specifying multiple axis_limits to keyword, it must be a list of same dimensions as resid_data'
               exit()
          else:
               axis_limits = [axis_limits]
          
               
     


     # This doesn't change throughout, so set it up now:
     possible_xunits = ['mjd', 'mjd0', 'year', 'orbphase', 'serial']
     
     # Set up conversion factors for possible yunits
     yconv = {'s':1e-6, 'ms':1e-3, 'us':1, 'ns':1e3}


# Set up the plot:
     fig = plt.figure(figsize=canvassize)

     ax = []
     for i_plot in np.arange(len(resid_data)):

          res_data = resid_data[i_plot]

          ax.append(fig.add_axes(axis_limits[i_plot]))
          ax[i_plot].xaxis.set_tick_params(labelsize=16)
          ax[i_plot].yaxis.set_tick_params(labelsize=16)



          if(possible_xunits.count(xunits[i_plot]) == 0):
               print "There is no value "+xunits[i_plot]+ \
                    " in the data that we can plot on x axis. Exiting."
               exit()


     # mjd0 will subtract the nearest 100 of the smallest mjd value from the 
     # mjd arrays
          if(xunits[i_plot]=='mjd0'):
              min_mjd = np.amin(res_data['mjd'])
              mjdint = np.floor(min_mjd)
              res_data['mjd0'] = res_data['mjd'] - mjdint

          if(xunits[i_plot]=='year'):
              res_data['year'] = \
                   np.array([mjd.mjdtoyear(m) for m in res_data['mjd']])

          # Convert units based on yunits keyword
          res_data['res'] *= yconv[yunits[i_plot]]
          res_data['reserr'] *= yconv[yunits[i_plot]]


#         xmax = np.amax(res_data[xunits[i_plot]])+0.003
#          ymin = np.amin(res_data['res'] - res_data['reserr'])
#          ymax = np.amax(res_data['res'] + res_data['reserr'])
#          xspan = abs(xmax - xmin)
#          yspan = abs(ymax - ymin)

#          if(xlim==None):
#               xlim=(xmin, xmax)
#          if(ylim==None):
#               ylim=(ymin, ymax)

          if (xlabel[i_plot]):          
               if(xunits[i_plot]=='serial'):               
                    ax[i_plot].set_xlabel('Serial number', fontsize=axislabelsize, labelpad=12)
               elif(xunits[i_plot]=='orbphase'):
                    ax[i_plot].set_xlabel('Orbital phase', fontsize=axislabelsize, labelpad=12)           
               elif(xunits[i_plot]=='mjd'):
                    ax[i_plot].set_xlabel('MJD', fontsize=axislabelsize, labelpad=12)
               elif(xunits[i_plot]=='mjd0'):
                    ax[i_plot].set_xlabel('MJD - {:d}'.format(int(mjdint)), 
                                  fontsize=axislabelsize, labelpad=12)
               elif(xunits[i_plot]=='year'):
                    ax[i_plot].set_xlabel('Year', fontsize=axislabelsize, labelpad=12)
                    xmajorFormatter = FormatStrFormatter('%4d')
                    ax[i_plot].xaxis.set_major_formatter(xmajorFormatter)


          if (ylabel[i_plot]):
               if(yunits[i_plot]=='s'):
                    ax[i_plot].set_ylabel('Residuals (s)', fontsize=axislabelsize, labelpad=6)
               if(yunits[i_plot]=='ms'):
                    ax[i_plot].set_ylabel('Residuals (ms)', fontsize=axislabelsize, labelpad=6)
               if(yunits[i_plot]=='us'):
                    ax[i_plot].set_ylabel('Residuals ($\mu$s)', fontsize=axislabelsize, labelpad=6)
               if(yunits[i_plot]=='ns'):
                    ax[i_plot].set_ylabel('Residuals (ns)', fontsize=axislabelsize, labelpad=6)


          ax[i_plot].tick_params(labelsize=ticklabelsize, pad=10)

          if(xticks[i_plot]):
               for tick in ax[i_plot].xaxis.get_major_ticks():
                    tick.label1On = True
                    tick.label2On = False #top ticks
          else:
               for tick in ax[i_plot].xaxis.get_major_ticks():
                    tick.label1On = False
                    tick.label2On = False 


          if(yticks[i_plot]):
               for tick in ax[i_plot].yaxis.get_major_ticks():
                    tick.label1On = True
                    tick.label2On = False # right ticks
          else:
               for tick in ax[i_plot].yaxis.get_major_ticks():
                    tick.label1On = False
                    tick.label2On = False


          if(gridlines!=None):
               for ycoord in gridlines:
                    ax[i_plot].axhline(ycoord, linestyle='--', color='black', \
                                    linewidth=0.4)

        

     # To keep track of min and max x and y
          x_min = []
          x_max = []
          y_min = []
          y_max = []
     # Set things up for plotting, especially if there are many infos 
     # and colours:
          if (res_data['info'] != None):
               # Find common info numbers between the command-line 
               # requested info numbers, and those in the info file itself, 
               # and only plot these:
               # If no info array given, just plot all of them
               if(info_plot==None): 
                    info_common = res_data['info_val']
               else:
                    # call this something else rather than info_plot since 
                    # list may have different lengths in each element, 
                    # depending on what user wants to plot
                    # first, cast as numpy array
                    info_this_plot=np.array(info_plot[i_plot])
                    # This will work but sorts the result... need way to 
                    # preserve original order...
                    info_common = \
                        np.intersect1d(np.unique(info_this_plot), 
                                       res_data['info_val'])
                    print 'info_common = ', info_common
               # Set up colours depending on number of info numbers          
               n_info = len(info_common)
               for i_info in np.arange(n_info):
                    info_condition = \
                        res_data['info']==info_common[i_info]
                    info_ind = np.where(info_condition)
                    res_x = res_data[xunits[i_plot]][info_ind]
               # np.where(res_data['info']==res_data['info_val'][i_info])
                    res_y = res_data['res'][info_ind] + \
                            float(i_info)*resoffset # zero by default
                    res_err = res_data['reserr'][info_ind]
                    if(binsize!=None):
                         res_x, res_y, res_err = bin_resids(res_x, res_y, 
                                                          res_err, 
                                                          binsize=binsize, 
                                                          weigh_x=True)

                    # Keep track of min and max x and y values
                    x_min.append(np.amin(res_x))#[res_x > xlim[0]]))
                    x_max.append(np.amax(res_x))#[res_x < xlim[1]]))
                    y_min.append(np.amin(res_y - res_err))#[res_y > ylim[0]]))
                    y_max.append(np.amax(res_y + res_err))#[res_y < ylim[1]]))
                    if (n_info==1):
                         clr = 'black'
                    else:
                         if (colour==None):
                             clr = cm.jet(float(i_info)/float(n_info-1)) 
                         else:
                             if (len(colour) == n_info):
                                 clr = colour[i_info]
                                 print 'COLOUR = ', clr
                             else:
                                 print 'Error: Must use same number of colours as info numbers.  Exiting'
                                 return
                    ax[i_plot].plot(res_x, res_y, sym, 
                                    markersize=symsize, 
                                    markeredgecolor=clr, 
                                    markerfacecolor=clr)
                    ax[i_plot].errorbar(res_x, res_y, yerr=res_err, 
                                        fmt=None, capsize=csize, 
                                        ecolor=clr)
          else:
               n_info = 1
               if (colour==None):
                   clr = 'black'
               else:
                   clr = colour
               res_x = res_data[xunits[i_plot]]
               res_y = res_data['res']
               res_err = res_data['reserr']
               if(binsize!=None):
                    res_x, res_y, res_err = bin_resids(res_x, res_y, 
                                                     res_err, 
                                                     binsize=binsize,
                                                     weigh_x=True)

               ax[i_plot].plot(res_x, res_y, sym, markersize=symsize, \
                            markeredgecolor=clr, markerfacecolor=clr)
               ax[i_plot].errorbar(res_x, res_y, yerr=res_err, \
                                fmt=None, capsize=csize, ecolor=clr)
               x_min.append(np.amin(res_x))#[res_x > xlim[0]]))
               x_max.append(np.amax(res_x))#[res_x < xlim[1]]))
               y_min.append(np.amin(res_y - res_err))#[res_y > ylim[0]]))
               y_max.append(np.amax(res_y + res_err))#[res_y < ylim[1]]))

     # Now set limits based on min and max *plotted* data:
          x_min = np.amin(np.array(x_min))
          x_max = np.amax(np.array(x_max))
          y_min = np.amin(np.array(y_min))
          y_max = np.amax(np.array(y_max))

          if(xlim==None):
              ax[i_plot].set_xlim(x_min, x_max)
          else:
              ax[i_plot].set_xlim(xlim)
          
          if(ylim==None):
              ax[i_plot].set_ylim(y_min, y_max)
          else:
              ax[i_plot].set_ylim(ylim)
               

     # Adjust limits to have about 5% breathing room of the plotted limits 
     # on either side:
          #print 'GOT TO HERE'
          x_lim = ax[i_plot].get_xlim()
          x_buffer = 0.025*(x_lim[1]-x_lim[0])
          ax[i_plot].set_xlim(x_lim[0]-x_buffer, x_lim[1]+x_buffer)
          y_lim = ax[i_plot].get_ylim()
          y_buffer = 0.05*(y_lim[1]-y_lim[0])
          ax[i_plot].set_ylim(y_lim[0]-y_buffer, y_lim[1]+y_buffer)

          #print 'ylim = ', y_lim
          # Figure text must be a list of tuples: [(x, y, text), (x, y, text), ...]
          if(figtext!=None):
               for txt in figtext:
                    ax[i_plot].text(txt[0], txt[1], txt[2], 
                                    fontsize=10, \
                                    horizontalalignment='center', \
                                    verticalalignment='center')





# Plot a single profile to an output file
def plot_prof(prof_data, canvassize=None, xticks=True, yticks=True, \
                   xlabel=True, ylabel=True, linecolour='black', \
                   overlap=False, xlim=None, ylim=None, figtext=None, \
                   hgrid=False, vgrid=False):

     # If we have just the one profile, cast it as a one-element list 
     # to make things nice and general
     if(type(prof_data) is dict):
          prof_data = [prof_data] 
     xmin = np.min(prof_data[0]['phase'])
     xmax = np.max(prof_data[0]['phase'])
     ymin = np.min(prof_data[0]['i'])
     ymax = np.max(prof_data[0]['i'])
     xspan = abs(xmax - xmin)
     yspan = abs(ymax - ymin)
     

# Set up the plot:
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

     if (xlabel):
          ax.set_xlabel('Pulse phase', fontsize=18)
     if (ylabel):
          ax.set_ylabel('Flux density (arbitrary units)', fontsize=18)
       
     if(xticks==False):
          for tick in ax.xaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False
          for tick in ax.xaxis.get_minor_ticks():
               tick.label1On = False
               tick.label2On = False


     if(yticks==False):
          for tick in ax.yaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False
               tick.tick1On = False
               tick.tick2On = False
          for tick in ax.yaxis.get_minor_ticks():
               tick.label1On = False
               tick.label2On = False
               tick.tick1On = False
               tick.tick2On = False


     if(hgrid):
          ax.yaxis.grid(linestyle='--', color='black', \
                                    linewidth=0.4)
     if(vgrid):
          ax.xaxis.grid(linestyle='--', color='black', \
                                    linewidth=0.4)          

     for i_prf in range(len(prof_data)):
          prf = prof_data[i_prf] 
          if (type(linecolour) is list):
               col = linecolour[i_prf]
          else:
               col = linecolour
#         ax.plot(prof_data['phase'], prof_data['i'], color=linecolour)
          ax.plot(prf['phase'], prf['i'], color=col)

# Figure text must be a list of tuples: [(x, y, text), (x, y, text), ...]
     if(figtext!=None):
          for txt in figtext:
               ax.text(txt[0], txt[1], txt[2], fontsize=10, \
                            horizontalalignment='left', \
                            verticalalignment='center')
     
    #plt.xlim(xmin, xmax)
   # plt.ylim(ymin - 0.01*yspan, ymax + 0.01*yspan)
    #  plt.grid()
#    if(xticks and yticks):
#         plt.ticklabel_format(axis='both')
#    elif(xticks and not yticks):
#         plt.ticklabel_format(axis='x')
#    if(not xticks):
#         ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useOffset=False)) 
#    else:
#         plt.ticklabel_format(axis='')
       
#    plt.tick_params(labelsize=16)
        

#    if(xlabel):
#         plt.xlabel('Phase', fontsize=18)
#    if(ylabel):
#         plt.ylabel('Flux density (arbitrary units)', fontsize=18)



def plot_widths(width_data, canvassize=None, msize=None, \
                     xval='year', yunits='phase', \
                     xticks=True, yticks=True, xlabel=True, ylabel=True, \
                     sym='o', colour=None, mf_colour=None, me_colour=None, \
                     csize=3, xlim=None, ylim=None, \
                     figtext=None, gridlines=None,
                     ticklabelsize=18, axislabelsize=18):

     possible_xval = ['mjd', 'mjd0', 'year', 'serial']
     if(possible_xval.count(xval) == 0):
          print "There is not value "+xval+ \
              " in the data that we can plot on x axis. Exiting."
          exit()

     possible_yunits = ['phase', 'deg', 'rad', 'cos_phi']
     if(possible_yunits.count(yunits) == 0):
          print yunits+ " is not an option for y axis units. Exiting."
          exit()

          

     # If we have just the one set of residuals, cast it as a one-element
     # list to make things nice and general
     if(type(width_data) is dict):
          width_data = [width_data] 

# mjd0 will subtract the nearest 100 of the smallest mjd value from the 
# mjd arrays
     if(xval=='mjd0'):
          min_mjd = np.array([])
          for width in width_data:
               min_mjd = np.append(min_mjd, np.amin(width['mjd']))
          mjdint = np.floor(np.amin(min_mjd))
          for width in width_data:
               width['mjd0'] = width['mjd'] - mjdint
      
     if(xval=='year'):
          for width in width_data:
               # date_out = [mjd.mjdtodate(m) for m in width['mjd']]
               width['year'] = [mjd.mjdtoyear(m) for m in width['mjd']]
               # width['year'] = [d.year + d.day/365. + \
               #                    d.hour/(365.*24.) + \
               #                    d.minute/(365.*24.*60.) + \
               #                    d.second/(365.*24.*60.*60.) \
               #                    for d in date_out]
     
# Set up plot limits now
     xmin = np.amin(width_data[0][xval])-0.003
     xmax = np.amax(width_data[0][xval])+0.003
     ymin = np.amin(width_data[0]['width'] - width_data[0]['werr'])
     ymax = np.amax(width_data[len(width_data)-1]['width'] + \
                         width_data[len(width_data)-1]['werr'])
     xspan = abs(xmax - xmin)
     yspan = abs(ymax - ymin)
          

# Set up the plot:
     fig = plt.figure(figsize=canvassize)
     ax = fig.add_axes([0.12, 0.14, 0.86, 0.83])
     ax.xaxis.set_tick_params(labelsize=ticklabelsize, pad=8)
     ax.yaxis.set_tick_params(labelsize=ticklabelsize, pad=8)
     if(xlim==None):
          ax.set_xlim(xmin - 0.01*xspan, xmax + 0.01*xspan)
     else:
          ax.set_xlim(xlim)
     if(ylim==None):
          ax.set_ylim(ymin - 0.01*yspan, ymax + 0.02*yspan)
     else:
          ax.set_ylim(ylim)

     if (xlabel):          
          if(xval=='serial'):               
               ax.set_xlabel('Serial number', fontsize=axislabelsize, labelpad=12)
          elif(xval=='mjd'):
               ax.set_xlabel('MJD', fontsize=axislabelsize, labelpad=12)
          elif(xval=='mjd0'):
               ax.set_xlabel('MJD - {:d}'.format(int(mjdint)), fontsize=axislabelsize, labelpad=12)
          elif(xval=='year'):
               # Set formatting for years so that they have %d formatting:
               xmajorFormatter = FormatStrFormatter('%d')
               ax.set_xlabel('Year', fontsize=axislabelsize, labelpad=12)
               ax.xaxis.set_major_formatter(xmajorFormatter)
               
     if (ylabel):
          ax.set_ylabel('Pulse width (degrees)', fontsize=axislabelsize)#, labelpad=8)

     if(not xticks):
          for tick in ax.xaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     if(not yticks):
          for tick in ax.yaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     

     for i_width in range(len(width_data)):
          width = width_data[i_width]
          # Get colours
          if(colour!=None):
              if (type(colour) is list):
                  clr = colour[i_width]
              else:
                  clr = colour
          else:
              # Set up automated colours
              if (len(width_data)==1):
                  clr = 'black'
              else:
                  clr = cm.gist_heat(float(i_width)/float(len(width_data))) 
              
          if (type(mf_colour) is list):
               mf_clr = mf_colour[i_width]
          else:
               mf_clr = clr
          if (type(me_colour) is list):
               me_clr = me_colour[i_width]
          else:
               me_clr = clr
          #ax.plot(res[xval], res['res'], 'o', markersize=msize, color=col)
          if(gridlines!=None):
               for ycoord in gridlines:
                    ax.axhline(ycoord, linestyle='--', color='black', \
                                    linewidth=0.4)
# Change to appropriate units: 
          if (yunits=='deg'):
               # Set formatting for degrees so that they have correct formatting:
               ymajorFormatter = FormatStrFormatter('%5.1f')
               ax.yaxis.set_major_formatter(ymajorFormatter)
               y_plot = width['width']*360.
               yerr_plot = width['werr']*360.
          elif (yunits=='rad'):
               y_plot = width['width']*2.*np.pi
               yerr_plot = width['werr']*2.*np.pi               
          elif (yunits=='cos_phi'):
               y_plot = cos(width['width']/2.)
               yerr_plot = cos(width['werr']/2.)
          else: # Just in units of phase as given in data file
               y_plot = width['width']
               yerr_plot = width['werr']

          if(i_width==0):
               xmin = np.amin(width[xval])
               xmax = np.amax(width[xval])
               ymin = np.amin(y_plot-yerr_plot)
               ymax = np.amax(y_plot+yerr_plot)               
          else:
               xmin = np.amin(np.append(width[xval], xmin))
               xmax = np.amax(np.append(width[xval], xmax))
               ymin = np.amin(np.append(y_plot-yerr_plot, ymin))
               ymax = np.amax(np.append(y_plot+yerr_plot, ymax))

          xspan = abs(xmax - xmin)
          yspan = abs(ymax - ymin)



# Overplot error bars.  Use fmt=None to tell it not to plot points:
          ax.plot(width[xval], y_plot, sym, color=clr, mfc=mf_clr, mec=me_clr)
          ax.errorbar(width[xval], y_plot, yerr=yerr_plot, \
                           capsize=csize, fmt=None, ecolor=clr, \
                           markersize=msize)
          
     if(xlim==None):
          ax.set_xlim(xmin - 0.025*xspan, xmax + 0.025*xspan)
     else:
          ax.set_xlim(xlim)
     if(ylim==None):
          ax.set_ylim(ymin - 0.1*yspan, ymax + 0.1*yspan)
     else:
          ax.set_ylim(ylim)

# Figure text must be a list of tuples: [(x, y, text), (x, y, text), ...]
     if(figtext!=None):
          for txt in figtext:
               ax.text(txt[0], txt[1], txt[2], fontsize=10, \
                            horizontalalignment='center', \
                            verticalalignment='center')


#     plt.savefig('test_widths.png')

     return ax
#     plt.plot(resid_data['mjd'], resid_data['res'], 'o')
#     plt.savefig(plot_file)




def plot_m1m2(m1m2_file, plot_pk=None, m1gr=None, m2gr=None, 
              m1gr_err=None, m2gr_err=None,
              m1m2_contour=None,
              pk_label_coord=None,
              plot_sin1=False, parfile=None, tempo_version='tempo1',
              xlim=None, ylim=None, plot_inset=False,
              colour=None, line_style=None, m1m2_pbdot_uncorr=None):

     
     n_plot = len(plot_pk)
# Ensure that we are plotting more than 1 parameter:
     if(n_plot < 1):
          print 'plot_m1m2:  Must plot at least one PK parameter. Exiting...'   
     if(colour == None):
          clr = [cm.jet(float(i_plot)/float(n_plot-1)) for i_plot in range(n_plot)]
     else:
          clr = colour
          
     if(line_style == None):
         line_style='-'
 
# Now, read input m1m2.dat-style file:     
     try:
          f_m1m2 = open(m1m2_file, 'r')
     except IOError as (errno, strerror):
          if (errno == 2):  # file not found
               print "IOError ({0}): File".format(errno), m1m2_file, "not found."
          else:
               print "IOError ({0}): {1}".format(errno, strerror)
          return
        
     print "File", m1m2_file, "open."

# Read file.  Ordering of parameters is m1, then m2 for limiting values of each 
# PK parameter, in the following order:
#    M1 M2(OMDOT_LOW OMDOT_HI GAMMA_LOW GAMMA_HI PBDOT_LOW PBDOT_HI R_LOW R_HI S_LOW S_HI)

     m1m2_data = np.array([line.split() for line in f_m1m2.readlines()], dtype=float)
     n_rows = m1m2_data.shape[0]
     n_cols = m1m2_data.shape[1]
     m1 = m1m2_data[:,0]
     m2 = {'omdot':m1m2_data[:, 1:3], 'gamma':m1m2_data[:, 3:5], 
           'pbdot':m1m2_data[:, 5:7], 'r':m1m2_data[:, 7:9], 
           's':m1m2_data[:, 9:11]}
     f_m1m2.close()

# Set up labels for each PK parameter:
     pk_label={'omdot':'$\\dot{\\omega}$', 'gamma':'$\\gamma$', 
               'pbdot':'$\\dot{P_b}$', 'r':'$r$', 's':'$s$'}

# Read in uncorrected pbdot file if given
     if(m1m2_pbdot_uncorr != None):
          try:
               f_m1m2 = open(m1m2_pbdot_uncorr, 'r')
          except IOError as (errno, strerror):
               if (errno == 2):  # file not found
                    print "IOError ({0}): File".format(errno), m1m2_pbdot_uncorr, "not found."
               else:
                    print "IOError ({0}): {1}".format(errno, strerror)
               return
        
          print "File", m1m2_pbdot_uncorr, "open."

     # Read file.  Ordering of parameters is m1, then m2 for limiting values of each 
     # PK parameter, in the following order:
     #    M1 M2(OMDOT_LOW OMDOT_HI GAMMA_LOW GAMMA_HI PBDOT_LOW PBDOT_HI R_LOW R_HI S_LOW S_HI)

          m1m2_data = np.array([line.split() for line in f_m1m2.readlines()], dtype=float)
          n_rows = m1m2_data.shape[0]
          n_cols = m1m2_data.shape[1]
          m1_pbdot_uncorr = m1m2_data[:,0]
          m2_pbdot_uncorr = {'pbdot_uncorr':m1m2_data[:, 5:7]}
          f_m1m2.close()

     # Read in contour data, if file is given.  Must be same .npy format as used for m2mtot grid:
     if(m1m2_contour!=None):
         load_array = np.load(m1m2_contour+'_params.npy')
         # for the purposes of this routine, only need the following 
         # things in p_out
         p_out = {'m2':load_array[0],
                  'mtot':load_array[1],
                  'm1':load_array[2],
                  'm1_prob':load_array[3]}
         p_out['norm_like'] = np.load(m1m2_contour+'_prob.npy')
         m1_contour_data = p_out['mtot']-p_out['m2']
    #     m1lim = (args.mtotlim[0] - args.m2lim[0], args.mtotlim[1] - args.m2lim[1])
    
         prob_intervals = [0.683, 0.954]
         weights=1.0
     
     # x_val = m1_contour_data, y_val = p_out['m2'], p_out['norm_like']
          # Create levels at which to plot contours at each of the above intervals.  
          # Will not assume going in that Z values are normalized to total volume of 1. 
         contour_level = get_prob_2D_levels(p_out['norm_like'], prob_intervals)
     #         if (norm==True):
     #             z_val = (contour_data*weights)/np.sum(contour_data*weights)
     #         else:
         z_val = p_out['norm_like']
     
# Now start plotting, in order given by command-line choices of PK parameters:
     fig = plt.figure(figsize=(11,8.5))
     ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])

     for i_pk in range(n_plot):
          for i_lim in range(2):  # 2 for low and hi pkk param values of m2
               ax.plot(m1, m2[plot_pk[i_pk]][:,i_lim], 
                       linestyle=line_style, color=clr[i_pk])

     if(m1m2_pbdot_uncorr != None):
          for i_lim in range(2):  # 2 for low and hi pkk param values of m2
               ax.plot(m1_pbdot_uncorr, m2_pbdot_uncorr['pbdot_uncorr'][:,i_lim], 
                       linestyle='dashed', color=clr[plot_pk.index('pbdot')])
         

# Get overall (x,y) ranges for plot.  Needed only if xlim and ylim are not 
# provided:
     if(xlim == None):
          # Base xlim on multiple of error bar
          # min x is max between m1 values and 10*error, and max x is min 
          # between them.
          xlim = (max([min(m1), m1gr - m1gr_err - 200.*m1gr_err]), 
                  min([max(m1), m1gr + m1gr_err + 200.*m1gr_err]))
     if(ylim == None):
          m2_all=[]
          for pk in plot_pk:
               m2_all.append(m2[pk])
          m2_min = np.amin(np.array(m2_all))
          m2_max = np.amax(np.array(m2_all))
          ylim = (max([m2_min, m2gr - m2gr_err - 200.*m2gr_err]), 
                  min([m2_max, m2gr + m2gr_err + 200.*m2gr_err]))
         
# Calculate curve representing sini = 1:
     if(plot_sin1):
         if (parfile==None):
             print 'Require par file to calculate mass function.'
             print 'Will proceed without plotting sin(i) restriction.'
         else:
             f_mass = get_mass_func(parfile, tempo_ver=tempo_version)
             m2_sin1 = np.linspace(ylim[0], ylim[1], num=128)
             m1_sin1 = ((m2_sin1**3.0)/(f_mass**2.0)) - m2_sin1

# Now plot area forbidden by sini>1, if requested:
             ax.plot(m1_sin1, m2_sin1, linestyle='dashed', color='black')
         # Fill in below these values



     ax.set_xlim(xlim)
     ax.set_ylim(ylim)
     ax.set_xlabel('Pulsar mass (M$_\\odot$)', fontsize=20.0)
     ax.set_ylabel('Companion mass (M$_\\odot$)', fontsize=20.0)
     ax.tick_params(labelsize=18)         
     
     
# If asked for, label PK params plotted at user-given coords
     if(pk_label_coord!=None):
          for i_pk in range(n_plot):
               ax.text(pk_label_coord[i_pk][0], pk_label_coord[i_pk][1], 
                       pk_label[plot_pk[i_pk]],
                       fontsize=18,
                       horizontalalignment='center',
                       verticalalignment='center')

    # Plot m1m2 contour on main plot if requested:
     if(m1m2_contour!=None):
 

     # Now plot the pdf data
     # Need to redo with regiular contour to get a line bordering it...
          ax.contour(m1_contour_data, p_out['m2'], z_val, levels=contour_level, 
                        colors=('0.0','none'), linewidths=[0.6, 0.0])
     # Need to do with contourf so that it gets filled
          ax.contourf(m1_contour_data, p_out['m2'], z_val, levels=contour_level, 
                        colors=('0.35','none'))
                      
                        
                        
     else:
         # Plot GR-derived masses, if given:
          if(m1gr!=None and m2gr!=None):
              ax.plot(m1gr, m2gr, 'o', markersize=2.8, color='black')
              ax.errorbar(m1gr, m2gr, xerr=m1gr_err, yerr=m2gr_err, ecolor='black')
    


# If asked for, plot inset with zoom into area around actual masses:
     if(plot_inset):
          xlim_inset = (m1gr - 8.0*m1gr_err, m1gr + 8.0*m1gr_err)
          ylim_inset = (m2gr - 8.0*m2gr_err, m2gr + 8.0*m2gr_err)
          # Plot little box in main plot:
          inset_box = plt.Rectangle((xlim_inset[0], ylim_inset[0]), xlim_inset[1]-xlim_inset[0], ylim_inset[1]-ylim_inset[0], 
                                      fc='none', ec='black', linestyle='dotted')
          plt.gca().add_patch(inset_box)
          
          # For now, put it in top right corner:
          ax_inset = fig.add_axes([0.55, 0.57, 0.35, 0.36])
          for i_pk in range(n_plot):
               for i_lim in range(2):  # 2 for low and hi pkk param values of m2
                    ax_inset.plot(m1, m2[plot_pk[i_pk]][:,i_lim], 
                                  linestyle='-', color=clr[i_pk])
          if(m1m2_pbdot_uncorr != None):
               for i_lim in range(2):  # 2 for low and hi pkk param values of m2
                    ax_inset.plot(m1_pbdot_uncorr, m2_pbdot_uncorr['pbdot_uncorr'][:,i_lim], 
                                  linestyle='dashed', color=clr[plot_pk.index('pbdot')])
         
          ax_inset.set_xlim(xlim_inset)
          ax_inset.set_ylim(ylim_inset)
#          ax_inset.set_xlabel('$m_1$')
#          ax_inset.set_ylabel('$m_2$')

          if(m1m2_contour!=None):
         # Need to redo with regiular contour to get a line bordering it...
              ax_inset.contour(m1_contour_data, p_out['m2'], z_val, levels=contour_level, 
                            colors=('0.0','none'), linewidths=[1.5, 0.0])
         # Need to do with contourf so that it gets filled
              ax_inset.contourf(m1_contour_data, p_out['m2'], z_val, levels=contour_level, 
                            colors=('0.6','none'))
                      
          else:
              # Plot GR-derived masses, if given:
               if(m1gr != None and m2gr != None):
                     ax_inset.plot(m1gr, m2gr, 'o', markersize=3.8, color='black')
                     ax_inset.errorbar(m1gr, m2gr, xerr=m1gr_err, 
                                           yerr=m2gr_err, ecolor='black')

          ax_inset.tick_params(labelsize=11)



     return


          
def plot_dmx(dmx_data, canvassize=None, msize=None, 
             xval='year', 
             xticks=True, yticks=True, xlabel=True, ylabel=True, 
             sym='o', symsize=3.4, 
             colour='black', mf_colour=None, me_colour=None,
             csize=2, xlim=None, ylim=None,
             figtext=None, gridlines=None):

     possible_xval = ['mjd', 'mjd0', 'year', 'serial']
     if(possible_xval.count(xval) == 0):
          print "There is not value "+xval+ \
              " in the data that we can plot on x axis. Exiting."
          exit()

     # If we have just the one set of residuals, cast it as a one-element
     # list to make things nice and general
     if(type(dmx_data) is dict):
          dmx_data = [dmx_data] 

# mjd0 will subtract the nearest 100 of the smallest mjd value from the 
# mjd arrays
     if(xval=='mjd0'):
          min_mjd = np.array([])
          for dmx in dmx_data:
               min_mjd = np.append(min_mjd, np.amin(dmx['mjd']))
          mjdint = np.floor(np.amin(min_mjd))
          for dmx in dmx_data:
               dmx['mjd0'] = dmx['mjd'] - mjdint
      
     if(xval=='year'):
          for dmx in dmx_data:
               # date_out = [mjd.mjdtodate(m) for m in width['mjd']]
               dmx['year'] = [mjd.mjdtoyear(m) for m in dmx['mjd']]

# Set up plot limits now
##     xmin = np.amin(dmx_data[0][xval])-0.003
##     xmax = np.amax(dmx_data[0][xval])+0.003
##     ymin = np.amin(dmx_data[0]['width'] - dmx_data[0]['werr'])
##     ymax = np.amax(dmx_data[len(dmx_data)-1]['width'] + \
##                         dmx_data[len(dmx_data)-1]['werr'])
##     xspan = abs(xmax - xmin)
##     yspan = abs(ymax - ymin)
          
# Set up the plot:
     fig = plt.figure(figsize=canvassize)
     ax = fig.add_axes([0.15, 0.1, 0.8, 0.85])
     ax.xaxis.set_tick_params(labelsize=16)
     ax.yaxis.set_tick_params(labelsize=16)
##     if(xlim==None):
##          ax.set_xlim(xmin - 0.01*xspan, xmax + 0.01*xspan)
##     else:
##          ax.set_xlim(xlim)
##     if(ylim==None):
##          ax.set_ylim(ymin - 0.01*yspan, ymax + 0.02*yspan)
##     else:
##          ax.set_ylim(ylim)

     if (xlabel):          
          if(xval=='serial'):               
               ax.set_xlabel('Serial number', fontsize=18)
          elif(xval=='mjd'):
               ax.set_xlabel('MJD', fontsize=18)
          elif(xval=='mjd0'):
               ax.set_xlabel('MJD - {:d}'.format(int(mjdint)), fontsize=18)
          elif(xval=='year'):
               # Set formatting for years so that they have %d formatting:
               xmajorFormatter = FormatStrFormatter('%d')
               ax.set_xlabel('Year', fontsize=18)
               ax.xaxis.set_major_formatter(xmajorFormatter)
               
     if (ylabel):
          ax.set_ylabel('Delta DM (pc cm$^{-3}$)', fontsize=18)


     if(not xticks):
          for tick in ax.xaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     if(not yticks):
          for tick in ax.yaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     for i_dmx in range(len(dmx_data)):
          dmx = dmx_data[i_dmx]
          # Get colours
          if (type(colour) is list):
               clr = colour[i_dmx]
          else:
               clr = colour
          if (type(mf_colour) is list):
               mf_clr = mf_colour[i_dmx]
          else:
               mf_clr = mf_colour
          if (type(me_colour) is list):
               me_clr = me_colour[i_dmx]
          else:
               me_clr = me_colour
          #ax.plot(res[xval], res['res'], 'o', markersize=msize, color=col)

          if(gridlines!=None):
               for ycoord in gridlines:
                    ax.axhline(ycoord, linestyle='--', color='black', \
                                    linewidth=0.4)

          # Change to appropriate units: 
          # Set formatting for degrees so that they have correct 
          # formatting:
          ymajorFormatter = FormatStrFormatter('%6.3f')
          ax.yaxis.set_major_formatter(ymajorFormatter)
          y_plot = dmx['delta_dm']
          yerr_plot = dmx['delta_dm_err']
          
          if(i_dmx==0):
               xmin = np.amin(dmx[xval])
               xmax = np.amax(dmx[xval])
               ymin = np.amin(y_plot-yerr_plot)
               ymax = np.amax(y_plot+yerr_plot)               
          else:
               xmin = np.amin(np.append(dmx[xval], xmin))
               xmax = np.amax(np.append(dmx[xval], xmax))
               ymin = np.amin(np.append(y_plot-yerr_plot, ymin))
               ymax = np.amax(np.append(y_plot+yerr_plot, ymax))

          xspan = abs(xmax - xmin)
          yspan = abs(ymax - ymin)

# Overplot error bars.  Use fmt=None to tell it not to plot points:
          ax.plot(dmx[xval], y_plot, marker=sym, markersize=symsize,
                  linestyle='None',
                  color=clr, mfc=mf_clr, mec=me_clr)
          ax.errorbar(dmx[xval], y_plot, yerr=yerr_plot, 
                      capsize=csize, fmt=None, ecolor=clr, 
                      markersize=msize)
     if(xlim==None):
          ax.set_xlim(xmin - 0.025*xspan, xmax + 0.025*xspan)
     else:
          ax.set_xlim(xlim)
     if(ylim==None):
          ax.set_ylim(ymin - 0.1*yspan, ymax + 0.1*yspan)
     else:
          ax.set_ylim(ylim)

# Figure text must be a list of tuples: [(x, y, text), (x, y, text), ...]
     if(figtext!=None):
          for txt in figtext:
               ax.text(txt[0], txt[1], txt[2], fontsize=10, \
                            horizontalalignment='center', \
                            verticalalignment='center')

     return
