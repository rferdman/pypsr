# Scripts for doing various pulsar-related plotting.  Based on passing data through dictionary objects

from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
from datetime import date, datetime, timedelta
import mjd


# Routine to plot the residuals against date.  For now, use MJD, but will 
# update it to plot against years.  Will give command line option for MJD 
# soon, though.
def plot_resid(resid_data, info_plot=None, canvassize=None, \
                    preres=False, xunits='mjd', yunits='us', \
                    xticks=True, yticks=True, xlabel=True, ylabel=True, \
                    sym='o', symsize=1.8, colour='black', csize=1.3, \
                    xlim=None, ylim=None, figtext=None, gridlines=None, \
                    mjdlim=None, yearlim=None, orbphaselim=None):

     possible_xunits = ['mjd', 'mjd0', 'year', 'orbphase', 'serial']
     if(possible_xunits.count(xunits) == 0):
          print "There is not value "+xunits+ \
              " in the data that we can plot on x axis. Exiting."
          exit()

# Set up x limits:

##     if(mjdlim==None):
##         if(yearlim!=None):
##               if(len(mjdlim) != 2):
##                    print 'Error: yearlim array must be length 2.'
##                    return
              # Check that the first element is smaller than the second:
##               if(yearlim[0] > yearlim[1]):
##                    tempyear=yearlim[0]
##                    yearlim[0]=yearlim[1]
##                    yearlim[1]=tempyear
##               mjdlim=[mjd.yeartomjd(yearlim[0]), mjd.yeartomjd(yearlim[1])]
##          else:
##               mjdlim=[np.amin(resid_data['mjd']), np.amax(resid_data['mjd'])]
##     else:
##          if(len(mjdlim) != 2):
##               print 'Error: mjdlim array must be length 2.'
##               return
          # Check that the first element is smaller than the second:
##          if(mjdlim[0] > mjdlim[1]):
##               tempmjd=mjdlim[0]
##               mjdlim[0]=mjdlim[1]
##               mjdlim[1]=tempmjd
##          if(yearlim!=None):
##               print 'Must choose either mjdlim OR yearlim.'
##               print 'Will use mjdlim.'
##               yearlim==None
##          if(xunits=='orbphase'):
##               print 'Chosen x coordinate is orbital phase.'
##               print 'Will not apply mjdlim.'
         
          
##     if(xunits=='mjd'):
##          xlim=(mjdlim[0], mjdlim[1])

# mjd0 will subtract the nearest 100 of the smallest mjd value from the 
# mjd arrays
     if(xunits=='mjd0'):
#          min_mjd = np.array([])
#          for res inresid_data:
#          min_mjd = np.append(min_mjd, np.amin(res['mjd']))
          min_mjd = np.amin(resid_data['mjd'])
          mjdint = np.floor(min_mjd)
          resid_data['mjd0'] = resid_data['mjd'] - mjdint
##          xlim=(mjdlim[0]-mjdint, mjdlim[1]-mjdint)
      
     if(xunits=='year'):
          #date_out = [mjd.mjdtodate(m) for m in resid_data['mjd']]
          resid_data['year'] = \
              np.array([mjd.mjdtoyear(m) for m in resid_data['mjd']])
##          xlim=(mjd.mjdtoyear(mjdlim[0]), mjd.mjdtoyear(mjdlim[1]))
     
##     if(xunits=='orbphase'):
##          if(orbphaselim==None):
##               xlim=(0., 1.)
##          else:
##               xlim=(orbphaselim[0], orbphaselim[1])

     xmin = np.amin(resid_data[xunits])-0.003
     xmax = np.amax(resid_data[xunits])+0.003
     ymin = np.amin(resid_data['res'] - resid_data['reserr'])
     ymax = np.amax(resid_data['res'] + resid_data['reserr'])
     xspan = abs(xmax - xmin)
     yspan = abs(ymax - ymin)

# Set up the plot:
     fig = plt.figure(figsize=canvassize)
     ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])
     ax.xaxis.set_tick_params(labelsize=16)
     ax.yaxis.set_tick_params(labelsize=16)
     if(xlim==None):
          xlim=(xmin, xmax)
#          ax.set_xlim(xmin, xmax)
#     else:
#          ax.set_xlim(xlim)
     if(ylim==None):
          ylim=(ymin, ymax)
#          ax.set_ylim(ymin - 0.01*yspan, ymax + 0.02*yspan)
#     else:
#          ax.set_ylim(ylim)

     if (xlabel):          
          if(xunits=='serial'):               
               ax.set_xlabel('Serial number', fontsize=18)
          elif(xunits=='orbphase'):
               ax.set_xlabel('Orbital phase', fontsize=18)           
          elif(xunits=='mjd'):
               ax.set_xlabel('MJD', fontsize=18)
          elif(xunits=='mjd0'):
               ax.set_xlabel('MJD - {:d}'.format(int(mjdint)), fontsize=18)
          elif(xunits=='year'):
               ax.set_xlabel('Year', fontsize=18)


     if (ylabel):
          ax.set_ylabel('Residuals ($\mu$s)', fontsize=18)

     if(not xticks):
          for tick in ax.xaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     if(not yticks):
          for tick in ax.yaxis.get_major_ticks():
               tick.label1On = False
               tick.label2On = False

     if(gridlines!=None):
          for ycoord in gridlines:
               ax.axhline(ycoord, linestyle='--', color='black', \
                               linewidth=0.4)

# To keep track of min and max x and y
     x_min = []
     x_max = []
     y_min = []
     y_max = []
# Set things up for plotting, especially if there are many infos and colours:
     if (resid_data['info'] != None):
          # Find common info numbers between the command-line requested info
          # numbers, and those in the info file itself, and only plot these:
          if(info_plot==None): # If no info array given, just plot all of them
               info_common = resid_data['info_val']
          else:
               info_plot=np.array(info_plot) # first, cast as numpy array
               info_common = \
                   np.intersect1d(np.unique(info_plot), resid_data['info_val'])
               print 'info_common = ', info_common
          # Set up colours depending on number of info numbers          
          n_info = len(info_common)
          for i_info in np.arange(n_info):
               info_condition = \
                   resid_data['info']==info_common[i_info]
               info_ind = np.where(info_condition)
#               print resid_data[xunits]
               res_x = resid_data[xunits][info_ind]
                  # np.where(resid_data['info']==resid_data['info_val'][i_info])
               res_y = resid_data['res'][info_ind]
               # Keep track of min and max x and y values
               x_min.append(np.amin(res_x[res_x > xlim[0]]))
               x_max.append(np.amax(res_x[res_x < xlim[1]]))
               y_min.append(np.amin(res_y[res_y > ylim[0]]))
               y_max.append(np.amax(res_y[res_y < ylim[1]]))
               res_err = resid_data['reserr'][info_ind]
               if (n_info==1):
                    clr = 'black'
               else:
                    clr = cm.jet(float(i_info)/float(n_info-1)) 
               ax.plot(res_x, res_y, sym, markersize=symsize, \
                            markeredgecolor=clr, markerfacecolor=clr)
               ax.errorbar(res_x, res_y, yerr=res_err, \
                                fmt=None, capsize=csize, ecolor=clr)
     else:
          n_info = 1
          clr = colour
          res_x = resid_data[xunits]
          res_y = resid_data['res']
          res_err = resid_data['reserr']
          #ax.plot(res[xunits], res['res'], 'o', markersize=msize, color=col)
# Overplot error bars.  Use fmt=None to tell it not to plot points:
###          ax.plot(resid_data[xunits], resid_data['res'], \
###                       sym, markersize=symsize, \
###                       color=clr)
###          ax.errorbar(resid_data[xunits], resid_data['res'], \
###                           yerr=resid_data['reserr'], \
###                           fmt=None, capsize=csize, ecolor=clr)         
          ax.plot(res_x, res_y, sym, markersize=symsize, \
                       markeredgecolor=clr, markerfacecolor=clr)
          ax.errorbar(res_x, res_y, yerr=res_err, \
                           fmt=None, capsize=csize, ecolor=clr)
          x_min.append(np.amin(res_x[res_x > xlim[0]]))
          x_max.append(np.amax(res_x[res_x < xlim[1]]))
          y_min.append(np.amin(res_y[res_y > ylim[0]]))
          y_max.append(np.amax(res_y[res_y < ylim[1]]))

# Now set limits based on min and max *plotted* data:
     x_min = np.amin(np.array(x_min))
     x_max = np.amax(np.array(x_max))
     y_min = np.amin(np.array(y_min))
     y_max = np.amax(np.array(y_max))
 
     ax.set_xlim(x_min, x_max)
     ax.set_ylim(y_min, y_max)
#     ax.set_xlim(np.amax(np.array([x_min, xlim[0]])), 
#                 np.amin(np.array([x_max, xlim[1]])))
#     ax.set_ylim(np.amax(np.array([y_min, ylim[0]])), 
#                 np.amin(np.array([y_max, ylim[1]])))

# Adjust limits to have about 5% breathing room of the plotted limits 
# on either side:
     xlim = ax.get_xlim()
     x_buffer = 0.025*(xlim[1]-xlim[0])
     ax.set_xlim(xlim[0]-x_buffer, xlim[1]+x_buffer)
     ylim = ax.get_ylim()
     y_buffer = 0.05*(ylim[1]-ylim[0])
     ax.set_ylim(ylim[0]-y_buffer, ylim[1]+y_buffer)

# Figure text must be a list of tuples: [(x, y, text), (x, y, text), ...]
     if(figtext!=None):
          for txt in figtext:
               ax.text(txt[0], txt[1], txt[2], fontsize=10, \
                            horizontalalignment='center', \
                            verticalalignment='center')



#     return ax
#     plt.plot(resid_data['mjd'], resid_data['res'], 'o')
#     plt.savefig(plot_file)

#     plt.savefig('test_resid.eps')








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
                            horizontalalignment='center', \
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
                     sym='o', colour='black', mf_colour=None, me_colour=None, \
                     csize=3, xlim=None, ylim=None, \
                     figtext=None, gridlines=None):

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
     ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])
     ax.xaxis.set_tick_params(labelsize=16)
     ax.yaxis.set_tick_params(labelsize=16)
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
          ax.set_ylabel('Pulse width (degrees)', fontsize=18)

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
          if (type(colour) is list):
               clr = colour[i_width]
          else:
               clr = colour
          if (type(mf_colour) is list):
               mf_clr = mf_colour[i_width]
          else:
               mf_clr = mf_colour
          if (type(me_colour) is list):
               me_clr = me_colour[i_width]
          else:
               me_clr = me_colour
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
              pk_label_coord=None,
              xlim=None, ylim=None, plot_inset=False):

     
     n_plot = len(plot_pk)
# Ensure that we are plotting more than 1 parameter:
     if(n_plot < 2):
          print 'plot_m1m2:  Must plot more than one PK parameter. Exiting...'
     clr = [cm.jet(float(i_plot)/float(n_plot-1)) for i_plot in range(n_plot)]

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

# Set up labels for each PK parameter:
     pk_label={'omdot':'$\.{\\omega}$', 'gamma':'$\\gamma$', 
               'pbdot':'$\.{P_b}$', 'r':'$r$', 's':'$s$'}

# Now start plotting, in order given by command-line choices of PK parameters:
     fig = plt.figure(figsize=(11,8.5))
     ax = fig.add_axes([0.12, 0.1, 0.8, 0.85])

     for i_pk in range(n_plot):
          for i_lim in range(2):  # 2 for low and hi pkk param values of m2
               ax.plot(m1, m2[plot_pk[i_pk]][:,i_lim], 
                       linestyle='-', color=clr[i_pk])

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
         

     ax.set_xlim(xlim)
     ax.set_ylim(ylim)
     ax.set_xlabel('Pulsar mass ($M_\\odot$)', fontsize=20.0)
     ax.set_ylabel('Companion mass ($M_\\odot$)', fontsize=20.0)
     ax.tick_params(labelsize=18)
# Plot GR-derived masses, if given:
     ax.plot(m1gr, m2gr, 'o', markersize=2.8, color='black')
     ax.errorbar(m1gr, m2gr, xerr=m1gr_err, yerr=m2gr_err, ecolor='black')
     
# If asked for, label PK params plotted at user-given coords
     if(pk_label_coord!=None):
          for i_pk in range(n_plot):
               ax.text(pk_label_coord[i_pk][0], pk_label_coord[i_pk][1], 
                       pk_label[plot_pk[i_pk]],
                       fontsize=18,
                       horizontalalignment='center',
                       verticalalignment='center')



# If asked for, plot inset with zoom into area around actual masses:
     if(plot_inset):
          xlim_inset = (m1gr - 8.0*m1gr_err, m1gr + 8.0*m1gr_err)
          ylim_inset = (m2gr - 8.0*m2gr_err, m2gr + 8.0*m2gr_err)
          # For now, put it in top right corner:
          ax_inset = fig.add_axes([0.55, 0.57, 0.35, 0.36])
          for i_pk in range(n_plot):
               for i_lim in range(2):  # 2 for low and hi pkk param values of m2
                    ax_inset.plot(m1, m2[plot_pk[i_pk]][:,i_lim], 
                                  linestyle='-', color=clr[i_pk])
          ax_inset.set_xlim(xlim_inset)
          ax_inset.set_ylim(ylim_inset)
#          ax_inset.set_xlabel('$m_1$')
#          ax_inset.set_ylabel('$m_2$')
     # Plot GR-derived masses, if given:
          ax_inset.plot(m1gr, m2gr, 'o', markersize=3.8, color='black')
          ax_inset.errorbar(m1gr, m2gr, xerr=m1gr_err, 
                            yerr=m2gr_err, ecolor='black')
          ax_inset.tick_params(labelsize=11)

          

               
     


     return
