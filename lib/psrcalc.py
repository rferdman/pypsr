from sys import exit
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as cst
from astropy import coordinates as coord
from scipy.optimize import newton
import swiftmonitor.utils as smu
from psrread import read_par


#c_light = 299792458.
tsun_us = 4.925490947  # solar constant in microsecs

#cst.c.value
def get_pkparams(ma=None, mb=None, ecc=None, pb=None, x=None):
    # assume pb in given in days, x=asini in lt-s, and masses in solar masses
    tsun_us=4.925490947
    tsun_sec = tsun_us * 10**(-6)
    pb *= 86400.0
    
    omdot = 3.0 * (tsun_sec)**(2./3.) * ((pb/(2.0*np.pi)))**(-5./3.) * (1.0/(1.0 - ecc*ecc))*(ma + mb)**(2.0/3.0)
    gamma = (tsun_sec)**(2./3.) * ((pb/(2.0*np.pi)))**(1./3.) * ecc*mb*(ma + 2.0*mb)/((ma+mb)**(4.0/3.0))
    ge = (1.0 + (73.0/24.0)*(ecc**2.0) + (37.0/96.0)*(ecc**4.0)) * ((1.0 - ecc**2.0)**(-7.0/2.0))
    pbdot = (-192.0*np.pi/5.0)*(tsun_sec**(5.0/3.0)) * ((pb/(2.0*np.pi)))**(-5./3.) * ge * ma*mb/((ma+mb)**(1.0/3.0))
    r = mb
    sini = (tsun_sec)**(-1./3.) * ((pb/(2.0*np.pi)))**(-2./3.) * x * ((ma+mb)**(2.0/3.0))/mb
    
    # convert to normal units
    omdot = np.degrees(omdot) * 86400.0 * 365.25   # deg/year
    gamma *= 1000.0    # ms
    
    pk = {'omdot':omdot, 'gamma':gamma, 'pbdot':pbdot, 'r':r, 'sini':sini}
    
    
    return pk
    

def get_mass_func(parfile, tempo_ver='tempo1'):
    # read in pb and asini from par file, adn calculation mass function:
#    if(tempo_ver='tempo2'):
#        tempo_ver = 'tempo2'
    params=read_par(parfile, file_format=tempo_ver)
    # asini is in (light) seconds.  Convert pb to secs as well.
    # tsun is in us to give f_mass in solar units.
    tsun = 4.925490947
    tsun_secs = tsun*10**(-6)
    pb_secs = float(params['pb'])*86400.0
    f_mass = 4*np.pi*np.pi*(float(params['a1'])**3.)/(tsun_secs*(pb_secs**2.))
#    f_mass = 4*np.pi*np.pi*(float(params['a1'])**3.)/(tsun*(float(params['pb'])**2.))
    print 'pb = ', float(params['pb']), ' = ', pb_secs, ' sec'
    print 'a1 = ', float(params['a1'])
    print 'fmass = ', f_mass
    return f_mass

    

def get_prec_period(ma=None, mb=None, ecc=None, pb=None, units='years'):
    possible_units = ['centuries', 'century', 'cen', 'years', 'year', 'y', 'days', 'day', 'd']    
    if (possible_units.count(units)==0):
        print 'Can only output precession period in days, years, or centuries'
        exit()
    
    # Calculated in days (pb is input in days), so convert if needed 
    if(units[0:1]=='c'):
        conv_div = u.year.to(u.day, 100.0)
    elif(units[0:1]=='y'):
        conv_div = u.year.to(u.day, 1.0)
        #conv_div = 365.242
    else:  # days
        conv_div = 1.0
            
 
    tsun_days = tsun_us / 10.0**6 / 86400.   # solar constant converted to units of days
    
    # angular precession rate in rad/days
    prec_rate_days_a = ((2.0*np.pi/pb)**(5./3.)) * \
        ((tsun_days)**(2./3.)) * (mb*(4.0*ma + 3.0*mb)/(2.0*((ma + mb)**(4./3.) ) ) )*\
        ((1.0 - ecc**2)**(-1))
    
    prec_period = 2.0 * np.pi / prec_rate_days_a  # In days.
    
    prec_period /= conv_div
    
    return prec_period

    
def get_pbdot_gal(l, b, d, derr=None, v0=220.0, v0err=None, R0=7.7, R0err=None, 
                verbose=False):
    
    # l is in degrees
    # b is in degrees
    # d is in kpcm -- only convert after calculating z, which is based on d
    #   being in kpc.
    # v0 is in km/s
    # R0 is in kpc
    
    
    # change quantities to consistent units (m, radians)
    l = u.deg.to(u.radian, l)
    b = u.deg.to(u.radian, b)
    
    v0 *= 1000.0
    if(v0err != None):
        v0_err *= 1000.0
    R0 = u.kpc.to(u.m, R0)
    if(R0err != None):
        R0err = u.kpc.to(u.m, R0err)
    
    # First calculate component orthogonal to plane:
    z = d * np.sin(b)
    if(derr != None):
        z_err = derr * np.sin(b)
        
    if(verbose):
        print 'z = ', z
    
    # The following is based on Nice and Taylor (1995)
    az_c = 1.09e-19 * ( ( (1.25 * z)/np.sqrt(z**2. + 0.0324) ) + (0.58 * z) )
    
    pbdot_pb_orth = -az_c * np.sin(b)
    if(verbose):
        print 'pbdot/pb (z) = ', pbdot_pb_orth
    
    # Next calculate component due to centripetal acceleration due to planar
    # rotation
    d = u.kpc.to(u.m, d)
    if(derr != None):
        derr = u.kpc.to(u.m, derr)
    
    beta = (d/R0)*np.cos(b) - np.cos(l)
    
    pbdot_pb_planar = -np.cos(b)*(v0**2./(cst.c.value*R0)) * \
                    (np.cos(l) + ( beta/(np.sin(l)**2. + beta**2.) ) )
    
    if(verbose):
        print 'pbdot/pb (planar) = ', pbdot_pb_planar
    
    # Put them together!
    pbdot_pb_gal = pbdot_pb_orth + pbdot_pb_planar
    
    return pbdot_pb_gal
    
        
    
def get_shklovskii(pm, d, pmerr=None, derr=None):
    
    # pm is in mas/year
    # d is in kpc
    
    # change quantities to consistent units (rad/s, m)
    
    pm = (u.mas.to(u.radian, pm)) / u.year.to(u.s, 1.0)
    d = u.kpc.to(u.m, d)
    
    pbdot_pb_shk = pm**2. * d/cst.c.value
    
    return pbdot_pb_shk
    
    
def get_pbdot_corr(ra, dec, pm, d, pb, v0=220.0, R0=7.7, 
                pmerr=None, derr=None, v0err=None, R0err=None,
                verbose=False):
    
    # ra and dec are given in degrees
    # pm is given in mas/year
    # d is given in kpc
    # pb given is seconds
    
    # convert coordinates to degrees:
#    ra = 360.0 * (ra[0] + ra[1]/60. + ra[2]/3600.)/24.0
#    if(dec[0] < 0):
#        sgn = -1.
#    else:
#        sgn = 1.
#    dec =  sgn * (np.abs(dec[0]) + dec[1]/60. + dec[2]/3600.)
    # convert ra and dec to galactic coordinates, in degrees

    # Now convert ra and dec to l and b 
    c = coord.ICRSCoordinates(ra=ra, dec=dec, unit=(u.degree, u.degree))
    l = c.galactic.l.degrees
    b = c.galactic.b.degrees
        
    
    # convert pb to seconds
    # pb = u.day.to(u.s, pb)
    
    
    # Calculate corrections.  Returns pbdot/pb in each case.
    
    # calculate galactic correction -- assume v0 and R0 are fine as 
    # the default values:
    if(verbose):
        print 'Calculating Galactic contribution:'
    pbdot_pb_gal = get_pbdot_gal(l, b, d, v0=v0, R0=R0, 
                                derr=derr, v0err=v0err, R0err=R0err,
                                verbose=verbose)
    if(verbose):
        print 'Galactic contribution:  ', pbdot_pb_gal * pb, ' s/s\n'
    
    # calculate the contribution from the Shklovskii effect:
    if(verbose):
        print 'Calculating Shklovskii contribution:'
    pbdot_pb_shk = get_shklovskii(pm, d, pmerr=pmerr, derr=derr)
    if(verbose):
        print 'Shklovskii contrbution: ', pbdot_pb_shk * pb, ' s/s\n'
    
    # get total correction -- should now be unitless:
    pbdot_corr = (pbdot_pb_gal + pbdot_pb_shk) * pb
    if(verbose):
        print 'Total pbdot correction: ', pbdot_corr, ' s/s\n'
    
    return pbdot_corr
    
    
def get_dist_pbdot(ra, dec, pm, xpbdot, pb, v0=222.0, R0=7.95, dist_guess=2.4,
                verbose=False):
    
    # pm is given in mas/year, convert to radians per second
    # d is given in kpc
    # xpbdot is dimensionless
    # v0 is in km/s, convert to kpc/s
    # R0 is in kpc
    
    # Following the above conversions, d will pop out in kpc.
    
    # convert velocity to kpc
    # No need to convert R0, which is already in R0
    v0 = u.km.to(u.kpc, v0)

    # convert pm to rad/s
    pm = (u.mas.to(u.radian, pm)) / u.year.to(u.s, 1.0)

    # Everything is not converted to good units.  
    
    dist = scipy.optimize.newton(dist_func, dist_guess, 
                                args=(l, b, v0, R0, pm, pb, xpbdot))
 
    return dist
 

# Use Newton's method to get distance from pbdot correction equation, given 
# xpbdot=observed-GR pbdot

def dist_func(d, l, b, v0, R0, pm, pb, xpbdot): 
 
    # All distances must be in kpc because of the az_c component
 
    # First, convert speed of light to kpc/s
    c_kpc = u.m.to(u.kpc,cst.c.value)
    
    # Calculate galactic component
    # First calculate component orthogonal to plane:
    z = d * np.sin(b)
            
    # The following is based on Nice and Taylor (1995)
    az_c = 1.09e-19 * ( ( (1.25 * z)/np.sqrt(z**2. + 0.0324) ) + (0.58 * z) )
    
    pbdot_pb_orth = -az_c * np.sin(b)
        
    beta = (d/R0)*np.cos(b) - np.cos(l)
    
    pbdot_pb_planar = -np.cos(b)*(v0**2./(c_kpc*R0)) * \
                    (np.cos(l) + ( beta/(np.sin(l)**2. + beta**2.) ) )
    
    # Put them together!
    pbdot_pb_gal = pbdot_pb_orth + pbdot_pb_planar
 
    # For Shklovskii component:
    pbdot_pb_shk = pm**2. * d/c_kpc
 
    pbdot_corr = (pbdot_pb_gal + pbdot_pb_shk) * pb
 
    # Final function is xpbdot - pbdot_corr = 0
    f = xpbdot - pbdot_corr
    
    return f
    
    
    # Can solve for distance analytically if we use low-b approximation, 
    # so that a_z term --> 0 (technically depends on other variables as well,
    # but in the case of 1756-2251, for example, this is the case).
    # Otherwise, will need to do root finding... which I'll do later.
    
    
    
    
    
    
    
    
    
    
# This routine will bin x,y data, and recalculate the error in y.
# Will default to 64 bins unless given a bin size
def bin_data(x, y, yerr=None, binsize=None, weigh_x=False, even_bins=False):
     
     # Default will be to start with first point, bin all points within binsize after that 
     # point, then find the next point after that, and do the same, and so on.
     if(binsize==None):
         n_bin = 64
         binsize = (np.amax(x) - np.amin(x))/float(n_bin)
     else:
         n_bin = np.int_((np.amax(x) - np.amin(x))/float(binsize))

     # Best to sort x, y, and y_err to start
     ind_sort = np.argsort(x)
     x = x[ind_sort]
     y = y[ind_sort]
     if(yerr != None):
         y_err = yerr[ind_sort]
     else:
         y_err = np.ones_like(y)
     # Have weight aarray ready
     weight = 1.0/(y_err**2.) 
  
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

             if(len(x[bin_inds]) > 0):
                 if(yerr !=None):
                     bin_weight = np.sum(weight[bin_inds])
                 else:
                     bin_weight = 1.
                 ybin_cur = np.sum(y[bin_inds]*weight[bin_inds])/bin_weight
                 ybin.append(ybin_cur)
                 ybinerr.append(np.sqrt(1.0/bin_weight))
                 if(weigh_x==True):
                     xbin.append(np.sum(x[bin_inds]*weight[bin_inds])/bin_weight)
             else:
                 ybin.append(0.)
                 ybinerr.append(0.)
          
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
             bin_weight = np.sum(weight[bin_inds])
             ybin.append(np.sum(y[bin_inds]*weight[bin_inds])/bin_weight)
             ybinerr.append(np.sqrt(1.0/bin_weight))
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
     
     if(yerr==None):
         bin_dict = {'xbin':xbin, 'ybin':ybin}
         return bin_dict
     else:
         bin_dict = {'xbin':xbin, 'ybin':ybin, 'ybinerr':ybinerr}
         return bin_dict
#         return xbin, ybin, ybinerr

def testfunc(x):
    f = x**2.0 - 9.0
    return f

# Given two neighbouring par files and an epoch in MJD (e.g. glitch epoch), 
# calculate the delta_nu (nuder=0) or delta_nudot (nuder=1) at the glitch epoch.
def fit_delta_nu_OLD(par1, par2, mjd, mjd_error=0, nuder=0):
    pars1 = smu.read_parfile(par1)
    pars2 = smu.read_parfile(par2)
    midt1 = pars1['TZRMJD'].value 
    midt2 = pars2['TZRMJD'].value 
#    ferr1 = np.ones_like(mjd)*pars1['F' + str(nuder)].error
#    ferr2 = np.ones_like(mjd)*pars2['F' + str(nuder)].error
#    for i in range(nuder+1, 12):
#        fdot_name = 'F' + str(i)
#        if fdot_name in pars1.keys() and  pars1[fdot_name].value != 0.0:
            #ferr1 = np.sqrt( ferr1**2 + ((pars1[fdot_name].error)*((mjd-midt1)*3600.*24)**(i-nuder))**2) 
#    ferr1 = np.ones_like(mjd)*pars1['F' + str(nuder)].error
#    ferr2 = np.ones_like(mjd)*pars2['F' + str(nuder)].error
#    ferr1 = pars1['F' + str(nuder)].error
#    ferr2 = pars2['F' + str(nuder)].error
    ferrs1 = np.zeros(12) #, len(mjd)))
    ferrs2 = np.zeros(12) #, len(mjd)))
    ferr_t1 = 0.0
    ferr_t2 = 0.0
    for i in range(nuder, 12):
        fdot_name = 'F' + str(i)
        if fdot_name in pars1.keys() and pars1[fdot_name].value != 0.0:
            ferrs1[i] = (pars1[fdot_name].error)*(((mjd-midt1)*3600.*24)**(i-nuder))
            if(i>nuder):
                ferr_t1 = ferr_t1 + (i-nuder)*pars1[fdot_name].value*(((mjd-midt1)*3600.*24)**(i-nuder-1))
        if fdot_name in pars2.keys() and  pars2[fdot_name].value != 0.0:
            ferrs2[i] = (pars2[fdot_name].error)*(((mjd-midt2)*3600.*24)**(i-nuder))
            if(i>nuder):
                ferr_t2 = ferr_t2 + (i-nuder)*pars2[fdot_name].value*(((mjd-midt2)*3600.*24)**(i-nuder-1))
        # include error on mjd if given
    ferr1 = np.sqrt((ferrs1**2).sum(0) + ((mjd_error*3600.*24*ferr_t1)**2)) 
    ferr2 = np.sqrt((ferrs2**2).sum(0) + ((mjd_error*3600.*24*ferr_t2)**2)) 
#    print 'ferrs1, ferrs2 = ', (ferrs1**2).sum(0) , (ferrs2**2).sum(0)
#    print 'ferr_t1, feerr_t2 = ', ((mjd_error*3600.*24*ferr_t1)**2), ((mjd_error*3600.*24*ferr_t2)**2)
#    print 'ferr1, ferr2 = ', ferr1, ferr2
##        if fdot_name in pars2.keys() and  pars2[fdot_name].value != 0.0:
##            ferr2 = np.sqrt( ferr2**2 + ((pars2[fdot_name].error)*((mjd-midt2)*3600.*24)**(i-nuder))**2) 
            
    f1 = smu.times2freqs(mjd, par1, nuder)
    f2 = smu.times2freqs(mjd, par2, nuder)
    
    
    
    deltaf = f2-f1
    deltaf_err = np.sqrt(ferr1**2.0 + ferr2**2.0)
    
    deltaf_over_f = deltaf/f1
    
    return deltaf, deltaf_err, deltaf_over_f




def calc_nu_err_at_epoch(par, mjd, mjd_off=0., nuder = 0):
    
    pars = smu.read_parfile(par)
    # t = np.linspace(pars['START].value, pars['FINISH'].value, 100)
    midt = pars['TZRMJD'].value 
    # nuderr = 0
    freqs = np.zeros(0)
    ferrs = np.zeros(0)
#    freqs = 0.
#    ferrs = 0.
    for i in range(nuder, 12):
        fdot_name = 'F' + str(i)
        if fdot_name in pars.keys() and  pars[fdot_name].value != 0.0:
            freqs = np.append(freqs, pars[fdot_name].value)
            ferrs = np.append(ferrs, pars[fdot_name].error)
    fs = smu.calc_freq(mjd, pars['PEPOCH'].value, *freqs)
    f_err = smu.calc_freq(mjd_off,  0, *ferrs)
#    f_err = smu.calc_freq(mjd, midt, *ferrs)
#    f_err +=  smu.calc_freq(mjd + mjd_off, midt, *freqs)-smu.calc_freq(mjd, midt, *freqs) 
    return f_err
#    plt.plot(t, fs)
#    plt.fill_between(t, fs-f_errs, fs+f_errs, alpha=0.5, color='k')
    
def fit_delta_nu(par1, par2, mjd, mjd_error=0, nuder=0):
    pars1 = smu.read_parfile(par1)
    pars2 = smu.read_parfile(par2)
    midt1 = pars1['TZRMJD'].value 
    midt2 = pars2['TZRMJD'].value 
            
    f1 = smu.times2freqs(mjd, par1, nuder)
    f2 = smu.times2freqs(mjd, par2, nuder)
    
    ferr1 = calc_nu_err_at_epoch(par1, mjd, mjd_off=mjd_error, nuder=nuder)
    ferr2 = calc_nu_err_at_epoch(par2, mjd, mjd_off=mjd_error, nuder=nuder)
    
    
    deltaf = f2-f1
    deltaf_err = np.sqrt(ferr1**2.0 + ferr2**2.0)
    
    deltaf_over_f = deltaf/f1
    delta_f_over_f_err = deltaf_over_f*np.sqrt((deltaf_err/deltaf)**2 + (ferr1/f1)**2)
    
    return deltaf, deltaf_err, deltaf_over_f, delta_f_over_f_err
