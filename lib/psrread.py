# Scripts for reading in various pulsar data  products and packaging them into dictionary objects

from sys import argv, exit
import os
import numpy as np
import matplotlib.pyplot as plt


duty_file='/Users/ferdman/Work/pulsar/data_tables/duty_lookup_table.dat'
duty_default = 0.55
tempo2_resid_format_file='/Users/ferdman/Work/pulsar/data_tables/tempo2_resid_format.dat'

def get_duty(psr_name):
     
     # Open duty cycle data file
     try:
          f_duty = open(duty_file,'r')
     except IOError as (errno, strerror):
          if (errno == 2):  # file not found
               print "IOError ({0}): File".format(errno), duty_file,"not found"
          else:
               print "IOError ({0}): {1}".format(errno, strerror)
               return
          
          print "File", duty_file, "open."
          
     # read in file
     lines = f_duty.readlines()
     f_duty.close()

     # Cycle through lines in file to find psr_name
     duty = -1
     for line in lines:
        cur_line = line.split()
        if(psr_name in cur_line[0]):
            duty = float(cur_line[1])
            print "Duty fraction for pulsar " + psr_name + ": ", duty
            break

     if(duty < 0):
         duty = duty_default
         print "Could not find pulsar " + psr_name + " in duty lookup table."
         print "Consider adding a new line to " + duty_file
         print "Will set duty to " + duty_default + " for no good reason."
         print ""

     return duty



# Routine to read print_resids-format ascii residual data file.  returns a 
# dictionary that contains mjd, res, res_err, and freq
def read_resid(resid_file, info_file=None, tempo2=False):     

    # Open residual file
    try:
        f_res = open(resid_file, 'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), resid_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
            return
        
    print "File", resid_file, "open."
        
# Test to see whether this residual file is from tempo2
    if (tempo2):
# Read columns in file into numpy arrays (change this to have more info, 
# like tempo1 output
##         mjd, freq, res, reserr = \
##             np.loadtxt(f_res, dtype='double', comments='#', 
##                        usecols=(0,1,2,3), unpack=True)
###         resp, res, reserr, preres, orbphase, mjd, obsfreq, nobs = \
###             np.loadtxt(f_res, dtype='double', comments='#', 
###                        usecols=(1,2,3,4,5,6,7,8), unpack=True)
# Read in lines of file, and cast each column into appropriate data types
         resid_arr = [line.split() for line in f_res.readlines()]
# Will assume that there are 8 columns in file as commented out above:
         n_rows = len(resid_arr)
         n_col = len(resid_arr[0])
         if(n_col != 9):
              print 'For tempo2 output there should be 9 columns.  Exiting...'
              return
         resp=np.array([resid_arr[i][1] for i in range(n_rows)], dtype=float)
         res=np.array([resid_arr[i][2] for i in range(n_rows)], dtype=float)
         reserr=np.array([resid_arr[i][3] for i in range(n_rows)], dtype=float)
         preres=np.array([resid_arr[i][4] for i in range(n_rows)], dtype=float)
         orbphase=np.array([resid_arr[i][5] for i in range(n_rows)], dtype=float)
         mjd=np.array([resid_arr[i][6] for i in range(n_rows)], dtype=float)
         obsfreq=np.array([resid_arr[i][7] for i in range(n_rows)], dtype=float)
         nobs=np.array([resid_arr[i][8] for i in range(n_rows)])
         
# Just make serial number to be in order of residuals printed by tempo2
         serial = [i for i in range(len(res))] 
# Default units for tempo2 is seconds -- make it microseconds here:
         res = res*1e6
         preres = preres*1e6
# For now, don't use info numbers
         info_num = None
         info_val = None
         info_instr = None
# Instead of info numbers, here, if user has used the '--info' flag, we will
# plot different nobs different colours
         info_val, info_indices = np.unique(nobs, return_index=True)
         print 'INFO_VAL = ', info_val
         info_num = nobs # just to make compatible with rest of code and plotter
       
# Otherwise, assume it is in c.tmp format from tempo1:    
    else:
         resp, res, reserr, preres, orbphase, mjd, serial, obsfreq = \
             np.loadtxt(f_res, dtype='double', comments='#', 
                        usecols=(1,2,3,4,5,6,8,10), unpack=True)

# resp:     postfit residual in periods
# res:      postfit residual in microseconds
# reserr:   residual error in microseconds
# preres:   prefit residual in microseconds
# orbphase: orbital phase
# mjd:      MJD
# serial:   Residual serial number (integer)
# obsfreq:  Observing frequency
         # if an info file is given, read it, assigning each first column 
         # to second column (and in plotting, to a colour)
         # Otherwise, fill it with NULLs:
         if(info_file != None and os.path.isfile(info_file)):
              try:
                   f_info = open(info_file, 'r')
              except IOError as (errno, strerror):
                   if (errno == 2):  # file not found
                        print "IOError ({0}): File".format(errno), \
                            info_file, "not found"
                   else:
                        print "IOError ({0}): {1}".format(errno, strerror)
                        return
        
              print "Info file", info_file, "open."
         
         # Assuming info.tmp format:
              info_arr = [line.split() for line in f_info.readlines()]
              info_num = np.array([info_arr[i][0] for i in range(len(info_arr))])
              # Test to see if there are two columns by checking number of 
              # elements in each line.  If there is one column, assign 
              # info_instr=NULL
              instrument=[]
              for info_line in info_arr:
                   if (len(info_line) == 2):
                        instrument.append(info_line[1])
                   else:
                        instrument.append(NULL)
              instrument=np.array(instrument)
         
         # Array for unique values of info number, and corresponding
         # second column, if it exists
              info_val, info_indices = np.unique(info_num, return_index=True)
              info_instr = instrument[info_indices]
                   
              print 'Found ', len(info_val), ' info numbers: '
              for i_info in np.arange(len(info_val)):
                   print '   ', info_val[i_info], '  ', info_instr[i_info]
         else:
              info_num = None
              info_val = None
              info_instr = None

# Just conert serial into an integer list:
    serial = [int(x) for x in serial]


# Create data packsg in dictionary form
    resid_data = {'resp':resp, 'res':res, 'reserr':reserr, 'preres':preres, \
                       'orbphase':orbphase, 'mjd':mjd, 'serial':serial, \
                       'obsfreq':obsfreq, \
                       'info':info_num, 'info_val':info_val, \
                       'info_instr':info_instr}

    
    f_res.close()

    

    return resid_data





# Routine that reads in asp-styles ascii profile
def read_asc_prof(prof_file):
    
    # Open ascii profile file
    try:
        f_prof = open(prof_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), prof_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
            return
        
    print "File", prof_file, "open."
      
    bin, I, Q, U, V = \
        np.loadtxt(f_prof, dtype='double', comments='#', usecols=(0,1,2,3,4), 
                   unpack=True)
      
    phase = bin/(len(bin)-1)

    prof_data = {'phase':phase, 'i':I, 'q':Q, 'u':U, 'v':V}

    f_prof.close()
    return prof_data



# Routine to read header of asp-style ascii file
def read_asc_header(prof_file):

# Open ascii profile file
    try:
        f_prof = open(prof_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), prof_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
            return
        

# Read header line
    line = f_prof.readline()
    line = line.split()

# Test to make sure that first character in header is '#'
    if(line.count('#') != 1 or line.index('#') != 0):
        print "Header line is not correct format.  Exiting."
        exit()

# Remove '#' character from header line list
    del line[0]

# Set up dictionary for each quantity of header line:

    header = {'imjd':float(line[0]), 'smjd':float(line[1]), \
              'period':float(line[2]), 'obsfreq':float(line[4]), \
              'dm':float(line[5]), 'nbins':int(line[6]), \
              'obscode':line[7], 'psrname':line[9], 'phase':float(line[10])}
        
    f_prof.close()
    return header


# Routine to read  ascii profile width data file.  Returns a 
# dictionary that contains mjd, and width, width error in units of 
# profile phase [0,1].
def read_widths(width_file, rfile=None):     

    # Open residual file
    try:
        f_width = open(width_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), width_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
            exit
        
    print "File", width_file, "open."
        
    mjd, width, werr = \
        np.loadtxt(f_width, dtype='double', comments='#', 
                   usecols=(0, 1, 2), unpack=True)        

# Create data packsg in dictionary form
    width_data = {'mjd':mjd, 'width':width, 'werr':werr}

# If we want to adjust widths by comparing means of a reference width file and 
# subtracting by difference of means, then user must provide a reference width 
# file (same format), and this adjustment will be made.
    if(rfile != None):
        rw = read_widths(rfile)
        width_data['width'] = width_data['width'] - \
            (np.average(width_data['width'], \
                             weights=1./(width_data['werr']**2)) - \
                  np.average(rw['width'], weights=1./(rw['werr']**2)))
        width_data['werr'] = np.sqrt(width_data['werr']**2 + \
                                (np.std(width_data['width']))**2 + \
                                (np.std(rw['width']))**2 )

  
    f_width.close()

    return width_data




def read_par(par_file,file_format='tempo1'):
    
# define an exception that will be used if the arguments are bad:
    class ArgError(Exception):
        """Exception called when something is wrong with one of the arguments
        passed to a function
    
        Attributes:
            name    -- argument name
            value   -- argument value given by user
            message -- explanation of the error
        """
        def __init__(self, name, value):
            self.name = name
            self.value = value
            self.output = "Argument \'"+value+"\' to \'"+name+ \
                "\' keyword is not valid."
    
# Check that input file_format is either tempo1 or tempo2
# can be given as:
#    "tempo1", "tempo2", "t1", "t2", "1", or "2", and case-insensitive
    try:
        # make sure it's al string, then convert to lower case
        file_format = str(file_format).lower()
        # list of acceptable file_format keyword arguments
        good_format = ['tempo', 'tempo1', 'tempo2', 't1', 't2', '1', '2']
        # determine if user-inputted argument is in list of acceptable arguments
        is_good_format = good_format.count(file_format)
        if (not is_good_format):
            raise ArgError('file_format',file_format)
    except ArgError as err:
        # print "name = ", err.name
        # print "value = ", err.value
        print err.output
        return 

# now make file_format either 'tempo1' or 'tempo2'
    t1_format = ['tempo', 'tempo1', 't1', '1']
    t2_format = ['tempo2', 't2', '2']
    if (t1_format.count(file_format)):
        file_format = 'tempo1'
    else:
        file_format = 'tempo2'

    print "Parameter file format: ", file_format 
    
# Open par file
    try:
        f_par = open(par_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), par_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
        return

    print "File", par_file, "open.  Format is "+file_format+"."

# Read file, putting each line as a separate string array element
    lines = f_par.readlines()

# Length of file
    f_len = len(lines)

# Here is a (somewhat exhaustive) list of possible par file parameters:
    

# Parse each line
#    param_name = ['psr', 'raj', 'decj', 'pmra', 'pmdec', 'px', 'f', \
    param_name = ['psr', 'raj', 'decj', 'pmra', 'pmdec', 'px', 'f', \
                       'pepoch', 'start', 'finish', 'dm', 'ephem', 'clk', \
                       'ntoa', 'tres', 'tzrmjd', 'tzrfrq', 'tzrsite', \
                       'nits', 'binary', 'pb', 'ecc', 'om', 't0', 'a1', \
                       'tasc', 'eps1', 'eps2', 'omdot', 'pbdot', 'gamma', \
                       'm2', 'sini', 'shapmax']
    param_type = ['s',   's',   's',    'f',    'f',     'f',  'f', \
                       'f',      'f',     'f',      'f',  's',    's', \
                       'd',    'f',    'f',      'f',      's',  \
                       'd',    's',      'f',  'f',   'f',  'f',  'f', \
                       'f',    'f',    'f',    'f',     'f',     'f', \
                       'f',   'f',   'f']
    f_val = ['0.0' for i in range(16)] # max 15 f derivatives
    dm_val = ['0.0' for i in range(8)] # max 7 dm derivatives
##    temp_params = []
# make a list as long as the param list and fill with empty strings:
    param_val = ['' for ip in range(len(param_name))]
# Write first two columns to separate temporary variables
    for cur_line in lines:
        field = cur_line.split()
# If first character is "C" then ignore
        if (field[0] == 'C'):
            pass
            # print "Line is commented out:\n"
            # print "  ", cur_line
        else:
# convert first column to lowercase
            param = field[0].lower()
            value = field[1]
            ## print "param = ", param, ", value = ", value
# change 'e' to 'ecc' if this is a tempo1 par file:
            if (file_format == 'tempo1' and param == 'e'): param = 'ecc' 
# change 'psrj' to 'psr' if file is in tempo2 format -- still not sure about 
# storing psr j vs b names...
            if (file_format == 'tempo2' and param == 'psrj'): param = 'psr' 
##            temp_params.append( (param, value) )
# If param if rot. frequency or a derivative...
            if (param[0] == 'f' and len(param) < 3):
# ...determine derivative...
                deriv = int(param[1:])
# ...and finally place value in appropriate element of f_val array:
                if (deriv >= 0 and deriv < len(f_val)):
                    f_val[deriv] = value
# If param if DM or derivative...
            elif (param[0:2] == 'dm' and len(param) < 5):
                if(len(param) == 2): # i.e. actual DM
# ...place actual DM in first element of dm_val array:
                    dm_val[0] = value
                else: # i.e. one of the derivatives
                    deriv = int(param[2:])
# ... place in appropriate place in dm_val array if derivative:
                    if (deriv > 0 and deriv < len(dm_val)):
                        dm_val[deriv] = value
                    else:
                        print "Error: DM derivative in par file is either", \
                            "zero or larger than the maximum of", len(dm_val)
                        return
# If we find the parameter that corresponds to the official list of names...
            elif (param_name.count(param)):
# ...then enter it into corresponding element in value array:
                param_val[param_name.index(param)] = value
            else:
                print "Unrecognized parameter "+param+".  Continuing for now..."

# Now insert f_val and dm_val lists into appropriate place in param_val array
    ## print "f_val =", f_val
    ## print "dm_val =", dm_val
    f_ind = param_name.index('f')
    del param_name[f_ind]
    del param_val[f_ind]
    for i_deriv in range(len(f_val)):
        f_name='f{0}'.format(i_deriv)
        param_name.insert(f_ind+i_deriv, f_name)
        param_val.insert(f_ind+i_deriv, f_val[i_deriv])

    dm_ind = param_name.index('dm')
    del param_name[dm_ind]
    del param_val[dm_ind]
    for i_deriv in range(len(dm_val)):
        if (i_deriv == 0):
            dm_name='dm'
        else:
            dm_name='dm{0}'.format(i_deriv)
        param_name.insert(dm_ind+i_deriv, dm_name)
        param_val.insert(dm_ind+i_deriv, dm_val[i_deriv])

#    print temp_params
#    params = dict(temp_params)
    params = dict(zip(param_name, param_val))
##    print "List of parameters in file "+par_file+":\n"
##    for p, v in params.iteritems():
##        print p, v

# Do comparison to expected fields in par file.  Map read-in parameter naems to 
# param names as we want them to be stored (If only have jname or bname, not
# sure how to get other), and make separate mappings for tempo1 and tempo2 
# formats
    

#    t1_params = ['psr', 'raj', 'decj', 'pmra', 'pmdec', 'px', 'f', \
#                     'pepoch', 'start', 'finish', 'dm', 'ephem', 'clk', \
#                     'ntoa', 'tres', 'tzrmjd', 'tzrfreq', 'tzrsite', 'nits' \
#                     'binary', 'pb', 'e', 'omega', 't0', 'a1', 'tasc', \
#                     'eps1', 'eps2', 'omdot', 'pbdot', 'gamma', 'm2', 'sini']


        

# Close par file
    f_par.close()

    print "File", par_file, "closed"
    
    return params
