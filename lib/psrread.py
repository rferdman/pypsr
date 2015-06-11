# Scripts for reading in various pulsar data  products and packaging them into dictionary objects

from sys import argv, exit
import os
import numpy as np
import matplotlib.pyplot as plt
from psrcalc import rad_to_cos_phi


duty_file='/Users/ferdman/Work/pulsar/data_tables/duty_lookup_table.dat'
duty_default = 0.55
tempo2_resid_format_file='/Users/ferdman/Work/pulsar/data_tables/tempo2_resid_format.dat'

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

def get_duty(psr_name):
     
     # Strip off leading 'J' or 'B' if there
     if(psr_name[0] == 'J' or psr_name[0] == 'B'):
          psr_name = psr_name[1:]

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
         print "Will set duty to ", duty_default, " for no good reason."
         print ""

     return duty



# Routine to read print_resids-format ascii residual data file.  returns a 
# dictionary that contains mjd, res, res_err, and freq
def read_resid(resid_file, tempo2=False, info_file=None, info_flag=None):     

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
         n_tempo2_cols = 9 # Standard, non-flag number of columns in tempo2 residual file
         file_lines = f_res.readlines()
         resid_arr = [line.split() for line in file_lines[0:len(file_lines)-2]]
         # Get rid of lines like "JUMP", "SKIP", "NO SKIP", etc.

# Will assume that there are 8 columns in file as commented out above:
         n_rows = len(resid_arr)
         print 'NROWS = ', n_rows
         min_n_col = np.amin(np.array([len(resid_arr[i]) for i in np.arange(n_rows)]))
         max_n_col = np.amax(np.array([len(resid_arr[i]) for i in np.arange(n_rows)]))
#         n_col = len(resid_arr[0])
         if(min_n_col < n_tempo2_cols):
              print 'MIN_N_COL = ', min_n_col
              print 'For tempo2 output there should be at least 9 columns.  Exiting...'
              min_ind = np.argmin(np.array([len(resid_arr[i]) for i in np.arange(n_rows)]))
              print 'line', min_ind
              print resid_arr[min_ind]
              return
         resp=np.array([resid_arr[i][1] for i in range(n_rows)], dtype=float)
         res=np.array([resid_arr[i][2] for i in range(n_rows)], dtype=float)
         reserr=np.array([resid_arr[i][3] for i in range(n_rows)], dtype=float)
         preres=np.array([resid_arr[i][4] for i in range(n_rows)], dtype=float)
         orbphase=np.array([resid_arr[i][5] for i in range(n_rows)], dtype=float)
         mjd=np.array([resid_arr[i][6] for i in range(n_rows)], dtype=float)
         obsfreq=np.array([resid_arr[i][7] for i in range(n_rows)], dtype=float)
         nobs=np.array([resid_arr[i][8] for i in range(n_rows)])
         # If there are more than 9 columns, the rest are flags from the TOA line
         # if(max_n_col > n_tempo2_cols):
              # Will handle each row separately since different rows may have different 
              # numbers of flags
         #flag_exists = []
         flag_d = [] # turn this and args into np.array() later
         #flag_arg = []
         for i_row in range(n_rows):
              n_col = len(resid_arr[i_row])
              flag_d_temp = {}
              if(n_col > n_tempo2_cols):
                   flag_data=resid_arr[i_row][9:n_col]
                   i_flag=0
                   while(i_flag < len(flag_data)):
                        if(flag_data[i_flag].startswith('-') and 
                           not flag_data[i_flag+1].startswith('-')):
                             flag_d_temp[flag_data[i_flag]] = flag_data[i_flag+1]
                             i_flag += 2
                        elif(flag_data[i_flag].startswith('-') and 
                             flag_data[i_flag+1].startswith('-') or
                             flag_data[i_flag].startswith('-') and 
                             i_flag==flag_data-1):
                             flag_d_temp[flag_data[i_flag]] = ''
                             i_flag += 1
                        else:
                             print 'Problem with flag section of TOA line: '
                             print '    '+' '.join(resid_arr[i_row])
                             
                             
                   #flag_exists_temp = True
                   ##flag_data=np.array(resid_arr[i_row][9:n_col])
                   ##n_flag_data = len(flag_data)
                   ##is_flag = [flag_data[i][0]=='-' for i in range(len(flag_data))]
                   ##flag_inds = np.where(is_flag)[0]
                   ##flag_name_temp = flag_data[flag_inds]
                   ##flag_arg_temp = np.empty_like(flag_name_temp)
                   ##for i_flag in range(len(flag_inds)):
                   ##     next_ind = flag_inds[i_flag]+1 
                   ##     if next_ind < n_flag_data:
                   ##          if (is_flag[next_ind] == False):
                   ##               flag_arg_temp[i_flag] = flag_data[next_ind]
                   ##          else:
                   ##               flag_arg_temp[i_flag] = ''
                   ##     else:
                   ##          flag_arg_temp[i_flag] = ''
              else:
                   #flag_exists_temp = False
                   flag_d_temp['no_flag'] = ''
                   ##flag_name_temp = ['']
                   ##flag_arg_temp = ['']
 
              flag_d.append(flag_d_temp)
              
              
              #flag_exists.append(flag_exists_temp) 
              #flag_arg.append(flag_arg_temp)

         ##flag_name = np.array(flag_name)
         ##flag_arg = np.array(flag_arg)
#         print 'flag_name = ', flag_name
#         print 'flag_arg  = ', flag_arg
     
     # Now we can select which values to use for info for colour coding when plotting.
     # Choose backend by default, unless there is no backend flag, then revert to telescope`

     #flag_data=np.array([resid_arr[i][9:n_col] for i in range(n_rows)])
              #flag_name=
         
# Just make serial number to be in order of residuals printed by tempo2
         serial = np.arange(len(res)) #[i for i in range(len(res))] 
# Default units for tempo2 is seconds -- make it microseconds here:
         res = res*1e6
         preres = preres*1e6
# For now, don't use info numbers
         info_id = None
         info_val = None
         info_instr = None
# If user has used the --infoflag option, then use given flag to colour-code
# residuals.
         if (info_flag != None):
# Check the flag names and if the requested flag is not present for all TOAs, 
# then give error and exit
              flag_exists = np.sum(np.array([flag_d[i_flag].has_key(info_flag) 
                                             for i_flag in range(len(flag_d))])) == len(flag_d)
              ##test_flags = []               
              ##for i_flag in range(np.shape(flag_name)[0]):
              ##     test_flags.append(np.sum(np.array([flag_name[i_flag][j_flag]==info_flag 
              ##                           for j_flag in range(len(flag_name[i_flag]))])) )
              ##flag_exists = np.sum(np.array(test_flags)) == np.shape(flag_name)[0]
              print 'FLAG_EXISTS = ',flag_exists
              # If the flag exists in every TOA line, then proceed with assigning values
              # to info variables.
              ##print 'FLAG_NAME = ', flag_name
              ##print 'FLAG_ARG = ', flag_arg
              if(flag_exists):
                   info_id = np.array([flag_d[i_flag][info_flag] for i_flag in range(len(flag_d))])
                   ##info_id = np.extract(flag_name==info_flag, flag_arg)
                   info_val, info_indices = np.unique(info_id, return_index=True)
              # Otherwise, give error message and quit. 
              else:
                   print 'Requested flag ' +info_flag+ ' for colour-coding does not exist'
                   print ' for some TOAs.  Exiting.'
                   exit()
                   

         else:
# Instead of info numbers, here, if user has used the '--info' flag, we will
# plot different nobs different colours by deafult
              info_val, info_indices = np.unique(nobs, return_index=True)
              info_id = nobs # just to make compatible with rest of code and plotter
         # print 'INFO_ID = ', info_id
         # print 'INFO_VAL = ', info_val
      
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
              info_id = np.array([info_arr[i][0] for i in range(len(info_arr))])
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
              info_val, info_indices = np.unique(info_id, return_index=True)
              info_instr = instrument[info_indices]
                   
              print 'Found ', len(info_val), ' info numbers: '
              for i_info in np.arange(len(info_val)):
                   print '   ', info_val[i_info], '  ', info_instr[i_info]
         else:
              info_id = None
              info_val = None
              info_instr = None

# Just conert serial into an integer list:
    serial = [int(x) for x in serial]


# Create data packsg in dictionary form
    resid_data = {'resp':resp, 'res':res, 'reserr':reserr, 'preres':preres, \
                       'orbphase':orbphase, 'mjd':mjd, 'serial':serial, \
                       'obsfreq':obsfreq, \
                       'info':info_id, 'info_val':info_val, \
                       'info_instr':info_instr}
    
    f_res.close()

    

    return resid_data





# Routine that reads in asp-styles ascii profile
def read_asc_prof(prof_file, ionly=False):
    
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
    
    if(ionly): # Only read Stokes I  
        bin, I = \
            np.loadtxt(f_prof, dtype='double', comments='#', usecols=(0,1), 
                   unpack=True)
        # Set Q, U, V to be zeroes in this case
        Q=np.zeros_like(I)
        U=np.zeros_like(I)
        V=np.zeros_like(I)
    else:
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
def read_widths(width_file, rfile=None, units_in='deg', units_out='deg'):     
            
    print "File", width_file, "open."
    
    possible_units=['radians', 'radian', 'rad', 'degrees', 'degree', \
        'deg', 'phase', 'cos_phi']
    if(possible_units.count(units_in)==0):
        print 'ERROR: Input width units must be either \'rad\', \'deg\', \'phase\', or \'cos_phi\'.'
        exit()
    if(possible_units.count(units_out)==0):
        print 'ERROR: Input width units must be either \'rad\', \'deg\', \'phase\', or \'cos_phi\'.'
        exit()
                
    if(units_in=='cos_phi'):
        print 'ERROR: \'cos_phi\' is not currently supported as an input width unit.'
        exit()            
 
    # Open width file
    try:
        f_width = open(width_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print 'IOError ({0}): File {1}'.format(errno), width_file, ' not found.'
        else:
            print 'IOError ({0}): {1}'.format(errno, strerror)
            exit()
        
    # Read header line
    header_line = f_width.readline().split()
    if(header_line[0] == '#'): # Indicates the existence of a header line
        header_line_exists = True
        psr_name = header_line[1] # Zeroth char is a '#'
        percent_height=header_line[2]
        phase_range = np.array([header_line[3], header_line[4]])
    else:
        header_line_exists = False
        print 'WARNING: Input width file does not have a header line.'
    # In any case, rewind file pointer position for loadtxt (not sure if needed)
    f_width.seek(0)
    
    # Read input file
    mjd, width, werr = \
        np.loadtxt(f_width, dtype='double', comments='#', 
        usecols=(0, 1, 2), unpack=True)        

    # Now sort out input/output unit conversions
    if(units_in=='phase'):
        if(units_out[0:3]=='rad'):
            width *= 2.0*np.pi
            werr  *= 2.0*np.pi
        elif(units_out[0:3]=='deg'):
            width *= 360.0
            werr  *= 360.0
        elif(units_out=='cos_phi'):
            phi = 0.5 * width * 2.0*np.pi 
            phi_err = 0.5 * werr * 2.0*np.pi
            width, werr = rad_to_cos_phi(phi, phi_err)
    elif(units_in[0:3]=='rad'):
        if(units_out=='phase'):
            width /= 2.0*np.pi
            werr  /= 2.0*np.pi
        elif(units_out[0:3]=='deg'):
            width *= 180.0/np.pi
            werr  *= 180.0/np.pi
        elif(units_out=='cos_phi'):
            phi = 0.5 * width 
            phi_err = 0.5 * werr
            width, werr = rad_to_cos_phi(phi, phi_err)
    elif(units_in[0:3]=='deg'):
        if(units_out=='phase'):
            width /= 360.0
            werr  /= 360.0
        elif(units_out[0:3]=='rad'):
            width *= np.pi/180.0
            werr  *= np.pi/180.0
        elif(units_out=='cos_phi'):
            phi = 0.5 * width * np.pi/180.0
            phi_err = 0.5 * werr * np.pi/180.0
            width, werr = rad_to_cos_phi(phi, phi_err)

# Create data packsg in dictionary form
    width_data = {'mjd':mjd, 'width':width, 'werr':werr}
    if(header_line_exists):
        width_data['psr'] = psr_name
        width_data['percent'] = percent_height
        width_data['phase'] = phase_range

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

    
    
def read_par(par_file, file_format='tempo1', return_tuple=False, verbose=False):
#    def __init__(self, par_file, file_format='tempo1'):

    maxf=16 # maximum number of frequency derivatives
    maxdm=8 # maximum number of dm derivatives

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

    if (verbose==True):
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

    if (verbose==True):
        print "File", par_file, "open.  Format is "+file_format+"."

# Read file, putting each line as a separate string array element
    lines = f_par.readlines()

# Length of file
    f_len = len(lines)

# Here is a (somewhat exhaustive) list of possible par file parameters:


# Parse each line
#    param_name = ['psr', 'raj', 'decj', 'pmra', 'pmdec', 'px', 'f', \
    param_name = ['psr', 'psrj', 'raj', 'decj', 'pmra', 'pmdec', 'px', 'f', \
                       'pepoch', 'start', 'finish', 'dm', 'ephem', 'clk', \
                       'ntoa', 'tres', 'tzrmjd', 'tzrfrq', 'tzrsite', \
                       'nits', 'binary', 'pb', 'ecc', 'om', 't0', 'a1', \
                       'tasc', 'eps1', 'eps2', 'omdot', 'pbdot', 'gamma', \
                       'm2', 'sini', 'shapmax', 'posepoch', 'dmepoch', 'ephver', \
                       'glep', 'glph', 'glf0', 'glf1', 'units', 'mode']
    param_type = ['s',   's',    's',   's',    'f',    'f',     'f',  'f', \
                       'f',      'f',     'f',      'f',  's',    's', \
                       'd',    'f',    'f',      'f',      's',  \
                       'd',    's',      'f',  'f',   'f',  'f',  'f', \
                       'f',    'f',    'f',    'f',     'f',     'f', \
                       'f',   'f',   'f',        'f',         'f',       'd', \
                       'f',    'f',    'f',    'f',    's',     'd']
    f_val = [0.0 for i in range(maxf)] # max 15 f derivatives
    f_err = [-1.0 for i in range(maxf)] # max 15 f derivatives
    f_fit = [False for i in range(maxf)] # max 15 f derivatives
    dm_val = [0.0 for i in range(maxdm)] # max 7 dm derivatives
    dm_err = [-1.0 for i in range(maxdm)] # max 7 dm derivatives
    dm_fit = [False for i in range(maxdm)] # max 7 dm derivatives
##    temp_params = []
# make a list as long as the param list and fill with empty strings:
    param_val = ['' for ip in range(len(param_name))]
    param_fit = [False for ip in range(len(param_name))]
    param_err = ['' for ip in range(len(param_name))]
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
            fit = False
            error = '-1'
            ################### NEED TO DISTINGUISH BETWEEN TEMPO1 ANND 2,AND ALSO GET UNCERTAINTIES!!! #######
            if(len(field) > 2):
                if(field[2] == '1'):
                    fit = True
                if(file_format == 'tempo1'):
                    if(len(field) == 4):
                        error = field[3]
                else:  # tempo2 format
                    if(len(field) == 3 and fit==False):
                        error = field[2]
                    elif(len(field) == 4):
                        error = field[3]
                        #print param, ' ', value, ' ', error, ' ', fit
#                else:
#                    fit = False
#            else:
#                fit=False
            ## print "param = ", param, ", value = ", value
# change 'e' to 'ecc' if this is a tempo1 par file:
            if (file_format == 'tempo1' and param == 'e'): param = 'ecc' 
# change 'psrj' to 'psr' if file is in tempo2 format -- still not sure about 
# storing psr j vs b names...
            if (file_format == 'tempo2' and param == 'psrj'): param = 'psr'

##            temp_params.append( (param, value) )
# If param if rot. frequency or a derivative...
            if (param[0].startswith('f') and len(param) < 4):  # < 4 to avoid the "finish" parameter
# ...determine derivative...
                deriv = int(param[1:])
# ...and finally place value in appropriate element of f_val array:
                if (deriv >= 0 and deriv < len(f_val)):
                    f_val[deriv] = float(value)
                    f_err[deriv] = float(error)
                    f_fit[deriv] = fit
                    #print param, ' ', f_val[deriv], ' ', f_err[deriv], ' ', f_fit[deriv]
# If param if DM or derivative...
            elif (param.startswith('dm') and len(param) <= 5):
                if(len(param) == 2): # i.e. actual DM
# ...place actual DM in first element of dm_val array:
                    dm_val[0] = float(value)
                    dm_err[0] = float(error)
                    dm_fit[0] = fit
                else: # i.e. one of the derivatives
                    if(param[2] != 'x'):
                        deriv = int(param[2:])
    # ... place in appropriate place in dm_val array if derivative:
                        if (deriv > 0 and deriv < len(dm_val)):
                            dm_val[deriv] = float(value)
                            dm_err[deriv] = float(error)
                            dm_fit[deriv] = fit
                        else:
                            print "Error: DM derivative in par file is either", \
                                "zero or larger than the maximum of", len(dm_val)
                            return

# HERE AT SOME POINT INCLUDE GLITCH PARAMETERS, IN AN ARRAY FORMAT AS WITH F AND DM

# If we find the parameter that corresponds to the official list of names...
            elif (param_name.count(param)):
# ...then enter it into corresponding element in value array:
#setattr(self, param, )
                param_ind = param_name.index(param)
                if (param_type[param_ind] == 'f'):
                    value = float(value)
                elif (param_type[param_ind] == 'd'):
                    value = int(value)
                param_val[param_ind] = value
                param_fit[param_ind] = fit
                param_err[param_ind] = float(error)
            else:
                print "Unrecognized parameter "+param+".  Continuing for now..."

# Now insert f_val and dm_val lists into appropriate place in param_val array, keeping 
    ## print "f_val =", f_val
    ## print "dm_val =", dm_val
    f_ind = param_name.index('f')
    param_val[f_ind] = np.array(f_val)
    param_fit[f_ind] = np.array(f_fit)
    param_err[f_ind] = np.array(f_err)
    
#    del param_name[f_ind]
#    del param_val[f_ind]
#    for i_deriv in range(len(f_val)):
#        f_name='f{0}'.format(i_deriv)
#        param_name.insert(f_ind+i_deriv, f_name)
#        param_val.insert(f_ind+i_deriv, f_val[i_deriv])

    dm_ind = param_name.index('dm')
    param_val[dm_ind] = np.array(dm_val)
    param_fit[dm_ind] = np.array(dm_fit)
    param_err[dm_ind] = np.array(dm_err)
#    del param_name[dm_ind]
#    del param_val[dm_ind]
#    for i_deriv in range(len(dm_val)):
#        if (i_deriv == 0):
#            dm_name='dm'
#        else:
#            dm_name='dm{0}'.format(i_deriv)
#        param_name.insert(dm_ind+i_deriv, dm_name)
#        param_val.insert(dm_ind+i_deriv, dm_val[i_deriv])



#    print temp_params
#    params = dict(temp_params)
    params = dict(zip(param_name, param_val))
#    params = dict(zip(param_name, param_val, param_err, param_fit))
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

    if (verbose==True):
        print "File", par_file, "closed"

#    return params
    if(return_tuple):
        return param_name, param_val, param_fit, param_err
    else:
        #print 'Par file:'
        #print params
        return params



#class f0(real):
#    def __init__(self, f0):
#        super(self).__init__(f0)
        
def read_dmx(dmx_file):     

    # Open residual file
    try:
        f_dmx = open(dmx_file,'r')
    except IOError as (errno, strerror):
        if (errno == 2):  # file not found
            print "IOError ({0}): File".format(errno), dmx_file, "not found"
        else:
            print "IOError ({0}): {1}".format(errno, strerror)
            exit
        
    print "File", dmx_file, "open."
        
    mjd, delta_dm, delta_dm_err = \
        np.loadtxt(f_dmx, dtype='double', comments='#', 
                   usecols=(0, 1, 2), unpack=True)        

# Create data packsg in dictionary form
    dmx_data = {'mjd':mjd, 'delta_dm':delta_dm, 'delta_dm_err':delta_dm_err}

    return dmx_data
