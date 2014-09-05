#!/usr/bin/python

import time
import math
from datetime import date, datetime, timedelta
import numpy as np

def mjdtodate(mjdin, dateformat=None):

    timestamp = (mjdin - 40587.)*86400.
#    mjdin = time.time()/86400. + 40587.
    
    if (type(timestamp) is list or type(timestamp) is np.ndarray):
        date_out = np.array([datetime.fromtimestamp(timestamp[i_time]) for i_time in range(len(timestamp))])
    else:
        date_out = datetime.fromtimestamp(timestamp)

    if (dateformat != None):
##    print dateout.strftime("%Y-%b-%d")  # year-month-date
        date_out = date_out.strftime(dateformat)  # year-month-date
       
    
    return date_out


# Returns current MJD
def mjd():

    mjd_out = time.time()/86400. + 40587.
 
    return mjd_out


# Returns MJD at (floating point) year given)
def datetomjd(year, month, day, hour, minutes, seconds):

    date_in = \
    datetime(year, month, day, hour, minutes, seconds)
    
    date_in = date_in.utctimetuple()

#    date_in.tm_isdst=0
    timestamp_in = time.mktime(date_in)

    mjd_out = (timestamp_in/86400.) + 40587.

    return mjd_out

# Input here is a decimal year, returns mjd after figuring out datetime
# object value 
def yeartomjd(year_in):
    year = int(math.floor(year_in))
    yday = (year_in - float(year))*365.


# use yday-1, since days are counted in datetime object starting from 1
    date_in = datetime(year, 1, 1) + timedelta(days=yday-1)

# Now use datetime objects to input for datetomjd:
    mjd_out = datetomjd(date_in.year, date_in.month, date_in.day, \
                        date_in.hour, date_in.minute, date_in.second)
    
    return mjd_out


def mjdtoyear(mjdin):

    dateout = mjdtodate(mjdin).utctimetuple()
    
    yearout = dateout.tm_year + dateout.tm_yday/365. + dateout.tm_hour/24./365. + \
        dateout.tm_min/60./24./365. + dateout.tm_sec/60./60./24./365.

    return yearout
