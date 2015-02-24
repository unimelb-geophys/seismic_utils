# -*- coding: utf-8 -*-

#=============================================================================
                            ###START PROGRAM###
#=============================================================================

"""
Edited on Thu Feb 24 13:21:33 2015
@author: AGOS UniMelb Team incl: Louis Moresi, Dan Sandiford, Abe Jones. 
@editor: Benjamin Boland
"""

#=============================================================================
                            ###START IMPORTS###
#=============================================================================
from __future__ import with_statement
from __future__ import division

from datetime import date
from math import *
from math import radians, cos, sin, asin, sqrt
from obspy.core import UTCDateTime, Stream, Trace, read, AttribDict
from obspy.core.util import NamedTemporaryFile
from obspy.mseed import util
from obspy.mseed import * 
from obspy.mseed.core import readMSEED, writeMSEED, isMSEED
from obspy.mseed.headers import clibmseed, ENCODINGS
from obspy.mseed.msstruct import _MSStruct
from obspy import read
from obspy import signal
from os import open
from scipy import signal

import calendar, copy, csv, datetime, glob, io, math, obspy, os, shutil, \
StringIO, time, unittest, urllib2, warnings

import matplotlib.pyplot as plt, numpy as np

#=============================================================================
                            ###END IMPORTS###
#=============================================================================

#=============================================================================
                            ###START FUNCTIONS###
#=============================================================================


def split_channels(Y, julian, network_name, name, date_station):
    """
    Function to import mini-seed files from working directory, split them based
    on the specified channels, then save the files with a new function name. 
    """
    
    channels = ['Z'] 
    
    for channel in channels:
    
        #read mseed files
        sts = read("%s.mseed" %(date_station))
        #select channel to split
        
        channel_split = sts.select(component="%s" %(channel))
        
        
        #check to see if there is already a directory with that name, if not, 
        #create one
        if not os.path.exists("networks/%s/%s" %(network_name, name)):\
        os.makedirs("networks/%s/%s" %(network_name, name))
    
        #save individual channels to correct directories and correct file name
        #for BUD structure for MSNOISE
        channel_split.write('networks/%s/%s/%s.%s..EN%s.%s.%s' %(network_name, name,\
        name, network_name, channel, Y, julian), format='MSEED')


    
def read_channels(date_station):
    
    sts = read("%s.mseed" %(date_station))
    print(sts)
    #sts.plot()



def spectrum(directory):
    """
    produces Power Density Spectrum for given cross-correlation waveforms. 
    Input takes directory where waveforms are being stored. 
    Outputs are the power density spectrum/ frequency spectrum plots

    """

    for data_folder in os.listdir(directory):
        plt.figure()
        for j in range(1,8):
        
            j = str(j).zfill(2)

            for i in range(1,11):
            
                i = str(i).zfill(2)
            
                try:
                    sts = read('/home/boland/Documents/MSNoise-master/msnoise/STACKS/%s/\
001_DAYS/ZZ/%s/2014-01-%s.MSEED'%(j, data_folder, i))
                    tr = sts[0]
                    wave = tr.data #this is how to extract a data array from a mseed file
                    fs = tr.stats.sampling_rate
                
                    #hour = str(hour).zfill(2) #create correct format for eqstring
                    f, Pxx_spec = signal.welch(wave, fs, 'flattop', 1024, scaling='spectrum')
                    #plt.semilogy(f, np.sqrt(Pxx_spec))
                    plt.title("Frequency Density Plot of: %s" %(data_folder))
                    plt.plot(f, np.sqrt(Pxx_spec))
                    plt.xlim([0, 1])
                    plt.ylim([0, 0.01])
            
                    plt.xlabel('frequency [Hz]')
                    plt.ylabel('Linear spectrum [V RMS]')
                    plt.show()

                except(IOError):
                    continue

def stations():
    
    """
    The following function returns a list of latitudes and longitudes, as 
    well as returning seismic station names. All of this is from a .csv file
    stored in the directory variable below. No inputs. Remember to filter values!
    """
    
    directory = "Stations/stations.csv"
    
    stations = io.open(directory, "rb") 
    stations = csv.reader(stations)#import in csv file format
    stations.next() #get rid of titles
    stations = [row for row in stations] #create array of data
    stations1 = [row[0] for row in stations] #take station name data only
    latitude = [row[1] for row in stations]
    longitude = [row[2] for row in stations]
    
    #stations, latitude, longitude = stations()
    #latitude = filter(None, latitude) #remove all null strings from list
    #longitude = filter(None, longitude) #remove all null strings from list
    #station_names = filter(None, stations) #remove all null strings from list  
    #latitude = [float(i) for i in latitude]   
    #longitude = [float(i) for i in longitude]  
    
    return stations1, latitude, longitude
    
 


def distance(latitude, longitude, station_names):
    """ 
    Function to return the distance (in kms) between two lat, long points
    using the Haversine Formula. 
    Input is a list of latitudes and longitudes, as well as the station names
    associated with them. Output is a list of station correlations names, 
    distances between those stations, and a concatenated list of both of these
    together. 
    """
    
    dlat = [0]*len(station_names)**2   #initialise difference in latitudes for
    dlong = [0]*len(station_names)**2
    dist = [0]*len(station_names)**2
    data_folder_string = [0]*len(station_names)**2
    
    counts = 0
    counts1 = 0
    
    for name in station_names:
        
        for game in range(0, len(station_names)):
            
            #name of two stations taking the difference between
            #gives data folder name for stacks folders format!
            data_folder_string[counts] = "_%s__%s"%(name, station_names[game]) 
            
            lat1 = radians(latitude[game])
            lat2 = radians(latitude[counts1])
            long1 = radians(longitude[game])
            long2 = radians(longitude[counts1])
            
            dlat[counts] = abs(lat1 - lat2) #gives lat difference for all combinations
            dlong[counts] = abs(long1 - long2)
            
            dlat1 = dlat[counts]; dlong1 = dlong[counts];
            
            a = sin(dlat1/2)**2 + cos(lat1) \
            * cos(lat2) * sin(dlong1/2)**2
            
            c = 2 * asin(sqrt(a))    
                        
            dist[counts] = 6367 * c
            
            counts += 1
            
        counts1 += 1 #used for keeping counts of number of first set of stations used. 
    
    #create .csv file of data_folder_string and dist (distances)
    #np.column stack, stacks two lists into columns as required for .csv file
    #stat_dist = np.column_stack((data_folder_string, dist))
    
    return data_folder_string, dist
            
    
def leap_year (year):
    """
    Simple function that returns 1 if the year entered is a leap year, 
    returns 0 if not. 
    """
    
    if (year % 100 == 0 ): # Gregorian fix
        if (year % 400 == 0 ):
            return (1)
        else:
            return (0)
    else:
        if (year % 4 == 0 ):
            return (1)
        else:
            return (0)
      


def calendar_date(year, doy):
    """
    Function that returns a string of the year, the month and the day of the
    month respectively. Input is the year and the day of the year or julian
    day. 
    Example use would be: 
    print(year, month, day = calendar_date(2014, 29))
    
    year = '2014'
    month = '1'
    day = '29'
    
    """
    
    if doy < 32: month = 1; day = doy
    
    elif doy < 60 + leap_year(year): month = 2; day = doy - 31
    
    else:
        if leap_year(year) == 0:
            doy += 1
            month = int((doy+31.39)/30.61)
            day = doy + 2 - (month-1)*30-int((month+1)*0.61)
    
    return str(year), str(month), str(day)


def merge_delete(Y,M,D,station_name, date_station):
    
    #merge individual station's hourly data into a whole day's 

    for wave in glob.glob("%s*.mseed" %(date_station)):
                                  #inputs should be of the same form as the
                                  #saved hour duration outputs from AGOS 
                                  #network
        os.system("cat "+wave+" >> %s.mseed" %(date_station)) 
                                             #output should have station name
                                             #and date stamp only. 
        
    #delete individual station's hourly data
    for s in glob.glob("%sT*.mseed" %(date_station)):
          os.remove(s)
          
          
   
def return_index(main_list, search_str): 
    """
    The following function returns the index of a given search string
    from a specified searched list. This can either be a partial search 
    or a full search string
    """    
    
    indices = [i for i, s in enumerate(main_list) if search_str in s]
    #indices = indices[0] #make indices an int NOT an array/list
    indices = indices[0]
    
    return indices
    

#=============================================================================
                            ###END FUNCTIONS###
#=============================================================================



#=============================================================================
                            ###START VARIABLES###
#=============================================================================

#set important variables
project_duration = 2 #set duration of the project measured in days

Y = "2014"                #start year
M = "01"                  #start month
D = "09"                  #start day of that month







#=============================================================================
                            ###CALL FOR START-DATE DURATION###
#=============================================================================

check = False  
while check == False:
    start_date = raw_input("\nEnter start date of your project: ")
    
    if "-" not in start_date:
        check = False
        print("\nError, input incorrect. Date format MUST BE dd-mm-yy\n")
        print("For example: 01-03-2011")
        #note WILL NOT work if someone's input DOES have - but is the wrong format
        
    else:
        start_date_list = start_date.split("-")
        year = start_date_list[2] #returns year as a string
        month = start_date_list[1]#returns month as a string
        day = start_date_list[0] #returns day of month as a string
    
        if int(day) not in range(1,32) or int(month) not in range(1,13):
            check = False
            print("\nError, input incorrect. Date format MUST BE dd-mm-yy.\n")
            print("For example 01-03-2011.\n")
            print("Day must be between 1 and 31.\n")
            print("Month must be be between 1 and 12.")
        
        elif int(year) not in range(1700, (date.today().year + 1)):
            check = False
            print("\nError, input incorrect. Date format MUST BE dd-mm-yy.\n")
            print("For example 01-03-2011.\n")
            print("Year must be between 1700 and %d" %(date.today().year + 1))       
        else:
            check = True
            print("\nYour inputs are:  year = %s, month = %s, day = %s." \
            %(year, month, day))
            
del check

#=============================================================================
                            ###CALL FOR OUTPUT DURATION###
#=============================================================================
months = ['Jan', 'Feb', 'Mar','Apr','May','Jun', 'Jul', 'Aug', 'Sep', 'Oct', \
                                                                  'Nov', 'Dec']
months_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

time_dur_types = ['minute', 'hour', 'day']
dur_secs = [60, 3600, 86400]

check = False  
while check == False:
    time_dur_string = raw_input("\nEnter output duration you require: ")
                   
    if time_dur_string not in time_dur_types:
        check = False
        print("\nError, input incorrect. Options are: minute, hour, day, month.")
        
    else:
        check = True
        print("\nYou have entered %s for your time duration output files." \
        %(time_dur_string))
        
        dur_secs_index = return_index(time_dur_types, time_dur_string)
        
        time_dur = dur_secs[dur_secs_index]
 
del check


#=============================================================================
                            ###CALL FOR DATA STRUCTURE###
#=============================================================================


#data_structure = {}
#data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
#data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
#data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
#data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

struct_types = ['Singles', 'BUD', 'SDS', 'IDDS', 'PDF', 'ANT']


check = False  
while check == False:
    data_struct = raw_input("\nEnter the data file structure you require: ")
                   
    if data_struct not in struct_types:
        check = False
        string_struct = ''
        for a in struct_types: string_struct += "'%s' " %(a) 
        print("\nError, input incorrect. Options are: %s" %(string_struct))
        
    else:
        check = True
        print("\nYou have entered %s as the data structure for your output files." \
        %(data_struct))



#INPUTS SO FAR: data_struct, time_dur and year, month, day
        


            
        


#import array of stringed names for each of the different stations (in melbourne)


#stations = io.open("Stations/stations.csv", "rb") 
#stations = csv.reader(stations)#import in csv file format
#stations.next() #get rid of titles
#stations = [row for row in stations] #create array of data
#stations = [row[0] for row in stations] #take station name data only
#secs_in_hour = 3600 #86400       #number of seconds in a day


#the following function gives the julian date (day of the year) given what the
   #current gregorian date is. 
#initialise julian date i.e. day of the year since january!
julian = datetime.datetime(int(Y), int(M), int(D), 00, 00, 00)\
.timetuple().tm_yday

#=============================================================================
                            ###END VARIABLES###
#=============================================================================




#=============================================================================
                            ###START FUNCTION CALLS###
#=============================================================================

#stations, latitude, longitude = stations()

for day in range(1, project_duration):
    #calendar_date function returns date from day of the year
    date = calendar_date(int(Y), julian);

    Y = date[0];
    M = date[1];
    #month = str(M).zfill(2) #create correct format for eqstring
    D = date[2];
    #day = str(D).zfill(3) #create correct format for eqstring
    date_stamp = "%s-%s-%-s" %(Y, M, D);
    dur = "%d" %(secs_in_hour)#string telling the program the duration 
                                  #of data you require
    mdur = (float(dur))/60
    mdur = math.ceil(mdur)
    mdur = int(mdur)
    
        
    for name in stations:
        
        date_station = "%s_%s" %(name, date_stamp);
        
        for hour in range(0, 24): 
            
            #set start date
            H = "%d" %(hour)          #start hour
            m = "0"                  #start minute
            
            stamp = "%s-%s-%sT%s:00" %(Y, M, D, str(H).zfill(2))
            
    #files are meant to be one hour long and seem to start at minute 00
    #For broken files there will need to be some additions to the script 
            #hour = str(hour).zfill(2) #create correct format for eqstring
            #downloading miniseed file for each individual station. 
            
            
            eqstring = 'http://agos1.kelunji.net/eqserver/eqwaveextractor?year=%s&\
month=%s&day=%s&hour=%s&minute=%s&duration=%s&servernum=0&conttrig=0&\
sitechoice=list&sitelist=+%s+&siteradius=&closesite=&radius=&latitude=&\
longitude=&fileformat=miniseed&getwave=Get+Waveform' \
%(Y, M, D, hour, m, \
str(mdur), name)
              
        
            fin_string = 'wget --user=eq --password=event55s "%s" -O %s_%s.mseed'\
%(eqstring, name, stamp)

            os.system(fin_string)
    
    #merge files into one one day file, then delete the hourlies.
    #by putting the function here, it's deleting the hourlies as soon as each
    #day file is CREATED!
        

        merge_delete(Y, M, D, name, date_station)
       
        try:
            split_channels(Y, julian, network_name, name, date_station)
        
        except(TypeError):
            continue
        
    #delete all multiplexed data, and all files that have 0 data
    #the remaining data with 0 byte size will have returned typeerror above
        for s in glob.glob("*.mseed"):
            os.remove(s)     
        
    if julian < 365:
        julian += 1; #go to the next day's data, unless it's the end of the year!
    else:
        julian = 1; Y = str(int(Y) + 1);
    
    #call the MSNOISE script for each day that passes in order to correctly process cross-correlations
    #os.system('/home/boland/Documents/MSNoise-master/msnoise/scripts/cron.sh') 


#=============================================================================
                            ###END FUNCTION CALLS###
#=============================================================================


#=============================================================================
                            ###START PROGRAM###
#=============================================================================