# -*- coding: utf-8 -*-
                            ###START PROGRAM###

"""
Edited on Thu Feb 12 13:21:33 2015
@author: AGOS UniMelb Team incl: Louis Moresi, Dan Sandiford etc. 
@editor: Benjamin Boland
"""
                            ###START MODULES###
from __future__ import with_statement
import math
import os
from os import open
import obspy
import urllib2, StringIO, csv, numpy as np,matplotlib.pyplot as plt,time
import io, glob, shutil, time
from obspy import read
import datetime, calendar
from funcs import *

                            ###END MODULES###


                            ###START VARIABLES###

#set important variables
project_duration = 2 #set duration of the project measured in days

Y = "2014"                #start year
M = "01"                  #start month
D = "09"                  #start day of that month

secs_in_hour = 3600 #86400       #number of seconds in a day

network_name = "AGOS" #please enter the name of your network of seismometers

## Function to download individual waveforms for correct MSNOISE data structure
    



# raw input:

#Y = raw_input("enter year, e.g. 2010: ")
#M = raw_input("enter Month, e.g. 02: ")
#D = raw_input("enter Day, e.g. 07: ")
#H = raw_input("enter hour (24 hr), e.g. 13: ")
#m = raw_input("enter minute, e.g. 00: ")
#dur = raw_input("enter duration in seconds, e.g. 180: ")


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
                            ###END VARIABLES###

                            ###START FUNCTION CALLS###
stations, latitude, longitude = stations()

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


                            ###END FUNCTION CALL###



                            ###END PROGRAM###