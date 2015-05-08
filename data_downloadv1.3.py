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

import sys
sys.path.append("/home/abe/anaconda/lib/python2.7/site-packages")

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
from obspy import read, read_inventory
from obspy import station
from obspy.station import stationxml

from obspy import read, read_inventory
from obspy import station
from obspy.station import stationxml


import calendar, copy, csv, datetime, glob, io, math, obspy, os, shutil, \
StringIO, time, unittest, urllib2, warnings

import matplotlib.pyplot as plt, numpy as np

#=============================================================================
                            ###END IMPORTS###
#=============================================================================

#=============================================================================
                            ###START FUNCTIONS###
#=============================================================================

#data_structure = {}
#data_structure['SDS'] = "YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
#data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
#data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
#data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

#def Singles_struct():

def merge_months():
    month_file = 'MONTH_TIMOR'
    folders_list = os.listdir(month_file)
    stations = []
    for folder in folders_list:
        for files in os.listdir('%s/%s'%(month_file, folder)):
            stations.append(files.split(".")[0])
        
        stations = list(np.unique(stations))
        print(stations)
        for stat in stations:
            for wave in glob.glob("%s/%s/%s*" %(month_file, folder, stat)):
                os.system("cat "+wave+" >> %s/%s/UM.%s.EHZ.mseed"\
                    %(month_file, folder, stat))
                    
                os.system('rm -r %s'%(wave))

        
    

def move_month(read_dir):
    file_name = 'tmp'
    path = '%s/%s' %(file_name, read_dir)
    month_file = 'MONTH_TIMOR'
    info = path.split(".")
    stat_name = info[0]
    date = info[1]
    print(date)
    year = date.split('-')[0]
    month = date.split('-')[1].zfill(2)
    day = date.split('-')[2]
    year_month = '%s-%s' %(year, month)
    print(year_month)
    
    if not os.path.exists(month_file):\
        os.makedirs(month_file)
    
    
    month_folder = '%s/%s' %(month_file, year_month)
    
    if not os.path.exists(month_folder):\
        os.makedirs(month_folder)
    
    orig_path = path
    new_path = '%s/%s' %(month_folder, read_dir)
    shutil.move(orig_path, new_path)



def split_downsample(read_dir):
    sample = 20.0 #new downsizing sample rate!
    try:	
        st = read(read_dir)
        #select only the Z channel
        st = st.select(component="Z")
        
        print(st)
	
        for trace in st[:]:
            if trace.stats.channel != "EHZ":
                st.remove(trace)
            else:
                continue
         	   #print(trace.stats.channel)
            trace.stats.sampling_rate = sample


        st.merge()
        tr = st[0]
        tr.stats.sampling_rate = sample
        st.write(read_dir, format = "MSEED")
        
        print(st)


    except Exception as e:
        print("\nOops, there was a problem\n")
        print(e)
        
def move(read_dir):
    file_name = 'tmp'
    if not os.path.exists(file_name):\
        os.makedirs(file_name)
        
    orig_path = read_dir
    new_path = '%s/%s' %(file_name, orig_path)
    
    shutil.move(orig_path, new_path)
    


def metadata(read_dir):
    

#currently reading an example miniseed file in as a stream
    try:
        st = read(read_dir)
        print(st)
        x = len(st)



    #CORRECTING FOR NETWORK CODE AND LOCATION

        network_change1 = 'UM'
        network_new_name1 = 'UM'#'BM'
        network_change2 = '01'
        network_new_name2 = 'UM'
        network_change3 = 'A'
        network_change4 = ''
        location_blank = ''
    

        for i in range(0, x):
            tr = st[i]
    
    # removes LOCATION so it is blank, as listed in the metadata files (regardless of what it was previously)
   
            tr.stats["location"] = location_blank
    
    # Changes BOREHOLE network codes from UM to BM and SURFACE network codes from 01 to UM 
            net = tr.stats["network"]
            if network_change1 in net:
                tr.stats["network"] = network_new_name1
            elif network_change2 or network_change3 or network_change4 in net:
                tr.stats["network"] = network_new_name2
            else:
                continue

    #CORRECTING BOREHOLE STATION NAMES

        serial_no_1 = 'A346'
        site_name_1 = 'LOYU'
        serial_no_2 = 'BD5E'
        site_name_2 = 'MOSU'
        serial_no_3 = 'BD70'
        site_name_3 = 'SGWU'
        serial_no_4 = 'BD91'
        site_name_4 = 'WILU'

    #Changes station name from serial number to station code

        for i in range(0, x):
            tr = st[i]
            stat = tr.stats["station"] 
            if serial_no_1 in stat:
                tr.stats["station"] = site_name_1
            elif serial_no_2 in stat:
                tr.stats["station"] = site_name_2
            elif serial_no_3 in stat:
                tr.stats["station"] = site_name_3
            elif serial_no_4 in stat:
                tr.stats["station"] = site_name_4
            else:
                continue
        
        # CHANGES TO CHANNEL CODE

    # (this is a bit messy at the moment since the wildcard feature seemed to be failing)
        
        channel_new_name_E = 'EHE'
        channel_new_name_N = 'EHN'
        channel_new_name_Z = 'EHZ'
        
        for i in range(0, x):
            tr = st[i]
            chan = tr.stats["channel"] 

    # Changes CHANNEL names from '**E', '**N', '**Z', (e.g. BHE, DHZ) to a consitant format of EHE, EHN, EHZ
    # EXCEPT FOR BOREHOLE STATIONS, which will maintain channel codes BHE, BHN, BHZ

            if 'DHE' in chan:
                tr.stats["channel"] = channel_new_name_E
            elif 'DHN' in chan:
                tr.stats["channel"] = channel_new_name_N
            elif 'DHZ' in chan:
                tr.stats["channel"] = channel_new_name_Z    
            elif 'ENE' in chan:
                tr.stats["channel"] = channel_new_name_E
            elif 'ENN' in chan:
                tr.stats["channel"] = channel_new_name_N
            elif 'ENZ' in chan:
                tr.stats["channel"] = channel_new_name_Z         
            else:
                continue        
        
# saves stream as a combination of edited traces       
            st[i] = tr

# ATTACH METADATA TO STREAM

#import and read metadata file from GitHub
        os.system('wget https://raw.githubusercontent.com/unimelb-geophys/metadata/master/UOM.xml')


        metadata = stationxml.read_StationXML('UOM.xml')
    #print(metadata)

    #st.write("")

    #attach metadata to stream
        st.attach_response(metadata)
        print("\nResponse attached successfully\n")
    #st.remove_response()

        print(st)
    
    #delete downloaded metadata file 
        os.system('rm -r UOM.xml')


        st.write(read_dir, format="MSEED")
    
    #check saved response file has only one trace and is only one channel!
        st = read(read_dir)
        print(st)
    
    except Exception as e:
        print("\nOops, there was a problem with attaching the metadata\n")
        print(e)
    
def structure(file_name, data_struct):
    file_name1 = file_name
    file_name = file_name.split(".")
    
    ###IF YOU WISH TO CHANGE YOUR NETWORK NAME, PLEASE CHANGE THE CODE BELOW!###
    network_name = "UM"

    channel = file_name[0]
    station_name = file_name[1]
    dates = file_name[2].split("-") 
    Y = dates[0]
    M = dates[1]
    D = dates[2]
    julian = datetime.datetime(int(Y), int(M), int(D), 00, 00, 00)\
    .timetuple().tm_yday

    if data_struct == "BUD":\
        
        end_file_name = "%s.%s..%s.%s.%s" %(station_name, network_name, \
        channel, year, julian)
        
        destination = "networks/%s/%s/%s" %(network_name, station_name,\
        end_file_name)
        
        if not os.path.exists("networks/%s/%s" %(network_name, station_name) ):\
        os.makedirs("networks/%s/%s" %(network_name, station_name))
        
        #shutil.move(destination, file_name)
        os.rename(file_name1, destination)

    #elif data_struct == "SDS":
    #    end_file_name = "%s.%s..%s.%s.%s" %(station_name, network_name, channel, \
    #    year, julian)
    #    pathway = "networks/%s/%s/%s" %(network_name, station_name, end_file_name)
    #    
    #elif data_struct == "IDDS":
    #    end_file_name = "%s.%s..%s.%s.%s" %(station_name, network_name, channel, \
    #    year, julian)
    #    pathway = "networks/%s/%s/%s" %(network_name, station_name, end_file_name)

    #elif data_struct == "PDF":
    #    end_file_name = "%s.%s..%s.%s.%s" %(station_name, network_name, channel, \
    #    year, julian)
    #    pathway = "networks/%s/%s" %(network_name, station_name) 

    else:
        a = 5
    
    
#def SDS_struct(file_name):
#    continue
    
#def IDDS_struct(file_name):
#    continue
    
#def PDF_struct(file_name):
#    continue
    
#def ANT_struct(file_name):
#    continue
    



#def data_structure(data_struct):
        
    #if data_struct == 'Singles':
    #    continue
    
    
    #elif data_struct == 'BUD'
    
    #elif data_struct == 'SDS'
    
    #elif data_struct == 'ANT'
    
    #elif data_struct == 'IDDS'

    #elif data_strcut == 'PDF'    

    

    
def starttime(year, month, day, hour, mins):
    
    start_time = "%s-%s-%s-%s:%s" %(year, month.zfill(2), day, \
                                    hour.zfill(2), mins.zfill(2))
    return start_time
    
    
def download(year, month, day, hour, mins, time_dur, station_name):
    
    """
    This function returns a downloaded file from our eq server.
    The parameters in order to get said file are: 
    1. duration - this is in seconds and set in the inputs part of the program. 
    2. start date - in format dd-mm-yy, the input takes Y as year, M as month
        and D as day of the month from this input. 
    3. 
    """
    
    # the durataiton of data you require
    dur = time_dur
    time_dur = (float(time_dur))/60
    time_dur = math.ceil(time_dur)
    time_dur = str(int(time_dur))
    
    #set start date
    start_time = starttime(year, month, day, hour, mins)
    
    
    year  =  int(year) #returns year as a string
    month =  int(month)#returns month as a string
    day   =  int(day) #returns day of month as a string
    hour  =  int(hour)#returns hour in the day as a string
    mins  =  int(mins)#returns minute in the hour
    
    call_eq_string = 'http://agos2.kelunji.net/eqserver/eqwaveextractor?year=%s&\
month=%s&day=%s&hour=%s&minute=%s&duration=%s&servernum=0&conttrig=0&\
sitechoice=list&sitelist=+%s+&siteradius=&closesite=&radius=&latitude=&\
longitude=&fileformat=miniseed&getwave=Get+Waveform' \
%(str(year), str(month), str(day), str(hour), str(mins), time_dur, station_name)
    
              
    output_file_name = "%s.%s.%s.mseed" %(station_name, start_time, str(dur))
    print(output_file_name)
    final_string = 'wget "%s" -O %s'\
    %(call_eq_string, output_file_name)

    os.system(final_string)
    
    return(output_file_name)

    #--user=eq --password=event55s

def split_channels(time_dur, station_name, output_file_name, data_struct):
    
    """
    Function to import mini-seed files from working directory, split them based
    on the specified channels, then save the files with a new function name. 
    """
    
    try:
        sts = read(output_file_name)
        sts = sts.merge()
    
    
    
    
    #channels = ['E','N','Z']         
        channels = ['Z']         
        network_name = "UM"
    
        for channel in channels:
               
                #select channel to split
        
            channel_split = sts.select(component="%s" %(channel))
        
        #save individual channels to correct directories and correct file name
        #for BUD structure for MSNOISE
        
            file_name = 'EN%s.%s' %(channel, output_file_name)
        
            try:
                channel_split.write(file_name, format='MSEED')
                structure(file_name, data_struct)
            
            
            ###if the data is 0 byte size, then just move file of correct name!
            except(NotImplementedError, TypeError):
                continue
            #if not os.path.exists("networks/%s/%s" %(network_name, station_name) ):\
            #os.makedirs("networks/%s/%s" %(network_name, station_name))
    except(NotImplementedError, TypeError):
        a=5
            #shutil.move(destination, file_name)
            #os.rename(file_name1, destination)
            
        ###IF THE DATA STRUCTURE IS DIFFERENT TO THE SINGLES OPTIONS, THEN 
        #CHANGE THE FILE STRUCTURE!

        
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

def stations2():
    
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
            
    
def leap_year(year):
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
        else:
            doy += 1
            month = int((doy+31.39)/30.61)
            day = doy + 2 - (month-1)*30-int((month+1)*0.61)
    
    return str(year), str(month), str(day)


def merge_delete(year, month, day, hour, mins, time_dur, station_name):
    
    #merge individual station's hourly data into a whole day's 
    
    #file_name = "%s.*.mseed" %(station_name)
    file_name = "%s.%s-%s-%s-*.mseed" %(station_name, year, month.zfill(2), day)

    day_file_name = "%s.%s-%s-%s.mseed" % (station_name, year, month.zfill(2), day)
    
    for wave in glob.glob(file_name):
                                  #inputs should be of the same form as the
                                  #saved hour duration outputs from AGOS 
                                  #network
    
        os.system("cat "+wave+" >> %s" %(day_file_name)) 
                                             #output should have station name
                                             #and date stamp only. 
        
    #delete individual station's hourly data
    for s in glob.glob(file_name):
          os.remove(s)
    
    return day_file_name
          
          
   
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

project_duration = 1194 #set duration of the project measured in days


#=============================================================================
                            ###ENTER DEFAULTS###
#=============================================================================
#list of all the station names stored on ours servers
all_stations = ['LOCU','CLIF','HODL','HOLS','KRAN','BRAD','MELU','MRDN',\
                'NARR','PEGM','PEGU','S88U','SOMU','WALM','ADE','DNL','MRAT']

#start date example: 2012-03-24-00-12
start_date = '2012-01-01-00-00'
start_date_list = start_date.split("-")
year = start_date_list[0] #returns year as a string
month = start_date_list[1]#returns month as a string
day = start_date_list[2] #returns day of month as a string
hour = start_date_list[3]
mins = start_date_list[4]

#example: 3600 seconds for an hour, or 60 seconds for a minute etc. 
#DO NOT INPUT minute, hour or day as a string for this input!
#set to 0 for raw inputs
time_dur = 86400

#default output file structure is singles. This means all files will just be
#save to the current working directory. 
data_struct = 'BUD'

#default station name is 'CLIF' for testing. You can either one string or
#a list of strings containing each station name. 
#find out how to get a list of station names from the XML file as 

#set stations to 0 for raw input
stations = all_stations#["HOLS"]##["HOLS", "NARR", "MRDN", "LOCU"] #['HOLS']  # # ['CLIF'] or ['CLIF', 'HODL'] or set as 0 for raw input!

split = 'Yes'


#=============================================================================
                            ###CALL FOR START-DATE###
#=============================================================================

check = False  
while check == False:
    
    if start_date == 0:
        
        start_date = raw_input("\nEnter start date of your project: ")
    else:
        break
    
    if "-" not in start_date:
        check = False
        print("\nError, input incorrect. Date format MUST BE yy-MM-dd-hh-mm\n")
        print("For example: 2012-03-24-00-12")
        #note WILL NOT work if someone's input DOES have - but is the wrong format
        
    else:
        start_date_list = start_date.split("-")
        year = start_date_list[0] #returns year as a string
        month = start_date_list[1]#returns month as a string
        day = start_date_list[2] #returns day of month as a string
        hour = start_date_list[3]
        mins = start_date_list[4]
    
        if int(day) not in range(1,32) or int(month) not in range(1,13) \
        or int(hour) and int(mins) not in range (0,60):
            check = False
            print("\nError, input incorrect. Date format MUST BE yy-MM-dd-hh-mm\n")
            print("For example: 2012-03-24-00-12")
            print("Day must be between 1 and 31.\n")
            print("Month must be be between 1 and 12.")
        
        elif int(year) not in range(1700, (date.today().year + 1)):
            check = False
            print("\nError, input incorrect. Date format MUST BE dd-mm-yy.\n")
            print("For example 01-03-2011.\n")
            print("Year must be between 1700 and %d" %(date.today().year + 1))       
        else:
            check = True
            print("\nYour inputs are:  year = %s, month = %s, day = %s \
hour = %s, minute = %s" \
            %(year, month, day, hour, mins))
            
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
    #break the while loop if defaults aren't 0s
    if time_dur == 0:
        time_dur_string = raw_input("\nEnter output duration you require: ")
    else:
        break
                  
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


struct_types = ['Singles', 'BUD', 'SDS', 'IDDS', 'PDF', 'ANT']


if data_struct == 0:
    
    check = False  
    while check == False:
        #if 'Singles' in data_struct:
        #    break
    
        #else:
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
            
else:
    data_struct = data_struct



#INPUTS SO FAR: data_struct, time_dur and year, month, day

#check if the duration is a day long, if so, make sure to download hourlies first
#then concatenate. 

#also get it to initiate a download if more than one hour duration, but less than a day. 


#=============================================================================
                            ###CALL FOR STATION NAMES###
#=============================================================================
stations_list = []
check = False  
check2 = False

###CODE STILL NEEDS TO CHECK IF THE STATION ENTERED IS BOTH IN THE LIST OF
###STATIONS FOR OUR NETWORK AND THE LIST BEING CREATED! 
if stations != 0:
    
    for i in stations:
        if i in all_stations:
            stations_list = stations
        else:
            print("\nThe station '%s' you entered is not in our network." %(i))
            print("\nIf you would like to add a new station to the list of\
 stations in this network, please change the all_stations list")
            quit()
        
    
else: 
    while check == False:
        stations = raw_input("\nPlease enter the station name you require: ")
        
        if stations in all_stations:
            if stations not in stations_list:
                stations_list.append(stations)
                print("\nYour current station/s to be outputed is/are: ")
                print(stations_list)
            else:
                print("\nPlease choose a station different from those entered already ")
                continue
        
            check = False
            check2 = False
            while check2 == False:
                
                YN = raw_input("\nDo you also require data from another station? (Y/N) ")
            
                
                if YN == 'Y': # or 'yes' or 'YES' or 'Yes' or 'y':
                    check = False
                    check2 = True
                    
                elif YN == 'N': # or 'no' or 'NO' or 'No' or 'n': 
                    check = True
                    check2 = True
                    
                else:
                    print("\nPlease enter either Y or N as your answer")
                    check = False
                    check2 = False
        else:
            check = False
            print("\nThe station you entered is not in our network.")
            print("\nIf you would like to add a new station to the list of\
 stations in this network, please change the all_stations list")
 
 


#=============================================================================
          ###CALL FOR DOWNLOAD CONCATENATION AND CHECK STRUCTURE###
#=============================================================================

project_duration = 1194




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

year = start_date_list[0] #returns year as a string
month = start_date_list[1]#returns month as a string
day = start_date_list[2] #returns day of month as a string
hour = start_date_list[3]
mins = start_date_list[4]


julian = datetime.datetime(int(year), int(month), int(day), 00, 00, 00)\
.timetuple().tm_yday


#calendar_date function returns date from day of the year


if time_dur == 86400:    
    
    for k in range(0, project_duration):
        year, month, day = calendar_date(int(year), julian)

        
        
        for station_name in stations_list:
        
            #reset time_dur to one hour or 3600 seconds in order to not time out!
            time_dur = 3600
            ###RESET hour to original in order for the function to role through###
            hour1 = hour
            print(hour)
            
            for i in range(0, 24):
        
                print("Hour currently being processed: %s" %(str(hour1)))
                
                download(year, month, day, hour1, mins, time_dur, station_name)
                #change type to int, add hour to the int, the change back to string in 
                #order for the download function to work!
                hour1 = int(hour1); hour1 += 1; hour1 = str(hour1)
        
            start_time = starttime(year, month, day, hour, mins)
        
            merge_delete(year, month, day, hour, mins, time_dur, station_name)
        
            input_directory = merge_delete(year, month, day, hour, mins, time_dur, \
            station_name)
        
            #split_channels(time_dur, station_name, input_directory, data_struct)
            

            metadata(input_directory)
            
            split_downsample(input_directory)
            
            move(input_directory)
            
            move_month(input_directory)

        
  #      day1 = int(day1); day1 += 1; day1 = str(day1)
        
        
        if julian < 365:
            julian += 1; #go to the next day's data, unless it's the end of the year!
            print("Currently Julian Day of the year is: %d" %(julian))
        else:
            julian = 1; year = str(int(year) + 1)

    
    merge_months()

        
        #metadata(input_directory)
        
        #response(input_directory)
                

        ###THIS ELIF IF STILL UNDER CONSTRUCTION###
elif int(float(time_dur) / 3600) > 1:
        
        for station_name in stations_list:
    
            download(year, month, day, hour, mins, time_dur, station_name)

    
else:
    
    for station_name in stations_list:
    
    
        download(year, month, day, hour, mins, time_dur, station_name)
    
        input_directory = download(year, month, day, hour, mins, time_dur, \
        station_name)
            
        split_channels(time_dur, station_name, input_directory, data_struct)



    
        



#if split != None:
    
#    quit()

#else:
    
#    while check == False:
                
#        split = raw_input("\nWould you like to split the channels of your data? also require data from another station? (Y/N) ")
            
#        if YN == 'Y': # or 'yes' or 'YES' or 'Yes' or 'y':
#            check = False
                    
#        elif YN == 'N': # or 'no' or 'NO' or 'No' or 'n': 
#            check = True
                    
#        else:
#            print("\nPlease enter either Y or N as your answer")
#            check = False
#            check2 = False 
        
        
        #metadata(input_directory)
    
        #response(input_directory)


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

#julian = datetime.datetime(int(Y), int(M), int(D), 00, 00, 00)\
#.timetuple().tm_yday

#=============================================================================
                            ###START ATTACH METADATA###
#=============================================================================


#=============================================================================
                            ###END ATTACH METADATA###
#=============================================================================


#=============================================================================
                            ###END PROGRAM###
#=============================================================================