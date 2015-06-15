
# coding: utf-8

# # Triggered 3 min download Scipt 

# A script with two primary functions: To locate potential events using the STA/LTA triggering function, and to then download a 3 minute file containing this event arrivalat surface and borehole seismic stations in our network.
# 
# **N.B. this is still very much a work in progress!**
# 
# ### Additions/changes to Ben's initial triggering code  ###
# 
# Will primarily run this script over data for a period of time collected at just one station, so the arrival times for a particular event at individual station is not as impotant as defining a window within which all station recorded the event arrival. For this purpose, the ouput of this script is to produce a string of times which represent 90sec prior to the origin time of events. This then becomes the start time (and date) for the second half of the script - the download phase. Also added an ouput txt file to be written containing this data.
# 
# ### Additions/changes to Dan's download code  ###
# Firstly, the input of a password into the call to eqserver can now be automated, which is nice. Input into theis function is changed from raw input into the output from the first half of theis script - as mentioned above. Modified some print functions so that while the script is runing, there is more feedback about which stage of the process it is at. 
# 
# ###Still in progress###
# Need to modify the output so that the folder structure is more user friendly. The code written to remove unneccessary files after the conversion process is flawed, and was deleting the converted miniseed files, so this has been disabled. Will need to make a new version of this.
# Still some errors in the borehole downloads if files are missing - need to write in an exception. Also need to add TRPU to the list. LOYU has 2 recorders there and this causes errors, need to program a way around this or sacrifice its data.

# In[20]:

# Triggered 3 min download Scipt 


### Usage: 


# coding: utf-8

get_ipython().magic(u'matplotlib inline')
import obspy
import glob
import os
import shutil
import math
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from obspy import read 
from obspy.signal.trigger import *
from __future__ import with_statement
from fabric.api import *
from fabric.contrib.console import confirm
from fabric.api import env, run
from fabric.tasks import execute
from obspy import UTCDateTime

# firstly combine everything above into one function. The inputs of this function will be 
#the traces and the outputs the event times





### DEFINE FUNCTIONS ###
def UTC_times(times, 
              trace, 
              diff_thres = 30.0,):
    
    """
    Function that takes the output of the obspy.signal.trigger function triggerOnset() and computes
    an average for points that do no differ by a certain difference threshold (diff_thres), and gives 
    back any individual points that differ by more than this threshold on either side.
    
    Finally the function converts this output into seconds and then adds those on to the starting time
    of the trace 
    """
    
    # set times values to seconds
    times = times / 100.0 
    
    #remove unwanted parts of times numpy array 
    times = times[:,0]
    
    event_times = []
    event = []

    
    for i in range(0, len(times) - 1):
    
        # check if two events in times array have a difference < diff_thres, 
        #if not, run average of those times, if so append that events to a 
        #new events_times list
        
        time_diff = times[i + 1] - times[i]
        
        if time_diff > diff_thres:
        
            if len(event):    
        
                event_times.append(sum(event) / len(event))
                event = []
        
            event_times.append(times[i + 1])

    
        else:
            event.append(times[i])
            
            

    start_time = trace.stats.starttime
    UTC_events = []

## Abe changed start_time here to be 90 seconds earlier than triggered event
# so that this can be used as the initial time for the download script    
    
    for i in event_times:
        event = (start_time-90) + i      
        UTC_events.append(event)      
    return UTC_events


def trigger_times(trace, 
                  STA        = 5, 
                  LTA        = 10, 
                  ratio      = 0.001,
                  quiet      = 0.001,           
                  trig_on    = 1500.0,
                  trig_off   = 600.0,
                  DIFF_THRES = 30.0,
                  show       = False):
    
    """
    Function that returns UTCdatetime objects of events from a given input trace. It uses carlSTATrig to 
    find the characteristic function of the trace, then plots the trigger findings using the plotTrigger and finally
    the function UTC_events is used to compute the event timings in UTC format from the ouput of the function
    triggerOnset.
    
    Input parameters are:
    
    trace      - must be an obspy Trace object.
    STA        - must be an int object which defines the short term average time window in seconds.
    LTA        - must be an int object which defines the long term average time window in seconds.
    ratio      - must be a float object, the lower the ratio number, the more sensitive the 
                 characteristic function's output.
    quiet      - must be a float object, the lower the quiet number, the more sensitive the 
                 characteristic function's output.
    DIFF_THRES - must be a float object which defines the minimum window difference between successive events. 
                 this is then used by the UTC_times function in order to determine whether a times output
                 is within a singular event or it differs enough to be called a seperate event. 

    """
    

    # first produce the STA/LTA characteristic function ctf associated with a given trace
    ctf = abs(1-carlSTATrig(tr, STA, LTA, ratio, quiet))

    if show == True:
        # use plotting the trigger functions using plotTrigger
        plotTrigger(tr, ctf, trig_on, trig_off, show=True) 
    
    
    # use triggerOnset function to calculate trigger on and off times for a given time frame
    times = triggerOnset(ctf, trig_on, trig_off, max_len = 60.0)
    
    if len(times) >= 1:

        UTC_events = UTC_times(times, trace = tr, diff_thres = DIFF_THRES)
    
        return UTC_events
    
# define a folder with many miniseed files in it

folder_name = './trigger_miniseeds'

# get a list of all miniseed files within the directory above

file_list = glob.glob('{}/*.MSEED'.format(folder_name))


all_events = []
for f in file_list:
# import miniseed file as an obspy Stream object, st
    st = read(f)
    tr = st[0]

    UTC_events = trigger_times(tr)

    if UTC_events != None and len(UTC_events) >= 1: 
        for event in UTC_events:

#Defining inputs for download script 

#            Y = event.year
#            M = event.month
#            D = event.day
#            H = event.hour
#            m = event.minute
#            s = event.second
#            dur = 180

#will here use download function to download data from all active sites (inc boreholes) at this moment in time            
                        
            event_list = [event]
            all_events.append(event_list)
            
            
f = open("download_times.txt","w") 
for i in all_events: 
    event_string = '{}'.format(i[0])
    f.write(event_string)
    f.write("\n")

f.close()  



### Download Script ###


## Will use values from above, but defaults selected here for testing
# rather than relying on keeping the event times in memory, will reread in these values from the produced txt file.

Y = 2014
M = 12
D = 03
H = 01
m = 00
dur = 180

####change the t = line back later when is automoated for a number ot dateTimes

#t = obspy.UTCDateTime("%s-%s%-sT%s:%s" %(Y, M, D, H, m))
t = obspy.UTCDateTime(Y, M, D, H, m)
tend = t + float(dur)


year = t.year
day = t.strftime('%j')
hour = t.hour 


#files are meant to be one hour long and seem to start at minute 00
#For broken files there will need to be some additions to the script 


hour = str(hour).zfill(2)


#Borehole sites -- still some errors with LOYU, and need to add TRPU
#bstats = ["BD5E", "BD70", "BD91", "A346"]
bstats = ["BD91", "BD70"]


bpaths = []
for item in bstats:
    #path_st = "../home/reftek/Archive/%s%d/%s*/2/%d" % (year,day,item,hour)
    path_st = "../home/reftek/Archive/%s%s/%s/2/%s*" % (year,day, item, hour)
    bpaths.append(path_st)




#the following works...
#scp root@agoslog.earthsci.unimelb.edu.au:../home/reftek/Archive/2014290/BD91/2/07* 07.reftek


### get borehole files ####

dir = './bh_reftek'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)

env.user = "root"

def get_wave():
    for i in range(len(bstats)):
        bh_st = bpaths[i]
        name = bstats[i]
        name_st = "bh_reftek/%s" %name
        get(bh_st, name_st) 

execute(get_wave, host="128.250.59.212")



# The simplest way of running UNIX command is to use os.system().
# os.system('command with args') passes the command and arguments to our system's shell. 
# By using this can actually run multiple commands at once and set up pipes and input/output redirections. :

# In[60]:

dir = './bh_mseed'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)


#have to add Aquitools to path for this session: (May be a way to make this permanent and cut down on computational time?)
#os.system('PATH=$PATH:/home/abe/AcquiTools-0.1.3/')
#Despite working when running this from soource in the same dir, ckreftek couldn't be added
#to path with os.system
#instead manually make a copy of the ckreftek executable in usr/anaconda/bin - the default path 
    
for i in range(len(bstats)):
    st_name = bstats[i]
    os.system('ckreftek -NUM -S%s bh_reftek/%s \n' % (st_name, st_name))
    print 'converting', ('bh_reftek/%s...' % (st_name))
print 'finished borehole downloads'


### Get surface stations ###

#this is the form of a wget request to the server
#     wget --user=eq --password=event55s "http://agos1.kelunji.net/eqserver/eqwaveextractor?   year=2014&month=10&day=20&hour=03&minute=00&duration=3&servernum=0&conttrig=0&sitechoice=list&sitelist=+BRAD+PEGM+PEGU+S88U+WALM+MRAT+SDAN+ADNL+ADE+SOMU+LOCU+AKRAN+HOLS+CLIF+NARR+AMRDN+HODL+UOM+MELU+&siteradius=&closesite=&radius=&latitude=&longitude=&fileformat=miniseed&getwave=Get+Waveform" -O agos.ms


print "Commencing surface station downloads..."

mdur = (float(dur))/60
mdur = math.ceil(mdur)
mdur = int(mdur)


# In[75]:

eqstring = "http://agos1.kelunji.net/eqserver/eqwaveextractor?year=%d&month=%d&day=%d&hour=%d&minute=%d&duration=%d&servernum=0&conttrig=0&sitechoice=list&sitelist=+SOMU+LOCU+KRAN+HOLS+CLIF+NARR+MRDN+HODL+UOM+MELU+&siteradius=&closesite=&radius=&latitude=&longitude=&fileformat=miniseed&getwave=Get+Waveform" % (t.year, t.month, t.day, t.hour, t.minute, mdur) 

# removed password, won't work with it 
# *** UPDATE: is now working with password
fin_string = 'wget --user=eq --password=event55s "%s" -O surface.mseed'  % eqstring
print 'downloading waveforms from eq server...'

#need wget to be installed 
os.system(fin_string)


### Appending ###

walk_dir = os.getcwd()
pattern = "EH*"
tracy = []

print "merging waveforms to single file..."

#Search directory for files containing "EH" and append these to trace

for root, dirs, files in os.walk(walk_dir):
    for filename in fnmatch.filter(files, pattern):
        st = obspy.read(os.path.join(root, filename))
        tracy.append(st)
        
tt = t.utctimetuple()
tstr = "%s-%s-%sUTC%s_%s" %(tt[0], tt[1], tt[2], tt[3], tt[4])        
        
for i in range(0,len(tracy)-1):
    if i == 0: 
        trrpt = tracy[0]
    trnew = trrpt + tracy[i+1]
    trrpt = trnew

# Indented these since weren't defined initially
    
    stb = trrpt.slice(t, tend)
    stbv = stb.select(channel="EHZ")
    stb.write("./bh_mseed/borehole.mseed",format='MSEED')
    
    sts = read("surface.mseed")
    merge = stbv + sts
    
    fmerge = merge.slice(t, tend)
    fmerge.filter("bandpass", freqmin=1.5, freqmax=40)
    fmerge.plot(outfile='%s.png' %tstr)
    fmerge.write("%s.mseed" %tstr,format='MSEED')
     


# In[21]:

#### Remove accessory files    
print "Removing accessory files"

test = './bh*'
r = glob.glob(test)
for i in r:
   shutil.rmtree(i)

test = './UM*'
r = glob.glob(test)
for i in r:
   shutil.rmtree(i)

try:
    os.remove("surface.mseed")
except OSError:
    pass

try:
    os.remove("reflog.log")
except OSError:
    pass

try:
    os.remove("report")
except OSError:
    pass    
    
    
print "Task complete"


# In[ ]:



