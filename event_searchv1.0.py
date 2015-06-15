import numpy as np
import matplotlib.pyplot as plt
from obspy import read 
from obspy.signal.trigger import *
import sys
import os, os.path
import glob
import datetime
from obspy.station import stationxml
import scipy



def trig_conditions(tr):
    
    """
    Function outputs a trigger value based on maximum of characteristic 
    function for a given waveform series. This means that if there is a 
    larger than normal earthquake that occurs during this waveform, it'll 
    take that as peak and the trigger detections won't be thrown off 
    because every detection is based on the maximum amplitude during 
    that imported wave-series. 
    
    The key is to find the optimum seismic waveform time-series length to 
    detect as many earthquakes as possible.
    
    Note that perc1 MUST be larger than perc2 and BOTH must be 0 < perc < 1
    """
    

    
    ##############
    #TRACE RELATED
    ##############
    avg_tr = np.average(tr)
    
    shifted_tr =  tr - avg_tr   
    
    abs_tr = abs(shifted_tr)
        
    sort_tr =  np.sort(abs_tr)
            
    tr_99 = sort_tr[int(0.999*len(sort_tr))]

    #95% noise to signal area ratio calculation
    area_big = tr_99 * int(0.95*len(sort_tr))
    
    area_small = np.sum(sort_tr[:int(0.95*len(sort_tr))])
    
        
    area_condition = area_small /  area_big; print area_condition
    
    # if signal is quite noisy, use recursive STALTA, else use classic STALTA
    
    df = tr.stats.sampling_rate
    
    
    if area_condition <= 0.1:
        
        print("CLEAN SIGNAL!\n")

        ctf = abs(1-carlSTATrig(tr.data, int(5 * df), int(10 * df), 0.0001, 0.0001))
        
        sort_ctf = np.sort(ctf)
    
        ctf_low = sort_ctf[int(0.80*len(sort_ctf))]
    
        ctf_99 = sort_ctf[int(0.999*len(sort_ctf))]
    
        max_trig = ctf_99
        
        trig_on = ctf_low + 0.4 * (max_trig - ctf_low)

        trig_off =  0.8 * trig_on
        
        if not max_trig > 1.25 * ctf_low:
            trig_on = 0; trig_off = 0;
        
        
    
    elif 0.1 < area_condition < 0.16:
        
        print("NOISY SIGNAL!\n")

        ctf = abs(1-recSTALTA(tr.data, int(5 * df), int(10 * df)))
        
        sort_ctf = np.sort(ctf)
    
        ctf_low = sort_ctf[int(0.90*len(sort_ctf))]
    
        ctf_99 = sort_ctf[int(0.999*len(sort_ctf))]
    
        max_trig = ctf_99
        
        trig_on = ctf_low + 0.4 * (max_trig - ctf_low)

        trig_off = 0.8 * trig_on
        
        if not max_trig > 1.25 * ctf_low:
            trig_on = 0; trig_off = 0;
    
    else:
        print("TRACE TOO NOISY TO DETECT SIGNAL!")
        trig_on = 0; trig_off = 0;
        ctf = abs(1-recSTALTA(tr.data, int(5 * df), int(10 * df)))

    



        

    
    
    ############
    #ctf RELATED
    ############

    
    #set conditional for noise. This is MINIMUM threshold.
    #ctf_avg = np.average(ctf)
    #ctf_mode =  scipy.stats.mode(ctf)
    #ctf_std =  np.std(ctf)
    

    
    
    #plt.figure()
    #plt.plot(sort_ctf)
    #plt.show()

    

    
    
    
    return trig_on, trig_off, ctf
    

def UTC_times(times, 
              trace, 
              diff_thres = 30.0):
    
    """
    Function that takes the output of the obspy.signal.trigger function 
    triggerOnset() and computes an average for points that do no differ 
    by a certain difference threshold (diff_thres), and gives back any 
    individual points that differ by more than this threshold on either 
    side.
    
    Finally the function converts this output into seconds and then 
    adds those on to the starting time of the trace. 
    """
    # set times values to seconds
    
    #AUTOMATE THIS SECTION!
    #CHECK THAT THIS IS CORRECT
    times = times / trace.stats.sampling_rate
    #remove unwanted parts of times numpy array 
    times = times[:,0]
    
    #remove the first instance of time because it is 
    #somehow always of the wrong format!
    #times = np.delete(times, 0)    
    
    event_times = []
    event = [times[0]]
    
    start_time = trace.stats.starttime
    
    #for item in times:
    #    print start_time + item

    for i in range(1, len(times)):
    
        # check if two events in times array have a difference < diff_thres, 
        #if not, run average of those times, if so append that events to a 
        #new events_times list
        
        #time_diff = times[i + 1] - times[i]
        
        time_diff = times[i] - times[i-1]

        #save info until events are far enough apart! 
        if time_diff < diff_thres:

            event.append(times[i])
            
        
        #raise conditional for if events are far enough apart! 
        else:

            event_start = event[0] - 2 #minus 5 seconds
            event_end = max(event) + 2 #add 5 seconds

            event_times.append([event_start, event_end])
            
            event = [] 
            
            event.append(times[i])

        #if event still contains something for any reason, add it to event times
    if len(event) > 0:            
        event_start = event[0] - 2 #minus 5 seconds
        event_end = max(event) + 2 #add 5 seconds
        event_times.append([event_start, event_end])
        event = [] 
            


        #if len(event_times) == 0 and len(event) > 0 or time_diff > diff_thres and len(event) > 0:
    
            #event_times.append(sum(event) / len(event))
            
        #    event_start = event[0] - 2 #minus 5 seconds
        #    event_end = event[-1] + 2 #add 5 seconds
            
        #    event_times.append([event_start, event_end])
            
        #    event = []
        
            #event_times.append(times[i])
    
   #     else:
   #         event.append(times[i])
    

    UTC_events = []

    #earthquake length threshold is 10 seconds and above!
    eq_len = 0#5.0

    for i in event_times:
        estart = start_time + i[0]
        eend = start_time + i[1]
        
        if eend - estart > eq_len:
            UTC_events.append([estart, eend])
    
    #UTC_events = np.unique(np.asarray(UTC_events))

    
    return UTC_events


def trigger_times(trace, 
                  STA        = 1, 
                  LTA        = 10, 
                  ratio      = 0.0001,
                  quiet      = 0.0001,           
                  trig_auto  = True,
                  perc1      = 0.2,
                  DIFF_THRES = 60.0,
                  show       = False):
    
    """
    Function that returns UTCdatetime objects of events from a given input 
    trace. It uses carlSTATrig to find the characteristic function of the 
    trace, then plots the trigger findings using the plotTrigger and finally
    the function UTC_events is used to compute the event timings in UTC 
    format from the ouput of the function triggerOnset.
    
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
    
    

    

    #plotTrigger(trace, cft, 1.2, 0.5) 


    if trig_auto:# and trig_on != None:
        trig_on, trig_off, ctf = trig_conditions(tr)
        
        # use triggerOnset function to calculate trigger on and off times for a given time frame
    times = triggerOnset(ctf, trig_on, trig_off, max_len = DIFF_THRES)
        
    if show:
            # use plotting the trigger functions using plotTrigger
        plotTrigger(tr, ctf, trig_on, trig_off, show=True) 
        
    if len(times) > 0:

        UTC_events = UTC_times(times, trace = tr, diff_thres = DIFF_THRES)
        
    return UTC_events

    #else:
        
    #    UTC_events = []
    #    return UTC_events

        


def save_txt(txt_fname, line):
    """
    Function that takes a string as a line input, and the name of a file,
    then opens the text file, appends the string to the next line and closes
    the file when done
    """
        
    #if os.path.isfile(txt_fname):
    #    print("File {} already exists ...\n".format(txt_fname))
    #    print("Adding line: {} to file: {}".format(line, txt_fname))
        
    #else:
    #    print('Creating new text file ... \n') 
    #    print("Adding line: {} to file: {}".format(line, txt_fname))

    try:
            
        txt_file = open(txt_fname, 'a+')
        
        line_string = "{}\n".format(line)

        txt_file.write(line_string)
        
        txt_file.close()
        
        
            #print("\nLine added to file {} correctly".format(txt_fname))
        
    except Exception as e:
        a = e
        
def paths_sort(path):
    """
    Function defined for customised sorting of the abs_paths list
    and will be used in conjunction with the sorted() built in python
    function in order to produce file paths in chronological order.
    """
    base_name = os.path.basename(path)
    
    stat_name = base_name.split('.')[0]    

    date = base_name.split('.')[1]
    
    try:
        date = datetime.datetime.strptime(date, '%Y-%m-%d')
        
        return date, stat_name
    except Exception as e:
        print(e)
        
def paths(folder_path, extension):
    """
    Function that returns a list of desired absolute paths called abs_paths
    of files that contains a given extension e.g. .txt should be entered as
    folder_path, txt. This function will run recursively through and find
    any and all files within this folder with that extension!
    """

    abs_paths = []
    
    for root, dirs, files in os.walk(folder_path):
        
        for f in files:
            
            fullpath = os.path.join(root, f)
            
            if os.path.splitext(fullpath)[1] == '.{}'.format(extension):
                
                abs_paths.append(fullpath)

    abs_paths = sorted(abs_paths, key=paths_sort)
       
    return abs_paths
    
    

def metadata(read_dir):
    

#currently reading an example miniseed file in as a stream
    try:
        st = read(read_dir)
        #print(st)
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

        metadata_path = 'UOM.xml'
        
        if not os.path.isfile(metadata_path):
            os.system('wget https://raw.githubusercontent.com/unimelb-geophys/metadata/master/UOM.xml')


        metadata = stationxml.read_StationXML(metadata_path)
        #print(metadata)

        #st.write("")

        #attach metadata to stream
        st.attach_response(metadata)
        #st.remove_response()

    
        #delete downloaded metadata file 
        #os.system('rm -r UOM.xml')


        #st.write(read_dir, format="MSEED")
    

    
    except Exception as e:
        print("\nOops, there was a problem with attaching the metadata\n")
        print(e)
        
    
    
    return st
    
    


######################
# SET INPUT PARAMETERS
######################

# set events output txt file name
event_fname = "event_outputs.txt"

# set error txt file name
error_fname = "errors.txt"

# set folder with all waveform files in it. Can be recursive! 
#folder_path = '/storage/ABE/Gippsland/MONTH_DOWNLOAD'

#folder_path = '/home/boland/Dropbox/University/UniMelb/Research/SIMULATIONS/Triggers/trigger_miniseeds'
folder_path = '/home/boland/Desktop/Link to SIMULATIONS/Triggers/data/small_borehole_quakes'

# set file extension to mseed

extension = 'm'

 
# set desired component e.g. E, N or Z
comp = 'Z' 
channel = 'EHZ'

################
# CALL FUNCTIONS
################

# get necessary absolute path string list for scanning all waveforms
abs_paths = paths(folder_path, extension)
#abs_paths = ['/home/boland/Dropbox/University/UniMelb/Research/SIMULATIONS/Triggers/Hawaii_low_mag_BHZ_IU_POHA.mseed']

#abs_paths = ['/home/boland/Dropbox/University/UniMelb/Research/SIMULATIONS/Triggers/chch_earthquake.mseed']

counts = 0
with open(event_fname, "a+") as event_file:
    
    # run for loop on all files!
    for f in abs_paths:
        
        print("Processing file: {}\n".format(os.path.basename(f)))
        
        try:
            #read and attached responce with metadata() function
        
            #
        
            #st = metadata(f)
        
            st = read(f)
            try: 
                #st = st.select(station = "HOLS")
                
                #st = st.select(component = comp)
                
                
                tr = st[0]
   
                
                UTC_events = trigger_times(tr, DIFF_THRES = 30.0, show = True)

                all_events = []


                if len(UTC_events) > 0: 
    
                    for event in UTC_events:
                        
                        event_start = event[0]; event_end = event[1]

                        event_string = "%s,%s,%s,%s\n" %(tr.stats.network, tr.stats.station, str(event_start), str(event_end))
                        
                        sliced_st = st.slice(event_start, event_end)
                        
                        #save event to new events folder
                        if not os.path.exists("events"): os.makedirs("events")
                        
                        new_file_path = "events/{}_event{}.mseed".format(os.path.basename(f).split(".")[0], counts)

                        sliced_st.write(new_file_path, format="MSEED")
                        
                        counts += 1

                        
                        #read in new event stream

                        all_events.append(event_string)
            
            
                    
                # save events to text file
                if len(all_events) > 0:
                    event_file.writelines(all_events)
                    #print(all_events)

            
            except Exception as error:
                print("There are no channels within the stream of type: {}".format(channel))
                print(error)

            
        except Exception as error:
        # save all errors to txt file
        #save_txt(error_fname, error)
            a=error
            
            
print("file closed")


# automate time difference between earthquakes to be larger or smaller based
# on earthquake size! ask Abe about some sort of equation. 

# set a lower limit on the trigger maximum because the signal will give 
# loads of noise if there isn't a single earthquake in the time series! 

# set event time min and time max +- a certain automated amount of time as list
# then go through and save mseed of each event labelling it 
# 0 through to whatever the event number is! 

#NORMALISE THESE POINTS ABOUT ZERO!

#MAKE SURE YOU ADD IN A MINIMUM EARTHQUAKE TIME OCCURRENCE!