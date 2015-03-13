import os, sys
sys.path.append("/home/abe/anaconda/lib/python2.7/site-packages")

from obspy import read, read_inventory
from obspy import station
from obspy.station import stationxml

directory = "/storage/AGOS/MONTH/2014-10"

#use file name produced from download script to read .mseed file from wherever this has been saved.
#alternatively I imagine eventually we could pull files down from seishub instead.


for files in os.listdir(directory):


#currently reading an example miniseed file in as a stream

    read_dir = "%s/%s" %(directory, files)
    st = read(read_dir)
    print(st)
    x = len(st)



    #CORRECTING FOR NETWORK CODE AND LOCATION

    network_change1 = 'UM'
    network_new_name1 = 'UM'#'BM'
    network_change2 = '01'
    network_new_name2 = 'UM'
    network_change3 = 'A'
    network_new_name3 = 'UM'
    network_change4 = ''
    network_new_name4 = 'UM'
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


    metadata = stationxml.read_StationXML("/storage/AGOS/XML/UOM.xml")#'UOM.xml')
    #print(metadata)

    #st.write("")

    print(st)
    #attach metadata to stream
    st.attach_response(metadata)
    #st.remove_response()


    #delete downloaded metadata file 
    os.system('rm -r UOM.xml')

    #CAN CHECK EVERYTHING IS WORKING

    st.write(read_dir, format="MSEED") 