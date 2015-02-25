#!/usr/bin/python -u
"""
[Advice: run this script using python with unbuffered output:
`python -u crosscorrelation.py`]
This script reads seismic waveform data from a set of stations, and
 coordinates, number of days, inter-
                  station distance etc.
- .stations.txt = general information on the stations: coordinates,
                  nb of cross-correlations in which it appears, total
                  nb of days it has been cross-correlated etc.
- .png          = figure showing all the cross-correlations (normalized to
                  unity), stacked as a function of inter-station distance.
"""

import obspy.core.trace
from obspy.signal import cornFreq2Paz
import obspy.xseed
import numpy as np
from obspy import read, read_inventory
from obspy import station
from obspy.station import stationxml


def get_waveforms(starttime = str() , duration, Surface = True, Borehole = True, metadata = True, filestruct =  Channel,Stat, Net, Sta, BUD, zerofill = True):

# ====================================================
# Sort out time from user input
# ====================================================

# ========================
# Data download
# ========================

# ========================================
# Attach metadata
# ========================================


#use file name produced from download script to read .mseed file from wherever this has been saved.
#alternatively I imagine eventually we could pull files down from seishub instead.

#currently reading an example miniseed file in as a stream
st = read("2014-10-17UTC8_13.mseed")
x = len(st)

#CORRECTING FOR NETWORK CODE AND LOCATION

network_change1 = 'UM'
network_new_name1 = 'BM'
network_change2 = '01'
network_new_name2 = 'UM'
location_blank = ''


for i in range(0, x):
    tr = st[i]
    
# removes LOCATION so it is blank, as listed in the metadata files (regardless of what it was previously)
   
    tr.stats["location"] = location_blank
    
# Changes BOREHOLE network codes from UM to BM and SURFACE network codes from 01 to UM 
    net = tr.stats["network"]
    if network_change1 in net:
        tr.stats["network"] = network_new_name1
    elif network_change2 in net:
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

#import and read metadata file, before attaching it to our stream
#next step is to pull from github or include option here to update metadata file.

metadata = stationxml.read_StationXML("UOM.xml")
#can use'("/path/to/UOM.xml")' to direct script to the correct location

st.attach_response(metadata)

#then save the file, or push to seishub, open in obspyck etc.
# ========================================
# Process: fill, buffer, highpass lowpass, remove response, spectral whitening
# ========================================

# ========================================
# Output
# ========================================

if filestruct == "Single":
    
    elif
        



