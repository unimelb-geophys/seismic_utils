# -*- coding: utf-8 -*-
"""
Created on Sat May  9 17:06:22 2015

@author: boland
"""

import os, sys
sys.path.append("/home/abe/anaconda/lib/python2.7/site-packages")

from obspy import UTCDateTime
import datetime as dt
from obspy import read, read_inventory
from obspy import station
from obspy.station import stationxml
from obspy.xseed import Parser
meta_dir = '/home/boland/Dropbox/University/UniMelb/AGOS/METADATA/metadata/UOM.dataless'

meta_dir = "/home/boland/Dropbox/University/UniMelb/AGOS/METADATA/metadata/UOM.xml"

metadata = stationxml.read_StationXML(meta_dir)
    
    
    
t0 = UTCDateTime("2014-10-01T00:00:00.0")  # date at time 00h00m00s

dirs = '/storage/ANT/PROGRAMS/ANT_OUTPUT/INPUT/DATA/2014-10/UM.HOLS.EHZ.mseed'

st = read(dirs, starttime=t0, endtime=t0 + dt.timedelta(minutes=45))
print(st)
#inv = Parser(meta_dir)


#paz = inv.getPAZ("EHZ", t0)

    #attach metadata to stream
st.attach_response(metadata)
st.remove_response()

#CAN CHECK EVERYTHING IS WORKING

