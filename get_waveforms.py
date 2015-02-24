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


def get_waveforms(starttime = str() , duration, Surface = True, Borehole = True, metadata = True, filestruct =  Channel,Stat, Net, Sta, BUD, ):

# ====================================================
# Sort out time from user input
# ====================================================


# ========================
# Data download
# ========================


# ========================================
# Attach metadata
# ========================================



# ========================================
# Output
# ========================================

if filestruct == "Single":
    
    elif
        



