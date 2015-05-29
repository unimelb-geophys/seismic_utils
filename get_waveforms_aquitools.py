
# coding: utf-8
from __future__ import with_statement
import math


## Function to download merge all waveforms
### Setup time:

# In[41]:

Y = raw_input("enter year, e.g. 2010: ")
M = raw_input("enter Month, e.g. 02: ")
D = raw_input("enter Day, e.g. 07: ")
H = raw_input("enter hour (24 hr), e.g. 13: ")
m = raw_input("enter minute, e.g. 00: ")
dur = raw_input("enter duration in seconds, e.g. 180: ")


# In[42]:

import obspy


# In[43]:

t = obspy.UTCDateTime("%s-%s%-sT%s:%s" %(Y, M, D, H, m))
tend = t + float(dur)


# In[44]:

tend


# In[45]:

#In script we will promt user for year, month, day hour minute
#little earthquake at Moe
#t = obspy.UTCDateTime("2014-10-17T08:13")
#tend = t + 180


# In[46]:

year = t.year
#type(year)
day = t.strftime('%j')
hour = t.hour 
#type(day)


# In[47]:

#files are meant to be one hour long and seem to start at minute 00
#For broken files there will need to be some additions to the script 


# In[48]:

hour = str(hour).zfill(2)


# In[49]:

print(hour)


# In[50]:

#Borehole sites
#bstats = ["BD5E", "BD70","BD91"]
bstats = ["BD70","BD91","A346"]

# In[51]:

bpaths = []
for item in bstats:
    #path_st = "../home/reftek/Archive/%s%d/%s*/2/%d" % (year,day,item,hour)
    path_st = "../home/reftek/Archive/%s%s/%s/2/%s*" % (year,day, item, hour)
    bpaths.append(path_st)


# In[52]:

bpaths


# In[53]:

type(bpaths[1])


# In[54]:

#the following works...
#scp root@agoslog.earthsci.unimelb.edu.au:../home/reftek/Archive/2014290/BD91/2/07* 07.reftek


### get borehole files

# In[55]:

import os
import shutil
dir = './bh_reftek'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)


# In[56]:


from fabric.api import *
from fabric.contrib.console import confirm
from fabric.api import env, run
from fabric.tasks import execute
from fabric.api import *

env.user = "root"


# In[57]:

def get_wave():
    for i in range(len(bstats)):
        bh_st = bpaths[i]
        name = bstats[i]
        name_st = "bh_reftek/%s" %name
        get(bh_st, name_st) 

execute(get_wave, host="128.250.59.212")


# In[58]:

#have to add Aquitools to path:
#/Users/dansandiford/Documents/programming/waveform_tools/AcquiTools-0.1.0/build/bin
#!./build/bin/ckreftek 090000000_0036EE80


# In[59]:

#!ckreftek -NUM ./bh_waves


# The simplest way of running UNIX command is to use os.system().
# os.system('command with args') passes the command and arguments to our system's shell. 
# By using this can actually run multiple commands at once and set up pipes and input/output redirections. :

# In[60]:

import os
import shutil
dir = './bh_mseed'
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)


# In[61]:

#this is writing 3 files for two of the stations, and only one for the other, don't know why
for i in range(len(bstats)):
    st_name = bstats[i]
    os.system('ckreftek -NUM -S%s bh_reftek/%s' % (st_name, st_name))


# In[62]:

#!ckreftek -NUM -SBD70 bh_reftek/*BD70


### read borehole files files

# In[63]:

from obspy import read
import fnmatch


# In[64]:

walk_dir = os.getcwd()


# In[65]:

pattern = 'BH*.m'
tracy = []
for root, dirs, files in os.walk(walk_dir):
    for filename in fnmatch.filter(files, pattern):
        print(filename)
        st = obspy.read(os.path.join(root, filename))
        tracy.append(st)
    


# In[66]:

for i in range(0,len(tracy)-1):
    if i == 0: 
        trrpt = tracy[0]
    trnew = trrpt + tracy[i+1]
    trrpt = trnew


# In[67]:

stb = trrpt.slice(t, tend)


# In[68]:

#stb.write("./bh_mseed/borehole.mseed",format='MSEED')


# In[69]:

#stb.plot()


# In[70]:

stbv = stb.select(channel="BHZ")


# In[71]:

#stbv.plot()


### Delete unnecessary directories

# In[72]:

import glob, os
test = './UM*'
r = glob.glob(test)
for i in r:
   shutil.rmtree(i)


### Get surface stations

# #this is the form of a wget request to the server
# 
#     wget --user=eq --password=event55s "http://agos1.kelunji.net/eqserver/eqwaveextractor?   year=2014&month=10&day=20&hour=03&minute=00&duration=3&servernum=0&conttrig=0&sitechoice=list&sitelist=+BRAD+PEGM+PEGU+S88U+WALM+MRAT+SDAN+ADNL+ADE+SOMU+LOCU+AKRAN+HOLS+CLIF+NARR+AMRDN+HODL+UOM+MELU+&siteradius=&closesite=&radius=&latitude=&longitude=&fileformat=miniseed&getwave=Get+Waveform" -O agos.ms

# In[73]:

mdur = (float(dur))/60
mdur = math.ceil(mdur)
mdur = int(mdur)


# In[74]:

#eqstring = "http://agos1.kelunji.net/eqserver/eqwaveextractor?   year=%d&month=%d&day=%d&hour=%d&minute=%d&duration=%d&servernum=0&conttrig=0&sitechoice=list&sitelist=+BRAD+PEGM+PEGU+S88U+WALM+MRAT+SDAN+ADNL+ADE+SOMU+LOCU+AKRAN+HOLS+CLIF+NARR+AMRDN+HODL+UOM+MELU+&siteradius=&closesite=&radius=&latitude=&longitude=&fileformat=miniseed&getwave=Get+Waveform" % (t.year, t.month,t.day, t.hour, t.minute, dur ) 


# In[75]:

eqstring = "http://agos1.kelunji.net/eqserver/eqwaveextractor?year=%d&month=%d&day=%d&hour=%d&minute=%d&duration=%d&servernum=0&conttrig=0&sitechoice=list&sitelist=+SOMU+LOCU+KRAN+HOLS+CLIF+NARR+MRDN+HODL+UOM+MELU+&siteradius=&closesite=&radius=&latitude=&longitude=&fileformat=miniseed&getwave=Get+Waveform" % (t.year, t.month, t.day, t.hour, t.minute, mdur) 


# In[76]:

print(eqstring)


# In[77]:
#removed password, won't work with it
fin_string = 'wget --user=eq --password= "%s" -O surface.mseed'  % eqstring


# In[78]:

#fin_string


# In[79]:

#need wget
os.system(fin_string)


# In[80]:

from obspy import read
sts = read("surface.mseed")
#sstf= obspy.read('surface.ms')
#sstf = trrpt.slice(t, tend)


# In[81]:

#sts.plot()


# In[82]:

merge = stbv + sts


# In[83]:


#merge.plot()


# In[84]:

tt = t.utctimetuple()


# In[85]:

tstr = "%s-%s-%sUTC%s_%s" %(tt[0], tt[1], tt[2], tt[3], tt[4])


# In[86]:

tstr


# In[87]:

fmerge = merge.slice(t, tend)
fmerge.filter("bandpass", freqmin=1.5, freqmax=40)

fmerge.plot(outfile='%s.png' %tstr)
fmerge.write("%s.mseed" %tstr,format='MSEED')


import matplotlib.pyplot as plt

for i in range(0, len(fmerge)):
    ax = plt.subplot(len(fmerge),1,i)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    tr = fmerge[i]
    d = tr.data
    t = tr.times()
    ax.plot(t, d)

plt.savefig('dirty.png', dpi = 300)


### remove more accessory stuff

# In[88]:

import glob, os
test = './bh*'
r = glob.glob(test)
for i in r:
   shutil.rmtree(i)


# In[89]:

os.remove("surface.mseed")


# In[ ]:

os.remove("reflog.log")

