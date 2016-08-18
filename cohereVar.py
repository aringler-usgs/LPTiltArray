#!/usr/bin/env python

###############################################################
# Coherence estimates for Max Rohde's Long-Period Tilt Paper
# Written by Max Rohde and Adam Ringler

import sys
import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
import numpy as np

from scipy import signal
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import pazToFreqResp
from itertools import combinations
from matplotlib.mlab import csd

# Setup fonts
font = {'family' : 'sans-serif',
        'weight' : 'medium',
        'size'   : 18}

mpl.rc('text',usetex=True)

def choptocommon(stream):
    #A function to chop the data to a common time window
    stimes = []
    etimes = []        

    for trace in stream:
        stimes.append(trace.stats.starttime)
        etimes.append(trace.stats.endtime)
    newstime = stimes[0]
    newetime = etimes[0]        

    for curstime in stimes:
        if debug:
            print(curstime)
        if curstime >= newstime:
            newstime = curstime        

    for curetime in etimes:
        if debug:        
            print(curetime)
        if curetime <= newetime:
            newetime = curetime        

    if debug:
        print(newstime)
        print(newetime)
        print(stream)
    for trace in stream:    
        trace.trim(starttime=newstime,endtime=newetime)
    if debug:
        print(stream)
    return stream        


minper = 50.
maxper=100.
color_list = ['r','b','g','c','m','y','k','chartreuse']
debug= False 
direction_list = [0,1,2]

sensorType = ["STS-2LG_A", "STS-2HG_B", "STS-2LG_A", "STS-2HG_B", "STS-2HG_A", "STS-2LG_B" ,"STS-2HG_A", "STS-2LG_B"] 
stas = ['XX_FBA1', 'XX_ENG4', 'XX_ENG7', 'XX_ENG5']
locs = ['00','10']

length = 80000        

overlap = 0.5


day_list = [
    (170,1466208000,1466244000,1466269200),

    (186,1467590400,1467597600,1467615600),

    (192,1468108800,1468130400,1468152000),

    (180,1467072000,1467079200,1467104400),
]

# For each day go through all components
for pidx, day_entry in enumerate(day_list):
    for direction in direction_list:
        day = day_entry[0]
        trim = True
        start_time = UTCDateTime(day_entry[2])
        end_time = UTCDateTime(day_entry[3])        
        # read in all the data
        st = Stream()
        for sta in stas:
            for loc in locs:
                st += read("/tr1/telemetry_days/" + sta + 
                "/2016/2016_" + str(day) + "/" + loc + "_BH" + str(direction) + ".512.seed")

      

        st.trim(start_time,end_time)
        st = choptocommon(st)
        fs = st[0].stats.sampling_rate

        # Get all the combinations and compute the coherence
        for comb in combinations(range(0,8),2):
            idx = comb[0]
            idx2 = comb[1]
            #Different sensors compare the coherence
            f, Cxy = signal.coherence(st[idx].data, st[idx2].data, fs, nperseg=length,noverlap=overlap)
            per = 1.0/f
            if 'meanCoh' not in vars():
                meanCoh = Cxy
            else:
                meanCoh = np.vstack((Cxy, meanCoh))
                    
        # We now have an array of means and stds for the coherences              
        stdCoh = np.std(meanCoh,axis=0)
        meanCoh = np.mean(meanCoh,axis=0)
        # Now we need to plot everything
        fig=plt.figure(1)
        if pidx > 0 and direction == 0:
            ax=plt.subplot(4,1,pidx+1,sharex=ax)
        elif direction == 0:
            ax=plt.subplot(4,1,pidx+1)
        ax.errorbar(1.0/f, meanCoh,stdCoh,linewidth=0.8, label='BH' + str(direction))
        ax.set_xscale('log')
        plt.xticks([])
        plt.xlim((10,500))
        plt.ylim((0.0,1.0))
        plt.text(12,.5,'Array '+ str(pidx+1),fontsize=12)
        plt.yticks([.2, .4, .6, .8])
        
        if pidx == 3:
            plt.xlabel('Period (s)')
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.7),
                fancybox=False, shadow=False, ncol=3)
            #plt.xticks([10., 100.])
            #plt.ylabel()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.05,
                 box.width, box.height * 0.95])
        fig.text(0.05, 0.5, 'Coherence ($\gamma^2$)', ha='center', va='center', rotation='vertical')

        del meanCoh
plt.xticks([10., 100.])        
plt.savefig('Coheren'+'.eps',format='eps',dpi=1000,bbox_inches='tight')
