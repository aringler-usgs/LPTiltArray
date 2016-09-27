#!/usr/bin/env python
###############################################################
# Intra-Sensor Coherence plots for Max Rohde's Long-Period Tilt Paper
# Written by Max Rohde and Adam Ringler
# Additional code for the location was provided by Emily Wolin

import pickle
import sys
import numpy as np
import matplotlib.cm as cm
from scipy import signal
from obspy.core import UTCDateTime, read, Stream

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import pylab

import seaborn as sns

sns.set(font_scale=2.5)

font = {'family' : 'sans-serif',
        'weight' : 'medium',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rc('axes',edgecolor='black')


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
        
# In meters E/W first, N/S second E,N = 0,0
distance = [([0.5,0.],[.5,.125],[.5,.25],[.5,.325],[.5,.4],[.5,.525],[.5,.65],[.5,.725]),
            ([0.,.5],[.125,.5],[.25,.5],[.325,.5],[.4,.5],[.525,.5],[.65,.5],[.725,.5]),
            ([.2,.1],[.4,.4],[.8,.9],[.6,.4],[.8,.1],[.6,.6],[.4,.6],[.1,.8]),
            ([.7,.4],[.8,.4],[.9,.4],[.8,.5],[.9,.5],[.7,.6],[.8,.6],[.9,.6])]



length= 80000        

overlap= 0.5
color_list = ['r','b','g','c','m','y','k','chartreuse']
debug= False   
direction_list = [0,1,2]

day_list = [
    (170,1466208000,1466244000,1466269200),

    (180,1467072000,1467079200,1467104400),

    (186,1467590400,1467597600,1467615600),

    (192,1468108800,1468130400,1468152000),
]
min_max_array= [(20.0,30.0),(30.0,50.0),(50.0,100.0),(100.0,500.0)]    
stas = ['XX_FBA1', 'XX_ENG4', 'XX_ENG7', 'XX_ENG5']
locs = ['00','10']


for idx, day_entry in enumerate(day_list):
    array_number = idx + 1
    print('On array: ' + str(array_number))
    for direction in direction_list:

        day = day_entry[0]
        trim = True
        start_time = UTCDateTime(day_entry[2])
        end_time = UTCDateTime(day_entry[3])    
        st = Stream()
        for sta in stas:
            for loc in locs:
                st += read("/tr1/telemetry_days/" + sta + 
                    "/2016/2016_" + str(day) + "/" + loc + "_BH" + str(direction) + ".512.seed")

        st.trim(start_time,end_time)        
        st = choptocommon(st)        

        delta = st[0].stats.delta
    

        fs = st[0].stats.sampling_rate
        # sampling freq of time series        

    

        for period_ranges in min_max_array:        

            coh_means = []
            coh_stds = []
            for n in range(0, 8):
                for m in range(0, 8):            
                    x = st[n].data
                    y = st[m].data
                    f, Cxy = signal.coherence(x, y, fs, nperseg=length,noverlap=overlap)
                    minper = period_ranges[0]
                    maxper = period_ranges[1]
                    per = 1.0/f
                    mic_coh = Cxy[(minper <= per) & (per <= maxper)]
                    micper = per[(minper <= per) & (per <= maxper)]
                    coh_mean = np.mean(mic_coh)
                    coh_std = np.std(mic_coh)            
                    coh_means.append(coh_mean)
                    coh_stds.append(coh_std)
            # We now have our array of coherences as well as distances
            print(distance[idx])
            x=[]
            y=[]
            for ele in distance[idx]:
                x.append(ele[0]*9)
                y.append(ele[1]*9)
            fig = plt.figure(figsize=(10,5))
            for n in range(0,8):
                sp = fig.add_subplot(2,4,n+1, aspect='equal')
                cohs = np.array(coh_means[n*8:(n+1)*8])
                s= plt.scatter(x, y, c=cohs, cmap=cm.magma, marker='s', vmin=0.,
                        vmax=1., s=100*cohs**2, linewidth=2)
                
                plt.scatter(x[n], y[n], c=1., cmap=cm.magma, marker='s', 
                    vmin=0., vmax=1., s=100*cohs[n]**2, linewidth=2, edgecolor='k')
                ax = plt.gca()
                ax.set_axis_bgcolor('white')
                sp.set_xlim(-1.1,8.9)
                sp.set_ylim(-1.1,8.9)
                plt.axvline(x=-1.1, linewidth=1, color='k')
                plt.axvline(x=8.9, linewidth=2, color='k')
                plt.axhline(y=8.9, linewidth=1, color='k')
                plt.axhline(y=-1.1, linewidth=2, color='k')
                sp.set_title('Coh. with sen. {0}'.format(n+1), fontsize=14,fontweight='bold')
                plt.xticks([])
                plt.yticks([])

            fig.subplots_adjust(right=0.8)
            cbar_ax=fig.add_axes([.85, 0.15, 0.05, 0.7])    
            cbar = fig.colorbar(s, ticks=[0, .25, .5, .75, 1], cax=cbar_ax)
            for l in cbar.ax.yaxis.get_ticklabels():
                l.set_weight("bold")
                l.set_fontsize(14)
            fig.suptitle('Coherence Array ' + str(array_number) + ' BH' + str(direction) + ' ' + str(int(minper)) + ' s to ' + str(int(maxper)) + ' s Period' , fontsize=16, fontweight='bold')    

            fig.savefig("WolinPltDay"+str(day)+'_GRID_Plot_BH'+ str(direction) + '_from_' + (str(start_time)[:-11]).replace(':','') + '_to_'+ (str(end_time)[:-11]).replace(':','') +'_'+ str(int(period_ranges[0])) + '_' + str(int(period_ranges[1])) +'.jpg',format='jpeg',dpi=400,bbox_inches='tight')
            fig.savefig("WolinPltDay"+str(day)+'_GRID_Plot_BH'+ str(direction) + '_from_' + (str(start_time)[:-11]).replace(':','') + '_to_'+ (str(end_time)[:-11]).replace(':','') +'_'+ str(int(period_ranges[0])) + '_' + str(int(period_ranges[1])) +'.eps',format='eps',dpi=800,bbox_inches='tight')
            plt.close()








