#!/usr/bin/env python
###############################################################
# Intra-Sensor Coherence plots for Max Rohde's Long-Period Tilt Paper
# Written by Max Rohde and Adam Ringler

import pickle
import numpy as np

from scipy import signal
from obspy.core import UTCDateTime, read, Stream

from matplotlib.pyplot import (figure,xlabel,ylabel,title,subplot,legend,savefig)
import matplotlib
from matplotlib import pylab
from matplotlib import pyplot
import seaborn as sns

sns.set(font_scale=2.5)

font = {'family' : 'sans-serif',
        'weight' : 'medium',
        'size'   : 16}

matplotlib.rc('font', **font)


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

                    subplot(111)
                    x = st[n].data
                    y = st[m].data
                    f, Cxy = signal.coherence(x, y, fs, nperseg=length,noverlap=overlap)
                    xlabel('Period [s]')
                    ylabel('Coherence')            

                    
                    matplotlib.pyplot.close('all')            

                    minper = period_ranges[0]
                    maxper = period_ranges[1]
                    per = 1.0/f
                    mic_coh = Cxy[(minper <= per) & (per <= maxper)]
                    micper = per[(minper <= per) & (per <= maxper)]
                    coh_mean = np.mean(mic_coh)
                    coh_std = np.std(mic_coh)            

                    coh_means.append(coh_mean)
                    coh_stds.append(coh_std)            

            grid_plots = figure()    
            titlelegend = 'BH' + str(direction)+ ' Array ' + str(array_number) +' ' + str(int(period_ranges[0])) + ' s '+ 'to ' + str(int(period_ranges[1])) + ' s'
            title(titlelegend,fontsize=24)
            my_matrix = np.matrix([coh_means[0:8*1:1],coh_means[8*1:8*2:1],coh_means[8*2:8*3:1],coh_means[8*3:8*4:1],coh_means[8*4:8*5:1],coh_means[8*5:8*6:1],coh_means[8*6:8*7:1],coh_means[8*7:8*8:1]])        
        

            ax = grid_plots.add_subplot(1,1,1)        

            names = range(1, 9)
            heat= sns.heatmap(my_matrix, vmin=0, vmax=1,cmap = pylab.cm.magma,linewidths=.5,annot = True, xticklabels= names, yticklabels = names)
            ylabel('Sensor Number',fontsize=22)    
            xlabel('Sensor Number',fontsize=22)
            pylab.yticks(rotation=0)        
            savefig("Day"+str(day)+'_GRID_Plot_BH'+ str(direction)+'_' + '_from_' + (str(start_time)[:-11]).replace('/',':') + '_to_'+ (str(end_time)[:-11]).replace('/',':') +'_'+ str(period_ranges[0]) + 's' + str(period_ranges[1]) + 's'+'.eps',format='eps',dpi=1000,bbox_inches='tight')
        
        

            matplotlib.pyplot.close('all')        
        
        
        













