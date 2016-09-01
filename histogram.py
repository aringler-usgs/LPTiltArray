#!/usr/bin/env python
import argparse
import sys
import math
import pickle
import numpy
import numpy as np
from scipy import signal
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show,xscale, xlim, ylim,scatter)
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import pazToFreqResp, PPSD
from obspy.signal.spectral_estimation import get_NLNM, get_NHNM
from matplotlib.mlab import csd
from math import pi
#from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show,xscale, xlim, ylim,scatter)
import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
import pandas as pd
from scipy import stats, integrate

import seaborn as sns
sns.set()




"""
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=12)"""


mpl.rc('font',weight='bold') 
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from itertools import combinations

direction_list=[0,1,2]


def getpaz(sensor):
    #Function to get the poles and zeros
    #Input the sensor type
    #Output is a paz object
    debugpaz = False
    if debugpaz:
        print 'Sensor type: ' + sensor
    if sensor == 'STS-2HG_A':
        paz = {'gain': 5.96806*10**7, 'zeros': [0, 0], 'poles': [-0.035647 - 0.036879j,  
            -0.035647 + 0.036879j, -251.33, -131.04 - 467.29j, -131.04 + 467.29j],
            'sensitivity': 3.3554432*10**10}
    elif sensor == 'STS-2HG_B':
        paz = {'gain': 5.96806*10**7, 'zeros': [0, 0], 'poles': [-0.035647 - 0.036879j,  
            -0.035647 + 0.036879j, -251.33, -131.04 - 467.29j, -131.04 + 467.29j],
            'sensitivity': 8.388608*10**9}
    elif sensor == 'STS-2LG_A':
        paz = {'gain': 5.96806*10**7, 'zeros': [0, 0], 'poles': [-0.035647 - 0.036879j,  
            -0.035647 + 0.036879j, -251.33, -131.04 - 467.29j, -131.04 + 467.29j],
            'sensitivity': 2.5165824*10**9}
    elif sensor == 'STS-2LG_B':
        paz = {'gain': 5.96806*10**7, 'zeros': [0, 0], 'poles': [-0.035647 - 0.036879j,  
            -0.035647 + 0.036879j, -251.33, -131.04 - 467.29j, -131.04 + 467.29j],
            'sensitivity': 6.291456*10**8}
    if debugpaz:
        print(paz)
    return paz        

def computeresp(resp,delta,lenfft):
    respval = pazToFreqResp(resp['poles'],resp['zeros'],resp['sensitivity']*resp['gain'],t_samp = delta, 
    nfft=lenfft,freq = False)
    respval = numpy.absolute(respval*numpy.conjugate(respval))
    respval = respval[1:]
    return respval        

def cp(tr1,tr2,lenfft,lenol,delta):
    sr = 1/delta
    cpval,fre = csd(tr1.data,tr2.data,NFFT=lenfft,Fs=sr,noverlap=lenol,scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre        

def choptocommon(stream):
    #A function to chop the data to a common time window
    debug = False
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



day_list = [
    (170,1466208000,1466244000,1466269200),

    (186,1467590400,1467597600,1467615600),

    (192,1468108800,1468130400,1468152000),

    (180,1467072000,1467079200,1467104400),
]
min_max_array= [(20.0,30.0),(30.0,50.0),(50.0,100.0),(100.0,500.0)]  



for period_ranges in min_max_array:

    minper = period_ranges[0]
    maxper = period_ranges[1]

    for direction in direction_list: 


        meanCoh =[]
        sensorType = ["STS-2LG_A", "STS-2HG_B", "STS-2LG_A", "STS-2HG_B", "STS-2HG_A", "STS-2LG_B" ,"STS-2HG_A", "STS-2LG_B"] 
        stas = ['XX_FBA1', 'XX_ENG4', 'XX_ENG7', 'XX_ENG5']
        locs = ['00','10']

        length= 80000        

        overlap= 0.5

        for pidx, day_entry in enumerate(day_list):

            ###############################################################################
            # PARAMETERS
            day = day_entry[0]
            trim = True
            start_time = UTCDateTime(day_entry[2])
            end_time = UTCDateTime(day_entry[3])
       

          
            st = Stream()
            for sta in stas:
                for loc in locs:
                    st += read("/tr1/telemetry_days/" + sta + 
                    "/2016/2016_" + str(day) + "/" + loc + "_BH" + str(direction) + ".512.seed")


            if trim:
                st.trim(start_time,end_time)        
            st = choptocommon(st)
            fs = st[0].stats.sampling_rate
            print(st)


            for comb in combinations(range(0,8),2):
                idx = comb[0]
                idx2 = comb[1]
                #Different sensors compare the coherence
                f, Cxy = signal.coherence(st[idx].data, st[idx2].data, fs, nperseg=length,noverlap=overlap)
                per = 1.0/f
                Cxy = Cxy[(minper <= per) & (per <= maxper)]
                meanCoh.append(np.mean(Cxy))
                #if 'meanCoh' not in vars():
                #    meanCoh = Cxy
                #else:
                #    meanCoh = np.vstack((Cxy, meanCoh))







        fig=plt.figure(1)

        titlelegend = 'BH' + str(direction) + ' ' + str(int(minper)) + ' s '+ 'to ' + str(int(maxper)) + ' s'
        title(titlelegend,fontsize=44)

        ax = sns.distplot(meanCoh,bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],kde=False,rug=True,color='lightslategrey')

        ax.set(xlim=(0, 1))
        ax.set(ylim=(0, 112))
        xticklabels = [0.0,0.2,0.4,0.6,0.8,1.0]
        yticklabels = [0.0,20.0,40.0,60.0,80.0,100.0]
        ax.set_xticklabels(xticklabels,fontsize=22)
        ax.set_yticklabels(yticklabels,fontsize=22)
        ax.tick_params(axis='both', which='major', pad=15)

 
        plt.savefig('Histogram_BH'+ str(direction)+'_' + str(minper) + '_' + str(maxper) +'.svg',format='svg',dpi=1000) 
        mpl.pyplot.close('all')

     
     
                    
                        
            ## We now have an array of means and stds for the coherences              
            #stdCoh = np.std(meanCoh,axis=0)
            #meanCoh = np.mean(meanCoh,axis=0)
            #fig=plt.figure(1)
            #if pidx > 0 and direction == 0:
                #ax=plt.subplot(4,1,pidx+1,sharex=ax)
            #elif direction == 0:
                #ax=plt.subplot(4,1,pidx+1)
            #ax.errorbar(1.0/f, meanCoh,stdCoh,linewidth=0.8, label='BH' + str(direction))
            #ax.set_xscale('log')
            #plt.xticks([])
            #plt.xlim((10,500))
            #plt.ylim((0.0,1.0))
            #plt.text(12,.5,'Array '+ str(pidx+1),fontsize=12)
            #plt.yticks([.2, .4, .6, .8])
            
            #if pidx == 3:
                #plt.xlabel('Period (s)')
                #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5),
                    #fancybox=False, shadow=False, ncol=3)
                ##plt.xticks([10., 100.])
                ##plt.ylabel()
            #box = ax.get_position()
            #ax.set_position([box.x0, box.y0 + box.height * 0.05,
                     #box.width, box.height * 0.95])
            #fig.text(0.05, 0.5, 'Coherence ($\gamma^2$)', ha='center', va='center', rotation='vertical')

            #del meanCoh
    #plt.xticks([10., 100.])        
    ##plt.tight_layout()                
    #plt.savefig('Coheren.jpg',format='jpeg')
