#!/usr/bin/env python
###############################################################
# PSD plots for Max Rohde's Long-Period Tilt Paper
# Written by Max Rohde and Adam Ringler


import sys

import numpy
import numpy as np
from scipy import signal
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import pazToFreqResp
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from matplotlib.mlab import csd
from math import pi
import matplotlib.pyplot as plt
import pylab
import matplotlib as mpl
direction_list = [0,1,2]


mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=14)



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

day_list = [
    (170,1466208000,1466244000,1466269200),

    (186,1467590400,1467597600,1467615600),

    (192,1468108800,1468130400,1468152000),

    (180,1467072000,1467079200,1467104400),
]

direction_list = [0,1,2]
NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()
color_list = ['r','b','g','c','m','y','k','chartreuse']
debug= True

sensorType = ["STS-2LG_A", "STS-2HG_B", "STS-2LG_A", "STS-2HG_B", "STS-2HG_A", "STS-2LG_B" ,"STS-2HG_A", "STS-2LG_B"] 
stas = ['XX_FBA1', 'XX_ENG4', 'XX_ENG7', 'XX_ENG5']
locs = ['00','10']

length= 80000        
overlap= 0.5
  
for pidx, day_entry in enumerate(day_list):
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

        instresp =[]
        delta = st[0].stats.delta
        
        length= 80000

        overlap= 0.5
        for idx in range(0,8):
            pazval = getpaz(sensorType[idx])
            instresp.append(computeresp(pazval,delta,length))

    

        st = choptocommon(st)
        fs = st[0].stats.sampling_rate


        for idx in range(0,8):
            #Different sensors compare the power
            (p, f) = cp(st[idx],st[idx],length,overlap,delta)
            psd = 10.*numpy.log10(((2*pi*f)**2)*p/instresp[idx])
            per = 1.0/f
            if 'meanPow' not in vars():
                meanPow = psd
            else:
                meanPow = np.vstack((psd, meanPow))
                    
        # We now have an array of means and stds for the power             
        stdPow = np.std(meanPow,axis=0)
        meanPow = np.mean(meanPow,axis=0)
        ############################################# First part
        fig=plt.figure(1)
        if direction == 0:
            ax=plt.subplot(3,1,1)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                    box.width, box.height * 0.8])
        elif direction == 1:
            ax=plt.subplot(3,1,2)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                box.width, box.height * 0.8])
        elif direction == 2:
            ax=plt.subplot(3,1,3)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                box.width, box.height * 0.8])
            plt.xlabel('Period (s)')
               
        ax.plot(1.0/f, meanPow,linewidth=0.8, label='Array ' + str(pidx+1))
        
        ax.set_xscale('log')
        plt.text(12,-185,'BH' + str(direction),fontsize=12)
        #plt.xticks([])
        if pidx == 3:
            ax.plot(NLNMper,NLNMpower,'k', label='NLNM', linewidth=2.0)
            ax.plot(NHNMper,NHNMpower,'k', linewidth=2.0)
        
        plt.text(12,-185,'BH' + str(direction),fontsize=12)
        plt.xlim((10,1000))
        plt.ylim((-190.,-150.0))
        plt.yticks([-185., -170., -155.])

            
        
            #plt.xticks([10., 100.])
            #plt.ylabel()
    
        fig.text(0.03, 0.5, 'Power  (dB rel. 1 $(m/s^2)^2/Hz$)', ha='center', va='center', rotation='vertical')
        del meanPow
        ############################################## Second part
        fig=plt.figure(2)
        if direction == 0:
            ax=plt.subplot(3,1,1)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                    box.width, box.height * 0.8])
        elif direction == 1:
            ax=plt.subplot(3,1,2)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                box.width, box.height * 0.8])
        elif direction == 2:
            ax=plt.subplot(3,1,3)
            if pidx == 3:  
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height * 0.3,
                box.width, box.height * 0.8])
            plt.xlabel('Period (s)')
               
        ax.plot(1.0/f, stdPow,linewidth=0.8, label='Array ' + str(pidx+1))
        
        ax.set_xscale('log')
        plt.text(12,-185,'BH' + str(direction),fontsize=12)
        #plt.xticks([])
        
        plt.text(12,5,'BH' + str(direction),fontsize=12)
        plt.xlim((10,1000))
        plt.ylim((0,8))
        plt.yticks([0.,4.,8.])

    
        fig.text(0.03, 0.5, 'Power  (dB rel. 1 $(m/s^2)^2/Hz$)', ha='center', va='center', rotation='vertical')
  
fig = plt.figure(1)        
ax=plt.subplot(3,1,3)        
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5),
        fancybox=False, shadow=False, ncol=5,fontsize=12) 
        
plt.xticks([10., 100., 1000.])                      
plt.savefig('Power.eps',format='eps',dpi=1000)

fig = plt.figure(2)
ax=plt.subplot(3,1,3)        
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.5),
        fancybox=False, shadow=False, ncol=5,fontsize=12) 
        
plt.xticks([10., 100., 1000.])                       
plt.savefig('STDPower'+'.eps',format='eps',dpi=1000)

