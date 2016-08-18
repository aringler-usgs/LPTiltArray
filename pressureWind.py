#!/usr/bin/env python

###############################################################
# Wind and Pressure for Max Rohde's Long-Period Tilt Paper
# Written by Max Rohde and Adam Ringler

import numpy
import numpy as np
from scipy import signal
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import pazToFreqResp, PPSD

from matplotlib.mlab import csd
import matplotlib.pyplot as plt
import matplotlib as mpl
direction_list = [0,1,2]
from scipy.signal import detrend

mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

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
    respval = np.absolute(respval*np.conjugate(respval))
    respval = respval[1:]
    return respval        

def cp(tr1,tr2,lenfft,lenol,delta):
    sr = 1/delta
    cpval,fre = csd(tr1.data, tr2.data, NFFT=lenfft, Fs=sr, noverlap=lenol, scale_by_freq=True)
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




curdis = []
meanCoh =[]
sensorType = ["STS-2LG_A", "STS-2HG_B", "STS-2LG_A", "STS-2HG_B", "STS-2HG_A", "STS-2LG_B" ,"STS-2HG_A", "STS-2LG_B"] 
stas = ['XX_FBA1', 'XX_ENG4', 'XX_ENG7', 'XX_ENG5']
locs = ['00','10']

length= 80000        
overlap= 0.5
trim = True
for pidx, day_entry in enumerate(day_list):

    day = day_entry[0]

    start_time = UTCDateTime(day_entry[2])
    end_time = UTCDateTime(day_entry[3])        

      
    st1 = Stream()
    st2 = Stream()
    for sta in stas:
        for loc in locs:
            st1 += read("/tr1/telemetry_days/" + sta + 
                "/2016/2016_" + str(day) + "/" + loc + "_BH1.512.seed")
            st2 += read("/tr1/telemetry_days/" + sta + 
                "/2016/2016_" + str(day) + "/" + loc + "_BH2.512.seed")

    if trim:
        st1.trim(start_time,end_time)
        st2.trim(start_time,end_time)       
    st1 = choptocommon(st1)
    st2 = choptocommon(st2)
    fs = st1[0].stats.sampling_rate

    for sidx, tr in enumerate(st1):
        tr.detrend()
        tr.taper(.05)
        tr.simulate(paz_remove=getpaz(sensorType[sidx]))

    for sidx, tr in enumerate(st2):
        tr.detrend()
        tr.taper(.05)
        tr.simulate(paz_remove=getpaz(sensorType[sidx]))

    # Do the estimates for the period bands
    bands = [(20.,30.),(30.,50.),(50.,100.),(100., 500.)]
    for nidx, band in enumerate(bands):
        
        st1temp = st1.copy()
        st2temp = st2.copy()
        st1temp.filter('bandpass', freqmin=1./band[1], freqmax=1./band[0])
        st2temp.filter('bandpass', freqmin=1./band[1], freqmax=1./band[0])
        st1temp.trim(starttime=st1temp[0].stats.starttime+3.*60.*60.)
        st2temp.trim(starttime=st2temp[0].stats.starttime+3.*60.*60.)

    
        print(st1temp)
        for idx in range(0,8):
            if 'meanData' not in vars():
                meanData = np.sqrt(st1temp[idx].data**2 + st2temp[idx].data**2)
            else:
                meanData = np.vstack((np.sqrt(st1temp[idx].data**2 + st2temp[idx].data**2), meanData))
        
        stdData = np.std(meanData,axis=0)
        meanData = np.mean(meanData,axis=0) 
       
        t1 = np.arange(0, st1temp[0].stats.npts / st1temp[0].stats.sampling_rate, st1temp[0].stats.delta)/(60.*60.)
        fig = plt.figure(1,figsize=(12,12))
        plt.subplots_adjust(hspace=0.001)
        if nidx == 0:
            ax=plt.subplot(6,1,nidx+1)
            ax.plot(t1,meanData*10**9 + stdData*10**9,color='.5')
            ax.plot(t1,meanData*10**9,color='k')

        else:
            ax3 = plt.subplot(6,1,nidx+1,sharex=ax)
            ax3.plot(t1,meanData*10**9 + stdData*10**9,color='.5')
            ax3.plot(t1,meanData*10**9,color='k')
        del meanData
        del stdData
        plt.xticks([])

        plt.xlim((0,max(t1)))
        if nidx == 0:
            plt.ylim((0,15))
            plt.yticks([2,7.5,13])
            plt.text(0.2,10,'Band-Pass ' + str(int(band[0])) + ' s to ' + str(int(band[1])) + ' s',fontsize=18)
        elif nidx == 1:
            plt.ylim((0,15))
            plt.yticks([2,7.5,13])
            plt.text(0.2,10,'Band-Pass ' + str(int(band[0])) + ' s to ' + str(int(band[1])) + ' s',fontsize=18)
        elif nidx == 2:
            plt.ylim((0,25))
            plt.yticks([3,12.5,22])
            plt.text(0.2,19,'Band-Pass ' + str(int(band[0])) + ' s to ' + str(int(band[1])) + ' s',fontsize=18)
            yyl = plt.ylabel('Velocity  $(nm/s)$',fontsize=18, labelpad=20)
            yyl.set_position((yyl.get_position()[0],1)) 
            yyl.set_verticalalignment('center')
        elif nidx == 3:
            plt.ylim((0,100))
            plt.yticks([10, 50, 90])
            plt.text(0.2,80,'Band-Pass ' + str(int(band[0])) + ' s to ' + str(int(band[1])) + ' s',fontsize=18)


    # Grab the environment data
    wind_data = "/tr1/telemetry_days/IU_ANMO/2016/2016_" + str(day) + "/50_LWS.512.seed"        
    pressure_data = "/tr1/telemetry_days/IU_ANMO/2016/2016_" + str(day) + "/31_LDO.512.seed"        
    st = Stream()
    st +=read(wind_data)
    st[0].data = np.require(st[0].data, dtype=np.float64)
    # In m/s
    st[0].data /= 10.
    st +=read(pressure_data)
    st[1].data = np.require(st[1].data, dtype=np.float64)
    #In Passcals
    st[1].data = 6.0*10**4 + 0.01227*st[1].data
    st = choptocommon(st)
    st.trim(start_time,end_time)
    st.trim(starttime=st[0].stats.starttime+3.*60.*60.)
    t2 = np.arange(0, st[0].stats.npts / st[0].stats.sampling_rate, st[0].stats.delta)/(60.*60.)
    ax4 = plt.subplot(6,1,5)
    ax4.plot(t2,st[0].data,color='k')
    plt.xlim((0,max(t2)))
    plt.ylim((0,10))
    plt.yticks([1, 5, 9])
    plt.ylabel('Speed $(m/s)$', labelpad=15,fontsize=18)
    xticklabels = ax4.get_xticklabels() 
    plt.setp(xticklabels, visible=False)
    ax2 = plt.subplot(6,1,6)
    plt.plot(t2,st[1].data,color='k') 
    plt.xlim((0,max(t2)))
    plt.ylabel('Press. (kPa)', labelpad=10,fontsize=18)
    plt.xlabel('Time (hours)',fontsize=18)

    plt.savefig('PressWind' + str(pidx) + '.eps',format='eps',dpi=1000)
    plt.clf()

