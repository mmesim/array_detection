# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:13:22 2016

@author: lisa linville
Utilities for LTX detections
"""
import numpy as np
import datetime
from math import pi, cos, radians
import geopy.distance as pydist
from numpy import median, absolute
from scipy.signal import spectrogram

#ANF event catalog parser stolen from detex: github.com/d-chambers/detex
import pandas as pd
import glob
import os
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def readANF(anfdir,lon1=-180,lon2=180,lat1=0,lat2=90,getPhases=False,UTC1='1960-01-01',
            UTC2='3000-01-01',Pcodes=['P','Pg'],Scodes=['S','Sg']):
    """Function to read the ANF directories as downloaded from the ANF Earthscope Website"""
    monthDirs=glob.glob(os.path.join(anfdir,'*'))
    Eve=pd.DataFrame()
    for month in monthDirs:
        utc1=UTCDateTime(UTC1).timestamp
        utc2=UTCDateTime(UTC2).timestamp
        #read files for each month
        
        
        dfOrigin=_readOrigin(glob.glob(os.path.join(month,'*.origin'))[0])
        dfOrigerr=readOrigerr(glob.glob(os.path.join(month,'*.origerr'))[0])
                #merge event files togther
        DF=pd.merge(dfOrigin,dfOrigerr)
            
        #discard all events outside area of interest
        DF=DF[(DF.Lat>lat1)&(DF.Lat<lat2)&(DF.Lon>lon1)&(DF.Lon<lon2)&(DF.time>utc1)&(DF.time<utc2)]
        
        if getPhases:
            dfAssoc=_readAssoc(glob.glob(os.path.join(month,'*.assoc'))[0])
            dfArrival=_readArrival(glob.glob(os.path.join(month,'*.arrival'))[0])
 
            #link associated phases with files
            DF=_linkPhases(DF,dfAssoc,dfArrival,Pcodes,Scodes)
        
        Eve=pd.concat([DF,Eve],ignore_index=True)
        Eve.reset_index(drop=True,inplace=True)
    return Eve

def readOrigerr(origerrFile):
    columnNames=['orid','sobs','smajax','sminax','strike','sdepth','conf']
    columnSpecs=[(0,8),(169,179),(179,188),(189,198),(199,205),(206,215),(225,230)]
    df=pd.read_fwf(origerrFile,colspecs=columnSpecs,header=None,names=columnNames)
    return df

def _readOrigin(originFile):
    columnNames=['Lat','Lon','depth','time','orid','evid',
     'jdate','nass','ndef','ndp','grn','srn','etype','review','depdp','dtype',
     'mb','mbid','ms','msid','ml','mlid','algo','auth','commid','lddate']
    columnSpecs=[(0,9),(10,20),(20,29),(30,47),(48,56),(57,65),(66,74),(75,79),(80,84),
                (85,89),(90,98),(99,107),(108,110),(111,115),(116,125),(126,128),(128,136),
                (136,144),(145,152),(153,161),(162,169),(170,178),(179,194),(195,210),
                (211,219),(220,237)]
    
    df=pd.read_fwf(originFile,colspecs=columnSpecs,header=None,names=columnNames)
    df['DateString']=[UTCDateTime(x).format_iris_web_service() for x in df.time]
    return df
    
def _readAssoc(assocFile):
    columnNames=['arid','orid','sta','phase','belief','delta']
    columnSpecs=[(0,8),(9,17),(18,24),(25,33),(34,38),(39,47)]
    df=pd.read_fwf(assocFile,colspecs=columnSpecs,header=None,names=columnNames)
    return df

def _readArrival(arrivalFile):
    columnNames=['sta','time','arid','stassid','iphase','amp','per','snr']
    columnSpecs=[(0,6),(7,24),(25,33),(43,51),(70,78),(136,146),(147,154),(168,178)]
    df=pd.read_fwf(arrivalFile,colspecs=columnSpecs,header=None,names=columnNames)
    return df
    
def _linkPhases(DF,dfAssoc,dfArrival,Pcodes,Scodes):
    DF['Picks']=[{} for x in range(len(DF))]
    for a in DF.iterrows():
        dfas=dfAssoc[dfAssoc.orid==a[1].orid] #DF associated with orid, should be one row
        dfas=dfas[dfas.phase.isin(Pcodes+Scodes)]
        dfas['time']=float()
        dfas['snr']=float
        for b in dfas.iterrows(): #get times from df arrival
            dfar=dfArrival[dfArrival.arid==b[1].arid]
            dfas.time[b[0]]=dfar.time.iloc[0]
            dfas.snr[b[0]]=dfar.snr.iloc[0]
        for sta in list(set(dfas.sta.values)):
            dfasSta=dfas[dfas.sta==sta]
            dfasP=dfasSta[dfasSta.phase.isin(Pcodes)]
            dfasS=dfasSta[dfasSta.phase.isin(Scodes)]
            tempdict={sta:[0,0]}
            if len(dfasP)>0:
                tempdict[sta][0]=dfasP.time.iloc[0]
            if len(dfasS)>0:
                tempdict[sta][1]=dfasS.time.iloc[0]
            DF.Picks[a[0]]=dict(DF.Picks[a[0]].items()+tempdict.items())
    return DF
    
def ANFtoTemplateKey(anfDF,temKeyName='TemplateKey_anf.csv',saveTempKey=True):   
    """Convert the dataframe created by the readANF function to a detex templatekey csv"""
    ds=[x.split('.')[0].replace(':','-') for x in anfDF.DateString]
    ts=[x.replace(':','-') for x in anfDF.DateString]
    contrib=['ANF']*len(anfDF)
    mtype=['ML']*len(anfDF)
    stakey=['StationKey.csv']*len(anfDF)
    df=pd.DataFrame()
    df['CONTRIBUTOR'],df['NAME'],df['TIME'],df['LAT'],df['LON'],df['DEPTH'],df['MTYPE'],df['MAG'],df['STATIONKEY']=contrib,ds,ts,anfDF.Lat.tolist(),anfDF.Lon.tolist(),anfDF.depth.tolist(),mtype,anfDF.ml.tolist(),stakey
    if saveTempKey:
        df.to_csv(temKeyName)
    return df
#**********************************
##functions for computing polygon area and running mean etc    

def mad(data, axis=None):
    '''use numpy to calculate median absolute deviation (MAD), more robust than std'''
    return median(absolute(data - median(data, axis)), axis)

def get_levels(rays, madtimes):
    """Set detection threshold based off 4 times MAD of input data"""
    a,b = np.shape(rays)
    temp = np.reshape(rays, [1,a*b])
    return mad(temp)*madtimes # orignal value is 3.00
    
def get_saturation_value(rays, masktimes):
    a,b = np.shape(rays)
    temp = np.reshape(rays, [1,a*b])
    return mad(temp)*masktimes # mask
    
#clean up the array at saturation value just below the detection threshold
def saturateArray(array, masktimes):
    """saturate array at 6*MAD """
    shigh = get_saturation_value(array, masktimes)
    sloww = shigh*-1
    junk=[]    
    for i in range(np.shape(array)[0]):
        fill = np.sum(array,axis=1)
    if np.sum(array[i][:]) >= 2.5*mad(fill):
        array[i][:]= np.median(array)        
    junk = np.where(array>=shigh)
    array[junk]=shigh
    junk = np.where(array<=sloww)
    array[junk]=sloww
    narray = [array[x][y]/np.max(array[x][:]) for x in range(len(array)) for y in range(len(np.transpose(array)))]
    array = np.reshape(narray,np.shape(array))
#    #if there are too many hi-amp value swaps subdue that channel
#    zcrosst = [(np.tanh(array[x][:-1]) > .5).sum() for x in range(len(array))]
#    junk= np.where(np.array(zcrosst) >= 350)
#    for each in junk[0]:
#        array[each][:] = array[each]/2
    return array
    
def PolygonArea(corners):
    earth_radius = 6371009 # in meters
    lat_dist = pi * earth_radius / 180.0
    x=corners[:,0]
    y=corners[:,1]
    yg = [lat * lat_dist for lat in y]
    xg = [lon * lat_dist * cos(radians(lat)) for lat, lon in zip(y,x)]
    corners=[yg,xg]
    area = 0.0
    for i in range(-1, len(xg)-1):
        area += xg[i] * (yg[i+1] - yg[i-1])
    return abs(area) / 2.0
 
def runningmean(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]
    
    
def get_k(x,y,triangles,ratio):
    """Threshold value for area based on input ratio times
    the average station triad area"""
    out = []
    for points in triangles:
        a,b,c = points
        d0 = pydist.vincenty([x[a],y[a]],[x[b],y[b]]).meters
        d1 = pydist.vincenty([x[b],y[b]],[x[c],y[c]]).meters
        d2 = pydist.vincenty([x[c],y[c]],[x[a],y[a]]).meters
        s=d0+d1+d2/2
        arear=np.sqrt(s*(s-d0)*(s-d1)*(s-d2))
        out.append(arear)
    k_value=np.median(out)*ratio
    return k_value
    
def get_edge_ratio(x,y,triangles,ratio):
    """Mask out long edges from the station mesh. Too stingy and you get holes
    which cause problems (spurrious detections when single station amplitudes 
    span the gap), too large and you preferentially get detections from the 
    long interstation distances."""
    out,outl = [],[]
    for points in triangles:
        a,b,c = points
        d0 = np.sqrt( (x[a] - x[b]) **2 + (y[a] - y[b])**2 )
        d1 = np.sqrt( (x[b] - x[c]) **2 + (y[b] - y[c])**2 )
        d2 = np.sqrt( (x[c] - x[a]) **2 + (y[c] - y[a])**2 )
        out.append(d0)
        out.append(d1)
        out.append(d2)
        d0 = pydist.vincenty([x[a],y[a]],[x[b],y[b]]).meters
        outl.append(d0)
        d1 = pydist.vincenty([x[b],y[b]],[x[c],y[c]]).meters
        outl.append(d1)
        d2 = pydist.vincenty([x[c],y[c]],[x[a],y[a]]).meters
        outl.append(d2)
    mask_length=np.median(out)*ratio
    #need the median edge in meters
    median_edge=np.median(outl)
    return mask_length,median_edge

def long_edges(x, y, triangles, ratio=2.5):
    olen,edgeL=get_edge_ratio(x,y,triangles,ratio)
    out = []
    for points in triangles:
        #print points
        a,b,c = points
        d0 = np.sqrt( (x[a] - x[b]) **2 + (y[a] - y[b])**2 )
        d1 = np.sqrt( (x[b] - x[c]) **2 + (y[b] - y[c])**2 )
        d2 = np.sqrt( (x[c] - x[a]) **2 + (y[c] - y[a])**2 )
        #d0 = pydist.vincenty([x[a],y[a]],[x[b],y[b]]).meters
        #d1 = pydist.vincenty([x[b],y[b]],[x[c],y[c]]).meters
        #d2 = pydist.vincenty([x[c],y[c]],[x[a],y[a]]).meters
        max_edge = max([d0, d1, d2])
        #print points, max_edge
        if max_edge > olen:
            out.append(True)
        else:
            out.append(False)
    return out,edgeL

def templatetimes(detectiontime,tlength,delta):
    vec = detectiontime-datetime.timedelta(seconds = tlength/delta)
    #end = tr.stats.endtime
    step = datetime.timedelta(seconds=1.0/delta)
    out = []
    while len(out)<tlength*2:
        out.append(vec)
        vec += step
    return out 
    
def gettvals(tr1,tr2,tr3):
    mlen= max(len(tr1),len(tr2),len(tr3))
    if len(tr1)==mlen:
        tr=tr1
    elif len(tr2)==mlen:
        tr=tr2
    else:
        tr=tr3
    vec = tr.stats.starttime
    vec=vec.datetime
    #end = tr.stats.endtime
    step = datetime.timedelta(seconds=tr.stats.delta)
    out = []
    while len(out)<len(tr.data):
        out.append(vec)
        vec += step
    return out
    
def getfvals(tt,Bwhite,nseconds,edgebuffer):
    vec = tt.datetime
    ed = tt+(nseconds+edgebuffer)
    step = datetime.timedelta(seconds=((nseconds+edgebuffer)/float(len(np.transpose(Bwhite)))))
    out = []
    while vec <= ed.datetime:
        out.append(vec)
        vec += step
    return out
    

    
##get catalog data (ANF right now only)
def getCatalogData(tt,nseconds,lo,ll):
    import geopy.distance as pydist
    localE= readANF('anfdir',UTC1=tt,UTC2=tt+nseconds, lon1=min(lo),lon2=max(lo),lat1=min(ll),lat2=max(ll))
    globalE= readANF('anfdir',UTC1=tt,UTC2=tt+nseconds)
    #fix for issue 011: remove events from global that overlap with local
    dg = globalE
    for i in range(len(localE)):
        dg = globalE[globalE.DateString != localE.DateString[i]]
        globalE =dg
    globalE = globalE.reset_index(drop=True)
    distarray,closesti=[],[]
    for event in range(len(localE)):
        for each in range(len(ll)):
            distarray.append(pydist.vincenty([localE.Lat[event],localE.Lon[event]],[ll[each],lo[each]]).meters)
        closesti.append(np.argmin(distarray))
        distarray = []
    return localE,globalE,closesti


def bestCentroid(detections,localev,centroids,localE,ctimes):
    import bisect
    from obspy import UTCDateTime
    centroid = np.empty([len(detections),2])
    atimes=[]
    for j in range(len(localE)):
        atimes.append(UTCDateTime(localE.DateString[j]).datetime)
       
    for each in range(len(detections)):
        
        if localev.count(detections[each]) ==1:
            localEi=bisect.bisect_left(atimes, (ctimes[detections[each]]))
            if localEi == 0:
                centroid[each][0]=localE.Lat[localEi]
                centroid[each][1]=localE.Lon[localEi]
            else:
                centroid[each][0]=localE.Lat[localEi-1]
                centroid[each][1]=localE.Lon[localEi-1]
        else:
            centroid[each][0]=centroids[detections[each]][1]
            centroid[each][1]=centroids[detections[each]][0]
        
    return centroid

   
def markType(detections,centroids,localev,localE,ctimes,doubles):
    cents= bestCentroid(detections,localev,centroids,localE,ctimes)    
    temp = np.empty([len(detections)])
    dtype = []
    for event in range(len(detections)):
            dtype.append('earthquake')
    return dtype,cents
    

    
def w_spec(szdata,deltaf,fftsize,freq1,freq2):
    '''return whitened spectrogram in decibles'''
    specgram = spectrogram(szdata,fs=deltaf,nperseg=fftsize,window=('hanning'),scaling='spectrum',noverlap = fftsize/2)
    specgram_cut = [specgram[2][index] for index,  i in enumerate(specgram[0]) if i >= freq1 and i <= freq2]
    sg = 10*np.log10(specgram_cut)
    bgs=[runningmean(sg[count,:],50) for count in range(len(sg))]
    endbgs = [np.median(bgs[count][-101:-50]) for count in range(len(bgs))]
    begbgs = [np.median(bgs[count][52:103]) for count in range(len(bgs))]
    tbags = bgs
    for i in range(len(bgs)):
        for j in range(-1,-51,-1):
            tbags[i][j] = endbgs[i]
        for k in range(1,51,1):
            tbags[i][k] = begbgs[i]
    Bwhite=sg-tbags 
    return Bwhite
    
def spec(szdata,deltaf,fftsize):
    '''return spectrogram in decibles'''
    specgram = spectrogram(szdata,fs=deltaf,nperseg=fftsize,window=('hanning'),scaling='spectrum')
    sg = 10*np.log10(specgram[2]) 
    return sg

def reviewer(filestring='2010_*'):
#    import pandas as pd
#    import glob
#    import os
#    import matplotlib.pyplot as plt
#    import matplotlib.image as mpimg
    plt.rcParams['figure.figsize'] = 18,12 #width,then height
    dlist=sorted(glob.glob(filestring), key=os.path.getmtime)
    for eachdir in dlist:
        os.chdir(eachdir)
        if os.path.isfile('rptable.pkl'):
            os.chdir('../')
        else:
            try:
                df = pd.read_pickle('picktable.pkl')
                df = df[df.Type !='blast']
                df = df.reset_index(drop=True)
                #lets you put in confidence as an array, 1 value for each station
                df = df.astype(object)
                imlist=sorted(glob.glob('image*.eps'))
                if len(df) == len(imlist):
                #check df eq len and number of im's are equal
                    count = 0
                    for im in imlist:
                        img = mpimg.imread(im)
                        plt.imshow(img)
                        plt.show()
                        conf = input()
                        df.Confidence[count]=conf
                        count=count+1
                    df.to_pickle('rptable.pkl')    
                    os.chdir('../')
                else:
                    print('length of pick table and number of images are not the same')
                    
                    
                    os.chdir('../')        
            except:
                print('no picktable for '+str(eachdir))
            
                os.chdir('../')

def cat_df(filestring='2010_*'):
    
    dlist=sorted(glob.glob(filestring), key=os.path.getmtime)
    os.chdir(dlist[0])
    df1=pd.read_pickle('rptable.pkl')
    df1 = df1[df1.Confidence != -1]
    os.chdir('../')
    dlist = dlist[1:]
    
    for alldirs in range(len(dlist)):
        os.chdir(dlist[alldirs])
        try:
            df=pd.read_pickle('rptable.pkl')
            df = df[df.Confidence != -1]
            df1 = pd.concat([df1,df])
            os.chdir('../')
        except:
            print('no picktable for '+dlist[alldirs])
            os.chrid('../')
    df1.sort_values(by='S1time', inplace=True)
    df = df1.reset_index(drop=True)
    
    df.to_html('reviewed_picks.html') 
