##!/uufs/chpc.utah.edu/common/home/u0149900/.conda/envs/anaconda_2.7.11/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:08:15 2016

@authors: linville.seis.utah.edu; dkilb@ucsd.edu

Scipt to generate array images in specified year_day directories.
From: Linville, L., K. Pankow, D. Kilb, and A. Velasco (2014), 
doi:10.1002/2014JB011529.

Parameters
------------------
wb : int, which basin number to import station/blast site lists from

ndays: int, how many days to process, in 2 hour blocks

thresholdv: float, what value above the avg area based on station space

levels: which contours to generate and base detections off

"""
def freqtor(yr, mo, dy, hr, mn, sc, duration, ndays, datdir, freq1, freq2,
            thresholdv, deltaf, masktimes, madtimes, time_thres, distance_thres):

#control plot behavior

	import matplotlib.pylab as plt
	plt.switch_backend("nbagg")
	plt.style.use('ggplot')
	plt.rcParams['figure.figsize'] = 18,12 #width,then height
	plt.rcParams['savefig.dpi'] = 80
	from obspy import UTCDateTime
	import numpy as np
	import matplotlib.dates as mdates 
	import matplotlib.tri as tri
	from obspy.signal.trigger import recursive_sta_lta as recSTALTA
	from obspy.signal.trigger import trigger_onset as triggerOnset
	import copy,os,bisect,scipy,datetime,itertools
	import pandas as pd
	#suppress the chained assignment warning
	pd.options.mode.chained_assignment = None
	from mpl_toolkits.basemap import Basemap
	from obspy.taup import TauPyModel as TauP
	model = TauP(model="iasp91")
	from obspy.geodetics import locations2degrees as loc2d
	import Util as Ut
	import geopy.distance as pydist
	from obspy.core import read
	#############################
	homedir=''
	wb = 5 #which basin # are we working on for station list import
	maketemplates = 1
	tlength = 7200 #nsamples on either side of detection time for template
	counter = datetime.date(int(yr),int(mo),int(dy)).timetuple().tm_yday
	edgebuffer = 00
	duration = duration +edgebuffer
	#ndays= 2 #however many days you want to generate images for
	dayat = int(dy)
	#set parameter values; k = area threshold for detections:
    #thresholdv= 2.0
    #deltaf = 250.0
	nseconds = 7200
	npts = int(deltaf*(nseconds+edgebuffer))
	fftsize=256
	overlap=4   
	hop = fftsize / overlap
	w = scipy.hanning(fftsize+1)[:-1] 
    #delta=250.0

	if duration == 86400:
		im = 12
	elif duration == 7200:
		im=1
	else:
		im=1

	#parse the datetime
	counter_3char = str(counter).zfill(3)
	datest =yr+str('-')+mo+str('-')+str(dayat)+str('T')+hr+str(':')+mn+str('.')+sc 
	tt = UTCDateTime(datest)

	#####################################################################
	# Now start making the detections, in 2 hour data chunks, 1 day at a time
	print (os.getcwd())
	for days in range(ndays):
		plt.close('all')
		print(str(tt))
		sacyear=str(tt.date.year)
		sacmonth= str(tt.date.month)
		sacday = str(tt.date.day)
		if len(sacmonth)==1:
			sacmonth=str(0)+sacmonth
		if len(sacday)==1:
			sacday = str(0)+sacday
		sacname = str(sacyear)+str(sacmonth)+str(sacday)
		sacdir = datdir+sacyear+sacmonth+sacday+'/'+sacname+'*.sac'
		#############################
	
		s = homedir+'basin%s/'%wb+yr+str('_')+counter_3char
		if not os.path.exists(s):
			os.makedirs(s)     
		sz = read(sacdir)
		sz.sort()
		sz.detrend()
		sz.trim(starttime=tt,endtime=tt+duration,pad=True,fill_value=000,nearest_sample=False)
		sz.filter('highpass',freq=1.0)

	
		alltimes=Ut.gettvals(sz[0],sz[1],sz[2])
		#############################
		#########################
		#%%
		nptsf = edgebuffer*deltaf
		blockette = 0
		d = {'Contributor': 'NA', 'Latitude': 'NA','Longitude': 'NA', 'S1': 'NA',
			 'S1time': 'NA', 'Magnitude': -999.00, 'mag_error': -999.00,'cent_er': -999.00,
			 'Confidence': 0,'S2':'NA','S3':'NA','S4': 'NA', 'S5': 'NA',
			 'S2time': 'NA','S3time': 'NA','S4time': 'NA','S5time': 'NA',
			 'Type': 'Event'}
		index = [0]; df1 = pd.DataFrame(data=d, index=index)   
		stations,latitudes,longitudes,distances=[],[],[],[]
		snames,latitudes,longitudes = [],[],[]
		for i in range(len(sz)):
			snames.append(str(sz[i].stats.station))
			latitudes.append(sz[i].stats.sac['stla'])
			longitudes.append(sz[i].stats.sac['stlo'])
		latmin=min(latitudes);lonmin=max(longitudes) 
		newlat= np.empty([len(snames)])
		newlon= np.empty([len(snames)])
		stations = copy.copy(snames)
		for i in range(len(snames)):
			reindex = stations.index(snames[i])
			newlat[i]=latitudes[reindex]
			newlon[i]=longitudes[reindex]
			distances.append(pydist.vincenty([newlat[i],newlon[i]],[latmin,lonmin]).meters)
		#####this is where maths happends and arrays are created
		for block in range(im):
			print(blockette,tt)
			ll,lo,stalist,vizray,dist=[],[],[],[],[]
			shorty = 0
			for z in range(len(snames)):
				szdata = sz[z].data[blockette:blockette+npts]   
			   # if len(szdata)==npts:
				vizray.append([])
				Bwhite = Ut.w_spec(szdata,deltaf,fftsize,freq1,freq2)
				vizray[shorty].append(np.sum(Bwhite[:,:],axis=0))
				ll.append(newlat[z])
				lo.append(newlon[z])
				dist.append(distances[z])
				stalist.append(snames[z])
				shorty=shorty+1
			rays = np.vstack(np.array(vizray))
			ix = np.where(np.isnan(rays))
			rays[ix] =0
			rayz=np.copy(rays)
			latudes=copy.copy(ll)
			longitudes=copy.copy(lo)
			slist=copy.copy(stalist)
			#sort the array orders by distance from lomin,latmin
			for i in range(len(slist)):
				junk=np.where(np.array(dist)==max(dist))
				rayz[i]=rays[junk[0][0]]
				ll[i]=latudes[junk[0][0]]
				lo[i]=longitudes[junk[0][0]]
				slist[i]=stalist[junk[0][0]]
				dist[junk[0][0]]=-9999999999
			timevector = Ut.getfvals(tt,Bwhite,nseconds,edgebuffer)

			#clean up the array 
			rayz = Ut.saturateArray(rayz, masktimes)
			ix = np.where(np.isnan(rayz))
			rayz[ix] =0
	   
			#determine which level to use as detections 4* MAD
			levels=[Ut.get_levels(rayz, madtimes)]
	   
			#get the ANF catalog events and get closest station
		
			localE,globalE,closesti=Ut.getCatalogData(tt,nseconds,lo,ll)


			#closesti = np.flipud(closesti) 
			#unstructured triangular mesh with stations as verticies, mask out the long edges
			triang = tri.Triangulation(lo, ll)
			mask,edgeL = Ut.long_edges(lo,ll, triang.triangles)
			triang.set_mask(mask)
			kval=Ut.get_k(lo,ll,triang.triangles,thresholdv)
		
	#%%
			#get contour areas by frame
			av,aa,xc,yc,centroids,ctimes,ctimesdate,junkx,junky=[],[],[],[],[],[],[],[],[]
			for each in range(len(rayz[0,:])):
	#                refiner = tri.UniformTriRefiner(triang)
	#                tri_refi, z_refi = refiner.refine_field(rayz[0:,each], subdiv=0)
				cs=plt.tricontour(triang,rayz[0:,each],mask=mask, levels=levels,colors='c', linewidths=[1.5])
				contour = cs.collections[0].get_paths()
				for alls in range(len(contour)):
					vs=contour[alls].vertices
					a = Ut.PolygonArea(vs)
					aa.append(a)
					x = vs[:,0]
					y = vs[:,1]
					points = np.array([x,y])
					points = points.transpose()
					sx = sy = sL = 0
					for i in range(len(points)):   # counts from 0 to len(points)-1
						x0, y0 = points[i - 1]     # in Python points[-1] is last element of points
						x1, y1 = points[i]
						L = ((x1 - x0)**2 + (y1 - y0)**2) ** 0.5
						sx += (x0 + x1)/2 * L
						sy += (y0 + y1)/2 * L
						sL += L
					xc.append(sx/sL)
					yc.append(sy/sL)
				if aa != []:
					idi = np.where(np.array(aa) > kval)
					filler = np.where(np.array(aa) <= kval)
					chained = itertools.chain.from_iterable(filler)
					chain = itertools.chain.from_iterable(idi)
					idi = list(chain)
					filler = list(chained)
					for alls in range(len(aa)):
						if aa[alls] > kval:
							centroids.append([xc[idi[0]],yc[idi[0]]])
							ctimes.append(timevector[each])
							ctimesdate.append(timevector[each])
							av.append(aa[idi[0]])
						else:
							centroids.append([0,0])
							ctimes.append(timevector[each])
							ctimesdate.append(timevector[each])
							av.append(0)
		
				aa,yc,xc=[],[],[]
	#%%     Filter peaks in av above threshold by time and distance to remove redundant.
			idxx,idx,regionals,localev=[],[],[],[]
			coordinatesz = np.transpose(centroids)
			avz=av
			abovek=np.where(np.array(avz)>0)
			idxx=abovek[0]
			iii=[]
			for i in range(len(abovek[0])-1):
				junk = ctimes[idxx[i+1]]-ctimes[idxx[i]]
				junk1 = centroids[idxx[i]]
				junk2 = centroids[idxx[i+1]]
				if junk.seconds < time_thres and pydist.vincenty(junk2,junk1).meters < distance_thres:
					iii.append(idxx[i+1])
		
			idxx=set(idxx)-set(iii)
			idxx=list(idxx)
			idxx.sort()
			idx=idxx
			ltxlocal,ltxlocalexist=[],[]
			ltxglobal=[]
			ltxglobalexist=[]
			doubles,localev = [],[]
			dit2=[]
	#%%           
	#if there are no picks but cataloged events exist, make null arrays           
			if len(idx) == 0 and len(globalE) >0:
				ltxglobalexist = np.ones(len(globalE))
			if len(idx) == 0 and len(localE) >0:
				ltxlocalexist = np.ones(len(localE))
	#try to match detections with known catalog events based on time and location
			if len(idx) > 0:
				distarray=[]
				dmin= np.zeros([5])
				dval= np.zeros([5])
				closestl = np.empty([len(idx),5])
				dvals = np.empty([len(idx),5])
				closestl=closestl.astype(np.int64)
				for i in range(len(idx)):
					#find distance to the 5 nearest stations and save them for plotting templates
					for each in range(len(ll)):
						distarray.append(pydist.vincenty([coordinatesz[1][idx[i]],coordinatesz[0][idx[i]]],[ll[each],lo[each]]).meters)
					for all5 in range(5):
						dmin[all5] =np.argmin(distarray)
						dmin=dmin.astype(np.int64)
						dval[all5]=distarray[dmin[all5]]
						distarray[dmin[all5]]= 9e10
					closestl[i][:]=dmin
					dvals[i][:]=dval
					dmin=np.zeros_like(dmin)                
					distarray=[]
					#get timeseries for this pick
					stg=slist[closestl[i][0]]
					timeindex=bisect.bisect_left(alltimes, ctimes[idx[i]])
					sss=sz.select(station=stg)
					av = sss[0].data[timeindex-tlength:timeindex+tlength]
					cf=recSTALTA(av, int(40), int(1200))
					peaks = triggerOnset(cf, 3, .2)
					#get rid of peaks that are way off LTX times
					peaksi=[]   
					for peak in peaks:
						peak=peak[0]
						junk=alltimes[timeindex]-alltimes[timeindex-tlength+peak]
						if abs(junk.seconds) >45:
							peaksi.append(i) 
						
					peaks= np.delete(peaks,peaksi,axis=0)
					#look for overlap with ANF global

					for j in range(len(globalE)):
						#get distance between stations and depth for theoretical ttime calc
						# the time buffers are somewhat arbitrary
						dep = globalE.depth[j]
						dit = loc2d(centroids[idx[i]][1],centroids[idx[i]][0], globalE.Lat[j],globalE.Lon[j])
						arrivals = model.get_travel_times(dep,dit,phase_list=['P'])
						#if no calculated tt but sta/lta peak
						if len(arrivals) == 0 and len(peaks)!=0:
							junk= UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]) - UTCDateTime(globalE.DateString[j])
							if junk > -40 and junk < 40:
								doubles.append(idx[i])
								ltxglobal.append(UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]))
								ltxglobalexist.append(0)
							else:
								ltxglobalexist.append(1)
						#if no calculated tt and no sta/lta peak use ltx time
						elif len(arrivals) == 0 and len(peaks)==0:
							junk= UTCDateTime(alltimes[timeindex]) - UTCDateTime(globalE.DateString[j])
							if junk > -40 and junk < 40:
								doubles.append(idx[i])
								ltxglobal.append(UTCDateTime(alltimes[timeindex]))
								ltxglobalexist.append(0)
							else:
								ltxglobalexist.append(1)
						#if there are calculated arrivals and sta/lta peak
						elif len(peaks) != 0:
							junk= UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]) - (UTCDateTime(globalE.DateString[j]) + datetime.timedelta(seconds = arrivals[0].time))
							if junk > -30 and junk < 30:
								doubles.append(idx[i])
								ltxglobal.append(UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]))
								ltxglobalexist.append(0)
							else:
								ltxglobalexist.append(1)
						#if there are calculated arrivals and no sta/lta peaks
						else:
						
							junk= UTCDateTime(alltimes[timeindex]) - (UTCDateTime(globalE.DateString[j]) + datetime.timedelta(seconds = arrivals[0].time))
							if junk > -60 and junk < 60:
								doubles.append(idx[i])
								ltxglobalexist.append(0)
							else:
								ltxglobalexist.append(1)
					#look for overlap with ANF local
					if len(localE) > 0 and len(peaks) != 0:
						for eachlocal in range(len(localE)):
							#junk= UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]) - UTCDateTime(localE.DateString[eachlocal])
							#took this out because faulty picks disassociated too many events
							#calculate with LTX pick time instead
							dep = localE.depth[eachlocal]
							dit = pydist.vincenty([centroids[idx[i]][1],centroids[idx[i]][0]], [localE.Lat[eachlocal],localE.Lon[eachlocal]]).meters
							junk= UTCDateTime(alltimes[timeindex]) - UTCDateTime(localE.DateString[eachlocal])
							if junk > -60 and junk < 60 and dit <2.0*edgeL:
								localev.append(idx[i])
								ltxlocal.append(UTCDateTime(alltimes[timeindex-tlength+peaks[0][0]]))
								ltxlocalexist.append(0)
							else:
								ltxlocalexist.append(1)
					if len(localE) > 0 and len(peaks) ==0:
						for eachlocal in range(len(localE)):
							dep = localE.depth[eachlocal]
							dit = pydist.vincenty([centroids[idx[i]][1],centroids[idx[i]][0]], [localE.Lat[eachlocal],localE.Lon[eachlocal]]).meters
							junk= UTCDateTime(alltimes[timeindex]) - UTCDateTime(localE.DateString[eachlocal])
							if junk > -60 and junk < 60 and dit <2.0*edgeL:
								localev.append(idx[i])
								ltxlocal.append(UTCDateTime(alltimes[timeindex]))
								ltxlocalexist.append(0)
							else:
								ltxlocalexist.append(1)
				#if it goes with a local- don't let it go with a double too
				dupe=[]
				for dl in range(len(doubles)):
					if localev.count(doubles[dl]) >0:
						dupe.append(doubles[dl])
				for repeats in range(len(dupe)):
					doubles.remove(dupe[repeats])
	#
				detections = []
				detections = set(idx)#-set(doubles)
				detections = list(detections)
				#or if there are more locals LTX detections than ANF locals, fix it
				#by assuming double pick on closest pair
				pdist=[]
				if len(localev) > len(localE):
					for i in range(len(localev)-1):
						pdist.append(localev[i+1]-localev[i])
						junk=np.where(pdist==min(pdist))
					localev.pop(junk[0][0]+1)
					#detections.remove(localev[junk[0][0]+1])
				detections.sort()
				idx = detections
				dtype,cents = Ut.markType(detections,centroids,localev,localE,ctimes,doubles)
				#get the nearest station also for cataloged events
				closestd = np.zeros([len(doubles)])
				distarray = np.zeros([len(ll)])
				for event in range(len(doubles)):
					for each in range(len(ll)):
						distarray[each]=pydist.vincenty([coordinatesz[1][doubles[event]],coordinatesz[0][doubles[event]]],[ll[each],lo[each]]).meters
				
					finder = np.argmin(distarray)
					closestd[event]=finder
					distarray[finder] = 9e10
					closestd=closestd.astype(np.int64)
				
				closestp = []                
				distarray = np.zeros([len(ll)])
				for event in range(len(localev)):
					for each in range(len(ll)):
						distarray[each]=pydist.vincenty([coordinatesz[1][localev[event]],coordinatesz[0][localev[event]]],[ll[each],lo[each]]).meters
				
					finder = np.argmin(distarray)
					closestp.append(finder)
					distarray[finder] = 9e10
				

	#%%#save templates from this round of picks to verify on closest station
			ss = str(tt)
			ss = ss[0:13] 
			if 'detections' in locals():
				index = range(len(detections))
			else: 
				index=[0]
				detections = []
			df = pd.DataFrame(data=d, index=index)
			if maketemplates == 1 and len(detections) > 0:
				ptimes,confidence = [],[]
				magi = np.zeros_like(dvals)
				dum=0
				for fi in range(len(detections)):
					if localev.count(detections[fi]) == 0:
						df.Contributor[fi]='LTX'
					else:
						df.Contributor[fi]='ANF,LTX'
						allmags = [localE.ms[dum],localE.mb[dum],localE.ml[dum]]
						df.Magnitude[fi]=np.max(allmags)
						dum = dum+1
					#df.Latitude[fi] = coordinatesz[1][detections[fi]]
					#df.Longitude[fi]=coordinatesz[0][detections[fi]]
					df.Latitude[fi] = cents[fi][0]
					df.Longitude[fi]=cents[fi][1]
					df.Type[fi] = dtype[fi]
					plt.cla()
					ax = plt.gca()
					timeindex=bisect.bisect_left(alltimes, (ctimes[detections[fi]]))
					sss = np.zeros([5,tlength*2])
					for stas in range(5):
						stg = slist[closestl[fi][stas]]
						tr = sz.select(station=stg)
						if ctimes[detections[fi]]-datetime.timedelta(seconds=80) < tt.datetime:
							sss[stas][tlength:]=tr[0].data[timeindex:timeindex+tlength] 
						elif ctimes[detections[fi]]+datetime.timedelta(seconds=80) > tt.datetime + datetime.timedelta(seconds=nseconds+edgebuffer):
							sss[stas][0:tlength]=tr[0].data[timeindex-tlength:timeindex]
						else:
							sss[stas][:]=tr[0].data[timeindex-tlength:timeindex+tlength]
					sss=np.nan_to_num(sss)
					stg=slist[closestl[0][0]]    
					#plt.figure(fi)
					peak=None
					plt.suptitle('nearest station:'+stg+' '+str(ctimes[detections[fi]])+'TYPE = '+dtype[fi])
					for plots in range(5):
						plt.subplot(5,1,plots+1)
						cf=recSTALTA(sss[plots][:], int(80), int(500))
						peaks = triggerOnset(cf, 3, .1)
						peaksi=[]
						dummy=0 
						for pk in peaks:
							endi=pk[1]
							peak=pk[0]
							mcdur=alltimes[timeindex-tlength+endi]-alltimes[timeindex-tlength+peak]
							mdur=mcdur.total_seconds()
							if alltimes[timeindex]>alltimes[timeindex-tlength+peak]:
								junk=alltimes[timeindex]-alltimes[timeindex-tlength+peak]
							else:
								junk=alltimes[timeindex-tlength+peak]-alltimes[timeindex]
							if (junk.seconds) >40:
								peaksi.append(dummy)
							dummy=dummy+1
						peaks= np.delete(peaks,peaksi,axis=0)                            
					
						sss[plots]=np.nan_to_num(sss[plots])
						#if the channel is blank underflow problems occur plotting station name
						sss = np.round(sss,decimals=10)
						plt.plot(Ut.templatetimes(alltimes[timeindex],tlength,deltaf),sss[plots][:],'black')
						plt.axvline(x=alltimes[timeindex])
						plt.text(alltimes[timeindex],0,slist[closestl[fi][plots]],color='red',fontsize=20)
						plt.axis('tight')
						for arc in range(len(peaks)):
							plt.axvline(x=alltimes[timeindex-tlength-10+peaks[arc][0]],color='orange')
							plt.axvline(x=alltimes[timeindex-tlength-10+peaks[arc][1]],color='purple')
					
						if len(peaks)>0:
							ptimes.append(UTCDateTime(alltimes[timeindex-tlength-10+peaks[0][0]]))
							confidence.append(len(peaks))
							magi[fi][plots]=(-2.25+2.32*np.log10(mdur)+0.0023*dvals[fi][plots]/1000)
							#magi[fi][plots]=(1.86*np.log10(mdur)-0.85)
						else:
							ptimes.append(UTCDateTime(alltimes[timeindex]))
							confidence.append(2)
				

					magi= np.round(magi,decimals=2)
					magii = pd.DataFrame(magi)
					magu= magii[magii != 0]
					if df.Contributor[fi]=='ANF,LTX':
						df.mag_error[fi]= np.round(np.max(allmags)-np.mean(magu,axis=1)[fi],decimals=2)                    
						df.Magnitude[fi]=str(str(df.Magnitude[fi])+','+str(np.round(np.mean(magu,axis=1)[fi],decimals=2)))
						df.cent_er[fi] = np.round(pydist.vincenty([coordinatesz[1][detections[fi]],
							coordinatesz[0][detections[fi]]],[cents[fi][0],cents[fi][1]]).meters/1000.00,decimals=2)
					else:
						df.Magnitude[fi]=np.round(np.mean(magu,axis=1)[fi],decimals=2)
					#ptimes = np.reshape(ptimes,[len(ptimes)/5,5])       
					df.S1[fi]= slist[closestl[fi][0]]
					df.S1time[fi] = ptimes[0]
					df.S2[fi]= slist[closestl[fi][1]]
					df.S2time[fi] = (ptimes[1])
					df.S3[fi]= slist[closestl[fi][2]]
					df.S3time[fi] = (ptimes[2])
					df.S4[fi]= slist[closestl[fi][3]]
					df.S4time[fi] = (ptimes[3])
					df.S5[fi]= slist[closestl[fi][4]]
					df.S5time[fi] = (ptimes[4])
					#df.Confidence[fi]= confidence[0]
					ptimes = []
					if dtype[fi]=='earthquake':
						svname=homedir+str(s)+"/image"+ss[11:13]+"_pick_"+str(fi+1)+".png"
						plt.savefig(svname,format='png')
					plt.clf()

	#%%
			df1 = [df1,df]
			df1= pd.concat(df1)

		################################################
	#%%     
		
			fig = plt.figure()
			plt.cla()
			ax = plt.gca()
			#plot it all
			for i in range(len(detections)):
			
				if localev.count(detections[i]) ==1:
					color='c'
				elif doubles.count(detections[i])==1:
					color='blue'
				else:
					color='white'
				if dtype[i] =='blast':
					facecolor='none'
				else: 
					facecolor = color
				plt.scatter(mdates.date2num(ctimes[detections[i]]),closestl[i][0],s=200,color=color,facecolor=facecolor)
		   
	#         
			for i in range(len(globalE)):
				plt.scatter(mdates.date2num(UTCDateTime(globalE.time[i])),1,s=100, color='b', alpha=.8)
			for i in range(len(localE)):
				plt.scatter(mdates.date2num(UTCDateTime(localE.time[i])),closesti[i],s=100,facecolor='c',edgecolor='grey')
			plt.imshow(np.flipud(rayz),extent = [mdates.date2num(tt.datetime), mdates.date2num((tt+nseconds+edgebuffer).datetime),  0, len(slist)],
						 aspect='auto',interpolation='nearest',cmap='bone',vmin=np.min(rayz)/2,vmax=np.max(rayz)*2)

			ax.set_adjustable('box-forced')
			ax.xaxis_date() 
			plt.yticks(np.arange(len(ll)))
			ax.set_yticklabels(slist)
			tdate = yr+'-'+mo+'-'+str(dayat).zfill(2)
			plt.title(tdate)
			ax.grid(color='black')
			ss = str(tt)
			ss = ss[0:13]
			kurs = "%s/"%s +"%s.png"%ss
			svpath=homedir+kurs
			plt.savefig(svpath, format='png')
			plt.close()
		
			#%%
			blockette = blockette+(npts-nptsf)
			tt = tt+nseconds
			detections=[]
			localev=[]
			doubles=[] 
		#############################
	
		svpath = homedir+'%s'%s+"/picktable.html"  
		df1.to_html(open(svpath, 'w'),index=False)
		svpath = homedir+'%s'%s+"/picktable.pkl"  
		df1.to_pickle(svpath)    
		dayat = dayat+1
		counter=counter+1
		counter_3char = str(counter).zfill(3)

		#############################
	
	if __name__ == '__main__':
		detection_function()
