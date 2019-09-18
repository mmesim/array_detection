#!/bin/bash
#Pseudo parallel run for frequency detector                                         #  
#This scripts generates N, where N is the number of days defined by the user,       #
#scripts and runs them in background                                                #
#Preferably run in chpc            -- see manual for instructions                   #
#################################################### MM 04/2019 #####################
date
# Decide how you want to run the detector (uncomment the appropriate line)
# answer=1 : run locally, use one core, and run all days
#answer=1
# answer=2 : run in chpc or locally using multiple cores, assign each day to a core
answer=2
#------------------------------------------------------------------------------------
#Define parameters 
#------------------------------------------------------------------------------------
#Data directory
datdir='/home/mesimeri/Documents/WORK/FREQUENCY_DETECTOR/SOURCE/'
#------------------------------------------------------------------------------------
#Starting time
yr=2016   # year 
mo=12     #month
dy=18     #day
hr=00     #hour
mn=00     #minute
sc=00    #second
#------------------------------------------------------------------------------------
#Duration and sampling rate
days=1          #total number of days
duration=86400  #one day in seconds
deltaf=250     #delta= 1/sampling rate
#------------------------------------------------------------------------------------
#Detection parameters
freq1=15                # minimum frequency to sum amplitudes
freq2=25                # maximum frequency to sum amplitudes
thresholdv=2.0           # area threshold
masktimes=10.0          # N times average interstation distance
madtimes=4.0            # N times the MAD
time_thres=20.0          #time detection threshold in sec
distance_thres=250000   # distance detection threshold in meters
#-------------------------------------------------------------------------------------
# Do not edit below this line
#-------------------------------------------------------------------------------------
case $answer in
 1)
 echo 'Work with one core'
 cat > run_freqtor.py <<EOF 
from detection_script import freqtor

freqtor(yr ='$yr', mo='$mo', dy='$dy', hr='$hr', mn='$mn', sc='$sc', duration=$duration, ndays=$days, datdir='$datdir', freq1=$freq1, freq2=$freq2, thresholdv=$thresholdv, deltaf=$deltaf, masktimes=$masktimes, madtimes=$madtimes, time_thres=$time_thres, distance_thres=$distance_thres)
EOF
chmod +x run_freqtor.py
python run_freqtor.py
 ;;
 2)
 echo 'Work with multiple cores'
 stop=`expr $dy + $days`
for day in $(seq $dy $stop);do
cat > run_freqtor_$day.py <<EOF
from detection_script import freqtor

freqtor(yr ='$yr', mo='$mo', dy='$day', hr='$hr', mn='$mn', sc='$sc', duration=$duration, ndays=1, datdir='$datdir', freq1=$freq1, freq2=$freq2, thresholdv=$thresholdv, deltaf=$deltaf, masktimes=$masktimes, madtimes=$madtimes, time_thres=$time_thres, distance_thres=$distance_thres)
EOF
chmod +x run_freqtor_$day.py
python run_freqtor_$day.py &
done
 ;;
 *)
 echo 'Please uncomment line 10 or 12 and try again!'
 exit
 ;;
esac

date
exit
