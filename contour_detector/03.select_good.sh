#!/bin/bash

workdir=$(pwd)

datadir=$workdir/basin5
tempdir=$workdir/temp
evtlst=$tempdir/evt.lst

## functions
function select_good () {
cd $datadir
ls > $evtlst

for tstamp in $(cat $evtlst); do
	cd $datadir/$tstamp
	echo $tstamp

	if [[ ! -d BAD ]]; then
		mkdir BAD
	fi

        nimg=$(ls image* | wc -l) 
        if [[ $nimg -lt 1 ]]; then
                continue
        fi  

	for tseg in $(ls 201*); do
		time=$(echo $tseg | awk -F"T" '{print $2}' | awk -F"." '{print $1}')
		echo "Time is $time!"
		open $tseg
		pid1=$(ps aux | grep -i 'preview' | awk '$11 != "grep"{print $2}')
		echo "Please press [ENTER]: $pid1"
		read pause
		nline=$(ls image$time* | wc -l)
		if [[ $nline -lt 1 ]]; then
			continue
		fi
        kill -9 $pid1
        for i in $(seq 1 $nline); do
			open image$time*_$i.png
            pid2=$(ps aux | grep -i 'preview' | awk '$11 != "grep"{print $2}')
#			pid2=$(ps aux | grep -i 'preview' | awk '$11 != "grep" && $2 != pid1 {print $2}' "pid1=$pid1")
			echo "Please input: 0 bad, 1 good: $i"
			read igood # igood 1 is good, 0 is bad
			if [[ $igood -ne 1 ]]; then
				mv image$time*_$i.png BAD
			fi
			kill -9 $pid2
		done
#		kill -9 $pid1
	done
done

return 0
}

## Exective
date

select_good

date

exit 0
