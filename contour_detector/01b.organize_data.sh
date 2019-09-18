#/bin/bash

#Path to SAC files - organized by station name
workdir='/uufs/chpc.utah.edu/common/home/koper-group1/mesimeri/FORGE_2019/SAC'
#Directory to save the data organized by date
mydata='/uufs/chpc.utah.edu/common/home/koper-group1/mesimeri/FORGE_2019/DATA'
###############################################################################

mkdir -p $mydata

for station in $workdir/* ; do
sta=`basename $station`
  for file in $workdir/$sta/* ; do
  dirname=`echo $file |basename $file|  cut  -c1,2,3,4,5,6,7,8`
  mkdir -p $mydata/$dirname
  cp -v  $file $mydata/$dirname
  done
done
exit
