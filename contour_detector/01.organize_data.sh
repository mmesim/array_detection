#!/bin/bash
## remove the two horizontal components from the raw data

workdir=$(cd ..;pwd)


## functions 

function reorg_dir() {
for stnm in $(seq 1 96); do
cd $workdir/$stnm
echo "$workdir/$stnm"
for file in $(ls *.sac); do
dir=${file:0:8}
if [[ ! -d $dir ]]; then
mkdir -p $workdir/$dir
fi
cp $file $workdir/$dir
done
done

return 0
}

function remove_data () {
for stnm in $(seq 1 96); do
cd $workdir/$stnm
echo "$workdir/$stnm"

rm -f *.EHE.FG.sac *.EHN.FG.sac

done

return 0
}


## exe
date

reorg_dir
remove_data

date

exit 0
