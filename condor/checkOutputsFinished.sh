#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Usage: ./go_plots.sh Jan15(the date)"
	echo
	exit 0
fi

folder=$1
storage="/eos/uscms/store/user/bfrancis/inputs_v7/signal_contamination_"

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

for x in `ls -d *_$folder | grep -v "SingleElectron" | grep -v "SingleMu"`
do

	[ "$(ls $storage*.root)" ] || continue

	[ -f $storage${x%_$folder}.root ] && continue

	if [ `ls $storage${x%_$folder}_*.root | wc -l` -eq `ls $x/filelist_* | wc -l` ]
	then
	    echo
	    echo hadd $storage${x%_$folder}.root $storage${x%_$folder}_\*.root
	    echo
	    hadd $storage${x%_$folder}.root $storage${x%_$folder}_*.root
	    if [ ! -f $storage${x%_$folder}.root ]
	    then
		echo Failed!
		exit 0
	    fi
	    rm $storage${x%_$folder}_*.root
	    rm -r $x
	fi

done
