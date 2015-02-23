#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Usage: ./go_plots.sh kSR2 (kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0)"
	echo
	exit 0
fi

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

PHOTON_REGION=$1

rm limitInputs_${PHOTON_REGION#?}.root
rm met_differences_${PHOTON_REGION#?}.root

cat makePlots_template.C | sed s:PHOTON_REGION:$PHOTON_REGION: > makePlots.C
root -b -q -l makePlots.C
rm makePlots.C

file=`date +"%b%d"`

mkdir $file
mv *.pdf $file
cd $file

tar -czf $file.tgz *
scp $file.tgz $HEP

mv $file.tgz ..
cd ..
rm -r $file