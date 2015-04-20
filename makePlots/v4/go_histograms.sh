#!/bin/bash

if [ $# -ne 2 ]; then
	echo
	echo "Usage: ./go_plots.sh kSR2 kSignal"
	echo "                     (kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot, kAny)"
	echo "                     (kSignal, kFake)"
	echo
	exit 0
fi

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

PHOTON_REGION=$1
PHOTON_MODE=$2

cat makeHistograms_template.C | sed s:PHOTON_REGION:$PHOTON_REGION: | sed s:PHOTON_MODE:$PHOTON_MODE: > makeHistograms.C
root -b -q -l makeHistograms.C 2>&1 | sed '/does not exist/d'
rm makeHistograms.C
