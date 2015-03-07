#!/bin/bash

if [ $# -ne 2 ]; then
	echo
	echo "Usage: ./go_plots.sh kSR2 false (kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot)"
	echo "                     (kSR1, kSR2, kCR1, kCR2, kCR2a, kCR0, kSigmaPlot)"
	echo "                     useSuperFakes"
	echo
	exit 0
fi

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

PHOTON_REGION=$1
USE_SUPER_FAKES=$2

cat makeHistograms_template.C | sed s:PHOTON_REGION:$PHOTON_REGION: | sed s:USE_SUPER_FAKES:$USE_SUPER_FAKES: > makeHistograms.C
root -b -q -l makeHistograms.C 2>&1 | sed '/does not exist/d'
rm makeHistograms.C
