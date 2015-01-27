#!/bin/bash

if [ $# -ne 1 ]; then
	echo
	echo "Usage: ./go_plots.sh kSR2 (kSR1, kSR2, kCR1, kCR2)"
	echo
	exit 0
fi

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

PHOTON_REGION=$1
ELE_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs_v2/SingleElectron.root
MUON_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs_v2/SingleMu.root

rm *.root

cat makePlots_template.C | sed s:ELE_FILE_TO_RUN:$ELE_FILE_TO_RUN: | sed s:MUON_FILE_TO_RUN:$MUON_FILE_TO_RUN: | sed s:PHOTON_REGION:$PHOTON_REGION: > makePlots.C
root -b -q -l makePlots.C 2>&1 | sed '/does not exist/d' | sed '/has been created/d'
rm makePlots.C