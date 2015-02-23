#!/bin/bash

WORK_DIR=`pwd`

PHOTON_REGION=$1
ELE_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs_v3/SingleElectron.root
MUON_FILE_TO_RUN=/eos/uscms/store/user/bfrancis/inputs_v3/SingleMu.root

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3

cd CMSSW_5_3_8_patch3/src/
eval `scramv1 runtime -sh`
mkdir acc/
cat $WORK_DIR/makeHistograms_template.C | sed s:ELE_FILE_TO_RUN:$ELE_FILE_TO_RUN: | sed s:MUON_FILE_TO_RUN:$MUON_FILE_TO_RUN: | sed s:PHOTON_REGION:$PHOTON_REGION: > makeHistograms.C

mv $WORK_DIR/rootRoutines.h .
mv $WORK_DIR/Binning.h .
mv $WORK_DIR/CreateHistograms.* .

root -b -q -l makeHistograms.C

mv *.root $WORK_DIR
mv acc/*.pdf $WORK_DIR

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
rm makeHistograms_template.C
