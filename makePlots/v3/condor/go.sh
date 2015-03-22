#!/bin/bash

WORK_DIR=`pwd`

PHOTON_REGION=$1
PHOTON_MODE=$2

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_8_patch3
cd CMSSW_5_3_8_patch3/src/
eval `scramv1 runtime -sh`

mkdir acc/
mv $WORK_DIR/makeHistograms.C .
mv $WORK_DIR/CreateHistograms_C.* .

sed -i "s/PHOTON_REGION/$PHOTON_REGION/g" makeHistograms.C
sed -i "s/PHOTON_MODE/$PHOTON_MODE/g" makeHistograms.C

root -b -q -l makeHistograms.C

mv limitInputs_bjj.root limitInputs_bjj_$PHOTON_REGION.root
mv limitInputs_jjj.root limitInputs_jjj_$PHOTON_REGION.root

mv *.root $WORK_DIR
mv acc/*.pdf $WORK_DIR

cd $WORK_DIR
rm -rf CMSSW_5_3_8_patch3/
rm makeHistograms_template.C
