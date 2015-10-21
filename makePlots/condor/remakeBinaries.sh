#!/bin/bash

WORK_DIR=`pwd`

cd ../
source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`

root -b -q -l compile.C

cp CreateHistograms_C.so $WORK_DIR
cp CreateHistograms_C.d $WORK_DIR

cd $WORK_DIR
