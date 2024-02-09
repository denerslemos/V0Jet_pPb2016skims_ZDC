#!/bin/bash

echo "Setup CMSSW (ROOT version)"
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/d/ddesouza/UIC/pPbMultAna/CMSSW_13_0_5/src/V0inUPC/V0Jet_pPb2016skims_ZDC
echo "Submit skim jobs at "
echo PWD: $PWD

./V0Jet_pPb2016skims_ZDC $1 $2 $3
