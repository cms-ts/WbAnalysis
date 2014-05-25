#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

root -l -q -b DataMCComp8.C+\(\"w_mt\",1,0\)

root -l -q -b DataMCComp8.C+\(\"w_first_bjet_pt\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_eta\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_mass\",1,0\)

root -l -q -b DataMCComp8.C+\(\"w_second_bjet_pt\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_second_bjet_eta\",1,0\)
root -l -q -b DataMCComp8.C+\(\"w_second_bjet_mass\",1,0\)

root -l -q -b DataMCComp8.C+\(\"w_mt\",1,1\)

root -l -q -b DataMCComp8.C+\(\"w_first_bjet_pt\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_eta\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_first_bjet_mass\",1,1\)

root -l -q -b DataMCComp8.C+\(\"w_second_bjet_pt\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_second_bjet_eta\",1,1\)
root -l -q -b DataMCComp8.C+\(\"w_second_bjet_mass\",1,1\)

exit
