#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C+\(\"w_mt_wenu\",1,$i,0\)
  root -l -q -b DataMCComp7.C+\(\"w_mt_wmnu\",1,$i,0\)

  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_pt\",1,$i,0\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta\",1,$i,0\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_mass\",1,$i,0\)
 
  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_pt\",1,$i,0\)
  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_eta\",1,$i,0\)
  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_mass\",1,$i,0\)
 
  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C+\(\"w_mt_wenu\",1,$i,1\)
  root -l -q -b DataMCComp7.C+\(\"w_mt_wmnu\",1,$i,1\)

  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_pt\",1,$i,1\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_eta\",1,$i,1\)
  root -l -q -b DataMCComp7.C+\(\"w_first_bjet_mass\",1,$i,1\)

  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_pt\",1,$i,1\)
  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_eta\",1,$i,1\)
  root -l -q -b DataMCComp7.C+\(\"w_second_bjet_mass\",1,$i,1\)
 
  i=$((i+1))
done

exit
