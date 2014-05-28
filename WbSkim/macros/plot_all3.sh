#!/bin/sh

d=0
if [ ! -z "$1" ]; then
  d=$1
fi

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp3.C+\($d,\"w_mt_wenu_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mt_wenu_bb\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mt_wmnu_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_mt_wmnu_bb\",1,$i\)

  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_pt_bb\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_eta_bb\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_mass_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_first_jet_mass_bb\",1,$i\)
  
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_pt_bb\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_eta_bb\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_mass_b\",1,$i\)
  root -l -q -b DataMCComp3.C+\($d,\"w_second_jet_mass_bb\",1,$i\)

  i=$((i+1))
done

exit
