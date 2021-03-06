#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp7.C+\(\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_first_jet_mass_b\",1,$i\)
 
  root -l -q -b DataMCComp7.C+\(\"w_second_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_second_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_second_jet_mass_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_dijet_pt_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_dijet_eta_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_dijet_mass_b\",1,$i\)
 
  root -l -q -b DataMCComp7.C+\(\"w_Ht_b\",1,$i\)

  root -l -q -b DataMCComp7.C+\(\"w_pt_W_wenu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_pt_W_wmnu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_eta_W_wenu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_eta_W_wmnu_b\",1,$i\)

  root -l -q -b DataMCComp7.C+\(\"w_delta_wenu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_delta_wmnu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_deltaR_wenu_b\",1,$i\)
  root -l -q -b DataMCComp7.C+\(\"w_deltaR_wmnu_b\",1,$i\)

  i=$((i+1))
done

exit
