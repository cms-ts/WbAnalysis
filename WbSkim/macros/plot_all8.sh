#!/bin/sh

export ROOT_HIST=0

unset PYTHIA8175DATA
unset G4NEUTRONXS

cd $CMS_PATH/slc6_amd64_gcc472/cms/cmssw/CMSSW_6_2_7
eval `scramv1 runtime -sh`
cd -

root -l -q -b DataMCComp8.C+\(\"w_first_jet_pt_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_first_jet_eta_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_first_jet_mass_b\",1\)

root -l -q -b DataMCComp8.C+\(\"w_second_jet_pt_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_second_jet_eta_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_second_jet_mass_b\",1\)

root -l -q -b DataMCComp8.C+\(\"w_dijet_pt_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_dijet_eta_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_dijet_mass_b\",1\)

root -l -q -b DataMCComp8.C+\(\"w_Ht_b\",1\)

root -l -q -b DataMCComp8.C+\(\"w_pt_W_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_eta_W_b\",1\)

root -l -q -b DataMCComp8.C+\(\"w_delta_b\",1\)
root -l -q -b DataMCComp8.C+\(\"w_deltaR_b\",1\)

exit
