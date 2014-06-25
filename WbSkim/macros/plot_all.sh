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

i=7
while [ $i -le 8 ]; do

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_wide\",1,$i,0,1\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b_wide\",1,$i,0,1\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb_wide\",1,$i,0,1\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_wide\",1,$i,0,1\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b_wide\",1,$i,0,1\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb_wide\",1,$i,0,1\)

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_b\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_bb\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_b\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_bb\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i,1,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_b\",1,$i,1,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_bb\",1,$i,1,3\)

  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_b\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_bb\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_b\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_bb\",1,$i,0,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i,1,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_b\",1,$i,1,3\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_bb\",1,$i,1,3\)

  i=$((i+1))
done

i=1
while [ $i -le 8 ]; do

  root -l -q -b DataMCComp.C+\($d,\"h_pu_weights\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"h_recoVTX\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_recoVTX\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"h_tracks\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_tracks\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_ele_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_ele_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_ele_iso\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_muon_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_muon_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_muon_iso\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_wide\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_wide\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_b_wide\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_b_wide\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_wide\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_wide\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b_wide\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b_wide\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb_wide\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb_wide\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wenu_blepton\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wmnu_blepton\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_ee\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_mm\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_ee\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_mm\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wenu\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wmnu\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wenu\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wmnu\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wenu_blepton_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wmnu_blepton_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_pt_Z_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_ee_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_mm_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wenu_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wmnu_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wenu_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wmnu_b\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_mass_ee_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_mm_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wenu_blepton_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mass_wmnu_blepton_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wenu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wmnu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wenu_2b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_delta_wmnu_2b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wenu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wmnu_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wenu_2b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_deltaR_wmnu_2b\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_jetmultiplicity\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_mass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_mass\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_mass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_mass_b\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_pt_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_eta_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_jet_mass_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_pt_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_eta_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_jet_mass_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_bjetmultiplicity\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_first_bjet_mass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_pt\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_eta\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_second_bjet_mass\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_MET\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_MET_sign_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_Ht\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_Ht_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_Ht_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_first_ele\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"h_scaleFactor_first_muon\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_jet\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_jet_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_jet_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_trk\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_trk_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_trk_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_SVTX_mass_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_zoom\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_mass\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_secondvtx_N_nomass\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_mass_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_mass_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP_nomass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_nomass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP_nomass_bb\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_nomass\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_nomass_b\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_JBP_nomass_bb\",1,$i\)

  root -l -q -b DataMCComp.C+\($d,\"w_BJP0\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP1\",1,$i\)
  root -l -q -b DataMCComp.C+\($d,\"w_BJP2\",1,$i\)

  i=$((i+1))
done

i=1
while [ $i -le 2 ]; do

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_wide\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b_wide\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb_wide\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_wide\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b_wide\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb_wide\",1,$i,1,0\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b_wide\",1,$i,0,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb_wide\",1,$i,0,2\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_wide\",1,$i,1,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b_wide\",1,$i,1,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb_wide\",1,$i,1,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_wide\",1,$i,1,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b_wide\",1,$i,1,2\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb_wide\",1,$i,1,2\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b\",1,$i,1,0\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb\",1,$i,1,0\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu\",1,$i,0,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b\",1,$i,0,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb\",1,$i,0,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu\",1,$i,0,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b\",1,$i,0,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb\",1,$i,0,4\)

  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu\",1,$i,1,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_b\",1,$i,1,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wenu_bb\",1,$i,1,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu\",1,$i,1,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_b\",1,$i,1,4\)
  root -l -q -b DataMCComp.C+\($d,\"w_mt_wmnu_bb\",1,$i,1,4\)

  i=$((i+1))
done

exit
