
#include "LumiInfo_v12.h"

TList *FileList;
TFile *Target;

void hmerge_(TDirectory *target, TList *sourcelist, double crossArray[]);

void hmerge(string path=".", string version="v01", string title="W") {

  FileList = new TList();
  double crossSections[10];

  if (title=="W"||title=="W_patgen"||title=="W_gen") {
    if (title=="W") {
      Target = TFile::Open((path + "/" + version + "/" + "Wj_merge.root").c_str(), "RECREATE");
      FileList->Add(TFile::Open((path + "/" + version + "/" + "Wj.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W1j.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W2j.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W3j.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W4j.root").c_str()));
    }
    if (title=="W_patgen") {
      Target = TFile::Open((path + "/" + version + "/" + "Wj_patgen_merge.root").c_str(), "RECREATE");
      FileList->Add(TFile::Open((path + "/" + version + "/" + "Wj_patgen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W1j_patgen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W2j_patgen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W3j_patgen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W4j_patgen.root").c_str()));
    }
    if (title=="W_gen") {
      Target = TFile::Open((path + "/" + version + "/" + "Wj_gen_merge.root").c_str(), "RECREATE");
      FileList->Add(TFile::Open((path + "/" + version + "/" + "Wj_gen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W1j_gen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W2j_gen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W3j_gen.root").c_str()));
      FileList->Add(TFile::Open((path + "/" + version + "/" + "W4j_gen.root").c_str()));
    }

    crossSections[0] = (Xsec_wj / Ngen_wj);
    crossSections[1] = (Xsec_w1j / Ngen_w1j) * (Xsec_wj / 30400.);
    crossSections[2] = (Xsec_w2j / Ngen_w2j) * (Xsec_wj / 30400.);
    crossSections[3] = (Xsec_w3j / Ngen_w3j) * (Xsec_wj / 30400.);
    crossSections[4] = (Xsec_w4j / Ngen_w4j) * (Xsec_wj / 30400.);

    crossSections[0] = crossSections[0] / (Xsec_wj / Ngen_wj);
    crossSections[1] = crossSections[1] / (Xsec_wj / Ngen_wj);
    crossSections[2] = crossSections[2] / (Xsec_wj / Ngen_wj);
    crossSections[3] = crossSections[3] / (Xsec_wj / Ngen_wj);
    crossSections[4] = crossSections[4] / (Xsec_wj / Ngen_wj);
  }

  if (title=="T") {
    Target = TFile::Open((path + "/" + version + "/" + "T_merge.root").c_str(), "RECREATE");
    FileList->Add(TFile::Open((path + "/" + version + "/" + "TBar_s.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "TBar_t.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "TBar_tW.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "T_s.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "T_t.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "T_tW.root").c_str()));

    crossSections[0] = (Xsec_tbar_s / Ngen_tbar_s);
    crossSections[1] = (Xsec_tbar_t / Ngen_tbar_t);
    crossSections[2] = (Xsec_tbar_tw / Ngen_tbar_tw);
    crossSections[3] = (Xsec_t_s / Ngen_t_s);
    crossSections[4] = (Xsec_t_t / Ngen_t_t);
    crossSections[5] = (Xsec_t_tw / Ngen_t_tw);

    crossSections[0] = crossSections[0] / (Xsec_t / Ngen_t);
    crossSections[1] = crossSections[1] / (Xsec_t / Ngen_t);
    crossSections[2] = crossSections[2] / (Xsec_t / Ngen_t);
    crossSections[3] = crossSections[3] / (Xsec_t / Ngen_t);
    crossSections[4] = crossSections[4] / (Xsec_t / Ngen_t);
    crossSections[5] = crossSections[5] / (Xsec_t / Ngen_t);
  }

  if (title=="TTbar") {
    Target = TFile::Open((path + "/" + version + "/" + "TTbar_merge.root").c_str(), "RECREATE");
    FileList->Add(TFile::Open((path + "/" + version + "/" + "TTbar_FullLept.root").c_str()));
    FileList->Add(TFile::Open((path + "/" + version + "/" + "TTbar_SemiLept.root").c_str()));

    crossSections[0] = (Xsec_tt_fl / Ngen_tt_fl);
    crossSections[1] = (Xsec_tt_sl / Ngen_tt_sl);

    crossSections[0] = crossSections[0] / (Xsec_tt / Ngen_tt);
    crossSections[1] = crossSections[1] / (Xsec_tt / Ngen_tt);
  }

  //cout << "going to call merging routine..." << endl;
  hmerge_( Target, FileList, crossSections );
  //cout << "done." << endl;

}

void hmerge_( TDirectory *target, TList *sourcelist, double crossArray[] ) {

  //cout << "Target path: " << target->GetPath() << endl;
  TString rpath( (char*)strstr( target->GetPath(), ":" ) );
  rpath.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();

  first_source->cd( rpath );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // loop over all keys in this directory
  TChain *globChain = 0;
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( rpath );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> merge it

      //cout << "Merging histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
      h1->Sumw2();

      // Scale by the cross-section factor
      h1->Scale(crossArray[0]);

      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );

      int q = 1; // This keeps track of which
                 // cross section factor to use
      while ( nextsource ) {
        // make sure we are at the correct directory level by cd'ing to rpath
        nextsource->cd( rpath );
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
           TH1 *h2 = (TH1*)key2->ReadObj();
           h2->Sumw2();

           // Scale by the cross-section factor before adding.
           h2->Scale(crossArray[q]);
           h1->Add( h2 );
           q++;
           delete h2;
        }

        nextsource = (TFile*)sourcelist->After( nextsource );
      }

    } else if ( obj->IsA()->InheritsFrom( "TTree" ) ) {

      // loop over all source files create a chain of Trees "globChain"
      const char* obj_name= obj->GetName();

      globChain = new TChain(obj_name);
      globChain->Add(first_source->GetName());
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {
        globChain->Add(nextsource->GetName());
        nextsource = (TFile*)sourcelist->After( nextsource );
      }

    } else if ( obj->IsA()->InheritsFrom( "TDirectory" ) ) {

      cout << "Found subdirectory " << obj->GetName() << endl;

      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      hmerge_( newdir, sourcelist,crossArray );

    } else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
           << obj->GetName() << " title: " << obj->GetTitle() << endl;

    }

    // now write the merged histogram (which is "in" obj) to the target file
    // note that this will just store obj in the current directory level,
    // which is not persistent until the complete directory itself is stored
    // by "target->Write()" below
    if ( obj ) {
      target->cd();

      //!!if the object is a tree, it is stored in globChain...
        if(obj->IsA()->InheritsFrom( "TTree" ))
          globChain->Merge(target->GetFile(),0,"keep");
        else
        obj->Write( key->GetName() );
    }

  } // while ( ( TKey *key = (TKey*)nextkey() ) )

  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);

}

