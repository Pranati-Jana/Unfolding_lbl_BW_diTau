//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  6 23:17:28 2024 by ROOT version 6.28/06
// from TTree output_tree/
// found on file: RootFile/ForDielectron/QEDSC__Costheatstar_GenReco_6thMarch.root
//////////////////////////////////////////////////////////

#ifndef ReadqedTree_h
#define ReadqedTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ReadqedTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         gen_Pt1;
   Float_t         gen_Eta1;
   Float_t         gen_Phi1;
   Float_t         gen_Pt2;
   Float_t         gen_Eta2;
   Float_t         gen_Phi2;
   Float_t         gen_diEle_M;
   Float_t         gen_diEle_Pt;
   Float_t         gen_diEle_Rapidity;
   Float_t         gen_diEle_Costhetastar;
   Int_t           run;
   Int_t           ls;
   Int_t           evtnb;
   Int_t           nEle;
   Int_t           eleCharge_1;
   Float_t         elePt_1;
   Float_t         eleEta_1;
   Float_t         elePhi_1;
   Float_t         eleSCEt_1;
   Float_t         eleSCEta_1;
   Float_t         eleSCPhi_1;
   Float_t         eleEoverP_1;
   Int_t           eleCharge_2;
   Float_t         elePt_2;
   Float_t         eleEta_2;
   Float_t         elePhi_2;
   Float_t         eleSCEt_2;
   Float_t         eleSCEta_2;
   Float_t         eleSCPhi_2;
   Float_t         eleEoverP_2;
   Float_t         vSum_M;
   Float_t         vSum_Energy;
   Float_t         vSum_Pt;
   Float_t         vSum_Eta;
   Float_t         vSum_Phi;
   Float_t         vSum_Rapidity;
   Float_t         ele_dpt;
   Float_t         ele_deta;
   Float_t         ele_dphi;
   Float_t         ele_acop;
   Int_t           ok_neuexcl;
   Int_t           ok_zdcexcl;
   Int_t           ok_chexcl;
   Int_t           ok_chexcl_extrk;
   Int_t           nExtrk;
   Int_t           ok_trigger;
   Int_t           ok_event;
   Float_t         zdc_energy_pos;
   Float_t         zdc_energy_neg;
   Float_t         costhetastar;
   Float_t         cos_photon_pair_helicity0;
   Float_t         cos_photon_pair_helicity1;

   // List of branches
   TBranch        *b_gen_Pt1;   //!
   TBranch        *b_gen_Eta1;   //!
   TBranch        *b_gen_Phi1;   //!
   TBranch        *b_gen_Pt2;   //!
   TBranch        *b_gen_Eta2;   //!
   TBranch        *b_gen_Phi2;   //!
   TBranch        *b_gen_diEle_M;   //!
   TBranch        *b_gen_diEle_Pt;   //!
   TBranch        *b_gen_diEle_Rapidity;   //!
   TBranch        *b_gen_diEle_Costhetastar;   //!
   TBranch        *b_run;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_evtnb;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge_1;   //!
   TBranch        *b_elePt_1;   //!
   TBranch        *b_eleEta_1;   //!
   TBranch        *b_elePhi_1;   //!
   TBranch        *b_eleSCEt_1;   //!
   TBranch        *b_eleSCEta_1;   //!
   TBranch        *b_eleSCPhi_1;   //!
   TBranch        *b_eleEoverP_1;   //!
   TBranch        *b_eleCharge_2;   //!
   TBranch        *b_elePt_2;   //!
   TBranch        *b_eleEta_2;   //!
   TBranch        *b_elePhi_2;   //!
   TBranch        *b_eleSCEt_2;   //!
   TBranch        *b_eleSCEta_2;   //!
   TBranch        *b_eleSCPhi_2;   //!
   TBranch        *b_eleEoverP_2;   //!
   TBranch        *b_vSum_ee_M;   //!
   TBranch        *b_vSum_ee_Energy;   //!
   TBranch        *b_vSum_ee_Pt;   //!
   TBranch        *b_vSum_ee_Eta;   //!
   TBranch        *b_vSum_ee_Phi;   //!
   TBranch        *b_vSum_ee_Rapidity;   //!
   TBranch        *b_ele_dpt;   //!
   TBranch        *b_ele_deta;   //!
   TBranch        *b_ele_dphi;   //!
   TBranch        *b_ele_acop;   //!
   TBranch        *b_ok_neuexcl;   //!
   TBranch        *b_ok_zdcexcl;   //!
   TBranch        *b_ok_chexcl;   //!
   TBranch        *b_ok_chexcl_extrk;   //!
   TBranch        *b_nExtrk;   //!
   TBranch        *b_ok_trigger;   //!
   TBranch        *b_ok_event;   //!
   TBranch        *b_zdc_energy_pos;   //!
   TBranch        *b_zdc_energy_neg;   //!
   TBranch        *b_costhetastar;   //!
   TBranch        *b_cos_photon_pair_helicity0;   //!
   TBranch        *b_cos_photon_pair_helicity1;   //!

   ReadqedTree(TTree *tree=0);
   virtual ~ReadqedTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadqedTree_cxx
ReadqedTree::ReadqedTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootFile/ForDielectron/QEDSC__Costheatstar_GenReco_6thMarch.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootFile/ForDielectron/QEDSC__Costheatstar_GenReco_6thMarch.root");
      }
      f->GetObject("output_tree",tree);

   }
   Init(tree);
}

ReadqedTree::~ReadqedTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadqedTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadqedTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ReadqedTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("gen_Pt1", &gen_Pt1, &b_gen_Pt1);
   fChain->SetBranchAddress("gen_Eta1", &gen_Eta1, &b_gen_Eta1);
   fChain->SetBranchAddress("gen_Phi1", &gen_Phi1, &b_gen_Phi1);
   fChain->SetBranchAddress("gen_Pt2", &gen_Pt2, &b_gen_Pt2);
   fChain->SetBranchAddress("gen_Eta2", &gen_Eta2, &b_gen_Eta2);
   fChain->SetBranchAddress("gen_Phi2", &gen_Phi2, &b_gen_Phi2);
   fChain->SetBranchAddress("gen_diEle_M", &gen_diEle_M, &b_gen_diEle_M);
   fChain->SetBranchAddress("gen_diEle_Pt", &gen_diEle_Pt, &b_gen_diEle_Pt);
   fChain->SetBranchAddress("gen_diEle_Rapidity", &gen_diEle_Rapidity, &b_gen_diEle_Rapidity);
   fChain->SetBranchAddress("gen_diEle_Costhetastar", &gen_diEle_Costhetastar, &b_gen_diEle_Costhetastar);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("evtnb", &evtnb, &b_evtnb);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge_1", &eleCharge_1, &b_eleCharge_1);
   fChain->SetBranchAddress("elePt_1", &elePt_1, &b_elePt_1);
   fChain->SetBranchAddress("eleEta_1", &eleEta_1, &b_eleEta_1);
   fChain->SetBranchAddress("elePhi_1", &elePhi_1, &b_elePhi_1);
   fChain->SetBranchAddress("eleSCEt_1", &eleSCEt_1, &b_eleSCEt_1);
   fChain->SetBranchAddress("eleSCEta_1", &eleSCEta_1, &b_eleSCEta_1);
   fChain->SetBranchAddress("eleSCPhi_1", &eleSCPhi_1, &b_eleSCPhi_1);
   fChain->SetBranchAddress("eleEoverP_1", &eleEoverP_1, &b_eleEoverP_1);
   fChain->SetBranchAddress("eleCharge_2", &eleCharge_2, &b_eleCharge_2);
   fChain->SetBranchAddress("elePt_2", &elePt_2, &b_elePt_2);
   fChain->SetBranchAddress("eleEta_2", &eleEta_2, &b_eleEta_2);
   fChain->SetBranchAddress("elePhi_2", &elePhi_2, &b_elePhi_2);
   fChain->SetBranchAddress("eleSCEt_2", &eleSCEt_2, &b_eleSCEt_2);
   fChain->SetBranchAddress("eleSCEta_2", &eleSCEta_2, &b_eleSCEta_2);
   fChain->SetBranchAddress("eleSCPhi_2", &eleSCPhi_2, &b_eleSCPhi_2);
   fChain->SetBranchAddress("eleEoverP_2", &eleEoverP_2, &b_eleEoverP_2);
   fChain->SetBranchAddress("vSum_M", &vSum_M, &b_vSum_ee_M);
   fChain->SetBranchAddress("vSum_Energy", &vSum_Energy, &b_vSum_ee_Energy);
   fChain->SetBranchAddress("vSum_Pt", &vSum_Pt, &b_vSum_ee_Pt);
   fChain->SetBranchAddress("vSum_Eta", &vSum_Eta, &b_vSum_ee_Eta);
   fChain->SetBranchAddress("vSum_Phi", &vSum_Phi, &b_vSum_ee_Phi);
   fChain->SetBranchAddress("vSum_Rapidity", &vSum_Rapidity, &b_vSum_ee_Rapidity);
   fChain->SetBranchAddress("ele_dpt", &ele_dpt, &b_ele_dpt);
   fChain->SetBranchAddress("ele_deta", &ele_deta, &b_ele_deta);
   fChain->SetBranchAddress("ele_dphi", &ele_dphi, &b_ele_dphi);
   fChain->SetBranchAddress("ele_acop", &ele_acop, &b_ele_acop);
   fChain->SetBranchAddress("ok_neuexcl", &ok_neuexcl, &b_ok_neuexcl);
   fChain->SetBranchAddress("ok_zdcexcl", &ok_zdcexcl, &b_ok_zdcexcl);
   fChain->SetBranchAddress("ok_chexcl", &ok_chexcl, &b_ok_chexcl);
   fChain->SetBranchAddress("ok_chexcl_extrk", &ok_chexcl_extrk, &b_ok_chexcl_extrk);
   fChain->SetBranchAddress("nExtrk", &nExtrk, &b_nExtrk);
   fChain->SetBranchAddress("ok_trigger", &ok_trigger, &b_ok_trigger);
   fChain->SetBranchAddress("ok_event", &ok_event, &b_ok_event);
   fChain->SetBranchAddress("zdc_energy_pos", &zdc_energy_pos, &b_zdc_energy_pos);
   fChain->SetBranchAddress("zdc_energy_neg", &zdc_energy_neg, &b_zdc_energy_neg);
   fChain->SetBranchAddress("costhetastar", &costhetastar, &b_costhetastar);
   fChain->SetBranchAddress("cos_photon_pair_helicity0", &cos_photon_pair_helicity0, &b_cos_photon_pair_helicity0);
   fChain->SetBranchAddress("cos_photon_pair_helicity1", &cos_photon_pair_helicity1, &b_cos_photon_pair_helicity1);
   Notify();
}

Bool_t ReadqedTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadqedTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadqedTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ReadqedTree_cxx
