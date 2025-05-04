#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
#include <fstream>
#include <sstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TROOT.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TStyle.h"
#include "TPaveStats.h"
#include "TText.h"
#include "THStack.h"
#include "TCut.h"
#include "TAxis.h"
#include "TChain.h"
#include "Hist.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"


#endif
#include "ReadTauTauTree.C"
//#include "ReadTree_TauTau_gammaUPC.C"
#include "tnp_muon_UPC_PbPb2018.h"

const int nSample = 1;
//const char *sample[nSample]={"Data","QED_SC","QED_SL"};
const char *Sample[nSample]={"TauTau"};

void make_canvas(TCanvas *&);
void make_canvas_ratio(TCanvas *&);
void make_hist(TH1D *&, Color_t , int );
void make_hist_ratio(TH1D *&,Color_t, int) ;
void make_Tpad(TPad *& pad,double,double);
double Width =700;
double Height = 800;
const double luminosity       = 1700;

const double nEventsGeneratedGUPC = 5270000; // 
const double xsecGeneratedGUPC    = 1060.8; // ub
double scaleFactors = 0.943*0.961237*0.9*0.99; //
const double nEventsGeneratedSC = 1000000;
const double xsecGeneratedSC = 570;

double lumiNormGUPC = xsecGeneratedGUPC*luminosity*scaleFactors/nEventsGeneratedGUPC;
double lumiNormSC = xsecGeneratedSC*luminosity*scaleFactors/nEventsGeneratedSC;
double crossSectionUPCGen = 564.89;
double generatedEvtUPCGen = 10000000;
double Other_UPCGen = (luminosity * crossSectionUPCGen*scaleFactors)/generatedEvtUPCGen;

void printOutput(TH1D *hist);

const double wt[nSample] = {lumiNormGUPC};
float getAcoBinSF(float acop);


void drawText(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kMagenta);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
const int nMassbins=4;
double Massbin[nMassbins]={5.0,8.0,11.0,16.0};
const int nMassbin= sizeof(Massbin)/sizeof(double) - 1;

const int nPtbins=18;

double Ptbin[nPtbins]={0.0,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,20.5};
//double Ptbin[nPtbins]={0.0,2.5,4.5,5.5,22.5};
const int nPtbin= sizeof(Ptbin)/sizeof(double) - 1;

TH1D* subtractBkg(TH1D *hdata, TH1D *hmumufsr, const char* name, int nbins, float xmin, float xmax);
TH1D* subtractInvBkg(TH1D *hdata, TH1D *hmumufsr, const char* name);
TH1D* getInvXSecHist(TH1D *hist, const char* name, float lumi);
TH1D* getXSecHist(TH1D *hist, const char* name, float lumi);
TCanvas* PlotHistsAndRatio(TCanvas* , TH1D* , TH1D*, TH1D*, double , double , double , double , double , double , const char *, bool); 
void makeUnfold_TauTau(){
  int itr = 1;
  //gROOT->LoadMacro("CMS_lumi.C");
  bool outOfFrame    = false;
  TH1::SetDefaultSumw2();
  gStyle->SetOptStat(0);
//  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
  TFile *outf= new TFile(Form("Test_TauTau.root"),"recreate");  
  TChain *TauTau[nSample];
  TH1D* hMu_Pt[nSample], *hEleMu_Mass[nSample], *hEleMu_vSumPt[nSample], *hAcoplanarity[nSample]; TH1F *hSF[nSample];
  TH1D* hTau_GenPt[nSample], *hMu_GenPt[nSample], *hEleMu_Gen_Mass[nSample], *hEleMu_Gen_vSumPt[nSample];
  TH2D* hEleMu_RecoGenPt[nSample],  *hEleMu_RecoGenInvmass[nSample], *hEleMu_RecoGenvSumPt[nSample],  *hEleMu_RecoGenTauPt[nSample], *hMu_PtEta[nSample];
  cout << "define hists"<< endl;
  int i =0;
  //for (int i = 0; i < nSample; i++){
    TauTau[i] = new TChain("output_tree");
//TauTau[i]->Add("RootFile/For_TauTau/aTau_0_ChFF_TUnfold_RunningHFThreshold6.4GeV_Crack_HEM_applied.root");
 //TauTau[i]->Add("RootFile/For_TauTau/gammaUPC_ChFF_Unfold_RunningHFThreshold6.4GeV_Crack_HEM_applied.root"); 
 TauTau[i]->Add("RootFile/For_TauTau/aTau0_RunningHFThreshold90Percent_Crack_HEM_applied.root"); 
 //TauTau[i]->Add("RootFile/For_TauTau/EleMu_MuMuFSR_TUnfold_ChFF.root");

    hMu_Pt[i]    = new TH1D(Form("hMu_Pt%s", Sample[i]),"",nPtbin,Ptbin);
    hTau_GenPt[i]         = new TH1D(Form("hTau_GenPt%s", Sample[i]),"",nPtbin,Ptbin);
    hMu_GenPt[i]   = new TH1D(Form("hMu_GenPt%s",Sample[i]),"",nPtbin,Ptbin);
    hEleMu_RecoGenPt[i]         = new TH2D(Form("hEleMu_RecoGenPt%s", Sample[i]),"",nPtbin,Ptbin,nPtbin,Ptbin);
    hMu_PtEta[i]         = new TH2D(Form("hMu_PtEta%s", Sample[i]),"",nPtbin,Ptbin,nPtbin,Ptbin);
   
    TH1D* hUnfoldMuPt  = new TH1D("hUnfoldMuPt", "Pt unfo MC",nPtbin,Ptbin); 
    TH1D* hUnfold_MuPt_data  = new TH1D("hUnfold_MuPt_data", "Pt unfo Data",nPtbin,Ptbin);  
    RooUnfoldResponse responseMuPt (hMu_Pt[i], hMu_GenPt[i], hEleMu_RecoGenPt[i]);
  
    cout << "file " << TauTau[i]->GetEntries()  << endl;
    ReadTauTauTree  TauTauR(TauTau[i]);
    //ReadTree_TauTau_gammaUPC TauTauR(TauTau[i]);
    TauTauR.fChain->SetBranchStatus("*",1);
    if (TauTauR.fChain == 0) return;
    Long64_t nentries = TauTauR.fChain->GetEntriesFast();
    cout << TauTau[i]->GetName() << "    " << nentries << endl;
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry_evt_R = TauTauR.LoadTree(jentry);
      if (ientry_evt_R < 0) break;
      nb = TauTauR.fChain->GetEntry(jentry);   nbytes += nb;   

    //   // "===================      Fill gen Info here ======================     
   //if(TauTauR.Gen_Mu_Pt < 0.0 || abs(TauTauR.Gen_Mu_Eta) > 3.0 || TauTauR.Gen_Ele_Pt < 0.0 || abs(TauTauR.Gen_Ele_Eta) > 3.0 ) continue;
 
   if(TauTauR.EleMuEvents == 0) continue;
   
   if(TauTauR.Tau_Pt2 <= 1.0 || abs(TauTauR.Tau_Eta2) >= 3 ) continue;
  
  // if((abs(TauTauR.Gen_Mu_Eta) < 1.2) && (TauTauR.Gen_Mu_Pt < 3.5)) continue;
    //   else if(((abs(TauTauR.Gen_Mu_Eta) > 1.2) && abs(TauTauR.Gen_Mu_Eta) < 2.4) && (TauTauR.Gen_Mu_Pt < 1.5)) continue;
  //if(TauTauR.hasMuon == 0) continue;

  double ET = TauTauR.Gen_Mu_Eta;
  double PT = TauTauR.Gen_Mu_Pt; 
    // double W = lumiNormGUPC;
  //double lumiNormUPCGen1 = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0) * Other_UPCGen;
   //double W = (466.087*429.51)/10000000;  
  //double W = tnp_weight_softid_upc_pbpb( PT, ET, 0) * tnp_weight_trigger_upc_pbpb( PT, ET, 0));
  double W = 1;
      //hTau_GenPt[i]->Fill(TauTauR.Tau_Pt2,W);
hMu_GenPt[i]->Fill(TauTauR.Gen_Mu_Pt,W);
hMu_PtEta[i]->Fill(TauTauR.Gen_Mu_Pt, TauTauR.Gen_Ele_Eta);
 
      // ===================   Apply analysis cuts on reco MC ======================           
bool isPassing = TauTauR.ok_chexcl_goodtracks == 1 && TauTauR.ok_neuexcl == 1 && TauTauR.ok_neuexcl_HFonly == 1 
   && abs(TauTauR.Mu_Eta) < 2.4
&& TauTauR.EleMu_acop >  0.01 &&  TauTauR.EleMu_acop < 0.15 ;
//bool isPassing = TauTauR.Check ==1;
double ET2 = TauTauR.Mu_Eta;
double PT2 = TauTauR.Mu_Pt;
    //  double Weight3 = tnp_weight_softid_upc_pbpb( PT2, ET2, 0);
    //  double Weight4 = tnp_weight_trigger_upc_pbpb( PT2, ET2, 0);
   //  double W2 = 1.0 ;
//double lumiNormUPCGen2 = tnp_weight_softid_upc_pbpb( PT2, ET2, 0) * tnp_weight_trigger_upc_pbpb( PT2, ET2, 0) * scaleFactors;
double W2 = tnp_weight_softid_upc_pbpb( PT2, ET2, 0) * tnp_weight_trigger_upc_pbpb( PT2, ET2, 0)*1;
//double W2 =1;		
 if(TauTauR.ok_trigger == 1) { //trigger 
	if(isPassing){    	
	hMu_Pt[i]->Fill(TauTauR.Mu_Pt,W2);
  //hEleMu_RecoGenTauPt[i]->Fill(TauTauR.Mu_Pt,TauTauR.Tau_Pt2,W2);
  hEleMu_RecoGenPt[i]->Fill(TauTauR.Mu_Pt,TauTauR.Gen_Mu_Pt,W2);
	
  //responseTauPt.Fill(TauTauR.Mu_Pt,TauTauR.Tau_Pt2,W2);
	responseMuPt.Fill(TauTauR.Mu_Pt,TauTauR.Gen_Mu_Pt,W2);
	} //isPassing analysis cuts
	else {
	//responseTauPt.Miss(TauTauR.Tau_Pt2,W2);
	responseMuPt.Miss(TauTauR.Gen_Mu_Pt,W2);
	}
      } //trigger
	else {
	//responseTauPt.Miss(TauTauR.Tau_Pt2,W2);
	responseMuPt.Miss(TauTauR.Gen_Mu_Pt,W2);
	}   
    } //entry
  
/*************************************************************/
TDecompSVD *svd= new TDecompSVD (responseMuPt.Mresponse());
auto singular_values = svd->GetSig();
//cout << "Test1" << endl;
svd->Print();
double singular_value_min;
      for (int i = 0; i < 5 ; i++)
      {
        //cout << "Test2" << endl;
        cout << singular_values[i] << endl;
      }

//cout << "condition number: "  << singular_values[0]/singular_value_min << endl;
  cout<< "======================================Response matrix========================="<<endl;
  // Closure test using GammaUPC MC
 // TFile *f2=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/For_TauTau/EleMu_histos_gUPC_signal.root","r");
  TFile *f2=new TFile("/Users/pranatijana/Documents/Unfolding/NewUnfold/RooUnfold/examples/RootFile/For_TauTau/EleMu_histos_gUPC_signal_0-1.5-2.5-3.5-4.5-5.5-22.5GeV.root","r");
  f2->cd("signal");
  TH1D *hMu_Pt_gUPC = (TH1D*)gDirectory->Get("hmu_pt");
  
  TCanvas* ccc = new TCanvas("cc","MuPt",50,50,600,670);
  hMu_Pt_gUPC->GetXaxis()->SetTitle("#mu pT");
  hMu_Pt_gUPC->Draw("hist");
  //ccc->Print("Plot/Plots_For_Paper_TauTau/Mu_Pt_gUPC.png");

RooUnfoldInvert  unfold(&responseMuPt, hMu_Pt[0]);
//RooUnfoldInvert  unfold(&responseMuPt, hMu_Pt_gUPC);
unfold.IncludeSystematics(1);
hUnfoldMuPt= (TH1D*) unfold.Hunfold()->Clone("hUnfoldMuPt");

cout << "Test1" << endl;
unfold.PrintTable (cout,hMu_GenPt[0]);

TCanvas* cc = new TCanvas("cc","MuPt",50,50,600,670);
hEleMu_RecoGenPt[0]->GetXaxis()->SetTitle("#mu pT");
cc->SetRightMargin(0.15);
cc->SetBottomMargin(0.15); 
hEleMu_RecoGenPt[0]->Draw("colz");
cc->Print("Plot/Plots_For_Paper_TauTau/Mu_Pt_RecoGen.png");

TCanvas* cc1 = new TCanvas("cc1","GenMuPtEta",50,50,600,670);
hMu_PtEta[0]->GetXaxis()->SetTitle("Gen #mu pT");
hMu_PtEta[0]->GetYaxis()->SetTitle(" Gen electron pT");
cc1->SetRightMargin(0.18);
cc1->SetRightMargin(0.15);
hMu_PtEta[0]->Draw("colz");
cc1->Print("Plot/Plots_For_Paper_TauTau/Mu_PtEtaGen.png");
          
TCanvas *c2 = new TCanvas();

TMatrixD responseMatrix3 = responseMuPt.Mresponse();

// Create a TH2D histogram to display the matrix
       TH2D *hist2D = new TH2D("responseMatrix", "Response Matrix (#mu)",
                        responseMatrix3.GetNcols(), 0, responseMatrix3.GetNcols(),
                        responseMatrix3.GetNrows(), 0, responseMatrix3.GetNrows());
// Fill the TH2D histogram with the matrix values
double total = 0.0;
for (int i = 0; i < responseMatrix3.GetNcols(); i++) {
    for (int j = 0; j < responseMatrix3.GetNrows(); j++) {
        total += responseMatrix3(i, j);
    }
}
for (int i = 0; i < responseMatrix3.GetNcols(); i++) {
    for (int j = 0; j < responseMatrix3.GetNrows(); j++) {
        double normVal = responseMatrix3(i, j) / total;
        hist2D->SetBinContent(i + 1, j + 1, responseMatrix3(i, j));
        //hist2D->SetBinContent(i + 1, j + 1, normVal);
   // cout << "Bincontent:" << hist2D->GetBinContent(i + 1, j + 1, responseMatrix3(i, j)); 
   }
}
hist2D->GetXaxis()->SetTitle("Reco bin index");
hist2D->GetYaxis()->SetTitle("Gen bin index");
c2->SetRightMargin(0.15);
c2->SetBottomMargin(0.15); 
hist2D->Draw("textcolz");
c2->SaveAs("Plot/Plots_For_Paper_TauTau/MuPtMatrix_wonormalized.pdf");

TCanvas*c5 =  new TCanvas();
TMatrixD  responseMatrix5 = responseMuPt.Mresponse(true);
responseMatrix5.Draw("textcolz");
 c5->SaveAs("Plot/Plots_For_Paper_TauTau/CovMatrixmu_UPCGen.png"); 

  /***************************  Start reading data and subtract background contributions ********************************/
 //TFile *f1=new TFile("RootFile/For_TauTau/EleMu_histos.root","r");
//TFile *f1=new TFile("RootFile/For_TauTau/EleMu_histos_data_0-1.5-2.5-3.5-4.5-5.5-22.5GeV.root","r");
TFile *f1=new TFile("RootFile/For_TauTau/EleMu_histos_gUPC_signal_0-1.5-2.5-3.5-4.5-5.5-6.5-7.5-8.5-9.5-10.5-11.5-12.5-13.5-14.5-15.5-16.5-20.5GeV.root","r");
 //TFile *f1=new TFile("RootFile/For_TauTau/EleMu_histos_Data_0-2.5-4.5-5.5-22.5GeV.root","r");
 f1->cd("data");
 TH1D *hMu_Pt_data    = (TH1D*)gDirectory->Get("hmu_pt");
 f1->cd("mumuFSR");
 TH1D *hMu_Pt_mumufsr    = (TH1D*)gDirectory->Get("hmu_pt");
 f1->cd("mumuFSR");
 TH1D* hMu_Pt_data_bkgSub       = subtractInvBkg(hMu_Pt_data, hMu_Pt_mumufsr, "hMu_Pt_data_bkgSub" );


// "===================     Unfold data ======================     
 //TH1D* hMu_Pt_data_bkgSub_rebinned = (TH1D*)hMu_Pt_data_bkgSub->Rebin(nPtbin, "hMu_Pt_data_bkgSub_rebinned", Ptbin);
  

  RooUnfoldInvert  unfold_data_MuPt(&responseMuPt, hMu_Pt_data_bkgSub); // unfold data

  unfold_data_MuPt.IncludeSystematics(1);
  hUnfold_MuPt_data= (TH1D*) unfold_data_MuPt.Hunfold()->Clone("hUnfold_MuPt_data");
  unfold_data_MuPt.PrintTable (cout, hMu_Pt_data_bkgSub);


/**************************Normalization*******************************/


/***********************************************************/
 // cout << " Pt Gen Tau" <<   " :" << hTau_GenPt[0]->Integral() << endl;  
  cout << "Pt Data bkgSub:" << ":"  << hMu_Pt_data_bkgSub->Integral() << endl;
  cout << "Pt Reco:" << ":"  << hMu_Pt[0]->Integral() << endl;
  cout << "Pt unfold data:" << ":"  << hUnfold_MuPt_data->Integral() << endl;
 // cout << "Pt Data bkgSub:" << ":"  << hMu_Pt_data_bkgSub->Integral() << endl;  
  
/****************************************Draw**********************************/
 /************************************TauPt*****************************/
  
/****************************************************************/
TH1D* hMu_GenPt_xSec              = new TH1D("hMu_GenPt_xSec", "Pt gen MC",nPtbin,Ptbin);
TH1D* hUnfoldMuPt_xSec            = new TH1D("hUnfoldMuPt_xSec", "Pt Unfo MC",nPtbin,Ptbin);
TH1D* hUnfold_MuPt_data_xSec      = new TH1D("hUnfold_MuPt_data_xSec", "Pt Unfo data",nPtbin,Ptbin);
TH1D* hMu_Pt_xSec                 = new TH1D("hMu_Pt_xSec", "Pt reco MC",nPtbin,Ptbin);
TH1D* hMu_Pt_gUPC_xSec            = new TH1D("hMu_Pt_gUPC_xSec", "Pt reco MC gUPC",nPtbin,Ptbin);
TH1D* hMu_Pt_data_bkgSub_xSec     = new TH1D("hMu_Pt_data_bkgSub_xSec", "Pt reco data ",nPtbin,Ptbin);

// hMu_GenPt[0]->Scale(564/hMu_GenPt[0]->Integral());//564 for UPCgen//629 for gUPC
// hUnfoldMuPt->Scale(564/hUnfoldMuPt->Integral());
// hUnfold_MuPt_data->Scale(542/hUnfold_MuPt_data->Integral());
// hMu_Pt[0]->Scale(0.075/hMu_Pt[0]->Integral());
// hMu_Pt_gUPC->Scale(1.058/hMu_Pt_gUPC->Integral());
// hMu_Pt_data_bkgSub->Scale(1.058/hMu_Pt_data_bkgSub->Integral());

//hMu_GenPt[0]->Scale(1700*564/20000000);//564 for UPCgen//629 for gUPC
//hUnfoldMuPt->Scale(1700*564/10000000);
//hUnfold_MuPt_data->Scale(542/hUnfold_MuPt_data->Integral());
//hMu_Pt[0]->Scale(0.075/hMu_Pt[0]->Integral());
//hMu_Pt_gUPC->Scale(1.058/hMu_Pt_gUPC->Integral());
//hMu_Pt_data_bkgSub->Scale(1.058/hMu_Pt_data_bkgSub->Integral());

/*
hMu_GenPt[0]->Scale(1/hMu_GenPt[0]->Integral());
hUnfoldMuPt->Scale(1/hUnfoldMuPt->Integral());
hUnfold_MuPt_data->Scale(1/hUnfold_MuPt_data->Integral());
hMu_Pt[0]->Scale(1/hMu_Pt[0]->Integral());
hMu_Pt_gUPC->Scale(1/hMu_Pt_gUPC->Integral());
hMu_Pt_data_bkgSub->Scale(1/hMu_Pt_data_bkgSub->Integral());
*/
// hMu_GenPt_xSec = getInvXSecHist(hMu_GenPt[0], "hMu_GenPt_xSec ",1);
// hUnfoldMuPt_xSec = getInvXSecHist(hUnfoldMuPt, "hUnfoldMuPt_xSec ",1);
// hUnfold_MuPt_data_xSec = getInvXSecHist(hUnfold_MuPt_data, "hUnfold_MuPt_data_xSec ",1);
// hMu_Pt_xSec = getInvXSecHist(hMu_Pt[0], "hMu_Pt_xSec ",1);
// hMu_Pt_gUPC_xSec = getInvXSecHist(hMu_Pt_gUPC, "hMu_Pt_gUPC_xSec ",1);
// hMu_Pt_data_bkgSub_xSec = getInvXSecHist(hMu_Pt_data_bkgSub, "hMu_Pt_data_bkgSub_xSec ",1);

hMu_GenPt_xSec = getXSecHist(hMu_GenPt[0], "hMu_GenPt_xSec ",1700*0.0612);
hUnfoldMuPt_xSec = getInvXSecHist(hUnfoldMuPt, "hUnfoldMuPt_xSec ",1700);
hUnfold_MuPt_data_xSec = getInvXSecHist(hUnfold_MuPt_data, "hUnfold_MuPt_data_xSec ",1700*0.0612);
hMu_Pt_xSec = getInvXSecHist(hMu_Pt[0], "hMu_Pt_xSec ",1700);
hMu_Pt_gUPC_xSec = getInvXSecHist(hMu_Pt_gUPC, "hMu_Pt_gUPC_xSec ",1700);
hMu_Pt_data_bkgSub_xSec = getInvXSecHist(hMu_Pt_data_bkgSub, "hMu_Pt_data_bkgSub_xSec ",1700);

/****************************************************************/
  cout << "After scaling:" << endl;
  cout << "Data" << endl;
  printOutput(hMu_Pt_data);
  cout << "MuMuFSR" << endl;
  printOutput(hMu_Pt_mumufsr);
  cout << "Data wo bkg" << endl;
  printOutput(hMu_Pt_data_bkgSub);
  cout << "Gen MC" << endl;
  printOutput(hMu_GenPt[0]);
  cout << "Unfolded data" << endl;
  printOutput(hUnfold_MuPt_data);
  cout << "Unfolded MC gUPC" << endl;
  printOutput(hUnfoldMuPt);
  cout << "Gen Xsec norm by bin:" << endl;
  printOutput(hMu_GenPt_xSec);
  cout << "Unfolded Data Xsec norm by bin:" << endl;
  printOutput(hUnfold_MuPt_data_xSec);
  cout << "Unfolded MC Xsec norm by bin:" << endl;
  printOutput(hUnfoldMuPt_xSec);

  outf->cd();
  hMu_Pt[0]->Write();
  hUnfoldMuPt->Write();
  hUnfold_MuPt_data->Write();
  hMu_GenPt[0]->Write();
  hMu_Pt_data_bkgSub->Write();
  hEleMu_RecoGenPt[0]->Write();
  hMu_Pt_gUPC->Write();
  hMu_GenPt_xSec->Write();
  hUnfoldMuPt_xSec->Write();
  hUnfold_MuPt_data_xSec->Write();
  hMu_Pt_xSec->Write();
  hMu_Pt_gUPC_xSec->Write();
  hMu_Pt_data_bkgSub_xSec->Write();
  outf->Close();

  
}

TH1D* subtractBkg(TH1D *hdata, TH1D *hmumufsr, const char* name, int nbins, float xmin, float xmax){
  TH1D *hans = new TH1D(name,"",nbins,xmin,xmax);
  //TH1D* hans = (TH1D*)hdata->Clone();
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hdata->GetNbinsX(); i++) {
      double binCont  = hdata->GetBinContent(i)-hmumufsr->GetBinContent(i);
      double binError = sqrt(pow(hdata->GetBinError(i),2) + pow(hmumufsr->GetBinError(i),2) );
      hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError);
      //cout << "******bin:" << i << " " << hdata->GetBinContent(i) << " " << hqed->GetBinContent(i) << " " << hcep->GetBinContent(i) << "   " <<  hans->GetBinContent(i) << "+/-" <<  binError<< endl;
      //double err=0;
      //for (int ivar=1; ivar<14; ivar++) err += pow(hvari[ivar]->GetBinContent(i)-hvari[0]->GetBinContent(i),2);
      //hans->SetBinError(i,sqrt(err+pow(hans->GetBinError(i),2)+pow(glob_syst,2)));
   }
   return hans;
}

TH1D* subtractInvBkg(TH1D *hdata, TH1D *hmumufsr, const char* name) {
    // Clone the data histogram to create the result histogram
    TH1D *hans = (TH1D*)hdata->Clone(name);
    // Add the histograms with unequal binning
    hdata->Add(hmumufsr, -1);
    // Set the result histogram values based on the modified data histogram
    //  
   for (int i = 1; i <= hdata->GetNbinsX(); i++) {
        double binCont  = hdata->GetBinContent(i);
        double binError = sqrt(pow(hdata->GetBinError(i),2) + pow(hmumufsr->GetBinError(i),2));
        //double binError = hdata->GetBinError(i);
        // If data bin content is zero, set the hans bin content to 0
        if (binCont == 0) {
            hans->SetBinContent(i, 0);
           hans->SetBinError(i, 0);
        } else {
            hans->SetBinContent(i, binCont);
            hans->SetBinError(i, binError);
        }
    }
    // Reset the original data histogram
    hdata->Add(hmumufsr);
    return hans;
}

TH1D* getInvXSecHist(TH1D *hist, const char* name, float luminosity){
  TH1D *hans = new TH1D(name,"",nPtbin,Ptbin);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {
    double binCont  = hist->GetBinContent(i)/(luminosity*hist->GetBinWidth(i));
    //double binCont  = hist->GetBinContent(i)/(luminosity);
    double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));
     // double binError = hist->GetBinError(i)/hist->GetBinWidth(i);
      if(binCont <=0){
       hans->SetBinContent(i, 0);
       hans->SetBinError(i,0);
       }
	else{
     hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError); 
    }
      cout << "binerror:" << hans->GetBinError(4) << endl;
      cout << "******bin:" << i << " " << hans->GetBinContent(i) << "+/-" <<  binError<< ": Real bin content: " << hist->GetBinContent(i) <<": real error:" << hist->GetBinError(i) << ":Bin center:" <<  hans->GetXaxis()->GetBinCenter(i) <<  endl;
}
return hans;
}

TH1D* getXSecHist(TH1D *hist, const char* name, float luminosity){
  TH1D *hans = new TH1D(name,"",nPtbin,Ptbin);
   cout << " Name of histogram :" << name << endl;
   for (int i=1; i<=hist->GetNbinsX(); i++) {
    double binCont  = hist->GetBinContent(i)*1700*564/(10000000*luminosity*hist->GetBinWidth(i));
    double binError = binCont * sqrt(pow(hist->GetBinError(i)/hist->GetBinContent(i),2));

      if(binCont <=0){
       hans->SetBinContent(i, 0);
       hans->SetBinError(i,0);
       }
	else{
     hans->SetBinContent(i, binCont);
      hans->SetBinError(i, binError); 
    }
      cout << "binerror:" << hans->GetBinError(4) << endl;
      cout << "******bin:" << i << " " << hans->GetBinContent(i) << "+/-" <<  binError<< ": Real bin content: " << hist->GetBinContent(i) <<": real error:" << hist->GetBinError(i) << ":Bin center:" <<  hans->GetXaxis()->GetBinCenter(i) <<  endl;
}
return hans;
}

void printOutput(TH1D *hist){
 for (int i=1; i<=hist->GetNbinsX(); i++) {
      cout << "******bin:" << i << " " << "bin width:" <<  hist->GetBinWidth(i) << " "<<hist->GetBinContent(i) << "+/-" <<  hist->GetBinError(i) << endl; 
    }
}
void make_canvas(TCanvas *& canvas){

  int W = 600;
  int H = 670;
  float T = 0.08;
  float B = 0.14; 
  float L = 0.18;
  float R = 0.04;

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L );
  canvas->SetRightMargin( R );
  canvas->SetTopMargin( T );
  canvas->SetBottomMargin( B );
  canvas->SetTickx(0);
  canvas->SetTicky(0);
}

void make_canvas_ratio(TCanvas *& can){

  int W = 700;
  int H = 600;
  float T = 0.08;
  float B = 0.14; 
  float L = 0.14;
  float R = 0.04;

  can->SetFillColor(0);
  can->SetBorderMode(0);
  can->SetFrameFillStyle(0);
  can->SetFrameBorderMode(0);
  can->SetLeftMargin( L );
  can->SetRightMargin( R );
  can->SetTopMargin( T );
  can->SetBottomMargin( B );
  can->SetTickx(0);
  can->SetTicky(0);
}


void make_hist(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.05, "XYZ");
 
  hist->SetMarkerColor(kcolor);
  hist->SetLineColor(kcolor);
  hist->SetMarkerStyle(kstyle);
 
}

void make_hist_ratio(TH1D *& hist, Color_t kcolor, int kstyle){
  // For the axis labels:
  hist->SetLabelColor(1, "XYZ");
  hist->SetLabelFont(42, "XYZ");
  hist->SetLabelOffset(0.007, "XYZ");
  hist->SetLabelSize(0.1, "XYZ");

  // For the axis titles:
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleSize(0.11, "XYZ");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(0.4);
  hist->SetMarkerColor(kcolor);
  hist->SetFillColor(kcolor);
  hist->SetMarkerStyle(kstyle);
}

void make_Tpad(TPad *& pad, double topmargin, double bottommargin){
pad->SetTopMargin(topmargin);
pad->SetBottomMargin(bottommargin);

}

