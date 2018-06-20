//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 23 15:51:14 2017 by ROOT version 6.08/07
// from TTree tree/tree
// found on file: DYJetsToLL_M-50_NoPU.root
//////////////////////////////////////////////////////////

#ifndef TreeAnalyzer_h
#define TreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.

class TreeAnalyzer {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t           gen1_pdgId;
    Double_t        gen1_pt;
    Double_t        gen1_eta;
    Double_t        gen1_phi;
    Double_t        gen1_m;
    Double_t        gen1_beta;
    Int_t           gen2_pdgId;
    Double_t        gen2_pt;
    Double_t        gen2_eta;
    Double_t        gen2_phi;
    Double_t        gen2_m;
    Double_t        gen2_beta;
    UShort_t        rpcHit_n;
    Bool_t          rpcHit_isBarrel[1000];   //[rpcHit_n]
    Bool_t          rpcHit_isIRPC[1000];   //[rpcHit_n]
    Short_t         rpcHit_sector[1000];   //[rpcHit_n]
    Short_t         rpcHit_station[1000];   //[rpcHit_n]
    Short_t         rpcHit_wheel[1000];   //[rpcHit_n]
    Short_t         rpcHit_layer[1000];   //[rpcHit_n]
    Short_t         rpcHit_disk[1000];   //[rpcHit_n]
    Short_t         rpcHit_ring[1000];   //[rpcHit_n]
    Double_t        rpcHit_x[1000];   //[rpcHit_n]
    Double_t        rpcHit_y[1000];   //[rpcHit_n]
    Double_t        rpcHit_z[1000];   //[rpcHit_n]
    Double_t        rpcHit_time[1000];   //[rpcHit_n]
    Double_t        rpcHit_timeErr[1000];   //[rpcHit_n]
    Double_t        rpcHit_lx[1000];   //[rpcHit_n]
    Double_t        rpcHit_ly[1000];   //[rpcHit_n]
    Short_t         rpcHit_bx[1000];   //[rpcHit_n]
    UShort_t        muon_n;
    Double_t        muon_pt[100];   //[muon_n]
    Double_t        muon_eta[100];   //[muon_n]
    Double_t        muon_phi[100];   //[muon_n]
    Short_t         muon_q[100];   //[muon_n]
    Bool_t          muon_isLoose[100];   //[muon_n]
    Bool_t          muon_isTight[100];   //[muon_n]
    Bool_t          muon_isRPC[100];   //[muon_n]
    Double_t        muon_time[100];   //[muon_n]
    Double_t        muon_RPCTime[100];   //[muon_n]
    Double_t        muon_RPCTimeNew[100];   //[muon_n]
    Double_t        muon_RPCBeta[100];   //[muon_n]
    UShort_t        muon_nRPC[100];   //[muon_n]
    UShort_t        muon_nIRPC[100];   //[muon_n]
    Double_t        muon_genDR[100];   //[muon_n]
    Short_t         muon_genPdgId[100];   //[muon_n]
    UShort_t        simDigi1_n; 
    Double_t        simDigi1_x[1000];   //[simDigi1_n]
    Double_t        simDigi1_y[1000];   //[simDigi1_n]
    Double_t        simDigi1_z[1000];   //[simDigi1_n]
    UShort_t        simDigi2_n;
    Double_t        simDigi2_x[1000];   //[simDigi2_n]
    Double_t        simDigi2_y[1000];   //[simDigi2_n]
    Double_t        simDigi2_z[1000];   //[simDigi2_n]

    // List of branches
    TBranch        *b_gen1_pdgId;   //!
    TBranch        *b_gen1_pt;   //!
    TBranch        *b_gen1_eta;   //!
    TBranch        *b_gen1_phi;   //!
    TBranch        *b_gen1_m;   //!
    TBranch        *b_gen1_beta;   //!
    TBranch        *b_gen2_pdgId;   //!
    TBranch        *b_gen2_pt;   //!
    TBranch        *b_gen2_eta;   //!
    TBranch        *b_gen2_phi;   //!
    TBranch        *b_gen2_m;   //!
    TBranch        *b_gen2_beta;   //!
    TBranch        *b_simDigi1_n;   //!
    TBranch        *b_simDigi1_isBarrel;   //!
    TBranch        *b_simDigi1_isIRPC;   //!
    TBranch        *b_simDigi1_sector;   //!
    TBranch        *b_simDigi1_station;   //!
    TBranch        *b_simDigi1_wheel;   //!
    TBranch        *b_simDigi1_layer;   //!
    TBranch        *b_simDigi1_disk;   //!
    TBranch        *b_simDigi1_ring;   //!
    TBranch        *b_simDigi1_x;   //!
    TBranch        *b_simDigi1_y;   //!
    TBranch        *b_simDigi1_z;   //!
    TBranch        *b_simDigi1_tof;   //!
    TBranch        *b_simDigi1_t0;   //!
    TBranch        *b_simDigi1_lx;   //!
    TBranch        *b_simDigi1_ly;   //!
    TBranch        *b_simDigi1_bx;   //!
    TBranch        *b_simDigi2_n;   //!
    TBranch        *b_simDigi2_isBarrel;   //!
    TBranch        *b_simDigi2_isIRPC;   //!
    TBranch        *b_simDigi2_sector;   //!
    TBranch        *b_simDigi2_station;   //!
    TBranch        *b_simDigi2_wheel;   //!
    TBranch        *b_simDigi2_layer;   //!
    TBranch        *b_simDigi2_disk;   //!
    TBranch        *b_simDigi2_ring;   //!
    TBranch        *b_simDigi2_x;   //!
    TBranch        *b_simDigi2_y;   //!
    TBranch        *b_simDigi2_z;   //!
    TBranch        *b_simDigi2_tof;   //!
    TBranch        *b_simDigi2_t0;   //!
    TBranch        *b_simDigi2_lx;   //!
    TBranch        *b_simDigi2_ly;   //!
    TBranch        *b_simDigi2_bx;   //!
    TBranch        *b_rpcHit_n;   //!
    TBranch        *b_rpcHit_isBarrel;   //!
    TBranch        *b_rpcHit_isIRPC;   //!
    TBranch        *b_rpcHit_sector;   //!
    TBranch        *b_rpcHit_station;   //!
    TBranch        *b_rpcHit_wheel;   //!
    TBranch        *b_rpcHit_layer;   //!
    TBranch        *b_rpcHit_disk;   //!
    TBranch        *b_rpcHit_ring;   //!
    TBranch        *b_rpcHit_x;   //!
    TBranch        *b_rpcHit_y;   //!
    TBranch        *b_rpcHit_z;   //!
    TBranch        *b_rpcHit_time;   //!
    TBranch        *b_rpcHit_timeErr;   //!
    TBranch        *b_rpcHit_lx;   //!
    TBranch        *b_rpcHit_ly;   //!
    TBranch        *b_rpcHit_bx;   //!
    TBranch        *b_muon_n;   //!
    TBranch        *b_muon_pt;   //!
    TBranch        *b_muon_eta;   //!
    TBranch        *b_muon_phi;   //!
    TBranch        *b_muon_q;   //!
    TBranch        *b_muon_isLoose;   //!
    TBranch        *b_muon_isTight;   //!
    TBranch        *b_muon_isRPC;   //!
    TBranch        *b_muon_time;   //!
    TBranch        *b_muon_RPCTime;   //!
    TBranch        *b_muon_RPCTimeNew;   //!
    TBranch        *b_muon_RPCBeta;   //!
    TBranch        *b_muon_nRPC;   //!
    TBranch        *b_muon_nIRPC;   //!
    TBranch        *b_muon_genDR;   //!
    TBranch        *b_muon_genPdgId;   //!

    TreeAnalyzer(TTree *tree=0);
    virtual ~TreeAnalyzer();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop(TFile* outFile);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    std::vector<std::vector<unsigned>> clusterHitsByGenP4s(const TLorentzVector p4s[]) const;
    std::vector<std::vector<unsigned>> clusterHitsByEtaPhi() const;
    std::vector<double> fitTrackBxConstrained(const std::vector<unsigned>& hits) const;
    std::vector<double> fitTrackSlope(const std::vector<unsigned>& hits) const;

};

#endif

#ifdef TreeAnalyzer_cxx
TreeAnalyzer::TreeAnalyzer(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DYJetsToLL_M-50_NoPU.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("DYJetsToLL_M-50_NoPU.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("DYJetsToLL_M-50_NoPU.root:/HSCPTree");
    dir->GetObject("tree",tree);

  }
  Init(tree);
}

TreeAnalyzer::~TreeAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t TreeAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t TreeAnalyzer::LoadTree(Long64_t entry)
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

void TreeAnalyzer::Init(TTree *tree)
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

  fChain->SetBranchAddress("gen1_pdgId", &gen1_pdgId, &b_gen1_pdgId);
  fChain->SetBranchAddress("gen1_pt", &gen1_pt, &b_gen1_pt);
  fChain->SetBranchAddress("gen1_eta", &gen1_eta, &b_gen1_eta);
  fChain->SetBranchAddress("gen1_phi", &gen1_phi, &b_gen1_phi);
  fChain->SetBranchAddress("gen1_m", &gen1_m, &b_gen1_m);
  fChain->SetBranchAddress("gen1_beta", &gen1_beta, &b_gen1_beta);
  fChain->SetBranchAddress("gen2_pdgId", &gen2_pdgId, &b_gen2_pdgId);
  fChain->SetBranchAddress("gen2_pt", &gen2_pt, &b_gen2_pt);
  fChain->SetBranchAddress("gen2_eta", &gen2_eta, &b_gen2_eta);
  fChain->SetBranchAddress("gen2_phi", &gen2_phi, &b_gen2_phi);
  fChain->SetBranchAddress("gen2_m", &gen2_m, &b_gen2_m);
  fChain->SetBranchAddress("gen2_beta", &gen2_beta, &b_gen2_beta);
  fChain->SetBranchAddress("simDigi1_n", &simDigi1_n, &b_simDigi1_n);
  fChain->SetBranchAddress("simDigi1_isBarrel", simDigi1_isBarrel, &b_simDigi1_isBarrel);
  fChain->SetBranchAddress("simDigi1_isIRPC", simDigi1_isIRPC, &b_simDigi1_isIRPC);
  fChain->SetBranchAddress("simDigi1_sector", simDigi1_sector, &b_simDigi1_sector);
  fChain->SetBranchAddress("simDigi1_station", simDigi1_station, &b_simDigi1_station);
  fChain->SetBranchAddress("simDigi1_wheel", simDigi1_wheel, &b_simDigi1_wheel);
  fChain->SetBranchAddress("simDigi1_layer", simDigi1_layer, &b_simDigi1_layer);
  fChain->SetBranchAddress("simDigi1_disk", simDigi1_disk, &b_simDigi1_disk);
  fChain->SetBranchAddress("simDigi1_ring", simDigi1_ring, &b_simDigi1_ring);
  fChain->SetBranchAddress("simDigi1_x", simDigi1_x, &b_simDigi1_x);
  fChain->SetBranchAddress("simDigi1_y", simDigi1_y, &b_simDigi1_y);
  fChain->SetBranchAddress("simDigi1_z", simDigi1_z, &b_simDigi1_z);
  fChain->SetBranchAddress("simDigi1_tof", simDigi1_tof, &b_simDigi1_tof);
  fChain->SetBranchAddress("simDigi1_t0", simDigi1_t0, &b_simDigi1_t0);
  fChain->SetBranchAddress("simDigi1_lx", simDigi1_lx, &b_simDigi1_lx);
  fChain->SetBranchAddress("simDigi1_ly", simDigi1_ly, &b_simDigi1_ly);
  fChain->SetBranchAddress("simDigi1_bx", simDigi1_bx, &b_simDigi1_bx);
  fChain->SetBranchAddress("simDigi2_n", &simDigi2_n, &b_simDigi2_n);
  fChain->SetBranchAddress("simDigi2_isBarrel", simDigi2_isBarrel, &b_simDigi2_isBarrel);
  fChain->SetBranchAddress("simDigi2_isIRPC", simDigi2_isIRPC, &b_simDigi2_isIRPC);
  fChain->SetBranchAddress("simDigi2_sector", simDigi2_sector, &b_simDigi2_sector);
  fChain->SetBranchAddress("simDigi2_station", simDigi2_station, &b_simDigi2_station);
  fChain->SetBranchAddress("simDigi2_wheel", simDigi2_wheel, &b_simDigi2_wheel);
  fChain->SetBranchAddress("simDigi2_layer", simDigi2_layer, &b_simDigi2_layer);
  fChain->SetBranchAddress("simDigi2_disk", simDigi2_disk, &b_simDigi2_disk);
  fChain->SetBranchAddress("simDigi2_ring", simDigi2_ring, &b_simDigi2_ring);
  fChain->SetBranchAddress("simDigi2_x", simDigi2_x, &b_simDigi2_x);
  fChain->SetBranchAddress("simDigi2_y", simDigi2_y, &b_simDigi2_y);
  fChain->SetBranchAddress("simDigi2_z", simDigi2_z, &b_simDigi2_z);
  fChain->SetBranchAddress("simDigi2_tof", simDigi2_tof, &b_simDigi2_tof);
  fChain->SetBranchAddress("simDigi2_t0", simDigi2_t0, &b_simDigi2_t0);
  fChain->SetBranchAddress("simDigi2_lx", simDigi2_lx, &b_simDigi2_lx);
  fChain->SetBranchAddress("simDigi2_ly", simDigi2_ly, &b_simDigi2_ly);
  fChain->SetBranchAddress("simDigi2_bx", simDigi2_bx, &b_simDigi2_bx);
  fChain->SetBranchAddress("rpcHit_n", &rpcHit_n, &b_rpcHit_n);
  fChain->SetBranchAddress("rpcHit_isBarrel", rpcHit_isBarrel, &b_rpcHit_isBarrel);
  fChain->SetBranchAddress("rpcHit_isIRPC", rpcHit_isIRPC, &b_rpcHit_isIRPC);
  fChain->SetBranchAddress("rpcHit_sector", rpcHit_sector, &b_rpcHit_sector);
  fChain->SetBranchAddress("rpcHit_station", rpcHit_station, &b_rpcHit_station);
  fChain->SetBranchAddress("rpcHit_wheel", rpcHit_wheel, &b_rpcHit_wheel);
  fChain->SetBranchAddress("rpcHit_layer", rpcHit_layer, &b_rpcHit_layer);
  fChain->SetBranchAddress("rpcHit_disk", rpcHit_disk, &b_rpcHit_disk);
  fChain->SetBranchAddress("rpcHit_ring", rpcHit_ring, &b_rpcHit_ring);
  fChain->SetBranchAddress("rpcHit_x", rpcHit_x, &b_rpcHit_x);
  fChain->SetBranchAddress("rpcHit_y", rpcHit_y, &b_rpcHit_y);
  fChain->SetBranchAddress("rpcHit_z", rpcHit_z, &b_rpcHit_z);
  fChain->SetBranchAddress("rpcHit_time", rpcHit_time, &b_rpcHit_time);
  fChain->SetBranchAddress("rpcHit_timeErr", rpcHit_timeErr, &b_rpcHit_timeErr);
  fChain->SetBranchAddress("rpcHit_lx", rpcHit_lx, &b_rpcHit_lx);
  fChain->SetBranchAddress("rpcHit_ly", rpcHit_ly, &b_rpcHit_ly);
  fChain->SetBranchAddress("rpcHit_bx", rpcHit_bx, &b_rpcHit_bx);
  fChain->SetBranchAddress("muon_n", &muon_n, &b_muon_n);
  fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
  fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
  fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
  fChain->SetBranchAddress("muon_q", muon_q, &b_muon_q);
  fChain->SetBranchAddress("muon_isLoose", muon_isLoose, &b_muon_isLoose);
  fChain->SetBranchAddress("muon_isTight", muon_isTight, &b_muon_isTight);
  fChain->SetBranchAddress("muon_isRPC", muon_isRPC, &b_muon_isRPC);
  fChain->SetBranchAddress("muon_time", muon_time, &b_muon_time);
  fChain->SetBranchAddress("muon_RPCTime", muon_RPCTime, &b_muon_RPCTime);
  fChain->SetBranchAddress("muon_RPCTimeNew", muon_RPCTimeNew, &b_muon_RPCTimeNew);
  fChain->SetBranchAddress("muon_RPCBeta", muon_RPCBeta, &b_muon_RPCBeta);
  fChain->SetBranchAddress("muon_nRPC", muon_nRPC, &b_muon_nRPC);
  fChain->SetBranchAddress("muon_nIRPC", muon_nIRPC, &b_muon_nIRPC);
  fChain->SetBranchAddress("muon_genDR", muon_genDR, &b_muon_genDR);
  fChain->SetBranchAddress("muon_genPdgId", muon_genPdgId, &b_muon_genPdgId);
  Notify();
}

Bool_t TreeAnalyzer::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void TreeAnalyzer::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
#endif // #ifdef TreeAnalyzer_cxx
