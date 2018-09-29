#define TreeAnalyzer_cxx
#include "TreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include <vector>
#include <iostream>
#include <algorithm>

const double muonMass = 0.1057;
const double speedOfLight = 29.979; // 30cm/ns
const int bxLo = -8, nBx = 17;

const unsigned minNhitCluster = 2;
enum class ClusterAlgo { GenMatch, Histogram };
enum class FitAlgo { BxContrained, FitSlope };
//ClusterAlgo clusterAlgo = ClusterAlgo::Histogram;
ClusterAlgo clusterAlgo = ClusterAlgo::GenMatch;
//FitAlgo fitAlgo = FitAlgo::FitSlope;
FitAlgo fitAlgo = FitAlgo::BxContrained;


double deltaPhi(const double phi1, const double phi2)
{
  double dphi = phi1 - phi2;
  if ( dphi < -TMath::TwoPi() ) dphi += TMath::TwoPi();
  else if ( dphi > TMath::TwoPi() ) dphi -= TMath::TwoPi();
  return dphi;
}

struct HEtaPhiBeta
{
  constexpr static int nbinsEta = 25, nbinsPhi = 50, nbinsBeta = 15;
  constexpr static double minEta = -2.5, maxEta = 2.5;
  constexpr static double minPhi = 0, maxPhi = 2*3.14159265358979323846L+1e-9;
  constexpr static double minBeta = -1, maxBeta = 2;
  constexpr static double binwEta  = (maxEta-minEta)/nbinsEta;
  constexpr static double binwPhi  = (maxPhi-minPhi)/nbinsPhi;
  constexpr static double binwBeta = (maxBeta-minBeta)/nbinsBeta;

  unsigned findBin(const double eta, const double phi, const double beta) const
  {
    const int ieta  = std::max(-1, std::min(nbinsEta, int(floor((eta-minEta)/binwEta))));
    const int iphi  = std::max(-1, std::min(nbinsPhi, int(floor((phi-minPhi)/binwPhi))));
    const int ibeta = std::max(-1, std::min(nbinsBeta, int(floor((beta-minBeta)/binwBeta))));

    const unsigned ibin = (ieta+1)
                        + (iphi+1)*(nbinsEta+2);
                        //+ (ibeta+1)*(nbinsEta+2)*(nbinsPhi+2);
    return ibin;
  }

  std::vector<double> findBinLowEdge(unsigned ibin) const
  {
    const int ieta = ibin % (nbinsEta+2) - 1;
    ibin /= nbinsEta+2;
    const int iphi = ibin % (nbinsPhi+2) - 1;
    //ibin /= nbinsPhi+2;
    //const int ibeta = ibin - 1;

    const double eta  = ieta*binwEta+minEta;
    const double phi  = iphi*binwPhi+minPhi;
    //const double beta = ibeta*binwBeta+minBeta;

    std::vector<double> res = {eta, phi};//, beta};
    return res;
  }


  void fill(const double eta, const double phi, const double beta, const unsigned idx)
  {
    const unsigned ibin = findBin(eta, phi, beta);
    auto itr = contents.find(ibin);
    if ( itr == contents.end() ) itr = contents.insert(std::make_pair(ibin, std::vector<unsigned>())).first;
    itr->second.push_back(idx);
  }

  std::map<unsigned, std::vector<unsigned>> contents;
};

std::vector<std::vector<unsigned>> TreeAnalyzer::clusterHitsByGenP4s(const TLorentzVector p4s[]) const
{
  // Cluster hits along the GenParticle's four momentum, just for initial testing.
  std::vector< std::vector<unsigned> > clusters;
  clusters.resize(2);
/*  
  for ( unsigned i=0; i<rpcHit_n; ++i ) {
    const TVector3 pos_rpc(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    const TVector3 pos_gem(gemSegment_x[i], gemSegment_y[i], gemSegment_z[i]);
    double minDR_rpc = 0.3;
    int match = -1;
    for ( unsigned j=0; j<2; ++j ) {
      const double dR_rpc = p4s[j].Vect().DeltaR(pos_rpc);
      if ( dR_rpc < minDR_rpc ) {
        minDR_rpc = dR_rpc;
        match = j;
      }
    }
    if ( match >= 0 ) {
      for ( unsigned k=0; k<2; ++k ) { 
        const double dEta_rpc = p4s[k].Eta()-pos_rpc.Eta();
        const double dPhi_rpc = p4s[k].Vect().DeltaPhi(pos_rpc);
    //    if ( std::abs(dEta) < 0.2 && std::abs(dPhi) < 0.02 ) clusters.at(match).push_back(i);
      }
    clusters.at(match).push_back(i);
    }
  }  
*/
  

  TVector3 pos_rpc;
  TVector3 pos_gem;
  TVector3 pos_csc;
  double minDR_rpc = 0.3;
  std::vector< std::vector<unsigned> > dR_rpc, dEta_rpc, dPhi_rpc, dEta_gem, dPhi_gem, dEta_csc, dPhi_csc;
  for ( auto iRPC=0; iRPC<rpcHit_n; ++iRPC ) {
    pos_rpc.SetXYZ(rpcHit_x[iRPC], rpcHit_y[iRPC], rpcHit_z[iRPC]);
    for ( auto jRPC=0; jRPC<2; ++jRPC ) {
      dR_rpc.at(jRPC).push_back(std::abs(p4s[jRPC].Vect().DeltaR(pos_rpc)));
      dEta_rpc.at(jRPC).push_back(std::abs(p4s[jRPC].Eta() - pos_rpc.Eta()));
      dPhi_rpc.at(jRPC).push_back(std::abs(p4s[jRPC].Vect().DeltaPhi(pos_rpc)));
      //if ( dR_rpc >= minDR_rpc ) continue;
      //if ( dR_gem < minDR_gem && dR_rpc < minDR_rpc && std::abs(dEta) < 0.2 && std::abs(dPhi) < 0.02 ) clusters.at(j).push_back(iRPC);
    }
  }
  for ( auto iGEM = 0; iGEM<gemSegment_n; ++iGEM ) {
    pos_gem.SetXYZ(gemSegment_x[iGEM], gemSegment_y[iGEM], gemSegment_z[iGEM]);
    for ( auto jGEM=0; jGEM<2; ++jGEM ) {
      dEta_gem.at(jGEM).push_back(std::abs(p4s[jGEM].Eta() - pos_gem.Eta()));
      dPhi_gem.at(jGEM).push_back(std::abs(p4s[jGEM].Vect().DeltaPhi(pos_gem)));
      
    }
  }
  for ( auto iCSC = 0; iCSC<cscSegment_n; ++iCSC ) {
    pos_csc.SetXYZ(cscSegment_x[iCSC], cscSegment_y[iCSC], cscSegment_z[iCSC]);
    for ( auto jCSC=0; jCSC<2; ++jCSC ) {
      dEta_csc.at(jCSC).push_back(std::abs(p4s[jCSC].Eta() - pos_csc.Eta()));
      dPhi_csc.at(jCSC).push_back(std::abs(p4s[jCSC].Vect().DeltaPhi(pos_csc)));
    }
  }
  for ( auto event=0; event < rpcHit_n; ++event ) {
    for ( auto j=0; j<2; ++j ) {
      if ( dR_rpc[event][j] >= minDR_rpc ) continue;
      //iRPC region cut
      if ( p4s[j].Eta() >= 1.8 && p4s[j].Eta() < 2.4 && (dEta_rpc[event][j] > 0.03 || dPhi_rpc[event][j] > 0.008) ) continue;
      //RPC and GEM
      if ( dR_rpc[event][j] < minDR_rpc && dPhi_gem[event][j] <= 0.005 && dEta_gem[event][j] <= 0.06 ) clusters.at(j).push_back(event);
      //RPC and CSC
      else if ( dR_rpc[event][j] < minDR_rpc && dPhi_csc[event][j] <= 0.008 && dEta_csc[event][j] <= 0.06 ) clusters.at(j).push_back(event);
    }
  }
  
  return clusters;
}

std::vector<std::vector<unsigned>> TreeAnalyzer::clusterHitsByEtaPhi() const
{
  // Hit clustering in 3D, eta-phi-beta space with vx=vy=vz=0 & bx=0 hypothesis
  std::vector<std::vector<unsigned>> clusters;
  // First step to fill "histograms"
  HEtaPhiBeta h;
  for ( unsigned i=0; i<rpcHit_n; ++i ) {
    const TVector3 pos(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    const double eta = pos.Eta(), phi = pos.Phi();
    const double ct = speedOfLight*rpcHit_time[i];
    const double beta = 1./(1+ct/pos.Mag());

    h.fill(eta, phi, beta, i);
  }

  for ( auto item : h.contents ) {
    if ( item.second.size() < minNhitCluster ) continue;
    clusters.push_back(item.second);
  }

  return clusters;
}

std::vector<double> TreeAnalyzer::fitTrackBxConstrained(const std::vector<unsigned>& hits) const
{
  // result: qual, beta, t0
  std::vector<double> result = {1e9, 0, 0, 0, 0};
  const unsigned n = hits.size();
  if ( n == 0 ) return result;

  double sumBeta = 0;
  double sumTimeDiff2 = 0;
  double firstPhi = 0; // this is necessary to shift all phi's to the new origin. without shifting, +2pi and -2pi will result unphysical large variance
  double sumEta = 0, sumEta2 = 0, sumDphi = 0, sumDphi2 = 0;

  unsigned nValidBeta = 0;
  for ( auto i : hits ) {
    const TVector3 pos(rpcHit_x[i], rpcHit_y[i], rpcHit_z[i]);
    const double ct = speedOfLight*rpcHit_time[i];
    const double ibeta = 1./(1+ct/pos.Mag());
    if ( ibeta > 2 or ibeta < -1 ) continue;
    sumBeta += ibeta;
    sumTimeDiff2 += rpcHit_time[i]*rpcHit_time[i];

    sumEta += pos.Eta();
    sumEta2 += pos.Eta()*pos.Eta();
    if ( nValidBeta == 0 )  firstPhi = pos.Phi();
    else {
      const double dphi = pos.Phi()-firstPhi;
      sumDphi += dphi;
      sumDphi2 += dphi*dphi;
    }


    ++nValidBeta;
  }
  const double dRErr2 = ((sumDphi2-sumDphi*sumDphi/nValidBeta) +
                         (sumEta2-sumEta*sumEta/nValidBeta))/nValidBeta;

  result = {dRErr2, sumBeta/nValidBeta, sumTimeDiff2/nValidBeta};

  return result;
}

std::vector<double> TreeAnalyzer::fitTrackSlope(const std::vector<unsigned>& hits) const
{
  std::vector<double> result = {1e9, 0, 1e9, 0, 0};
  const unsigned n = hits.size();
  if ( n <= 2 ) return result;

  double sx = 0, sy = 0, sxy = 0, sxx = 0, syy = 0;
  for ( auto i : hits ) {
    const double r2 = rpcHit_x[i]*rpcHit_x[i] + rpcHit_y[i]*rpcHit_y[i] + rpcHit_z[i]*rpcHit_z[i];
    const double r = std::sqrt(r2);
    const double t = rpcHit_time[i];
    sx  += r;
    sy  += t;
    sxy += r*t;
    sxx += r2;
    syy += t*t;
  }
  const double ssxy = sxy-sx*sy/n;
  const double ssxx = sxx-sx*sx/n;
  const double ssyy = syy-sy*sy/n;
  const double b = ssxy/ssxx;
  const double s = std::sqrt((ssyy - b*ssxy)/(n-2));
  const double bStdErr = s/std::sqrt(ssxx);

  const double t0 = (sy-b*sx)/n; // is "a" in the original code
  const double t0Err = s*std::sqrt(1./n + sx*sx*ssxx/n/n); // is "aStdErr" in the original code
  const double beta = 1./(b*speedOfLight+1.);
  const double betaErr = speedOfLight*bStdErr/((b*speedOfLight+1)*(b*speedOfLight+1));

  const int nbx = std::round(t0/25.);
  const double bxPull = t0 - nbx*25;
  //const double bxPull = t0/25.;
  const double RelativeSlopeErr = bStdErr/b;
  result = {bxPull, beta, betaErr, t0, t0Err};
  //result = {bxPull, beta, RelativeSlopeErr, t0, t0Err};

  return result;
}

void TreeAnalyzer::Loop(TFile* fout)
{
  int oneHSCP = 0, twoHSCP = 0, nEvent = 0;
  double oneTotEff = 0.0;
  double twoTotEff = 0.0;
  fout->cd();
  TTree* tree = new TTree("tree", "tree");

  TLorentzVector out_gens_p4[2];
  int out_gens_pdgId[2];
  tree->Branch("gen1_p4", "TLorentzVector", &out_gens_p4[0]);
  tree->Branch("gen2_p4", "TLorentzVector", &out_gens_p4[1]);
  tree->Branch("gen1_pdgId", &out_gens_pdgId[0], "gen1_pdgId/I");
  tree->Branch("gen2_pdgId", &out_gens_pdgId[1], "gen2_pdgId/I");

  unsigned out_muons_n;
  TLorentzVector out_muons_p4[3];
  int out_muons_q[3];
  tree->Branch("muon1_p4", "TLorentzVector", &out_muons_p4[0]);
  tree->Branch("muon2_p4", "TLorentzVector", &out_muons_p4[1]);
  tree->Branch("muon3_p4", "TLorentzVector", &out_muons_p4[2]);
  tree->Branch("muon1_q", &out_muons_q[0], "muon1_q/I");
  tree->Branch("muon2_q", &out_muons_q[1], "muon2_q/I");
  tree->Branch("muon3_q", &out_muons_q[2], "muon3_q/I");

  TVector3 out_rpc_p3[2];
  TVector3 out_gem_p3[2];
  TVector3 out_csc_p3[2];
  tree->Branch("rpc1_clust_p3", "TVector3", &out_rpc_p3[0]);
  tree->Branch("rpc2_clust_p3", "TVector3", &out_rpc_p3[1]);
  tree->Branch("gem1_clust_p3", "TVector3", &out_gem_p3[0]);
  tree->Branch("gem2_clust_p3", "TVector3", &out_gem_p3[1]);
  tree->Branch("csc1_clust_p3", "TVector3", &out_csc_p3[0]);
  tree->Branch("csc2_clust_p3", "TVector3", &out_csc_p3[1]);

  float out_fit_quals[2];
  float out_fit_betas[2];
  unsigned out_fit_nhits[2];
  tree->Branch("fit_qual1", &out_fit_quals[0], "fit_qual1/F");
  tree->Branch("fit_qual2", &out_fit_quals[1], "fit_qual2/F");
  tree->Branch("fit_beta1", &out_fit_betas[0], "fit_beta1/F");
  tree->Branch("fit_beta2", &out_fit_betas[1], "fit_beta2/F");
  tree->Branch("fit_nhit1", &out_fit_nhits[0], "fit_nhit1/i");
  tree->Branch("fit_nhit2", &out_fit_nhits[1], "fit_nhit2/i");

  TH1* hDeta_1 = new TH1F("hDeta_1", "positive HSCP deta(gen-sim)", 100, -0.2, 0.2);
  TH1* hDphi_1 = new TH1F("hDphi_1", "positive HSCP dphi", 100, -0.03, 0.03);  
  tree->Branch("hDeta_1", &hDeta_1, "TH1F");
  tree->Branch("hDphi_1", &hDphi_1, "TH1F");

  TH1* hDeta_2 = new TH1F("hDeta_2", "negative HSCP deta(gen-sim)", 100, -0.2, 0.2);
  TH1* hDphi_2 = new TH1F("hDphi_2", "negative HSCP dphi", 100, -0.03, 0.03);
  tree->Branch("hDeta_2", &hDeta_2, "TH1F");
  tree->Branch("hDphi_2", &hDphi_2, "TH1F");

  TH1* hDeta_rpcHit_1 = new TH1F("hDeta_rpcHit_1", "positive HSCP deta(gen-rpcHit)", 100, -0.2, 0.2);
  TH1* hDphi_rpcHit_1 = new TH1F("hDphi_rpcHit_1", "positive HSCP dphi(gen-rpcHit)", 100, -0.02, 0.02);
  tree->Branch("hDeta_rpcHit_1", &hDeta_rpcHit_1, "TH1F");
  tree->Branch("hDphi_rpcHit_1", &hDphi_rpcHit_1, "TH1F");

  TH1* hDeta_rpcHit_2 = new TH1F("hDeta_rpcHit_2", "negative HSCP deta(gen-rpcHit)", 100, -0.2, 0.2);
  TH1* hDphi_rpcHit_2 = new TH1F("hDphi_rpcHit_2", "negative HSCP dphi(gen-rpcHit)", 100, -0.02, 0.02);
  tree->Branch("hDeta_rpcHit_2", &hDeta_rpcHit_2, "TH1F");
  tree->Branch("hDphi_rpcHit_2", &hDphi_rpcHit_2, "TH1F");

  TH1* hDeta_gemSeg_1 = new TH1F("hDeta_gemSeg_1", "positive HSCP deta(gen-gemSegment)", 100, -3.0, 3.0);
  TH1* hDphi_gemSeg_1 = new TH1F("hDphi_gemSeg_1", "positive HSCP dphi(gen-gemSegment)", 100, -0.5, 0.5);
  tree->Branch("hDeta_gemSeg_1", &hDeta_gemSeg_1, "TH1F");
  tree->Branch("hDphi_gemSeg_1", &hDphi_gemSeg_1, "TH1F");

  TH1* hDeta_gemSeg_2 = new TH1F("hDeta_gemSeg_2", "negative HSCP deta(gen-gemSegment)", 100, -3.0, 3.0);
  TH1* hDphi_gemSeg_2 = new TH1F("hDphi_gemSeg_2", "negative HSCP dphi(gen-gemSegment)", 100, -0.5, 0.5);
  tree->Branch("hDeta_gemSeg_2", &hDeta_gemSeg_2, "TH1F");
  tree->Branch("hDphi_gemSeg_2", &hDphi_gemSeg_2, "TH1F");

  TH1* hDeta_cscSeg_1 = new TH1F("hDeta_cscSeg_1", "positive HSCP deta(gen-cscSegment)", 100, -3.0, 3.0);
  TH1* hDphi_cscSeg_1 = new TH1F("hDphi_cscSeg_1", "positive HSCP dphi(gen-cscSegment)", 100, -0.5, 0.5);
  tree->Branch("hDeta_cscSeg_1", &hDeta_cscSeg_1, "TH1F");
  tree->Branch("hDphi_cscSeg_1", &hDphi_cscSeg_1, "TH1F");

  TH1* hDeta_cscSeg_2 = new TH1F("hDeta_cscSeg_2", "negative HSCP deta(gen-cscSegment)", 100, -3.0, 3.0);
  TH1* hDphi_cscSeg_2 = new TH1F("hDphi_cscSeg_2", "negative HSCP dphi(gen-cscSegment)", 100, -0.5, 0.5);
  tree->Branch("hDeta_cscSeg_2", &hDeta_cscSeg_2, "TH1F");
  tree->Branch("hDphi_cscSeg_2", &hDphi_cscSeg_2, "TH1F");

  //x-axis VS y-axis
  TH2* hDetaVSDphi_rpcHit_1 = new TH2F("hDetaVSDphi_rpcHit_1","positive HSCP deta VS dphi (gen-rpcHit)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_rpcHit_2 = new TH2F("hDetaVSDphi_rpcHit_2","negative HSCP deta VS dphi (gen-rpcHit)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_rpcHit_1", &hDetaVSDphi_rpcHit_1, "TH2F");
  tree->Branch("hDetaVSDphi_rpcHit_2", &hDetaVSDphi_rpcHit_2, "TH2F");
  hDetaVSDphi_rpcHit_1->SetOption("COLZ");
  hDetaVSDphi_rpcHit_2->SetOption("COLZ");

  TH2* hDetaVSDphi_rpcHit_iRPC_1 = new TH2F("hDetaVSDphi_rpcHit_iRPC_1","positive HSCP deta VS dphi (gen-rpcHit_iRPC)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_rpcHit_iRPC_2 = new TH2F("hDetaVSDphi_rpcHit_iRPC_2","negative HSCP deta VS dphi (gen-rpcHit_iRPC)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_rpcHit_iRPC_1", &hDetaVSDphi_rpcHit_iRPC_1, "TH2F");
  tree->Branch("hDetaVSDphi_rpcHit_iRPC_2", &hDetaVSDphi_rpcHit_iRPC_2, "TH2F");
  hDetaVSDphi_rpcHit_iRPC_1->SetOption("COLZ");
  hDetaVSDphi_rpcHit_iRPC_2->SetOption("COLZ");

  TH2* hDetaVSDphi_gemSeg_1 = new TH2F("hDetaVSDphi_gemSeg_1","positive HSCP deta VS dphi (gen-gemSeg)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_gemSeg_2 = new TH2F("hDetaVSDphi_gemSeg_2","negative HSCP deta VS dphi (gen-gemSeg)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_gemSeg_1", &hDetaVSDphi_gemSeg_1, "TH2F");
  tree->Branch("hDetaVSDphi_gemSeg_2", &hDetaVSDphi_gemSeg_2, "TH2F");
  hDetaVSDphi_gemSeg_1->SetOption("COLZ");
  hDetaVSDphi_gemSeg_2->SetOption("COLZ");

  TH2* hDetaVSDphi_cscSeg_1 = new TH2F("hDetaVSDphi_cscSeg_1","positive HSCP deta VS dphi (gen-cscSeg)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_cscSeg_2 = new TH2F("hDetaVSDphi_cscSeg_2","negative HSCP deta VS dphi (gen-cscSeg)", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_cscSeg_1", &hDetaVSDphi_cscSeg_1, "TH2F");
  tree->Branch("hDetaVSDphi_cscSeg_2", &hDetaVSDphi_cscSeg_2, "TH2F");
  hDetaVSDphi_cscSeg_1->SetOption("COLZ");
  hDetaVSDphi_cscSeg_2->SetOption("COLZ");

  //x-axis VS y-axis, bx=0
  TH2* hDetaVSDphi_bx0_rpcHit_iRPC_1 = new TH2F("hDetaVSDphi_bx0_rpcHit_iRPC_1","positive HSCP deta VS dphi (gen-rpcHit_iRPC), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_bx0_rpcHit_iRPC_2 = new TH2F("hDetaVSDphi_bx0_rpcHit_iRPC_2","negative HSCP deta VS dphi (gen-rpcHit_iRPC), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_bx0_rpcHit_iRPC_1", &hDetaVSDphi_bx0_rpcHit_iRPC_1, "TH2F");
  tree->Branch("hDetaVSDphi_bx0_rpcHit_iRPC_2", &hDetaVSDphi_bx0_rpcHit_iRPC_2, "TH2F");
  hDetaVSDphi_bx0_rpcHit_iRPC_1->SetOption("COLZ");
  hDetaVSDphi_bx0_rpcHit_iRPC_2->SetOption("COLZ");

  TH2* hDetaVSDphi_bx0_gemSeg_1 = new TH2F("hDetaVSDphi_bx0_gemSeg_1","positive HSCP deta VS dphi (gen-gemSeg), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_bx0_gemSeg_2 = new TH2F("hDetaVSDphi_bx0_gemSeg_2","negative HSCP deta VS dphi (gen-gemSeg), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_bx0_gemSeg_1", &hDetaVSDphi_bx0_gemSeg_1, "TH2F");
  tree->Branch("hDetaVSDphi_bx0_gemSeg_2", &hDetaVSDphi_bx0_gemSeg_2, "TH2F");
  hDetaVSDphi_bx0_gemSeg_1->SetOption("COLZ");
  hDetaVSDphi_bx0_gemSeg_2->SetOption("COLZ");

  TH2* hDetaVSDphi_bx0_cscSeg_1 = new TH2F("hDetaVSDphi_bx0_cscSeg_1","positive HSCP deta VS dphi (gen-cscSeg), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  TH2* hDetaVSDphi_bx0_cscSeg_2 = new TH2F("hDetaVSDphi_bx0_cscSeg_2","negative HSCP deta VS dphi (gen-cscSeg), bx0", 100, -0.1, 0.1, 100, -0.1, 0.1);
  tree->Branch("hDetaVSDphi_bx0_cscSeg_1", &hDetaVSDphi_bx0_cscSeg_1, "TH2F");
  tree->Branch("hDetaVSDphi_bx0_cscSeg_2", &hDetaVSDphi_bx0_cscSeg_2, "TH2F");
  hDetaVSDphi_bx0_cscSeg_1->SetOption("COLZ");
  hDetaVSDphi_bx0_cscSeg_2->SetOption("COLZ");

  float out_t0[2];
  tree->Branch("t0_1", &out_t0[0], "t0_1/F");
  tree->Branch("t0_2", &out_t0[1], "t0_2/F");

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if ( (10000*(jentry+1)/nentries) % 100 == 0 ) {
      printf("Processing %lld/%lld, %.0f %% done...\r",  jentry+1, nentries, 100.*jentry/nentries);
    }

    // Initialize variables
    for ( unsigned i=0; i<2; ++i ) {
      out_gens_p4[i].SetXYZT(0,0,0,0);
      out_gens_pdgId[i] = 0;
    }
    out_muons_n = 0;
    for ( unsigned i=0; i<3; ++i ) {
      out_muons_p4[i].SetXYZT(0,0,0,0);
      out_muons_q[i] = 0;
    }
    for ( unsigned i=0; i<2; ++i ) {
      out_rpc_p3[i].SetXYZ(0,0,0);
      out_gem_p3[i].SetXYZ(0,0,0);
      out_csc_p3[i].SetXYZ(0,0,0);
    }
    for ( unsigned i=0; i<2; ++i ) {
      out_fit_quals[i] = 1e9;
      out_fit_betas[i] = 0;
      out_fit_nhits[i] = 0;
      out_t0[i] = 1e9;
    }

    //pdgId, pt
    if ( gen1_pdgId == 0 or gen2_pdgId == 0 ) continue;
    //if ( gen1_pt < 100 or gen2_pt < 100 ) continue;
    // Fill particles
    if ( gen1_pdgId != 0 ) {
      out_gens_p4[0].SetPtEtaPhiM(gen1_pt, gen1_eta, gen1_phi, gen1_m);
      out_gens_pdgId[0] = gen1_pdgId;
    }
    if ( gen2_pdgId != 0 ) {
      out_gens_p4[1].SetPtEtaPhiM(gen2_pt, gen2_eta, gen2_phi, gen2_m);
      out_gens_pdgId[1] = gen2_pdgId;
    }
    std::vector<unsigned> muonIdxs;
    for ( unsigned i=0; i<muon_n; ++i ) {
      if ( muon_pt[i] < 30 or std::abs(muon_eta[i]) > 2.4 ) continue;
      muonIdxs.push_back(i);
    }
    std::sort(muonIdxs.begin(), muonIdxs.end(), [&](unsigned i, unsigned j){return muon_pt[i] > muon_pt[j];});
    for ( unsigned i=0, n=std::min(3ul, muonIdxs.size()); i<n; ++i ) {
      const auto ii = muonIdxs[i];
      out_muons_p4[i].SetPtEtaPhiM(muon_pt[ii], muon_eta[ii], muon_phi[ii], muonMass);
      out_muons_q[i] = muon_q[i];
    }

    // Cluster hits and do the fitting
    std::vector<std::vector<unsigned>> hitClusters;
    if      ( clusterAlgo == ClusterAlgo::GenMatch  ) hitClusters = clusterHitsByGenP4s(out_gens_p4);
    else if ( clusterAlgo == ClusterAlgo::Histogram ) hitClusters = clusterHitsByEtaPhi();
    //std::sort(hitClusters.begin(), hitClusters.end(),
              //[](const std::vector<unsigned>& a, const std::vector<unsigned>& b){return a.size() > b.size();});
    for ( unsigned i=0, n=std::min(2ul, hitClusters.size()); i<n; ++i ) {
      std::vector<double> res;
      if      ( fitAlgo == FitAlgo::BxContrained ) res = fitTrackBxConstrained(hitClusters[i]);
      else if ( fitAlgo == FitAlgo::FitSlope     ) res = fitTrackSlope(hitClusters[i]);

      for ( unsigned j : hitClusters[i] ) {
        //Fill clustered hits
        out_rpc_p3[i].SetXYZ(rpcHit_x[j], rpcHit_y[j], rpcHit_z[j]);
        out_gem_p3[i].SetXYZ(gemSegment_x[j], gemSegment_y[j], gemSegment_z[j]);
        out_csc_p3[i].SetXYZ(cscSegment_x[j], cscSegment_y[j], cscSegment_z[j]);
        //For find cuts
	const TVector3 pos_rpcHit(rpcHit_x[j], rpcHit_y[j], rpcHit_z[j]);
        const TVector3 pos_gem(gemSegment_x[j], gemSegment_y[j], gemSegment_z[j]);
        const TVector3 pos_csc(cscSegment_x[j], cscSegment_y[j], cscSegment_z[j]);
        const double ct_rpc = speedOfLight*rpcHit_time[j];
        const double ibeta_rpc = 1./(1+ct_rpc/pos_rpcHit.Mag());
        const double ct_gem = speedOfLight*gemSegment_time[j];
        const double ibeta_gem = 1./(1+ct_gem/pos_gem.Mag());
        const double ct_csc = speedOfLight*cscSegment_time[j];
        const double ibeta_csc = 1./(1+ct_csc/pos_csc.Mag());
        if( i == 0 ) {
          hDeta_rpcHit_1->Fill(out_gens_p4[0].Eta()-pos_rpcHit.Eta());
          hDphi_rpcHit_1->Fill(out_gens_p4[0].Phi()-pos_rpcHit.Phi());
          hDetaVSDphi_rpcHit_1->Fill( out_gens_p4[0].Eta()-pos_rpcHit.Eta(), out_gens_p4[0].Phi()-pos_rpcHit.Phi() );
          if( std::abs( out_gens_p4[0].Eta() ) > 1.8 and std::abs( out_gens_p4[0].Eta() ) < 2.5 ) {
            hDetaVSDphi_rpcHit_iRPC_1->Fill( out_gens_p4[0].Eta()-pos_rpcHit.Eta(), out_gens_p4[0].Phi()-pos_rpcHit.Phi() );
            if( ibeta_rpc >= -1 and ibeta_rpc <= 2 ) hDetaVSDphi_bx0_rpcHit_iRPC_1->Fill( out_gens_p4[0].Eta()-pos_rpcHit.Eta(), out_gens_p4[0].Phi()-pos_rpcHit.Phi() );
          }

          if( gemSegment_n > 0 ) {
            if ( std::abs(out_gens_p4[0].Vect().DeltaPhi(pos_gem)) < 0.1 )hDeta_gemSeg_1->Fill(out_gens_p4[0].Eta()-pos_gem.Eta());
            hDphi_gemSeg_1->Fill(out_gens_p4[0].Vect().DeltaPhi(pos_gem));
            hDetaVSDphi_gemSeg_1->Fill( out_gens_p4[0].Eta()-pos_gem.Eta(), out_gens_p4[0].Vect().DeltaPhi(pos_gem) );
            if( ibeta_gem >= -1 and ibeta_gem <= 2 ) hDetaVSDphi_bx0_gemSeg_1->Fill( out_gens_p4[0].Eta()-pos_gem.Eta(), out_gens_p4[0].Vect().DeltaPhi(pos_gem) );
          }
          if( cscSegment_n > 0 ) {
            if ( std::abs(out_gens_p4[0].Vect().DeltaPhi(pos_csc)) < 0.1 )hDeta_cscSeg_1->Fill(out_gens_p4[0].Eta()-pos_csc.Eta());
            hDphi_cscSeg_1->Fill(out_gens_p4[0].Vect().DeltaPhi(pos_csc));
            hDetaVSDphi_cscSeg_1->Fill( out_gens_p4[0].Eta()-pos_csc.Eta(), out_gens_p4[0].Vect().DeltaPhi(pos_csc) );
            if( ibeta_csc >= -1 and ibeta_csc <= 2 ) hDetaVSDphi_bx0_cscSeg_1->Fill( out_gens_p4[0].Eta()-pos_csc.Eta(), out_gens_p4[0].Vect().DeltaPhi(pos_csc) );
          }
          
        }
        else if ( i == 1 ) {
          hDeta_rpcHit_2->Fill(out_gens_p4[1].Eta()-pos_rpcHit.Eta());
          hDphi_rpcHit_2->Fill(out_gens_p4[1].Phi()-pos_rpcHit.Phi());
          hDetaVSDphi_rpcHit_2->Fill( out_gens_p4[1].Eta()-pos_rpcHit.Eta(), out_gens_p4[1].Phi()-pos_rpcHit.Phi() );
          if( std::abs( out_gens_p4[1].Eta() ) > 1.8 and std::abs( out_gens_p4[1].Eta() ) < 2.5 ) {
            hDetaVSDphi_rpcHit_iRPC_2->Fill( out_gens_p4[1].Eta()-pos_rpcHit.Eta(), out_gens_p4[1].Phi()-pos_rpcHit.Phi() );
            if( ibeta_rpc >= -1 and ibeta_rpc <= 2 ) hDetaVSDphi_bx0_rpcHit_iRPC_2->Fill( out_gens_p4[1].Eta()-pos_rpcHit.Eta(), out_gens_p4[1].Phi()-pos_rpcHit.Phi() );
          }

          if( gemSegment_n > 0 ) {
            if ( std::abs(out_gens_p4[1].Vect().DeltaPhi(pos_gem)) < 0.1 )hDeta_gemSeg_2->Fill(out_gens_p4[1].Eta()-pos_gem.Eta());
            hDphi_gemSeg_2->Fill(out_gens_p4[1].Vect().DeltaPhi(pos_gem));
            hDetaVSDphi_gemSeg_2->Fill( out_gens_p4[1].Eta()-pos_gem.Eta(), out_gens_p4[1].Vect().DeltaPhi(pos_gem) );
            if( ibeta_gem >= -1 and ibeta_gem <= 2 ) hDetaVSDphi_bx0_gemSeg_2->Fill( out_gens_p4[1].Eta()-pos_gem.Eta(), out_gens_p4[1].Vect().DeltaPhi(pos_gem) );
          }
          if( cscSegment_n > 0 ) {
            if ( std::abs(out_gens_p4[1].Vect().DeltaPhi(pos_csc)) < 0.1 )hDeta_cscSeg_2->Fill(out_gens_p4[1].Eta()-pos_csc.Eta());
            hDphi_cscSeg_2->Fill(out_gens_p4[1].Vect().DeltaPhi(pos_csc));
            hDetaVSDphi_cscSeg_2->Fill( out_gens_p4[1].Eta()-pos_csc.Eta(), out_gens_p4[1].Vect().DeltaPhi(pos_csc) );
            if( ibeta_csc >= -1 and ibeta_csc <= 2 ) hDetaVSDphi_bx0_cscSeg_2->Fill( out_gens_p4[1].Eta()-pos_csc.Eta(), out_gens_p4[1].Vect().DeltaPhi(pos_csc) );
          }
        }
      }

      if ( fitAlgo == FitAlgo::BxContrained ) out_fit_quals[i] = res[0]; //bxErr2, using at BxConstrained algo
      else if ( fitAlgo == FitAlgo::FitSlope ) out_fit_quals[i] = res[2]; //betaErr, using at FitSlope algo 
      out_fit_betas[i] = res[1];
      out_fit_nhits[i] = hitClusters[i].size();
      //out_t0[i] = res[3];

      /////////simDigi///////////////
      const TVector3 pos1(simDigi1_x[i], simDigi1_y[i], simDigi1_z[i]);

      hDeta_1->Fill(out_gens_p4[0].Eta()-pos1.Eta());
      hDphi_1->Fill(out_gens_p4[0].Phi()-pos1.Phi());
        
      const TVector3 pos2(simDigi2_x[i], simDigi2_y[i], simDigi2_z[i]);

      hDeta_2->Fill(out_gens_p4[1].Eta()-pos2.Eta());
      hDphi_2->Fill(out_gens_p4[1].Phi()-pos2.Phi());

    }
        
    if ( out_fit_quals[0] < 1e9 and out_fit_quals[1] < 1e9) ++twoHSCP;
    if ( out_fit_quals[0] < 1e9 or out_fit_quals[1] < 1e9) ++oneHSCP;
    ++nEvent;
    //if ( out_fit_quals[0] >= 1e9 or out_fit_quals[1] >= 1e9 ) continue;

    tree->Fill();
    
    
  }
  cout << "Processing " << nentries << "/" << nentries << "\n"; // Just to print last event
  oneTotEff = (double)oneHSCP/nEvent;
  twoTotEff = (double)twoHSCP/nEvent;
  cout << "Total Efficiency of single HSCP : " << oneHSCP << "/" << nEvent << " = " << oneTotEff << endl;
  cout << "Total Efficiency of pair HSCP : " << twoHSCP << "/" << nEvent << " = " << twoTotEff << endl;
  fout->Write();

  

}
