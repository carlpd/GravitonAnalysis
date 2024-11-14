#define Cfunctions_cxx
#include <stdio.h>
#include <ROOT/RVec.hxx>
#include "CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"

using Vec2_t = const ROOT::VecOps::RVec<float>;
using Vec_t = const ROOT::VecOps::RVec<int>;

using Vectlv_t = const ROOT::VecOps::RVec<TLorentzVector>;
#include "Cfunctions.h"

void help(){

  printf("Library of handy functions to be used with RDataFrame\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("isOS(const ROOT::VecOps::RVec<int>& chlep)\n");
  printf("\t Checks if pair of leptons has opposite sign. Returns bool\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("isSF(Vec_t& fllep)\n");
  printf("\t Checks if pair of leptons has same flavour (i.e. electron, muon, tau etc.). Returns bool\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e)\n");
  printf("\t Computes invariant mass of leptons. Input can be any size, function will compute the total invariant mass of all objects.  Returns float\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi)\n");
  printf("\t Computes the stransverse mass [Ref.: https://gitlab.cern.ch/atlas-phys-susy-wg/CalcGenericMT2].  Returns float\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e)\n");
  printf("\t Computes the cos(theta)* of two leptons.  Returns float\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e)\n");
  printf("\t Computes the difference in phi between two leptons.  Returns float\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e)\n");
  printf("\t Computes the difference in phi between missing transverse energy vector and the highest pT lepton.  Returns float\n");
  printf("%.*s", 20, "=================");
  printf("\n");
  printf("checkPt(Vec2_t& pt, float cut1, float cut2)\n");
  printf("\t Checks if leptons passes some certain pT threshold set by the inputs.  Returns bool\n");
}

void readMCCat(std::string infile){
  std::ifstream myfile (infile);
  TString tstr;
  TString descr;
  TString dsid;
  TString cate;
  std::string line;
  

  
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
        {
          //std::cout<<"line = "<<line<<std::endl;
          tstr = line;
          TObjArray *tx = tstr.Tokenize(" ");
          if(tx->GetEntries() >= 3){
            descr = ((((TObjString *)tx->At(0))->String()));
            dsid = ((((TObjString *)tx->At(1))->String()));
            cate = ((((TObjString *)tx->At(2))->String()));

            pairs.push_back(std::make_pair(descr.Data(),cate.Data()));

            // std::cout<<"dsid = "<<dsid<<std::endl;
            // std::cout<<"descr = "<<descr<<std::endl;
            // std::cout<<"cate = "<<cate<<std::endl;
            
            MCcat[dsid.Atoi()] = pairs.back();
          }
        
        }
      myfile.close();
    }
  //std::cout<<MCcat[305580].first<<"   " << MCcat[305580].second<<std::endl;
}

std::pair<const char *, const char *> getMCMetaData(int dsid){

  if ( MCcat.find(dsid) == MCcat.end() ) {
    //std::cout<<"Reading"<<std::endl;
    readMCCat("/storage/shared/software/Input/Signal_samples_13TeV.txt");
    readMCCat("/storage/shared/software/Input/Background_samples_13TeV.txt");
  }// else{
  //   
  // }
  //std::cout<<"For "<<dsid<<" found "<< MCcat[dsid].first << " and " << MCcat[dsid].second << std::endl;
  if ( MCcat.find(dsid) == MCcat.end() ) {
    std::cout<<"Could not find infor for dsid "<<dsid<<std::endl;
  }
  
  
  return MCcat[dsid];

}

std::pair<const char *, const char *> getDataMetaData(int dsid){

  TString d = "data";
  std::pair <const char *,const char *> p = std::make_pair(d.Data(),d.Data());
  return p;

}

bool myfilter(float x) {
   return x > 5;
}


double compareEvents(Int_t rnum_orig, const ROOT::VecOps::RVec<ULong64_t>& rnum_wgt, const ROOT::VecOps::RVec<Double_t>& value) {

  //const auto n_orig = rnum_orig.size();
  const auto n_wgt  = rnum_wgt.size();
  
  //for (size_t i=0; i < n_orig; ++i) {
  for (size_t j=0; j < n_wgt; ++j) {
    if(rnum_orig == rnum_wgt[j]){
      return value[j];
      
    }
  }  
  //}
  
  return -999;
}

int isOS(const ROOT::VecOps::RVec<int>& chlep) {
  if(chlep[0]*chlep[1] < 0)return 1;
  return 0;
}


int isSF(Vec_t& fllep) {
    if(fllep[0] == fllep[1])return 1;
    return 0;
}

//// Taken from https://root-forum.cern.ch/t/drop-columns-from-rdataframe/30910/2
std::vector<std::string> DropColumns(std::vector<std::string> &&good_cols)
{
   // your blacklist
   static const std::vector<std::string> blacklist = {"useless", "columns"};
   // a lambda that checks if `s` is in the blacklist
   auto is_blacklisted = [&blacklist](const std::string &s)  { return std::find(blacklist.begin(), blacklist.end(), s) != blacklist.end(); };

   // removing elements from std::vectors is not pretty, see https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
   good_cols.erase(std::remove_if(good_cols.begin(), good_cols.end(), is_blacklisted), good_cols.end());
   
   return good_cols;
}

/// A wrapper for ROOT's InvariantMass function that takes std::vector instead of RVecs
template <typename T>
T InvariantMassStdVector(std::vector<T>& pt, std::vector<T>& eta, std::vector<T>& phi, std::vector<T>& energy)
{
   assert(pt.size() == 2 && eta.size() == 2 && phi.size() == 2 && energy.size() == 2);

   // We adopt the memory here, no copy
   ROOT::RVec<float> rvPt(pt);
   ROOT::RVec<float> rvEta(eta);
   ROOT::RVec<float> rvPhi(phi);
   ROOT::RVec<float> rvMass(energy);

   return ComputeInvariantMass(rvPt, rvEta, rvPhi, rvMass);
}

float ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {
  TLorentzVector total;
  const auto size = pt.size();
  std::vector<TLorentzVector> tlv(size);
  for (size_t i=0; i < size; ++i) {
    tlv[i].SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
    if(i == 0)total = tlv[i];
    else total += tlv.back();		
  }
  //p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
  //p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
  return (total).M();
}

float getVar(Vec2_t& var, int idx){
    const auto size = var.size();
    if(idx > size){
    std::cout<<"Can not ask for idx "<<idx<<" when there are only "<< size <<" object(s)"<<std::endl;
    return -1;       
  }
    return var[idx];  
}

int countBJets(Vec2_t& pt, Vec2_t& tagger, float ptcut, int WP){
    btagWP[0]=-9999;
    btagWP[85]=0.1758;
    btagWP[77]=0.6459;
    btagWP[70]=0.8244;
    btagWP[60]=0.9349;
    int nbjet = 0;
    const auto size = pt.size();
    for (size_t i=0; i < size; ++i) {
        if(pt[i] > ptcut and tagger[i] > btagWP[WP])nbjet+=1;
    }
    return nbjet;
}

int countJets(Vec2_t& pt, float ptcut){
    int njet = 0;
    const auto size = pt.size();
    for (size_t i=0; i < size; ++i) {
        if(pt[i] > ptcut)njet+=1;
    }
    return njet;
}

int getTypeTimesCharge(Vec_t& chlep, Vec_t& fllep, int idx){
    const auto sizech = chlep.size();
    const auto sizefl = fllep.size();
    
    if(idx > sizech or idx > sizefl){
    std::cout<<"Can not ask for idx "<<idx<<" when there are only "<< sizech <<" object(s)"<<std::endl;
    return -1;       
  }
    if(sizech != sizefl){
    std::cout<<"Different size of fl "<<sizefl<<" and charge "<< sizech <<std::endl;
    return -1;       
  }
    return chlep[idx]*fllep[idx];  
}


float calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return ComputeMT2(p1,p2,met,0.,0.).Compute();
}

float ptllboost(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return (met+p1+p2).Pt();
}

float costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return TMath::ATan(fabs(p1.Eta()-p2.Eta())/2.);
}

float deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return p1.DeltaPhi(p2);
}

float deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if(p1.Pt() > p2.Pt()){
        return p1.DeltaPhi(met);
    }else{
        return p2.DeltaPhi(met);
    }
}

bool checkPt(Vec2_t& pt, float cut1, float cut2){
    if((pt[0] > cut1 && pt[1] > cut2) || (pt[1] > cut1 && pt[0] > cut2))return kTRUE;
    return kFALSE;
}

