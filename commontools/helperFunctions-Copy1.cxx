#define helperFunctions_cxx


#include <ROOT/RVec.hxx>
#include "CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"

using VecF_t = const ROOT::RVec<float>&;
using VecD_t = const ROOT::RVec<double>&;
using VecI_t = const ROOT::RVec<int>&;
using VecUI_t = const ROOT::RVec<UInt_t>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;

#include "helperFunctions.h"

ROOT::RVec<double> wgt_vec;
ROOT::RVec<ULong64_t> evn_vec;



bool myfilter(float x) {
   return x > 5;
}

auto sum = [](int a, int b) {
        return a + b;
    };

std::vector<float> getTaggerProb(ROOT::VecOps::RVec<UInt_t> key, ROOT::VecOps::RVec<UInt_t> index, ROOT::VecOps::RVec<float> tagger){
  std::vector< float > result;
  for (UInt_t i = 0; i < key.size(); i++){
    result.push_back(tagger[index[i]]);
  }
  return result;
}


std::vector<int> checkJVT(ROOT::VecOps::RVec<char> NNJVT){
  std::vector< int > result;
  for (UInt_t i = 0; i < NNJVT.size(); i++){
    result.push_back(int(NNJVT[i]));
  }
  return result;
}

Double_t getScaleFactor(UInt_t this_did, VecI_t dids, VecD_t sfs){
    std::cout<<"this dsid = "<<this_did<<std::endl;
    for(unsigned int i=0; i<dids.size(); i++){
        std::cout<<"dids.at("<<i<<") = "<< dids.at(i) << std::endl;
        if(dids.at(i) == this_did){
           return sfs.at(i);
        }
    }
    return -1.0;
}

Double_t getSumOfWeights(ROOT::VecOps::RVec<std::string> name, VecD_t& sow, std::string key){
    for(unsigned int i=0; i<name.size(); i++){
        if(name.at(i) == key){
           return sow.at(i);
        }
    }
    return 1.0;
}

int checkTriggerMatch(VecF_t& pt, VecB_t trigmatch, Float_t th){
  int n_match = 0;
  for(unsigned int i=0; i<pt.size(); i++)
    {
      if(trigmatch[i] && pt[i] > th)n_match += 1;
    }
  return n_match;
}

std::vector<int> getUniqueDSIDs(VecI_t alldsid){
  //using namespace ROOT::VecOps;
  std::vector<int> unique_dsid;
  //auto v_1 = Combinations(alldsid,1);
  std::cout<<"Length of alldsid = "<<alldsid.size()<<std::endl;
  for(unsigned int i=0; i<alldsid.size(); i++){
    int key = alldsid.at(i);
    if(std::count(unique_dsid.begin(), unique_dsid.end(), key))
      continue;
    std::cout<<"Found a new dsid "<<key<<std::endl;
    unique_dsid.push_back(key);  
  }
  return unique_dsid;
}

Double_t getMetRel(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi){
  TLorentzVector l;
  Double_t min_dphi_lep_met = 9999;
  TLorentzVector met;
  met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
  for(unsigned int i=0; i<pt.size(); i++)
    {
      l.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
      Double_t dphi = fabs(l.DeltaPhi(met));
      if(dphi < min_dphi_lep_met){
	min_dphi_lep_met = dphi;
      }
    }
  return (min_dphi_lep_met < M_PI/2.0) ? (met_et)*sin(min_dphi_lep_met) : (met_et);
 
}



std::pair <int,int> num_bl_sg_lep(VecF_t& pt, VecF_t& eta, VecI_t& fllep, VecB_t passOR, VecB_t passLOOSE, VecB_t passMEDIUM, VecB_t passBL, VecF_t z0sinth, VecB_t ISO, VecF_t d0sig, VecF_t passTIGHT){
  int nbl = 0;
  int nsg = 0;
  std::pair <int,int> nlep;
  for(unsigned int i=0; i<fllep.size(); i++)
    {
      if(pt[i] < 9)continue;
      if((fllep[i] == 1 && fabs(eta[i])>2.47) || ((fllep[i] == 2 && fabs(eta[i])>2.6)))continue;
      if(!passOR[i])continue;
      if((fllep[i] == 1 && (!passLOOSE[i] || !passBL[i])) || (fllep[i] == 2 && !passMEDIUM[i]))continue;
      if(fabs(z0sinth[i])>0.5)continue;
      
      nbl += 1;
      
      if((fllep[i] == 1 && !passTIGHT[i]))continue;
      if(!ISO[i])continue;
      if((fllep[i] == 1 && fabs(d0sig[i])>5) || ((fllep[i] == 2 && fabs(d0sig[i])>3)))continue;
      
      nsg += 1;
    }
  nlep = std::make_pair(nbl,nsg);
  return nlep;
}


std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi){
  double diff = 10000000000.0;
  /**
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Wlep1 = -999;
  */
  double Zmass = -1.0;
  double Wmass = -1.0;
  bool foundSFOS = false;
  std::pair <double,double> masses;
  for(unsigned int i=0; i<chlep.size(); i++)
    {
      for(unsigned int j=i+1; j<chlep.size(); j++)
	{
	  //Opposite-Sign
	  if(chlep[i]*chlep[j]<0)
	    {
	      //Same-Flavor
	      if(abs(fllep[i])==abs(fllep[j]))
		{
		  TLorentzVector p1;
		  p1.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
		  TLorentzVector p2;
		  p2.SetPtEtaPhiM(pt[j], eta[j], phi[j], e[j]);
		  double mass = (p1+p2).M();
		  double massdiff = fabs(mass-91187.6);
		  if(massdiff<diff)
		    {
		      diff=massdiff;
		      Zmass=mass;
		      Zlep1 = i;
		      Zlep2 = j;
		      foundSFOS = true;
		    }
		}
	    }
	}

    }
  
  if(foundSFOS){
    TLorentzVector met;
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if((Zlep1==0 && Zlep2==1) || (Zlep1==1 && Zlep2==0) ) Wlep1=2;
    else if((Zlep1==0 && Zlep2==2) || (Zlep1==2 && Zlep2==0) ) Wlep1=1;
    else if((Zlep1==1 && Zlep2==2) || (Zlep1==2 && Zlep2==1) ) Wlep1=0;
    
    TLorentzVector lepW;
    lepW.SetPtEtaPhiM(pt[Wlep1], eta[Wlep1], phi[Wlep1], e[Wlep1]);
    double wlepMetphi = lepW.DeltaPhi(met);
    Wmass = sqrt(2*lepW.Pt()*met.Pt()*(1-cos(wlepMetphi)));
  }
  masses = std::make_pair(Zmass,Wmass);
    
  return masses;
}

Float_t getLumiSF(Int_t randrnum){
    Float_t lumi15 = 3219.56;
    Float_t lumi16 = 32988.1;
    Float_t lumi17 = 44307.4;
    Float_t lumi18 = 58450.1;
    if(randrnum < 320000)return lumi15+lumi16;
    else if(randrnum > 320000 && randrnum < 348000)return lumi17;
    else if(randrnum > 348000)return lumi18;
    else{std::cout<<"ERROR \t RandomRunnumber "<<randrnum<<" has no period attached"<<std::endl;}
    return 1.0;
}

double getSF(VecF_t& sf){
  const auto size = sf.size();
  double scalef = 1.0;
  for (size_t i=0; i < size; ++i) {
    scalef *= sf[i];
  }
  return scalef;
}


bool isOS(const ROOT::VecOps::RVec<int>& chlep) {
  if(chlep[0]*chlep[1] < 0)return kTRUE;
  return kFALSE;
}

// bool isTriggerMatched(const ROOT::VecOps::RVec<int>& isTM) {
//   std::vector<int> tm_vec; 
//   const auto size = isTM.size();
//    for (size_t i=0; i < size; ++i) {
//      if(isTM[i])tm_vec.push_back(0);
//      else tm_vec.push_back(1);
//    }
   
//   return kFALSE;
// }

int flavourComp(VecI_t& fllep){//, int nlep) {
  //const auto size = fllep.size();
  int nlep = 2;
  if(nlep==3){
    //std::cout<<"ERROR \t Vector must be at least 3 long!"<<std::endl;
    if(fllep[0] == 1 && fllep[1] == 1 && fllep[2] == 1)return 0;
    if(fllep[0] == 1 && fllep[1] == 1 && fllep[2] == 2)return 1;
    if(fllep[0] == 1 && fllep[1] == 2 && fllep[2] == 2)return 2;
    if(fllep[0] == 2 && fllep[1] == 2 && fllep[2] == 2)return 3;
    if(fllep[0] == 2 && fllep[1] == 2 && fllep[2] == 1)return 4;
    if(fllep[0] == 2 && fllep[1] == 1 && fllep[2] == 1)return 5;
  }else if(nlep==2){
    if(fllep[0] == 1 && fllep[1] == 1)return 6;                                                                                                                                                                                                                    
    if(fllep[0] == 2 && fllep[1] == 2)return 7;                                                                                                                                                                                                                    
    if((fllep[0] == 1 && fllep[1] == 2) || (fllep[0] == 2 && fllep[1] == 1))return 8;                                                                                                                                                                                               
  }
  return -1;
}

bool deltaRlepjet(float lpt, float leta, float lphi, float le, VecF_t& jpt, VecF_t& jeta, VecF_t& jphi, VecF_t& je){
  const auto njet = int(jpt.size());
  TLorentzVector p2;
  TLorentzVector p1;
  double deltaR;
  double mindr = 9999;
  //for (size_t i=0; i < nlep; ++i) {
  p1.SetPtEtaPhiM(lpt, leta, lphi, le);
  for (int j=0; j < njet; ++j) {
    p2.SetPtEtaPhiM(jpt[j], jeta[j], jphi[j], je[j]);
    deltaR = p1.DeltaR(p2);
    if(deltaR < mindr){
      mindr = deltaR;
    }
    //}
  }
  //  std::cout<<"mindr = "<<mindr<<std::endl;
  return mindr;  
}

bool isSF(VecI_t& fllep) {
    if(fllep[0] == fllep[1])return kTRUE;
    return kFALSE;
}

bool isEE(VecI_t& fllep) {
    if(fllep[0] == fllep[1] && fllep[0] == 1)return kTRUE;
    return kFALSE;
}

bool isMM(VecI_t& fllep) {
    if(fllep[0] == fllep[1] && fllep[0] == 2)return kTRUE;
    return kFALSE;
}

float ComputeInvariantMass2L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass, int i, int j) {
  TLorentzVector p1;
  TLorentzVector p2;
  
  p1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
  p2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
  return (p1 + p2).M();
}

float ComputeInvariantMass4L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass) {
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p3;
  TLorentzVector p4;
  
  p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
  p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
  p3.SetPtEtaPhiM(pt[2], eta[2], phi[2], mass[2]);
  p4.SetPtEtaPhiM(pt[3], eta[3], phi[3], mass[3]);
  return (p1 + p2+p3+p4).M();
}

ROOT::RVec<float> getVector(VecF_t& inp1, Float_t m1){
  ROOT::RVec<float> ret_vec;
  const auto ninp1 = int(inp1.size());
  for (int j=0; j < ninp1; ++j) {
    if(m1)ret_vec.push_back(fabs(inp1.at(j))*m1);
    else ret_vec.push_back(inp1.at(j));
  }
  return ret_vec;
}
/*Hva gjor denne*/
ROOT::RVec<float> getVector(VecF_t& inp1, VecF_t& inp2, Float_t m1 , Float_t m2){
  ROOT::RVec<float> ret_vec;
  const auto ninp1 = int(inp1.size());
  for (int j=0; j < ninp1; ++j) {
    if(m1)ret_vec.push_back(fabs(inp1.at(j))*m1);
    else ret_vec.push_back(inp1.at(j));
  }
  const auto ninp2 = int(inp2.size());
  for (int j=0; j < ninp2; ++j) {
    if(m2)ret_vec.push_back(fabs(inp2.at(j))*m2);
    else ret_vec.push_back(inp2.at(j));
  }
  return ret_vec;
}


float calcMT2(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, Float_t met_et, Float_t met_px, Float_t met_py, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("calcMT2::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }
  
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector met;
  p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], m[idx1]);
  p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], m[idx2]);
  met.SetPtEtaPhiM(met_et, 0.0, TMath::ATan(met_py/met_px), 0.0);
  return ComputeMT2(p1,p2,met,0.,0.).Compute();
}

float ptllboost(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

    TLorentzVector p1;
    TLorentzVector p2;
    TLorentzVector met;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    return (met+p1+p2).Pt();
}

float costhetastar(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return TMath::ATan(fabs(p1.Eta()-p2.Eta())/2.);
}

float deltaPhi_ll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e) {

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
    p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);
    return p1.DeltaPhi(p2);
}

float deltaPhi_metl(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

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

float deltaPhi_metll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi) {

  if(pt.size() < 2){
    return -999;
  }
  
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector dil;
  TLorentzVector met;
  p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], e[0]);
  p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], e[1]);

  dil = (TLorentzVector)(p1+p2);
  met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
  return dil.DeltaPhi(met);

}

bool checkPt(VecF_t& pt, float cut1, float cut2){
    if((pt[0] > cut1 && pt[1] > cut2) || (pt[1] > cut1 && pt[0] > cut2))return kTRUE;
    return kFALSE;
}


double getWeight(ULong64_t evnum){
    /**
    dfevnum = wgtdf.Filter(ROOT::Form("EventNumber == %llu",evnum));
    std::cout<<"Found "<<dfevnum.Count()<<" elements"<<std::endl;
    */
    
    // std::cout<<"Size of vector = "<<evn_vec.size()<<std::endl;
    for(int i = 0; i<evn_vec.size(); i++){
        // std::cout<<"Checking "<<evn_vec.at(i) <<" against "<<evnum<<std::endl;
        if(evn_vec.at(i) == evnum){
            return wgt_vec.at(i);
        }        
    }  
    
    return -999;
}

// void writeToFile(ROOT::RDF::RNode df, std::string treename, std::string outfilename, std::vector<std::string> &&good_cols){
    
//     df.Snapshot(treename,outfilename,good_cols);
    
//     return;
// }



