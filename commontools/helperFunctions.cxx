#define helperFunctions_cxx


#include <ROOT/RVec.hxx>
#include "CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"

using VecF_t = const ROOT::RVec<float>&;
using VecD_t = const ROOT::RVec<double>&;
using VecI_t = const ROOT::RVec<int>&;
using VecUI_t = const ROOT::RVec<UInt_t>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;
using VecF_t = const ROOT::RVec<float>&;

#include "helperFunctions.h"

ROOT::RVec<double> wgt_vec;
ROOT::RVec<ULong64_t> evn_vec;


// To read data correctly
void setGRL(int vb)                                                                                                                                                                                                       
{

  std::ifstream input("/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data22_13p6TeV/20230207/data22_13p6TeV.periodAllYear_DetStatus-v109-pro28-04_MERGED_PHYS_StandardGRL_All_Good_25ns_ignore_TRIGLAR.xml");
  bool startLB = false;
  bool firstRun = true;
  std::vector<UInt_t> thislbvec;
  UInt_t thisrun;
  for( std::string line; getline( input, line ); )
    {
      TString l = line;
      // <Run>430536</Run>
      if(l.Contains("<Run>") && l.Contains("</Run>")){
        if(!firstRun){
          //GRL_lb[thisrun] = thislbvec;
          GRL_lb.insert( std::make_pair(thisrun, thislbvec) );
          thislbvec.clear();
          //firstRun = false;
        }
        int st = l.First('>');
        int en = 6;//l.Last("Run");
        if(vb)printf("Index start = %i, stop = %i\n",st,en);
        TString runnum( l(st+1,en) );
        thisrun = (UInt_t)runnum.Atoi();
        
        firstRun = false;
      }else if(l.Contains("<LBRange")){
        TObjArray *tx = l.Tokenize("\"");
        UInt_t start = (UInt_t)(((TObjString *)tx->At(1))->String()).Atoi();
        UInt_t end   = (UInt_t)(((TObjString *)tx->At(3))->String()).Atoi();
        if(vb)printf("Run %i has LB start = %i and end = %i\n",thisrun,start,end);
        for(UInt_t i = start; i<=end; i++){
          //printf("Adding %i to run %i\n",i,thisrun);
          thislbvec.push_back(i);
        }
        
      }
      
    }
  // Need to get the last run too!
  GRL_lb.insert( std::make_pair(thisrun, thislbvec) );
  thislbvec.clear();
}


Int_t checkLB(UInt_t rnum, UInt_t lb, int vb){
  //time_t start = time(NULL);
  //GRLlist.Summary(kTRUE);
  //setGRL();
  //vb = 1;
  if(vb){
    printf("Checking LB %i\n",lb);
    for (auto x = GRL_lb.begin(); x != GRL_lb.end(); x++){
      if(x->first != rnum)continue;
      printf("Run %i has following LBs:\n",x->first);
      for(UInt_t j= 0; j< (x->second).size(); j++){
        printf("%i  ",(x->second).at(j));
      }
    }
  }
  int isGRL = 0;//GRLlist.HasRunLumiBlock(rnum,lb);
  if(GRL_lb.find(rnum) != GRL_lb.end()) {
    if ( std::find(GRL_lb[rnum].begin(), GRL_lb[rnum].end(), lb) != GRL_lb[rnum].end() ){
      isGRL = 1;
    }
  }
  //time_t end = time(NULL);
  //elapsed_time["checkLB"] += difftime(end, start);
  n_calls["checkLB"] += 1;
  return isGRL;
}

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

ROOT::VecOps::RVec<bool>  checkJVT(ROOT::VecOps::RVec<char> NNJVT){
  time_t start = time(NULL);
  ROOT::VecOps::RVec<bool> result;
  for (UInt_t i = 0; i < NNJVT.size(); i++){
    result.push_back(bool(NNJVT[i]));
  }
  time_t end = time(NULL);
  elapsed_time["getTaggerProb"] += difftime(end, start);
  n_calls["getTaggerProb"] += 1;
  return result;
}




void readMetaData(std::string infile){
  lumi = 1000.;
  metadata_rdf = ROOT::RDF::FromCSV(infile);
  // std::cout<<"Number of entris in metadata ttree is "<<metadata_rdf.Count().GetValue()<<std::endl;
}

Double_t getScalefactor(Long64_t dsid){
  
  time_t start = time(NULL);
  //TString dsid_str  = Form("%lld",dsid);
  if(scalefactors.find(dsid) != scalefactors.end()){
    time_t end = time(NULL);
    elapsed_time["getScaleFactor"] += difftime(end, start);
    n_calls["getScaleFactor"] += 1;
    return scalefactors[dsid];
  }
  //std::cout<<"getScalefactor::before"<<std::endl;
  auto rdf = metadata_rdf.Filter(Form("dataset_number == %lld",dsid));
  Double_t xsec = rdf.Take<Double_t>("crossSection")->at(0);
  Double_t kfac = rdf.Take<Double_t>("kFactor")->at(0);
  Double_t filtereff = rdf.Take<Double_t>("genFiltEff")->at(0);
  Double_t sow = rdf.Take<Double_t>("sumofweights")->at(0);
  
  if(sow >0){
    scalefactors.insert( std::make_pair(dsid, (xsec*kfac*filtereff*lumi)/sow) );
    time_t end = time(NULL);
    elapsed_time["getScaleFactor"] += difftime(end, start);
    n_calls["getScaleFactor"] += 1;
    return scalefactors[dsid];   
  }else{
    scalefactors.insert( std::make_pair(dsid, -1.0 ));
    time_t end = time(NULL);
    elapsed_time["getScaleFactor"] += difftime(end, start);
    n_calls["getScaleFactor"] += 1;
    return -1;
  }
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

ROOT::RVec<float> getChargeList(VecF_t ChEl, VecF_t ChMu){
  ROOT::RVec<float> ChLep;
  for(unsigned int i=0; i<ChEl.size(); i++){
    ChLep.push_back(ChEl[i]);
  }
  for(unsigned int j=0; j<ChMu.size(); j++){
    ChLep.push_back(ChMu[j]);
  }
  return(ChLep);
}

ROOT::RVec<int> getFlavourList(VecF_t ChEl, VecF_t ChMu){
  ROOT::RVec<int> FlLep;
  for(unsigned int i=0; i<ChEl.size(); i++){
    if(ChEl[i]==(-1) || ChEl[i]==1){
      FlLep.push_back(1);
    }
  }
  for(unsigned int j=0; j<ChMu.size(); j++){
    if(ChMu[j]==(-1) || ChMu[j]==1){
      FlLep.push_back(2);
    }
  }
  return(FlLep);
}

int getActualIDFromFlavourList(int FlavourListObject){
  int flavour;
  if(FlavourListObject==1){
    flavour=getElecPdgID();
  }
  if(FlavourListObject==2){
    flavour=getMuonPdgID();
  }
  return flavour;
}

ROOT::RVec<float> getPtList(VecF_t PtEl, VecF_t PtMu){
  ROOT::RVec<float> PtLep;
  for(unsigned int i=0;i<PtEl.size(); i++){
    PtLep.push_back(PtEl[i]);
  }
  for(unsigned int j=0;j<PtMu.size(); j++){
    PtLep.push_back(PtMu[j]);
  }
  return(PtLep);

}
ROOT::RVec<float> sortedPtList(ROOT::RVec<float> PtLep){
  ROOT::RVec<float> SortedPt =Sort(PtLep);
  SortedPt=Reverse(SortedPt);
  return SortedPt;
}
ROOT::RVec<float> getEtaList(VecF_t EtaEl, VecF_t EtaMu){
  ROOT::RVec<float> EtaLep;
  for(unsigned int i=0;i<EtaEl.size(); i++){
    EtaLep.push_back(EtaEl[i]);
  }
  for(unsigned int j=0;j<EtaMu.size(); j++){
    EtaLep.push_back(EtaMu[j]);
  }
  return(EtaLep);
}

ROOT::RVec<float> getPhiList(VecF_t PhiEl, VecF_t PhiMu){
  ROOT::RVec<float> PhiLep;
  for(unsigned int i=0;i<PhiEl.size(); i++){
    PhiLep.push_back(PhiEl[i]);
  }
  for(unsigned int j=0;j<PhiMu.size(); j++){
    PhiLep.push_back(PhiMu[j]);
  }
  return(PhiLep);
}

ROOT::RVec<float> getMList(VecF_t ChEl, VecF_t ChMu, float EM, float MM){
  ROOT::RVec<float> MLep;
  for(unsigned int i=0; i<ChEl.size(); i++){
    if(ChEl[i]==(-1) || ChEl[i]==1){
      MLep.push_back(EM);
    }
  }
  for(unsigned int j=0; j<ChMu.size(); j++){
    if(ChMu[j]==(-1) || ChMu[j]==1){
      MLep.push_back(MM);
    }
  }
  return(MLep);
}


//std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi){
std::pair <double,double> getLeptonsFromZ(VecF_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e){
  double diff = 10000000000.0;
  double diff2 = 10000000000.0;
  
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Zlep3 = -999;
  int Zlep4 = -999;
  double Zmass = -1.0;
  double Zmass2 = -1.0;
  bool foundSFOS = false;
  bool foundSFOS2 =false;
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
      /*
      if(mass<5000){
        continue; // Filtering out the mll<5GeV
      }
      */
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
    for(unsigned int i=0; i<chlep.size(); i++)
      {if(i==Zlep1){
        continue;
      }
      if(i==Zlep2){
        continue;
      }
      for(unsigned int j=i+1; j<chlep.size(); j++){
        if(j==Zlep1){
          continue;
        }
        if(j==Zlep2){
          continue;
        }
        if(chlep[i]*chlep[j]<0){
          if(abs(fllep[i])==abs(fllep[j])){
            TLorentzVector p3;
            p3.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
            TLorentzVector p4;
            p4.SetPtEtaPhiM(pt[j], eta[j], phi[j], e[j]);
            double mass2=(p3+p4).M();
            /*
            if(mass2<5000){
              continue; //To filter out the mll<5GeV
            }
            */
            double massdiff2 = fabs(mass2-91187.6);
            if(massdiff2<diff2){
              diff2=massdiff2;
              Zmass2=mass2;
              Zlep3=i;
              Zlep4=j;
              foundSFOS2=true;
            }
      
          }
        }


      }}
    /*
    TLorentzVector met;
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if((Zlep1==0 && Zlep2==1) || (Zlep1==1 && Zlep2==0) ) Wlep1=2;
    else if((Zlep1==0 && Zlep2==2) || (Zlep1==2 && Zlep2==0) ) Wlep1=1;
    else if((Zlep1==1 && Zlep2==2) || (Zlep1==2 && Zlep2==1) ) Wlep1=0;
    
    TLorentzVector lepW;
    lepW.SetPtEtaPhiM(pt[Wlep1], eta[Wlep1], phi[Wlep1], e[Wlep1]);
    double wlepMetphi = lepW.DeltaPhi(met);
    Wmass = sqrt(2*lepW.Pt()*met.Pt()*(1-cos(wlepMetphi)));
    */
  }
  masses = std::make_pair(Zmass,Zmass2);
    
  return masses;
}

bool CloseToZ(double mll, double diffOK){
  double m_Z=91187.6;
  double mdif;
  mdif=abs(mll-m_Z);
  bool ac;
  if(mdif<=diffOK){
    ac=true;
  }
  if(mdif>diffOK){
    ac=false;
  }
  return ac;
}

 ROOT::RVec<int> getLeptonsPairsFromZ(VecF_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e){
  double diff = 10000000000.0;
  double diff2 = 10000000000.0;
  
  int Zlep1 = -99;
  int Zlep2 = -99;
  int Zlep3 = -999;
  int Zlep4 = -999;
  double Zmass = -1.0;
  double Zmass2 = -1.0;
  bool foundSFOS = false;
  bool foundSFOS2 =false;
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
      /*
      if(mass<5000){
        continue; // Filtering out the mll<5GeV
      }
      */
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
    for(unsigned int i=0; i<chlep.size(); i++)
      {if(i==Zlep1){
        continue;
      }
      if(i==Zlep2){
        continue;
      }
      for(unsigned int j=i+1; j<chlep.size(); j++){
        if(j==Zlep1){
          continue;
        }
        if(j==Zlep2){
          continue;
        }
        if(chlep[i]*chlep[j]<0){
          if(abs(fllep[i])==abs(fllep[j])){
            TLorentzVector p3;
            p3.SetPtEtaPhiM(pt[i], eta[i], phi[i], e[i]);
            TLorentzVector p4;
            p4.SetPtEtaPhiM(pt[j], eta[j], phi[j], e[j]);
            double mass2=(p3+p4).M();
            /*
            if(mass2<5000){
              continue; //To filter out the mll<5GeV
            }
            */
            double massdiff2 = fabs(mass2-91187.6);
            if(massdiff2<diff2){
              diff2=massdiff2;
              Zmass2=mass2;
              Zlep3=i;
              Zlep4=j;
              foundSFOS2=true;
            }
      
          }
        }


      }}
    /*
    TLorentzVector met;
    met.SetPtEtaPhiM(met_et, 0.0, met_phi, 0.0);
    
    if((Zlep1==0 && Zlep2==1) || (Zlep1==1 && Zlep2==0) ) Wlep1=2;
    else if((Zlep1==0 && Zlep2==2) || (Zlep1==2 && Zlep2==0) ) Wlep1=1;
    else if((Zlep1==1 && Zlep2==2) || (Zlep1==2 && Zlep2==1) ) Wlep1=0;
    
    TLorentzVector lepW;
    lepW.SetPtEtaPhiM(pt[Wlep1], eta[Wlep1], phi[Wlep1], e[Wlep1]);
    double wlepMetphi = lepW.DeltaPhi(met);
    Wmass = sqrt(2*lepW.Pt()*met.Pt()*(1-cos(wlepMetphi)));
    */
  }
  ROOT::RVec<int> pairlist;
  pairlist.push_back(Zlep1);
  pairlist.push_back(Zlep2);
  pairlist.push_back(Zlep3);
  pairlist.push_back(Zlep4);
  //std::cout<<"pairlist"<<pairlist.at(0)<<pairlist.at(3)<<std::endl;
  return pairlist;
}

bool isFirstPairElectrons(VecI_t pairlist, VecI_t fllist){
  bool isFirstElectron;
  int firstlep=fllist[pairlist[0]];
  int secondlep=fllist[pairlist[1]];
  if(firstlep != secondlep){
    printf("ERROR The first leptonpair does not have the same flavour");
  }
  if(firstlep==1){
    isFirstElectron=true;
  }
  return isFirstElectron;
}
float getInvariantMass_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2){
  TLorentzVector p1;
	p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
  TLorentzVector p2;
  p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
  double mass = (p1+p2).M();
  return mass;
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

float ComputeInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass) {
  TLorentzVector p1;
  TLorentzVector p2;
  
  p1.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
  p2.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
  return (p1 + p2).M();
}

float ComputeInvariantMass2L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("calcMT2::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }
  TLorentzVector p1;
  TLorentzVector p2;
  
  p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], mass[idx1]);
  p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], mass[idx2]);
  return (p1 + p2).M();
}
float ComputeInvariantMass4L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass) {
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p3;
  TLorentzVector p4;

  const auto size = int(pt.size());
  if (size!=4){
    printf("ComputeInvariantMass4L::ERROR \t Particles are %i and not 4\n",size);
    return -1;
  }
  
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
/*
float deltaPhi_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("deltaR_ll::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
    p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
    return p1.DeltaPhi(p2);
}
*/
float deltaPhi_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("deltaR_ll::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }

    TLorentzVector p1;
    TLorentzVector p2;
    p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
    p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
    return p1.DeltaPhi(p2);
}
/*
float findE(VecF_t& pt, int idx){
  p1.
}
*/
float deltaR_ll_4(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, int idx1, int idx2) {

  const auto size = int(pt.size());
  if(idx1 > size || idx2 > size){
    printf("deltaR_ll::ERROR \t Indices %i and %i are higher than size of vector %i\n",idx1,idx2,size);
    return -1;
  }

  TLorentzVector p1;
  TLorentzVector p2;
  p1.SetPtEtaPhiM(pt[idx1], eta[idx1], phi[idx1], e[idx1]);
  p2.SetPtEtaPhiM(pt[idx2], eta[idx2], phi[idx2], e[idx2]);
  //std::cout<<"deltaR"<<p1.DeltaR(p2)<<std::endl;
  return p1.DeltaR(p2);
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
    for(unsigned int i = 0; i<evn_vec.size(); i++){
        // std::cout<<"Checking "<<evn_vec.at(i) <<" against "<<evnum<<std::endl;
        if(evn_vec.at(i) == evnum){
            return wgt_vec.at(i);
        }        
    }  
    
    return -999;
}

ROOT::VecOps::RVec<bool> checkIsolation(VecF_t& iso1, VecF_t& iso2, VecF_t& iso3, VecF_t& pt, float cutval){
  ROOT::VecOps::RVec<bool> result;
  for (UInt_t i = 0; i < iso1.size(); i++){
    //if(((TMath::Max(iso1.at(i),iso2.at(i)) + 0.4*iso3.at(i))/pt.at(i)) < cutval){
    result.emplace_back(((TMath::Max(iso1.at(i),iso2.at(i)) + 0.4*iso3.at(i))/pt.at(i)) < cutval);
    //}else{
    //result.push_back(false);
    //} 
  }
  return result;
}

bool isCleanEvent(UInt_t larFlag, UInt_t tileFlag, UInt_t sctFlag, UInt_t coreFlag){


  bool isCleanEvent = !( ((( larFlag >> 28 ) & 0xF) == 2) || ((( tileFlag >> 28 ) & 0xF) == 2) || ((( sctFlag >> 28 ) & 0xF) == 2) || ( ( coreFlag >> 18 ) & 0x1 ));
  
  return isCleanEvent;
  
}

float ComputeBornMass(VecF_t& px, VecF_t& py, VecF_t& pz, VecF_t& e, Long64_t dsid) {
  TLorentzVector p1;
  TLorentzVector p2;

  if(!((dsid >= 700320 && dsid <= 700328) || (dsid >= 700615 && dsid <= 700623)))return 0;
  
  const auto ninpt = int(px.size());
  if(ninpt > 2){
    std::cout<<"ERROR \t Number of born leptons is "<<ninpt<<std::endl;
    for (int j=0; j < ninpt; ++j) {
      std::cout<<"(px,py,pz,e) = ("<<px[j]<<","<<py[j]<<","<<pz[j]<<","<<e[j]<<")"<<std::endl;
    }
    
  }
  
  p1.SetPxPyPzE(px[0], py[0], pz[0], e[0]);
  p2.SetPxPyPzE(px[1], py[1], pz[1], e[1]);
 
  return (p1 + p2).M();
}

ROOT::VecOps::RVec<bool> JetLepOR(VecF_t& lep_pt, VecF_t& lep_eta, VecF_t& lep_phi, VecF_t& jet_pt, VecF_t& jet_eta, VecF_t& jet_phi, VecF_t& jet_m, float cutval, int isEl){
  TLorentzVector lep;
  TLorentzVector jet;
  ROOT::VecOps::RVec<bool> result;
  Float_t mass = 0.511 ? isEl : 105.66;
  bool goodjet;
  for(unsigned int j=0; j<jet_pt.size(); j++){
    goodjet = true;
    jet.SetPtEtaPhiM(jet_pt[j], jet_eta[j], jet_phi[j], jet_m[j]);
    for(unsigned int i=0; i<lep_pt.size(); i++){
      lep.SetPtEtaPhiM(lep_pt[i], lep_eta[i], lep_phi[i], mass);
      if(jet.DeltaR(lep) < cutval){
	goodjet = false;
	break;
      }
    }
     result.emplace_back(goodjet);
  }
  return result;
}


double mthFunction(double m4l){
  double xa=140.0;
  double xb=190.0;
  double ya=12;
  double yb=50;

  double mth;
  if(m4l<=xa){
    mth=ya;
  }
  if(m4l>=xb){
    mth=yb;
  }
  if(m4l>xa && m4l<xb){
    double c=(yb-ya)/(xb-xa);
    double d=ya-c*xa;
    mth=c*m4l+d;
  }
  return mth;
}

bool MassesInThreshold(double mll1, double mll2, double m4l){
  bool m1=false;
  if(mll1>40 && mll1<106){
    m1=true;
  }
  bool m2=false;
  double mth=mthFunction(m4l);
  if(mll2>mth && mll2<115){
    m2=true;
  }
  if(m1 && m2){
    return true;
  }
  else{
    return false;
  }
}

std::pair <double,double> unpairedSFOSmasses(ROOT::RVec<int> lepPairs, VecF_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e){
  int a=lepPairs[0];
  int b=lepPairs[1];
  int c=lepPairs[2];
  int d=lepPairs[3];
  bool acpair=false;
  bool adpair=false;
  bool bcpair=false;
  bool bdpair=false;
  if(chlep[a]*chlep[c]<0){
    if((abs(fllep[a])==abs(fllep[c]))){
      acpair=true;
    }
  }
  if(chlep[a]*chlep[d]<0){
    if((abs(fllep[a])==abs(fllep[d]))){
      adpair=true;
    }
  }
  if(acpair && adpair){
    printf("c and d is the same, both work on a");
  }
  if(acpair==false && adpair==false){
    printf("neither c or d work on a");
  }
  
  if(chlep[b]*chlep[c]<0){
    if((abs(fllep[b])==abs(fllep[c]))){
      bcpair=true;
    }
  }
  if(chlep[b]*chlep[d]<0){
    if((abs(fllep[b])==abs(fllep[d]))){
      bdpair=true;
    }
  }
  if(bcpair && bdpair){
    printf("c and d is the same, both work on b");
  }
  if(bcpair==false && bdpair==false){
    printf("neither c or d work on b");
  }
  TLorentzVector p1;
  p1.SetPtEtaPhiM(pt[a], eta[a], phi[a], e[a]);
  TLorentzVector p2;
  p2.SetPtEtaPhiM(pt[b], eta[b], phi[b], e[b]);
  TLorentzVector p3;
  p3.SetPtEtaPhiM(pt[c], eta[c], phi[c], e[c]);
  TLorentzVector p4;
  p4.SetPtEtaPhiM(pt[d], eta[d], phi[d], e[d]);
  double mass1;
  double mass2;
  if(acpair){
    mass1=(p1+p3).M();
  }
  if(adpair){
    mass1=(p1+p4).M();
  }
  if(bcpair){
    mass2=(p2+p3).M();
  }
  if(bdpair){
    mass2=(p2+p4).M();
  }
  std::pair <double,double> masses;
  masses = std::make_pair(mass1,mass2);
  return masses;
}
// void writeToFile(ROOT::RDF::RNode df, std::string treename, std::string outfilename, std::vector<std::string> &&good_cols){
    
//     df.Snapshot(treename,outfilename,good_cols);
    
//     return;
// }


