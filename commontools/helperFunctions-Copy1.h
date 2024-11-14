#ifndef helperFunctions_h
#define helperFunctions_h

#include "ROOT/RDF/RInterface.hxx"
#include <ROOT/RDataSource.hxx>
#include <ROOT/RCsvDS.hxx>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <ctime>
#include <dirent.h>
#include "TLorentzVector.h"
#include "TParameter.h"
#include <ROOT/RVec.hxx>
using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;
using VecB_t = const ROOT::VecOps::RVec<bool>;
using VecD_t = const ROOT::RVec<double>&;
using VecUI_t = const ROOT::RVec<UInt_t>&;


//TCanvas c("c","x hist");
std::string progressBar;
float nEvents;
int everyN;
ROOT::RDataFrame metadata_rdf(0);


clock_t c_begin = 0;
clock_t c_end = 0;
// void updateHistogram(TH1D &h_){		     
//   c.cd();
//   h_.Draw();
//   c.Update();    
// };

Int_t getElecPdgID(){return 11;}
Int_t getMuonPdgID(){return 13;}

Float_t getElecMass(){return 0.511;}
Float_t getMuonMass(){return 105.66;}
/*
void readMetaData(std::string infile){
    metadata_rdf = ROOT::RDF::FromCSV(infile);
    std::cout<<"Number of entris in metadata ttree is "<<metadata_rdf.Count().GetValue()<<std::endl;
}
*/
void dumpMetaData(){
    if(!metadata_rdf.Count().GetValue()){
       std::cout<<"ERROR \t RDataFrame for metadata has not been loaded. Call ROOT.readMetaData(in_csv_file)."<<std::endl;
       return;
    }
    auto colNames = metadata_rdf.GetColumnNames();
    // Print columns' names
    for (auto &&colName : colNames){
       auto colType = metadata_rdf.GetColumnType(colName);
       // Print column type
       std::cout << "Column " << colName << " has type " << colType << std::endl;
    }
}


Double_t getXsec(Long64_t dsid){
    //std::cout<<"getXsec::before"<<std::endl;
    auto var = metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Take<Double_t>("crossSection");
    if(var->size() <= 0){
        std::cout<<"Did not find xsec for "<<dsid<<". Returning -1"<<std::endl;
        return -1;
    }
    return var->at(0);
    //return metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Max("crossSection").GetValue();
}
Double_t getKfac(Long64_t dsid){
    //std::cout<<"getKfac::before"<<std::endl;
    auto var = metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Take<Double_t>("kFactor");
    if(var->size() <= 0){
        std::cout<<"Did not find kFactor for "<<dsid<<". Returning -1"<<std::endl;
        return -1;
    }
    return var->at(0);
    //return metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Max("kFactor").GetValue();
}
Double_t getFilterEff(Long64_t dsid){
    //std::cout<<"getFilterEff::before"<<std::endl;
    auto var = metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Take<Double_t>("genFiltEff");
    if(var->size() <= 0){
        std::cout<<"Did not find genFiltEff for "<<dsid<<". Returning -1"<<std::endl;
        return -1;
    }
    return var->at(0);
    //return metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Max("genFiltEff").GetValue();
}
Double_t getSOW(Long64_t dsid){
    //std::cout<<"getSOW::before"<<std::endl;
    auto var = metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Take<Double_t>("sumofweights");
    if(var->size() <= 0){
        std::cout<<"Did not find sumofweights for "<<dsid<<". Returning -1"<<std::endl;
        return -1;
    }
    return var->at(0);
    //return metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Max("sumofweights").GetValue();
}
// Double_t getScalefactor(Double_t xsec, Double_t kfac, Double_t filtereff, Double_t sow, Double_t lumi){
//     //std::cout<<"getScalefactor::before"<<std::endl;
//     if(sow >0){
//        return (xsec*kfac*filtereff*lumi)/sow;   
//     }else{
//        return -1;
//     }
// }
Double_t getScalefactor(Long64_t dsid, Double_t lumi){
  //std::cout<<"getScalefactor::before"<<std::endl;
  auto rdf = metadata_rdf.Filter(Form("dataset_number == %lld",dsid));
  Double_t xsec = rdf.Take<Double_t>("crossSection")->at(0);
  Double_t kfac = rdf.Take<Double_t>("kFactor")->at(0);
  Double_t filtereff = rdf.Take<Double_t>("genFiltEff")->at(0);
  Double_t sow = rdf.Take<Double_t>("sumofweights")->at(0);
  if(sow >0){
    return (xsec*kfac*filtereff*lumi)/sow;   
  }else{
    return -1;
  }
}

std::string getCategory(Long64_t dsid){
    auto category = metadata_rdf.Filter(Form("dataset_number == %lld",dsid)).Take<std::string>("category");
    //
    if(category->size() <= 0){
        std::cout<<"Did not find category for "<<dsid<<". Returning UNKNOWN"<<std::endl;
        return "UNKNOWN";
    }
    return category->at(0);
}


void setRunParameters(float n, int en){
  nEvents = n;
  everyN = en;
}

void printProgressBar(TH1D &h_){
  double elapsed_secs = 0;
  double ev_per_sec = 0;
  double time_left = 0.0;
  double minutesRemainder = 0.0;
  double secondsRemainder = 0.0;
  int hours = 0;
  int minutes = 0;
  int seconds = 0;
  double all_nev = progressBar.size()*everyN;
  c_end = clock();
  elapsed_secs = double(c_end - c_begin) / CLOCKS_PER_SEC;
  if(elapsed_secs > 0){
    ev_per_sec = ((double)everyN)/elapsed_secs;
    if(ev_per_sec > 0){
      time_left = (nEvents-all_nev)/ev_per_sec;
      time_left = time_left/(60.*60);
      hours = time_left;
      minutesRemainder = (time_left - hours) * 60;
      minutes = minutesRemainder;
      secondsRemainder = (minutesRemainder - minutes) * 60;
      seconds = secondsRemainder;
    }
    std::cout<<"Events/sec = "<<ev_per_sec<<"). "<<"Estimated time left = "<<hours<<"h"<<minutes<<"m"<<seconds<<"s"<<std::endl;
  }
  std::mutex barMutex; // Only one thread at a time can lock a mutex. Let's use this to avoid concurrent printing.
  // Magic numbers that yield good progress bars for nSlots = 1,2,4,8
  const auto barWidth = nEvents / everyN;
  std::lock_guard<std::mutex> l(barMutex); // lock_guard locks the mutex at construction, releases it at destruction
  progressBar.push_back('#');
  // re-print the line with the progress bar
  std::cout << "\r[" << std::left << std::setw(barWidth) << progressBar << ']' << std::flush;
  c_begin = clock();
}

int Zlep1 = -1;
int Zlep2 = -1;
int Wlep1 = -1;

std::vector<int> checkJVT(ROOT::VecOps::RVec<char> NNJVT);
std::vector<float> getTaggerProb(ROOT::VecOps::RVec<UInt_t> key, ROOT::VecOps::RVec<UInt_t> index, ROOT::VecOps::RVec<float> tagger);

std::pair <double,double> getLeptonsFromZ(VecI_t chlep, VecI_t& fllep, VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
bool  myfilter(float x);
bool  isOS(const ROOT::VecOps::RVec<int>& chlep);
bool  isSF(VecI_t& fllep);
bool isEE(VecI_t& fllep);
bool isMM(VecI_t& fllep);
double getSF(VecF_t& sf);
int flavourComp(VecI_t& fllep);//, int nlep);
float ComputeInvariantMass2L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& massfloat ComputeInvariantMass2L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass, int i, int j));
float ComputeInvariantMass4L(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& mass);
float calcMT2(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& m, Float_t met_et, Float_t met_px, Float_t met_py, int idx1=0, int idx2=1);
float ptllboost(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
float costhetastar(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
float deltaPhi_ll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e);
float deltaPhi_metl(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
bool  checkPt(VecF_t& pt, float cut1, float cut2);
float deltaPhi_metll(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
std::pair <int,int> num_bl_sg_lep(VecF_t& pt, VecF_t& eta, VecI_t& fllep, VecB_t passOR, VecB_t passLOOSE, VecB_t passMEDIUM, VecB_t passBL, VecF_t z0sinth, VecB_t ISO, VecF_t d0sig, VecF_t passTIGHT);
Double_t getMetRel(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecF_t& e, Float_t met_et, Float_t met_phi);
int getZlep1(){return Zlep1;}
int getZlep2(){return Zlep2;}
int getWlep1(){return Wlep1;}
Float_t getLumiSF(Int_t randrnum);
int checkTriggerMatch(VecF_t& pt, VecB_t trigmatch, Float_t th);
double getWeight(ULong64_t evnum);
ROOT::RDF::RResultPtr<ROOT::VecOps::RVec<double>> getWeightVec(ROOT::RDF::RNode wgtdf); 
ROOT::RDF::RNode defineNewColumn(ROOT::RDF::RNode wgtdf, ROOT::RDF::RNode df);
//void writeToFile(ROOT::RDF::RNode df, std::string treename, std::string outfilename, std::vector<std::string> &&good_cols);
//Double_t getSumOfWeights(ROOT::RDF::RNode df, int dsid);
std::vector<int>  getUniqueDSIDs(VecI_t alldsid);
ROOT::RVec<float> getVector(VecF_t& inp1, VecF_t& inp2, Float_t m1 = 0.0 , Float_t m2 = 0.0);
ROOT::RVec<float> getVector(VecF_t& inp1, Float_t m1);
Double_t getSumOfWeights(ROOT::VecOps::RVec<std::string> name, VecD_t& sow, std::string key = "PHYSLITEKernel");
Double_t getScaleFactor(UInt_t this_did, VecI_t dids, VecD_t sfs);
//std::vector<std::string> DropColumns(std::vector<std::string> &&good_cols, const std::vector<std::string> blacklist);
#endif
