#ifndef Cfunctions_h
#define Cfunctions_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TString.h"
#include <ROOT/RVec.hxx>
#include <map>
#include <utility>
using Vec2_t = const ROOT::VecOps::RVec<float>;
using Vec_t = const ROOT::VecOps::RVec<int>;


std::map<int,float> btagWP;
std::map<int,std::pair<TString,TString>> MCcat;
std::vector< std::pair<TString,TString> > pairs;

template <typename T>

T InvariantMassStdVector(std::vector<T>& pt, std::vector<T>& eta, std::vector<T>& phi, std::vector<T>& energy);
void  help();
bool  myfilter(float x);
int  isOS(const ROOT::VecOps::RVec<int>& chlep);
int  isSF(Vec_t& fllep);
float ComputeInvariantMass(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float calcMT2(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
float ptllboost(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
float costhetastar(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float deltaPhi_ll(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e);
float deltaPhi_metl(Vec2_t& pt, Vec2_t& eta, Vec2_t& phi, Vec2_t& e, Float_t met_et, Float_t met_phi);
bool  checkPt(Vec2_t& pt, float cut1, float cut2);
std::vector<std::string> DropColumns(std::vector<std::string> &&good_cols);
float getVar(Vec2_t& var, int idx);
int getTypeTimesCharge(Vec_t& chlep, Vec_t& fllep, int idx);
int countBJets(Vec2_t& pt, Vec2_t& tagger, float ptcut, int WP = 85);
int countJets(Vec2_t& pt, float ptcut);
std::pair<const char *,const char *> getMCMetaData(int dsid);
void readMCCat(std::string infile = "./Signal_samples_13TeV.txt");
std::pair<const char *, const char *> getDataMetaData(int dsid);
double compareEvents(Int_t rnum_orig, const ROOT::VecOps::RVec<ULong64_t>& rnum_wgt, const ROOT::VecOps::RVec<Double_t>& value);
#endif
