#ifndef example_h
#define example_h

#include "ROOT/RDF/RInterface.hxx"
#include <ROOT/RDataSource.hxx>
#include <ROOT/RCsvDS.hxx>
#include <iostream>
#include "TLorentzVector.h"
#include "TParameter.h"
#include <ROOT/RVec.hxx>
using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;

VecF_t ExampleInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecI_t& charge, VecF_t& e);
#endif
