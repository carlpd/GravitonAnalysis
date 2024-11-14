#define helperFunctions_cxx


#include <ROOT/RVec.hxx>
#include "TLorentzVector.h"

using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;

// Compute the invariant mass for oppositely charged pairs of objects
// Note that this function "sees" one event at a time, so the entries in each vector are for a single event
ROOT::RVec<float> ExampleInvariantMass(VecF_t& pt, VecF_t& eta, VecF_t& phi, VecI_t& charge, float mass) {
  // Construct oppositely charged pairs. 
  // This is done by constructing index pairs where each index is the same length as the vectors
  // and then checking whether they are oppositely charged. Selected index pairs are then used later
  // to select the correct values from the other vectors, to calculate the mass
  std::vector<std::pair<unsigned int, unsigned int> > indices{};
  unsigned int length{(unsigned int)charge.size()};
  for (unsigned int i=0; i<length; ++i) {
    for (unsigned int j=i+1; j<length; ++j) {
      if (charge[i]!=charge[j]) indices.emplace_back(std::make_pair(i,j));
    }
  }
  // Calculate mass for each pair
  ROOT::RVec<float> masses{};
  for (unsigned int k=0; k<indices.size(); ++k) {
    TLorentzVector p1;
    TLorentzVector p2;
    unsigned int index1{indices[k].first};
    unsigned int index2{indices[k].second};
    p1.SetPtEtaPhiM(pt[index1], eta[index1], phi[index1], mass);
    p2.SetPtEtaPhiM(pt[index2], eta[index2], phi[index2], mass);
    masses.push_back((p1 + p2).M());
  }
  return masses;
}



