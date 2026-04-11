/**
 * CBDsim ROOT 검증: librootIO 로드 후 트리 점검
 *   export ROOTIO_LIB=$PWD/../rootIO/librootIO.so
 *   root -l -b -e '.x verify_root.C("test_0.root")' -q
 */
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TSystem.h"

namespace {
bool TryLoadRootIO() {
  if (const char* env = gSystem->Getenv("ROOTIO_LIB")) {
    if (gSystem->Load(env) == 0) {
      std::cout << "[verify_root] Loaded ROOTIO_LIB=" << env << std::endl;
      return true;
    }
  }
  const char* candidates[] = {"../rootIO/librootIO.so", "./librootIO.so", "../../rootIO/librootIO.so"};
  for (const char* path : candidates) {
    if (gSystem->AccessPathName(path)) continue;
    if (gSystem->Load(path) == 0) {
      std::cout << "[verify_root] Loaded " << path << std::endl;
      return true;
    }
  }
  std::cerr << "[verify_root] Set ROOTIO_LIB to build/rootIO/librootIO.so\n";
  return false;
}
TH1* GetHTemp() {
  if (auto* h = dynamic_cast<TH1*>(gDirectory->Get("htemp"))) return h;
  return dynamic_cast<TH1*>(gROOT->FindObject("htemp"));
}
}  // namespace

void verify_root(const char* filename = "test_0.root", const char* treename = "CBDsim") {
  if (!TryLoadRootIO()) return;
  TFile f(filename);
  if (!f.IsOpen() || f.IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }
  f.cd();
  TTree* t = nullptr;
  f.GetObject(treename, t);
  if (!t) {
    std::cerr << "No tree " << treename << std::endl;
    return;
  }
  std::cout << "Entries=" << t->GetEntries() << std::endl;
  t->Print();
  t->Scan("primaryEkin:hasSiPMTrig0:hasSiPMTrig1:siPMPhotonSumTrig0:siPMPhotonSumTrig1", "", "", 5);
  t->Scan("primaryEkin:primaryVx:primaryVy:primaryVz", "", "", 5);
  f.cd();
  t->Draw("towerT1.triggerNum", "", "goff");
  {
    TH1* h = GetHTemp();
    if (h) std::cout << "towerT1.triggerNum entries=" << h->GetEntries() << std::endl;
  }
  f.cd();
  t->Draw("Edeps.triggerNum", "", "goff");
  {
    TH1* h = GetHTemp();
    if (h) std::cout << "Edeps.triggerNum entries=" << h->GetEntries() << std::endl;
  }
  t->Show(0);
}
