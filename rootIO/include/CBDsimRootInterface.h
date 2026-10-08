#ifndef CBDsimRootInterface_h
#define CBDsimRootInterface_h 1

#include "CBDsimInterface.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

class CBDsimRootInterface {
public:
  CBDsimRootInterface(const std::string& filename);
  ~CBDsimRootInterface();

  void fill(const CBDsimInterface::CBDsimEventData* evt);
  void GetChain(const std::string& treename);
  void read(CBDsimInterface::CBDsimEventData& evt);
  void create();
  void set();
  void write();
  void close();

  TTree* getTree();
  unsigned int entries() { return fTree->GetEntries(); }
  unsigned int numEvt() { return fNumEvt; }

private:
  void init();
  void PrepareChain();
  
  TChain* fChain = nullptr;
  TFile* fFile = nullptr;
  TTree* fTree = nullptr;
  std::string fFilename;
  CBDsimInterface::CBDsimEventData* fEventData = nullptr;
  unsigned int fNumEvt;
};

#endif
