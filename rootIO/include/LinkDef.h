#include <vector>
#include <map>
#include <utility>

#if defined(__CLING__) || defined(__CINT__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class std::pair<float,float>+;
#pragma link C++ class std::map<std::pair<float,float>, int>+;
#pragma link C++ class std::map<std::pair<float,float>, float>+;

#pragma link C++ struct CBDsimInterface::CBDsimSiPMData+;
#pragma link C++ struct CBDsimInterface::CBDsimPhoton+;
#pragma link C++ struct CBDsimInterface::CBDsimPhysicalevent+;
#pragma link C++ struct CBDsimInterface::CBDsimTowerData+;
#pragma link C++ struct CBDsimInterface::CBDsimEdepData+;
#pragma link C++ struct CBDsimInterface::CBDsimGenData+;
#pragma link C++ struct CBDsimInterface::CBDsimEventData+;

#endif
