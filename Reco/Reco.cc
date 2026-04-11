#include "CBDsimInterface.h"
#include "CBDsimRootInterface.h"

#include <iostream>

int main(int argc, char* argv[]) {
  std::string filenum = std::string(argv[1]);
  std::string filename = std::string(argv[2]);

  CBDsimRootInterface* drInterface = new CBDsimRootInterface(filename+"_"+filenum+".root");
  drInterface->set();

  unsigned int entries = drInterface->entries();
  while (drInterface->numEvt() < entries) {
    CBDsimInterface::CBDsimEventData evt;
    drInterface->read(evt);

    auto dumpTower = [](const char* label, const CBDsimInterface::CBDsimTowerData& tower) {
      std::cout << label << " triggerNum = " << tower.triggerNum << std::endl;
      for (const auto& sipm : tower.SiPMs) {
        std::cout << "  SiPM num = " << sipm.SiPMnum << " | Count = " << sipm.count << std::endl;
      }
    };
    dumpTower("towerT1 (T1)", evt.towerT1);
    dumpTower("towerT2 (T2)", evt.towerT2);

    for (auto edepItr = evt.Edeps.begin(); edepItr != evt.Edeps.end(); ++edepItr) {
      auto edep = *edepItr;

      // do something on the Edeps
      std::cout << "Trigger num = " << edep.triggerNum << " | Edep = " << edep.Edep << " (MeV)" << std::endl;
    }
  } // event loop

  return 0;
}
