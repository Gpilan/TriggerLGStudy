#ifndef CBDsimInterface_h
#define CBDsimInterface_h 1

#include <vector>
#include <utility>
#include <map>

class CBDsimInterface {
public:
  CBDsimInterface() {};
  ~CBDsimInterface() {};

  typedef std::pair<float,float> hitRange;
  typedef std::pair<int,int> hitXY;
  typedef std::map<hitRange, int> CBDsimTimeStruct;
  typedef std::map<hitRange, float> CBDsimWaveForm;
  typedef std::map<hitRange, int> CBDsimWavlenSpectrum;

  struct CBDsimSiPMData {
    CBDsimSiPMData() {};
    virtual ~CBDsimSiPMData() {};

    int count;
    int SiPMnum;
    int x;
    int y;
    CBDsimTimeStruct timeStruct;
    CBDsimWaveForm waveForm;
    CBDsimWavlenSpectrum wavlenSpectrum;
    /** Flattened time bins (ns edges) + counts for ROOT Draw; mirrors `timeStruct`. */
    std::vector<float> timeBinEdgeLow;
    std::vector<float> timeBinEdgeHigh;
    std::vector<int> timeBinCounts;
  };

  struct CBDsimPhysicalevent{
    CBDsimPhysicalevent() {};
    virtual ~CBDsimPhysicalevent() {};

    double x;
    double y;
    double z;
    int physicalID;
    int particleID;
    double energy;
  };

  struct CBDsimTowerData {
    CBDsimTowerData() {};
    virtual ~CBDsimTowerData() {};

    int numx;
    int numy;
    /** Legacy: tower index; proto: trigger index (0=T1, 1=T2). */
    int triggerNum;
    std::vector<CBDsimSiPMData> SiPMs;
  };

  struct CBDsimEdepData {
    CBDsimEdepData() {};
    virtual ~CBDsimEdepData() {};

    float Edep;
    int triggerNum;
  };

  struct CBDsimPhoton {
    CBDsimPhoton() {};
    virtual ~CBDsimPhoton() {};

    int opticalPhotonNumber;

  };

  struct CBDsimGenData {
    CBDsimGenData() {};
    virtual ~CBDsimGenData() {};

    float E;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float vt;
    int pdgId;
  };


  struct CBDsimEventData {
    CBDsimEventData() {};
    virtual ~CBDsimEventData() {};

    int event_number;
    /** Primary kinetic energy (MeV). */
    float primaryEkin = 0.f;
    /** Vertex position (mm). */
    float primaryVx = 0.f;
    float primaryVy = 0.f;
    float primaryVz = 0.f;
    /** Primary momentum direction (unit vector). */
    float primaryDirX = 0.f;
    float primaryDirY = 0.f;
    float primaryDirZ = 1.f;

    /** Per-trigger (0=T1, 1=T2) SiPM summary for quick Draw / cut. Legacy: only trig0 may fill. */
    int siPMPhotonSumTrig0 = 0;
    int siPMPhotonSumTrig1 = 0;
    int nSiPMChannelsTrig0 = 0;
    int nSiPMChannelsTrig1 = 0;
    int hasSiPMTrig0 = 0;
    int hasSiPMTrig1 = 0;
    /**
     * Sum of SiPM timeBinCounts over all channels in that trigger (same bin edges as first channel).
     * Optical photon arrival times are binned in SD; this is the merged histogram per event.
     */
    std::vector<float> timeMergedEdgeLowTrig0;
    std::vector<float> timeMergedEdgeHighTrig0;
    std::vector<int> timeMergedCountsTrig0;
    std::vector<float> timeMergedEdgeLowTrig1;
    std::vector<float> timeMergedEdgeHighTrig1;
    std::vector<int> timeMergedCountsTrig1;

    /** Trigger 1 (proto: T1, triggerNum 0). SiPM time bins: towerT1.SiPMs[].timeBinCounts / timeStruct. */
    CBDsimTowerData towerT1;
    /** Trigger 2 (proto: T2, triggerNum 1). */
    CBDsimTowerData towerT2;
    std::vector<CBDsimEdepData> Edeps;
    std::vector<CBDsimGenData> GenPtcs;
    std::vector<CBDsimPhoton> opticalPhotons;
    std::vector<CBDsimPhysicalevent> totPhysicals;
  };

};

#endif
