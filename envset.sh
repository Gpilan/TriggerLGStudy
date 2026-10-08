source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
source /cvmfs/geant4.cern.ch/geant4/11.4/x86_64-el9-gcc11-optdeb-MT/CMake-setup.sh

# Select the installed Geant4 config and its external Xerces dependency.
export Geant4_DIR=/cvmfs/geant4.cern.ch/geant4/11.4/x86_64-el9-gcc11-optdeb-MT/lib64/cmake/Geant4
export CMAKE_PREFIX_PATH="/cvmfs/sft.cern.ch/lcg/releases/XercesC/3.3.0-ff50f/x86_64-el9-gcc11-opt:${CMAKE_PREFIX_PATH:-}"
export LD_LIBRARY_PATH="/cvmfs/sft.cern.ch/lcg/releases/XercesC/3.3.0-ff50f/x86_64-el9-gcc11-opt/lib:${LD_LIBRARY_PATH:-}"

# Override dataset variables inherited from LCG with the Geant4 11.4 paths.
export G4NEUTRONHPDATA="$GEANT4_DATA_DIR/G4NDL4.7.1"
export G4LEDATA="$GEANT4_DATA_DIR/G4EMLOW8.8"
export G4LEVELGAMMADATA="$GEANT4_DATA_DIR/PhotonEvaporation6.1.2"
export G4RADIOACTIVEDATA="$GEANT4_DATA_DIR/RadioactiveDecay6.1.2"
export G4PARTICLEXSDATA="$GEANT4_DATA_DIR/G4PARTICLEXS4.2"
export G4PIIDATA="$GEANT4_DATA_DIR/G4PII1.3"
export G4REALSURFACEDATA="$GEANT4_DATA_DIR/RealSurface2.2"
export G4SAIDXSDATA="$GEANT4_DATA_DIR/G4SAIDDATA2.0"
export G4ABLADATA="$GEANT4_DATA_DIR/G4ABLA3.3"
export G4INCLDATA="$GEANT4_DATA_DIR/G4INCL1.3"
export G4ENSDFSTATEDATA="$GEANT4_DATA_DIR/G4ENSDFSTATE3.0"
export G4CHANNELINGDATA="$GEANT4_DATA_DIR/G4CHANNELING2.0"
export G4PARTICLEHPDATA="$GEANT4_DATA_DIR/G4TENDL1.4"
export G4NUDEXLIBDATA="$GEANT4_DATA_DIR/G4NUDEXLIB1.0"
export G4URRPTDATA="$GEANT4_DATA_DIR/G4URRPT1.1"

export TRIGGER_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CBDsim_STATS="${TRIGGER_ROOT}/build/CBDsim/stats"
