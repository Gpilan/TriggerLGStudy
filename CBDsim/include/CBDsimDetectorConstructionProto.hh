#ifndef CBDsimDetectorConstructionProto_h
#define CBDsimDetectorConstructionProto_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "CBDsimMaterials.hh"

class G4LogicalVolume;

/**
 * Prototype geometry (trigger/LG study): plastic scint slab + tessellated light guide.
 * Beam along +z crosses the thin (5 mm) z extent; air envelope + Al foil
 * on scint; LG on the 40x5 mm face (normal -y), extending to -y; SiPM at LG tip (axis y).
 * readout face (same SD naming as legacy back plane). Units: mm.
 *
 * Two triggers share protoAirEnvLog, placed inside one outer air volume (protoOuterAirEnvLog)
 * sized to contain both (T2 is Rz(+90 deg)); optional safety margin on that outer box.
 * Two placements of shared protoAirEnvLog; env centers at +/-(hzEnv + gap/2). Scint centers match env centers.
 * Effective gap between scint inner faces = 2*kEnvZHalfSep - 2*kHzThin (large vs 1 mm — required for LG in mother).
 *
 * Switch in CBDsim.cc: use this class instead of CBDsimDetectorConstruction.
 */
class CBDsimDetectorConstructionProto : public G4VUserDetectorConstruction {
public:
  CBDsimDetectorConstructionProto();
  ~CBDsimDetectorConstructionProto() override;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

private:
  void DefineMaterials();
  G4Material* FindMaterial(const G4String& name) { return fMaterials->GetMaterial(name); }
  G4OpticalSurface* FindSurface(const G4String& name) { return fMaterials->GetOpticalSurface(name); }

  CBDsimMaterials* fMaterials;
  G4VisAttributes* fVisWorld;
  G4VisAttributes* fVisScint;
  G4VisAttributes* fVisLG;
  G4VisAttributes* fVisFoil;
  G4VisAttributes* fVisSiPM;
  G4LogicalVolume* fProtoWaferLog = nullptr;
};

#endif
