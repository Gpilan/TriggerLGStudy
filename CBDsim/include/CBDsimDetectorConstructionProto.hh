#ifndef CBDsimDetectorConstructionProto_h
#define CBDsimDetectorConstructionProto_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "CBDsimMaterials.hh"

class G4LogicalVolume;

/**
 * Prototype geometry (trigger/LG study): plastic scint slab + tessellated light guide.
 * Beam along +z crosses the thin (5 mm) z extent; world is vacuum (G4_Galactic); Al foil
 * on scint; LG on the 40x5 mm face (normal -y), extending to -y; SiPM at LG tip (axis y).
 * readout face (same SD naming as legacy back plane). Units: mm.
 *
 * Only the world logical volume is G4_Galactic; two trigger assemblies (shared LVs, copy 0/1) are placed
 * directly in the world (T2 is Rz(+90 deg)). Along world z, the gap between the two scintillator slabs
 * is kTrig12TileFaceGapZ (1 mm); center separation is 2*kHzThin + that gap so tiles do not overlap.
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
