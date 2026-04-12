#include "CBDsimMaterials.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <algorithm>
#include <cmath>

CBDsimMaterials* CBDsimMaterials::fInstance = 0;

CBDsimMaterials::CBDsimMaterials() {
  fNistMan = G4NistManager::Instance();
  CreateMaterials();
}

CBDsimMaterials::~CBDsimMaterials() {}

CBDsimMaterials* CBDsimMaterials::GetInstance() {
  if (fInstance==0) fInstance = new CBDsimMaterials();

  return fInstance;
}

G4Material* CBDsimMaterials::GetMaterial(const G4String matName) {
  G4Material* mat = fNistMan->FindOrBuildMaterial(matName);

  if (!mat) mat = G4Material::GetMaterial(matName);
  if (!mat) {
    std::ostringstream o;
    o << "Material " << matName << " not found!";
    G4Exception("CBDsimMaterials::GetMaterial","",FatalException,o.str().c_str());
  }

  return mat;
}

G4OpticalSurface* CBDsimMaterials::GetOpticalSurface(const G4String surfName) {
  if (surfName=="SiPMSurf") return fSiPMSurf;
  else if (surfName=="FilterSurf") return fFilterSurf;
  else if (surfName=="AluminumSurf") return fAluminumSurf;
  else if (surfName=="AirSurf") return fAirSurf;
  else {
    std::ostringstream o;
    o << "OpticalSurface " << surfName << " not found!";
    G4Exception("CBDsimMaterials::GetOpticalSurface","",FatalException,o.str().c_str());
  }

  return nullptr;
}

void CBDsimMaterials::CreateMaterials() {
  fNistMan->FindOrBuildMaterial("G4_Galactic");
  fNistMan->FindOrBuildMaterial("G4_AIR");

  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z=1., a=1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z=6., a=12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z=7., a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z=8., a=16.00*g/mole);
  G4Element* F  = new G4Element("Fluorine",symbol="F" , z=9., a=18.9984*g/mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.085*g/mole);
  G4Element* Y  = new G4Element("Yttrium" ,symbol="Y" , z= 39., a= 88.905*g/mole);
  G4Element* Ce = new G4Element("Cerium"  ,symbol="Ce", z= 58., a= 140.116*g/mole);
  G4Element* Lu = new G4Element("Lutetium",symbol="Lu", z= 71., a= 174.9668*g/mole);
  G4Element* W  = new G4Element("Tungsten",symbol="W" , z=74.,a=183.8400*g/mole);
  G4Element* Pb = new G4Element("Lead"    ,symbol="Pb", z= 82., a=207.2*g/mole);
  G4Element* Al = new G4Element("Aluminum",symbol="Al",z=13,a=26.982*g/mole);



  fAl = new G4Material("Aluminum",z=13,a=26.982*g/mole,density=2.70*g/cm3);
  fCu = new G4Material("Copper", z=29., a=63.546*g/mole, density=8.96*g/cm3);
  // Unique name: NIST/FindOrBuildMaterial("Silicon") can resolve to a different G4Material without our MPT.
  fSi = new G4Material("SiPM_WaferSilicon", z=14., a=28.09*g/mole, density=2.33*g/cm3);

  fVacuum = G4Material::GetMaterial("G4_Galactic");
  fAir = G4Material::GetMaterial("G4_AIR");

  fPWO = new G4Material("PWO",density=8.28*g/cm3, ncomponents=3);
  fPWO -> AddElement(Pb,natoms=1);
  fPWO -> AddElement(W, natoms=1);
  fPWO -> AddElement(O, natoms=4);

  fLSO = new G4Material("LSO", density=7.40*g/cm3, ncomponents=4);
  fLSO->AddElement(Lu, 76.25*perCent);
  fLSO->AddElement(Ce, 0.2*perCent);
  fLSO->AddElement(Si, 6.12*perCent);
  fLSO->AddElement(O , 17.43*perCent);

  fLYSO = new G4Material("LYSO", density=7.25*g/cm3, ncomponents=5);
  fLYSO->AddElement(Lu, 73.856*perCent);
  fLYSO->AddElement(Y , 1.975*perCent);
  fLYSO->AddElement(Si, 6.240*perCent);
  fLYSO->AddElement(O , 17.773*perCent);
  fLYSO->AddElement(Ce, 0.156*perCent);

  fFluoPoly = new G4Material("FluorinatedPolymer", density=1.43*g/cm3, ncomponents=2);
  fFluoPoly->AddElement(C, 2);
  fFluoPoly->AddElement(F, 2);

  fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
  fGlass->AddElement(C, 91.533*perCent);
  fGlass->AddElement(H, 8.467*perCent);

  fPS = new G4Material("Polystyrene", density=1.05*g/cm3, ncomponents=2);
  fPS->AddElement(C, natoms=8);
  fPS->AddElement(H, natoms=8);

  fPMMA = new G4Material("PMMA", density= 1.19*g/cm3, ncomponents=3);
  fPMMA->AddElement(C, natoms=5);
  fPMMA->AddElement(H, natoms=8);
  fPMMA->AddElement(O, natoms=2);

  fProtoLGMatchScint = new G4Material("ProtoLG_MatchScint", density= 1.19*g/cm3, ncomponents=3);
  fProtoLGMatchScint->AddElement(C, natoms=5);
  fProtoLGMatchScint->AddElement(H, natoms=8);
  fProtoLGMatchScint->AddElement(O, natoms=2);

  fGelatin = new G4Material("Gelatin", density=1.27*g/cm3, ncomponents=4);
  fGelatin->AddElement(C, natoms=102);
  fGelatin->AddElement(H, natoms=151);
  fGelatin->AddElement(N, natoms=31);
  fGelatin->AddElement(O, natoms=39);


  ///--- Material property tables for absorber materials ---
  G4MaterialPropertiesTable* mpPWO;
  G4MaterialPropertiesTable* mpLSO;
  G4MaterialPropertiesTable* mpLYSO;

///--- Material property tables for fiber materials ---
  G4MaterialPropertiesTable* mpAir;
  G4MaterialPropertiesTable* mpPS;
  G4MaterialPropertiesTable* mpPMMA;
  G4MaterialPropertiesTable* mpFluoPoly;
  G4MaterialPropertiesTable* mpGlass;
  G4MaterialPropertiesTable* mpSiPM;
  G4MaterialPropertiesTable* mpFilter;
  G4MaterialPropertiesTable* mpFilterSurf;
  G4MaterialPropertiesTable* mpAl;
  G4MaterialPropertiesTable* mpAlSurf;
  G4MaterialPropertiesTable* mpAirSurf;

  G4double opEn[] = { // from 900nm to 300nm with 25nm step
    1.37760*eV, 1.41696*eV, 1.45864*eV, 1.50284*eV, 1.54980*eV, 1.59980*eV, 1.65312*eV, 1.71013*eV,
    1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV, 2.15625*eV, 2.25426*eV, 2.36160*eV,
    2.47968*eV, 2.61019*eV, 2.75520*eV, 2.91728*eV, 3.09960*eV, 3.30625*eV, 3.54241*eV, 3.81490*eV, 4.13281*eV
  };

  const G4int nEnt = sizeof(opEn) / sizeof(G4double);
//////http://nuclear.korea.ac.kr/~lamps/geant4/G4MaterialPropertiesIndex_8hh_source.html

G4double RI_Air[nEnt]; std::fill_n(RI_Air,nEnt,1.0);
mpAir = new G4MaterialPropertiesTable();
mpAir->AddProperty("RINDEX",opEn,RI_Air,nEnt);
fAir->SetMaterialPropertiesTable(mpAir);

//G4double refl_Air[nEnt]; std::fill_n(refl_Air,nEnt,0.);
mpAirSurf = new G4MaterialPropertiesTable();
//mpAirSurf->AddProperty("REFLECTIVITY",opEn,refl_Air,nEnt);
fAirSurf = new G4OpticalSurface("AirSurf", glisur, polished, dielectric_dielectric);

fAirSurf->SetMaterialPropertiesTable(mpAirSurf);


/*
  mpSiPM = new G4MaterialPropertiesTable();
  mpSiPM->AddProperty("REFLECTIVITY",opEn,refl_SiPM,nEnt);
  mpSiPM->AddProperty("EFFICIENCY",opEn,eff_SiPM,nEnt);
  fSiPMSurf = new G4OpticalSurface("SiPMSurf",glisur,polished,dielectric_metal);
  fSiPMSurf->SetMaterialPropertiesTable(mpSiPM);
*/
// G4double RI_Al[nEnt]={//https://refractiveindex.info/?shelf=3d&book=metals&page=aluminium
//   2.0412,2.2672,2.4942,2.6798,2.6717,2.5128,2.2317,1.9126,
//   1.6664,1.4831,1.3456,1.2366,1.1419,1.0530,0.96521,0.87691,
//   0.78834,0.70111,0.61722,0.53817,0.46582,0.40059,0.34288,0.29204,0.24765};
// G4double refl_Al[nEnt]; std::fill_n(refl_Al, nEnt,1.0);
G4double refl_Al[nEnt]={//https://eng.libretexts.org/Bookshelves/Materials_Science/Supplemental_Modules_(Materials_Science)/Optical_Properties/Metallic_Reflection
  0.894, 0.887, 0.880, 0.873, 0.866, 0.876, 0.887, 0.897,
  0.908, 0.910, 0.912, 0.914, 0.916, 0.918, 0.920, 0.921,
  0.923, 0.925, 0.927, 0.928, 0.930, 0.928, 0.925, 0.923, 0.92};
mpAlSurf = new G4MaterialPropertiesTable();
// G4double refl_Al[nEnt]; std::fill_n(refl_Al,nEnt,1);
mpAlSurf->AddProperty("REFLECTIVITY", opEn, refl_Al, nEnt);
//mpAlSurf->AddProperty("RINDEX",opEn,RI_Al,nEnt);
fAluminumSurf = new G4OpticalSurface("AluminumSurf", unified, polished, dielectric_metal);
fAluminumSurf->SetMaterialPropertiesTable(mpAlSurf);


/*G4double RI_Al[nEnt]={//https://refractiveindex.info/?shelf=3d&book=metals&page=aluminium
  2.0412,2.2672,2.4942,2.6798,2.6717,2.5128,2.2317,1.9126,
  1.6664,1.4831,1.3456,1.2366,1.1419,1.0530,0.96521,0.87691,
  0.78834,0.70111,0.61722,0.53817,0.46582,0.40059,0.34288,0.29204,0.24765};*/
//G4double Abslength_Al[nEnt]; std::fill_n(Abslength_Al,nEnt,100*cm);

//};
// G4double RI_Al[nEnt]; std::fill_n(RI_Al,nEnt,2);
mpAl = new G4MaterialPropertiesTable();
// mpAl->AddProperty("RINDEX",opEn,RI_Al,nEnt);
//mpAl->AddProperty("ABSLENGTH",opEn,Abslength_Al,nEnt);
//mpAl->AddProperty("ABSLENGTH",opEn,abslength_Al,nEnt);
// mpAl->AddProperty("REFLECTIVITY",opEn,refl_Al,nEnt);
fAl->SetMaterialPropertiesTable(mpAl);


  //--- LSO ---
G4double RI_LSO[nEnt]; std::fill_n(RI_LSO, nEnt, 1.82);
G4double Abslength_LSO[nEnt]; std::fill_n(Abslength_LSO, nEnt, 1.14*cm);
G4double scintFast_LSO[nEnt]={//이거 근거 더 찾기
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.03, 0.07, 0.12, 0.16, 0.32,
  0.50, 0.68, 0.92, 0.97, 0.50, 0.07, 0.00, 0.00, 0.00};
mpLSO = new G4MaterialPropertiesTable();
mpLSO->AddProperty("RINDEX",opEn,RI_LSO,nEnt);
mpLSO->AddProperty("ABSLENGTH",opEn,Abslength_LSO,nEnt);
mpLSO->AddConstProperty("SCINTILLATIONYIELD",30000./MeV);
mpLSO->AddConstProperty("RESOLUTIONSCALE",1.0);//이것도
mpLSO->AddProperty("SCINTILLATIONCOMPONENT1",opEn,scintFast_LSO,nEnt,true);
mpLSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1",40*ns,true); // true: linear interpolation 03.25
fLSO->SetMaterialPropertiesTable(mpLSO);

//--- LYSO ---
G4double RI_LYSO[nEnt]; std::fill_n(RI_LYSO, nEnt, 1.81);
G4double Abslength_LYSO[nEnt]; std::fill_n(Abslength_LYSO, nEnt, 42*cm);
G4double scintFast_LYSO[nEnt] = {
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
  0.00, 0.00, 0.00, 0.03, 0.07, 0.12, 0.16, 0.32,
  0.50, 0.68, 0.92, 0.97, 0.50, 0.07, 0.00, 0.00, 0.00
};
mpLYSO = new G4MaterialPropertiesTable();
mpLYSO->AddProperty("RINDEX",opEn,RI_LYSO,nEnt);
mpLYSO->AddProperty("ABSLENGTH",opEn,Abslength_LYSO,nEnt);
mpLYSO->AddProperty("SCINTILLATIONCOMPONENT1",opEn,scintFast_LYSO,nEnt,true);
mpLYSO->AddConstProperty("SCINTILLATIONYIELD",33200./MeV);
mpLYSO->AddConstProperty("RESOLUTIONSCALE",2.0);
mpLYSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1",36*ns,true);
fLYSO->SetMaterialPropertiesTable(mpLYSO);

//--- PbWO4 ---
G4double RI_PWO[nEnt]; std::fill_n(RI_PWO,nEnt,2.16);
G4double Abslength_PWO[nEnt]; std::fill_n(Abslength_PWO,nEnt,0.92*cm);
//G4double scintFast_PWO[nEnt]={};
//G4double scintSlow_PWO[nEnt]={};
mpPWO = new G4MaterialPropertiesTable();
mpPWO->AddProperty("RINDEX",opEn,RI_PWO,nEnt);
mpPWO->AddProperty("ABSLENGTH",opEn,Abslength_PWO,nEnt);
//mpPWO->AddProperty("SCINTILLATIONCOMPONENT1",opEn,scintFast_PWO,nEnt);
//mpPWO->AddProperty("SCINTILLATIONCOMPONENT2",opEn,scintSlow_PWO,nEnt);
mpPWO->AddConstProperty("SCINTILLATIONYIELD",400/MeV);
mpPWO->AddConstProperty("RESOLUTIONSCALE",1.0);
mpPWO->AddConstProperty("SCINTILLATIONTIMECONSTANT1",6*ns,true);
mpPWO->AddConstProperty("SCINTILLATIONTIMECONSTANT2",30*ns,true);
fPWO->SetMaterialPropertiesTable(mpPWO);


  G4double RI_PMMA[nEnt] = {
    1.48329, 1.48355, 1.48392, 1.48434, 1.48467, 1.48515, 1.48569, 1.48628,
    1.48677, 1.48749, 1.48831, 1.48899, 1.49000, 1.49119, 1.49219, 1.49372,
    1.49552, 1.49766, 1.49953, 1.50252, 1.50519, 1.51000, 1.51518, 1.52182, 1.53055
  };
  G4double AbsLen_PMMA[nEnt] = {
    0.414*m, 0.543*m, 0.965*m, 2.171*m, 2.171*m, 3.341*m, 4.343*m, 1.448*m,
    4.343*m, 14.48*m, 21.71*m, 8.686*m, 28.95*m, 54.29*m, 43.43*m, 48.25*m,
    54.29*m, 48.25*m, 43.43*m, 28.95*m, 21.71*m, 4.343*m, 2.171*m, 0.869*m, 0.434*m
  };
  mpPMMA = new G4MaterialPropertiesTable();
  mpPMMA->AddProperty("RINDEX",opEn,RI_PMMA,nEnt);
  mpPMMA->AddProperty("ABSLENGTH",opEn,AbsLen_PMMA,nEnt);
  fPMMA->SetMaterialPropertiesTable(mpPMMA);

  G4double RI_FluoPoly[nEnt]; std::fill_n(RI_FluoPoly, nEnt, 1.42);
  mpFluoPoly = new G4MaterialPropertiesTable();
  mpFluoPoly->AddProperty("RINDEX",opEn,RI_FluoPoly,nEnt);
  fFluoPoly->SetMaterialPropertiesTable(mpFluoPoly);

  G4double RI_PS[nEnt] = {
    1.57483, 1.57568, 1.57644, 1.57726, 1.57817, 1.57916, 1.58026, 1.58148,
    1.58284, 1.58435, 1.58605, 1.58796, 1.59013, 1.59328, 1.59621, 1.59960,
    1.60251, 1.60824, 1.61229, 1.62032, 1.62858, 1.63886, 1.65191, 1.66888, 1.69165
  };
  G4double AbsLen_PS[nEnt] = {
    2.714*m, 3.102*m, 3.619*m, 4.343*m, 5.791*m, 7.896*m, 4.343*m, 7.896*m,
    5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m, 28.95*m, 21.71*m, 14.48*m,
    12.41*m, 8.686*m, 7.238*m, 1.200*m, 0.200*m, 0.500*m, 0.200*m, 0.100*m, 0.100*m
  };
  G4double scintFast_PS[nEnt] = {
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.07, 0.13,
    0.33, 0.63, 1.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00
  };
  mpPS = new G4MaterialPropertiesTable();
  mpPS->AddProperty("RINDEX",opEn,RI_PS,nEnt);
  mpPS->AddProperty("ABSLENGTH",opEn,AbsLen_PS,nEnt);
  mpPS->AddProperty("SCINTILLATIONCOMPONENT1",opEn,scintFast_PS,nEnt,true);
  mpPS->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  mpPS->AddConstProperty("RESOLUTIONSCALE",1.0);
  mpPS->AddConstProperty("SCINTILLATIONTIMECONSTANT1",2.8*ns,true);
  fPS->SetMaterialPropertiesTable(mpPS);
  fPS->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // Legacy proto LG bulk (same as PMMA optically); current geometry does not place an LG — material kept for reuse/tests.
  G4MaterialPropertiesTable* mpProtoLGMatchScint = new G4MaterialPropertiesTable();
  mpProtoLGMatchScint->AddProperty("RINDEX", opEn, RI_PS, nEnt);
  mpProtoLGMatchScint->AddProperty("ABSLENGTH", opEn, AbsLen_PMMA, nEnt);
  fProtoLGMatchScint->SetMaterialPropertiesTable(mpProtoLGMatchScint);

  // G4double RI_Glass[nEnt]; std::fill_n(RI_Glass, nEnt, 1.81);
  G4double RI_Glass[nEnt]; std::fill_n(RI_Glass, nEnt, 1.52);
  G4double Abslength_Glass[nEnt]; std::fill_n(Abslength_Glass, nEnt, 420.*cm);
  mpGlass = new G4MaterialPropertiesTable();
  mpGlass->AddProperty("RINDEX",opEn,RI_Glass,nEnt);
  mpGlass->AddProperty("ABSLENGTH",opEn,Abslength_Glass,nEnt);
  fGlass->SetMaterialPropertiesTable(mpGlass);

  // SiPM wafer bulk: RINDEX + short ABSLENGTH so optical photons enter Silicon and absorb in-volume.
  // (Skin-surface EFFICIENCY on dielectric_metal kills photons at the boundary without a step in Si,
  //  so SiPMSD::ProcessHits never ran.) Detection is modeled by bulk absorption, not surface QE.
  G4MaterialPropertiesTable* mpSiBulk = new G4MaterialPropertiesTable();
  G4double RI_Si[nEnt];
  std::fill_n(RI_Si, nEnt, 3.8);
  G4double Abslength_Si[nEnt];
  std::fill_n(Abslength_Si, nEnt, 50. * um);
  mpSiBulk->AddProperty("RINDEX", opEn, RI_Si, nEnt);
  mpSiBulk->AddProperty("ABSLENGTH", opEn, Abslength_Si, nEnt);
  fSi->SetMaterialPropertiesTable(mpSiBulk);

  G4double refl_SiPM[nEnt];
  std::fill_n(refl_SiPM, nEnt, 0.);

  // 표면 EFFICIENCY(QE): 나중에 복구 시 참고. 지금은 mpSiBulk 체내 흡수로 SD 스텝을 쓰므로 0.
  // 복구 방법은 다를 수 있음(경계 검출만 켜면 SD ProcessHits와 안 맞을 수 있어 Stepping/별도 HC 등 검토).
  /*
  G4double eff_SiPM[nEnt] = {
    0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10,
    0.11, 0.13, 0.15, 0.17, 0.19, 0.20, 0.22, 0.23,
    0.24, 0.25, 0.24, 0.23, 0.21, 0.20, 0.17, 0.14, 0.10
  };
  */
  G4double eff_SiPM[nEnt];
  std::fill_n(eff_SiPM, nEnt, 0.);

  mpSiPM = new G4MaterialPropertiesTable();
  mpSiPM->AddProperty("REFLECTIVITY",opEn,refl_SiPM,nEnt);
  mpSiPM->AddProperty("EFFICIENCY",opEn,eff_SiPM,nEnt);
  fSiPMSurf = new G4OpticalSurface("SiPMSurf",glisur,polished,dielectric_metal);
  fSiPMSurf->SetMaterialPropertiesTable(mpSiPM);

  G4double filterEff[nEnt] = {
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.900, 0.734,
    0.568, 0.402, 0.296, 0.070, 0.000, 0.000, 0.000, 0.000, 0.000
  };
  G4double filterRef[nEnt]; std::fill_n(filterRef,nEnt,0.);
  G4double RI_gel[nEnt]; std::fill_n(RI_gel,nEnt,1.52);
  mpFilter = new G4MaterialPropertiesTable();
  mpFilter->AddProperty("RINDEX",opEn,RI_gel,nEnt);
  fGelatin->SetMaterialPropertiesTable(mpFilter);
  mpFilterSurf = new G4MaterialPropertiesTable();
  mpFilterSurf->AddProperty("TRANSMITTANCE",opEn,filterEff,nEnt);
  mpFilterSurf->AddProperty("REFLECTIVITY",opEn,filterRef,nEnt);
  fFilterSurf = new G4OpticalSurface("FilterSurf",glisur,polished,dielectric_dielectric);
  fFilterSurf->SetMaterialPropertiesTable(mpFilterSurf);
}
