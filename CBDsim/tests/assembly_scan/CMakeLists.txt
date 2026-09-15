cmake_minimum_required(VERSION 3.16)
project(AssemblyScanValidation LANGUAGES CXX)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
add_executable(assembly_scan_test scan_test.cc ../../src/CBDsimPrimaryGeneratorAction.cc
  ../../src/CBDsimDetectorConstructionProto.cc ../../src/CBDsimMaterials.cc
  ../../src/CBDsimSiPMSD.cc ../../src/CBDsimSiPMHit.cc ../../src/CBDsimOpticalDiagnostics.cc)
target_include_directories(assembly_scan_test PRIVATE ../../include ../../../rootIO/include)
target_compile_features(assembly_scan_test PRIVATE cxx_std_17)
target_link_libraries(assembly_scan_test PRIVATE ${Geant4_LIBRARIES})
