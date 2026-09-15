cmake_minimum_required(VERSION 3.16)
project(PrimaryGeneratorValidation LANGUAGES CXX)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
set(PROTO_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/CBDsimDetectorConstructionProto.cc" CACHE FILEPATH "Production geometry source (override for before regression)")
add_executable(primary_generator_test primary_test.cc ../../src/CBDsimPrimaryGeneratorAction.cc ${PROTO_SOURCE}
  ../../src/CBDsimMaterials.cc ../../src/CBDsimSiPMSD.cc ../../src/CBDsimSiPMHit.cc
  ../../src/CBDsimOpticalDiagnostics.cc)
target_include_directories(primary_generator_test PRIVATE ../../include ../../../rootIO/include)
target_compile_features(primary_generator_test PRIVATE cxx_std_17)
target_link_libraries(primary_generator_test PRIVATE ${Geant4_LIBRARIES})
