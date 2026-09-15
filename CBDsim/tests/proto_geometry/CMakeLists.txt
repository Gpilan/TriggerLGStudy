cmake_minimum_required(VERSION 3.16)
project(ProtoGeometryValidation LANGUAGES CXX)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
set(PROTO_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../src/CBDsimDetectorConstructionProto.cc" CACHE FILEPATH "Production geometry source (override for before regression)")
add_executable(proto_geometry_test geometry_test.cc ${PROTO_SOURCE}
  ../../src/CBDsimMaterials.cc ../../src/CBDsimSiPMSD.cc ../../src/CBDsimSiPMHit.cc
  ../../src/CBDsimOpticalDiagnostics.cc)
target_include_directories(proto_geometry_test PRIVATE ../../include ../../../rootIO/include)
target_compile_features(proto_geometry_test PRIVATE cxx_std_17)
target_link_libraries(proto_geometry_test PRIVATE ${Geant4_LIBRARIES})
