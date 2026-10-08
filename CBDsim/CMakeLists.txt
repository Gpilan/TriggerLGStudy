#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(CBDsim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
find_package(Geant4 11.4.0 EXACT REQUIRED ui_all vis_all)
else()
find_package(Geant4 11.4.0 EXACT REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find HepMC (required package)
#
#set(HEPMC_INCLUDE_DIR "$ENV{HEPMC_INCLUDE_DIR}")
#set(HEPMC_LIBRARIES "$ENV{HEPMC_LIBRARIES}")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
#                    ${HEPMC_INCLUDE_DIR}
                   )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(CBDsim CBDsim.cc ${sources} ${headers})
option(CBDsim_USE_BOUNDARY_ADAPTER "Use the version-specific Trigger boundary adapter" ON)
target_compile_definitions(CBDsim PRIVATE
  CBDsim_USE_BOUNDARY_ADAPTER=$<BOOL:${CBDsim_USE_BOUNDARY_ADAPTER}>)
target_compile_features(CBDsim PUBLIC cxx_std_17)
target_link_libraries(
  CBDsim
  ${Geant4_LIBRARIES}
  rootIO
#  ${HEPMC_LIBRARIES}
)

#----------------------------------------------------------------------------
# Copy interactive UI/visualization assets and ROOT inspection tools.
# Batch run macros are generated per run by scripts/frozen_baseline.py.
#
set(CBDsim_SCRIPTS
  icons.mac
  gui.mac
  run.png
  init.mac
  init_vis.mac
  vis.mac
  verify_root.C
  verify_root_run.sh
)

foreach(_script ${CBDsim_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

# ROOT 통계 산출물 보관 (빌드 트리에 stats/ 유지)
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/stats")
configure_file(
  ${PROJECT_SOURCE_DIR}/stats/README.txt
  ${PROJECT_BINARY_DIR}/stats/README.txt
  COPYONLY
)

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
#add_custom_target(WGR16 DEPENDS exampleB)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CBDsim DESTINATION bin)
