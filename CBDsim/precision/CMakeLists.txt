# Shared implementation comes from the parent CBDsim target's source list.
# Only this executable compiles the precision observer.
add_executable(CBDsimPrecision EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../CBDsim.cc ${sources} ${headers})
target_compile_features(CBDsimPrecision PRIVATE cxx_std_17)
target_compile_definitions(CBDsimPrecision PRIVATE CBDsim_PRECISION_TRACE=1)
target_include_directories(CBDsimPrecision PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(CBDsimPrecision ${Geant4_LIBRARIES} rootIO)
