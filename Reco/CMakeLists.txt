cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(Reco)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

add_executable(Reco Reco.cc ${sources} ${headers})
target_compile_features(Reco PUBLIC cxx_std_17)
target_link_libraries(
  Reco
  rootIO
  ${CMAKE_DL_LIBS}
)

set(
  Reco_SCRIPTS
)

foreach(_script ${Reco_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

install(TARGETS Reco DESTINATION bin)
