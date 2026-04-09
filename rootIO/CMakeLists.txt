cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(rootIO)

find_package(ROOT REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(DR_DIR "$ENV{DR_DIR}")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
file(GLOB sources
  ${PROJECT_SOURCE_DIR}/src/*.cc
)
file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*Interface.h
)

include(${ROOT_USE_FILE})
ROOT_GENERATE_DICTIONARY(G__rootIO ${headers} LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/LinkDef.h)
add_library(rootIO SHARED ${sources} G__rootIO.cxx)
target_include_directories(rootIO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})
target_compile_features(rootIO PUBLIC cxx_std_17)
target_link_libraries(rootIO ${ROOT_LIBRARIES})

install(TARGETS rootIO LIBRARY DESTINATION lib)
