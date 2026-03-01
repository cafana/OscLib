cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

project(OscLib VERSION 1.2.0)

find_package(cetmodules REQUIRED)
include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

# default to 'release' build if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

option(OscLib_USE_STAN "Build with Stan dependency?" OFF)
option(OscLib_USE_STAN_THREADS "Use threading with Stan?" OFF)

find_package(ROOT REQUIRED EXPORT)
find_package(GSL REQUIRED EXPORT)

cet_make_library(LIBRARY_NAME OscLib_options INTERFACE NO_SOURCE)

if(OscLib_USE_STAN)
  target_compile_definitions(OscLib_options INTERFACE OscLib_STAN _REENTRANT)
  if(OscLib_USE_STAN_THREADS)
    target_compile_definitions(OscLib_options INTERFACE STAN_THREADS)
  endif()
  find_package(stan_math REQUIRED MODULE EXPORT)
  target_link_libraries(OscLib_options INTERFACE stan_math::stan_math)
endif()

find_package(Eigen3 REQUIRED EXPORT)

cmake_policy(SET CMP0167 NEW) #don't use builtin FindBoost.cmake
find_package(Boost REQUIRED EXPORT)

target_link_libraries(OscLib_options INTERFACE Eigen3::Eigen GSL::gsl ROOT::Core ROOT::Hist)

add_subdirectory(OscLib)

cet_cmake_config()
