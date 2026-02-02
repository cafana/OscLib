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

option(OscLib_USE_STAN "Build with Stan dependency?" OFF)
option(OscLib_USE_STAN_THREADS "Use threading with Stan?" OFF)

find_package(ROOT REQUIRED)
find_package(Eigen3 REQUIRED)

cmake_policy(SET CMP0167 NEW) #don't use builtin FindBoost.cmake
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

if(OscLib_USE_STAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOscLib_STAN -D_REENTRANT")
  if(TBB_ONEAPI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTBB_INTERFACE_NEW")
  endif()
  if(OscLib_USE_STAN_THREADS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTAN_THREADS")
  endif()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(stan_math REQUIRED)
  find_package(SUNDIALS REQUIRED)
  find_package(TBB REQUIRED)
endif()

add_subdirectory(OscLib)

cet_cmake_config()
