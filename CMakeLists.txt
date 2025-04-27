cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

find_package(cetmodules REQUIRED)
project(osclib VERSION 0.27)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR 
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

option(STAN "Build with Stan dependency")

find_package(ROOT REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
find_package(GSL REQUIRED)

if(STAN)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOSCLIB_STAN -D_REENTRANT")
  if(TBB_ONEAPI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTBB_INTERFACE_NEW")
  endif()
  if(NOT "$ENV{USE_STAN_THREADS}" STREQUAL "")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSTAN_THREADS")
  endif
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  find_package(stan_math REQUIRED)
  find_package(TBB REQUIRED)
endif()

add_subdirectory(OscLib)

cet_cmake_config()
