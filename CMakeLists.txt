cmake_minimum_required (VERSION 3.20 FATAL_ERROR)

project(OscLib VERSION 01.00)

find_package(cetmodules REQUIRED)

include(CetCMakeEnv)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -DOSCLIB_STAN
)

cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
find_package(ROOT REQUIRED)
# find_package(Stan)
# find_package()

add_subdirectory(OscLib)

cet_cmake_config()