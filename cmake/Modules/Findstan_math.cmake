####  Findstan_math.cmake:
####     Find the stan_math library.
####     Assumes you have a $STAN_MATH_INC environment variable set.
####
####   Original author:  J. Wolcott <jwolcott@fnal.gov>
####   Date:             August 2020

set(STAN_MATH_DIR $ENV{STAN_MATH_DIR})
message(STATUS "Trying stan_math directory: ${STAN_MATH_DIR}")

set(STAN_MATH_VERSION $ENV{STAN_MATH_VERSION})
if(NOT STAN_MATH_VERSION)
	if (EXISTS ${STAN_MATH_DIR}/ups/stan_math.table)
		file(READ ${STAN_MATH_DIR}/ups/stan_math.table table_contents)
		if(${table_contents} MATCHES "VERSION=([^ \t\r\n]+)")
			set(STAN_MATH_VERSION ${CMAKE_MATCH_1})
		endif()
	endif()
endif()
message(STATUS "Found stan_math version: ${STAN_MATH_VERSION}")

set(STAN_MATH_INC $ENV{STAN_MATH_INC})
if(NOT STAN_MATH_INC)
	if (EXISTS ${STAN_MATH_DIR}/include)
		set(STAN_MATH_INC ${STAN_MATH_DIR}/include)
	endif()
endif()
message(STATUS "Found stan_math include dir: ${STAN_MATH_INC}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(stan_math
		REQUIRED_VARS STAN_MATH_DIR STAN_MATH_VERSION STAN_MATH_INC
		VERSION_VAR STAN_MATH_VERSION)

cmake_policy(SET CMP0167 NEW) #don't use builtin FindBoost.cmake
if (SPACK)
  find_package(Eigen3 REQUIRED EXPORT)
  find_package(Boost REQUIRED EXPORT)
  find_package(SUNDIALS REQUIRED EXPORT)
  find_package(TBB REQUIRED EXPORT)
else()
  find_package(Eigen3 REQUIRED EXPORT)
  find_package(Boost REQUIRED MODULE EXPORT)
  find_package(Sundials REQUIRED MODULE EXPORT)
  find_package(TBB REQUIRED MODULE EXPORT)
endif()

add_library(stan_math::include INTERFACE IMPORTED)

if(DEFINED ENV{STAN_MATH_LOCAL})
  target_include_directories(stan_math::include INTERFACE
    $<BUILD_INTERFACE:${STAN_MATH_INC}>)
else()
  target_include_directories(stan_math::include INTERFACE ${STAN_MATH_INC})
endif()

target_link_libraries(stan_math::include INTERFACE Boost::boost Eigen3::Eigen)

add_library(stan_math::stan_math INTERFACE IMPORTED)

target_link_libraries(stan_math::stan_math INTERFACE stan_math::include TBB::tbb)

if(SPACK)
  target_link_libraries(stan_math::stan_math INTERFACE SUNDIALS::kinsol_shared)
else()
  target_link_libraries(stan_math::stan_math INTERFACE Sundials::all)
endif()

