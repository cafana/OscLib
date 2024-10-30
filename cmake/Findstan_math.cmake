#[================================================================[.rst:
Findstan_math
----------
  find stan_math

#]================================================================]

#MESSAGE("stan_math_INC is $ENV{STAN_MATH_INC}")
if (stan_math_FOUND)
  set(_cet_stan_math_config_mode CONFIG_MODE)
else()
  unset(_cet_stan_math_config_mode)
  find_file(_cet_stan_math_h NAMES stan/math.hpp HINTS ENV STAN_MATH_INC)
  #MESSAGE("${_cet_stan_math_h}")
  if (_cet_stan_math_h)
    get_filename_component(_cet_stan_math_include_dir "${_cet_stan_math_h}" PATH)
    #MESSAGE("${_cet_stan_math_include_dir}")
    if (_cet_stan_math_include_dir STREQUAL "/")
      unset(_cet_stan_math_include_dir)
    endif()
  endif()
  if (EXISTS "${_cet_stan_math_include_dir}")
    set(stan_math_FOUND TRUE)
    #MESSAGE("FOUND - stan_math ${stan_math_FOUND}")
    get_filename_component(_cet_stan_math_dir "${_cet_stan_math_include_dir}" PATH)
    #MESSAGE("${_cet_stan_math_dir}")
    if (_cet_stan_math_dir STREQUAL "/")
      unset(_cet_stan_math_dir)
    endif()
    set(stan_math_INCLUDE_DIRS "$ENV{STAN_MATH_INC}")
    #MESSAGE(VERBOSE "DIR set to ${stan_math_INCLUDE_DIRS}")
  endif()
endif()

if (stan_math_FOUND
    AND NOT TARGET stan_math::stan_math)
  add_library(stan_math::stan_math INTERFACE IMPORTED)
  set_target_properties(stan_math::stan_math PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${stan_math_INCLUDE_DIRS}"
    INTERFACE_COMPILE_DEFINITIONS "${STAN_MATH_DEFINITIONS}"
  )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(stan_math ${_cet_stan_math_config_mode}
  REQUIRED_VARS stan_math_FOUND
  stan_math_INCLUDE_DIRS)

unset(_cet_stan_math_FIND_REQUIRED)
unset(_cet_stan_math_config_mode)
unset(_cet_stan_math_dir)
unset(_cet_stan_math_include_dir)
unset(_cet_stan_math_h CACHE)