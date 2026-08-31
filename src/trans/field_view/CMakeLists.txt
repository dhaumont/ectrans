# (C) Copyright 2026- ECMWF.
# (C) Copyright 2026- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


function(generate_backend_sources)
  set (options)
  set (oneValueArgs BACKEND DESTINATION PLATFORM OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  set(platform ${_PAR_PLATFORM})

  file(MAKE_DIRECTORY ${destination}/field_view)

  # The routines dir_trans_field_view, inv_trans_field_view and ectrans_field_view_internal_util_mod
  # are compatible with all plaftorms. They are precision dependent, because they make function calls
  # to dir_trans_view_ctl_mod and inv_trans_view_ctl_mod.
  # Platform specific code only implemented for CPU version for the time beeing.

  if (${platform} STREQUAL "cpu")
      list( APPEND files
        ectrans_field_view_internal_util_mod.F90
        dir_trans_field_view.F90
        inv_trans_field_view.F90
        ${platform}/internal/dir_trans_view_ctl_mod.F90
        ${platform}/internal/inv_trans_view_ctl_mod.F90
        ${platform}/internal/ftdir_view_ctl_mod.F90
        ${platform}/internal/ltdir_view_mod.F90
        ${platform}/internal/ltinv_view_mod.F90
        ${platform}/internal/trgl_view_mod.F90
        ${platform}/internal/trltog_view_mod.F90
        ${platform}/internal/updsp_view_mod.F90
        ${platform}/internal/fsc_view_mod.F90
        ${platform}/internal/ftinv_view_ctl_mod.F90
        ${platform}/internal/ltdir_view_ctl_mod.F90
        ${platform}/internal/ltinv_view_ctl_mod.F90
        ${platform}/internal/prfi1b_view_mod.F90
        ${platform}/internal/trgtol_view_mod.F90
        ${platform}/internal/updspb_view_mod.F90)
  endif()

  set(outfiles)

  foreach(file_i ${files})
    set(outfile "${destination}/${file_i}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)

  # For GPU, a fallaback based on inv_trans and dir_trans is used to implement
  # inv_trans_field_view and dir_trans_field_view.
  # As this fallback is precision agnostic, ectrans_field_view_internal_util_mod
  # is kept precision agnostic as well.
  if (${platform} STREQUAL "gpu")
      list( APPEND outfiles
        ectrans_field_view_internal_util_mod.F90
        gpu/dir_trans_field_view.F90
        gpu/inv_trans_field_view.F90)
  endif()

  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)

set (platforms)
if (HAVE_CPU)
    list(APPEND platforms cpu)
endif()
if (HAVE_GPU)
    list(APPEND platforms gpu)
endif()
foreach(platform ${platforms})
  foreach( prec dp sp)
    if( HAVE_${prec})
      if (${platform} MATCHES "gpu")
        set (platform_tag "gpu_")
      else()
        set (platform_tag "")
      endif()
      generate_backend_sources( BACKEND ${prec}
        PLATFORM ${platform}
        OUTPUT ectrans_field_view_${platform_tag}${prec}_src
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_${platform_tag}${prec})
      ecbuild_add_library(
          TARGET           ectrans_field_view_${platform_tag}${prec}
          LINKER_LANGUAGE  Fortran
          SOURCES          ${ectrans_field_view_${platform_tag}${prec}_src}
                           ectrans_field_view_mod.F90
          PUBLIC_LIBS      trans_${platform_tag}${prec} fiat parkind_${prec}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/field_view/include/ectrans>
                           $<INSTALL_INTERFACE:include/ectrans> )

      ecbuild_target_fortran_module_directory(
          TARGET            ectrans_field_view_${platform_tag}${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans_field_view_${platform_tag}${prec}
          INSTALL_DIRECTORY module/ectrans_field_view_${platform_tag}${prec} )
    endif()
  endforeach()
endforeach()

list( APPEND field_view_includes
  include/ectrans/dir_trans_field_view.h
  include/ectrans/inv_trans_field_view.h
)
install(
  FILES        ${field_view_includes}
  DESTINATION  include/ectrans
)
