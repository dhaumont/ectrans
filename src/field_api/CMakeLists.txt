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
  set (oneValueArgs BACKEND DESTINATION OUTPUT)
  set (multiValueArgs)

  cmake_parse_arguments(_PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(backend ${_PAR_BACKEND})
  set(destination ${_PAR_DESTINATION})
  file(MAKE_DIRECTORY ${destination}/field_view)

  list( APPEND files 
        ectrans_field_api_mod.F90
        dir_trans_field_api.F90
        inv_trans_field_api.F90)
  set(outfiles)
  foreach(file_i ${files})
    set(outfile "${destination}/${file_i}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)

set (platforms)
if (HAVE_CPU)
    list(APPEND platforms cpu)
endif()
# Note: ectrans_field_api_gpu_${prec} depends on ectrans_field_view_gpu_${prec} which
# is not built (see src/trans/field_view/CMakeLists.txt). GPU executables should link
# against the CPU variant ectrans_field_api_${prec} (which requires HAVE_CPU to be
# enabled) until a proper GPU implementation of field_view exists.
foreach(platform ${platforms})
  foreach( prec dp sp)
    if( HAVE_${prec})
      if (${platform} MATCHES "gpu")
        set (platform_tag "gpu_")
      else()
        set (platform_tag "")
      endif()
      generate_backend_sources( BACKEND ${prec} OUTPUT ectrans_field_api_${prec}_src  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_${prec})
      ecbuild_add_library(
          TARGET           ectrans_field_api_${platform_tag}${prec}
          LINKER_LANGUAGE  Fortran
          SOURCES          ${ectrans_field_api_${prec}_src}
          PUBLIC_LIBS      ectrans_field_view_${platform_tag}${prec} field_api_${prec} fiat parkind_${prec}
          PUBLIC_INCLUDES  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/field_api/include/ectrans>
                           $<INSTALL_INTERFACE:include/ectrans> )

      ecbuild_target_fortran_module_directory(
          TARGET            ectrans_field_api_${platform_tag}${prec}
          MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/ectrans_field_api_${platform_tag}${prec}
          INSTALL_DIRECTORY module/ectrans_field_api_${platform_tag}${prec} )
    endif()
  endforeach()
endforeach()

list( APPEND field_api_includes 
  include/ectrans/dir_trans_field_api.h
  include/ectrans/inv_trans_field_api.h
)
install(
  FILES        ${field_api_includes}
  DESTINATION  include/ectrans
)
