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
  file(MAKE_DIRECTORY ${destination}/field_api_${backend})
  file(MAKE_DIRECTORY ${destination}/field_api_${backend}/trans)


  list( APPEND files
  field_api_ectrans_module.F90
  trans/dir_trans_field_api.F90
  trans/inv_trans_field_api.F90
  )


  set(outfiles)
  foreach(file_i ${files})
    get_filename_component(outfile_name    ${file_i} NAME)
    get_filename_component(outfile_name_we ${file_i} NAME_WE)
    get_filename_component(outfile_ext     ${file_i} EXT)
    get_filename_component(outfile_dir     ${file_i} DIRECTORY)
    #set(outfile "${destination}/${file_i}")
    set(outfile "${destination}/${outfile_name_we}_${backend}${outfile_ext}")
    ecbuild_debug("Generate ${outfile}")
    generate_file(BACKEND ${backend} INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${file_i} OUTPUT ${outfile})
    list(APPEND outfiles ${outfile})
  endforeach(file_i)
  set(${_PAR_OUTPUT} ${outfiles} PARENT_SCOPE)
endfunction(generate_backend_sources)


# loop over precisions
foreach( prec dp )
  if( HAVE_${prec} )
   generate_backend_sources( BACKEND ${prec} 
                             OUTPUT ectrans_field_api_${prec}_src  
                             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/generated/ectrans_field_api_${prec})
 
    ecbuild_add_library(
       TARGET           ectrans_field_api_${prec}
       LINKER_LANGUAGE  Fortran
       SOURCES          ${ectrans_field_api_${prec}_src}
       PUBLIC_LIBS      field_api_${prec} fiat trans_${prec} parkind_${prec}
       PUBLIC_INCLUDES   $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/trans/include/ectrans>
       )
  endif()
endforeach()



