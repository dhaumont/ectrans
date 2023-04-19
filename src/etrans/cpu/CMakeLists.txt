# Create library etrans_${prec}.so

if( HAVE_SINGLE_PRECISION )
	list( APPEND precs sp )
endif()

if( HAVE_DOUBLE_PRECISION )
	list( APPEND precs dp )
endif()

ecbuild_list_add_pattern( LIST etrans_src
                          GLOB
                                external/*.F90
								module/*.F90
								aux/*.F90
                          QUIET
                        )#
				
foreach ( prec IN LISTS precs )

	ecbuild_add_library( TARGET etrans_${prec}

	  PUBLIC_INCLUDES
				$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include/etrans>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
				$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}/trans_${prec}>
				$<INSTALL_INTERFACE:include>
				$<INSTALL_INTERFACE:include/ectrans>

	  PUBLIC_LIBS  
				fiat
				parkind_${prec}

	  SOURCES
				${etrans_src}

      INSTALL_HEADERS_LIST
				${etrans_hdr}
				
      HEADER_DESTINATION include/etrans
	  #TYPE             STATIC
    )


	if( HAVE_FFTW )
	  target_link_libraries( etrans_${prec} PRIVATE ${FFTW_LIBRARIES} )
	  target_include_directories( etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
	  target_compile_definitions( etrans_${prec} PRIVATE WITH_FFTW )
	endif()

    ectrans_target_fortran_module_directory(
      TARGET            etrans_${prec}
      MODULE_DIRECTORY  ${CMAKE_BINARY_DIR}/module/etrans_${prec}
      INSTALL_DIRECTORY module/etrans_${prec}
    )

    if( HAVE_OMP )
      target_link_libraries( etrans_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    endif()
	
endforeach()
