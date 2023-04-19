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
                        )

ecbuild_list_add_pattern( LIST etrans_hdr
                          GLOB
                                interface/*.h
                          QUIET
                        )
						
foreach ( prec IN LISTS precs )

	ecbuild_add_library( TARGET etrans_${prec}

	  PUBLIC_INCLUDES
				$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interface>
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


	if( ectrans_HAVE_FFTW )
	  target_link_libraries( etrans_${prec} PRIVATE ${FFTW_LIBRARIES} )
	  target_include_directories( etrans_${prec} PRIVATE ${FFTW_INCLUDE_DIRS} )
	  target_compile_definitions( etrans_${prec} PRIVATE WITH_FFTW )
	endif()

    if( prec STREQUAL sp )
        target_compile_definitions( etrans_${prec} PUBLIC TRANS_SINGLE PARKINDTRANS_SINGLE )
    endif()

    set_target_properties( etrans_${prec} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec} )
    target_include_directories( etrans_${prec} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/module/etrans_${prec}> )

    target_include_directories(  etrans_${prec} PUBLIC $<INSTALL_INTERFACE:module/etrans_${prec}> )
    install( DIRECTORY ${CMAKE_BINARY_DIR}/module/etrans_${prec}/
	  	     DESTINATION module/etrans_${prec}
		     COMPONENT modules )
	
    ecbuild_add_executable(TARGET  aatestprog_${prec}
                           SOURCES programs/aatestprog.F90
                           LIBS
							 etrans_${prec}
                             trans_${prec}
                             fiat
                             parkind_${prec}
						   LINKER_LANGUAGE Fortran
                          )	
	
	target_link_libraries( aatestprog_${prec} ${OpenMP_Fortran_FLAGS} OpenMP::OpenMP_Fortran )
	
endforeach()
