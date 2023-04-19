foreach( prec sp )

  if( HAVE_${prec} )

    ecbuild_add_executable(TARGET  aatestprog_${prec}
                           SOURCES aatestprog.F90
                           LIBS
							 etrans_${prec}
                             trans_${prec}
                             fiat
                             parkind_${prec}
						   LINKER_LANGUAGE Fortran
                          )	
	target_link_libraries( aatestprog_${prec} ${OpenMP_Fortran_FLAGS} OpenMP::OpenMP_Fortran )

  endif()
	
endforeach()
						  