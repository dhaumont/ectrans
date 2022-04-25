# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

## Assemble sources

ecbuild_list_add_pattern( LIST trans_src
                          GLOB
                                algor/*.F90
                                internal/*.F90
                                external/*.F90
                                internal/sharedmem.c
                          QUIET
                        )

ecbuild_list_exclude_pattern( LIST trans_src REGEX dilatation_mod.F90 )


set( HAVE_dp ${HAVE_DOUBLE_PRECISION} )
set( HAVE_sp ${HAVE_SINGLE_PRECISION} )

set( IFS_ACC_Fortran_LIBRARIES OpenACC::OpenACC_Fortran PARENT_SCOPE)

foreach( prec sp dp )
  if( HAVE_${prec} )

    ecbuild_add_library( 
      TARGET           trans_gpu_${prec}
      LINKER_LANGUAGE  Fortran
      SOURCES          ${trans_src}
                       algor/external/fourier/destroy_plan_fftr.cpp
                       algor/external/fourier/create_plan_fftr.cpp
                       algor/external/fourier/execute_plan_fftr.cpp
      TYPE             STATIC
      PUBLIC_INCLUDES  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ectrans>
                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/algor/interface>
                       $<INSTALL_INTERFACE:include/ectrans>
                       $<INSTALL_INTERFACE:include>
      PRIVATE_INCLUDES  ${MPI_Fortran_INCLUDE_PATH}
      PUBLIC_LIBS      parkind_${prec}
                       fiat
      PRIVATE_LIBS     ${LAPACK_LIBRARIES}
    )

    ectrans_target_fortran_module_directory(
      TARGET            trans_gpu_${prec}
      MODULE_DIRECTORY  ${PROJECT_BINARY_DIR}/module/trans_gpu_${prec}
      INSTALL_DIRECTORY module/trans_gpu_${prec}
    )

    #if( HAVE_OMP )
    #  target_link_libraries( trans_gpu_${prec} PRIVATE OpenMP::OpenMP_Fortran )
    #endif()

    target_compile_options( trans_gpu_${prec} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp > )    # add specific compile options here !!!

    #target_link_libraries( trans_gpu_${prec} PRIVATE OpenACC::OpenACC_Fortran )

    if( prec STREQUAL sp )
		target_compile_definitions( trans_gpu_${prec} PUBLIC TRANS_SINGLE PARKINDTRANS_SINGLE )
    endif()
	
    target_compile_definitions( trans_gpu_${prec} PUBLIC USE_CUDA_AWARE_MPI_FT )

  endif()
endforeach()

## precision-independent GPU library with CUDA/roc kernels
ecbuild_add_library( TARGET gpu
  TYPE STATIC
  SOURCES
            algor/module/rocblasSgemmBatched.cpp
            algor/module/rocblasDgemmBatched.cpp
            #algor/module/IPC_Alltoall.cu
            PRIVATE_INCLUDES  ${MPI_C_INCLUDE_PATH}
)

## Install trans interface

file( GLOB trans_interface interface/* )
install(
  FILES        ${trans_interface}
  DESTINATION  include/ectrans
)
