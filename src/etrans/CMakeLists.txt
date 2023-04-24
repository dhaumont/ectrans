if( HAVE_CPU)
  add_subdirectory( cpu )
endif()

if( HAVE_GPU )
  add_subdirectory( gpu )
endif()

## Install trans interface

file( GLOB etrans_interface include/etrans/* )
install(
  FILES        ${etrans_interface}
  DESTINATION  include/etrans
)

#add_subdirectory( programs )
