module ectrans_field_api_helper

use field_module, only:field_1rb, field_2rb, field_3rb, field_4rb
use field_factory_module
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

type wrapped_fields

! Set of fields for spectral transform

  class (field_3rb), pointer :: spscalars     ! spectral scalar fields
  class (field_2rb), pointer :: spscalars2    ! spectral surfacic scalar fields
  class (field_2rb), pointer :: spvor, spdiv  ! spectral vorticity and divergence

  class (field_3rb), pointer :: vor, div      ! grid-point vorticity and divergence
  class (field_3rb), pointer :: u, v          ! grid-point u and v fields  
  class (field_3rb), pointer :: udm, vdm      ! grid-point u and derivatives

  class (field_4rb), pointer :: scalars       ! grid-point scalar fields
  class (field_4rb), pointer :: scalars_ew    ! grid-point scalar fields derivatives ew
  class (field_4rb), pointer :: scalars_ns    ! grid-point scalar fields derivatives ns

  class (field_3rb), pointer :: scalars2      ! grid-point surfacic scalar fields
  class (field_3rb), pointer :: scalars2_ew   ! grid-point surfacic scalar fields derivatives ew
  class (field_3rb), pointer :: scalars2_ns   ! grid-point surfacic scalar fields derivatives ns
end type wrapped_fields

type fields_lists

! List of field lists that will be used as parameter to inv_trans_field_api and dir_trans_field_api

  type (field_basic_ptr), allocatable:: u (:), v (:)                ! grid-point u and v fields
  type (field_basic_ptr), allocatable:: scalar (:)                  ! grid-point scalar fields
  type (field_basic_ptr), allocatable:: spvor (:), spdiv (:)        ! spectral vorticity and divergence
  type (field_basic_ptr), allocatable:: vor (:), div (:)            ! grid-point vorticity and diverence
  type (field_basic_ptr), allocatable:: spscalar (:)                ! spectral scalar fields
  type (field_basic_ptr), allocatable:: udm (:), vdm (:)            ! grid-point u and derivatives
  type (field_basic_ptr), allocatable:: scalardm (:), scalardl (:)  ! grid space scalar derivatives
  end type fields_lists
  
contains

subroutine wrap_benchmark_fields(ywflds, lvordiv, lscders, luvders,& 
                               & sp3d, spc2, zgmv, zgmvs, zgp2,    & 
                               & jbegin_uv,jend_uv,                & 
                               & jbegin_sc,jend_sc,                & 
                               & jbegin_scder_ns, jend_scder_ns,   & 
                               & jbegin_scder_ew, jend_scder_ew,   & 
                               & jbegin_uder_ew, jend_uder_ew,     & 
                               & jbegin_vder_ew, jend_vder_ew)

  ! Wrap the arrays given as input in field API objects

    type(wrapped_fields), intent(inout) :: ywflds
    logical, intent(in) :: lvordiv
    logical, intent(in) :: lscders
    logical, intent(in) :: luvders
    real(kind=jprb), intent(in) :: sp3d(:,:,:)
    real(kind=jprb), intent(in) :: spc2(:,:)
    real(kind=jprb), intent(in) :: zgmv(:,:,:,:)
    real(kind=jprb), intent(in) :: zgmvs(:,:,:)
    real(kind=jprb), intent(in) :: zgp2 (:,:,:)
    integer(kind=jpim), intent(in) :: jbegin_uv
    integer(kind=jpim), intent(in) :: jend_uv
    integer(kind=jpim), intent(in) :: jbegin_sc
    integer(kind=jpim), intent(in) :: jend_sc
    integer(kind=jpim), intent(in) :: jbegin_scder_ns
    integer(kind=jpim), intent(in) :: jend_scder_ns
    integer(kind=jpim), intent(in) :: jbegin_scder_ew
    integer(kind=jpim), intent(in) :: jend_scder_ew
    integer(kind=jpim), intent(in) :: jbegin_uder_ew
    integer(kind=jpim), intent(in) :: jend_uder_ew
    integer(kind=jpim), intent(in) :: jbegin_vder_ew
    integer(kind=jpim), intent(in) :: jend_vder_ew
  
  if (lvordiv) then
      if (jbegin_uv>0 )      call field_new(ywflds%u, data=zgmv(:,:,jbegin_uv,:))
      if (jend_uv>0 )        call field_new(ywflds%v, data=zgmv(:,:,jend_uv,:))
  
      ! In the benchmark, vorticity is not computed
      !call field_new(ywflds%vor,         data=zgmv(:,:,jbegin_uv,:))
      !call field_new(ywflds%div,         data=zgmv(:,:,jend_uv,:))
  endif
  
  ! spectral vector fields
  if (size(sp3d,3) >=1 ) call field_new(ywflds%spvor,      data=sp3d(:,:,1))
  if (size(sp3d,3) >=2 ) call field_new(ywflds%spdiv,      data=sp3d(:,:,2))
  
  ! spectral scalar fields
  if (size(sp3d,3) >=3 ) call field_new(ywflds%spscalars,  data=sp3d(:,:,3:))
  if (size(spc2,2) >=1 ) call field_new(ywflds%spscalars2, data=spc2(:,:))
  
  ! spectral surfacic scalar fields
  if (size(zgmvs,2)>=1)  call field_new(ywflds%scalars2,   data=zgmvs(:,1:1,:))

  ! grid-point vector derivatives
  if (luvders) then
    call field_new(ywflds%udm, data=zgmv(:,:,jbegin_uder_ew,:))
    call field_new(ywflds%vdm, data=zgmv(:,:,jend_uder_ew,:))
  endif

  ! grid-point scalars fields
  if (jend_sc>0 .and. jend_sc>=jbegin_sc ) call field_new(ywflds%scalars,  data=zgmv(:,:,jbegin_sc:jend_sc,:))

  ! grid-point scalars derivatives fields
  if (lscders) then    
    if (jend_scder_ew>0 .and. jend_scder_ew>=jbegin_scder_ew ) call field_new(ywflds%scalars_ew,  data=zgmv(:,:,jbegin_scder_ew:jend_scder_ew,:))
    if (jend_scder_ns>0 .and. jend_scder_ns>=jbegin_scder_ns ) call field_new(ywflds%scalars_ns,  data=zgmv(:,:,jbegin_scder_ns:jend_scder_ns,:))
  
    ! grid-point surfacic scalars derivatives fields
    if (size(zgmvs,2)>=2)     call field_new(ywflds%scalars2_ns, data=zgmvs(:,2:2,:))
    if (size(zgmvs,2)>=3)     call field_new(ywflds%scalars2_ew, data=zgmvs(:,3:3,:))      
  endif
end subroutine wrap_benchmark_fields

subroutine create_fields_lists(ywflds,ylf, nbsetlev,nbsetsc2)

  ! Create field lists in ylf from field API objects in ywflds

  type(wrapped_fields), intent(in) :: ywflds       !input fields api objects
  type(fields_lists), intent(inout), target :: ylf ! output field lists
  integer(kind=jpim), intent(in) :: nbsetlev(:)    ! 'b-set' for vector fields
  integer(kind=jpim), intent(in) :: nbsetsc2(:)    ! 'b-set' for surfacic fields
  
  if(associated(ywflds%spvor)) ylf%spvor=[b(ywflds%spvor,'sp_vor')]
  
  if(associated(ywflds%spdiv)) ylf%spdiv= [b(ywflds%spdiv,'spdiv')]
  
  if(associated(ywflds%u)) ylf%u = [b(ywflds%u,'u',nbsetlev)]
  if(associated(ywflds%v)) ylf%v = [b(ywflds%v,'v',nbsetlev)]
      
  if(associated(ywflds%udm)) ylf%udm=[b(ywflds%udm,'u_fdm', nbsetlev)]
  if(associated(ywflds%vdm)) ylf%vdm=[b(ywflds%vdm,'v_fdm', nbsetlev)]
    
  if(associated(ywflds%vor))  ylf%vor = [b(ywflds%vor,'vor', nbsetlev)]
  if(associated(ywflds%div))  ylf%div = [b(ywflds%div,'div', nbsetlev)]
  
  if (associated(ywflds%spscalars) .and. associated(ywflds%spscalars2) ) then
    ylf%spscalar = [b(ywflds%spscalars,'sp_scalars'), b(ywflds%spscalars2,'sp_scalars2')]
  else if (associated(ywflds%spscalars)) then
    ylf%spscalar = [b(ywflds%spscalars,'sp_scalars')]    
  else if (associated(ywflds%spscalars2)) then
    ylf%spscalar = [b(ywflds%spscalars2,'sp_scalars2')]  
  endif
    
  if (associated(ywflds%scalars) .and. associated(ywflds%scalars2) ) then
    ylf%scalar = [b(ywflds%scalars,'scalars', nbsetlev), b(ywflds%scalars2,'scalars2', nbsetsc2)]
  else if (associated(ywflds%scalars)) then
    ylf%scalar = [b(ywflds%scalars,'scalars', nbsetlev)]    
  else if (associated(ywflds%scalars2)) then
    ylf%scalar = [b(ywflds%scalars2,'scalars2', nbsetsc2)]  
  endif
  
  if (associated(ywflds%scalars_ns) .and. associated(ywflds%scalars2_ns) ) then
    ylf%scalardm = [b(ywflds%scalars_ns,'scalars_ns', nbsetlev), b(ywflds%scalars2_ns,'scalars2_ns', nbsetsc2)]
  else if (associated(ywflds%scalars_ns)) then
    ylf%scalardm = [b(ywflds%scalars_ns,'scalars_ns', nbsetlev)]
  else if (associated(ywflds%scalars2_ns)) then
    ylf%scalardm = [b(ywflds%scalars2_ns,'scalars2_ns', nbsetsc2)]  
  endif
  
  if (associated(ywflds%scalars_ew) .and. associated(ywflds%scalars2_ew) ) then
    ylf%scalardl = [b(ywflds%scalars_ew,'scalars_ew', nbsetlev), b(ywflds%scalars2_ew,'scalars2_ew', nbsetsc2)]
  else if (associated(ywflds%scalars_ew)) then
    ylf%scalardl = [b(ywflds%scalars_ew,'scalars_ew', nbsetlev)]    
  else if (associated(ywflds%scalars2_ew)) then
    ylf%scalardl = [b(ywflds%scalars2_ew,'scalars2_ew', nbsetsc2)]  
  endif
 end subroutine create_fields_lists

 subroutine delete_wrapped_fields(ywflds)
  
  ! Delete  all fields in ywflds

  type(wrapped_fields), intent(inout) :: ywflds
   
  if(associated(ywflds%spvor)) call field_delete(ywflds%spvor)
  if(associated(ywflds%spdiv)) call field_delete(ywflds%spdiv)
  if(associated(ywflds%spscalars)) call field_delete(ywflds%spscalars)
  if(associated(ywflds%spscalars2)) call field_delete(ywflds%spscalars2)

  if(associated(ywflds%u)) call field_delete(ywflds%u)
  if(associated(ywflds%v)) call field_delete(ywflds%v)
  if(associated(ywflds%udm)) call field_delete(ywflds%udm)
  if(associated(ywflds%vdm)) call field_delete(ywflds%vdm)
  if(associated(ywflds%scalars)) call field_delete(ywflds%scalars)
  if(associated(ywflds%scalars_ew)) call field_delete(ywflds%scalars_ew)
  if(associated(ywflds%scalars_ns)) call field_delete(ywflds%scalars_ns)
  if(associated(ywflds%vor)) call field_delete(ywflds%vor)
  if(associated(ywflds%div)) call field_delete(ywflds%div)

  if(associated(ywflds%scalars2)) call field_delete(ywflds%scalars2)
  if(associated(ywflds%scalars2_ew)) call field_delete(ywflds%scalars2_ew)
  if(associated(ywflds%scalars2_ns)) call field_delete(ywflds%scalars2_ns) 
end subroutine delete_wrapped_fields

subroutine delete_fields_lists(yfl)

  ! Delete  all field lists in yfl

  type(fields_lists), intent(inout) ::yfl
  if (allocated(yfl%u)) deallocate(yfl%u)
  if (allocated(yfl%v)) deallocate(yfl%v)
  if (allocated(yfl%scalar)) deallocate(yfl%scalar)
  if (allocated(yfl%spscalar)) deallocate(yfl%spscalar)
  if (allocated(yfl%spvor)) deallocate(yfl%spvor)
  if (allocated(yfl%spdiv)) deallocate(yfl%spdiv)
  if (allocated(yfl%vor)) deallocate(yfl%vor)
  if (allocated(yfl%div)) deallocate(yfl%div)
  if (allocated(yfl%udm)) deallocate(yfl%udm)
  if (allocated(yfl%vdm)) deallocate(yfl%vdm) 
  if (allocated(yfl%scalardm)) deallocate(yfl%scalardm)
  if (allocated(yfl%scalardl)) deallocate(yfl%scalardl)
end subroutine delete_fields_lists

subroutine nullify_wrapped_fields(ywflds)

    ! Nullify all pointers in ywflds

  type(wrapped_fields), intent(inout) :: ywflds
    
  nullify(ywflds%spvor)
  nullify(ywflds%spdiv)
  nullify(ywflds%spscalars)
  nullify(ywflds%spscalars2)

  nullify(ywflds%u)
  nullify(ywflds%v)
  nullify(ywflds%udm)
  nullify(ywflds%vdm)
  nullify(ywflds%scalars)
  nullify(ywflds%scalars_ew)
  nullify(ywflds%scalars_ns)
  nullify(ywflds%vor)
  nullify(ywflds%div)

  nullify(ywflds%scalars2)
  nullify(ywflds%scalars2_ew)
  nullify(ywflds%scalars2_ns)
end subroutine nullify_wrapped_fields


subroutine output_wrapped_fields(nout, ywflds)
  
  ! output the adress of all data fields in ywflds

  type(wrapped_fields), intent(in) :: ywflds
  integer(kind=jpim), intent(in) :: nout

  write(nout,*) "ywflds%spvor", loc(ywflds%spvor)
  write(nout,*) "ywflds%spdiv", loc(ywflds%spdiv)
  write(nout,*) "ywflds%spscalars", loc(ywflds%spscalars)
  write(nout,*) "ywflds%spscalars2", loc(ywflds%spscalars2)

  write(nout,*) "ywflds%u", loc(ywflds%u)
  write(nout,*) "ywflds%v", loc(ywflds%v)
  write(nout,*) "ywflds%udm", loc(ywflds%udm)
  write(nout,*) "ywflds%vdm", loc(ywflds%vdm)
  write(nout,*) "ywflds%scalars", loc(ywflds%scalars)
  write(nout,*) "ywflds%scalars_ew", loc(ywflds%scalars_ew)
  write(nout,*) "ywflds%scalars_ns", loc(ywflds%scalars_ns)
  write(nout,*) "ywflds%vor", loc(ywflds%vor)
  write(nout,*) "ywflds%div", loc(ywflds%div)

  write(nout,*) "ywflds%scalars2", loc(ywflds%scalars2)
  write(nout,*) "ywflds%scalars2_ew", loc(ywflds%scalars2_ew)
  write(nout,*) "ywflds%scalars2_ns", loc(ywflds%scalars2_ns)
end subroutine output_wrapped_fields



subroutine output_fields_lists(nout,yfl)

  ! output the size of all field lists in yfl

  integer(kind=jpim), intent(in) :: nout
  type(fields_lists), intent(in) :: yfl
  
  if (allocated(yfl%u)) write(nout,*) "yfl%u", size(yfl%u)  
  if (allocated(yfl%v)) write(nout,*) "yfl%v", size(yfl%v)
  if (allocated(yfl%scalar)) write(nout,*) "yfl%scalar", size(yfl%scalar)
  if (allocated(yfl%spscalar)) write(nout,*) "yfl%spscalar", size(yfl%spscalar)
  if (allocated(yfl%spvor)) write(nout,*) "yfl%spvor", size(yfl%spvor)
  if (allocated(yfl%spdiv)) write(nout,*) "yfl%spdiv", size(yfl%spdiv)
  if (allocated(yfl%vor)) write(nout,*) "yfl%vor", size(yfl%vor)
  if (allocated(yfl%div)) write(nout,*) "yfl%div", size(yfl%div)
  if (allocated(yfl%udm)) write(nout,*) "yfl%udm", size(yfl%udm)
  if (allocated(yfl%vdm)) write(nout,*) "yfl%vdm", size(yfl%vdm)
  if (allocated(yfl%scalardm)) write(nout,*) "yfl%scalardm", size(yfl%scalardm)
  if (allocated(yfl%scalardl)) write(nout,*) "yfl%scalardl", size(yfl%scalardl)
end subroutine output_fields_lists


end module ectrans_field_api_helper