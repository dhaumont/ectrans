module ectrans_field_api_helper

use field_module, only:field_1rb, field_2rb, field_3rb, field_4rb
use field_factory_module
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

type wrapped_fields
  class (field_3rb), pointer :: spscalars
  class (field_2rb), pointer :: spscalars2
  class (field_2rb), pointer :: spvor, spdiv

  class (field_3rb), pointer :: vor, div
  class (field_3rb), pointer :: u, v
  class (field_3rb), pointer :: udm, vdm

  class (field_4rb), pointer :: scalars
  class (field_4rb), pointer :: scalars_ew
  class (field_4rb), pointer :: scalars_ns

  class (field_3rb), pointer :: scalars2
  class (field_3rb), pointer :: scalars2_ew
  class (field_3rb), pointer :: scalars2_ns
end type wrapped_fields

type fields_lists
  type (field_basic_ptr), allocatable:: u (:), v (:)
  type (field_basic_ptr), allocatable:: scalar (:)
  type (field_basic_ptr), allocatable:: spvor (:), spdiv (:)
  type (field_basic_ptr), allocatable:: vor (:), div (:)
  type (field_basic_ptr), allocatable:: spscalar (:)
  type (field_basic_ptr), allocatable:: udm (:), vdm (:)
  type (field_basic_ptr), allocatable:: scalardm (:), scalardl (:)  
  end type fields_lists
  
contains

subroutine output_wrapped_fields(nout, ywflds)
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

subroutine nullify_wrapped_fields(ywflds)
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

subroutine wrap_fields(ywflds, lvordiv, lscders, luvders,& 
                    sp3d, spc2, zgmv, zgmvs, zgp2,& 
                    &jbegin_uv,jend_uv,& 
                    &jbegin_sc,jend_sc,& 
                    &jbegin_scder_ns, jend_scder_ns,& 
                    &jbegin_scder_ew, jend_scder_ew,& 
                    &jbegin_uder_ew, jend_uder_ew,& 
                    &jbegin_vder_ew, jend_vder_ew)

  type(wrapped_fields), intent(inout) :: ywflds
  logical :: lvordiv
  logical :: lscders
  logical :: luvders
  real(kind=jprb), intent(in) :: sp3d(:,:,:)
  real(kind=jprb), intent(in) :: spc2(:,:)
  real(kind=jprb), intent(in) :: zgmv(:,:,:,:)
  real(kind=jprb), intent(in) :: zgmvs(:,:,:)
  real(kind=jprb), intent(in) :: zgp2 (:,:,:)
  integer(kind=jpim) :: jbegin_uv
  integer(kind=jpim) :: jend_uv
  integer(kind=jpim) :: jbegin_sc
  integer(kind=jpim) :: jend_sc
  integer(kind=jpim) :: jbegin_scder_ns
  integer(kind=jpim) :: jend_scder_ns
  integer(kind=jpim) :: jbegin_scder_ew
  integer(kind=jpim) :: jend_scder_ew
  integer(kind=jpim) :: jbegin_uder_ew
  integer(kind=jpim) :: jend_uder_ew
  integer(kind=jpim) :: jbegin_vder_ew
  integer(kind=jpim) :: jend_vder_ew
    
  call delete_wrapped_fields(ywflds)
  call nullify_wrapped_fields(ywflds)

  
  if (lvordiv) then
    if (jbegin_uv>0 )      call field_new(ywflds%u,         data=zgmv(:,:,jbegin_uv,:))
    if (jend_uv>0 )        call field_new(ywflds%v,         data=zgmv(:,:,jend_uv,:))
  !if (jbegin_uv>0 )                    call field_new(ywflds%vor,         data=zgmv(:,:,jbegin_uv,:))
  !if (jend_uv>0 )                      call field_new(ywflds%div,         data=zgmv(:,:,jend_uv,:))

  endif
  
  if (size(sp3d,3) >=1 )   call field_new(ywflds%spvor,        data=sp3d(:,:,1))
  if (size(sp3d,3) >=2 )   call field_new(ywflds%spdiv,        data=sp3d(:,:,2))
  if (size(sp3d,3) >=3 )   call field_new(ywflds%spscalars,    data=sp3d(:,:,3:))
  if (size(spc2,2) >=1 ) call field_new(ywflds%spscalars2,  data=spc2(:,:))
  
  if (size(zgmvs,2)>=1)    call field_new(ywflds%scalars2,    data=zgmvs(:,1:1,:))

  if (luvders) then
    call field_new(ywflds%udm,         data=zgmv(:,:,jbegin_uder_ew,:))
    call field_new(ywflds%vdm,         data=zgmv(:,:,jend_uder_ew,:))
  endif

  if (jend_sc>0 .and. jend_sc>=jbegin_sc )  call field_new(ywflds%scalars,     data=zgmv(:,:,jbegin_sc:jend_sc,:))
  if (lscders) then    
    if (jend_scder_ew>0 .and. jend_scder_ew>=jbegin_scder_ew ) call field_new(ywflds%scalars_ew,  data=zgmv(:,:,jbegin_scder_ew:jend_scder_ew,:))
    if (jend_scder_ns>0 .and. jend_scder_ns>=jbegin_scder_ns ) call field_new(ywflds%scalars_ns,  data=zgmv(:,:,jbegin_scder_ns:jend_scder_ns,:))
  
    if (size(zgmvs,2)>=2)     call field_new(ywflds%scalars2_ns, data=zgmvs(:,2:2,:))
    if (size(zgmvs,2)>=3)     call field_new(ywflds%scalars2_ew, data=zgmvs(:,3:3,:))      
  endif
end subroutine wrap_fields

subroutine delete_wrapped_fields(ywflds)
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

subroutine output_fields_lists(nout,yfl)
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

subroutine create_fields_lists(ywflds,ylf, nbsetlev,nbsetsc)
  type(wrapped_fields), intent(in) :: ywflds
  type(fields_lists), intent(inout), target :: ylf
  integer(kind=jpim) :: nbsetlev(:)
  integer(kind=jpim) :: nbsetsc(:)
    
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
    ylf%scalar = [b(ywflds%scalars,'scalars', nbsetlev), b(ywflds%scalars2,'scalars2', nbsetsc)]
  else if (associated(ywflds%scalars)) then
    ylf%scalar = [b(ywflds%scalars,'scalars', nbsetlev)]    
  else if (associated(ywflds%scalars2)) then
    ylf%scalar = [b(ywflds%scalars2,'scalars2', nbsetsc)]  
  endif
  
  if (associated(ywflds%scalars_ns) .and. associated(ywflds%scalars2_ns) ) then
    ylf%scalardm = [b(ywflds%scalars_ns,'scalars_ns', nbsetlev), b(ywflds%scalars2_ns,'scalars2_ns', nbsetsc)]
  else if (associated(ywflds%scalars_ns)) then
    ylf%scalardm = [b(ywflds%scalars_ns,'scalars_ns', nbsetlev)]
  else if (associated(ywflds%scalars2_ns)) then
    ylf%scalardm = [b(ywflds%scalars2_ns,'scalars2_ns', nbsetsc)]  
  endif
  
  if (associated(ywflds%scalars_ew) .and. associated(ywflds%scalars2_ew) ) then
    ylf%scalardl = [b(ywflds%scalars_ew,'scalars_ew', nbsetlev), b(ywflds%scalars2_ew,'scalars2_ew', nbsetsc)]
  else if (associated(ywflds%scalars_ew)) then
    ylf%scalardl = [b(ywflds%scalars_ew,'scalars_ew', nbsetlev)]    
  else if (associated(ywflds%scalars2_ew)) then
    ylf%scalardl = [b(ywflds%scalars2_ew,'scalars2_ew', nbsetsc)]  
  endif

 end subroutine create_fields_lists

end module ectrans_field_api_helper