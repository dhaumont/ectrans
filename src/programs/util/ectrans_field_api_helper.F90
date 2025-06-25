module ectrans_field_api_helper

use field_module, only:field_1rb, field_2rb, field_3rb, field_4rb
use field_factory_module
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

type wrapped_fields
  class (field_3rb), pointer :: f_spscalars
  class (field_2rb), pointer :: f_spscalars2
  class (field_2rb), pointer :: f_spvor, f_spdiv

  class (field_3rb), pointer :: f_vor, f_div
  class (field_3rb), pointer :: f_u, f_v
  class (field_3rb), pointer :: f_udm, f_vdm

  class (field_4rb), pointer :: f_scalars
  class (field_4rb), pointer :: f_scalars_ew
  class (field_4rb), pointer :: f_scalars_ns

  class (field_3rb), pointer :: f_scalars2
  class (field_3rb), pointer :: f_scalars2_ew
  class (field_3rb), pointer :: f_scalars2_ns
end type wrapped_fields

type fields_lists
  type (field_basic_ptr), allocatable:: alloc_u (:), alloc_v (:)
  type (field_basic_ptr), allocatable:: alloc_scalar (:)
  type (field_basic_ptr), allocatable:: alloc_spvor (:), alloc_spdiv (:)
  type (field_basic_ptr), allocatable:: alloc_vor (:), alloc_div (:)
  type (field_basic_ptr), allocatable:: alloc_spscalar (:)
  type (field_basic_ptr), allocatable:: alloc_udm (:), alloc_vdm (:)
  type (field_basic_ptr), allocatable:: alloc_scalardm (:), alloc_scalardl (:)  

  type (field_basic_ptr), pointer:: u (:), v (:)
  type (field_basic_ptr), pointer:: scalar (:)
  type (field_basic_ptr), pointer:: spvor (:), spdiv (:)
  type (field_basic_ptr), pointer:: vor (:), div (:)
  type (field_basic_ptr), pointer:: spscalar (:)
  type (field_basic_ptr), pointer:: udm (:), vdm (:)
  type (field_basic_ptr), pointer:: scalardm (:), scalardl (:)  
  end type fields_lists
  
contains

subroutine output_wrapped_fields(nout, wflds)
  type(wrapped_fields), intent(in) :: wflds
  integer(kind=jpim), intent(in) :: nout

  write(nout,*) "wflds%f_spvor", loc(wflds%f_spvor)
  write(nout,*) "wflds%f_spdiv", loc(wflds%f_spdiv)
  write(nout,*) "wflds%f_spscalars", loc(wflds%f_spscalars)
  write(nout,*) "wflds%f_spscalars2", loc(wflds%f_spscalars2)

  write(nout,*) "wflds%f_u", loc(wflds%f_u)
  write(nout,*) "wflds%f_v", loc(wflds%f_v)
  write(nout,*) "wflds%f_udm", loc(wflds%f_udm)
  write(nout,*) "wflds%f_vdm", loc(wflds%f_vdm)
  write(nout,*) "wflds%f_scalars", loc(wflds%f_scalars)
  write(nout,*) "wflds%f_scalars_ew", loc(wflds%f_scalars_ew)
  write(nout,*) "wflds%f_scalars_ns", loc(wflds%f_scalars_ns)
  write(nout,*) "wflds%f_vor", loc(wflds%f_vor)
  write(nout,*) "wflds%f_div", loc(wflds%f_div)

  write(nout,*) "wflds%f_scalars2", loc(wflds%f_scalars2)
  write(nout,*) "wflds%f_scalars2_ew", loc(wflds%f_scalars2_ew)
  write(nout,*) "wflds%f_scalars2_ns", loc(wflds%f_scalars2_ns)
end subroutine output_wrapped_fields

subroutine nullify_wrapped_fields(flds)
  type(wrapped_fields), intent(inout) :: flds
    
  nullify(flds%f_spvor)
  nullify(flds%f_spdiv)
  nullify(flds%f_spscalars)
  nullify(flds%f_spscalars2)

  nullify(flds%f_u)
  nullify(flds%f_v)
  nullify(flds%f_udm)
  nullify(flds%f_vdm)
  nullify(flds%f_scalars)
  nullify(flds%f_scalars_ew)
  nullify(flds%f_scalars_ns)
  nullify(flds%f_vor)
  nullify(flds%f_div)

  nullify(flds%f_scalars2)
  nullify(flds%f_scalars2_ew)
  nullify(flds%f_scalars2_ns)
end subroutine nullify_wrapped_fields

subroutine wrap_fields(flds, lvordiv, lscders, luvders,& 
                    sp3d, spc2, zgmv, zgmvs, zgp2,& 
                    &jbegin_uv,jend_uv,& 
                    &jbegin_sc,jend_sc,& 
                    &jbegin_scder_ns, jend_scder_ns,& 
                    &jbegin_scder_ew, jend_scder_ew,& 
                    &jbegin_uder_ew, jend_uder_ew,& 
                    &jbegin_vder_ew, jend_vder_ew)

  type(wrapped_fields), intent(inout) :: flds
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
    
  call nullify_wrapped_fields(flds)

  
  if (lvordiv) then
    if (jbegin_uv>0 )      call field_new(flds%f_u,         data=zgmv(:,:,jbegin_uv,:))
    if (jend_uv>0 )        call field_new(flds%f_v,         data=zgmv(:,:,jend_uv,:))
  !if (jbegin_uv>0 )                    call field_new(flds%f_vor,         data=zgmv(:,:,jbegin_uv,:))
  !if (jend_uv>0 )                      call field_new(flds%f_div,         data=zgmv(:,:,jend_uv,:))

  endif
  
  if (size(sp3d,3) >=1 )   call field_new(flds%f_spvor,        data=sp3d(:,:,1))
  if (size(sp3d,3) >=2 )   call field_new(flds%f_spdiv,        data=sp3d(:,:,2))
  if (size(sp3d,3) >=3 )   call field_new(flds%f_spscalars,    data=sp3d(:,:,3:))
  if (size(spc2,2) >=1 ) call field_new(flds%f_spscalars2,  data=spc2(:,:))
  
  if (size(zgmvs,2)>=1)    call field_new(flds%f_scalars2,    data=zgmvs(:,1:1,:))

  if (luvders) then
    call field_new(flds%f_udm,         data=zgmv(:,:,jbegin_uder_ew,:))
    call field_new(flds%f_vdm,         data=zgmv(:,:,jend_uder_ew,:))
  endif

  if (jend_sc>0 .and. jend_sc>=jbegin_sc )  call field_new(flds%f_scalars,     data=zgmv(:,:,jbegin_sc:jend_sc,:))
  if (lscders) then    
    if (jend_scder_ew>0 .and. jend_scder_ew>=jbegin_scder_ew ) call field_new(flds%f_scalars_ew,  data=zgmv(:,:,jbegin_scder_ew:jend_scder_ew,:))
    if (jend_scder_ns>0 .and. jend_scder_ns>=jbegin_scder_ns ) call field_new(flds%f_scalars_ns,  data=zgmv(:,:,jbegin_scder_ns:jend_scder_ns,:))
  
    if (size(zgmvs,2)>=2)     call field_new(flds%f_scalars2_ns, data=zgmvs(:,2:2,:))
    if (size(zgmvs,2)>=3)     call field_new(flds%f_scalars2_ew, data=zgmvs(:,3:3,:))      
  endif
end subroutine wrap_fields

subroutine delete_wrapped_fields(wflds)
  type(wrapped_fields), intent(inout) :: wflds
   
  if(associated(wflds%f_spvor)) call field_delete(wflds%f_spvor)
  if(associated(wflds%f_spdiv)) call field_delete(wflds%f_spdiv)
  if(associated(wflds%f_spscalars)) call field_delete(wflds%f_spscalars)
  if(associated(wflds%f_spscalars2)) call field_delete(wflds%f_spscalars2)

  if(associated(wflds%f_u)) call field_delete(wflds%f_u)
  if(associated(wflds%f_v)) call field_delete(wflds%f_v)
  if(associated(wflds%f_udm)) call field_delete(wflds%f_udm)
  if(associated(wflds%f_vdm)) call field_delete(wflds%f_vdm)
  if(associated(wflds%f_scalars)) call field_delete(wflds%f_scalars)
  if(associated(wflds%f_scalars_ew)) call field_delete(wflds%f_scalars_ew)
  if(associated(wflds%f_scalars_ns)) call field_delete(wflds%f_scalars_ns)
  if(associated(wflds%f_vor)) call field_delete(wflds%f_vor)
  if(associated(wflds%f_div)) call field_delete(wflds%f_div)

  if(associated(wflds%f_scalars2)) call field_delete(wflds%f_scalars2)
  if(associated(wflds%f_scalars2_ew)) call field_delete(wflds%f_scalars2_ew)
  if(associated(wflds%f_scalars2_ns)) call field_delete(wflds%f_scalars2_ns) 
end subroutine delete_wrapped_fields

subroutine nullify_fields_lists(fl)
  type(fields_lists), intent(inout) :: fl
    
  nullify(fl%u)
  nullify(fl%v)
  nullify(fl%scalar)
  nullify(fl%spscalar)
  nullify(fl%spvor)
  nullify(fl%spdiv)
  nullify(fl%vor)
  nullify(fl%div)
  nullify(fl%udm)
  nullify(fl%vdm)
  nullify(fl%scalardm)
  nullify(fl%scalardl)
  
end subroutine nullify_fields_lists

subroutine delete_fields_lists(fl)
  type(fields_lists), intent(inout) ::fl
  if (associated(fl%u)) deallocate(fl%u)
  if (associated(fl%v)) deallocate(fl%v)
  if (associated(fl%scalar)) deallocate(fl%scalar)
  if (associated(fl%spscalar)) deallocate(fl%spscalar)
  if (associated(fl%spvor)) deallocate(fl%spvor)
  if (associated(fl%spdiv)) deallocate(fl%spdiv)
  if (associated(fl%vor)) deallocate(fl%vor)
  if (associated(fl%div)) deallocate(fl%div)
  if (associated(fl%udm)) deallocate(fl%udm)
  if (associated(fl%vdm)) deallocate(fl%vdm) 
  if (associated(fl%scalardm)) deallocate(fl%scalardm)
  if (associated(fl%scalardl)) deallocate(fl%scalardl)
  call nullify_fields_lists(fl) 
end subroutine delete_fields_lists

subroutine output_fields_lists(nout,fl)
  integer(kind=jpim), intent(in) :: nout
  type(fields_lists), intent(in) :: fl
  
  if (associated(fl%u)) write(nout,*) "fl%u", size(fl%u)  
  if (associated(fl%v)) write(nout,*) "fl%v", size(fl%v)
  if (associated(fl%scalar)) write(nout,*) "fl%scalar", size(fl%scalar)
  if (associated(fl%spscalar)) write(nout,*) "fl%spscalar", size(fl%spscalar)
  if (associated(fl%spvor)) write(nout,*) "fl%spvor", size(fl%spvor)
  if (associated(fl%spdiv)) write(nout,*) "fl%spdiv", size(fl%spdiv)
  if (associated(fl%vor)) write(nout,*) "fl%vor", size(fl%vor)
  if (associated(fl%div)) write(nout,*) "fl%div", size(fl%div)
  if (associated(fl%udm)) write(nout,*) "fl%udm", size(fl%udm)
  if (associated(fl%vdm)) write(nout,*) "fl%vdm", size(fl%vdm)
  if (associated(fl%scalardm)) write(nout,*) "fl%scalardm", size(fl%scalardm)
  if (associated(fl%scalardl)) write(nout,*) "fl%scalardl", size(fl%scalardl)
end subroutine output_fields_lists

subroutine create_fields_lists(flds,ylf, nbsetlev,nbsetsc)
  type(wrapped_fields), intent(in) :: flds
  type(fields_lists), intent(inout), target :: ylf
  integer(kind=jpim) :: nbsetlev(:)
  integer(kind=jpim) :: nbsetsc(:)
    
  if(associated(flds%f_spvor)) ylf%alloc_spvor=[b(flds%f_spvor,'sp_vor')]
  
  if(associated(flds%f_spdiv)) ylf%alloc_spdiv= [b(flds%f_spdiv,'spdiv')]
  
  if(associated(flds%f_u)) ylf%alloc_u = [b(flds%f_u,'u',nbsetlev)]
  if(associated(flds%f_v)) ylf%alloc_v = [b(flds%f_v,'v',nbsetlev)]
      
  if(associated(flds%f_udm)) ylf%alloc_udm=[b(flds%f_udm,'u_fdm', nbsetlev)]
  if(associated(flds%f_vdm)) ylf%alloc_vdm=[b(flds%f_vdm,'v_fdm', nbsetlev)]
    
  if(associated(flds%f_vor))  ylf%alloc_vor = [b(flds%f_vor,'vor', nbsetlev)]
  if(associated(flds%f_div))  ylf%alloc_div = [b(flds%f_div,'div', nbsetlev)]
  
  if (associated(flds%f_spscalars) .and. associated(flds%f_spscalars2) ) then
    ylf%alloc_spscalar = [b(flds%f_spscalars,'sp_scalars'), b(flds%f_spscalars2,'sp_scalars2')]
  else if (associated(flds%f_spscalars)) then
    ylf%alloc_spscalar = [b(flds%f_spscalars,'sp_scalars')]    
  else if (associated(flds%f_spscalars2)) then
    ylf%alloc_spscalar = [b(flds%f_spscalars2,'sp_scalars2')]  
  endif
    
  if (associated(flds%f_scalars) .and. associated(flds%f_scalars2) ) then
    ylf%alloc_scalar = [b(flds%f_scalars,'scalars', nbsetlev), b(flds%f_scalars2,'scalars2', nbsetsc)]
  else if (associated(flds%f_scalars)) then
    ylf%alloc_scalar = [b(flds%f_scalars,'scalars', nbsetlev)]    
  else if (associated(flds%f_scalars2)) then
    ylf%alloc_scalar = [b(flds%f_scalars2,'scalars2', nbsetsc)]  
  endif
  
  if (associated(flds%f_scalars_ns) .and. associated(flds%f_scalars2_ns) ) then
    ylf%alloc_scalardm = [b(flds%f_scalars_ns,'scalars_ns', nbsetlev), b(flds%f_scalars2_ns,'scalars2_ns', nbsetsc)]
  else if (associated(flds%f_scalars_ns)) then
    ylf%alloc_scalardm = [b(flds%f_scalars_ns,'scalars_ns', nbsetlev)]
  else if (associated(flds%f_scalars2_ns)) then
    ylf%alloc_scalardm = [b(flds%f_scalars2_ns,'scalars2_ns', nbsetsc)]  
  endif
  
  if (associated(flds%f_scalars_ew) .and. associated(flds%f_scalars2_ew) ) then
    ylf%alloc_scalardl = [b(flds%f_scalars_ew,'scalars_ew', nbsetlev), b(flds%f_scalars2_ew,'scalars2_ew', nbsetsc)]
  else if (associated(flds%f_scalars_ew)) then
    ylf%alloc_scalardl = [b(flds%f_scalars_ew,'scalars_ew', nbsetlev)]    
  else if (associated(flds%f_scalars2_ew)) then
    ylf%alloc_scalardl = [b(flds%f_scalars2_ew,'scalars2_ew', nbsetsc)]  
  endif

  if (allocated(ylf%alloc_u)) ylf%u=>ylf%alloc_u
  if (allocated(ylf%alloc_v)) ylf%v=>ylf%alloc_v
  if (allocated(ylf%alloc_scalar)) ylf%scalar=>ylf%alloc_scalar
  if (allocated(ylf%alloc_spscalar)) ylf%spscalar=>ylf%alloc_spscalar
  if (allocated(ylf%alloc_spvor)) ylf%spvor=>ylf%alloc_spvor
  if (allocated(ylf%alloc_spdiv)) ylf%spdiv=>ylf%alloc_spdiv
  if (allocated(ylf%alloc_vor)) ylf%vor=>ylf%alloc_vor
  if (allocated(ylf%alloc_div)) ylf%div=>ylf%alloc_div
  if (allocated(ylf%alloc_udm)) ylf%udm=>ylf%alloc_udm
  if (allocated(ylf%alloc_vdm)) ylf%vdm=>ylf%alloc_vdm 
  if (allocated(ylf%alloc_scalardm)) ylf%scalardm=>ylf%alloc_scalardm
  if (allocated(ylf%alloc_scalardl)) ylf%scalardl=>ylf%alloc_scalardl
end subroutine create_fields_lists

end module ectrans_field_api_helper