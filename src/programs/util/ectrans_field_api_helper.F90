module ectrans_field_api_helper

USE FIELD_MODULE, ONLY:FIELD_1RB, FIELD_2RB, FIELD_3RB, FIELD_4RB
USE FIELD_FACTORY_MODULE
use parkind1, only: jpim, jprb, jprd
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"

implicit none

TYPE WRAPPED_FIELDS
  CLASS (FIELD_3RB), POINTER :: F_SPSCALARS
  CLASS (FIELD_2RB), POINTER :: F_SPSCALARS2
  CLASS (FIELD_2RB), POINTER :: F_SPVOR, F_SPDIV

  CLASS (FIELD_3RB), POINTER :: F_VOR, F_DIV
  CLASS (FIELD_3RB), POINTER :: F_U, F_V
  CLASS (FIELD_3RB), POINTER :: F_UDM, F_VDM

  CLASS (FIELD_4RB), POINTER :: F_SCALARS
  CLASS (FIELD_4RB), POINTER :: F_SCALARS_EW
  CLASS (FIELD_4RB), POINTER :: F_SCALARS_NS

  CLASS (FIELD_3RB), POINTER :: F_SCALARS2
  CLASS (FIELD_3RB), POINTER :: F_SCALARS2_EW
  CLASS (FIELD_3RB), POINTER :: F_SCALARS2_NS
END TYPE WRAPPED_FIELDS

TYPE FIELDS_LISTS
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_U (:), ALLOC_V (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_SCALAR (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_SPVOR (:), ALLOC_SPDIV (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_VOR (:), ALLOC_DIV (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_SPSCALAR (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_UDM (:), ALLOC_VDM (:)
  TYPE (FIELD_BASIC_PTR), ALLOCATABLE:: ALLOC_SCALARDM (:), ALLOC_SCALARDL (:)  

  TYPE (FIELD_BASIC_PTR), POINTER:: U (:), V (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SCALAR (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SPVOR (:), SPDIV (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: VOR (:), DIV (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SPSCALAR (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: UDM (:), VDM (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SCALARDM (:), SCALARDL (:)  
  END TYPE FIELDS_LISTS
  
contains

subroutine output_wrapped_fields(nout, wflds)
  type(WRAPPED_FIELDS), INTENT(IN) :: wflds
  integer(kind=jpim), INTENT(IN) :: nout

  write(nout,*) "wflds%F_SPVOR", LOC(wflds%F_SPVOR)
  write(nout,*) "wflds%F_SPDIV", LOC(wflds%F_SPDIV)
  write(nout,*) "wflds%F_SPSCALARS", LOC(wflds%F_SPSCALARS)
  write(nout,*) "wflds%F_SPSCALARS2", LOC(wflds%F_SPSCALARS2)

  write(nout,*) "wflds%F_U", LOC(wflds%F_U)
  write(nout,*) "wflds%F_V", LOC(wflds%F_V)
  write(nout,*) "wflds%F_UDM", LOC(wflds%F_UDM)
  write(nout,*) "wflds%F_VDM", LOC(wflds%F_VDM)
  write(nout,*) "wflds%F_SCALARS", LOC(wflds%F_SCALARS)
  write(nout,*) "wflds%F_SCALARS_EW", LOC(wflds%F_SCALARS_EW)
  write(nout,*) "wflds%F_SCALARS_NS", LOC(wflds%F_SCALARS_NS)
  write(nout,*) "wflds%F_VOR", LOC(wflds%F_VOR)
  write(nout,*) "wflds%F_DIV", LOC(wflds%F_DIV)

  write(nout,*) "wflds%F_SCALARS2", LOC(wflds%F_SCALARS2)
  write(nout,*) "wflds%F_SCALARS2_EW", LOC(wflds%F_SCALARS2_EW)
  write(nout,*) "wflds%F_SCALARS2_NS", LOC(wflds%F_SCALARS2_NS)
end subroutine output_wrapped_fields

subroutine nullify_wrapped_fields(flds)
  type(WRAPPED_FIELDS), INTENT(INOUT) :: flds
    
  NULLIFY(flds%F_SPVOR)
  NULLIFY(flds%F_SPDIV)
  NULLIFY(flds%F_SPSCALARS)
  NULLIFY(flds%F_SPSCALARS2)

  NULLIFY(flds%F_U)
  NULLIFY(flds%F_V)
  NULLIFY(flds%F_UDM)
  NULLIFY(flds%F_VDM)
  NULLIFY(flds%F_SCALARS)
  NULLIFY(flds%F_SCALARS_EW)
  NULLIFY(flds%F_SCALARS_NS)
  NULLIFY(flds%F_VOR)
  NULLIFY(flds%F_DIV)

  NULLIFY(flds%F_SCALARS2)
  NULLIFY(flds%F_SCALARS2_EW)
  NULLIFY(flds%F_SCALARS2_NS)
end subroutine nullify_wrapped_fields

subroutine wrap_fields(flds, lvordiv, lscders, luvders,& 
                    sp3d, spc2, zgmv, zgmvs, zgp2,& 
                    &jbegin_uv,jend_uv,& 
                    &jbegin_sc,jend_sc,& 
                    &jbegin_scder_NS, jend_scder_NS,& 
                    &jbegin_scder_EW, jend_scder_EW,& 
                    &jbegin_uder_EW, jend_uder_EW,& 
                    &jbegin_vder_EW, jend_vder_EW)

  type(WRAPPED_FIELDS), INTENT(INOUT) :: flds
  logical :: lvordiv
  logical :: lscders
  logical :: luvders
  real(kind=jprb), INTENT(IN) :: sp3d(:,:,:)
  real(kind=jprb), INTENT(IN) :: spc2(:,:)
  real(kind=jprb), INTENT(IN) :: zgmv(:,:,:,:)
  real(kind=jprb), INTENT(IN) :: zgmvs(:,:,:)
  real(kind=jprb), INTENT(IN) :: zgp2 (:,:,:)
  integer(kind=jpim) :: jbegin_uv
  integer(kind=jpim) :: jend_uv
  integer(kind=jpim) :: jbegin_sc
  integer(kind=jpim) :: jend_sc
  integer(kind=jpim) :: jbegin_scder_NS
  integer(kind=jpim) :: jend_scder_NS
  integer(kind=jpim) :: jbegin_scder_EW
  integer(kind=jpim) :: jend_scder_EW
  integer(kind=jpim) :: jbegin_uder_EW
  integer(kind=jpim) :: jend_uder_EW
  integer(kind=jpim) :: jbegin_vder_EW
  integer(kind=jpim) :: jend_vder_EW
    
  call nullify_wrapped_fields(flds)

  
  if (lvordiv) then
    IF (jbegin_uv>0 )      CALL FIELD_NEW(flds%F_U,         DATA=zgmv(:,:,jbegin_uv,:))
    IF (jend_uv>0 )        CALL FIELD_NEW(flds%F_V,         DATA=zgmv(:,:,jend_uv,:))
  !IF (jbegin_uv>0 )                    CALL FIELD_NEW(flds%F_VOR,         DATA=zgmv(:,:,jbegin_uv,:))
  !IF (jend_uv>0 )                      CALL FIELD_NEW(flds%F_DIV,         DATA=zgmv(:,:,jend_uv,:))

  endif
  
  IF (SIZE(sp3d,3) >=1 )   CALL FIELD_NEW(flds%F_SPVOR,        DATA=sp3d(:,:,1))
  IF (SIZE(sp3d,3) >=2 )   CALL FIELD_NEW(flds%F_SPDIV,        DATA=sp3d(:,:,2))
  IF (SIZE(sp3d,3) >=3 )   CALL FIELD_NEW(flds%F_SPSCALARS,    DATA=sp3d(:,:,3:))
  IF (SIZE(spc2,2) >=1 ) CALL FIELD_NEW(flds%F_SPSCALARS2,  DATA=spc2(:,:))
  
  IF (SIZE(zgmvs,2)>=1)    CALL FIELD_NEW(flds%F_SCALARS2,    DATA=zgmvs(:,1:1,:))

  if (luvders) then
    CALL FIELD_NEW(flds%F_UDM,         DATA=zgmv(:,:,jbegin_uder_EW,:))
    CALL FIELD_NEW(flds%F_VDM,         DATA=zgmv(:,:,jend_uder_EW,:))
  endif

  IF (jend_sc>0 .AND. jend_sc>=jbegin_sc )  CALL FIELD_NEW(flds%F_SCALARS,     DATA=zgmv(:,:,jbegin_sc:jend_sc,:))
  if (lscders) then    
    IF (jend_scder_EW>0 .AND. jend_scder_EW>=jbegin_scder_EW ) CALL FIELD_NEW(flds%F_SCALARS_EW,  DATA=zgmv(:,:,jbegin_scder_EW:jend_scder_EW,:))
    IF (jend_scder_NS>0 .AND. jend_scder_NS>=jbegin_scder_NS ) CALL FIELD_NEW(flds%F_SCALARS_NS,  DATA=zgmv(:,:,jbegin_scder_NS:jend_scder_NS,:))
  
    IF (SIZE(ZGMVS,2)>=2)     CALL FIELD_NEW(flds%F_SCALARS2_NS, DATA=zgmvs(:,2:2,:))
    IF (SIZE(ZGMVS,2)>=3)     CALL FIELD_NEW(flds%F_SCALARS2_EW, DATA=zgmvs(:,3:3,:))      
  endif
end subroutine wrap_fields

subroutine delete_wrapped_fields(wflds)
  type(WRAPPED_FIELDS), INTENT(INOUT) :: wflds
   
  IF(ASSOCIATED(wflds%F_SPVOR)) CALL FIELD_DELETE(wflds%F_SPVOR)
  IF(ASSOCIATED(wflds%F_SPDIV)) CALL FIELD_DELETE(wflds%F_SPDIV)
  IF(ASSOCIATED(wflds%F_SPSCALARS)) CALL FIELD_DELETE(wflds%F_SPSCALARS)
  IF(ASSOCIATED(wflds%F_SPSCALARS2)) CALL FIELD_DELETE(wflds%F_SPSCALARS2)

  IF(ASSOCIATED(wflds%F_U)) CALL FIELD_DELETE(wflds%F_U)
  IF(ASSOCIATED(wflds%F_V)) CALL FIELD_DELETE(wflds%F_V)
  IF(ASSOCIATED(wflds%F_UDM)) CALL FIELD_DELETE(wflds%F_UDM)
  IF(ASSOCIATED(wflds%F_VDM)) CALL FIELD_DELETE(wflds%F_VDM)
  IF(ASSOCIATED(wflds%F_SCALARS)) CALL FIELD_DELETE(wflds%F_SCALARS)
  IF(ASSOCIATED(wflds%F_SCALARS_EW)) CALL FIELD_DELETE(wflds%F_SCALARS_EW)
  IF(ASSOCIATED(wflds%F_SCALARS_NS)) CALL FIELD_DELETE(wflds%F_SCALARS_NS)
  IF(ASSOCIATED(wflds%F_VOR)) CALL FIELD_DELETE(wflds%F_VOR)
  IF(ASSOCIATED(wflds%F_DIV)) CALL FIELD_DELETE(wflds%F_DIV)

  IF(ASSOCIATED(wflds%F_SCALARS2)) CALL FIELD_DELETE(wflds%F_SCALARS2)
  IF(ASSOCIATED(wflds%F_SCALARS2_EW)) CALL FIELD_DELETE(wflds%F_SCALARS2_EW)
  IF(ASSOCIATED(wflds%F_SCALARS2_NS)) CALL FIELD_DELETE(wflds%F_SCALARS2_NS) 
end subroutine delete_wrapped_fields

subroutine nullify_fields_lists(fl)
  type(FIELDS_LISTS), INTENT(INOUT) :: fl
    
  NULLIFY(fl%U)
  NULLIFY(fl%V)
  NULLIFY(fl%SCALAR)
  NULLIFY(fl%SPSCALAR)
  NULLIFY(fl%SPVOR)
  NULLIFY(fl%SPDIV)
  NULLIFY(fl%VOR)
  NULLIFY(fl%DIV)
  NULLIFY(fl%UDM)
  NULLIFY(fl%VDM)
  NULLIFY(fl%SCALARDM)
  NULLIFY(fl%SCALARDL)
  
end subroutine nullify_fields_lists

subroutine delete_fields_lists(fl)
  type(FIELDS_LISTS), INTENT(INOUT) ::fl
  IF (ASSOCIATED(fl%U)) DEALLOCATE(fl%U)
  IF (ASSOCIATED(fl%V)) DEALLOCATE(fl%V)
  IF (ASSOCIATED(fl%SCALAR)) DEALLOCATE(fl%SCALAR)
  IF (ASSOCIATED(fl%SPSCALAR)) DEALLOCATE(fl%SPSCALAR)
  IF (ASSOCIATED(fl%SPVOR)) DEALLOCATE(fl%SPVOR)
  IF (ASSOCIATED(fl%SPDIV)) DEALLOCATE(fl%SPDIV)
  IF (ASSOCIATED(fl%VOR)) DEALLOCATE(fl%VOR)
  IF (ASSOCIATED(fl%DIV)) DEALLOCATE(fl%DIV)
  IF (ASSOCIATED(fl%UDM)) DEALLOCATE(fl%UDM)
  IF (ASSOCIATED(fl%VDM)) DEALLOCATE(fl%VDM) 
  IF (ASSOCIATED(fl%SCALARDM)) DEALLOCATE(fl%SCALARDM)
  IF (ASSOCIATED(fl%SCALARDL)) DEALLOCATE(fl%SCALARDL)
  call nullify_fields_lists(fl) 
end subroutine delete_fields_lists

subroutine output_fields_lists(nout,fl)
  integer(kind=jpim), INTENT(IN) :: nout
  type(FIELDS_LISTS), INTENT(IN) :: fl
  
  IF (ASSOCIATED(fl%U)) write(nout,*) "fl%U", SIZE(fl%U)  
  IF (ASSOCIATED(fl%V)) write(nout,*) "fl%V", SIZE(fl%V)
  IF (ASSOCIATED(fl%SCALAR)) write(nout,*) "fl%SCALAR", SIZE(fl%SCALAR)
  IF (ASSOCIATED(fl%SPSCALAR)) write(nout,*) "fl%SPSCALAR", SIZE(fl%SPSCALAR)
  IF (ASSOCIATED(fl%SPVOR)) write(nout,*) "fl%SPVOR", SIZE(fl%SPVOR)
  IF (ASSOCIATED(fl%SPDIV)) write(nout,*) "fl%SPDIV", SIZE(fl%SPDIV)
  IF (ASSOCIATED(fl%VOR)) write(nout,*) "fl%VOR", SIZE(fl%VOR)
  IF (ASSOCIATED(fl%DIV)) write(nout,*) "fl%DIV", SIZE(fl%DIV)
  IF (ASSOCIATED(fl%UDM)) write(nout,*) "fl%UDM", SIZE(fl%UDM)
  IF (ASSOCIATED(fl%VDM)) write(nout,*) "fl%VDM", SIZE(fl%VDM)
  IF (ASSOCIATED(fl%SCALARDM)) write(nout,*) "fl%SCALARDM", SIZE(fl%SCALARDM)
  IF (ASSOCIATED(fl%SCALARDL)) write(nout,*) "fl%SCALARDL", SIZE(fl%SCALARDL)
end subroutine output_fields_lists

subroutine create_fields_lists(flds,ylf, NBSETLEV,NBSETSC)
  type(WRAPPED_FIELDS), INTENT(IN) :: flds
  type(FIELDS_LISTS), INTENT(INOUT), TARGET :: ylf
  INTEGER(KIND=JPIM) :: NBSETLEV(:)
  INTEGER(KIND=JPIM) :: NBSETSC(:)
    
  IF(ASSOCIATED(flds%F_SPVOR)) ylf%ALLOC_SPVOR=[B(flds%F_SPVOR,'SP_VOR')]
  
  IF(ASSOCIATED(flds%F_SPDIV)) ylf%ALLOC_SPDIV= [B(flds%F_SPDIV,'SPDIV')]
  
  IF(ASSOCIATED(flds%F_U)) ylf%ALLOC_U = [B(flds%F_U,'U',NBSETLEV)]
  IF(ASSOCIATED(flds%F_V)) ylf%ALLOC_V = [B(flds%F_V,'V',NBSETLEV)]
      
  IF(ASSOCIATED(flds%F_UDM)) ylf%ALLOC_UDM=[B(flds%F_UDM,'U_FDM', NBSETLEV)]
  IF(ASSOCIATED(flds%F_VDM)) ylf%ALLOC_VDM=[B(flds%F_VDM,'V_FDM', NBSETLEV)]
    
  IF(ASSOCIATED(flds%F_VOR))  ylf%ALLOC_VOR = [B(flds%F_VOR,'VOR', NBSETLEV)]
  IF(ASSOCIATED(flds%F_DIV))  ylf%ALLOC_DIV = [B(flds%F_DIV,'DIV', NBSETLEV)]
  
  IF (ASSOCIATED(flds%F_SPSCALARS) .AND. ASSOCIATED(flds%F_SPSCALARS2) ) THEN
    ylf%ALLOC_SPSCALAR = [B(flds%F_SPSCALARS,'SP_SCALARS'), B(flds%F_SPSCALARS2,'SP_SCALARS2')]
  ELSE IF (ASSOCIATED(flds%F_SPSCALARS)) THEN
    ylf%ALLOC_SPSCALAR = [B(flds%F_SPSCALARS,'SP_SCALARS')]    
  ELSE IF (ASSOCIATED(flds%F_SPSCALARS2)) THEN
    ylf%ALLOC_SPSCALAR = [B(flds%F_SPSCALARS2,'SP_SCALARS2')]  
  ENDIF
    
  IF (ASSOCIATED(flds%F_SCALARS) .AND. ASSOCIATED(flds%F_SCALARS2) ) THEN
    ylf%ALLOC_SCALAR = [B(flds%F_SCALARS,'SCALARS', NBSETLEV), B(flds%F_SCALARS2,'SCALARS2', NBSETSC)]
  ELSE IF (ASSOCIATED(flds%F_SCALARS)) THEN
    ylf%ALLOC_SCALAR = [B(flds%F_SCALARS,'SCALARS', NBSETLEV)]    
  ELSE IF (ASSOCIATED(flds%F_SCALARS2)) THEN
    ylf%ALLOC_SCALAR = [B(flds%F_SCALARS2,'SCALARS2', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(flds%F_SCALARS_NS) .AND. ASSOCIATED(flds%F_SCALARS2_NS) ) THEN
    ylf%ALLOC_SCALARDM = [B(flds%F_SCALARS_NS,'SCALARS_NS', NBSETLEV), B(flds%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]
  ELSE IF (ASSOCIATED(flds%F_SCALARS_NS)) THEN
    ylf%ALLOC_SCALARDM = [B(flds%F_SCALARS_NS,'SCALARS_NS', NBSETLEV)]
  ELSE IF (ASSOCIATED(flds%F_SCALARS2_NS)) THEN
    ylf%ALLOC_SCALARDM = [B(flds%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(flds%F_SCALARS_EW) .AND. ASSOCIATED(flds%F_SCALARS2_EW) ) THEN
    ylf%ALLOC_SCALARDL = [B(flds%F_SCALARS_EW,'SCALARS_EW', NBSETLEV), B(flds%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]
  ELSE IF (ASSOCIATED(flds%F_SCALARS_EW)) THEN
    ylf%ALLOC_SCALARDL = [B(flds%F_SCALARS_EW,'SCALARS_EW', NBSETLEV)]    
  ELSE IF (ASSOCIATED(flds%F_SCALARS2_EW)) THEN
    ylf%ALLOC_SCALARDL = [B(flds%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]  
  ENDIF

  IF (ALLOCATED(ylf%ALLOC_U)) ylf%U=>ylf%ALLOC_U
  IF (ALLOCATED(ylf%ALLOC_V)) ylf%V=>ylf%ALLOC_V
  IF (ALLOCATED(ylf%ALLOC_SCALAR)) ylf%SCALAR=>ylf%ALLOC_SCALAR
  IF (ALLOCATED(ylf%ALLOC_SPSCALAR)) ylf%SPSCALAR=>ylf%ALLOC_SPSCALAR
  IF (ALLOCATED(ylf%ALLOC_SPVOR)) ylf%SPVOR=>ylf%ALLOC_SPVOR
  IF (ALLOCATED(ylf%ALLOC_SPDIV)) ylf%SPDIV=>ylf%ALLOC_SPDIV
  IF (ALLOCATED(ylf%ALLOC_VOR)) ylf%VOR=>ylf%ALLOC_VOR
  IF (ALLOCATED(ylf%ALLOC_DIV)) ylf%DIV=>ylf%ALLOC_DIV
  IF (ALLOCATED(ylf%ALLOC_UDM)) ylf%UDM=>ylf%ALLOC_UDM
  IF (ALLOCATED(ylf%ALLOC_VDM)) ylf%VDM=>ylf%ALLOC_VDM 
  IF (ALLOCATED(ylf%ALLOC_SCALARDM)) ylf%SCALARDM=>ylf%ALLOC_SCALARDM
  IF (ALLOCATED(ylf%ALLOC_SCALARDL)) ylf%SCALARDL=>ylf%ALLOC_SCALARDL

end subroutine create_fields_lists


end module ectrans_field_api_helper