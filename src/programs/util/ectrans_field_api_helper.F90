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

subroutine output_wrapped_fields(nout, wrap_fields)
  type(WRAPPED_FIELDS), INTENT(IN) :: wrap_fields
  integer(kind=jpim), INTENT(IN) :: nout

  write(nout,*) "wrap_fields%F_SPVOR", LOC(wrap_fields%F_SPVOR)
  write(nout,*) "wrap_fields%F_SPDIV", LOC(wrap_fields%F_SPDIV)
  write(nout,*) "wrap_fields%F_SPSCALARS", LOC(wrap_fields%F_SPSCALARS)
  write(nout,*) "wrap_fields%F_SPSCALARS2", LOC(wrap_fields%F_SPSCALARS2)

  write(nout,*) "wrap_fields%F_U", LOC(wrap_fields%F_U)
  write(nout,*) "wrap_fields%F_V", LOC(wrap_fields%F_V)
  write(nout,*) "wrap_fields%F_UDM", LOC(wrap_fields%F_UDM)
  write(nout,*) "wrap_fields%F_VDM", LOC(wrap_fields%F_VDM)
  write(nout,*) "wrap_fields%F_SCALARS", LOC(wrap_fields%F_SCALARS)
  write(nout,*) "wrap_fields%F_SCALARS_EW", LOC(wrap_fields%F_SCALARS_EW)
  write(nout,*) "wrap_fields%F_SCALARS_NS", LOC(wrap_fields%F_SCALARS_NS)
  write(nout,*) "wrap_fields%F_VOR", LOC(wrap_fields%F_VOR)
  write(nout,*) "wrap_fields%F_DIV", LOC(wrap_fields%F_DIV)

  write(nout,*) "wrap_fields%F_SCALARS2", LOC(wrap_fields%F_SCALARS2)
  write(nout,*) "wrap_fields%F_SCALARS2_EW", LOC(wrap_fields%F_SCALARS2_EW)
  write(nout,*) "wrap_fields%F_SCALARS2_NS", LOC(wrap_fields%F_SCALARS2_NS)
end subroutine output_wrapped_fields

subroutine nullify_wrapped_fields(wf)
  type(WRAPPED_FIELDS), INTENT(INOUT) :: wf
    
  NULLIFY(wf%F_SPVOR)
  NULLIFY(wf%F_SPDIV)
  NULLIFY(wf%F_SPSCALARS)
  NULLIFY(wf%F_SPSCALARS2)

  NULLIFY(wf%F_U)
  NULLIFY(wf%F_V)
  NULLIFY(wf%F_UDM)
  NULLIFY(wf%F_VDM)
  NULLIFY(wf%F_SCALARS)
  NULLIFY(wf%F_SCALARS_EW)
  NULLIFY(wf%F_SCALARS_NS)
  NULLIFY(wf%F_VOR)
  NULLIFY(wf%F_DIV)

  NULLIFY(wf%F_SCALARS2)
  NULLIFY(wf%F_SCALARS2_EW)
  NULLIFY(wf%F_SCALARS2_NS)
end subroutine nullify_wrapped_fields

subroutine wrap_fields(wf, lvordiv, lscders, luvders,& 
                    sp3d, spc2, zgmv, zgmvs, zgp2,& 
                    &jbegin_uv,jend_uv,& 
                    &jbegin_sc,jend_sc,& 
                    &jbegin_scder_NS, jend_scder_NS,& 
                    &jbegin_scder_EW, jend_scder_EW,& 
                    &jbegin_uder_EW, jend_uder_EW,& 
                    &jbegin_vder_EW, jend_vder_EW)

  type(WRAPPED_FIELDS), INTENT(INOUT) :: wf
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
    
  call nullify_wrapped_fields(wf)

  
  if (lvordiv) then
    IF (jbegin_uv>0 )      CALL FIELD_NEW(wf%F_U,         DATA=zgmv(:,:,jbegin_uv,:))
    IF (jend_uv>0 )        CALL FIELD_NEW(wf%F_V,         DATA=zgmv(:,:,jend_uv,:))
  !IF (jbegin_uv>0 )                    CALL FIELD_NEW(wf%F_VOR,         DATA=zgmv(:,:,jbegin_uv,:))
  !IF (jend_uv>0 )                      CALL FIELD_NEW(wf%F_DIV,         DATA=zgmv(:,:,jend_uv,:))

  endif
  
  IF (SIZE(sp3d,3) >=1 )   CALL FIELD_NEW(wf%F_SPVOR,        DATA=sp3d(:,:,1))
  IF (SIZE(sp3d,3) >=2 )   CALL FIELD_NEW(wf%F_SPDIV,        DATA=sp3d(:,:,2))
  IF (SIZE(sp3d,3) >=3 )   CALL FIELD_NEW(wf%F_SPSCALARS,    DATA=sp3d(:,:,3:))
  IF (SIZE(spc2,2) >=1 ) CALL FIELD_NEW(wf%F_SPSCALARS2,  DATA=spc2(:,:))
  
  IF (SIZE(zgmvs,2)>=1)    CALL FIELD_NEW(wf%F_SCALARS2,    DATA=zgmvs(:,1:1,:))

  if (luvders) then
    CALL FIELD_NEW(wf%F_UDM,         DATA=zgmv(:,:,jbegin_uder_EW,:))
    CALL FIELD_NEW(wf%F_VDM,         DATA=zgmv(:,:,jend_uder_EW,:))
  endif

  IF (jend_sc>0 .AND. jend_sc>=jbegin_sc )  CALL FIELD_NEW(wf%F_SCALARS,     DATA=zgmv(:,:,jbegin_sc:jend_sc,:))
  if (lscders) then    
    IF (jend_scder_EW>0 .AND. jend_scder_EW>=jbegin_scder_EW ) CALL FIELD_NEW(wf%F_SCALARS_EW,  DATA=zgmv(:,:,jbegin_scder_EW:jend_scder_EW,:))
    IF (jend_scder_NS>0 .AND. jend_scder_NS>=jbegin_scder_NS ) CALL FIELD_NEW(wf%F_SCALARS_NS,  DATA=zgmv(:,:,jbegin_scder_NS:jend_scder_NS,:))
  
    IF (SIZE(ZGMVS,2)>=2)     CALL FIELD_NEW(wf%F_SCALARS2_NS, DATA=zgmvs(:,2:2,:))
    IF (SIZE(ZGMVS,2)>=3)     CALL FIELD_NEW(wf%F_SCALARS2_EW, DATA=zgmvs(:,3:3,:))      
  endif
end subroutine wrap_fields

subroutine delete_wrapped_fields(wrap_fields)
  type(WRAPPED_FIELDS), INTENT(INOUT) :: wrap_fields
   
  IF(ASSOCIATED(wrap_fields%F_SPVOR)) CALL FIELD_DELETE(wrap_fields%F_SPVOR)
  IF(ASSOCIATED(wrap_fields%F_SPDIV)) CALL FIELD_DELETE(wrap_fields%F_SPDIV)
  IF(ASSOCIATED(wrap_fields%F_SPSCALARS)) CALL FIELD_DELETE(wrap_fields%F_SPSCALARS)
  IF(ASSOCIATED(wrap_fields%F_SPSCALARS2)) CALL FIELD_DELETE(wrap_fields%F_SPSCALARS2)

  IF(ASSOCIATED(wrap_fields%F_U)) CALL FIELD_DELETE(wrap_fields%F_U)
  IF(ASSOCIATED(wrap_fields%F_V)) CALL FIELD_DELETE(wrap_fields%F_V)
  IF(ASSOCIATED(wrap_fields%F_UDM)) CALL FIELD_DELETE(wrap_fields%F_UDM)
  IF(ASSOCIATED(wrap_fields%F_VDM)) CALL FIELD_DELETE(wrap_fields%F_VDM)
  IF(ASSOCIATED(wrap_fields%F_SCALARS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS)
  IF(ASSOCIATED(wrap_fields%F_SCALARS_EW)) CALL FIELD_DELETE(wrap_fields%F_SCALARS_EW)
  IF(ASSOCIATED(wrap_fields%F_SCALARS_NS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS_NS)
  IF(ASSOCIATED(wrap_fields%F_VOR)) CALL FIELD_DELETE(wrap_fields%F_VOR)
  IF(ASSOCIATED(wrap_fields%F_DIV)) CALL FIELD_DELETE(wrap_fields%F_DIV)

  IF(ASSOCIATED(wrap_fields%F_SCALARS2)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2)
  IF(ASSOCIATED(wrap_fields%F_SCALARS2_EW)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2_EW)
  IF(ASSOCIATED(wrap_fields%F_SCALARS2_NS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2_NS) 
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

subroutine create_fields_lists(wf,ylf, NBSETLEV,NBSETSC)
  type(WRAPPED_FIELDS), INTENT(IN) :: wf
  type(FIELDS_LISTS), INTENT(INOUT), TARGET :: ylf
  INTEGER(KIND=JPIM) :: NBSETLEV(:)
  INTEGER(KIND=JPIM) :: NBSETSC(:)
    
  IF(ASSOCIATED(wf%F_SPVOR)) ylf%ALLOC_SPVOR=[B(wf%F_SPVOR,'SP_VOR')]
  
  IF(ASSOCIATED(wf%F_SPDIV)) ylf%ALLOC_SPDIV= [B(wf%F_SPDIV,'SPDIV')]
  
  IF(ASSOCIATED(wf%F_U)) ylf%ALLOC_U = [B(wf%F_U,'U',NBSETLEV)]
  IF(ASSOCIATED(wf%F_V)) ylf%ALLOC_V = [B(wf%F_V,'V',NBSETLEV)]
      
  IF(ASSOCIATED(wf%F_UDM)) ylf%ALLOC_UDM=[B(wf%F_UDM,'U_FDM', NBSETLEV)]
  IF(ASSOCIATED(wf%F_VDM)) ylf%ALLOC_VDM=[B(wf%F_VDM,'V_FDM', NBSETLEV)]
    
  IF(ASSOCIATED(wf%F_VOR))  ylf%ALLOC_VOR = [B(wf%F_VOR,'VOR', NBSETLEV)]
  IF(ASSOCIATED(wf%F_DIV))  ylf%ALLOC_DIV = [B(wf%F_DIV,'DIV', NBSETLEV)]
  
  IF (ASSOCIATED(wf%F_SPSCALARS) .AND. ASSOCIATED(wf%F_SPSCALARS2) ) THEN
    ylf%ALLOC_SPSCALAR = [B(wf%F_SPSCALARS,'SP_SCALARS'), B(wf%F_SPSCALARS2,'SP_SCALARS2')]
  ELSE IF (ASSOCIATED(wf%F_SPSCALARS)) THEN
    ylf%ALLOC_SPSCALAR = [B(wf%F_SPSCALARS,'SP_SCALARS')]    
  ELSE IF (ASSOCIATED(wf%F_SPSCALARS2)) THEN
    ylf%ALLOC_SPSCALAR = [B(wf%F_SPSCALARS2,'SP_SCALARS2')]  
  ENDIF
    
  IF (ASSOCIATED(wf%F_SCALARS) .AND. ASSOCIATED(wf%F_SCALARS2) ) THEN
    ylf%ALLOC_SCALAR = [B(wf%F_SCALARS,'SCALARS', NBSETLEV), B(wf%F_SCALARS2,'SCALARS2', NBSETSC)]
  ELSE IF (ASSOCIATED(wf%F_SCALARS)) THEN
    ylf%ALLOC_SCALAR = [B(wf%F_SCALARS,'SCALARS', NBSETLEV)]    
  ELSE IF (ASSOCIATED(wf%F_SCALARS2)) THEN
    ylf%ALLOC_SCALAR = [B(wf%F_SCALARS2,'SCALARS2', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(wf%F_SCALARS_NS) .AND. ASSOCIATED(wf%F_SCALARS2_NS) ) THEN
    ylf%ALLOC_SCALARDM = [B(wf%F_SCALARS_NS,'SCALARS_NS', NBSETLEV), B(wf%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]
  ELSE IF (ASSOCIATED(wf%F_SCALARS_NS)) THEN
    ylf%ALLOC_SCALARDM = [B(wf%F_SCALARS_NS,'SCALARS_NS', NBSETLEV)]
  ELSE IF (ASSOCIATED(wf%F_SCALARS2_NS)) THEN
    ylf%ALLOC_SCALARDM = [B(wf%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(wf%F_SCALARS_EW) .AND. ASSOCIATED(wf%F_SCALARS2_EW) ) THEN
    ylf%ALLOC_SCALARDL = [B(wf%F_SCALARS_EW,'SCALARS_EW', NBSETLEV), B(wf%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]
  ELSE IF (ASSOCIATED(wf%F_SCALARS_EW)) THEN
    ylf%ALLOC_SCALARDL = [B(wf%F_SCALARS_EW,'SCALARS_EW', NBSETLEV)]    
  ELSE IF (ASSOCIATED(wf%F_SCALARS2_EW)) THEN
    ylf%ALLOC_SCALARDL = [B(wf%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]  
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