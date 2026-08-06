! (C) Copyright 2026- ECMWF.
! (C) Copyright 2026- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIR_TRANS_FIELD_VIEW(KRESOL,                       &
                              & YDGPSCALAR, YDGPU, YDGPV,     &
                              & YDSPSCALAR, YDSPVOR, YDSPDIV)

USE ISO_FORTRAN_ENV, ONLY : INT32
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW


USE YOMHOOK  ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1 ,ONLY : JPRB, JPIM
USE ECTRANS_FIELD_VIEW_MOD, ONLY: FIELD_VIEW_GET_IVSET_PTR
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, IVSET_PTR, &
                                              & GET_NPROMA, GET_NFLD, GET_NSPEC2, GET_NBLK

USE DIR_TRANS_VIEW_CTL_MOD, ONLY: DIR_TRANS_VIEW_CTL

IMPLICIT NONE

! Arguments
INTEGER(KIND=INT32), INTENT(IN)      :: KRESOL
TYPE(FIELD_VIEW),   INTENT(IN)      :: YDGPSCALAR(:), YDGPU(:), YDGPV(:)
TYPE(FIELD_VIEW),   INTENT(INOUT)   :: YDSPSCALAR(:), YDSPVOR(:), YDSPDIV(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: NPROMA, NBLK

! Local variables

! List of FIELD_VIEW: intermediate representation of fields to facilitate copy to temporary arrays
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVSCALAR(:)
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVVOR(:), YLSPVDIV(:)

TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVSCALAR(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVU(:),YLGVV(:)

REAL(KIND=JPRB):: S
INTEGER(KIND=JPIM) :: C

! b-set for dir-trans
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETUV_LIST(:)
TYPE(IVSET_PTR), ALLOCATABLE :: IVSETSC_LIST(:)

INTEGER(KIND=JPIM) :: IFLDSPVOR
INTEGER(KIND=JPIM) :: IUVG
INTEGER(KIND=JPIM) :: JFLD      ! Field counter

INTEGER(KIND=JPIM) :: KFLEVG

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_VIEW',0,ZHOOK_HANDLE)

NPROMA              = GET_NPROMA(YDGPU, YDGPV, YDGPSCALAR)
NBLK                = GET_NBLK(YDGPU, YDGPV, YDGPSCALAR)

! 1. Vector fields transformation to spectral space

! Preliminary checks

! Do we have vector fields?
IF (SIZE(YDGPU) > 0) THEN

  IF ((SIZE(YDGPU)/= SIZE(YDGPV)).OR.(SIZE(YDGPU)/= SIZE(YDSPDIV)).OR.(SIZE(YDGPU)/= SIZE(YDSPVOR))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] The vector arrays have inconsitent sizes: YDGPU, YDGPV, YDSPDIV, YDSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW
  IFLDSPVOR = LS_COUNT(YDSPVOR)

  ALLOCATE(YLSPVVOR(IFLDSPVOR))
  ALLOCATE(YLSPVDIV(IFLDSPVOR))

  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  ALLOCATE(YLGVU(LG_COUNT(YDGPU)))
  ALLOCATE(YLGVV(LG_COUNT(YDGPV)))
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] inconsistent number of field_view for vectors")
  ENDIF
  KFLEVG = SIZE (YLGVU) / SIZE (YDGPU)    
  IUVG = SIZE(YDGPU)
    
  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields.
  ALLOCATE(IVSETUV_LIST(IUVG))
  DO JFLD=1,IUVG
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPVOR(JFLD), IVSETUV_LIST(JFLD)%PTR)
  ENDDO

  C = LG(YDGPU, YLGVU, IVSETUV_LIST)
  C = LG(YDGPV, YLGVV, IVSETUV_LIST)  
  
ENDIF

! 2. scalar fields transformation

! Preliminary checks

! Do we have scalar fields?
IF (SIZE(YDSPSCALAR) > 0 ) THEN
  IF ((SIZE(YDSPSCALAR)/= SIZE(YDGPSCALAR)))  CALL ABOR1("[DIR_TRANS_FIELD_API] Inconsistent size &
                                                         & for YDSPSCALAR and YDGPSCALAR")

  ! Convert list of spectral scalar fields of any dimension into a list of 2d fields
  ALLOCATE(YLGVSCALAR(LG_COUNT(YDGPSCALAR)))

  ALLOCATE(YLSPVSCALAR(LS_COUNT(YDSPSCALAR)))

  ! count the number of fields present on the processor
  C = LS(YDSPSCALAR, YLSPVSCALAR)

  ! For LG we need the ivset of each grid point field,
  ! so we extract a matching list from the spectral fields
  ALLOCATE(IVSETSC_LIST(SIZE(YDSPSCALAR)))
  DO JFLD=1,SIZE(YDSPSCALAR)
    CALL FIELD_VIEW_GET_IVSET_PTR(YDSPSCALAR(JFLD), IVSETSC_LIST(JFLD)%PTR)
  ENDDO

  ! Copy list of scalar fields into temporary arrays (2d copy thanks to field_view)
  C = LG(YDGPSCALAR, YLGVSCALAR, IVSETSC_LIST)

ENDIF
IF (.NOT. ALLOCATED(YLSPVVOR))       ALLOCATE(YLSPVVOR(0))
IF (.NOT. ALLOCATED(YLSPVDIV))       ALLOCATE(YLSPVDIV(0))
IF (.NOT. ALLOCATED(YLSPVSCALAR))    ALLOCATE(YLSPVSCALAR(0)) 
IF (.NOT. ALLOCATED(YLGVU))          ALLOCATE(YLGVU(0))
IF (.NOT. ALLOCATED(YLGVV))          ALLOCATE(YLGVV(0))
IF (.NOT. ALLOCATED(YLGVSCALAR))     ALLOCATE(YLGVSCALAR(0))

CALL DIR_TRANS_VIEW_CTL(NPROMA,NBLK,YLSPVVOR, YLSPVDIV, YLSPVSCALAR,YLGVU, YLGVV, YLGVSCALAR)
! 3. CALL DIR_TRANS using the regular interface and the temporary arrays

! 5. Final cleanup

! delete temporary arrays

IF (ALLOCATED(IVSETUV_LIST))  DEALLOCATE(IVSETUV_LIST)
IF (ALLOCATED(IVSETSC_LIST))  DEALLOCATE(IVSETSC_LIST)

! delete FIELD_VIEWS
IF (ALLOCATED(YLSPVVOR))    DEALLOCATE(YLSPVVOR)
IF (ALLOCATED(YLSPVDIV))    DEALLOCATE(YLSPVDIV)
IF (ALLOCATED(YLSPVSCALAR)) DEALLOCATE(YLSPVSCALAR)
IF (ALLOCATED(YLGVU))       DEALLOCATE(YLGVU)
IF (ALLOCATED(YLGVV))       DEALLOCATE(YLGVV)
IF (ALLOCATED(YLGVSCALAR))  DEALLOCATE(YLGVSCALAR)


IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_VIEW',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_FIELD_VIEW
