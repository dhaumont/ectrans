! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIR_TRANS_FIELD_API(YDFSPVOR,YDFSPDIV,YDFSPSCALAR, &
                             & YDFU, YDFV, YDFSCALAR, &
                             & KSPEC, KPROMA, KGPBLKS, KGPTOT, KFLEVG, KFLEVL,&
                             & LDACC)


!**** *DIR_TRANS_FIELD_API* - Field API interface to direct spectral transform

!     Purpose.
!     --------
!        Allow to call DIR_TRANS with a list of fields from field API

!**   Interface.
!     ----------
!     CALL DIR_TRANS_FIELD_API(...)

!     Explicit arguments :
!     --------------------
!      output
!       YDFSPVOR(:)    - List of spectral vector fields (vorticity)
!       YDFSPDIV(:)    - List of spectral vector fields (divergence)
!       YDFSPSCALAR(:) - List of spectral scalar fields
!      input
!       YDFU(:)        - List of grid-point vector fields (u)
!       YDFV(:)        - List of grid-point vector fields (v)
!       YDFSCALAR(:)   - List of grid-point scalar fields
!       KSPEC          - Number of spectral coefficients
!       KPROMA         - Blocking factor
!       KGPBLKS        - Number of blocks
!       KGPTOT         - Number of total grid points
!       KFLEVG         - Number of levels
!       KFLEVL         - Number of local levels
!       LDACC          - Field data on device

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE FIELD_API_BASIC_TYPE_MOD, ONLY: FIELD_BASIC_PTR
USE FIELD_API_ECTRANS_MOD
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)

INTEGER(KIND=JPIM), INTENT(IN) ::KSPEC
INTEGER(KIND=JPIM), INTENT(IN) ::KPROMA
INTEGER(KIND=JPIM), INTENT(IN) ::KGPBLKS
INTEGER(KIND=JPIM), INTENT(IN) ::KGPTOT
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVG
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVL
LOGICAL, INTENT(IN), OPTIONAL  :: LDACC

! Local variables

! List of FIELD_VIEW: intermediate representation of fields to facilitate copy to temporary arrays
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVSCALAR(:)

! Temporary arrays for dir_trans
REAL(KIND=JPRB),POINTER :: ZPSPVOR(:,:),ZPSPDIV(:,:)  ! spectral vector fields (out)
REAL(KIND=JPRB),POINTER :: ZPSPSC2(:,:)               ! spectral scalar fields(out)
REAL(KIND=JPRB),POINTER :: ZPGPUV(:,:,:,:)            ! grid vector fields (in)
REAL(KIND=JPRB),POINTER :: ZPGP2(:,:,:)               ! grid scalar fields (in)

! b-set for dir-trans
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETUV(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETSC2(:)

INTEGER(KIND=JPIM) :: ISPUV
INTEGER(KIND=JPIM) :: IFLDXG
INTEGER(KIND=JPIM) :: IFLDXL
INTEGER(KIND=JPIM) :: IFLDXGUV
INTEGER(KIND=JPIM) :: IFLDXLUV
INTEGER(KIND=JPIM) :: IFLDSPVOR
INTEGER(KIND=JPIM) :: IFLDSPSC
INTEGER(KIND=JPIM) :: IUVG
INTEGER(KIND=JPIM) :: ISCDIM
INTEGER(KIND=JPIM) :: IUVDIM
INTEGER(KIND=JPIM) :: ID
INTEGER(KIND=JPIM) :: IOFFSET
INTEGER(KIND=JPIM) :: JLEV      ! Level counter
INTEGER(KIND=JPIM) :: JFLD      ! Field counter

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "dir_trans.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',0,ZHOOK_HANDLE)

ISPUV = 0
IFLDXG  = 0
IFLDXL = 0
IFLDXGUV = 0
IFLDXLUV = 0
IFLDSPVOR= 0
IFLDSPSC= 0
IUVG = 0
JFLD  = 0
ISCDIM = 0
IUVDIM = 0
ID = 0
IOFFSET = 0
JLEV = 0
JFLD = 0

! 1. Vector fields transformation to spectral space

! Preliminary checks
IF (PRESENT(YDFU) .NEQV. PRESENT(YDFV)) CALL ABOR1("[DIR_TRANS_FIELD_API] YDFU and YDFV must be provided together")

! Do we have vector fields?
IF (PRESENT(YDFU)) THEN

  IF ((SIZE(YDFU)/= SIZE(YDFV)).OR.(SIZE(YDFU)/= SIZE(YDFSPDIV)).OR.(SIZE(YDFU)/= SIZE(YDFSPVOR))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] The vector arrays have inconsitent sizes: YDFU, YDFV, YDFSPDIV, YDFSPVOR")
  ENDIF

  ! Convert list of spectral vector fields into a list of 2d FIELD_VIEW
  YLSPVVOR = LS(YDFSPVOR, LDACC, .FALSE.)
  YLSPVDIV = LS(YDFSPDIV, LDACC, .FALSE.)

  ! Convert list of grid-point vector fields into a list of 2d FIELD_VIEW
  YLGVU = LG(YDFU, LDACC, .TRUE.)
  YLGVV = LG(YDFV, LDACC, .TRUE.)

  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] inconsistent number of field_view for vectors")
  ENDIF
   IF (((SIZE (YLGVU) / SIZE (YDFU)) /= KFLEVG) .OR. ((SIZE (YLSPVVOR) / SIZE (YDFSPVOR)) /= KFLEVL)) THEN
     CALL ABOR1("[DIR_TRANS_FIELD_API] inconsistent kflevg or kflevl")
  ENDIF

  IUVG = SIZE(YDFU)
  ISPUV = SIZE(YDFSPVOR)

  IUVDIM = 2

  IFLDSPVOR = SIZE(YLSPVVOR)
  ! allocate temporary vector field arrays in spectral space
  ALLOCATE(ZPSPVOR(IFLDSPVOR,KSPEC))
  ALLOCATE(ZPSPDIV(IFLDSPVOR,KSPEC))

  ! allocate temporary vector field array in grid space
  ALLOCATE(ZPGPUV(KPROMA,KFLEVG, IUVG * IUVDIM,KGPBLKS))

  ! allocate 'b-set' for vector fields
  ALLOCATE(IVSETUV(KFLEVG))


! temporary copies on gpu
   if ( LDACC ) THEN
    !$ACC ENTER DATA CREATE(ZPSPVOR,ZPSPDIV,ZPGPUV)
  ENDIF


  IOFFSET = 0

  ! Copy list of 2d views of grid point vector fields into temporary arrays
  IF (LDACC) THEN
    !$ACC KERNELS PRESENT(ZPGPUV,YLGVU(:)%VIEW%P,YLGVV(:)%VIEW%P) FIRSTPRIVATE(IUVG, KFLEVG,IOFFSET) PRIVATE(JFLD,JLEV,ID)
    DO JFLD=1,IUVG
      DO JLEV=1,KFLEVG
        ID = JLEV + (JFLD -1) * KFLEVG
        ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:) = YLGVU(ID)%VIEW%P(:,:)
        ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:) = YLGVV(ID)%VIEW%P(:,:)
      ENDDO
    ENDDO
    !$ACC END KERNELS
  ELSE
      DO JFLD=1,IUVG
        DO JLEV=1,KFLEVG
          ID = JLEV + (JFLD -1) * KFLEVG
          ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:) = YLGVU(ID)%VIEW%P(:,:)
          ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:) = YLGVV(ID)%VIEW%P(:,:)
        ENDDO
      ENDDO
  ENDIF

  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
      ID = JLEV + (JFLD -1) * KFLEVG
      IF (JFLD .EQ. 1) IVSETUV(JLEV) = YLGVU(ID)%IVSET
      IF (IVSETUV(JLEV) .NE. YLGVU(ID)%IVSET)  CALL ABOR1("[DIR_TRANS_FIELD_API] ivsetuv inconsistent with ylgvu%ivset")
    ENDDO
  ENDDO
ELSE
  ! No vector field provided
  IUVG = 0
  ZPGPUV=>NULL()
  ZPSPVOR=>NULL()
  ZPSPDIV=>NULL()
ENDIF

! 2. scalar fields transformation

! Preliminary checks
IF (PRESENT(YDFSPSCALAR) .NEQV. PRESENT(YDFSCALAR))  CALL ABOR1("[DIR_TRANS_FIELD_API]  YDFSPSCALAR and YDFSCALAR must be provided together")

! Do we have scalar fields?
IF (PRESENT(YDFSPSCALAR)) THEN
  IF ((SIZE(YDFSPSCALAR)/= SIZE(YDFSCALAR)))  CALL ABOR1("[DIR_TRANS_FIELD_API] Inconsistent size for YDFSPSCALAR and YDFSCALAR")

  ! Convert list of spectral scalar fields of any dimension into a list of 2d fields
  YLGVSCALAR = LG(YDFSCALAR, LDACC,.TRUE.)
  YLSPVSCALAR = LS(YDFSPSCALAR, LDACC,.FALSE.)

  IFLDSPSC = SIZE(YLSPVSCALAR)
  IFLDXG = SIZE(YLGVSCALAR)

  ! count the number of fields present on the processor
  IFLDXL = 0
  DO JFLD = 1, SIZE(YLGVSCALAR)
    IF (ASSOCIATED(YLSPVSCALAR(JFLD)%VIEW%P)) IFLDXL = IFLDXL + 1
  END DO

  ISCDIM = 1

   ! Allocate temporary scalar field array in spectral space
  ALLOCATE(ZPSPSC2(IFLDXL,KSPEC))

  ! Allocate temporary scalar field array in grid space
  ALLOCATE(ZPGP2(KPROMA,IFLDXG * ISCDIM,KGPBLKS))

  ! allocate 'b-set' for scalar fields
  ALLOCATE(IVSETSC2(IFLDXG))

  ! temporary copies on gpu
   if ( LDACC ) THEN
    !$ACC ENTER DATA CREATE(ZPSPSC2,ZPGP2)
  ENDIF


  ! Copy list of scalar fields into temporary arrays (2d copy thanks to field_view)

  IF (LDACC) THEN
    !$ACC KERNELS PRESENT(ZPGP2,YLGVSCALAR(:)%VIEW%P) FIRSTPRIVATE(IFLDXG)  PRIVATE(JFLD)
    DO JFLD=1, IFLDXG
      ZPGP2(:,JFLD,:) = YLGVSCALAR(JFLD)%VIEW%P(:,:)
    ENDDO
    !$ACC END KERNELS
  ELSE
    DO JFLD=1, IFLDXG
      ZPGP2(:,JFLD,:) = YLGVSCALAR(JFLD)%VIEW%P(:,:)
    ENDDO
  ENDIF

  DO JFLD=1, IFLDXG
    IVSETSC2(JFLD) = YLGVSCALAR(JFLD)%IVSET
  ENDDO

ELSE
  !No scalar field provided
  ISPUV = 0
  IFLDXG = 0
  ZPGP2=>NULL()
  ZPSPSC2=>NULL()
ENDIF

! 3. CALL DIR_TRANS using the regular interface and the temporary arrays

! We have to perform separated calls for nvfortran
IF (ASSOCIATED(ZPGP2) .AND. ASSOCIATED(ZPGPUV)) THEN
	CALL DIR_TRANS(PSPVOR = ZPSPVOR,PSPDIV = ZPSPDIV,PGPUV = ZPGPUV,KVSETUV = IVSETUV, &
	             & PSPSC2 = ZPSPSC2,PGP2 = ZPGP2, KVSETSC2 = IVSETSC2, &
	             & KPROMA = KPROMA)
ELSE IF (ASSOCIATED(ZPGP2)) THEN
	CALL DIR_TRANS(PSPSC2 = ZPSPSC2,PGP2 = ZPGP2, KVSETSC2 = IVSETSC2, &
	             & KPROMA = KPROMA)
ELSE IF (ASSOCIATED(ZPGPUV)) THEN
	CALL DIR_TRANS(PSPVOR = ZPSPVOR,PSPDIV = ZPSPDIV,PGPUV = ZPGPUV,KVSETUV = IVSETUV, &
	             & KPROMA = KPROMA)
ENDIF

! 4. Copy back temporary array data into spectral fields

! copy spectral vorticity and divergence
IF (IUVG>0) THEN
  IF (LDACC) THEN
    !$ACC KERNELS PRESENT(ZPSPVOR, ZPSPDIV,YLSPVVOR(:)%VIEW%P, YLSPVDIV(:)%VIEW%P) FIRSTPRIVATE(IFLDSPVOR) PRIVATE(JFLD)
    DO JFLD=1,IFLDSPVOR
      YLSPVVOR(JFLD)%VIEW%P(:) = ZPSPVOR(JFLD,:)
      YLSPVDIV(JFLD)%VIEW%P(:) = ZPSPDIV(JFLD,:)
    ENDDO
    !$ACC END KERNELS
  ELSE
    DO JFLD=1,IFLDSPVOR
      YLSPVVOR(JFLD)%VIEW%P(:) = ZPSPVOR(JFLD,:)
      YLSPVDIV(JFLD)%VIEW%P(:) = ZPSPDIV(JFLD,:)
    ENDDO
  ENDIF
ENDIF

 ! copy spectral scalar fields

 IF (IFLDSPSC > 0) THEN
   ID = 1
   IF (LDACC) THEN
      !$ACC KERNELS PRESENT(ZPSPSC2,YLSPVSCALAR(:)%VIEW%P) FIRSTPRIVATE(IFLDSPSC) PRIVATE(JFLD,ID)
      DO JFLD = 1, IFLDSPSC
          IF (ASSOCIATED(YLSPVSCALAR(JFLD)%VIEW%P)) THEN
            YLSPVSCALAR(JFLD)%VIEW%P(:) = ZPSPSC2(ID,:)
            ID = ID + 1
          ENDIF
      ENDDO
      !$ACC END KERNELS
   ELSE
      DO JFLD = 1, IFLDSPSC
        IF (ASSOCIATED(YLSPVSCALAR(JFLD)%VIEW%P)) THEN
          YLSPVSCALAR(JFLD)%VIEW%P(:) = ZPSPSC2(ID,:)
          ID = ID + 1
        ENDIF
      ENDDO
   ENDIF
ENDIF
! 5. Final cleanup

! delete temporary arrays

 IF ( LDACC ) THEN
  IF (ASSOCIATED(ZPSPVOR))  THEN
      !$ACC EXIT DATA DELETE(ZPSPVOR)
  ENDIF
  IF (ASSOCIATED(ZPSPDIV)) THEN
   !$ACC EXIT DATA DELETE(ZPSPDIV)
  ENDIF
  IF (ASSOCIATED(ZPGPUV)) THEN
     !$ACC EXIT DATA DELETE(ZPGPUV)
  ENDIF
  IF (ASSOCIATED(ZPSPSC2)) THEN
   !$ACC EXIT DATA DELETE(ZPSPSC2)
  ENDIF
  IF (ASSOCIATED(ZPGP2)) THEN
     !$ACC EXIT DATA DELETE(ZPGP2)
  ENDIF
ENDIF

IF (ASSOCIATED(ZPSPVOR)) DEALLOCATE(ZPSPVOR)
IF (ASSOCIATED(ZPSPDIV)) DEALLOCATE(ZPSPDIV)
IF (ASSOCIATED(ZPSPSC2)) DEALLOCATE(ZPSPSC2)
IF (ASSOCIATED(ZPGPUV))  DEALLOCATE(ZPGPUV)
IF (ASSOCIATED(ZPGP2))   DEALLOCATE(ZPGP2)
IF (ALLOCATED(IVSETUV))  DEALLOCATE(IVSETUV)
IF (ALLOCATED(IVSETSC2)) DEALLOCATE(IVSETSC2)

! delete FIELD_VIEWS
IF (ALLOCATED(YLSPVVOR))    DEALLOCATE(YLSPVVOR)
IF (ALLOCATED(YLSPVDIV))    DEALLOCATE(YLSPVDIV)
IF (ALLOCATED(YLSPVSCALAR)) DEALLOCATE(YLSPVSCALAR)
IF (ALLOCATED(YLGVU))       DEALLOCATE(YLGVU)
IF (ALLOCATED(YLGVV))       DEALLOCATE(YLGVV)
IF (ALLOCATED(YLGVSCALAR))  DEALLOCATE(YLGVSCALAR)

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE DIR_TRANS_FIELD_API


