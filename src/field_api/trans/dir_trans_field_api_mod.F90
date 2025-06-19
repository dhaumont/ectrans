MODULE DIR_TRANS_FIELD_API_MOD
USE FIELD_API_ECTRANS_MOD
USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE


CONTAINS

SUBROUTINE DIR_TRANS_FIELD_API(YDFSPVOR,YDFSPDIV,YDFSPSCALAR, &
                             & YDFU, YDFV, YDFVOR,YDFDIV,YDFSCALAR, &
                             & KSPEC, KPROMA, KGPBLKS, KGPTOT, KFLEVG, KFLEVL,&
                             & LDACC, LDVERBOSE)


TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPVOR(:), YDFSPDIV(:)        !Spectral vector fields : vorticity and divergence fields (in)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSPSCALAR(:)                !Spectral scalar fields (in)

TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFU(:),YDFV(:)                 !Grid vector fields     (out)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFVOR(:),YDFDIV(:)             !Grid vector fields :vorticity and divergence     (out)
TYPE(FIELD_BASIC_PTR),INTENT(IN), OPTIONAL  :: YDFSCALAR(:)                  !Grid scalar fields     (out)

INTEGER(KIND=JPIM), INTENT(IN) ::KSPEC
INTEGER(KIND=JPIM), INTENT(IN) ::KPROMA
INTEGER(KIND=JPIM), INTENT(IN) ::KGPBLKS
INTEGER(KIND=JPIM), INTENT(IN) ::KGPTOT
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVG
INTEGER(KIND=JPIM), INTENT(IN) :: KFLEVL
LOGICAL, INTENT(IN), OPTIONAL  :: LDACC
LOGICAL, INTENT(IN), OPTIONAL  :: LDVERBOSE

! Local variables

! Temporary arrays:

TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE  :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVVOR(:),YLGVDIV(:)
TYPE(GRID_VIEW), ALLOCATABLE  :: YLGVSCALAR(:)

REAL(KIND=JPRB),ALLOCATABLE :: ZPSPVOR(:,:),ZPSPDIV(:,:)          ! Spectral vector fields (in)
REAL(KIND=JPRB),ALLOCATABLE :: ZPSPSC2(:,:)                       ! Spectral surface scalar fields(in)
REAL(KIND=JPRB),ALLOCATABLE :: ZPGPUV(:,:,:,:)                    ! Grid vector fields (out)
REAL(KIND=JPRB),ALLOCATABLE :: ZPGP2(:,:,:)                       ! Grid surface scalar fields (out)

INTEGER(KIND=JPIM)          :: ISPUV                              ! Number of input spectral vector fields
INTEGER(KIND=JPIM)          :: IFLDXG
INTEGER(KIND=JPIM)          :: IFLDXL
INTEGER(KIND=JPIM)          :: IFLDXGUV
INTEGER(KIND=JPIM)          :: IFLDXLUV
INTEGER(KIND=JPIM)          :: IUVG                               ! Number of output vector fields
INTEGER(KIND=JPIM)          :: ISCDIM                             ! Size of output scalar fields array
INTEGER(KIND=JPIM)          :: IUVDIM                             ! Size of output vector fields array
INTEGER(KIND=JPIM)          :: ID,IOFFSET,JLEV
INTEGER(KIND=JPIM)          :: IEND


INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETUV(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IVSETSC2(:)

INTEGER(KIND=JPIM)          :: JFLD                                 ! field counter
LOGICAL          :: LLSCDERS                                        ! indicating if derivatives of scalar variables are req.
LOGICAL          :: LLVORGP                                         ! indicating if grid-point vorticity is req.
LOGICAL          :: LLDIVGP                                         ! indicating if grid-point divergence is req.
LOGICAL          :: LLUVDER                                         ! indicating if E-W derivatives of u and v are req.
LOGICAL          :: LLVERBOSE                                       ! indicating if verbose output is req.
#include "dir_trans.h"
#include "abor1.intfb.h"

LLVERBOSE = .FALSE.
IF (PRESENT(LDVERBOSE))  LLVERBOSE = LDVERBOSE

! 1. VECTOR FIELDS TRANSFORMATION

! Check if all provided vector field information is consistent
IF (PRESENT(YDFU) .NEQV. PRESENT(YDFV)) CALL ABOR1("[ECTRANS_API] IYDU/YDFV")
IF (PRESENT(YDFU) .NEQV. PRESENT(YDFSPVOR)) CALL ABOR1("[ECTRANS_API] IYDU/YDFVOR")
IF (PRESENT(YDFU) .NEQV. PRESENT(YDFSPDIV)) CALL ABOR1("[ECTRANS_API] IYDU/YDFDIV")

IUVDIM = 0
! Do we have vector fields?

IF (PRESENT(YDFU)) THEN

  IF ((SIZE(YDFU)/= SIZE(YDFV)).OR.(SIZE(YDFU)/= SIZE(YDFSPDIV)).OR.(SIZE(YDFU)/= SIZE(YDFSPVOR))) THEN
     CALL ABOR1("[ECTRANS_API] INVALID LIST SIZES:")
  ENDIF
  YLSPVVOR = LS (YDFSPVOR, LDACC)
  YLSPVDIV = LS (YDFSPDIV, LDACC)

  YLGVU = LG (YDFU, LDACC)
  YLGVV = LG (YDFV, LDACC)
  IF ((SIZE (YLGVU) /= SIZE (YLGVV)) .OR. (SIZE (YLSPVVOR) /= SIZE (YLSPVDIV))) THEN
     CALL ABOR1("[ECTRANS_API] INCONSISTENT NUMBER OF FIELD_VIEW FOR VECTORS:")
  ENDIF
   IF (((SIZE (YLGVU) / SIZE (YDFU)) /= KFLEVG) .OR. ((SIZE (YLSPVVOR) / SIZE (YDFSPVOR)) /= KFLEVL)) THEN
     CALL ABOR1("[ECTRANS_API] INCONSISTENT KFLEVG OR KFLEVL")
  ENDIF

  IUVG = SIZE(YDFU)           ! Number of output  vector fields
  ISPUV = SIZE(YDFSPVOR)

  IUVDIM = 2

  IF (PRESENT(YDFDIV)) THEN
     LLDIVGP = .TRUE.
     IUVDIM = 5
     YLGVDIV = LG (YDFDIV, LDACC)
  ENDIF

  IF (PRESENT(YDFVOR)) THEN
     LLVORGP = .TRUE.
     IUVDIM = 6
     YLGVVOR = LG (YDFVOR, LDACC)
  ENDIF
   
  ! Allocate vector field input in spectral space
  ALLOCATE(ZPSPVOR(SIZE(YLSPVVOR),KSPEC))
  ALLOCATE(ZPSPDIV(SIZE(YLSPVDIV),KSPEC))

  ! Allocate vector field output in grid space
  ALLOCATE(ZPGPUV(KPROMA,KFLEVG, IUVG * IUVDIM,KGPBLKS))
  ALLOCATE(IVSETUV(KFLEVG))

  ! Copy from fields to temporary arrays (1D copy thanks to FIELD VIEW)

  DO JFLD=1,SIZE(YLSPVVOR)
     ZPSPVOR(JFLD,:) = YLSPVVOR(JFLD)%VIEW%P(:)
     ZPSPDIV(JFLD,:) = YLSPVDIV(JFLD)%VIEW%P(:)
  ENDDO


ELSE
  ! YDFU is not provided, we do not have to compute the corresponding vector output
  ISPUV = 0
  ALLOCATE(ZPGPUV(0,ISPUV,2,0),ZPSPVOR(0,0),ZPSPDIV(0,0))
ENDIF

! 2. SCALAR FIELDS TRANSFORMATION

! Check if all provided scalar field information is consistent
IF (PRESENT(YDFSPSCALAR) .NEQV. PRESENT(YDFSCALAR))  CALL ABOR1("[ECTRANS_API] GRID/SPEC")

IF (PRESENT(YDFSPSCALAR)) THEN

  IF ((SIZE(YDFSPSCALAR)/= SIZE(YDFSCALAR)))  CALL ABOR1("[ECTRANS_API] INCONSISTENT NUMBER OF FIELD_VIEW FOR YDFSCALAR")

  YLGVSCALAR = LG (YDFSCALAR, LDACC)
  YLSPVSCALAR = LS (YDFSPSCALAR, LDACC)

  IFLDXG = SIZE(YLGVSCALAR) ! Number of output scalar fields in grid space

  IFLDXL = 0  ! Number of input scalar fields in spectral space
  DO JFLD = 1, SIZE(YLGVSCALAR) ! Number of output scalar fields in grid space
  IF (ASSOCIATED(YLSPVSCALAR(JFLD)%VIEW%P)) IFLDXL = IFLDXL + 1
  END DO 

  ISCDIM = 1

   ! Allocate scalar field input in spectral space
  ALLOCATE(ZPSPSC2(IFLDXL,KSPEC))

  ! Allocate scalar field output in grid space
  ALLOCATE(ZPGP2(KPROMA,IFLDXG * ISCDIM,KGPBLKS))
  ALLOCATE(IVSETSC2(IFLDXG))

 ! Copy scalar spectral fields to temporary arrays (1D copy thanks to FIELD VIEW)
  ID = 1
  DO JFLD = 1, SIZE(YLSPVSCALAR) ! Number of output scalar fields in grid space
     IF (ASSOCIATED(YLSPVSCALAR(JFLD)%VIEW%P)) THEN
         ZPSPSC2(ID,:) = YLSPVSCALAR(JFLD)%VIEW%P(:)
         ID = ID + 1 
     ENDIF
  ENDDO

 DO JFLD=1, IFLDXG
    IVSETSC2(JFLD) = YLGVSCALAR(JFLD)%IVSET
 ENDDO


DO JFLD=1,IUVG
  DO JLEV=1,KFLEVG
     ID = JLEV + (JFLD -1) * KFLEVG
     IF (JFLD .EQ. 1) IVSETUV(JLEV) = YLGVU(ID)%IVSET
     IF (IVSETUV(JLEV) .NE. YLGVU(ID)%IVSET)  CALL ABOR1("[ECTRANS_API] IVSETUV INCONSISTENT WITH YLGVU%IVSET")
     IF (IVSETUV(JLEV) .NE. YLGVV(ID)%IVSET)  CALL ABOR1("[ECTRANS_API] IVSETUV INCONSISTENT WITH YLGVV%IVSET")
  ENDDO
ENDDO
!
ELSE
  IFLDXG = 0
  ALLOCATE(ZPGP2(0,IFLDXG,0),ZPSPSC2(0,1))
ENDIF


! 3. CALL DIR_TRANS

	CALL DIR_TRANS(PSPVOR = ZPSPVOR,PSPDIV = ZPSPDIV,PGPUV = ZPGPUV,KVSETUV = IVSETUV, &
	             & PSPSC2 = ZPSPSC2,PGP2 = ZPGP2, KVSETSC2 = IVSETSC2, &
	             & KPROMA = KPROMA)
! 4. Copy back data to fields

! Remove garbage at the end of arrays
  IEND = KGPTOT - KPROMA * (KGPBLKS - 1)
  ZPGPUV (IEND+1:, :, :, KGPBLKS) = 0
  ZPGP2 (IEND+1:, :, KGPBLKS) = 0

  ! Copy vector fields back from temporary vector arrays
IOFFSET = 0

IF (LLVORGP) THEN
  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
     ID = JLEV + (JFLD -1) * KFLEVG
     YLGVVOR(ID)%VIEW%P(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
    ENDDO
  ENDDO
  IOFFSET = IOFFSET + 1
ENDIF

IF (LLDIVGP) THEN
  DO JFLD=1,IUVG
    DO JLEV=1,KFLEVG
     ID = JLEV + (JFLD -1) * KFLEVG
     YLGVDIV(ID)%VIEW%P(:,:) = ZPGPUV(:, JLEV,JFLD+IOFFSET*IUVG,:)
    ENDDO
  ENDDO
  IOFFSET = IOFFSET + 1
ENDIF

DO JFLD=1,IUVG
  DO JLEV=1,KFLEVG
     ID = JLEV + (JFLD -1) * KFLEVG
     YLGVU(ID)%VIEW%P(:,:) =  ZPGPUV(:,JLEV,JFLD+IOFFSET*IUVG,:)
     YLGVV(ID)%VIEW%P(:,:) =  ZPGPUV(:,JLEV,JFLD+(IOFFSET+1)*IUVG,:)
  ENDDO
ENDDO

! Copy scalar fields back from temporary scalar arrays
DO JFLD=1, IFLDXG
    YLGVSCALAR(JFLD)%VIEW%P(:,:) = ZPGP2(:,JFLD,:)
ENDDO


DEALLOCATE(ZPSPVOR)
DEALLOCATE(ZPSPDIV)
DEALLOCATE(ZPSPSC2)
DEALLOCATE(ZPGPUV)
DEALLOCATE(ZPGP2)
DEALLOCATE(IVSETUV)
DEALLOCATE(IVSETSC2)

END SUBROUTINE DIR_TRANS_FIELD_API

END MODULE DIR_TRANS_FIELD_API_MOD

