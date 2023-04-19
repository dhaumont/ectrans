MODULE ELTDIR_MOD
CONTAINS
SUBROUTINE ELTDIR(KF_FS,KF_UV,KF_SCALARS,KLED2,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPM_DIM         ,ONLY : R
USE TPM_DISTR       ,ONLY : D
USE TPMALD_DIM      ,ONLY : RALD
USE TPM_GEN          ,ONLY : LALLOPERM2
USE ELTDATA_MOD      ,ONLY : ZFFT_PERM, ZVODI_PERM



USE EPRFI2_MOD      ,ONLY : EPRFI2
USE ELEDIR_MOD      ,ONLY : ELEDIR
USE EUVTVD_MOD
USE EUPDSP_MOD      ,ONLY : EUPDSP
USE EUVTVD_COMM_MOD 
USE EXTPER_MOD      ,ONLY : EXTPER

USE TPM_DISTR       ,ONLY : D_NUMP
USE TPM_DIM         ,ONLY : R_NDGL

!
!**** *ELTDIR* - Control of Direct Legendre transform step

!     Purpose.
!     --------
!        Tranform from Fourier space to spectral space, compute
!        vorticity and divergence.

!**   Interface.
!     ----------
!        *CALL* *ELTDIR(...)*

!        Explicit arguments :
!        --------------------  IM     - zonal wavenumber
!                              JM  - local zonal wavenumber

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!         EPRFI2      - prepares the Fourier work arrays for model variables
!         ELEDIR      - direct Legendre transform
!         EUVTVD      -
!         EUPDSP      - updating of spectral arrays (fields)
!         EUVTVD_COMM -
!         EXTPER      -


!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-24
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
!        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified 94-04-06 R. El khatib Full-POS implementation
!        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
!                             instead of u,v->vor,div
!        MPP Group : 95-10-01 Support for Distributed Memory version
!        K. YESSAD (AUGUST 1996):
!               - Legendre transforms for transmission coefficients.
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!            01-03-14 G. Radnoti aladin version
!     01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2

REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPMEANU(:)
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPMEANV(:)

REAL(KIND=JPRB), POINTER    :: ZFFT(:,:,:)
INTEGER(KIND=JPIM) :: IINDEX(2*KF_FS), JF, JDIM
INTEGER(KIND=JPIM) :: IM
INTEGER(KIND=JPIM) :: JM
INTEGER(KIND=JPIM) :: IUS,IUE,IVS,IVE,IVORS,IVORE,IDIVS,IDIVE,IFC
REAL(KIND=JPRB), POINTER    :: ZVODI(:,:,:)
INTEGER(KIND=JPIM) :: JGL, IJR, IJI

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',0,ZHOOK_HANDLE)


!*     1.    PREPARE WORK ARRAYS.
!            --------------------

IF (ALLOCATED (ZFFT_PERM)) THEN
  IF ((UBOUND (ZFFT_PERM, 1) /= RALD%NDGLSUR+R%NNOEXTZG) &
& .OR. (UBOUND (ZFFT_PERM, 2) /= D%NUMP) &
& .OR. (UBOUND (ZFFT_PERM, 3) < KLED2)) THEN
    !$acc exit data delete (ZFFT_PERM)
    DEALLOCATE (ZFFT_PERM)
  ENDIF
ENDIF

IF (.NOT. ALLOCATED (ZFFT_PERM)) THEN
  ALLOCATE (ZFFT_PERM (RALD%NDGLSUR+R%NNOEXTZG,D%NUMP,KLED2))
  !$acc enter data create (ZFFT_PERM)
ENDIF


IF (ALLOCATED (ZVODI_PERM)) THEN
  IF ((UBOUND (ZVODI_PERM, 1) /= RALD%NDGLSUR+R%NNOEXTZG) &
& .OR. (UBOUND (ZVODI_PERM, 2) /= D%NUMP) &
& .OR. (UBOUND (ZVODI_PERM, 3) < MAX(4*KF_UV,1))) THEN
    !$acc exit data delete (ZVODI_PERM)
    DEALLOCATE (ZVODI_PERM)
  ENDIF
ENDIF

IF (.NOT. ALLOCATED (ZVODI_PERM)) THEN
  ALLOCATE (ZVODI_PERM (RALD%NDGLSUR+R%NNOEXTZG,D%NUMP,MAX(4*KF_UV,1)))
  !$acc enter data create (ZVODI_PERM)
ENDIF

ZFFT => ZFFT_PERM (:,:,1:KLED2)
ZVODI => ZVODI_PERM (:,:,1:MAX(4*KF_UV,1))

!$acc kernels present (ZVODI, ZFFT)
ZVODI = 0._JPRB
ZFFT = 0._JPRB
!$acc end kernels


IFC = 2 * KF_FS

CALL EPRFI2(KF_FS,ZFFT)

!*     2.    PERIODICIZATION IN Y DIRECTION
!            ------------------------------

IF(R%NNOEXTZG>0) THEN
  CALL ABOR1 ('ELTDIR: BIPERIODICIZATION NOT SUPPORTED')
ENDIF

!*     3.    DIRECT LEGENDRE TRANSFORM.
!            --------------------------

CALL ELEDIR(IFC,KLED2,ZFFT)


!*     4.    COMPUTE VORTICITY AND DIVERGENCE.
!            ---------------------------------

IF( KF_UV > 0 ) THEN
  IUS = 1
  IUE = 2*KF_UV
  IVS = 2*KF_UV+1
  IVE = 4*KF_UV
  IVORS = 1
  IVORE = 2*KF_UV
  IDIVS = 2*KF_UV+1
  IDIVE = 4*KF_UV

  CALL EUVTVD(KF_UV,ZFFT(:,:,IUS:IUE),ZFFT(:,:,IVS:IVE),&
   & ZVODI(:,:,IVORS:IVORE),ZVODI(:,:,IDIVS:IDIVE))

!*     5.    COMMUNICATION OF MEAN WIND
!            --------------------------


  DO JM=1,D%NUMP
    IM = D%MYMS(JM)

    CALL EUVTVD_COMM(IM,JM,KF_UV,KFLDPTRUV,ZFFT(:,JM,IUS:IUE), &
     & ZFFT(:,JM,IVS:IVE),ZVODI(:,JM,IVORS:IVORE),ZVODI(:,JM,IDIVS:IDIVE), &
     & PSPMEANU,PSPMEANV)

  ENDDO

ENDIF


!*       6.    UPDATE SPECTRAL ARRAYS.
!              -----------------------

CALL EUPDSP(KF_UV,KF_SCALARS,ZFFT,ZVODI, &
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)

IF (.NOT. LALLOPERM2) THEN
  !$acc exit data delete (ZFFT_PERM, ZVODI_PERM)
  DEALLOCATE (ZFFT_PERM)
  DEALLOCATE (ZVODI_PERM)
ENDIF

IF (LHOOK) CALL DR_HOOK('ELTDIR_MOD:ELTDIR',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE ELTDIR
END MODULE ELTDIR_MOD
