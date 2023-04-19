MODULE EPRFI1B_MOD
CONTAINS
SUBROUTINE EPRFI1B(PFFT,PSPEC,KFIELDS,KFLDPTR)

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!USE TPM_DIM
USE TPM_DISTR
USE TPMALD_DISTR    ,ONLY : DALD, DALD_NESM0, DALD_NCPL2M
!
!**** *PRFI1* - Prepare spectral fields for inverse Legendre transform

!     Purpose.
!     --------
!        To extract the spectral fields for a specific zonal wavenumber
!        and put them in an order suitable for the inverse Legendre           .
!        tranforms.The ordering is from NSMAX to KM for better conditioning.
!        Elements 1,2 and NLCM(KM)+1 are zeroed in preparation for computing
!        u,v and derivatives in spectral space.

!**   Interface.
!     ----------
!        *CALL* *PRFI1B(...)*

!        Explicit arguments :  KM     - zonal wavenumber
!        ------------------    PIA    - spectral components for transform
!                              PSPEC  - spectral array
!                              KFIELDS  - number of fields

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From PRFI1B in IFS CY22R1
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)  :: PFFT(:,:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, J, JFLD, ILCM, IOFF,IFLD
INTEGER(KIND=JPIM) :: IM, JM, MAX_NCPL2M
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
!              --------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',0,ZHOOK_HANDLE)

!$acc data present (PFFT, PSPEC)

!$acc kernels default(none)
PFFT = 0._JPRB
!$acc end kernels

IF(PRESENT(KFLDPTR)) THEN
  ! TODO 
  DO JFLD=1,KFIELDS
    IR = 2*(JFLD-1)+1
    II = IR+1
    IFLD = KFLDPTR(JFLD)
    DO JM = 1, D%NUMP
      IM   = D%MYMS(JM)
      ILCM = DALD%NCPL2M(IM)
      IOFF = DALD%NESM0(IM)
      DO J=1,ILCM,2
        INM = IOFF+(J-1)*2
        PFFT(J  ,JM,IR) = PSPEC(IFLD,INM  )
        PFFT(J+1,JM,IR) = PSPEC(IFLD,INM+1)
        PFFT(J  ,JM,II) = PSPEC(IFLD,INM+2)
        PFFT(J+1,JM,II) = PSPEC(IFLD,INM+3)
      ENDDO
    ENDDO
  ENDDO
ELSE
  MAX_NCPL2M = MAXVAL (DALD_NCPL2M)
  !$ACC parallel loop collapse(3) &
  !$ACC& present(D_MYMS,DALD_NCPL2M,DALD_NESM0) &
  !$ACC& present(PFFT,PSPEC) &
  !$ACC& private(IR,II,IM,ILCM,IOFF,INM) default(none)
  DO JFLD=1,KFIELDS
    DO JM = 1, D_NUMP
      DO J=1,MAX_NCPL2M,2
       IR = 2*(JFLD-1)+1
       II = IR+1
       IM   = D_MYMS(JM)
       ILCM = DALD_NCPL2M(IM)
       if (J > ILCM) CYCLE
       IOFF = DALD_NESM0(IM)
       INM = IOFF+(J-1)*2
       PFFT(J  ,JM,IR) = PSPEC(JFLD,INM  )
       PFFT(J+1,JM,IR) = PSPEC(JFLD,INM+1)
       PFFT(J  ,JM,II) = PSPEC(JFLD,INM+2)
       PFFT(J+1,JM,II) = PSPEC(JFLD,INM+3)
      ENDDO
    ENDDO
  ENDDO
ENDIF

!$acc end data


IF (LHOOK) CALL DR_HOOK('EPRFI1B_MOD:EPRFI1B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI1B
END MODULE EPRFI1B_MOD
