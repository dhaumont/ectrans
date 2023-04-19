MODULE EUPDSPB_MOD
CONTAINS
SUBROUTINE EUPDSPB(KFIELD,POA,PSPEC,KFLDPTR)

!**** *EUPDSPB* - Update spectral arrays after direct Legendre transform

!     Purpose.
!     --------
!        To update spectral arrays for a fixed zonal wave-number
!         from values in POA.

!**   Interface.
!     ----------
!        CALL EUPDSPB(....)

!        Explicit arguments :  
!        --------------------  KFIELD  - number of fields
!                              POA - work array
!                              PSPEC - spectral array

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-02-02
!        D. Giard : 93-03-19 truncations NSMAX and NTMAX (see NOTE)
!        R. El Khatib : 94-08-02 Replace number of fields by indexes of the
!                       first and last field
!        L. Isaksen : 95-06-06 Reordering of spectral arrays
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TPMALD_DISTR    ,ONLY : DALD, DALD_NESM0, DALD_NCPL2M
USE TPM_DISTR       ,ONLY : D, D_MYMS
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
REAL(KIND=JPRB)   ,INTENT(IN)  :: POA(:,:,:)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PSPEC(:,:)
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

INTEGER(KIND=JPIM) :: II, INM, IR, JFLD, JN,IFLD, JM, IM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------

!*       1.    UPDATE SPECTRAL FIELDS.
!              -----------------------
IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',0,ZHOOK_HANDLE)

!$ACC data present (POA, PSPEC)

IF(PRESENT(KFLDPTR)) THEN
  
  ! TODO

  DO JM = 1, D%NUMP
    IM = D%MYMS(JM)
  
    DO JN=1,DALD%NCPL2M(IM),2
      INM=DALD%NESM0(IM)+(JN-1)*2
      DO JFLD=1,KFIELD
        IR= 2*JFLD-1
        II=IR+1
        IFLD = KFLDPTR(JFLD)
        PSPEC(IFLD,INM)    =POA(JN  ,JM,IR)
        PSPEC(IFLD,INM+1)  =POA(JN+1,JM,IR)
        PSPEC(IFLD,INM+2)  =POA(JN  ,JM,II)
        PSPEC(IFLD,INM+3)  =POA(JN+1,JM,II)
      ENDDO
    ENDDO

  ENDDO

ELSE
  !$ACC parallel loop
  DO JM = 1, D%NUMP
    IM = D_MYMS(JM)
  
    !$ACC loop
    DO JN=1,DALD_NCPL2M(IM),2
      INM=DALD_NESM0(IM)+(JN-1)*2
      !$ACC loop
      DO JFLD=1,KFIELD
        IR= 2*JFLD-1
        II=IR+1
        PSPEC(JFLD,INM)    =POA(JN  ,JM,IR)
        PSPEC(JFLD,INM+1)  =POA(JN+1,JM,IR)
        PSPEC(JFLD,INM+2)  =POA(JN  ,JM,II)
        PSPEC(JFLD,INM+3)  =POA(JN+1,JM,II)
      ENDDO
    ENDDO
  
  ENDDO

ENDIF

!$ACC end data

IF (LHOOK) CALL DR_HOOK('EUPDSPB_MOD:EUPDSPB',1,ZHOOK_HANDLE)

END SUBROUTINE EUPDSPB
END MODULE EUPDSPB_MOD
