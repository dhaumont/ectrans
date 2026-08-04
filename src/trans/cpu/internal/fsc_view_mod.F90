! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FSC_VIEW_MOD
CONTAINS
SUBROUTINE FSC_VIEW(KGL,KF_UV,KF_SCALARS,KF_SCDERS,&
 & YLSPVUV,YLSPVSCALAR,&
 & YLSPVUV_EW,YLSPVSCALAR_NS,YLSPVSCALAR_EW)

!**** *FSC_VIEW - Division by a*cos(theta), east-west derivatives

!     Purpose.
!     --------
!        In Fourier space divide u and v and all north-south
!        derivatives by a*cos(theta). Also compute east-west derivatives
!        of u,v,thermodynamic, passiv scalar variables and surface
!        pressure.

!**   Interface.
!     ----------
!        CALL FSC_VIEW(..)
!        Explicit arguments :  PUV     - u and v
!        --------------------  PSCALAR - scalar valued varaibles
!                              PNSDERS - N-S derivative of S.V.V.
!                              PEWDERS - E-W derivative of S.V.V.
!                              PUVDERS - E-W derivative of u and v
!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03 (From SC2FSC_VIEW)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_TRANS       ,ONLY : LUVDER, LATLON
USE TPM_DISTR       ,ONLY : D, MYSETW
USE TPM_FIELDS      ,ONLY : F
USE TPM_GEOMETRY    ,ONLY : G
USE TPM_FLT                ,ONLY: S
!

IMPLICIT NONE

TYPE SPEC_VIEW
  ! Spectral field view
  REAL(KIND=JPRB),POINTER  :: P(:) => NULL()
  CHARACTER(LEN=:),POINTER :: NAME => NULL()
END TYPE


INTEGER(KIND=JPIM) , INTENT(IN) :: KGL,KF_UV,KF_SCALARS,KF_SCDERS
TYPE(SPEC_VIEW) :: YLSPVUV(:),YLSPVSCALAR(:)
TYPE(SPEC_VIEW) :: YLSPVUV_EW(:),YLSPVSCALAR_NS(:), YLSPVSCALAR_EW(:)

REAL(KIND=JPRB) :: ZACHTE,ZMUL, ZACHTE2, ZSHIFT, ZPI
REAL(KIND=JPRB) :: ZAMP, ZPHASE
INTEGER(KIND=JPIM) :: IMEN,ISTAGTF


INTEGER(KIND=JPIM) :: JLON,JF,IGLG,II,IR,JM

!     ------------------------------------------------------------------

IGLG    = D%NPTRLS(MYSETW)+KGL-1
ZACHTE  = REAL(F%RACTHE(IGLG),JPRB)
IMEN    = G%NMEN(IGLG)
ISTAGTF = D%NSTAGTF(KGL)
ZACHTE2 = REAL(F%RACTHE(IGLG),JPRB)

IF( LATLON.AND.S%LDLL ) THEN
  ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
  ZACHTE2 = 1._JPRB
  ZACHTE  = REAL(F%RACTHE2(IGLG),JPRB)
  
  ! apply shift for (even) lat-lon output grid
  IF( S%LSHIFTLL ) THEN
    ZSHIFT = ZPI/REAL(G%NLOEN(IGLG),JPRB)

    DO JF=1,KF_SCALARS
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        
        ! calculate amplitude and add phase shift then reconstruct A,B
        ZAMP = SQRT(YLSPVSCALAR(JF)%P(IR)**2 + YLSPVSCALAR(JF)%P(II)**2)
        ZPHASE = ATAN2(YLSPVSCALAR(JF)%P(II),YLSPVSCALAR(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
        
        YLSPVSCALAR(JF)%P(IR) = ZAMP*COS(ZPHASE)
        YLSPVSCALAR(JF)%P(II) = ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
    IF(KF_SCDERS > 0)THEN
      DO JF=1,KF_SCALARS
        DO JM=0,IMEN
          IR = ISTAGTF+2*JM+1
          II = IR+1          
          ! calculate amplitude and phase shift and reconstruct A,B
          ZAMP = SQRT(YLSPVSCALAR_NS(JF)%P(IR)**2 + YLSPVSCALAR_NS(JF)%P(II)**2)
          ZPHASE = ATAN2(YLSPVSCALAR_NS(JF)%P(II),YLSPVSCALAR_NS(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
          YLSPVSCALAR_NS(JF)%P(IR) = ZAMP*COS(ZPHASE)
          YLSPVSCALAR_NS(JF)%P(II) = ZAMP*SIN(ZPHASE)
        ENDDO
      ENDDO
    ENDIF
    DO JF=1,2*KF_UV
      DO JM=0,IMEN
        IR = ISTAGTF+2*JM+1
        II = IR+1
        ! calculate amplitude and phase shift and reconstruct A,B
        ZAMP = SQRT(YLSPVUV(JF)%P(IR)**2 + YLSPVUV(JF)%P(II)**2)
        ZPHASE = ATAN2(YLSPVUV(JF)%P(II),YLSPVUV(JF)%P(IR)) + REAL(JM,JPRB)*ZSHIFT
        YLSPVUV(JF)%P(IR) =  ZAMP*COS(ZPHASE)
        YLSPVUV(JF)%P(II) =  ZAMP*SIN(ZPHASE)
      ENDDO
    ENDDO
  ENDIF
ENDIF
  
  !     ------------------------------------------------------------------
  
!*       1.    DIVIDE U V AND N-S DERIVATIVES BY A*COS(THETA)
!              ----------------------------------------------

  
!*       1.1      U AND V.

IF(KF_UV > 0) THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,2*KF_UV
      YLSPVUV(JF)%P(JLON) = YLSPVUV(JF)%P(JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!*      1.2      N-S DERIVATIVES

IF(KF_SCDERS > 0)THEN
  DO JLON=ISTAGTF+1,ISTAGTF+2*(IMEN+1)
    DO JF=1,KF_SCALARS
      YLSPVSCALAR_NS(JF)%P(JLON) = YLSPVSCALAR_NS(JF)%P(JLON)*ZACHTE2
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       2.    EAST-WEST DERIVATIVES
!              ---------------------

!*       2.1      U AND V.

IF(LUVDER)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,2*KF_UV
      YLSPVUV_EW(JF)%P(IR) = -YLSPVUV(JF)%P(II)*ZMUL
      YLSPVUV_EW(JF)%P(II) =  YLSPVUV(JF)%P(IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!*       2.2     SCALAR VARIABLES

IF(KF_SCDERS > 0)THEN
  DO JM=0,IMEN
    IR = ISTAGTF+2*JM+1
    II = IR+1
    ZMUL = ZACHTE*JM
    DO JF=1,KF_SCALARS
      YLSPVSCALAR_EW(JF)%P(IR) = -YLSPVSCALAR(JF)%P(II)*ZMUL
      YLSPVSCALAR_EW(JF)%P(II) =  YLSPVSCALAR(JF)%P(IR)*ZMUL
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE FSC_VIEW
END MODULE FSC_VIEW_MOD
