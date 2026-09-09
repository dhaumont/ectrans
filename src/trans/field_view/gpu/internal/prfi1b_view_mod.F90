! (C) Copyright 2000- ECMWF.
! (C) Copyright 2000- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI1B_VIEW_MOD
  CONTAINS
  SUBROUTINE PRFI1B_VIEW(PIA,YDSP)
  
  USE PARKIND1,        ONLY: JPIM, JPRB
  USE TPM_DIM,         ONLY: R
  USE TPM_DISTR,       ONLY: D
  USE ABORT_TRANS_MOD, ONLY: ABORT_TRANS
  USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY: SPEC_VIEW

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
  !        *CALL* *PRFI1B_VIEW(...)*
  
  !        Explicit arguments :  KM     - zonal wavenumber
  !        ------------------    PIA    - spectral components for transform
  !                              YDSP    - spectral arrays
    
  
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
  !        Original : 00-02-01 From PRFI1B_VIEW in IFS CY22R1
  
  !     ------------------------------------------------------------------
  
  IMPLICIT NONE
    
  INTEGER(KIND=JPIM) :: KM,KMLOC
  TYPE(SPEC_VIEW), INTENT(IN) :: YDSP(:)
  REAL(KIND=JPRB)   ,INTENT(INOUT)  :: PIA(:,:,:)
      
  !     LOCAL INTEGER SCALARS
  INTEGER(KIND=JPIM) :: INM, IR, JN, JFLD, IASM0, IFIELDS
  INTEGER, PARAMETER :: ISIZE=32
  REAL(KIND=8) :: ZTEMP1(ISIZE+1,ISIZE)
  REAL(KIND=8) :: ZTEMP2(ISIZE+1,ISIZE)

INTEGER :: JNSEC,JFLDSEC,JNSECEND,JFLDSECEND
  !     ------------------------------------------------------------------
  
  !*       1.    EXTRACT FIELDS FROM SPECTRAL ARRAYS.
  !              --------------------------------------------------

  ASSOCIATE(D_NUMP=>D%NUMP, D_MYMS=>D%MYMS, D_NASM0=>D%NASM0, R_NSMAX=>R%NSMAX)
IFIELDS = SIZE(YDSP)
#ifdef ACCGPU
  !$ACC DATA PRESENT(D,D_NUMP,R,R_NSMAX,D_MYMS,D_NASM0,PIA) COPYIN(YDSP) ASYNC(1)
#endif
#ifdef OMPGPU
  !$OMP TARGET DATA MAP(PRESENT,ALLOC:D,D_NUMP,R,R_NSMAX,D_MYMS,D_NASM0,PIA,YDSP)
#endif

  
    !loop over wavenumber


!$ACC PARALLEL PRESENT(YDSP,PIA) DEFAULT(NONE) VECTOR_LENGTH(32) &
!$ACC& FIRSTPRIVATE(IFIELDS) &
!$ACC& PRIVATE(KMLOC,JN,JFLD,KM,IASM0,INM,JNSEC,JFLDSEC,JNSECEND,JFLDSECEND,ZTEMP1,ZTEMP2)
!$ACC LOOP GANG COLLAPSE(3) PRIVATE(ZTEMP1,ZTEMP2)
  DO KMLOC=1,D_NUMP
    DO JN=0,R_NSMAX+3, ISIZE
      DO JFLD=1,IFIELDS, ISIZE
        !$ACC CACHE(ZTEMP1,ZTEMP2)
        JNSECEND   = MIN(R_NSMAX+3, JN+ISIZE-1)
        JFLDSECEND = MIN(IFIELDS,   JFLD+ISIZE-1)
        KM = D_MYMS(KMLOC)
        IASM0 = D_NASM0(KM)

        !$ACC LOOP VECTOR COLLAPSE(2) PRIVATE(JNSEC,JFLDSEC,INM)
          DO JNSEC=JN,JNSECEND
            DO JFLDSEC=JFLD,JFLDSECEND
              IF (JNSEC > 1 .AND. JNSEC <= R_NSMAX+2-KM) THEN
                  INM = IASM0+((R_NSMAX+2-JNSEC)-KM)*2
                  ZTEMP1(JNSEC-JN+1,JFLDSEC-JFLD+1) = YDSP(JFLDSEC)%P(INM)
                  ZTEMP2(JNSEC-JN+1,JFLDSEC-JFLD+1) = YDSP(JFLDSEC)%P(INM+1)
              ENDIF
            ENDDO
          ENDDO

          !$ACC LOOP VECTOR COLLAPSE(2) PRIVATE(JNSEC,JFLDSEC)
          DO JFLDSEC=JFLD,JFLDSECEND
            DO JNSEC=JN,JNSECEND
              IF (JNSEC+1 <= UBOUND(PIA,2)) THEN
                IF (JNSEC <= 1) THEN
                    PIA(2*JFLDSEC-1,JNSEC+1,KMLOC) = 0.0_JPRB
                    PIA(2*JFLDSEC  ,JNSEC+1,KMLOC) = 0.0_JPRB
                ELSEIF (JNSEC <= R_NSMAX+2-KM) THEN
                    PIA(2*JFLDSEC-1,JNSEC+1,KMLOC) = ZTEMP1(JNSEC-JN+1,JFLDSEC-JFLD+1)
                    PIA(2*JFLDSEC  ,JNSEC+1,KMLOC) = ZTEMP2(JNSEC-JN+1,JFLDSEC-JFLD+1)
                ELSEIF (JNSEC <= R_NSMAX+3-KM) THEN
                    PIA(2*JFLDSEC-1,JNSEC+1,KMLOC) = 0.0_JPRB
                    PIA(2*JFLDSEC  ,JNSEC+1,KMLOC) = 0.0_JPRB
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
  ENDDO

 !$ACC END PARALLEL

#ifdef ACCGPU
  !$ACC END DATA
#endif
#ifdef OMPGPU
  !$OMP END TARGET DATA
#endif

  END ASSOCIATE

  !     ------------------------------------------------------------------

  END SUBROUTINE PRFI1B_VIEW
END MODULE PRFI1B_VIEW_MOD
