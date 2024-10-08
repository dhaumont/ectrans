! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE DIR_TRANS_FIELD_API(YD_FIELDS_PGP,YD_FIELDS_PSP, KVSETSC, &
                               &LDVORGP,  LDDIVGP, &
                               &LDUVDER,   LDLATLON,  KPROMA,KRESOL,      &
                               &LACC)


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE FIELD_ACCESS_MODULE
USE FIELD_FACTORY_MODULE
USE FIELD_MODULE
USE FIELD_API_INTERFACE

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR, NOUT
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP, LATLON, &
     &                      NF_SC2, NF_SC3A, NF_SC3B,        &
     &                      NGPBLKS, NPROMA
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DIR_TRANS_CTL_MOD ,ONLY : DIR_TRANS_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

TYPE (FIELDS_PGP), INTENT(IN) :: YD_FIELDS_PGP
TYPE (FIELDS_PSP), INTENT(OUT) :: YD_FIELDS_PSP
TYPE (FIELDS_KVSETC), INTENT(IN) :: KVSETSC

LOGICAL   ,OPTIONAL, INTENT(IN) :: LDVORGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDDIVGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDUVDER
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDLATLON

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPROMA

LOGICAL, INTENT(IN), OPTIONAL :: LACC
! LOCAL VARIABLES

!TEMPORARY INPUT
REAL(KIND=JPRB),POINTER :: PGP(:,:,:)
REAL(KIND=JPRB),POINTER :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),POINTER :: PGP_3D(:,:,:,:)
REAL(KIND=JPRB),POINTER :: PGP_2D(:,:,:)

!TEMPORARY OUTPUT
REAL(KIND=JPRB)    ,POINTER :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,POINTER :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,POINTER :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,POINTER :: PSPSC_3D(:,:,:)
REAL(KIND=JPRB)    ,POINTER :: PSPSC_2D(:,:,:)


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#include "dir_trans.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!CALL DIR_TRANS(PSPVOR, PSPDIV,    PSPSCALAR, PSPSC3A,  PSPSC3B,          &
 !             &PSPSC2, LDLATLON,  KPROMA,    KVSETUV,  KVSETSC,          &
  !            &KRESOL, KVSETSC3A, KVSETSC3B, KVSETSC2, PGP,              &
   !           &PGPUV,  PGP3A,     PGP3B,     PGP2 )                        

IF (LHOOK) CALL DR_HOOK('DIR_TRANS_FIELD_API',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------
!endif INTERFACE

END SUBROUTINE DIR_TRANS_FIELD_API
