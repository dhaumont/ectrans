! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
SUBROUTINE INV_TRANS_FIELD_API(YD_PSPVOR,   YD_PSPDIV,   YD_PSPSCALAR,   &
                              &YD_PSPSC_3D, YD_PSPSC_2D, FSPGL_PROC,     &
                              &LDSCDERS,    LDVORGP,     LDDIVGP,        &
                              &LDUVDER,     LDLATLON,    KPROMA,         &
                              &KVSETUV,     KVSETSC,     KRESOL,         &
                              &KVSETSC_3D,  KVSETSC_2D,  YD_PGP,         &
                              &YD_PGPUV,    YD_PGP_3D,   YD_PGP_2D,      &
                              &LACC)                                      

                                      

    

!**** *INV_TRANS* - Inverse spectral transform.

!     Purpose.
!     --------
!        Interface routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_FIELD_API(...)

!     Explicit arguments : All arguments are optional.
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (input)
!     PSPDIV(:,:) - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     PSPSC3A(:,:,:) - alternative to use of PSPSCALAR, see PGP3A below (input)

!     PSPSC2(:,:)  - alternative to use of PSPSCALAR, see PGP2 below (input)
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     LDSCDERS    - indicating if derivatives of scalar variables are req.
!     LDVORGP     - indicating if grid-point vorticity is req.
!     LDDIVGP     - indicating if grid-point divergence is req.
!     LDUVDER     - indicating if E-W derivatives of u and v are req.
!     LDLATLON   - indicating if regular lat-lon output requested
!     KPROMA      - required blocking factor for gridpoint output
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     KVSETSC3A(:) - as KVESETSC for PSPSC3A (distribution on first dimension)

!     KVSETSC2(:) - as KVESETSC for PSPSC2 (distribution on first dimension)
!     KRESOL   - resolution tag  which is required ,default is the
!                first defined resulution (input)
!     PGP(:,:,:) - gridpoint fields (output)
!                  PGP need to  dimensioned (NPROMA,IF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, IF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       vorticity     : IF_UV_G fields (if psvor present and LDVORGP)
!       divergence    : IF_UV_G fields (if psvor present and LDDIVGP)
!       u             : IF_UV_G fields (if psvor present)
!       v             : IF_UV_G fields (if psvor present)
!       scalar fields : IF_SCALARS_G fields (if pspscalar present)
!       N-S derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!       E-W derivative of u : IF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of v : IF_UV_G fields (if psvor present and and LDUVDER)
!       E-W derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!
!       Here IF_UV_G is the GLOBAL number of u/v fields as given by the length
!       of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!       IF_SCALARS_G is the GLOBAL number of scalar fields as giben by the
!       length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!       'b-set' split

!     As an alternative to using PGP you can also use a combination of the
!     following arrays. The reason for introducing these alternative ways
!     of calling INV_TRANS is to avoid uneccessary copies where your data
!     structures don't fit in to the 'PSPVOR,PSPDIV, PSPSCALAR, PGP' layout.
!     The use of any of these precludes the use of PGP and vice versa.

!     PGPUV(:,:,:,:) - the 'u-v' related grid-point variables in the order
!                      described for PGP. The second dimension of PGPUV should
!                      be the same as the "global" first dimension of
!                      PSPVOR,PSPDIV (in the IFS this is the number of levels)
!                      PGPUV need to be dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (u,v,vor,div ...)
!     PGP3A(:,:,:,:) - grid-point array directly connected with PSPSC3A
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3A if no derivatives, 3 times that with der.)
!     PGP2(:,:,:)    - grid-point array directly connected with PSPSC2
!                      dimensioned(NPROMA,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC2 if no derivatives, 3 times that with der.)
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  LTINV_CTL   - control of Legendre transform
!                 FTINV_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        26-02-03 Mats Hamrud & Gabor Radnoti : modified condition for scalar fields
!                                               and derivatives (IF_SCALARS_G)

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE FIELD_ACCESS_MODULE
USE FIELD_FACTORY_MODULE
USE FIELD_MODULE

USE YOMHOOK           ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments

CLASS (FIELD_2RB), POINTER, OPTIONAL, INTENT(IN) :: YD_PSPVOR
CLASS (FIELD_2RB), POINTER, OPTIONAL, INTENT(IN) :: YD_PSPDIV
CLASS (FIELD_2RB), POINTER, OPTIONAL, INTENT(IN) :: YD_PSPSCALAR
CLASS (FIELD_3RB), POINTER, OPTIONAL, INTENT(IN) :: YD_PSPSC_3D(:)
CLASS (FIELD_2RB), POINTER, OPTIONAL, INTENT(IN) :: YD_PSPSC_2D(:)

LOGICAL   ,OPTIONAL, INTENT(IN) :: LDSCDERS
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDVORGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDDIVGP
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDUVDER
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDLATLON
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC_3D(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC_2D(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC
CLASS (FIELD_3RB), POINTER, OPTIONAL, INTENT(OUT) :: YD_PGP
CLASS (FIELD_4RB), POINTER, OPTIONAL, INTENT(OUT) :: YD_PGPUV
CLASS (FIELD_3RB), POINTER, OPTIONAL, INTENT(OUT) :: YD_PGP_3D(:)
CLASS (FIELD_2RB), POINTER, OPTIONAL, INTENT(OUT) :: YD_PGP_2D(:)
LOGICAL, INTENT(IN), OPTIONAL :: LACC

!IN
REAL(KIND=JPRB), POINTER  :: PSPVOR(:,:)
REAL(KIND=JPRB), POINTER  :: PSPDIV(:,:)
REAL(KIND=JPRB), POINTER  :: PSPSCALAR(:,:)
REAL(KIND=JPRB), POINTER  :: PSPSC_3D(:,:,:)
REAL(KIND=JPRB), POINTER  :: PSPSC_2D(:,:)

!OUT
REAL(KIND=JPRB),POINTER :: PGP(:,:,:)
REAL(KIND=JPRB),POINTER :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),POINTER :: PGP_3D(:,:,:,:)
REAL(KIND=JPRB),POINTER :: PGP_2D(:,:,:)

REAL(KIND=JPHOOK), POINTER, CONTIGUOUS :: ZHOOK_HANDLE
#include "inv_trans.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INV_TRANS_FIELD_API',0,ZHOOK_HANDLE)

PSPVOR => NULL ()
PSPDIV => NULL ()
PSPSCALAR => NULL ()
PSPSC_3D => NULL ()
PSPSC_2D => NULL ()

if (LACC) THEN
  IF(PRESENT(YD_PSPVOR)) PSPVOR => GET_HOST_DATA_RDONLY (YD_PSPVOR)
  IF(PRESENT(YD_PSPDIV)) PSPDIV => GET_HOST_DATA_RDONLY (YD_PSPDIV)
  IF(PRESENT(YD_PSPSCALAR)) PSPSCALAR => GET_HOST_DATA_RDONLY (YD_PSPSCALAR)
  !IF(PRESENT(YD_PSPSC3A)) PSPSC3A => GET_HOST_DATA_RDONLY (YD_PSPSC3A)
  !IF(PRESENT(YD_PSPSC2)) PSPSC2 => GET_HOST_DATA_RDONLY (YD_PSPSC2)

  IF(PRESENT(YD_PGP)) PGP => GET_HOST_DATA_RDWR (YD_PGP)
  IF(PRESENT(YD_PGPUV)) PGPUV => GET_HOST_DATA_RDWR (YD_PGPUV)
  !IF(PRESENT(YD_PGP3A)) PGP3A => GET_HOST_DATA_RDWR (YD_PGP3A)
  !IF(PRESENT(YD_PGP2)) PGP2 => GET_HOST_DATA_RDWR (YD_PGP2)

ELSE
  IF(PRESENT(YD_PSPVOR)) PSPVOR => GET_DEVICE_DATA_RDONLY (YD_PSPVOR)
  IF(PRESENT(YD_PSPDIV)) PSPDIV => GET_DEVICE_DATA_RDONLY (YD_PSPDIV)
  IF(PRESENT(YD_PSPSCALAR)) PSPSCALAR => GET_DEVICE_DATA_RDONLY (YD_PSPSCALAR)
  !IF(PRESENT(YD_PSPSC3A)) PSPSC3A => GET_DEVICE_DATA_RDONLY (YD_PSPSC3A)
  !IF(PRESENT(YD_PSPSC2)) PSPSC2 => GET_DEVICE_DATA_RDONLY (YD_PSPSC2)

  IF(PRESENT(YD_PGP)) PGP => GET_DEVICE_DATA_RDWR (YD_PGP)
  IF(PRESENT(YD_PGPUV)) PGPUV => GET_DEVICE_DATA_RDWR (YD_PGPUV)
  !IF(PRESENT(YD_PGP3A)) PGP3A => GET_DEVICE_DATA_RDWR (YD_PGP3A)
  !IF(PRESENT(YD_PGP2)) PGP2 => GET_DEVICE_DATA_RDWR (YD_PGP2)
ENDIF

CALL INV_TRANS(PSPVOR=PSPVOR,       PSPDIV=PSPDIV,                       &
              &PSPSCALAR=PSPSCALAR, PSPSC3A=PSPSC_3D,                    &
              &PSPSC2=PSPSC_2D,     FSPGL_PROC=FSPGL_PROC,               &
              &LDSCDERS=LDSCDERS,   LDVORGP=LDVORGP,                     &
              &LDDIVGP=LDDIVGP,     LDUVDER=LDUVDER,                     &
              &LDLATLON=LDLATLON,   KPROMA=KPROMA,                       &
              &KVSETUV=KVSETUV,     KVSETSC=KVSETSC,                     &
              &KRESOL=KRESOL,       KVSETSC3A=KVSETSC_3D,                &
              &KVSETSC2=KVSETSC_2D, PGP=PGP,                             &
              &PGPUV=PGPUV,         PGP3A=PGP_3D,                        &
              &PGP2=PGP_2D)                                               

                                                 

                      

IF (LHOOK) CALL DR_HOOK('INV_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_FIELD_API