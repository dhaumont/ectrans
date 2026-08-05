! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_VIEW_CTL_MOD
CONTAINS
SUBROUTINE INV_TRANS_VIEW_CTL(&
                            & YLSPVVOR, YLSPVDIV, YLSPVSCALAR, &  
                            & YLGVU,YLGVV,&
                            & YLGVVOR,YLGVDIV,&
                            & YLGVSCALAR,&
                            & YLGVU_EW,YLGVV_EW,&
                            & YLGVSCALAR_NS, YLGVSCALAR_EW,&
                            & FSPGL_PROC)

!**** *INV_TRANS_VIEW_CTL* - Control routine for inverse spectral transform.

!     Purpose.
!     --------
!        Control routine for the inverse spectral transform

!**   Interface.
!     ----------
!     CALL INV_TRANS_VIEW_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_OUT_LT    - total number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
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
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition
!     PGP(:,:,:)  - gridpoint fields (output)

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):

!       vorticity     : KF_UV_G fields
!       divergence    : KF_UV_G fields
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields
!       N-S derivative of scalar fields : KF_SCALARS_G fields
!       E-W derivative of u : KF_UV_G fields
!       E-W derivative of v : KF_UV_G fields
!       E-W derivative of scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTINV_CTL   - control of Legendre transform
!                 FTINV_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE TPM_GEN         ,ONLY : NPROMATR
USE TPM_TRANS       ,ONLY : LDIVGP, LSCDERS, LUVDER, LVORGP

USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT
USE LTINV_VIEW_CTL_MOD   ,ONLY : LTINV_VIEW_CTL
USE FTINV_VIEW_CTL_MOD   ,ONLY : FTINV_VIEW_CTL
USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW

IMPLICIT NONE


! Declaration of arguments
TYPE(SPEC_VIEW) :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW) :: YLSPVSCALAR(:)

TYPE(GRID_VIEW) :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW) :: YLGVVOR(:),YLGVDIV(:)
TYPE(GRID_VIEW) :: YLGVSCALAR(:)

TYPE(GRID_VIEW) :: YLGVU_EW(:),YLGVV_EW(:)
TYPE(GRID_VIEW) :: YLGVSCALAR_NS(:), YLGVSCALAR_EW(:)


EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC


INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP,IF_OUT_LT
INTEGER(KIND=JPIM) :: IF_SCDERS,IF_UV_PAR
INTEGER(KIND=JPIM) :: IF_SC2_G,IF_SC3A_G2,IF_SC3A_G3,IF_SC3B_G2,IF_SC3B_G3
! :TODO:
! total number of fields in grid point space

IF_UV = SIZE(YLGVU) 
IF_SCALARS = SIZE(YLGVSCALAR)
IF_SCDERS = SIZE(YLGVSCALAR_NS) 
LVORGP = SIZE(YLGVVOR) > 0
LDIVGP = SIZE(YLGVDIV) > 0
LUVDER = SIZE(YLGVU_EW) > 0 .OR. SIZE(YLGVV_EW) > 0
LSCDERS = SIZE(YLGVSCALAR_NS) > 0 .OR. SIZE(YLGVSCALAR_EW) > 0

WRITE(6,*) "INV_TRANS_view_ctl: IF_UV = ",IF_UV," IF_SCALARS = ",IF_SCALARS," IF_SCDERS = ",IF_SCDERS
WRITE(6,*) "LVORGP = ",LVORGP," LDIVGP = ",LDIVGP," LUVDER = ",LUVDER," LSCDERS = ",LSCDERS
IF_GP = SIZE(YLGVU) + SIZE(YLGVV) + SIZE(YLGVSCALAR) + SIZE(YLGVVOR) + SIZE(YLGVDIV) +&
        &SIZE(YLGVU_EW) + SIZE(YLGVV_EW) + SIZE(YLGVSCALAR_NS) + SIZE(YLGVSCALAR_EW)


!total number of fields coming out from inverse LT
IF_OUT_LT = 2*IF_UV + IF_SCALARS+IF_SCDERS

IF(IF_UV > 0 .AND. LVORGP) THEN
  IF_OUT_LT = IF_OUT_LT+IF_UV
ENDIF
IF(IF_UV > 0 .AND. LDIVGP) THEN
  IF_OUT_LT = IF_OUT_LT+IF_UV
ENDIF
IF_FS = IF_OUT_LT+IF_SCDERS
IF(IF_UV > 0 .AND. LUVDER) THEN
  IF_FS = IF_FS+2*IF_UV
ENDIF


!     ------------------------------------------------------------------

! Perform transform

IF(NPROMATR > 0) THEN
  CALL ABOR1('INV_TRANS_VIEW NPROMATR not implemented')
ELSE

  ! No splitting of fields, transform done in one go
  CALL LTINV_VIEW_CTL(IF_OUT_LT,&
   &YLSPVVOR=YLSPVVOR,YLSPVDIV=YLSPVDIV,YLSPVSCALAR=YLSPVSCALAR,&
   &FSPGL_PROC=FSPGL_PROC)

  CALL FTINV_VIEW_CTL(IF_GP, IF_FS, IF_OUT_LT, &
   & YLGVU=YLGVU,YLGVV=YLGVV,YLGVSCALAR=YLGVSCALAR,&
   & YLGVVOR=YLGVVOR,YLGVDIV=YLGVDIV,&
   & YLGVU_EW=YLGVU_EW,YLGVV_EW=YLGVV_EW, &
   & YLGVSCALAR_NS=YLGVSCALAR_NS,YLGVSCALAR_EW=YLGVSCALAR_EW)
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_VIEW_CTL
END MODULE INV_TRANS_VIEW_CTL_MOD
