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
SUBROUTINE INV_TRANS_VIEW_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_OUT_LT,&
 & KF_UV,KF_SCALARS,KF_SCDERS,&
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
!SE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW, LS_COUNT, LG_COUNT, LS, LG, IVSET_PTR, &
 !                                             & GET_NPROMA, GET_NFLD, GET_NSPEC2, GET_NBLK
!

IMPLICIT NONE
TYPE SPEC_VIEW
  ! Spectral field view
  REAL(KIND=JPRB),POINTER  :: P(:) => NULL()
  CHARACTER(LEN=:),POINTER :: NAME => NULL()
END TYPE

TYPE GRID_VIEW
  ! Grid point field view
  REAL(KIND=JPRB),POINTER   :: P(:,:) => NULL()
  CHARACTER(LEN=:), POINTER :: NAME   => NULL()
  INTEGER(KIND=JPIM)        :: IVSET
END TYPE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS

TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVVOR(:), YLSPVDIV(:)
TYPE(SPEC_VIEW), ALLOCATABLE :: YLSPVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU(:),YLGVV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVVOR(:),YLGVDIV(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR(:)

TYPE(GRID_VIEW), ALLOCATABLE :: YLGVU_EW(:),YLGVV_EW(:)
TYPE(GRID_VIEW), ALLOCATABLE :: YLGVSCALAR_NS(:), YLGVSCALAR_EW(:)

EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC


! Local variables

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_SCDERS,IF_OUT_LT
INTEGER(KIND=JPIM) :: IOFFD,IOFFU,IOFFV,IOFFUVD,IOFFSC,IOFFSCNS,IOFFSCEW,IOFF,IF_GPB

!     ------------------------------------------------------------------

! Perform transform

IF_GPB = 2*KF_UV_G+KF_SCALARS_G
IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN

  CALL ABOR1('INV_TRANS_VIEW NPROMATR not implemented')

ELSE

  ! No splitting of fields, transform done in one go
  CALL LTINV_VIEW_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS, &
   &YLSPVVOR=YLSPVVOR,YLSPVDIV=YLSPVDIV,YLSPVSCALAR=YLSPVSCALAR,&
   &FSPGL_PROC=FSPGL_PROC)

  CALL FTINV_VIEW_CTL(KF_UV_G,KF_SCALARS_G,&
   & KF_UV,KF_SCALARS,KF_SCDERS,KF_GP,KF_FS,KF_OUT_LT,&
   & YLGVU=YLGVU,YLGVV=YLGVV,YLGVSCALAR=YLGVSCALAR,&
   & YLGVU_EW=YLGVU_EW,YLGVV_EW=YLGVV_EW, &
   & YLGVSCALAR_NS=YLGVSCALAR_NS,YLGVSCALAR_EW=YLGVSCALAR_EW
   )
ENDIF

!     ------------------------------------------------------------------

END SUBROUTINE INV_TRANS_VIEW_CTL
END MODULE INV_TRANS_VIEW_CTL_MOD
