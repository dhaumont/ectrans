! (C) Copyright 2001- ECMWF.
! (C) Copyright 2001- Meteo-France.
! (C) Copyright 2022- NVIDIA.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE INV_TRANS_VIEW_CTL_MOD
CONTAINS
  SUBROUTINE INV_TRANS_VIEW_CTL(KPROMA,KGPBLKS, &
                            & YDSPVVOR, YDSPVDIV, YDSPVSCALAR, &
                            & YDGVU,YDGVV,&
                            & YDGVVOR,YDGVDIV,&
                            & YDGVSCALAR,&
                            & YDGVU_EW,YDGVV_EW,&
                            & YDGVSCALAR_NS, YDGVSCALAR_EW,&
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
    
    !     KF_OUT_LT    - total number of fields coming out from inverse LT
    !     IF_UV        - local number of spectral u-v fields
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

    !       vorticity     : IF_UV_G fields
    !       divergence    : IF_UV_G fields
    !       u             : IF_UV_G fields
    !       v             : IF_UV_G fields
    !       scalar fields : IF_SCALARS_G fields
    !       N-S derivative of scalar fields : IF_SCALARS_G fields
    !       E-W derivative of u : IF_UV_G fields
    !       E-W derivative of v : IF_UV_G fields
    !       E-W derivative of scalar fields : IF_SCALARS_G fields

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


    USE PARKIND_ECTRANS,        ONLY: JPIM, JPRB, JPRBT, JPRD
    USE TPM_GEN,                ONLY: NPROMATR
    USE TPM_DISTR,              ONLY : MYSETV
    USE TPM_TRANS,              ONLY: LDIVGP, LSCDERS, LUVDER, LVORGP, GROWING_ALLOCATION, NPROMA, NGPBLKS
    USE ABORT_TRANS_MOD,        ONLY: ABORT_TRANS
    USE BUFFERED_ALLOCATOR_MOD, ONLY: BUFFERED_ALLOCATOR, MAKE_BUFFERED_ALLOCATOR, &
      &                               INSTANTIATE_ALLOCATOR
    USE TRMTOL_MOD,             ONLY: PREPARE_TRMTOL, TRMTOL_HANDLE, TRMTOL
    USE LTINV_VIEW_MOD,         ONLY: PREPARE_LTINV, LTINV_HANDLE, LTINV_VIEW
    USE TRMTOL_PACK_UNPACK,     ONLY: TRMTOL_PACK_HANDLE, TRMTOL_UNPACK_HANDLE, &
      &                               PREPARE_TRMTOL_PACK, PREPARE_TRMTOL_UNPACK, TRMTOL_PACK, &
      &                               TRMTOL_UNPACK
    USE FSC_MOD,                ONLY: FSC
    USE FTINV_MOD,              ONLY: FTINV_HANDLE, PREPARE_FTINV, FTINV
    USE TRLTOG_VIEW_MOD,        ONLY: TRLTOG_HANDLE, PREPARE_TRLTOG, TRLTOG_VIEW
    USE ECTRANS_FIELD_VIEW_INTERNAL_UTIL_MOD, ONLY : SPEC_VIEW, GRID_VIEW

    IMPLICIT NONE

    #include "fspgl_intf.h"
    ! Declaration of arguments
    INTEGER(KIND=JPIM) :: KPROMA, KGPBLKS
    TYPE(SPEC_VIEW) :: YDSPVVOR(:), YDSPVDIV(:)
    TYPE(SPEC_VIEW) :: YDSPVSCALAR(:)

    TYPE(GRID_VIEW) :: YDGVU(:),YDGVV(:)
    TYPE(GRID_VIEW) :: YDGVVOR(:),YDGVDIV(:)
    TYPE(GRID_VIEW) :: YDGVSCALAR(:)

    TYPE(GRID_VIEW) :: YDGVU_EW(:),YDGVV_EW(:)
    TYPE(GRID_VIEW) :: YDGVSCALAR_NS(:), YDGVSCALAR_EW(:)

    PROCEDURE (FSPGL_INTF), POINTER, OPTIONAL, INTENT(IN)  :: FSPGL_PROC

    ! Local variables

    REAL(KIND=JPRB), POINTER :: FOUBUF(:), FOUBUF_IN(:)
    REAL(KIND=JPRBT), POINTER :: PREEL_REAL(:), PREEL_COMPLEX(:)
    REAL(KIND=JPRBT), POINTER :: ZOUTS(:), ZOUTA(:)
    REAL(KIND=JPRD), POINTER :: ZOUTS0(:), ZOUTA0(:)
    INTEGER(KIND=JPIM) :: KUV_OFFSET, KSCALARS_OFFSET, KSCALARS_NSDER_OFFSET, &
        & KUV_EWDER_OFFSET, KSCALARS_EWDER_OFFSET
    INTEGER(KIND=JPIM) :: IF_LEG, IF_FOURIER
    INTEGER(KIND=JPIM) :: IF_GP
    INTEGER(KIND=JPIM) :: IFIRST
    TYPE(GRID_VIEW),ALLOCATABLE :: YLGP(:)
    INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP,IF_OUT_LT
    INTEGER(KIND=JPIM) :: IF_SCDERS,IF_SCDERS_G,IF_UV_PAR
    
    TYPE(BUFFERED_ALLOCATOR) :: ALLOCATOR
    TYPE(LTINV_HANDLE) :: HLTINV
    TYPE(TRMTOL_PACK_HANDLE) :: HTRMTOL_PACK
    TYPE(TRMTOL_HANDLE) :: HTRMTOL
    TYPE(TRMTOL_UNPACK_HANDLE) :: HTRMTOL_UNPACK
    TYPE(FTINV_HANDLE) :: HFTINV
    TYPE(TRLTOG_HANDLE) :: HTRLTOG
    
    INTEGER(KIND=JPIM) :: I,J
    INTEGER(KIND=JPIM), ALLOCATABLE :: IVSET(:)

    REAL(KIND=JPRB), POINTER :: ZZ2 (:,:)


    INTEGER :: M, N, P
    CHARACTER*1 :: CLENV
    !     ------------------------------------------------------------------

    NPROMA = KPROMA
    NGPBLKS = KGPBLKS

    YLGP = [YDGVVOR,YDGVDIV,YDGVU,YDGVV,YDGVSCALAR,YDGVSCALAR_NS,YDGVU_EW,YDGVV_EW,YDGVSCALAR_EW]

IF_FOURIER = 5043

M = 5043
P = 28480
N = 145238400
ALLOCATE (PREEL_REAL (N))

!$ACC ENTER DATA CREATE (PREEL_REAL)

CLENV = '0'

CALL GETENV ('ECTRANS_USE_DEVICE', CLENV)

IF (CLENV == '1') THEN

  DO J = 1, SIZE (YLGP)
    ZZ2 => YLGP (J)%P
    !$ACC HOST_DATA USE_DEVICE (ZZ2)
    YLGP (J)%P => ZZ2
    !$ACC END HOST_DATA
  ENDDO

ENDIF

!$ACC KERNELS PRESENT (PREEL_REAL)
PREEL_REAL = 0
!$ACC END KERNELS

CALL TRLTOG_VIEW(PREEL_REAL,IF_FOURIER,YLGP)

!$ACC EXIT DATA DELETE (PREEL_REAL)
DEALLOCATE (PREEL_REAL)

  END SUBROUTINE INV_TRANS_VIEW_CTL
END MODULE INV_TRANS_VIEW_CTL_MOD
