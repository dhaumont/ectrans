MODULE EFTDIR_MOD
CONTAINS
SUBROUTINE EFTDIR(PREEL, KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRB
USE PARKIND_ECTRANS, ONLY : JPRBT

USE TPM_DISTR       ,ONLY : D
USE TPM_DIM         ,ONLY : R
#ifdef HAVE_CUFFT
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD
use cudafor
#endif

USE TPM_FFTH        ,ONLY : CREATE_PLAN_FFT, EXECUTE_PLAN_FFT
USE ISO_C_BINDING

!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IRLEN,ICLEN
!INTEGER(KIND=JPIM) :: IPLAN_R2C
TYPE(C_PTR) :: IPLAN_R2C
REAL(KIND=JPRBT)   :: ZSCAL
!REAL(KIND=JPRBT), ALLOCATABLE :: ZGTF2(:,:), ZGTF3(:,:)
!INTEGER(KIND=JPIM) :: IDIM1, IDIM2, LOT,NX

integer :: istat

!     ------------------------------------------------------------------

IRLEN=R%NDLON+R%NNOEXTZG
ICLEN=D%NLENGTF/D%NDGL_FS

!write (*,*) __FILE__, __LINE__
!write (*,*) 'IRLEN = ',IRLEN,'; ICLEN = ',ICLEN
!write (*,*) 'SHAPE(PREEL) = ',SHAPE(PREEL)
!write (*,*) 'D%NLENGTF = ',D%NLENGTF,'; D%NDGL_FS = ',D%NDGL_FS
!write (*,*) 'KFIELDS = ',KFIELDS
!call flush(6)


!NX=128
!ICLEN=NX+2
!LOT=4
!ALLOCATE(ZGTF2(LOT,ICLEN))
!ALLOCATE(ZGTF3(LOT,ICLEN))
!write (*,*) __FILE__, __LINE__

!!$ACC ENTER DATA CREATE(ZGTF2,ZGTF3)
!!$ACC KERNELS
!ZGTF2(:,:) = 0._JPRBT
!ZGTF3(:,:) = 0._JPRBT
!!$ACC END KERNELS

write (*,*) __FILE__, __LINE__; call flush(6)
CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,IRLEN,KFIELDS*D%NDGL_FS)
write (*,*) __FILE__, __LINE__; call flush(6)
CALL EXECUTE_PLAN_FFT(-1,IRLEN,PREEL(1,1),PREEL(1,1),IPLAN_R2C)
write (*,*) __FILE__, __LINE__; call flush(6)

!!$ACC EXIT DATA DELETE(ZGTF2,ZGTF3)

#ifdef HAVE_CUFFT
CALL CREATE_PLAN_FFT (IPLAN_R2C, -1, KN=IRLEN, KLOT=KFIELDS*D%NDGL_FS, &
                    & KISTRIDE=1, KIDIST=ICLEN, KOSTRIDE=1, KODIST=ICLEN/2)
!$acc host_data use_device(PREEL)
CALL EXECUTE_PLAN_FFTC_INPLACE (IPLAN_R2C, -1, PREEL (1, 1))
!$acc end host_data
#endif

ZSCAL = 1._JPRB / REAL (R%NDLON, JPRB)

!$acc kernels present (PREEL) copyin (ZSCAL)
PREEL = ZSCAL * PREEL
!$acc end kernels



!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD
