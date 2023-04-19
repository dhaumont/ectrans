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
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPM_DIM         ,ONLY : R
USE CUDA_DEVICE_MOD
use cudafor

!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
REAL(KIND=JPRB), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IRLEN,ICLEN
INTEGER(KIND=JPIM) :: IPLAN_R2C
REAL(KIND=JPRBT)   :: ZSCAL

integer :: istat

!     ------------------------------------------------------------------

IRLEN=R%NDLON+R%NNOEXTZG
ICLEN=D%NLENGTF/D%NDGL_FS

!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)

CALL CREATE_PLAN_FFT (IPLAN_R2C, -1, KN=IRLEN, KLOT=KFIELDS*D%NDGL_FS, &
                    & KISTRIDE=1, KIDIST=ICLEN, KOSTRIDE=1, KODIST=ICLEN/2)

!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)

!$acc host_data use_device(PREEL)
!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)
CALL EXECUTE_PLAN_FFTC_INPLACE (IPLAN_R2C, -1, PREEL (1, 1))
!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)
!$acc end host_data

!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)

ZSCAL = 1._JPRB / REAL (R%NDLON, JPRB)

!$acc kernels present (PREEL) copyin (ZSCAL)
PREEL = ZSCAL * PREEL
!$acc end kernels

!write (0,*) __FILE__, __LINE__,'; cudaDeviceSynchronize returns ',cudaDeviceSynchronize(); call flush(0)

!     ------------------------------------------------------------------

END SUBROUTINE EFTDIR
END MODULE EFTDIR_MOD
