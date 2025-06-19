module dir_trans_field_api_mod
use field_api_ectrans_mod
use parkind1  ,only : jpim     ,jprb

implicit none


contains

subroutine dir_trans_field_api(ydfspvor,ydfspdiv,ydfspscalar, &
                             & ydfu, ydfv, ydfvor,ydfdiv,ydfscalar, &
                             & kspec, kproma, kgpblks, kgptot, kflevg, kflevl,&
                             & ldacc, ldverbose)


type(field_basic_ptr),intent(in), optional  :: ydfspvor(:), ydfspdiv(:)        !spectral vector fields : vorticity and divergence fields (in)
type(field_basic_ptr),intent(in), optional  :: ydfspscalar(:)                !spectral scalar fields (in)

type(field_basic_ptr),intent(in), optional  :: ydfu(:),ydfv(:)                 !grid vector fields     (out)
type(field_basic_ptr),intent(in), optional  :: ydfvor(:),ydfdiv(:)             !grid vector fields :vorticity and divergence     (out)
type(field_basic_ptr),intent(in), optional  :: ydfscalar(:)                  !grid scalar fields     (out)

integer(kind=jpim), intent(in) ::kspec
integer(kind=jpim), intent(in) ::kproma
integer(kind=jpim), intent(in) ::kgpblks
integer(kind=jpim), intent(in) ::kgptot
integer(kind=jpim), intent(in) :: kflevg
integer(kind=jpim), intent(in) :: kflevl
logical, intent(in), optional  :: ldacc
logical, intent(in), optional  :: ldverbose

! local variables

! temporary arrays:

type(spec_view), allocatable  :: ylspvvor(:), ylspvdiv(:)
type(spec_view), allocatable  :: ylspvscalar(:)

type(grid_view), allocatable  :: ylgvu(:),ylgvv(:)
type(grid_view), allocatable  :: ylgvvor(:),ylgvdiv(:)
type(grid_view), allocatable  :: ylgvscalar(:)

real(kind=jprb),allocatable :: zpspvor(:,:),zpspdiv(:,:)          ! spectral vector fields (in)
real(kind=jprb),allocatable :: zpspsc2(:,:)                       ! spectral surface scalar fields(in)
real(kind=jprb),allocatable :: zpgpuv(:,:,:,:)                    ! grid vector fields (out)
real(kind=jprb),allocatable :: zpgp2(:,:,:)                       ! grid surface scalar fields (out)

integer(kind=jpim)          :: ispuv                              ! number of input spectral vector fields
integer(kind=jpim)          :: ifldxg
integer(kind=jpim)          :: ifldxl
integer(kind=jpim)          :: ifldxguv
integer(kind=jpim)          :: ifldxluv
integer(kind=jpim)          :: iuvg                               ! number of output vector fields
integer(kind=jpim)          :: iscdim                             ! size of output scalar fields array
integer(kind=jpim)          :: iuvdim                             ! size of output vector fields array
integer(kind=jpim)          :: id,ioffset,jlev
integer(kind=jpim)          :: iend


integer(kind=jpim),allocatable :: ivsetuv(:)
integer(kind=jpim),allocatable :: ivsetsc2(:)

integer(kind=jpim)          :: jfld                                 ! field counter
logical          :: llscders                                        ! indicating if derivatives of scalar variables are req.
logical          :: llvorgp                                         ! indicating if grid-point vorticity is req.
logical          :: lldivgp                                         ! indicating if grid-point divergence is req.
logical          :: lluvder                                         ! indicating if e-w derivatives of u and v are req.
logical          :: llverbose                                       ! indicating if verbose output is req.
#include "dir_trans.h"
#include "abor1.intfb.h"

llverbose = .false.
if (present(ldverbose))  llverbose = ldverbose

! 1. vector fields transformation

! check if all provided vector field information is consistent
if (present(ydfu) .neqv. present(ydfv)) call abor1("[ectrans_api] iydu/ydfv")
if (present(ydfu) .neqv. present(ydfspvor)) call abor1("[ectrans_api] iydu/ydfvor")
if (present(ydfu) .neqv. present(ydfspdiv)) call abor1("[ectrans_api] iydu/ydfdiv")

iuvdim = 0
! do we have vector fields?

if (present(ydfu)) then

  if ((size(ydfu)/= size(ydfv)).or.(size(ydfu)/= size(ydfspdiv)).or.(size(ydfu)/= size(ydfspvor))) then
     call abor1("[ectrans_api] invalid list sizes:")
  endif
  ylspvvor = ls (ydfspvor, ldacc)
  ylspvdiv = ls (ydfspdiv, ldacc)

  ylgvu = lg (ydfu, ldacc)
  ylgvv = lg (ydfv, ldacc)
  if ((size (ylgvu) /= size (ylgvv)) .or. (size (ylspvvor) /= size (ylspvdiv))) then
     call abor1("[ectrans_api] inconsistent number of field_view for vectors:")
  endif
   if (((size (ylgvu) / size (ydfu)) /= kflevg) .or. ((size (ylspvvor) / size (ydfspvor)) /= kflevl)) then
     call abor1("[ectrans_api] inconsistent kflevg or kflevl")
  endif

  iuvg = size(ydfu)           ! number of output  vector fields
  ispuv = size(ydfspvor)

  iuvdim = 2

  if (present(ydfdiv)) then
     lldivgp = .true.
     iuvdim = 5
     ylgvdiv = lg (ydfdiv, ldacc)
  endif

  if (present(ydfvor)) then
     llvorgp = .true.
     iuvdim = 6
     ylgvvor = lg (ydfvor, ldacc)
  endif
   
  ! allocate vector field input in spectral space
  allocate(zpspvor(size(ylspvvor),kspec))
  allocate(zpspdiv(size(ylspvdiv),kspec))

  ! allocate vector field output in grid space
  allocate(zpgpuv(kproma,kflevg, iuvg * iuvdim,kgpblks))
  allocate(ivsetuv(kflevg))

  ! copy from fields to temporary arrays (1d copy thanks to field view)

  do jfld=1,size(ylspvvor)
     zpspvor(jfld,:) = ylspvvor(jfld)%view%p(:)
     zpspdiv(jfld,:) = ylspvdiv(jfld)%view%p(:)
  enddo


else
  ! ydfu is not provided, we do not have to compute the corresponding vector output
  ispuv = 0
  allocate(zpgpuv(0,ispuv,2,0),zpspvor(0,0),zpspdiv(0,0))
endif

! 2. scalar fields transformation

! check if all provided scalar field information is consistent
if (present(ydfspscalar) .neqv. present(ydfscalar))  call abor1("[ectrans_api] grid/spec")

if (present(ydfspscalar)) then

  if ((size(ydfspscalar)/= size(ydfscalar)))  call abor1("[ectrans_api] inconsistent number of field_view for ydfscalar")

  ylgvscalar = lg (ydfscalar, ldacc)
  ylspvscalar = ls (ydfspscalar, ldacc)

  ifldxg = size(ylgvscalar) ! number of output scalar fields in grid space

  ifldxl = 0  ! number of input scalar fields in spectral space
  do jfld = 1, size(ylgvscalar) ! number of output scalar fields in grid space
  if (associated(ylspvscalar(jfld)%view%p)) ifldxl = ifldxl + 1
  end do 

  iscdim = 1

   ! allocate scalar field input in spectral space
  allocate(zpspsc2(ifldxl,kspec))

  ! allocate scalar field output in grid space
  allocate(zpgp2(kproma,ifldxg * iscdim,kgpblks))
  allocate(ivsetsc2(ifldxg))

 ! copy scalar spectral fields to temporary arrays (1d copy thanks to field view)
  id = 1
  do jfld = 1, size(ylspvscalar) ! number of output scalar fields in grid space
     if (associated(ylspvscalar(jfld)%view%p)) then
         zpspsc2(id,:) = ylspvscalar(jfld)%view%p(:)
         id = id + 1 
     endif
  enddo

 do jfld=1, ifldxg
    ivsetsc2(jfld) = ylgvscalar(jfld)%ivset
 enddo


do jfld=1,iuvg
  do jlev=1,kflevg
     id = jlev + (jfld -1) * kflevg
     if (jfld .eq. 1) ivsetuv(jlev) = ylgvu(id)%ivset
     if (ivsetuv(jlev) .ne. ylgvu(id)%ivset)  call abor1("[ectrans_api] ivsetuv inconsistent with ylgvu%ivset")
     if (ivsetuv(jlev) .ne. ylgvv(id)%ivset)  call abor1("[ectrans_api] ivsetuv inconsistent with ylgvv%ivset")
  enddo
enddo
!
else
  ifldxg = 0
  allocate(zpgp2(0,ifldxg,0),zpspsc2(0,1))
endif


! 3. call dir_trans

	call dir_trans(pspvor = zpspvor,pspdiv = zpspdiv,pgpuv = zpgpuv,kvsetuv = ivsetuv, &
	             & pspsc2 = zpspsc2,pgp2 = zpgp2, kvsetsc2 = ivsetsc2, &
	             & kproma = kproma)
! 4. copy back data to fields

! remove garbage at the end of arrays
  iend = kgptot - kproma * (kgpblks - 1)
  zpgpuv (iend+1:, :, :, kgpblks) = 0
  zpgp2 (iend+1:, :, kgpblks) = 0

  ! copy vector fields back from temporary vector arrays
ioffset = 0

if (llvorgp) then
  do jfld=1,iuvg
    do jlev=1,kflevg
     id = jlev + (jfld -1) * kflevg
     ylgvvor(id)%view%p(:,:) = zpgpuv(:, jlev,jfld+ioffset*iuvg,:)
    enddo
  enddo
  ioffset = ioffset + 1
endif

if (lldivgp) then
  do jfld=1,iuvg
    do jlev=1,kflevg
     id = jlev + (jfld -1) * kflevg
     ylgvdiv(id)%view%p(:,:) = zpgpuv(:, jlev,jfld+ioffset*iuvg,:)
    enddo
  enddo
  ioffset = ioffset + 1
endif

do jfld=1,iuvg
  do jlev=1,kflevg
     id = jlev + (jfld -1) * kflevg
     ylgvu(id)%view%p(:,:) =  zpgpuv(:,jlev,jfld+ioffset*iuvg,:)
     ylgvv(id)%view%p(:,:) =  zpgpuv(:,jlev,jfld+(ioffset+1)*iuvg,:)
  enddo
enddo

! copy scalar fields back from temporary scalar arrays
do jfld=1, ifldxg
    ylgvscalar(jfld)%view%p(:,:) = zpgp2(:,jfld,:)
enddo


deallocate(zpspvor)
deallocate(zpspdiv)
deallocate(zpspsc2)
deallocate(zpgpuv)
deallocate(zpgp2)
deallocate(ivsetuv)
deallocate(ivsetsc2)

end subroutine dir_trans_field_api

end module dir_trans_field_api_mod

