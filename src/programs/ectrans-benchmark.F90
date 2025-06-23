! (C) Copyright 2014- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

program ectrans_benchmark

#ifdef USE_PINNED
#define PINNED_TAG , pinned
#else
#define PINNED_TAG
#endif

!
! Spectral transform test
!
! This test performs spectral to real and real to spectral transforms repeated in
! timed loop.
!
! 1) One "surface" field is always transformed:
!      zspsc2(1,1:nspec2) <-> zgmvs(1:nproma,1:1,1:ngbplk)
!
! 2) A Multiple "3d" fields are transformed and can be disabled with "--nfld 0"
!
!      zspsc3a(1:nlev,1:nspec2,1:nfld) <-> zgp3a(1:nproma,1:nlev,1:nfld,1:ngpblk)
!
! 3) Optionally a "3d" vorticity/divergence field is transformed to uv (wind) and
!   can be enabled with "--vordiv"
!
!      zspvor(1:nlev,1:nspec2) / zspdiv(1:nlev,1:nspec2) <-> zgpuv(1:nproma,1:nlev,1:2,1:ngpblk)
!
! 4) Optionally scalar derivatives can be computed for the fields described in 1) and 2)
!    This must be enabled with "--scders"
!
! 5) Optionally uv East-West derivate can be computed from vorticity/divergence.
!    This must be enabled with "--vordiv --uvders"
!
!
! Authors : George Mozdzynski
!           Willem Deconinck
!           Ioan Hadade
!           Sam Hatfield
!

use parkind1, only: jpim, jprb, jprd
use oml_mod ,only : oml_max_threads
use mpl_module
use yomgstats, only: jpmaxstat, gstats_lstats => lstats
use yomhook, only : dr_hook_init
USE FIELD_MODULE, ONLY:FIELD_1RB, FIELD_2RB, FIELD_3RB, FIELD_4RB
USE FIELD_FACTORY_MODULE
#include "field_basic_type_ptr.h"
#include "field_api_ectrans.h"
implicit none

TYPE WRAPPED_FIELDS
  CLASS (FIELD_3RB), POINTER :: F_SPSCALARS
  CLASS (FIELD_2RB), POINTER :: F_SPSCALARS2
  CLASS (FIELD_2RB), POINTER :: F_SPVOR, F_SPDIV

  CLASS (FIELD_3RB), POINTER :: F_VOR, F_DIV
  CLASS (FIELD_3RB), POINTER :: F_U, F_V
  CLASS (FIELD_3RB), POINTER :: F_UDM, F_VDM

  CLASS (FIELD_4RB), POINTER :: F_SCALARS
  CLASS (FIELD_4RB), POINTER :: F_SCALARS_EW
  CLASS (FIELD_4RB), POINTER :: F_SCALARS_NS

  CLASS (FIELD_3RB), POINTER :: F_SCALARS2
  CLASS (FIELD_3RB), POINTER :: F_SCALARS2_EW
  CLASS (FIELD_3RB), POINTER :: F_SCALARS2_NS
END TYPE WRAPPED_FIELDS

TYPE FIELDS_LISTS
  TYPE (FIELD_BASIC_PTR), POINTER:: U (:), V (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SCALAR (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SPVOR (:), SPDIV (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: VOR (:), DIV (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SPSCALAR (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: UDM (:), VDM (:)
  TYPE (FIELD_BASIC_PTR), POINTER:: SCALARDM (:), SCALARDL (:)
END TYPE FIELDS_LISTS

! Number of points in top/bottom latitudes
integer(kind=jpim), parameter :: min_octa_points = 20

integer(kind=jpim) :: istack, getstackusage
real(kind=jprd), dimension(1) :: zmaxerr(5), zerr(5)
real(kind=jprd) :: zmaxerrg

! Output unit numbers
integer(kind=jpim), parameter :: nerr     = 0 ! Unit number for STDERR
integer(kind=jpim), parameter :: nout     = 6 ! Unit number for STDOUT
integer(kind=jpim), parameter :: noutdump = 7 ! Unit number for field output

! Default parameters
integer(kind=jpim) :: nsmax   = 79  ! Spectral truncation
integer(kind=jpim) :: iters   = 10  ! Number of iterations for transform test
integer(kind=jpim) :: nfld    = 1   ! Number of scalar fields 
integer(kind=jpim) :: nlev    = 1   ! Number of vertical levels
integer(kind=jpim) :: iters_warmup = 3 ! Number of warm up steps (for which timing statistics should be ignored)

integer(kind=jpim) :: nflevg
integer(kind=jpim) :: ndgl ! Number of latitudes
integer(kind=jpim) :: nspec2
integer(kind=jpim) :: ngptot
integer(kind=jpim) :: ngptotg
integer(kind=jpim) :: ifld
integer(kind=jpim) :: jroc
integer(kind=jpim) :: jb
integer(kind=jpim) :: nspec2g
integer(kind=jpim) :: i
integer(kind=jpim) :: ja
integer(kind=jpim) :: ib
integer(kind=jpim) :: jprtrv

integer(kind=jpim), allocatable :: nloen(:), nprcids(:)
integer(kind=jpim) :: myproc, jj
integer :: jstep

real(kind=jprd) :: ztinit, ztloop, timef, ztstepmax, ztstepmin, ztstepavg, ztstepmed
real(kind=jprd) :: ztstepmax1, ztstepmin1, ztstepavg1, ztstepmed1
real(kind=jprd) :: ztstepmax2, ztstepmin2, ztstepavg2, ztstepmed2
real(kind=jprd), allocatable :: ztstep(:), ztstep1(:), ztstep2(:)

real(kind=jprb), allocatable :: znormsp(:), znormsp1(:), znormdiv(:), znormdiv1(:)
real(kind=jprb), allocatable :: znormvor(:), znormvor1(:), znormt(:), znormt1(:)
real(kind=jprd) :: zaveave(0:jpmaxstat)

! Grid-point space data structures
real(kind=jprb), allocatable, target PINNED_TAG :: zgmv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), allocatable, target PINNED_TAG :: zgmvs  (:,:,:)   ! Single level fields at t and t-dt
real(kind=jprb), pointer :: zgp3a (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgpuv   (:,:,:,:) ! Multilevel fields at t and t-dt
real(kind=jprb), pointer :: zgp2 (:,:,:) ! Single level fields at t and t-dt

! Spectral space data structures
real(kind=jprb), allocatable, target PINNED_TAG :: sp3d(:,:,:)
real(kind=jprb), pointer :: zspvor(:,:) => null()
real(kind=jprb), pointer :: zspdiv(:,:) => null()
real(kind=jprb), pointer :: zspsc3a(:,:,:) => null()
real(kind=jprb), allocatable PINNED_TAG :: zspsc2(:,:)

TYPE(WRAPPED_FIELDS), ALLOCATABLE :: wf
TYPE(FIELDS_LISTS), ALLOCATABLE :: ylf

logical :: lstack = .false. ! Output stack info
logical :: luserpnm = .false.
logical :: lkeeprpnm = .false.
logical :: luseflt = .false. ! Use fast legendre transforms
logical :: ltrace_stats = .false.
logical :: lstats_omp = .false.
logical :: lstats_comms = .false.
logical :: lstats_mpl = .false.
logical :: lstats = .true. ! gstats statistics
logical :: lbarrier_stats = .false.
logical :: lbarrier_stats2 = .false.
logical :: ldetailed_stats = .false.
logical :: lstats_alloc = .false.
logical :: lsyncstats = .false.
logical :: lstatscpu = .false.
logical :: lstats_mem = .false.
logical :: lxml_stats = .false.
logical :: lvordiv = .false.
logical :: lscders = .false.
logical :: luvders = .false.
logical :: lprint_norms = .false. ! Calculate and print spectral norms
logical :: lmeminfo = .false. ! Show information from FIAT routine ec_meminfo at the end
logical :: LLACC = .false.
logical :: lfield_api = .false.

integer(kind=jpim) :: nstats_mem = 0
integer(kind=jpim) :: ntrace_stats = 0
integer(kind=jpim) :: nprnt_stats = 1
integer(kind=jpim) :: nopt_mem_tr = 0

character*64 :: clfile

! The multiplier of the machine epsilon used as a tolerance for correctness checking
! ncheck = 0 (the default) means that correctness checking is disabled
integer(kind=jpim) :: ncheck = 0

logical :: lmpoff = .false. ! Message passing switch

! Verbosity level (0 or 1)
integer :: verbosity = 0

real(kind=jprd) :: zra = 6371229._jprd

integer(kind=jpim) :: nmax_resol = 37 ! Max number of resolutions
integer(kind=jpim) :: npromatr = 0 ! nproma for trans lib

integer(kind=jpim) :: nproc ! Number of procs
integer(kind=jpim) :: nthread
integer(kind=jpim) :: nprgpns ! Grid-point decomp
integer(kind=jpim) :: nprgpew ! Grid-point decomp
integer(kind=jpim) :: nprtrv = 0 ! Spectral decomp
integer(kind=jpim) :: nprtrw = 0 ! Spectral decomp
integer(kind=jpim) :: mysetv
integer(kind=jpim) :: mysetw
integer(kind=jpim) :: mp_type = 2 ! Message passing type
integer(kind=jpim) :: mbx_size = 150000000 ! Mailbox size

integer(kind=jpim), allocatable :: numll(:), ivset(:)
integer(kind=jpim) :: ivsetsc(1)

integer(kind=jpim) :: nflevl

! sumpini
integer(kind=jpim) :: isqr
logical :: lsync_trans = .true. ! Activate barrier sync
logical :: leq_regions = .true. ! Eq regions flag


integer(kind=jpim) :: nproma = 0
integer(kind=jpim) :: ngpblks
! locals
integer(kind=jpim) :: iprtrv
integer(kind=jpim) :: iprtrw
integer(kind=jpim) :: iprused, ilevpp, irest, ilev, jlev

integer(kind=jpim) :: ndimgmv  = 0 ! Third dim. of gmv "(nproma,nflevg,ndimgmv,ngpblks)"
integer(kind=jpim) :: ndimgmvs = 0 ! Second dim. gmvs "(nproma,ndimgmvs,ngpblks)"

integer(kind=jpim) :: jbegin_uv = 0
integer(kind=jpim) :: jend_uv   = 0
integer(kind=jpim) :: jbegin_sc = 0
integer(kind=jpim) :: jend_sc   = 0
integer(kind=jpim) :: jbegin_scder_NS = 0
integer(kind=jpim) :: jend_scder_NS = 0
integer(kind=jpim) :: jbegin_scder_EW = 0
integer(kind=jpim) :: jend_scder_EW = 0
integer(kind=jpim) :: jbegin_uder_EW = 0
integer(kind=jpim) :: jend_uder_EW = 0
integer(kind=jpim) :: jbegin_vder_EW = 0
integer(kind=jpim) :: jend_vder_EW = 0
integer(kind=jpim) :: iend = 0 

TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: U (:), V (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: SCALAR (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: SPVOR (:), SPDIV (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: VOR (:), DIV (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: SPSCALAR (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: UDM (:), VDM (:)
TYPE (FIELD_BASIC_PTR), ALLOCATABLE,TARGET:: SCALARDM (:), SCALARDL (:)

logical :: ldump_values = .false.
logical :: ldump_crc = .false.

integer, external :: ec_mpirank
logical :: luse_mpi = .true.

character(len=16) :: cgrid = ''

integer(kind=jpim) :: ierr

real(kind=jprb), allocatable :: global_field(:,:)

!===================================================================================================
#include "fspgl_intf.h"
#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "inv_trans_field_api.h"
#include "trans_inq.h"
#include "gath_grid.h"
#include "specnorm.h"
#include "abor1.intfb.h"
#include "gstats_setup.intfb.h"
#include "ec_meminfo.intfb.h"
#include "trans_end.h"

!===================================================================================================

luse_mpi = detect_mpirun()

! Setup
call get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, lscders, luvders, &
  & luseflt, nopt_mem_tr, nproma, verbosity, ldump_values, ldump_crc, lprint_norms, lmeminfo, nprtrv, nprtrw, ncheck,lfield_api)
if (cgrid == '') cgrid = cubic_octahedral_gaussian_grid(nsmax)
call parse_grid(cgrid, ndgl, nloen)
nflevg = nlev

!===================================================================================================

if (luse_mpi) then
  call mpl_init(ldinfo=(verbosity>=1))
  nproc  = mpl_nproc()
  myproc = mpl_myrank()
else
  nproc = 1
  myproc = 1
  mpl_comm = -1
  lsync_trans = .false.
endif
nthread = oml_max_threads()

call dr_hook_init()

!===================================================================================================

if( lstats ) call gstats(0,0)
ztinit = timef()

! only output to stdout on pe 1
if (nproc > 1) then
  if (myproc /= 1) then
    open(unit=nout, file='/dev/null')
  endif
endif

if (ldetailed_stats) then
  lstats_omp    = .true.
  lstats_comms  = .true.
  lstats_mpl    = .true.
  lstatscpu     = .true.
  nprnt_stats   = nproc
!  lstats_mem   = .true.
!  lstats_alloc = .true.
endif

!===================================================================================================

allocate(nprcids(nproc))
do jj = 1, nproc
  nprcids(jj) = jj
enddo

if (nproc <= 1) then
  lmpoff = .true.
endif

! Compute nprgpns and nprgpew
! This version selects most square-like distribution
! These will change if leq_regions=.true.
if (nproc == 0) nproc = 1
isqr = int(sqrt(real(nproc,jprb)))
do ja = isqr, nproc
  ib = nproc/ja
  if (ja*ib == nproc) then
    nprgpns = max(ja,ib)
    nprgpew = min(ja,ib)
    exit
  endif
enddo

! Compute nprtrv and nprtrw if not provided on the command line
if (nprtrv > 0 .or. nprtrw > 0) then
  if (nprtrv == 0) nprtrv = nproc/nprtrw
  if (nprtrw == 0) nprtrw = nproc/nprtrv
  if (nprtrw*nprtrv /= nproc) call abor1('ectrans_benchmark:nprtrw*nprtrv /= nproc')
else
  do jprtrv = 4, nproc
    nprtrv = jprtrv
    nprtrw = nproc/nprtrv
    if (nprtrv*nprtrw /= nproc) cycle
    if (nprtrv > nprtrw) exit
  enddo
  ! Go for approx square partition for backup
  if (nprtrv*nprtrw /= nproc .or. nprtrv > nprtrw) then
    isqr = int(sqrt(real(nproc,jprb)))
    do ja = isqr, nproc
      ib = nproc/ja
      if (ja*ib == nproc) then
        nprtrw = max(ja, ib)
        nprtrv = min(ja, ib)
        exit
      endif
    enddo
  endif
endif

! Create communicators for mpi groups
if (.not.lmpoff) then
  call mpl_groups_create(nprtrw, nprtrv)
endif

if (lmpoff) then
  mysetw = (myproc - 1)/nprtrv + 1
  mysetv = mod(myproc - 1, nprtrv) + 1
else
  call mpl_cart_coords(myproc, mysetw, mysetv)

  ! Just checking for now...
  iprtrv = mod(myproc - 1, nprtrv) + 1
  iprtrw = (myproc - 1)/nprtrv + 1
  if (iprtrv /= mysetv .or. iprtrw /= mysetw) then
    call abor1('ectrans_benchmark:inconsistency when computing mysetw and mysetv')
  endif
endif

if (.not. lmpoff) then
  call mpl_buffer_method(kmp_type=mp_type, kmbx_size=mbx_size, kprocids=nprcids, ldinfo=(verbosity>=1))
endif

! Determine the number of levels attributed to each member of the V set
allocate(numll(nprtrv))
iprused = min(nflevg+1, nprtrv)
ilevpp = nflevg/nprtrv
irest = nflevg -ilevpp*nprtrv
do jroc = 1, nprtrv
  if (jroc <= irest) then
    numll(jroc) = ilevpp+1
  else
    numll(jroc) = ilevpp
  endif
enddo

nflevl = numll(mysetv)

ivsetsc(1) = iprused
ifld = 0

!===================================================================================================
! Setup gstats
!===================================================================================================

if (lstats) then
  call gstats_setup(nproc, myproc, nprcids,                                            &
    & lstats, lstatscpu, lsyncstats, ldetailed_stats, lbarrier_stats, lbarrier_stats2, &
    & lstats_omp, lstats_comms, lstats_mem, nstats_mem, lstats_alloc,                  &
    & ltrace_stats, ntrace_stats, nprnt_stats, lxml_stats)
  call gstats_psut

  ! Assign labels to GSTATS regions
  call gstats_labels
endif

!===================================================================================================
! Call ecTrans setup routines
!===================================================================================================

if (verbosity >= 1) write(nout,'(a)')'======= Setup ecTrans ======='

call gstats(1, 0)
call setup_trans0(kout=nout, kerr=nerr, kprintlev=merge(2, 0, verbosity == 1),                 &
  &               kmax_resol=nmax_resol, kpromatr=npromatr, kprgpns=nprgpns, kprgpew=nprgpew,  &
  &               kprtrw=nprtrw, ldsync_trans=lsync_trans,                                     &
  &               ldeq_regions=leq_regions, prad=zra, ldalloperm=.true., ldmpoff=.not.luse_mpi,&
  &               kopt_memory_tr=nopt_mem_tr)
call gstats(1, 1)

call gstats(2, 0)
! IFS spectral fields are dimensioned NFLEVL, Nils !!
call set_ectrans_gpu_nflev(nflevl)
  ! We pass nflevl via environment variable in order not to change API
  ! In long run, ectrans should grow its internal buffers automatically
call setup_trans(ksmax=nsmax, kdgl=ndgl, kloen=nloen, ldsplit=.true.,       &
  &              lduserpnm=luserpnm, ldkeeprpnm=lkeeprpnm, &
  &              lduseflt=luseflt)
call gstats(2, 1)

call trans_inq(kspec2=nspec2, kspec2g=nspec2g, kgptot=ngptot, kgptotg=ngptotg)

if (nproma == 0) then ! no blocking (default when not specified)
  nproma = ngptot
endif

! Calculate number of NPROMA blocks
ngpblks = (ngptot - 1)/nproma+1

!===================================================================================================
! Print information before starting
!===================================================================================================

! Print configuration details
if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a)')'======= Start of runtime parameters ======='
  write(nout,'(" ")')
  write(nout,'("nsmax      ",i0)') nsmax
  write(nout,'("grid       ",a)') trim(cgrid)
  write(nout,'("ndgl       ",i0)') ndgl
  write(nout,'("nproc      ",i0)') nproc
  write(nout,'("nthread    ",i0)') nthread
  write(nout,'("nprgpns    ",i0)') nprgpns
  write(nout,'("nprgpew    ",i0)') nprgpew
  write(nout,'("nprtrw     ",i0)') nprtrw
  write(nout,'("nprtrv     ",i0)') nprtrv
  write(nout,'("ngptot     ",i0)') ngptot
  write(nout,'("ngptotg    ",i0)') ngptotg
  write(nout,'("nfld       ",i0)') nfld
  write(nout,'("nlev       ",i0)') nlev
  write(nout,'("nproma     ",i0)') nproma
  write(nout,'("ngpblks    ",i0)') ngpblks
  write(nout,'("nspec2     ",i0)') nspec2
  write(nout,'("nspec2g    ",i0)') nspec2g
  write(nout,'("luseflt    ",l1)') luseflt
  write(nout,'("nopt_mem_tr",i0)') nopt_mem_tr
  write(nout,'("lvordiv    ",l1)') lvordiv
  write(nout,'("lscders    ",l1)') lscders
  write(nout,'("luvders    ",l1)') luvders
  write(nout,'("lfield_api ",l1)') lfield_api
  write(nout,'(" ")')
  write(nout,'(a)') '======= End of runtime parameters ======='
  write(nout,'(" ")')
end if

!===================================================================================================
! Allocate and Initialize spectral arrays
!===================================================================================================

! Allocate spectral arrays
! Try to mimick IFS layout as much as possible
nullify(zspvor)
nullify(zspdiv)
nullify(zspsc3a)
allocate(sp3d(nflevl,nspec2,2+nfld))
allocate(zspsc2(1,nspec2))

call initialize_spectral_arrays(nsmax, zspsc2, sp3d)

! Point convenience variables to storage variable sp3d
zspvor  => sp3d(:,:,1)
zspdiv  => sp3d(:,:,2)
zspsc3a => sp3d(:,:,3:3+(nfld-1))

!===================================================================================================
! Allocate gridpoint arrays
!===================================================================================================

allocate(ivset(nflevg))

! Compute spectral distribution
ilev = 0
do jb = 1, nprtrv
  do jlev=1, numll(jb)
    ilev = ilev + 1
    ivset(ilev) = jb
  enddo
enddo

! Allocate grid-point arrays
if (lvordiv) then
  jbegin_uv = 1
  jend_uv = 2  
endif
if (luvders) then
  jbegin_uder_EW  = jend_uv + 1
  jend_uder_EW    = jbegin_uder_EW + 1
  jbegin_vder_EW  = jend_uder_EW + 1
  jend_vder_EW    = jbegin_vder_EW + 1
else
  jbegin_uder_EW = jend_uv
  jend_uder_EW   = jend_uv
  jbegin_vder_EW = jend_uv
  jend_vder_EW   = jend_uv
endif

jbegin_sc = jend_vder_EW + 1
jend_sc   = jend_vder_EW + nfld

if (lscders) then
  ndimgmvs = 3
  jbegin_scder_NS = jend_sc + 1
  jend_scder_NS   = jend_sc + nfld
  jbegin_scder_EW = jend_scder_NS + 1
  jend_scder_EW   = jend_scder_NS + nfld
else
  ndimgmvs = 1
  jbegin_scder_NS = jend_sc
  jend_scder_NS   = jend_sc
  jbegin_scder_EW = jend_sc
  jend_scder_EW   = jend_sc
endif

ndimgmv = jend_scder_EW

allocate(zgmv(nproma,nflevg,ndimgmv,ngpblks))
allocate(zgmvs(nproma,ndimgmvs,ngpblks))

zgpuv => zgmv(:,:,1:jend_vder_EW,:)
zgp3a => zgmv(:,:,jbegin_sc:jend_scder_EW,:)
zgp2  => zgmvs(:,:,:)

write(nout,*) "jbegin_uv", jbegin_uv
write(nout,*) "jend_uv", jend_uv
write(nout,*) "jbegin_sc", jbegin_sc
write(nout,*) "jend_sc", jend_sc
write(nout,*) "jbegin_scder_NS", jbegin_scder_NS
write(nout,*) "jend_scder_NS", jend_scder_NS
write(nout,*) "jbegin_scder_EW", jbegin_scder_EW
write(nout,*) "jend_scder_EW", jend_scder_EW
write(nout,*) "jbegin_uder_EW", jbegin_uder_EW
write(nout,*) "jend_uder_EW", jend_uder_EW
write(nout,*) "jbegin_vder_EW", jbegin_vder_EW
write(nout,*) "jend_vder_EW", jend_vder_EW

write(nout,*) "sp3d", shape(sp3d)
write(nout,*) "zspsc2", shape(zspsc2)
write(nout,*) "zgmv", shape(zgmv)
write(nout,*) "zgmvs", shape(zgmvs)
write(nout,*) "zgp2", shape(zgp2)

if (lfield_api) then
wf = wrap_fields(lvordiv, lscders, luvders, &
                & sp3d, zspsc2, zgmv, zgmvs, zgp2, &
                & jbegin_uv,jend_uv, &
                & jbegin_sc,jend_sc, &
                & jbegin_scder_NS, jend_scder_NS, &
                & jbegin_scder_EW, jend_scder_EW, &
                & jbegin_uder_EW, jend_uder_EW, &
                & jbegin_vder_EW, jend_vder_EW)
call output_wrapped_file(wf)

ylf = create_fields_lists(wf,ivset,ivsetsc)
call output_fields_lists(ylf)
endif
!===================================================================================================
! Allocate norm arrays
!===================================================================================================

if (lprint_norms .or. ncheck > 0) then
  allocate(znormsp(1))
  allocate(znormsp1(1))
  allocate(znormvor(nflevg))
  allocate(znormvor1(nflevg))
  allocate(znormdiv(nflevg))
  allocate(znormdiv1(nflevg))
  allocate(znormt(nflevg))
  allocate(znormt1(nflevg))

  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor1, kvset=ivset(1:nflevg))
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv1, kvset=ivset(1:nflevg))
  if (nfld > 0) then
    call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt1,   kvset=ivset(1:nflevg))
  endif
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp1,  kvset=ivsetsc)

  if (verbosity >= 1 .and. myproc == 1) then
    do ifld = 1, nflevg
      write(nout,'("norm zspvor( ",i4,",:)   = ",f20.15)') ifld, znormvor1(ifld)
      write(nout,'("0x",Z16.16)') znormvor1(ifld)
    enddo
    do ifld = 1, nflevg
      write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15)') ifld, znormdiv1(ifld)
      write(nout,'("0x",Z16.16)') znormdiv1(ifld)
    enddo
    if (nfld > 0) then
      do ifld = 1, nflevg
        write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15)') ifld, znormt1(ifld)
        write(nout,'("0x",Z16.16)') znormt1(ifld)
      enddo
    endif
    do ifld = 1, 1
      write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15)') ifld, znormsp1(ifld)
      write(nout,'("0x",Z16.16)') znormsp1(ifld)
    enddo
  endif
endif

!===================================================================================================
! Setup timers
!===================================================================================================

ztinit = (timef() - ztinit)/1000.0_jprd

if (verbosity >= 0 .and. myproc == 1) then
  write(nout,'(" ")')
  write(nout,'(a,i0,a,f9.2,a)') "ectrans_benchmark initialisation, on ",nproc,&
                                & " tasks, took ",ztinit," sec"
  write(nout,'(" ")')
endif

if (iters <= 0) call abor1('ectrans_benchmark:iters <= 0')

allocate(ztstep(iters+iters_warmup))
allocate(ztstep1(iters+iters_warmup))
allocate(ztstep2(iters+iters_warmup))

if (verbosity >= 1 .and. myproc == 1) then
  write(nout,'(a)') '======= Start of spectral transforms  ======='
  write(nout,'(" ")')
endif



!===================================================================================================
! Do spectral transform loop
!===================================================================================================

gstats_lstats = .false.

write(nout,'(a,i0,a,i0,a)') 'Running for ', iters, ' iterations with ', iters_warmup, &
  & ' extra warm-up iterations'
write(nout,'(" ")')

do jstep = 1, iters+iters_warmup
  if (jstep == iters_warmup + 1) then
    gstats_lstats = .true.
    ztloop = timef()
  endif

  call gstats(3,0)
  ztstep(jstep) = timef()

  !=================================================================================================
  ! Do inverse transform
  !=================================================================================================

  ztstep1(jstep) = timef()
  call gstats(4,0)

  if (ldump_crc) then
    zgmv(:,:,:,:) = 0
    zgmvs(:,:,:) = 0
  endif

  if (lfield_api) then
    
    CALL INV_TRANS_FIELD_API (YDFSPVOR=ylf%SPVOR, YDFSPDIV=ylf%SPDIV, YDFSPSCALAR=ylf%SPSCALAR, &
                            & YDFU=ylf%U, YDFV=ylf%V, YDFSCALAR=ylf%SCALAR, &
                            & YDFUDM=ylf%UDM, YDFVDM=ylf%VDM, &
                            & YDFSCALARDM=ylf%SCALARDM, YDFSCALARDL=ylf%SCALARDL, &
                            & YDFVOR=ylf%VOR, YDFDIV=ylf%DIV, &
                            & KSPEC=NSPEC2, KPROMA=NPROMA, KGPBLKS=NGPBLKS, KGPTOT=NGPTOT, KFLEVG=NFLEVG, KFLEVL=NFLEVL,& 
                            & LDACC=LLACC, LDVERBOSE =.TRUE.)
else
  if (lvordiv) then
    call inv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspvor=zspvor,                     & ! spectral vorticity
       & pspdiv=zspdiv,                     & ! spectral divergence
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  &
       & ldvorgp=.false.,                   & ! no gridpoint vorticity
       & lddivgp=.false.,                   & ! no gridpoint divergence
       & lduvder=luvders,                   &
       & kvsetuv=ivset,                     &
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgpuv=zgpuv,                       &
       & pgp3a=zgp3a)
  else
    call inv_trans(kresol=1, kproma=nproma, &
       & pspsc2=zspsc2,                     & ! spectral surface pressure
       & pspsc3a=zspsc3a,                   & ! spectral scalars
       & ldscders=lscders,                  & ! scalar derivatives
       & kvsetsc2=ivsetsc,                  &
       & kvsetsc3a=ivset,                   &
       & pgp2=zgp2,                         &
       & pgp3a=zgp3a)
  endif  
endif

if (ldump_crc) then  
    if (lfield_api) then        
      write (clfile,'(A)') 'inv_trans_field_api.txt'
    else
      write (clfile,'(A)') 'inv_trans.txt'
      ! Remove trash at end of last block    
      IEND = NGPTOT - NPROMA * (NGPBLKS - 1)
      !iend = nproma * ngpblks - ngptot      
      zgmvs (iend+1:, :, ngpblks) = 0
    endif
    
    call dump_crc(clfile, iter=jstep,zgmv=zgmv, zgmvs=zgmvs)
endif

  call gstats(4,1)

  ztstep1(jstep) = (timef() - ztstep1(jstep))/1000.0_jprd

  !=================================================================================================
  ! While in grid point space, dump the values to disk, for debugging only
  !=================================================================================================

  if (ldump_values .and. mod(jstep,10) == 1) then
    if (myproc == 1) then
      allocate(global_field(ngptotg,1))
    endif
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp2(:,1:1,:), 's', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,1,:), 'u', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgpuv(:,nflevg:nflevg,2,:), 'v', noutdump)
    call dump_gridpoint_field(jstep, myproc, nproma, global_field, zgp3a(:,nflevg:nflevg,1,:), 't', noutdump)
    if (myproc == 1) then
      deallocate(global_field)
    endif
  endif

  !=================================================================================================
  ! Do direct transform
  !=================================================================================================

  ztstep2(jstep) = timef()

  call gstats(5,0)
  if (lvordiv) then
    call dir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgmvs(:,1:1,:),                &
      & pgpuv=zgpuv(:,:,1:2,:),             &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspvor=zspvor,                      &
      & pspdiv=zspdiv,                      &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetuv=ivset,                      &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset)
  else
    call dir_trans(kresol=1, kproma=nproma, &
      & pgp2=zgmvs(:,1:1,:),                &
      & pgp3a=zgp3a(:,:,1:nfld,:),          &
      & pspsc2=zspsc2,                      &
      & pspsc3a=zspsc3a,                    &
      & kvsetsc2=ivsetsc,                   &
      & kvsetsc3a=ivset)
  endif
  if (ldump_crc) call dump_crc("dir_trans.txt", iter=jstep,sp3d=sp3d,zspc2=zspsc2)
  call gstats(5,1)

  ztstep2(jstep) = (timef() - ztstep2(jstep))/1000.0_jprd

  ztstep(jstep) = (timef() - ztstep(jstep))/1000.0_jprd
  
  !=================================================================================================
  ! Print norms
  !=================================================================================================

  if (lprint_norms) then
    call gstats(6,0)
    call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc(1:1))
    call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset(1:nflevg))
    call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset(1:nflevg))
    if (nfld > 0) then
      call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset(1:nflevg))
    endif

    ! Surface pressure
    if (myproc == 1) then
      zmaxerr(:) = -999.0
      do ifld = 1, 1
        zerr(1) = abs(znormsp1(ifld)/znormsp(ifld) - 1.0_jprb)
        zmaxerr(1) = max(zmaxerr(1), zerr(1))
      enddo
      ! Divergence
      do ifld = 1, nflevg
        zerr(2) = abs(znormdiv1(ifld)/znormdiv(ifld) - 1.0_jprb)
        zmaxerr(2) = max(zmaxerr(2), zerr(2))
      enddo
      ! Vorticity
      do ifld = 1, nflevg
        zerr(3) = abs(znormvor1(ifld)/znormvor(ifld) - 1.0_jprb)
        zmaxerr(3) = max(zmaxerr(3),zerr(3))
      enddo
      ! Temperature
      if (nfld > 0) then
        do ifld = 1, nflevg
          zerr(4) = abs(znormt1(ifld)/znormt(ifld) - 1.0_jprb)
          zmaxerr(4) = max(zmaxerr(4), zerr(4))
        enddo
        write(nout,'("time step ",i6," took", f8.4," | zspvor max err=",e10.3,&
                    & " | zspdiv max err=",e10.3," | zspsc3a max err=",e10.3," | zspsc2 max err=",e10.3)') &
                    &  jstep, ztstep(jstep), zmaxerr(3), zmaxerr(2), zmaxerr(4), zmaxerr(1)
      else
        write(nout,'("time step ",i6," took", f8.4," | zspvor max err=",e10.3,&
                    & " | zspdiv max err=",e10.3," | zspsc2 max err=",e10.3)') &
                    &  jstep, ztstep(jstep), zmaxerr(3), zmaxerr(2), zmaxerr(1)
      endif
    endif
    call gstats(6,1)
  else
    write(nout,'("Time step ",i6," took", f8.4)') jstep, ztstep(jstep)
  endif
  call gstats(3,1)
enddo

!===================================================================================================

ztloop = (timef() - ztloop)/1000.0_jprd

write(nout,'(" ")')
write(nout,'(a)') '======= End of spectral transforms  ======='
write(nout,'(" ")')

if (lprint_norms .or. ncheck > 0) then
  call specnorm(pspec=zspvor(1:nflevl,:),    pnorm=znormvor, kvset=ivset)
  call specnorm(pspec=zspdiv(1:nflevl,:),    pnorm=znormdiv, kvset=ivset)
  if (nfld > 0) then
    call specnorm(pspec=zspsc3a(1:nflevl,:,1), pnorm=znormt,   kvset=ivset)
  endif
  call specnorm(pspec=zspsc2(1:1,:),         pnorm=znormsp,  kvset=ivsetsc)

  if (myproc == 1) then
    zmaxerr(:) = -999.0
    do ifld = 1, nflevg
      zerr(3) = abs(real(znormvor1(ifld),kind=jprd)/real(znormvor(ifld),kind=jprd) - 1.0_jprd)
      zmaxerr(3) = max(zmaxerr(3), zerr(3))
      if (verbosity >= 1) then
        write(nout,'("norm zspvor( ",i4,")     = ",f20.15,"        error = ",e10.3)') ifld, znormvor(ifld), zerr(3)
        write(nout,'("0x",Z16.16)') znormvor(ifld)
      endif
    enddo
    do ifld = 1, nflevg
      zerr(2) = abs(real(znormdiv1(ifld),kind=jprd)/real(znormdiv(ifld),kind=jprd) - 1.0d0)
      zmaxerr(2) = max(zmaxerr(2),zerr(2))
      if (verbosity >= 1) then
        write(nout,'("norm zspdiv( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormdiv(ifld), zerr(2)
        write(nout,'("0x",Z16.16)') znormdiv(ifld)
      endif
    enddo
    if (nfld > 0) then
      do ifld = 1, nflevg
        zerr(4) = abs(real(znormt1(ifld),kind=jprd)/real(znormt(ifld),kind=jprd) - 1.0d0)
        zmaxerr(4) = max(zmaxerr(4), zerr(4))
        if (verbosity >= 1) then
          write(nout,'("norm zspsc3a(",i4,",:,1) = ",f20.15,"        error = ",e10.3)') ifld, znormt(ifld), zerr(4)
          write(nout,'("0x",Z16.16)') znormt(ifld)
        endif
      enddo
    endif
    do ifld = 1, 1
      zerr(1) = abs(real(znormsp1(ifld),kind=jprd)/real(znormsp(ifld),kind=jprd) - 1.0d0)
      zmaxerr(1) = max(zmaxerr(1), zerr(1))
      if (verbosity >= 1) then
        write(nout,'("norm zspsc2( ",i4,",:)   = ",f20.15,"        error = ",e10.3)') ifld, znormsp(ifld), zerr(1)
        write(nout,'("0x",Z16.16)') znormsp(ifld)
      endif
    enddo

    ! maximum error across all fields
    if (nfld > 0) then
      zmaxerrg = max(zmaxerr(1), zmaxerr(2), zmaxerr(3), zmaxerr(4))
    else
      zmaxerrg = max(zmaxerr(1), zmaxerr(2), zmaxerr(3))
    endif

    if (verbosity >= 1) write(nout,*)
    write(nout,'("max error zspvor(1:nlev,:)    = ",e10.3)') zmaxerr(3)
    write(nout,'("max error zspdiv(1:nlev,:)    = ",e10.3)') zmaxerr(2)
    if (nfld > 0) write(nout,'("max error zspsc3a(1:nlev,:,1) = ",e10.3)') zmaxerr(4)
    write(nout,'("max error zspsc2(1:1,:)       = ",e10.3)') zmaxerr(1)
    write(nout,*)
    write(nout,'("max error combined =          = ",e10.3)') zmaxerrg
    write(nout,*)
  endif
  if (ncheck > 0) then
    ierr = 0
    if (myproc == 1) then
      ! If the maximum spectral norm error across all fields is greater than 100 times the machine
      ! epsilon, fail the test
      if (zmaxerrg > real(ncheck, jprb) * epsilon(1.0_jprb)) then
        write(nout, '(a)') '*******************************'
        write(nout, '(a)') 'Correctness test failed'
        write(nout, '(a,1e9.2)') 'Maximum spectral norm error = ', zmaxerrg
        write(nout, '(a,1e9.2)') 'Error tolerance = ', real(ncheck, jprb) * epsilon(1.0_jprb)
        write(nout, '(a)') '*******************************'
        ierr = 1
      endif
    endif

    ! Root rank broadcasts the correctness checker result to the other ranks
    if (luse_mpi) then
      call mpl_broadcast(ierr,kroot=1,ktag=1)
    endif

    ! Halt if correctness checker failed
    if (ierr == 1) then
      error stop
    endif
  endif
endif

!===================================================================================================
! Calculate timings
!===================================================================================================

ztstepavg = sum(ztstep(iters_warmup+1:))
ztstepmin = minval(ztstep(iters_warmup+1:))
ztstepmax = maxval(ztstep(iters_warmup+1:))
ztstepavg1 = sum(ztstep1(iters_warmup+1:))
ztstepmin1 = minval(ztstep1(iters_warmup+1:))
ztstepmax1 = maxval(ztstep1(iters_warmup+1:))
ztstepavg2 = sum(ztstep2(iters_warmup+1:))
ztstepmin2 = minval(ztstep2(iters_warmup+1:))
ztstepmax2 = maxval(ztstep2(iters_warmup+1:))

if (luse_mpi) then
  call mpl_allreduce(ztloop,     'sum', ldreprod=.false.)
  call mpl_allreduce(ztstep,     'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg,  'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax,  'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin,  'min', ldreprod=.false.)

  call mpl_allreduce(ztstep1,    'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg1, 'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax1, 'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin1, 'min', ldreprod=.false.)

  call mpl_allreduce(ztstep2,    'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepavg2, 'sum', ldreprod=.false.)
  call mpl_allreduce(ztstepmax2, 'max', ldreprod=.false.)
  call mpl_allreduce(ztstepmin2, 'min', ldreprod=.false.)
endif

ztstepavg = (ztstepavg/real(nproc,jprb))/real(iters,jprd)
ztloop = ztloop/real(nproc,jprd)
ztstep(:) = ztstep(:)/real(nproc,jprd)
ztstepmed = get_median(ztstep(iters_warmup+1:))

ztstepavg1 = (ztstepavg1/real(nproc,jprb))/real(iters,jprd)
ztstep1(:) = ztstep1(:)/real(nproc,jprd)
ztstepmed1 = get_median(ztstep1(iters_warmup+1:))

ztstepavg2 = (ztstepavg2/real(nproc,jprb))/real(iters,jprd)
ztstep2(:) = ztstep2(:)/real(nproc,jprd)
ztstepmed2 = get_median(ztstep2(iters_warmup+1:))

write(nout,'(a)') '======= Start of time step stats ======='
write(nout,'(" ")')
write(nout,'("Inverse transforms")')
write(nout,'("------------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg1
write(nout,'("min  (s): ",f8.4)') ztstepmin1
write(nout,'("max  (s): ",f8.4)') ztstepmax1
write(nout,'("med  (s): ",f8.4)') ztstepmed1
write(nout,'(" ")')
write(nout,'("Direct transforms")')
write(nout,'("-----------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg2
write(nout,'("min  (s): ",f8.4)') ztstepmin2
write(nout,'("max  (s): ",f8.4)') ztstepmax2
write(nout,'("med  (s): ",f8.4)') ztstepmed2
write(nout,'(" ")')
write(nout,'("Inverse-direct transforms")')
write(nout,'("-------------------------")')
write(nout,'("avg  (s): ",f8.4)') ztstepavg
write(nout,'("min  (s): ",f8.4)') ztstepmin
write(nout,'("max  (s): ",f8.4)') ztstepmax
write(nout,'("med  (s): ",f8.4)') ztstepmed
write(nout,'("loop (s): ",f8.4)') ztloop
write(nout,'(" ")')
write(nout,'(a)') '======= End of time step stats ======='
write(nout,'(" ")')

if (lstack) then
  ! Gather stack usage statistics
  istack = getstackusage()
  if (myproc == 1) then
    print 9000, istack
    9000 format("Stack utilisation information",/,&
         &"=============================",//,&
         &"Task           size(bytes)",/,&
         &"====           ===========",//,&
         &"   1",11x,i10)

    do i = 2, nproc
      call mpl_recv(istack, ksource=nprcids(i), ktag=i, cdstring='ectrans_benchmark:')
      print '(i4,11x,i10)', i, istack
    enddo
  else
    call mpl_send(istack, kdest=nprcids(1), ktag=myproc, cdstring='ectrans_benchmark:')
  endif
endif


!===================================================================================================
! Cleanup
!===================================================================================================

deallocate(zgmv)
deallocate(zgmvs)

if (lfield_api) then
  call delete_wrapped_fields(wf)
  call delete_fields_lists(ylf)
endif
!===================================================================================================

if (lstats) then
  call gstats(0,1)
  call gstats_print(nout, zaveave, jpmaxstat)
endif

if (lmeminfo) then
  write(nout,*)
  call ec_meminfo(nout, "", mpl_comm, kbarr=1, kiotask=-1, &
      & kcall=1)
endif

call trans_end

!===================================================================================================
! Finalize MPI
!===================================================================================================

if (luse_mpi) then
  call mpl_end(ldmeminfo=.false.)
endif

!===================================================================================================
! Close file
!===================================================================================================

if (nproc > 1) then
  if (myproc /= 1) then
    close(unit=nout)
  endif
endif

!===================================================================================================

contains

!===================================================================================================

subroutine parse_grid(cgrid,ndgl,nloen)

  character(len=*) :: cgrid
  integer, intent(inout) :: ndgl
  integer, intent(inout), allocatable :: nloen(:)
  integer :: ios
  integer :: gaussian_number
  read(cgrid(2:len_trim(cgrid)),*,IOSTAT=ios) gaussian_number
  if (ios==0) then
    ndgl = 2 * gaussian_number
    allocate(nloen(ndgl))
    if (cgrid(1:1) == 'F') then ! Regular Gaussian grid
      nloen(:) = gaussian_number * 4
      return
    endif
    if (cgrid(1:1) == 'O') then ! Octahedral Gaussian grid
      do i = 1, ndgl / 2
        nloen(i) = 20 + 4 * (i - 1)
        nloen(ndgl - i + 1) = nloen(i)
      end do
      return
    endif
  endif
  call parsing_failed("ERROR: Unsupported grid specified: "// trim(cgrid))

end subroutine

!===================================================================================================

subroutine str2int(str, int, stat)

  character(len=*), intent(in) :: str
  integer, intent(out) :: int
  integer, intent(out) :: stat
  read(str, *, iostat=stat) int

end subroutine str2int

!===================================================================================================

function get_int_value(cname, iarg) result(value)

  integer :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg
  character(len=128) :: carg
  integer :: stat

  carg = get_str_value(cname, iarg)
  call str2int(carg, value, stat)

  if (stat /= 0) then
    call parsing_failed("Invalid argument for " // trim(cname) // ": " // trim(carg))
  end if

end function

!===================================================================================================

function get_str_value(cname, iarg) result(value)

  character(len=128) :: value
  character(len=*), intent(in) :: cname
  integer, intent(inout) :: iarg

  iarg = iarg + 1
  call get_command_argument(iarg, value)

  if (value == "") then
    call parsing_failed("Invalid argument for " // trim(cname) // ": no value provided")
  end if

end function

!===================================================================================================

subroutine print_help(unit)

  integer, optional :: unit
  integer :: nout = 6
  if (present(unit)) then
    nout = unit
  endif

  write(nout, "(a)") ""

  if (jprb == jprd) then
    write(nout, "(a)") "NAME    ectrans-benchmark-" // VERSION // "-dp"
  else
    write(nout, "(a)") "NAME    ectrans-benchmark-" // VERSION // "-sp"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "DESCRIPTION"
  write(nout, "(a)") "        This program tests ecTrans by transforming fields back and forth&
    & between spectral "
  if (jprb == jprd) then
    write(nout, "(a)") "        space and grid-point space (double-precision version)"
  else
    write(nout, "(a)") "        space and grid-point space (single-precision version)"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "USAGE"
  if (jprb == jprd) then
    write(nout, "(a)") "        ectrans-benchmark-" // VERSION // "-dp [options]"
  else
    write(nout, "(a)") "        ectrans-benchmark-" // VERSION // "-sp [options]"
  end if
  write(nout, "(a)") ""

  write(nout, "(a)") "OPTIONS"
  write(nout, "(a)") "    -h, --help          Print this message"
  write(nout, "(a)") "    -v                  Run with verbose output"
  write(nout, "(a)") "    -t, --truncation T  Run with this triangular spectral truncation&
    & (default = 79)"
  write(nout, "(a)") "    -g, --grid GRID     Run with this grid. Possible values: O<N>, F<N>"
  write(nout, "(a)") "                        If not specified, O<N> is used with N=truncation+1&
    & (cubic relation)"
  write(nout, "(a)") "    -n, --niter NITER   Run for this many inverse/direct transform&
    & iterations (default = 10)"
  write(nout, "(a)") "    --niter-warmup      Number of warm up iterations,&
    & for which timing statistics should be ignored (default = 3)"
  write(nout, "(a)") "    -f, --nfld NFLD     Number of scalar fields (default = 1)"
  write(nout, "(a)") "    -l, --nlev NLEV     Number of vertical levels (default = 1)"
  write(nout, "(a)") "    --vordiv            Also transform vorticity-divergence to wind"
  write(nout, "(a)") "    --scders            Compute scalar derivatives (default off)"
  write(nout, "(a)") "    --uvders            Compute uv East-West derivatives (default off). Only&
    & when also --vordiv is given"
  write(nout, "(a)") "    --flt               Run with fast Legendre transforms (default off)"
  write(nout, "(a)") "    --nproma NPROMA     Run with NPROMA (default no blocking: NPROMA=ngptot)"
  write(nout, "(a)") "    --norms             Calculate and print spectral norms of transformed&
    & fields"
  write(nout, "(a)") "                        The computation of spectral norms will skew overall&
    & timings"
  write(nout, "(a)") "    --meminfo           Show diagnostic information from FIAT's ec_meminfo&
    & subroutine on memory usage, thread-binding etc."
  write(nout, "(a)") "    --nprtrv            Size of V set in spectral decomposition"
  write(nout, "(a)") "    --nprtrw            Size of W set in spectral decomposition"
  write(nout, "(a)") "    --field-api         Use field API interface"
  write(nout, "(a)") "    -c, --check VALUE   The multiplier of the machine epsilon used as a&
   & tolerance for correctness checking"
  write(nout, "(a)") ""
  write(nout, "(a)") "DEBUGGING"
  write(nout, "(a)") "    --dump-values       Output gridpoint fields in unformatted binary file"
  write(nout, "(a)") ""

end subroutine print_help

!===================================================================================================

subroutine parsing_failed(message)

  character(len=*), intent(in) :: message
  if (luse_mpi) call mpl_init(ldinfo=.false.)
  if (ec_mpirank() == 0) then
    write(nerr,"(a)") trim(message)
    call print_help(unit=nerr)
  endif
  if (luse_mpi) call mpl_end(ldmeminfo=.false.)
  stop

end subroutine

!===================================================================================================

subroutine get_command_line_arguments(nsmax, cgrid, iters, iters_warmup, nfld, nlev, lvordiv, lscders, luvders, &
  &                                   luseflt, nopt_mem_tr, nproma, verbosity, ldump_values, ldump_crc, lprint_norms, &
  &                                   lmeminfo, nprtrv, nprtrw, ncheck,lfield_api)

#ifdef _OPENACC
  use openacc, only: acc_init, acc_get_device_type
#endif

  integer, intent(inout) :: nsmax           ! Spectral truncation
  character(len=16), intent(inout) :: cgrid ! Spectral truncation
  integer, intent(inout) :: iters           ! Number of iterations for transform test
  integer, intent(inout) :: iters_warmup    ! Number of iterations for transform test
  integer, intent(inout) :: nfld            ! Number of scalar fields
  integer, intent(inout) :: nlev            ! Number of vertical levels
  logical, intent(inout) :: lvordiv         ! Also transform vorticity/divergence
  logical, intent(inout) :: lscders         ! Compute scalar derivatives
  logical, intent(inout) :: luvders         ! Compute uv East-West derivatives
  logical, intent(inout) :: luseflt         ! Use fast Legendre transforms
  integer, intent(inout) :: nopt_mem_tr     ! Use of heap or stack memory for ZCOMBUF arrays in transposition arrays (0 for heap, 1 for stack)
  integer, intent(inout) :: nproma          ! NPROMA
  integer, intent(inout) :: verbosity       ! Level of verbosity
  logical, intent(inout) :: ldump_values    ! Dump values of grid point fields for debugging
  logical, intent(inout) :: ldump_crc   !
  logical, intent(inout) :: lprint_norms    ! Calculate and print spectral norms of fields
  logical, intent(inout) :: lmeminfo        ! Show information from FIAT ec_meminfo routine at the
                                            ! end
  integer, intent(inout) :: nprtrv          ! Size of V set (spectral decomposition)
  integer, intent(inout) :: nprtrw          ! Size of W set (spectral decomposition)
  integer, intent(inout) :: ncheck          ! The multiplier of the machine epsilon used as a
                                            ! tolerance for correctness checking

  logical, intent(inout) :: lfield_api
  character(len=128) :: carg          ! Storage variable for command line arguments
  integer            :: iarg = 1      ! Argument index

#ifdef _OPENACC
  call acc_init(acc_get_device_type())
#endif

  do while (iarg <= command_argument_count())
    call get_command_argument(iarg, carg)

    select case(carg)
      ! Parse help argument
      case('-h', '--help')
        if (luse_mpi) call mpl_init(ldinfo=.false.)
        if (ec_mpirank()==0) call print_help()
        if (luse_mpi) call mpl_end(ldmeminfo=.false.)
        stop
      ! Parse verbosity argument
      case('-v')
        verbosity = 1
      ! Parse number of iterations argument
      case('-n', '--niter')
        iters = get_int_value('-n', iarg)
        if (iters < 1) then
          call parsing_failed("Invalid argument for -n: must be > 0")
        end if
      case('--niter-warmup')
        iters_warmup = get_int_value('--niter-warmup', iarg)
        if (iters_warmup < 0) then
          call parsing_failed("Invalid argument for --niter-warmup: must be >= 0")
        end if
      ! Parse spectral truncation argument
      case('-t', '--truncation')
        nsmax = get_int_value('-t', iarg)
        if (nsmax < 1) then
          call parsing_failed("Invalid argument for -t: must be > 0")
        end if
      case('-g', '--grid'); cgrid = get_str_value('-g', iarg)
      case('-f', '--nfld'); nfld = get_int_value('-f', iarg)
      case('-l', '--nlev'); nlev = get_int_value('-l', iarg)
      case('--vordiv'); lvordiv = .True.
      case('--scders'); lscders = .True.
      case('--uvders'); luvders = .True.
      case('--flt'); luseflt = .True.
      case('--mem-tr'); nopt_mem_tr = get_int_value('--mem-tr', iarg)
      case('--nproma'); nproma = get_int_value('--nproma', iarg)
      case('--dump-values'); ldump_values = .true.
      case('--dump-crc'); ldump_crc = .true.
      case('--norms'); lprint_norms = .true.
      case('--meminfo'); lmeminfo = .true.
      case('--nprtrv'); nprtrv = get_int_value('--nprtrv', iarg)
      case('--nprtrw'); nprtrw = get_int_value('--nprtrw', iarg)
      case('-c', '--check'); ncheck = get_int_value('-c', iarg)
      case('--field-api'); lfield_api = .True.
              
      case default
        call parsing_failed("Unrecognised argument: " // trim(carg))

    end select
    iarg = iarg + 1
  end do

  if (.not. lvordiv) then
    luvders = .false.
  endif

end subroutine get_command_line_arguments

!===================================================================================================

function cubic_octahedral_gaussian_grid(nsmax) result(cgrid)

  character(len=16) :: cgrid
  integer, intent(in) :: nsmax
  write(cgrid,'(a,i0)') 'O',nsmax+1

end function


!===================================================================================================

function get_median(vec) result(median)

  real(kind=jprd), intent(in) :: vec(:)
  real(kind=jprd) :: median

  real(kind=jprd) :: vec_sorted(size(vec))
  real(kind=jprd) :: x

  integer :: i, j, n

  n = size(vec)

  ! Sort in ascending order
  vec_sorted = vec
  do i = 2, n
    x = vec_sorted(i)
    j = i - 1
    do while (j >= 1)
      if (vec_sorted(j) <= x) exit
      vec_sorted(j + 1) = vec_sorted(j)
      j = j - 1
    end do
    vec_sorted(j + 1) = x
  end do

  ! Calculate median according to if there is an even or odd number of elements
  if (mod(n, 2) == 0) then
    median = (vec_sorted(n/2) + vec_sorted(n/2+1))/2.0_jprd
  else
    median = vec_sorted((n+1)/2)
  endif

end function get_median

!===================================================================================================

subroutine initialize_spectral_arrays(nsmax, zsp, sp3d)

  integer,         intent(in)    :: nsmax       ! Spectral truncation
  real(kind=jprb), intent(inout) :: zsp(:,:)    ! Surface pressure
  real(kind=jprb), intent(inout) :: sp3d(:,:,:) ! 3D fields

  integer(kind=jpim) :: nflevl
  integer(kind=jpim) :: nfield

  integer :: i, j

  nflevl = size(sp3d, 1)
  nfield = size(sp3d, 3)

  ! First initialize surface pressure
  call initialize_2d_spectral_field(nsmax, zsp(1,:))

  ! Then initialize all of the 3D fields
  do i = 1, nflevl
    do j = 1, nfield
      call initialize_2d_spectral_field(nsmax, sp3d(i,:,j))
    end do
  end do

end subroutine initialize_spectral_arrays

!===================================================================================================

subroutine initialize_2d_spectral_field(nsmax, field)

  integer,         intent(in)    :: nsmax    ! Spectral truncation
  real(kind=jprb), intent(inout) :: field(:) ! Field to initialize

  integer :: index, num_my_zon_wns
  integer, allocatable :: my_zon_wns(:), nasm0(:)

  ! Choose a spherical harmonic to initialize arrays
  integer :: m_num = 4  ! Zonal wavenumber
  integer :: l_num = 19  ! Total wavenumber

  ! First initialise all spectral coefficients to zero
  field(:) = 0.0

  ! Get zonal wavenumbers this rank is responsible for
  call trans_inq(knump=num_my_zon_wns)
  allocate(my_zon_wns(num_my_zon_wns))
  call trans_inq(kmyms=my_zon_wns)

  ! If rank is responsible for the chosen zonal wavenumber...
  if (any(my_zon_wns == m_num) ) then
    ! Get array of spectral array addresses (this maps (m, n=m) to array index)
    allocate(nasm0(0:nsmax))
    call trans_inq(kasm0=nasm0)

    ! Find out local array index of chosen spherical harmonic
    index = nasm0(m_num) + 2 * (l_num - m_num) + 1

    ! Set just that element to a constant value
    field(index) = 1.0
  else
    return
  end if

end subroutine initialize_2d_spectral_field

!===================================================================================================

subroutine dump_gridpoint_field(jstep, myproc, nproma, gfld, fld, fldchar, noutdump)

  ! Dump a 2d field to a binary file.

  integer(kind=jpim), intent(in) :: jstep ! Time step, used for naming file
  integer(kind=jpim), intent(in) :: myproc ! MPI rank, used for naming file
  integer(kind=jpim), intent(in) :: nproma ! Size of nproma
  real(kind=jprb)   , intent(inout) :: gfld(:,:) ! 2d global field
  real(kind=jprb)   , intent(in) :: fld(:,:,:) ! 3d local field
  character         , intent(in) :: fldchar ! Single character field identifier
  integer(kind=jpim), intent(in) :: noutdump ! Tnit number for output file

  character(len=10) :: filename = "x.xxxx.dat"

  if (myproc == 1) then
    write(filename(1:1),'(a1)') fldchar
    write(filename(3:6),'(i4.4)') jstep
    open(noutdump,file=filename,form='unformatted')
  endif
  do ilev=1,size(fld,2)
    call gath_grid(gfld(:,:),nproma,1,(/1/),1,fld(:,ilev:ilev,:))
    if (myproc == 1) write(unit=noutdump) gfld(:,1)
  enddo
  if (myproc == 1) then
    close(noutdump)
  endif

end subroutine dump_gridpoint_field

!===================================================================================================

subroutine dump_crc(filename,iter, zgmv, zgmvs,sp3d,zspc2)
  real(kind=jprb), optional :: zgmv   (:,:,:,:) 
  real(kind=jprb), optional :: zgmvs   (:,:,:) 
  real(kind=jprb), optional :: sp3d  (:,:,:) 
  real(kind=jprb), optional :: zspc2   (:,:) 
  INTEGER(KIND = JPIM):: iter
  INTEGER*8 :: ICRC
  INTEGER(KIND = JPIM):: JLEV, JFLD
  character(len=*):: filename
  LOGICAL:: EXIST = .False.
  
  IF (iter>1)  INQUIRE(FILE = filename, EXIST = EXIST)
  IF (EXIST) THEN
    OPEN(10, FILE = filename, STATUS="OLD", POSITION="APPEND", ACTION="WRITE")
  ELSE
    OPEN(10, FILE = filename, ACTION="WRITE")
  END IF

  WRITE(10,*) "===================="
  WRITE(10,*) "iteration", iter  
  WRITE(10,*) "===================="

  IF (PRESENT(ZGMV)) THEN
    ICRC = 0
    DO JFLD = 1, SIZE (zgmv, 3)
      DO JLEV = 1, SIZE (zgmv, 2)
        CALL CRC64 (zgmv (:, JLEV, JFLD, :), INT (SIZE (zgmv (:, JLEV, JFLD, :)) * KIND (zgmv), 8), ICRC)
        WRITE (10, '(A," (",I0,", ",I0,") = ",Z16.16)') "zgmv", JLEV, JFLD, ICRC 
      ENDDO
    ENDDO
    ENDIF
  
  IF (PRESENT(ZGMVS)) THEN
  ICRC = 0
  DO JFLD = 1, SIZE (zgmvs, 2)
       CALL CRC64 (zgmvs (:, JFLD, :), INT (SIZE (zgmvs (:, JFLD, :)) * KIND (zgmvs), 8), ICRC)
       WRITE (10, '(A," (",I0,") = ",Z16.16)') "zgmvs", JFLD, ICRC            
  ENDDO
  ENDIF
  IF (PRESENT(sp3d)) THEN
  ICRC = 0
  DO JFLD = 1, SIZE (sp3d, 3)
    DO JLEV = 1, SIZE (sp3d, 1)
       CALL CRC64 (sp3d (JLEV,:, JFLD), INT (SIZE (sp3d (JLEV, :, JFLD)) * KIND (sp3d), 8), ICRC)
       WRITE (10, '(A," (",I0,", ",I0,") = ",Z16.16)') "sp3d", JLEV, JFLD, ICRC 
    ENDDO
  ENDDO
ENDIF
IF (PRESENT(zspc2)) THEN
  ICRC = 0
  DO JFLD = 1, SIZE (zspc2, 1)
       CALL CRC64 (zspc2 (JFLD, :), INT (SIZE (zspc2 (JFLD, :)) * KIND (zspc2), 8), ICRC)
       WRITE (10, '(A," (",I0,") = ",Z16.16)') "zspc2", JFLD, ICRC     
  ENDDO
ENDIF
  CLOSE(10)
end subroutine dump_crc

!===================================================================================================

function detect_mpirun() result(lmpi_required)
  use ec_env_mod, only : ec_putenv
  logical :: lmpi_required
  integer :: ilen
  integer, parameter :: nvars = 4
  character(len=32), dimension(nvars) :: cmpirun_detect
  character(len=4) :: clenv
  integer :: ivar

  ! Environment variables that are set when mpirun, srun, aprun, ... are used
  cmpirun_detect(1) = 'OMPI_COMM_WORLD_SIZE'  ! openmpi
  cmpirun_detect(2) = 'ALPS_APP_PE'           ! cray pe
  cmpirun_detect(3) = 'PMI_SIZE'              ! intel
  cmpirun_detect(4) = 'SLURM_NTASKS'          ! slurm

  lmpi_required = .false.
  do ivar = 1, nvars
    call get_environment_variable(name=trim(cmpirun_detect(ivar)), length=ilen)
    if (ilen > 0) then
      lmpi_required = .true.
      exit ! break
    endif
  enddo

  call get_environment_variable(name="ECTRANS_USE_MPI", value=clenv, length=ilen )
  if (ilen > 0) then
      lmpi_required = .true.
      if( trim(clenv) == "0" .or. trim(clenv) == "OFF" .or. trim(CLENV) == "off" .or. trim(clenv) == "F" ) then
        lmpi_required = .false.
      endif
      call ec_putenv("DR_HOOK_ASSERT_MPI_INITIALIZED=0", overwrite=.true.)
  endif
end function

!===================================================================================================

! Assign GSTATS labels to the main regions of ecTrans
subroutine gstats_labels

  call gstats_label(0,   '   ', 'PROGRAM        - Total')
  call gstats_label(1,   '   ', 'SETUP_TRANS0   - Setup ecTrans')
  call gstats_label(2,   '   ', 'SETUP_TRANS    - Setup ecTrans handle')
  call gstats_label(3,   '   ', 'TIME STEP      - Time step')
  call gstats_label(4,   '   ', 'INV_TRANS      - Inverse transform')
  call gstats_label(5,   '   ', 'DIR_TRANS      - Direct transform')
  call gstats_label(6,   '   ', 'NORMS          - Norm comp. (optional)')
  call gstats_label(102, '   ', 'LTINV_CTL      - Inv. Legendre transform')
  call gstats_label(103, '   ', 'LTDIR_CTL      - Dir. Legendre transform')
  call gstats_label(106, '   ', 'FTDIR_CTL      - Dir. Fourier transform')
  call gstats_label(107, '   ', 'FTINV_CTL      - Inv. Fourier transform')
  call gstats_label(140, '   ', 'SULEG          - Comp. of Leg. poly.')
  call gstats_label(152, '   ', 'LTINV_CTL      - M to L transposition')
  call gstats_label(153, '   ', 'LTDIR_CTL      - L to M transposition')
  call gstats_label(157, '   ', 'FTINV_CTL      - L to G transposition')
  call gstats_label(158, '   ', 'FTDIR_CTL      - G to L transposition')
  call gstats_label(400, '   ', 'GSTATS         - GSTATS itself')

end subroutine gstats_labels

!===================================================================================================

subroutine set_ectrans_gpu_nflev(kflev)
  use ec_env_mod, only : ec_putenv
  integer(kind=jpim), intent(in) :: kflev
  character(len=32) :: ECTRANS_GPU_NFLEV
  write(ECTRANS_GPU_NFLEV,'(A,I0)') "ECTRANS_GPU_NFLEV=",kflev
  call ec_putenv(ECTRANS_GPU_NFLEV, overwrite=.true.)
end subroutine

subroutine output_wrapped_file(wrap_fields)
  type(WRAPPED_FIELDS), INTENT(IN) :: wrap_fields
  write(nout,*) "wrap_fields%F_SPVOR", LOC(wrap_fields%F_SPVOR)
  write(nout,*) "wrap_fields%F_SPDIV", LOC(wrap_fields%F_SPDIV)
  write(nout,*) "wrap_fields%F_SPSCALARS", LOC(wrap_fields%F_SPSCALARS)
  write(nout,*) "wrap_fields%F_SPSCALARS2", LOC(wrap_fields%F_SPSCALARS2)

  write(nout,*) "wrap_fields%F_U", LOC(wrap_fields%F_U)
  write(nout,*) "wrap_fields%F_V", LOC(wrap_fields%F_V)
  write(nout,*) "wrap_fields%F_UDM", LOC(wrap_fields%F_UDM)
  write(nout,*) "wrap_fields%F_VDM", LOC(wrap_fields%F_VDM)
  write(nout,*) "wrap_fields%F_SCALARS", LOC(wrap_fields%F_SCALARS)
  write(nout,*) "wrap_fields%F_SCALARS_EW", LOC(wrap_fields%F_SCALARS_EW)
  write(nout,*) "wrap_fields%F_SCALARS_NS", LOC(wrap_fields%F_SCALARS_NS)
  write(nout,*) "wrap_fields%F_VOR", LOC(wrap_fields%F_VOR)
  write(nout,*) "wrap_fields%F_DIV", LOC(wrap_fields%F_DIV)

  write(nout,*) "wrap_fields%F_SCALARS2", LOC(wrap_fields%F_SCALARS2)
  write(nout,*) "wrap_fields%F_SCALARS2_EW", LOC(wrap_fields%F_SCALARS2_EW)
  write(nout,*) "wrap_fields%F_SCALARS2_NS", LOC(wrap_fields%F_SCALARS2_NS)
end subroutine output_wrapped_file

function allocate_wrapped_fields()
  type(WRAPPED_FIELDS), ALLOCATABLE :: allocate_wrapped_fields
  ALLOCATE(allocate_wrapped_fields)
  
  NULLIFY(allocate_wrapped_fields%F_SPVOR)
  NULLIFY(allocate_wrapped_fields%F_SPDIV)
  NULLIFY(allocate_wrapped_fields%F_SPSCALARS)
  NULLIFY(allocate_wrapped_fields%F_SPSCALARS2)

  NULLIFY(allocate_wrapped_fields%F_U)
  NULLIFY(allocate_wrapped_fields%F_V)
  NULLIFY(allocate_wrapped_fields%F_UDM)
  NULLIFY(allocate_wrapped_fields%F_VDM)
  NULLIFY(allocate_wrapped_fields%F_SCALARS)
  NULLIFY(allocate_wrapped_fields%F_SCALARS_EW)
  NULLIFY(allocate_wrapped_fields%F_SCALARS_NS)
  NULLIFY(allocate_wrapped_fields%F_VOR)
  NULLIFY(allocate_wrapped_fields%F_DIV)

  NULLIFY(allocate_wrapped_fields%F_SCALARS2)
  NULLIFY(allocate_wrapped_fields%F_SCALARS2_EW)
  NULLIFY(allocate_wrapped_fields%F_SCALARS2_NS)
end function allocate_wrapped_fields

function wrap_fields(lvordiv, lscders, luvders,& 
                    sp3d, spc2, zgmv, zgmvs, zgp2,& 
                    &jbegin_uv,jend_uv,& 
                    &jbegin_sc,jend_sc,& 
                    &jbegin_scder_NS, jend_scder_NS,& 
                    &jbegin_scder_EW, jend_scder_EW,& 
                    &jbegin_uder_EW, jend_uder_EW,& 
                    &jbegin_vder_EW, jend_vder_EW)

  type(WRAPPED_FIELDS), ALLOCATABLE :: wrap_fields
  logical :: lvordiv
  logical :: lscders
  logical :: luvders
  real(kind=jprb), INTENT(IN) :: sp3d(:,:,:)
  real(kind=jprb), INTENT(IN) :: spc2(:,:)
  real(kind=jprb), INTENT(IN) :: zgmv(:,:,:,:)
  real(kind=jprb), INTENT(IN) :: zgmvs(:,:,:)
  real(kind=jprb), INTENT(IN) :: zgp2 (:,:,:)
  integer(kind=jpim) :: jbegin_uv
  integer(kind=jpim) :: jend_uv
  integer(kind=jpim) :: jbegin_sc
  integer(kind=jpim) :: jend_sc
  integer(kind=jpim) :: jbegin_scder_NS
  integer(kind=jpim) :: jend_scder_NS
  integer(kind=jpim) :: jbegin_scder_EW
  integer(kind=jpim) :: jend_scder_EW
  integer(kind=jpim) :: jbegin_uder_EW
  integer(kind=jpim) :: jend_uder_EW
  integer(kind=jpim) :: jbegin_vder_EW
  integer(kind=jpim) :: jend_vder_EW
    
  wrap_fields = allocate_wrapped_fields()

  
  if (lvordiv) then
    IF (jbegin_uv>0 )      CALL FIELD_NEW(wrap_fields%F_U,         DATA=zgmv(:,:,jbegin_uv,:))
    IF (jend_uv>0 )        CALL FIELD_NEW(wrap_fields%F_V,         DATA=zgmv(:,:,jend_uv,:))
  !IF (jbegin_uv>0 )                    CALL FIELD_NEW(wrap_fields%F_VOR,         DATA=zgmv(:,:,jbegin_uv,:))
  !IF (jend_uv>0 )                      CALL FIELD_NEW(wrap_fields%F_DIV,         DATA=zgmv(:,:,jend_uv,:))

  endif
  
  IF (SIZE(sp3d,3) >=1 )   CALL FIELD_NEW(wrap_fields%F_SPVOR,        DATA=sp3d(:,:,1))
  IF (SIZE(sp3d,3) >=2 )   CALL FIELD_NEW(wrap_fields%F_SPDIV,        DATA=sp3d(:,:,2))
  IF (SIZE(sp3d,3) >=3 )   CALL FIELD_NEW(wrap_fields%F_SPSCALARS,    DATA=sp3d(:,:,3:))
  IF (SIZE(zspsc2,2) >=1 ) CALL FIELD_NEW(wrap_fields%F_SPSCALARS2,  DATA=zspsc2(:,:))
  
  IF (SIZE(ZGMVS,2)>=1)    CALL FIELD_NEW(wrap_fields%F_SCALARS2,    DATA=zgmvs(:,1:1,:))

  if (luvders) then
    CALL FIELD_NEW(wrap_fields%F_UDM,         DATA=zgmv(:,:,jbegin_uder_EW,:))
    CALL FIELD_NEW(wrap_fields%F_VDM,         DATA=zgmv(:,:,jbegin_vder_EW,:))
  endif

  IF (jend_sc>0 .AND. jend_sc>=jbegin_sc )  CALL FIELD_NEW(wrap_fields%F_SCALARS,     DATA=zgmv(:,:,jbegin_sc:jend_sc,:))
  if (lscders) then    
    IF (jend_scder_EW>0 .AND. jend_scder_EW>=jbegin_scder_EW ) CALL FIELD_NEW(wrap_fields%F_SCALARS_EW,  DATA=zgmv(:,:,jbegin_scder_EW:jend_scder_EW,:))
    IF (jend_scder_NS>0 .AND. jend_scder_NS>=jbegin_scder_NS ) CALL FIELD_NEW(wrap_fields%F_SCALARS_NS,  DATA=zgmv(:,:,jbegin_scder_NS:jend_scder_NS,:))
  
    IF (SIZE(ZGMVS,2)>=2)     CALL FIELD_NEW(wrap_fields%F_SCALARS2_NS, DATA=zgmvs(:,2:2,:))
    IF (SIZE(ZGMVS,2)>=3)     CALL FIELD_NEW(wrap_fields%F_SCALARS2_EW, DATA=zgmvs(:,3:3,:))      
  endif
end function wrap_fields

subroutine delete_wrapped_fields(wrap_fields)
  type(WRAPPED_FIELDS), INTENT(INOUT) :: wrap_fields
   
  IF(ASSOCIATED(wrap_fields%F_SPVOR)) CALL FIELD_DELETE(wrap_fields%F_SPVOR)
  IF(ASSOCIATED(wrap_fields%F_SPDIV)) CALL FIELD_DELETE(wrap_fields%F_SPDIV)
  IF(ASSOCIATED(wrap_fields%F_SPSCALARS)) CALL FIELD_DELETE(wrap_fields%F_SPSCALARS)
  IF(ASSOCIATED(wrap_fields%F_SPSCALARS2)) CALL FIELD_DELETE(wrap_fields%F_SPSCALARS2)

  IF(ASSOCIATED(wrap_fields%F_U)) CALL FIELD_DELETE(wrap_fields%F_U)
  IF(ASSOCIATED(wrap_fields%F_V)) CALL FIELD_DELETE(wrap_fields%F_V)
  IF(ASSOCIATED(wrap_fields%F_UDM)) CALL FIELD_DELETE(wrap_fields%F_UDM)
  IF(ASSOCIATED(wrap_fields%F_VDM)) CALL FIELD_DELETE(wrap_fields%F_VDM)
  IF(ASSOCIATED(wrap_fields%F_SCALARS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS)
  IF(ASSOCIATED(wrap_fields%F_SCALARS_EW)) CALL FIELD_DELETE(wrap_fields%F_SCALARS_EW)
  IF(ASSOCIATED(wrap_fields%F_SCALARS_NS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS_NS)
  IF(ASSOCIATED(wrap_fields%F_VOR)) CALL FIELD_DELETE(wrap_fields%F_VOR)
  IF(ASSOCIATED(wrap_fields%F_DIV)) CALL FIELD_DELETE(wrap_fields%F_DIV)

  IF(ASSOCIATED(wrap_fields%F_SCALARS2)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2)
  IF(ASSOCIATED(wrap_fields%F_SCALARS2_EW)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2_EW)
  IF(ASSOCIATED(wrap_fields%F_SCALARS2_NS)) CALL FIELD_DELETE(wrap_fields%F_SCALARS2_NS)
  
end subroutine delete_wrapped_fields

function allocate_fields_lists()
  type(FIELDS_LISTS), ALLOCATABLE :: allocate_fields_lists
  ALLOCATE(allocate_fields_lists)
  
  NULLIFY(allocate_fields_lists%U)
  NULLIFY(allocate_fields_lists%V)
  NULLIFY(allocate_fields_lists%SCALAR)
  NULLIFY(allocate_fields_lists%SPSCALAR)
  NULLIFY(allocate_fields_lists%SPVOR)
  NULLIFY(allocate_fields_lists%SPDIV)
  NULLIFY(allocate_fields_lists%VOR)
  NULLIFY(allocate_fields_lists%DIV)
  NULLIFY(allocate_fields_lists%UDM)
  NULLIFY(allocate_fields_lists%VDM)
  NULLIFY(allocate_fields_lists%SCALARDM)
  NULLIFY(allocate_fields_lists%SCALARDL)
end function allocate_fields_lists

subroutine delete_fields_lists(fl)
  type(FIELDS_LISTS), INTENT(INOUT) ::fl
  IF (ASSOCIATED(fl%U)) DEALLOCATE(fl%U)
  IF (ASSOCIATED(fl%V)) DEALLOCATE(fl%V)
  IF (ASSOCIATED(fl%SCALAR)) DEALLOCATE(fl%SCALAR)
  IF (ASSOCIATED(fl%SPSCALAR)) DEALLOCATE(fl%SPSCALAR)
  IF (ASSOCIATED(fl%SPVOR)) DEALLOCATE(fl%SPVOR)
  IF (ASSOCIATED(fl%SPDIV)) DEALLOCATE(fl%SPDIV)
  IF (ASSOCIATED(fl%VOR)) DEALLOCATE(fl%VOR)
  IF (ASSOCIATED(fl%DIV)) DEALLOCATE(fl%DIV)
  IF (ASSOCIATED(fl%UDM)) DEALLOCATE(fl%UDM)
  IF (ASSOCIATED(fl%VDM)) DEALLOCATE(fl%VDM) 
  IF (ASSOCIATED(fl%SCALARDM)) DEALLOCATE(fl%SCALARDM)
  IF (ASSOCIATED(fl%SCALARDL)) DEALLOCATE(fl%SCALARDL)
end subroutine delete_fields_lists

subroutine output_fields_lists(fl)
  type(FIELDS_LISTS), INTENT(IN) :: fl
  
  IF (ASSOCIATED(fl%U)) write(nout,*) "fl%U", SIZE(fl%U)  
  IF (ASSOCIATED(fl%V)) write(nout,*) "fl%V", SIZE(fl%V)
  IF (ASSOCIATED(fl%SCALAR)) write(nout,*) "fl%SCALAR", SIZE(fl%SCALAR)
  IF (ASSOCIATED(fl%SPSCALAR)) write(nout,*) "fl%SPSCALAR", SIZE(fl%SPSCALAR)
  IF (ASSOCIATED(fl%SPVOR)) write(nout,*) "fl%SPVOR", SIZE(fl%SPVOR)
  IF (ASSOCIATED(fl%SPDIV)) write(nout,*) "fl%SPDIV", SIZE(fl%SPDIV)
  IF (ASSOCIATED(fl%VOR)) write(nout,*) "fl%VOR", SIZE(fl%VOR)
  IF (ASSOCIATED(fl%DIV)) write(nout,*) "fl%DIV", SIZE(fl%DIV)
  IF (ASSOCIATED(fl%UDM)) write(nout,*) "fl%UDM", SIZE(fl%UDM)
  IF (ASSOCIATED(fl%VDM)) write(nout,*) "fl%VDM", SIZE(fl%VDM)
  IF (ASSOCIATED(fl%SCALARDM)) write(nout,*) "fl%SCALARDM", SIZE(fl%SCALARDM)
  IF (ASSOCIATED(fl%SCALARDL)) write(nout,*) "fl%SCALARDL", SIZE(fl%SCALARDL)
end subroutine output_fields_lists

function create_fields_lists(wf,NBSETLEV,NBSETSC)
  type(WRAPPED_FIELDS), INTENT(IN), ALLOCATABLE :: wf
  INTEGER(KIND=JPIM) :: NBSETLEV(:)
  INTEGER(KIND=JPIM) :: NBSETSC(:)
  
  
  type(FIELDS_LISTS), ALLOCATABLE ::create_fields_lists
  create_fields_lists = allocate_fields_lists()

  ASSOCIATE(ylf=>create_fields_lists,f=>wf) 

  IF(ASSOCIATED(f%F_SPVOR)) SPVOR=[B(f%F_SPVOR,'SP_VOR')]
  
  IF(ASSOCIATED(f%F_SPDIV)) SPDIV= [B(f%F_SPDIV,'SPDIV')]
  
  IF(ASSOCIATED(f%F_U)) U = [B(f%F_U,'U',NBSETLEV)]
  IF(ASSOCIATED(f%F_V)) V = [B(f%F_V,'V',NBSETLEV)]
      
  IF(ASSOCIATED(f%F_UDM)) UDM=[B(f%F_UDM,'U_FDM', NBSETLEV)]
  IF(ASSOCIATED(f%F_VDM)) VDM=[B(f%F_VDM,'V_FDM', NBSETLEV)]
    
  IF(ASSOCIATED(f%F_VOR))  VOR = [B(f%F_VOR,'VOR', NBSETLEV)]
  IF(ASSOCIATED(f%F_DIV))  DIV = [B(f%F_DIV,'DIV', NBSETLEV)]
  
  IF (ASSOCIATED(f%F_SPSCALARS) .AND. ASSOCIATED(f%F_SPSCALARS2) ) THEN
      SPSCALAR = [B(f%F_SPSCALARS,'SP_SCALARS'), B(f%F_SPSCALARS2,'SP_SCALARS2')]
  ELSE IF (ASSOCIATED(f%F_SPSCALARS)) THEN
      SPSCALAR = [B(f%F_SPSCALARS,'SP_SCALARS')]    
  ELSE IF (ASSOCIATED(f%F_SPSCALARS2)) THEN
    SPSCALAR = [B(f%F_SPSCALARS2,'SP_SCALARS2')]  
  ENDIF
    
  IF (ASSOCIATED(f%F_SCALARS) .AND. ASSOCIATED(f%F_SCALARS2) ) THEN
    SCALAR = [B(f%F_SCALARS,'SCALARS', NBSETLEV), B(f%F_SCALARS2,'SCALARS2', NBSETSC)]
  ELSE IF (ASSOCIATED(f%F_SCALARS)) THEN
    SCALAR = [B(f%F_SCALARS,'SCALARS', NBSETLEV)]    
  ELSE IF (ASSOCIATED(f%F_SCALARS2)) THEN
    SCALAR = [B(f%F_SCALARS2,'SCALARS2', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(f%F_SCALARS_NS) .AND. ASSOCIATED(f%F_SCALARS2_NS) ) THEN
    SCALARDM = [B(f%F_SCALARS_NS,'SCALARS_NS', NBSETLEV), B(f%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]
  ELSE IF (ASSOCIATED(f%F_SCALARS_NS)) THEN
    SCALARDM = [B(f%F_SCALARS_NS,'SCALARS_NS', NBSETLEV)]
  ELSE IF (ASSOCIATED(f%F_SCALARS2_NS)) THEN
    SCALARDM = [B(f%F_SCALARS2_NS,'SCALARS2_NS', NBSETSC)]  
  ENDIF
  
  IF (ASSOCIATED(f%F_SCALARS_EW) .AND. ASSOCIATED(f%F_SCALARS2_EW) ) THEN
    SCALARDL = [B(f%F_SCALARS_EW,'SCALARS_EW', NBSETLEV), B(f%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]
  ELSE IF (ASSOCIATED(f%F_SCALARS_EW)) THEN
    SCALARDL = [B(f%F_SCALARS_EW,'SCALARS_EW', NBSETLEV)]    
  ELSE IF (ASSOCIATED(f%F_SCALARS2_EW)) THEN
    SCALARDL = [B(f%F_SCALARS2_EW,'SCALARS2_EW', NBSETSC)]  
  ENDIF

  IF (ALLOCATED(U)) ylf%U=>U
  IF (ALLOCATED(V)) ylf%V=>V
  IF (ALLOCATED(SCALAR)) ylf%SCALAR=>SCALAR
  IF (ALLOCATED(SPSCALAR)) ylf%SPSCALAR=>SPSCALAR
  IF (ALLOCATED(SPVOR)) ylf%SPVOR=>SPVOR
  IF (ALLOCATED(SPDIV)) ylf%SPDIV=>SPDIV
  IF (ALLOCATED(VOR)) ylf%VOR=>VOR
  IF (ALLOCATED(DIV)) ylf%DIV=>DIV
  IF (ALLOCATED(UDM)) ylf%UDM=>UDM
  IF (ALLOCATED(VDM)) ylf%VDM=>VDM 
  IF (ALLOCATED(SCALARDM)) ylf%SCALARDM=>SCALARDM
  IF (ALLOCATED(SCALARDL)) ylf%SCALARDL=>SCALARDL

END ASSOCIATE

end function create_fields_lists

end program ectrans_benchmark

!===================================================================================================

