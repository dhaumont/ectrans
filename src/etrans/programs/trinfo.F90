PROGRAM TRINFO

USE PARKIND1, ONLY : JPIM, JPRB
USE XRD_GETOPTIONS
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT
USE SUEMPLAT_MOD
USE SUMPLAT_MOD
USE EQ_REGIONS_MOD

IMPLICIT NONE

CHARACTER (LEN=128) :: CLFA
INTEGER (KIND=JPIM) :: INBARP, INBARI, IREP
CHARACTER (LEN=*), PARAMETER :: CLNOMC = 'cadre'
INTEGER (KIND=JPIM), PARAMETER :: ILUN = 77

INTEGER (KIND=JPIM) :: ITYPTR, ITRONC, INLATI, INXLON, INIVER
INTEGER (KIND=JPIM), ALLOCATABLE :: INLOPA (:), INOZPA (:)
REAL (KIND=JPRB) :: ZSLAPO, ZCLOPO, ZSLOPO, ZCODIL, ZREFER
REAL (KIND=JPRB), ALLOCATABLE :: ZSINLA (:), ZAHYBR (:), ZBHYBR (:)
LOGICAL :: LLGARD

LOGICAL :: LELAM
INTEGER (KIND=JPIM) :: NSMAX, NMSMAX, NFLEVG, NDGLG, NDLON
INTEGER (KIND=JPIM), ALLOCATABLE :: NLOENG (:)

INTEGER (KIND=JPIM) :: NPRINTLEV
INTEGER (KIND=JPIM) :: NPROC, NPRGPNS, NPRGPEW, NDGUXG
LOGICAL :: LEQ_REGIONS, LSPLIT, LLWEIGHTED_DISTR

INTEGER (KIND=JPIM), ALLOCATABLE :: NFRSTLAT (:), NLSTLAT (:), NPROCAGP (:), &
                                  & NPTRLAT (:), NPTRFRSTLAT (:), NPTRLSTLAT (:)
LOGICAL, ALLOCATABLE :: LSPLITLAT (:)

INTEGER (KIND=JPIM) :: IFRSTLOFF, IPTRFLOFF, IMEDIAP, IRESTM
REAL (KIND=JPRB) :: ZWEIGHT (1), ZMEDIAP

INTEGER (KIND=JPIM) :: JGL

TYPE (EQ_REGIONS_T) :: YLER

CALL INITOPTIONS (KOPTMIN=0)

CALL GETOPTION ('--fa-file', CLFA, MND = .TRUE., USE = "FA file")

NPRINTLEV = 1
CALL GETOPTION ('--nprintlev', NPRINTLEV)

NPRGPNS = 1
CALL GETOPTION ('--nprgpns', NPRGPNS)

NPRGPEW = 1
CALL GETOPTION ('--nprgpew', NPRGPEW)

NPROC = NPRGPNS * NPRGPEW
CALL GETOPTION ('--nproc', NPROC)

CALL GETOPTION ('--leq_regions', LEQ_REGIONS)

CALL GETOPTION ('--lsplit', LSPLIT)

CALL CHECKOPTIONS ()

CALL FAITOU (IREP, ILUN, .TRUE., CLFA, 'OLD', .TRUE., .FALSE., 0_JPIM, INBARP, INBARI, CLNOMC)

ALLOCATE (INLOPA (FA%JPXPAH), INOZPA (FA%JPXIND), &
        & ZSINLA (FA%JPXGEO), ZAHYBR (0:FA%JPXNIV), ZBHYBR (0:FA%JPXNIV))

CALL FACIES (CLNOMC, ITYPTR, ZSLAPO, ZCLOPO, ZSLOPO, &
&            ZCODIL, ITRONC, INLATI, INXLON, INLOPA, &
&            INOZPA, ZSINLA, INIVER, ZREFER, ZAHYBR, &
&            ZBHYBR, LLGARD)


CALL FAIRME (IREP, ILUN, 'KEEP')

LELAM = ITYPTR < 0
NSMAX = ITRONC
IF (LELAM) THEN
  NMSMAX = - ITYPTR
  NDGUXG = INLOPA (6)
ELSE
  NMSMAX = NSMAX
  NDGUXG = 0
ENDIF

NFLEVG = INIVER
NDGLG  = INLATI
NDLON  = INXLON

WRITE (*, *) " NPROC       = ", NPROC
WRITE (*, *) " NPRGPNS     = ", NPRGPNS
WRITE (*, *) " NPRGPEW     = ", NPRGPEW
WRITE (*, *) " LELAM       = ", LELAM
WRITE (*, *) " NDLON       = ", NDLON
WRITE (*, *) " NDGLG       = ", NDGLG
WRITE (*, *) " NDGUXG      = ", NDGUXG
WRITE (*, *) " NFLEVG      = ", NFLEVG
WRITE (*, *) " NSMAX       = ", NSMAX
WRITE (*, *) " NMSMAX      = ", NMSMAX
WRITE (*, *) " LSPLIT      = ", LSPLIT
WRITE (*, *) " LEQ_REGIONS = ", LEQ_REGIONS

ALLOCATE (NLOENG (NDGLG))

IF (LELAM) THEN
  NLOENG = NDLON
ELSE
  DO JGL = 1, (NDGLG+1)/2
    NLOENG (JGL) = INLOPA (JGL)
  ENDDO
  DO JGL = (NDGLG+1)/2+1, NDGLG
    NLOENG (JGL) = INLOPA (NDGLG-JGL+1)
  ENDDO
ENDIF

WRITE (*, *) " NLOENG = ", NLOENG

CALL EQ_REGIONS_SAVE (YLER)

IF (LEQ_REGIONS) THEN
  ALLOCATE (N_REGIONS (NPROC+2))
  N_REGIONS = 0
  CALL EQ_REGIONS (NPROC)
ELSE
  N_REGIONS_NS = NPRGPNS
  ALLOCATE (N_REGIONS (N_REGIONS_NS))
  N_REGIONS    = NPRGPEW
  N_REGIONS_EW = NPRGPEW
ENDIF

WRITE (*, *) " N_REGIONS    = ", N_REGIONS
WRITE (*, *) " N_REGIONS_NS = ", N_REGIONS_NS
WRITE (*, *) " N_REGIONS_EW = ", N_REGIONS_EW

ALLOCATE (NFRSTLAT (N_REGIONS_NS), NLSTLAT (N_REGIONS_NS), NPROCAGP (N_REGIONS_NS), &
        & NPTRLAT (NDGLG), NPTRFRSTLAT (N_REGIONS_NS), NPTRLSTLAT (N_REGIONS_NS),   &
        & LSPLITLAT (NDGLG))


LLWEIGHTED_DISTR = .FALSE.

IF (LELAM) THEN

  CALL SUEMPLAT (KDGL               = NDGLG              , &
               & KPROC              = NPROC              , &
               & KPROCA             = N_REGIONS_NS       , &
               & KMYSETA            = 1_JPIM             , &
               & LDSPLIT            = LSPLIT             , &
               & LDEQ_REGIONS       = LEQ_REGIONS        , &
               & KFRSTLAT           = NFRSTLAT           , &
               & KLSTLAT            = NLSTLAT            , &
               & KFRSTLOFF          = IFRSTLOFF          , &
               & KPTRLAT            = NPTRLAT            , &
               & KPTRFRSTLAT        = NPTRFRSTLAT        , &
               & KPTRLSTLAT         = NPTRLSTLAT         , &
               & KPTRFLOFF          = IPTRFLOFF          , &
               & PWEIGHT            = ZWEIGHT            , &
               & LDWEIGHTED_DISTR   = LLWEIGHTED_DISTR   , &
               & PMEDIAP            = ZMEDIAP            , &
               & KPROCAGP           = NPROCAGP           , &
               & KMEDIAP            = IMEDIAP            , &
               & KRESTM             = IRESTM             , &
               & LDSPLITLAT         = LSPLITLAT          , &
               & KMYPROC            = 1_JPIM             , &
               & KLOEN              = NLOENG             , &
               & KDGUX              = NDGLG) ! Should be NDGUXG+INT((RDISTR_E*(NDGLG-NDGUXG)))

ELSE
  
  CALL SUMPLAT (KDGL               = NDGLG              , &
              & KPROC              = NPROC              , &
              & KPROCA             = N_REGIONS_NS       , &
              & KMYSETA            = 1_JPIM             , &
              & LDSPLIT            = LSPLIT             , &
              & LDEQ_REGIONS       = LEQ_REGIONS        , &
              & KFRSTLAT           = NFRSTLAT           , &
              & KLSTLAT            = NLSTLAT            , &
              & KFRSTLOFF          = IFRSTLOFF          , &
              & KPTRLAT            = NPTRLAT            , &
              & KPTRFRSTLAT        = NPTRFRSTLAT        , &
              & KPTRLSTLAT         = NPTRLSTLAT         , &
              & KPTRFLOFF          = IPTRFLOFF          , &
              & PWEIGHT            = ZWEIGHT            , &
              & LDWEIGHTED_DISTR   = LLWEIGHTED_DISTR   , &
              & PMEDIAP            = ZMEDIAP            , &
              & KPROCAGP           = NPROCAGP           , &
              & KMEDIAP            = IMEDIAP            , &
              & KRESTM             = IRESTM             , &
              & LDSPLITLAT         = LSPLITLAT          , &
              & KMYPROC            = 1_JPIM             , &
              & KLOEN              = NLOENG             )

ENDIF

WRITE (*, *) " NFRSTLAT    = ", NFRSTLAT
WRITE (*, *) " NLSTLAT     = ", NLSTLAT
WRITE (*, *) " NPTRLAT     = ", NPTRLAT
WRITE (*, *) " NPTRFRSTLAT = ", NPTRFRSTLAT
WRITE (*, *) " NPTRLSTLAT  = ", NPTRLSTLAT

IF (LEQ_REGIONS) THEN
  DEALLOCATE (N_REGIONS)
  NULLIFY (N_REGIONS)
  CALL EQ_REGIONS_LOAD (YLER)
ENDIF

END PROGRAM TRINFO

