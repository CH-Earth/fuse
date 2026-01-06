MODULE multiforce

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Nans Addor to enable distributed modeling, 9/2016
 ! Modified by Cyril Thébault to allow different metrics as objective function, 2024
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype

 USE multiforce_types, only: TDATA, VDATA, ADATA, FDATA

 implicit none
 private

 public :: forcefile

 public :: ncid_forc, ncid_var

 public :: nForce, nInput

 public :: timDat, valDat, aValid
 public :: AFORCE, CFORCE, MFORCE
 public :: ancilF, ancilF_3d
 public :: gForce, gForce_3d
 
 public :: date_start_input, date_end_input
 public :: numtim_in, numtim_sim, numtim_sub, numtim_sub_cur
 public :: itim_in, itim_sim, itim_sub
 public :: sim_beg, sim_end, eval_beg, eval_end
 public :: istart, jdayRef
 public :: deltim

 public :: SUB_PERIODS_FLAG, GRID_FLAG

 public :: startSpat2, nSpat1, nSpat2
 public :: xlon, ylat, latitude, longitude
 public :: latUnits, lonUnits, timeUnits

 public :: time_steps, julian_day_input

 public :: NUMPSET, name_psets

 public :: vname_iy, vname_im, vname_id, vname_ih, vname_imin, vname_dsec, vname_dtime

 public :: vname_aprecip, vname_potevap, vname_airtemp, vname_q, vname_spechum, vname_airpres, vname_swdown
 public :: ilook_aprecip, ilook_potevap, ilook_airtemp, ilook_q, ilook_spechum, ilook_airpres, ilook_swdown

 public :: ivarid_iy, ivarid_im, ivarid_id, ivarid_ih, ivarid_imin, ivarid_dsec
 public :: ivarid_ppt, ivarid_temp, ivarid_pet, ivarid_q

 public :: amult_ppt, amult_pet, amult_q

 public ::  NA_VALUE, NA_VALUE_SP

 SAVE
 
 ! general
 INTEGER(I4B),PARAMETER                  :: STRLEN=256 ! length of the character string
 
 ! time data structures
 TYPE(tData)                             :: timDat     ! model time structure
 
 ! response data structures
 TYPE(vData)                             :: valDat     ! validation structure
 TYPE(vData), DIMENSION(:,:,:), POINTER  :: aValid     ! all model validation data
 
 ! forcing data structures
 TYPE(FDATA), DIMENSION(:), POINTER      :: AFORCE     ! all model forcing data
 TYPE(FDATA), DIMENSION(:), POINTER      :: CFORCE     ! COPY of model forcing data
 TYPE(FDATA)                             :: MFORCE     ! model forcing data for a single time step
 TYPE(aData), DIMENSION(:,:), POINTER    :: ancilF     ! ancillary forcing data for the 2-d grid
 TYPE(fData), DIMENSION(:,:), POINTER    :: gForce     ! model forcing data for a 2-d grid
 TYPE(fData), DIMENSION(:,:,:), POINTER  :: gForce_3d  ! model forcing data for a 3-d grid (time as 3rd dimension)
 TYPE(aData), DIMENSION(:,:,:), POINTER  :: ancilF_3d  ! ancillary forcing data for the 3-d grid

 ! NetCDF

 CHARACTER(len=StrLen)                 :: forcefile    = 'undefined'   ! name of forcing file

 INTEGER(i4b), PARAMETER               :: nForce = 7                   ! number of forcing variables
 INTEGER(i4b)                          :: nInput = 3                   ! number of variable to retrieve from input file

 INTEGER(i4b)                          :: ncid_forc = -1               ! NetCDF forcing file ID
 INTEGER(i4b), DIMENSION(nForce)       :: ncid_var                     ! NetCDF forcing variable ID

 ! timing information - note that numtim_in >= numtim_sim >= numtim_sub
 
 CHARACTER(len=20)                     :: date_start_input            ! date start input time series
 CHARACTER(len=20)                     :: date_end_input              ! date end input time series

 INTEGER(i4b)                          :: numtim_in = -1              ! number of time steps of input (atmospheric forcing)
 INTEGER(i4b)                          :: numtim_sim = -1             ! number of time steps of FUSE simulations (including spin-up)
 INTEGER(i4b)                          :: numtim_sub = -1             ! number of time steps of subperiod (will be kept in memory)
 INTEGER(i4b)                          :: numtim_sub_cur = -1         ! number of time steps of current subperiod (allows for the last subperiod to be shorter)
 INTEGER(i4b)                          :: itim_in = -1                ! indice within numtim_in
 INTEGER(i4b)                          :: itim_sim = -1               ! indice within numtim_sim
 INTEGER(i4b)                          :: itim_sub = -1               ! indice within numtim_sub

 INTEGER(i4b)                          :: sim_beg = -1                ! index for the start of the simulation in fuse_metric
 INTEGER(i4b)                          :: sim_end = -1                ! index for the end of the simulation in fuse_metric
 INTEGER(i4b)                          :: eval_beg = -1               ! index for the start of evaluation period
 INTEGER(i4b)                          :: eval_end = -1               ! index for the end of the inference period

 INTEGER(i4b)                          :: istart = -1                 ! index for start of inference period (in reduced array)
 REAL(sp)                              :: jdayRef                     ! reference time (days)
 REAL(sp)                              :: deltim = -1._dp             ! length of time step (days)

 LOGICAL(LGT)                          :: SUB_PERIODS_FLAG            ! .true. if subperiods are used to run FUSE
 LOGICAL(LGT)                          :: GRID_FLAG                   ! spatial flag .true. if grid

 ! dimension information
 
 INTEGER(i4b)                          :: startSpat2 = -1             ! number of points in 1st spatial dimension
 INTEGER(i4b)                          :: nSpat1 = -1                 ! number of points in 1st spatial dimension
 INTEGER(i4b)                          :: nSpat2 = -1                 ! number of points in 2nd spatial dimension
 REAL(sp)                              :: xlon                        ! longitude (degrees) for PET computation
 REAL(sp)                              :: ylat                        ! latitude (degrees) for PET computation
 REAL(sp),dimension(:),allocatable     :: latitude                    ! latitude (degrees)
 REAL(sp),dimension(:),allocatable     :: longitude                   ! longitude (degrees)
 CHARACTER(len=strLen)                 :: latUnits                    ! units string for latitude
 CHARACTER(len=strLen)                 :: lonUnits                    ! units string for longitude
 CHARACTER(len=strLen)                 :: timeUnits                   ! units string for time

 REAL(sp),dimension(:),allocatable     :: time_steps                  ! time steps (days)
 REAL(sp),dimension(:),allocatable     :: julian_day_input            ! time steps (julian days)

 INTEGER(I4B)                          :: NUMPSET                     ! number of parameter sets
 CHARACTER(len=strLen),dimension(:),allocatable   :: name_psets       ! name of parameter sets

 ! name of time variables
 CHARACTER(len=StrLen)                 :: vname_iy     = 'undefined'   ! name of variable for year
 CHARACTER(len=StrLen)                 :: vname_im     = 'undefined'   ! name of variable for month
 CHARACTER(len=StrLen)                 :: vname_id     = 'undefined'   ! name of variable for day
 CHARACTER(len=StrLen)                 :: vname_ih     = 'undefined'   ! name of variable for hour
 CHARACTER(len=StrLen)                 :: vname_imin   = 'undefined'   ! name of variable for minute
 CHARACTER(len=StrLen)                 :: vname_dsec   = 'undefined'   ! name of variable for second
 CHARACTER(len=StrLen)                 :: vname_dtime  = 'undefined'   ! name of variable for time

 ! forcing variable names
 CHARACTER(len=StrLen)                 :: vname_aprecip = 'undefined'  ! variable name: precipitation
 CHARACTER(len=StrLen)                 :: vname_potevap = 'undefined'  ! variable name: potential ET
 CHARACTER(len=StrLen)                 :: vname_airtemp = 'undefined'  ! variable name: temperature
 CHARACTER(len=StrLen)                 :: vname_q       = 'undefined'  ! variable name: observed runoff
 CHARACTER(len=StrLen)                 :: vname_spechum = 'undefined'  ! variable name: specific humidity
 CHARACTER(len=StrLen)                 :: vname_airpres = 'undefined'  ! variable name: surface pressure
 CHARACTER(len=StrLen)                 :: vname_swdown  = 'undefined'  ! variable name: downward shortwave radiation

 ! indices for forcing variables
 INTEGER(i4b),PARAMETER                :: ilook_aprecip = 1  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_potevap = 2  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_airtemp = 3  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_q       = 4  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_spechum = 5  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_airpres = 6  ! named element in lCheck
 INTEGER(i4b),PARAMETER                :: ilook_swdown  = 7  ! named element in lCheck

 ! indices for time data (only used in ASCII files)
 INTEGER(i4b)                          :: ivarid_iy   = -1   ! variable ID for year
 INTEGER(i4b)                          :: ivarid_im   = -1   ! variable ID for month
 INTEGER(i4b)                          :: ivarid_id   = -1   ! variable ID for day
 INTEGER(i4b)                          :: ivarid_ih   = -1   ! variable ID for hour
 INTEGER(i4b)                          :: ivarid_imin = -1   ! variable ID for minute
 INTEGER(i4b)                          :: ivarid_dsec = -1   ! variable ID for second

 ! indices for variables
 INTEGER(i4b)                          :: ivarid_ppt  = -1   ! variable ID for precipitation
 INTEGER(i4b)                          :: ivarid_temp = -1   ! variable ID for temperature
 INTEGER(i4b)                          :: ivarid_pet  = -1   ! variable ID for potential ET
 INTEGER(i4b)                          :: ivarid_q    = -1   ! variable ID for runoff

  ! multipliers for variables to convert fluxes to mm/day
 REAL(sp)                              :: amult_ppt = -1._dp ! convert precipitation to mm/day
 REAL(sp)                              :: amult_pet = -1._dp ! convert potential ET to mm/day
 REAL(sp)                              :: amult_q   = -1._dp ! convert runoff to mm/day

 ! missing values
 INTEGER(I4B),PARAMETER                :: NA_VALUE    = -9999     ! integer designating missing values - TODO: retrieve from NetCDF file
 REAL(SP),PARAMETER                    :: NA_VALUE_SP = -9999._sp ! integer designating missing values - TODO: retrieve from NetCDF file

END MODULE multiforce
