module info_types

 use nrtype

 use multiparam_types, only: par_id

 private
 public :: cli_options
 public :: time_info
 public :: space_info
 public :: fuse_info

 ! --------------------------------------------------------------------------------------
 
 type :: mpi_info
   logical(lgt) :: enabled = .false.
   integer(i4b) :: rank    = 0
   integer(i4b) :: nproc   = 1
 end type mpi_info

 ! -------------------------------------------------------------------------------------

 ! options for the command-line interface

 type :: cli_options
   character(len=:), allocatable :: tag              ! string to add to output file
   character(len=:), allocatable :: control_file
   character(len=:), allocatable :: domain_id
   character(len=:), allocatable :: runmode          ! def/idx/opt/sce
   character(len=:), allocatable :: sets_file        ! for idx,opt
   integer(i4b)                  :: indx = -1        ! for idx
   character(len=:), allocatable :: restart_freq     ! y/m/d/e/never
   logical(lgt)                  :: show_version = .false.
   logical(lgt)                  :: show_help    = .false.
   character(len=:), allocatable :: param_name(:)    ! list of parameter names
   real(wp), allocatable         :: param_value(:)   ! list of parameter values
 end type cli_options

 ! -------------------------------------------------------------------------------------

 type :: space_info

   ! observation dimension
   integer(i4b) :: nobs = 1

   ! HRU and stream segment dimensions
   integer(i4b) :: n_hru = 1
   integer(i4b) :: n_seg = 1

   ! global dimensions (full forcing file)
   integer(i4b) :: nx_global = 1
   integer(i4b) :: ny_global = 1

   ! local dimensions (after MPI split)
   integer(i4b) :: nx_local = 1
   integer(i4b) :: ny_local = 1

   ! decomposition along y dimension
   integer(i4b) :: y_start_global = 1
   integer(i4b) :: y_end_global   = 1

   ! coordinate geometry
   logical(lgt) :: is_gridded = .false.
   logical(lgt) :: is_clinear = .false.
 
 end type space_info

 ! -------------------------------------------------------------------------------------
 
 type :: time_info
 
   ! forcing axis (global)
   integer(i4b) :: nt_global = 0
  
   ! simulation & evaluation indices into forcing time axis
   integer(i4b) :: sim_beg  = 1
   integer(i4b) :: sim_end  = 1
   integer(i4b) :: eval_beg = 1
   integer(i4b) :: eval_end = 1
  
   ! derived lengths
   integer(i4b) :: nt_sim   = 0
  
   ! subperiod / windowing
   logical(lgt) :: use_subperiods = .false.
   integer(i4b) :: nt_window      = 0         ! (= numtim_sub)
   integer(i4b) :: nt_window_cur  = 0         ! runtime: current window length
  
   ! bookkeeping for time axis
   character(len=:), allocatable :: units
   real(wp)       :: jdate_ref = 0._wp
   real(wp), allocatable :: jdate(:)          ! julian day for each forcing record

   real(wp), allocatable :: time_steps(:)     ! time since reference time (transferred to output)
   real(wp), allocatable :: time_bounds(:,:)  ! (2, nt)

   real(wp)              :: deltim_days       ! forcing time step in units of days

 end type time_info

 ! -------------------------------------------------------------------------------------

 type :: snow_info
   integer(i4b) :: n_bands = 0
 end type snow_info

 ! -------------------------------------------------------------------------------------
 
 ! --- hydromet_vars used in file_info
 
 type :: hydromet_vars ! NVAR_HYDROMET
  character(len=64), allocatable :: name(:) 
  character(len=64), allocatable :: units(:)
  integer(i4b),      allocatable :: varid(:)
  integer(i4b),      allocatable :: ndims(:)
  real(wp),          allocatable :: multiplier(:)
  real(wp),          allocatable :: fill_value(:)
  logical(lgt),      allocatable :: has_fill_value(:)
 end type

 ! ---

 type :: file_info

   ! directories
   character(len=:), allocatable :: setngs_path
   character(len=:), allocatable :: input_path
   character(len=:), allocatable :: output_path
  
   ! settings filenames
   character(len=:), allocatable :: constraints
   character(len=:), allocatable :: mod_numerix
   character(len=:), allocatable :: m_decisions
  
   ! domain-derived input suffixes
   character(len=:), allocatable :: suffix_hydromet
   character(len=:), allocatable :: suffix_elev_bands
  
   ! actual input filenames for this domain (derived once dom_id known)
   character(len=512) :: hydromet_file   ! dom_id//suffix_hydromet
   character(len=512) :: elevbands_file  ! dom_id//suffix_elev_bands
  
   ! output base name + concrete outputs
   character(len=512) :: fname_tempry
   character(len=512) :: fname_netcdf_hmet
   character(len=512) :: fname_netcdf_runs
   character(len=512) :: fname_netcdf_para

   ! NetCDF hydromet file info
   integer(i4b)                  :: ncid_hydromet = -9999  ! NetCDF file ID for hydromet data

   character(len=:), allocatable :: time_name          ! name of coordinate variables
   character(len=:), allocatable :: latitude_name      ! name of coordinate variables
   character(len=:), allocatable :: longitude_name     ! name of coordinate variables
   
   character(len=:), allocatable :: precip_name        ! name of hydromet variables   
   character(len=:), allocatable :: temp_name          ! name of hydromet variables   
   character(len=:), allocatable :: pet_name           ! name of hydromet variables  
   character(len=:), allocatable :: qobs_name          ! name of hydromet variables

   type(hydromet_vars)           :: hydromet           ! name/varid table

 end type file_info

 
 ! -------------------------------------------------------------------------------------

 type :: mizu_info

   character(len=:), allocatable :: namelist_path       ! namelist path
   character(len=:), allocatable :: namelist_file       ! namelist file

   character(len=:), allocatable :: methods             ! string of integers defining methods

   real(wp)                      :: dt = 3600._wp       ! routing time step (s)

 end type mizu_info

 ! -------------------------------------------------------------------------------------

 type :: topo_info

   ! Hydrofabric path/filenames
   character(len=:), allocatable :: hfabric_path       ! hydrofabric path
   character(len=:), allocatable :: hfabric_file       ! hydrofabric file
   character(len=:), allocatable :: hfabric_newfile    ! hydrofabric file (new)
  
   ! NetCDF dimensions
   character(len=:), allocatable :: dname_hru          ! dimension name: hru ID
   character(len=:), allocatable :: dname_seg          ! dimension name: segment
  
   ! NetCDF variable names
   character(len=:), allocatable :: varname_HRUid      ! variable name: HRU ID
   character(len=:), allocatable :: varname_segId      ! variable name: segment ID
   character(len=:), allocatable :: varname_hruSegId   ! variable name: ID of segment in HRU
   character(len=:), allocatable :: varname_downSegId  ! variable name: downstream segment ID
  
   character(len=:), allocatable :: varname_area       ! variable name: HRU area
   character(len=:), allocatable :: varname_slope      ! variable name: segment slope
   character(len=:), allocatable :: varname_length     ! variable name: segment length

   ! Network topology
   integer(i4b)                  :: idSegOut = -9999   ! ID of outlet segment
   integer(i4b)                  :: ixSegOut = -9999   ! index of outlet segment

 end type topo_info

 ! -------------------------------------------------------------------------------------

 type :: remap_info

   ! Remapping filename
   character(len=:), allocatable :: remap_file         ! remapping file

   ! NetCDF dimensions
   character(len=:), allocatable :: dname_hru          ! name of dimension of river network HRU ID
   character(len=:), allocatable :: dname_data         ! name of dimension of runoff HRU overlapping with river network HRU

   ! NetCDF variable names
   character(len=:), allocatable :: vname_hruid        ! name of variable containing ID of river network HRU
   character(len=:), allocatable :: vname_weight       ! name of variable contating areal weights of runoff HRUs within each river network HRU
   character(len=:), allocatable :: vname_num_qhru     ! name of variable containing numbers of runoff HRUs within each river network HRU
   character(len=:), allocatable :: vname_i_index      ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
   character(len=:), allocatable :: vname_j_index      ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)

 end type remap_info

 ! -------------------------------------------------------------------------------------


 type :: run_config

  ! provenance
  character(len=512) :: file_manager_file = ""

  ! CLI options
  type(cli_options)  :: cli_opts

  ! model selection
  character(len=:), allocatable :: fmodel_id

  ! model information
  integer(i4b)       :: nState = -9999
  integer(i4b)       :: nParam = -9999
 
  ! number of output variables
  integer(i4b)       :: nOutput

  ! list of output variables
  character(len=strLen), allocatable :: outvar_names(:)

  ! list of model parameters
  type(par_id), allocatable  :: listParam(:)

  ! run flags
  logical(lgt) :: write_timeseries = .true.

  ! requested time windows (strings as read from filemanager)
  character(len=:), allocatable :: date_start_sim
  character(len=:), allocatable :: date_end_sim
  character(len=:), allocatable :: date_start_eval
  character(len=:), allocatable :: date_end_eval
  character(len=:), allocatable :: numtim_sub_str

  ! parsed / derived values (optional convenience)
  integer(i4b) :: numtim_sub = -9999      ! parsed from numtim_sub_str

  ! output dimension for number of parameter sets
  integer(i4b) :: nSets

  ! calibration metrics and metric transformations
  character(len=:), allocatable :: metric
  character(len=:), allocatable :: transfo

  ! SCE settings (store as numeric types)
  integer(i4b)      :: maxn  = -9999
  integer(i4b)      :: kstop = -9999
  real(wp)          :: pcento = -9999._wp

  ! store raw strings too
  character(len=20) :: maxn_str  = ""
  character(len=20) :: kstop_str = ""
  character(len=20) :: pcento_str = ""

 end type run_config

 ! -------------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------------
 
 type :: fuse_info
   type(mpi_info)   :: mpi
   type(space_info) :: space
   type(time_info)  :: time
   type(snow_info)  :: snow
   type(mizu_info)  :: mrout
   type(topo_info)  :: ntopo
   type(remap_info) :: remap
   type(file_info)  :: files
   type(run_config) :: config
 end type fuse_info

end module info_types
