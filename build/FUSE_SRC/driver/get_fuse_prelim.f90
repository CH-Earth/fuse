module get_fuse_prelim_MODULE

  USE nrtype
  USE data_types, only: cli_options, PARATT 

  implicit none

  private
  public :: get_fuse_prelim

contains

  subroutine get_fuse_prelim(opts, APAR, BL, BU, err, message)

  ! access subroutines
  use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_strerror
  USE fuse_fileManager,  only: fuse_SetDirsUndPhiles        ! sets directories and filenames
  USE selectmodl_module, only: selectmodl                   ! reads model control file
  USE force_info_module, only: force_info                   ! get forcing info for NetCDF files
  USE get_gforce_module, only: read_ginfo                   ! get dimension lengths from the NetCDF file
  USE get_mbands_module, only: GET_MBANDS_INFO              ! get elevation bands for snow modeling 
  USE GET_TIME_INDICES_MODULE, only: GET_TIME_INDICES       ! get time indices
  USE get_gforce_module, only: get_varid                    ! list of var ids
  USE DEF_PARAMS_MODULE, only: DEF_PARAMS                   ! define model parameters
  USE DEF_OUTPUT_MODULE, only: DEF_OUTPUT                   ! define model output
  USE getpar_str_module                                     ! extracts parameter metadata
  USE par_insert_module                                     ! inserts model parameters
  
  ! shared data
  USE fuse_fileManager, only: SETNGS_PATH,MBANDS_INFO,MBANDS_NC, &
                              OUTPUT_PATH,FORCINGINFO,INPUT_PATH,&
                              FMODEL_ID,&
                              suffix_forcing,suffix_elev_bands,&
                              numtim_sub_str,&
                              KSTOP_str, MAXN_str, PCENTO_str
  USE model_defn, only: FNAME_TEMPRY, FNAME_NETCDF_RUNS, FNAME_NETCDF_PARA
  USE multiforce, only: ncid_forc, forcefile, GRID_FLAG, SUB_PERIODS_FLAG
  USE multiForce, only: AFORCE, gForce, gForce_3d, ancilF, ancilF_3d, aValid
  USE multiForce, only: nSpat1, nSpat2, numtim_sub 
  USE multiforce, only: NUMPSET, numtim_sim
  USE multiState, only: gState, gState_3d 
  USE multiBands, only: MBANDS_VAR_4d, N_BANDS
  USE multiparam, only: LPARAM, NUMPAR
  USE multiparam, only: KSTOP, MAXN, PCENTO
  USE multiRoute, only: aRoute, AROUTE_3d 
  implicit none
  ! input
  type(cli_options)   , intent(in)                  :: opts            ! command line interface options
  ! output
  real(sp)            , intent(out) , allocatable   :: aPar(:)         ! parameter vector
  real(sp)            , intent(out) , allocatable   :: BL(:), BU(:)    ! parameter bounds
  integer(i4b)        , intent(out)                 :: err             ! error code
  character(len=1024) , intent(out)                 :: message         ! error message
  ! ----- internal -----------------------------------------------------------------------
  INTEGER(I4B)                                      :: IPAR            ! parameter index
  INTEGER(I4B)                                      :: NMOD            ! number of models
  CHARACTER(LEN=1024)                               :: ELEV_BANDS_NC   ! name of NetCDF file for elevation bands
  TYPE(PARATT)                                      :: PARAM_META      ! parameter metadata (model parameters)
  CHARACTER(LEN=64)                                 :: TAG             ! tag for output file
  CHARACTER(LEN=1024)                               :: CMESSAGE        ! error message
  ! ---------------------------------------------------------------------------------------
  associate(&
            run_mode => opts%runmode,      &   ! FUSE run mode
            ffm_file => opts%control_file, &   ! FUSE file manager file
            dom_id   => opts%domain_id     )   ! Domain ID
  ! ---------------------------------------------------------------------------------------
  err=0; message='get_fuse_prelim/'

  ! ----- set paths and file names --------------------------------------------------------
  
  ! set directories and filenames for control files
  call fuse_SetDirsUndPhiles(fuseFileManagerIn=ffm_file,err=err,message=cmessage)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! define name of forcing info and elevation band file
  forcefile= trim(dom_id)//suffix_forcing
  ELEV_BANDS_NC=trim(dom_id)//suffix_elev_bands

  ! define tag
  tag = ""; if(allocated(opts%tag)) tag = trim(opts%tag)

  ! temporary file name
  FNAME_TEMPRY = TRIM(OUTPUT_PATH)//TRIM(dom_id)//'_'//TRIM(FMODEL_ID)//'_'//trim(tag)

  ! files to which model run and parameter set will be saved
  FNAME_NETCDF_RUNS = trim(FNAME_TEMPRY)//'_runs_'//trim(run_mode)//'.nc'
  FNAME_NETCDF_PARA = trim(FNAME_TEMPRY)//'_para_'//trim(run_mode)//'.nc'

  ! convert characters to integer
  READ (MAXN_STR,*) MAXN     ! maximum number of trials before optimization is terminated
  READ (KSTOP_STR,*) KSTOP   ! number of shuffling loops the value must change by PCENTO (MAX=9)
  READ (PCENTO_STR,*) PCENTO ! the percentage

  PRINT *, 'Variables defined based on domain name:'
  PRINT *, 'forcefile:', TRIM(forcefile)
  PRINT *, 'ELEV_BANDS_NC:', TRIM(ELEV_BANDS_NC)
  
  ! ----- read information on numerical decisions, forcing files, and grid info -----------
  
  ! defines method/parameters used for numerical solution based on numerix file
  CALL GETNUMERIX(ERR,CMESSAGE)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! get forcing info from the txt file, ?? including NA_VALUE ??
  call force_info(err,cmessage)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  print *, 'Open forcing file:', trim(INPUT_PATH)//trim(forcefile)
  
  ! open NetCDF forcing file
  err = nf90_open(trim(INPUT_PATH)//trim(forcefile), nf90_nowrite, ncid_forc)
  if (err/=0)then; message=trim(message)//' nf90_open failed: '//trim(nf90_strerror(err)); return; endif
  PRINT *, 'NCID_FORC is', ncid_forc
  
  ! get the grid info (spatial and temporal dimensions) from the NetCDF file
  call read_ginfo(ncid_forc,err,cmessage)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! determine period over which to run and evaluate FUSE and their associated indices
  CALL GET_TIME_INDICES()
  
  ! check time indices are OK
  IF((.NOT.GRID_FLAG).AND.SUB_PERIODS_FLAG)THEN
    write(*,*) 'Error: in catchment mode:'
    write(*,*) 'FUSE must run over entire time series at once'
    write(*,*) 'Please set numtim_sub to -9999 in the filemanager (', trim(ffm_file),').'
    stop 1
  endif
 
  ! get elevation band info, in particular N_BANDS
  CALL GET_MBANDS_INFO(ELEV_BANDS_NC,err,cmessage) ! read band data from NetCDF file
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! Get NetCDF ID for each variable of the forcing file
  ! NOTE: populates data structures in multiforce
  call get_varID(ncid_forc, err, cmessage)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif

  ! ----- define characteristics of the current model -------------------------------------

  ! Define model attributes (valid for all models)
  CALL UNIQUEMODL(NMOD)            ! get nmod unique models
  CALL GETPARMETA(ERR,CMESSAGE)    ! read parameter metadata (parameter bounds etc.)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! Identify a single model
  CALL SELECTMODL(FMODEL_ID,ERR=ERR,MESSAGE=CMESSAGE)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! Define list of states and parameters for the current model
  CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
  CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
  CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam
  
  ! Compute derived model parameters (bucket sizes, etc.)
  CALL PAR_DERIVE(ERR,CMESSAGE)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif

  ! ----- initialize parameters, statistics, and output -----------------------------------

  ! get number of parameter sets
  ! will be used to define the parameter set dimension of the NetCDF files
  select case(opts%runmode)
  case('def', 'idx', 'opt'); NUMPSET=1
  case('sce');               NUMPSET=1.2*MAXN ! using 1.2MAXN since the final number of parameter sets produced by SCE is unknown
  end select

  CALL DEF_PARAMS(NUMPSET)                ! define model parameters (initial CREATE)
  CALL DEF_SSTATS()                            ! define summary statistics (REDEF)
  CALL DEF_OUTPUT(nSpat1,nSpat2,N_BANDS,numtim_sim)    ! define model output time series (REDEF)
 
  ! get parameter bounds and random numbers
  ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR))
 
  DO IPAR=1,NUMPAR
   CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
   BL(IPAR)   = PARAM_META%PARLOW  ! lower boundary
   BU(IPAR)   = PARAM_META%PARUPP  ! upper boundary
   APAR(IPAR) = PARAM_META%PARDEF  ! using default parameter values
  END DO

  ! ----- allocate space for time series, grids, and states -------------------------------

  ! allocate space for the basin/grid-average time series
  allocate(aForce(numtim_sub),aRoute(numtim_sub),stat=err)
  if(err/=0)then; message=trim(message)//'unable to allocate space for basin-average time series [aForce,aRoute]'; return; endif

  ! allocate space for the forcing grid and states
  allocate(ancilF(nspat1,nspat2), gForce(nspat1,nspat2), gState(nspat1,nspat2), stat=err)
  if(err/=0)then; message=trim(message)//'unable to allocate space for forcing grid GFORCE'; return; endif

  ! allocate space for the forcing grid and states with a time dimension - only for subperiod
  allocate(AROUTE_3d(nspat1,nspat2,numtim_sub), gState_3d(nspat1,nspat2,numtim_sub+1),gForce_3d(nspat1,nspat2,numtim_sub),aValid(nspat1,nspat2,numtim_sub),stat=err)
  if(err/=0)then; message=trim(message)//'unable to allocate space for 3d structure'; return; endif

  ! allocate space for elevation bands
  allocate(MBANDS_VAR_4d(nspat1,nspat2,N_BANDS,numtim_sub+1),stat=err)
  if(err/=0)then; message=trim(message)//'unable to allocate space for elevation bands'; return; endif

  end associate

  end subroutine get_fuse_prelim

end module get_fuse_prelim_MODULE
