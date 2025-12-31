module setup_model_definition_MODULE

  USE nrtype
  USE data_types, only: cli_options 
  USE data_types, only: domain_type 
  USE data_types, only: PARATT 

  implicit none

  private
  public :: setup_model_definition

contains

  subroutine setup_model_definition(opts, domain, APAR, BL, BU, err, message)

  ! access subroutines
  use uniquemodl_module, only: uniquemodl                   ! Defines unique strings for all FUSE models
  use GETPARMETA_module, only: GETPARMETA                   ! Reads parameter metadata from the parameter constraints file
  use selectmodl_module, only: selectmodl                   ! reads model control file
  use ASSIGN_STT_module, only: ASSIGN_STT                   ! state definitions:     data are stored in module model_defn
  use ASSIGN_FLX_module, only: ASSIGN_FLX                   ! flux definitions:      data are stored in module model_defn
  use ASSIGN_PAR_module, only: ASSIGN_PAR                   ! parameter definitions: data are stored in module multiparam
  use PAR_DERIVE_module, only: PAR_DERIVE                   ! Compute derived model parameters (bucket sizes, etc.)
  USE DEF_SSTATS_MODULE, only: DEF_SSTATS                   ! define summary statistics
  USE DEF_PARAMS_MODULE, only: DEF_PARAMS                   ! define model parameters
  USE DEF_OUTPUT_MODULE, only: DEF_OUTPUT                   ! define model output
  USE getpar_str_module, only: GETPAR_STR                   ! extracts parameter metadata

  ! data stored in legacy modules
  USE model_defn, only: NSTATE             ! number of state variables
  USE multiparam, only: NUMPAR             ! number of paramters for the current model
  USE multiparam, only: LPARAM             ! list of model parameters
  USE multiparam, only: MAXN               ! maximum number of function evaluations in SCE -- used for NUMPSET
  USE multiforce, only: NUMPSET            ! number of model parameter sets
  
  implicit none
  
  ! input
  type(cli_options)   , intent(in)                  :: opts            ! command line interface options
  type(domain_type)   , intent(inout)               :: domain          ! the domain structure that stores "everything"
  
  ! output
  real(sp)            , intent(out) , allocatable   :: aPar(:)         ! parameter vector
  real(sp)            , intent(out) , allocatable   :: BL(:), BU(:)    ! parameter bounds
  integer(i4b)        , intent(out)                 :: err             ! error code
  character(len=1024) , intent(out)                 :: message         ! error message
  
  ! ----- internal -----------------------------------------------------------------------
  INTEGER(I4B)                                      :: IPAR            ! parameter index
  INTEGER(I4B)                                      :: NMOD            ! number of models
  TYPE(PARATT)                                      :: PARAM_META      ! parameter metadata (model parameters)
  CHARACTER(LEN=1024)                               :: CMESSAGE        ! error message
  ! ----- output dimensions --------------------------------------------------------------
  integer(i4b) :: nx, ny, nt, nb, nSet, nPar
  ! ---------------------------------------------------------------------------------------
  associate(fmodel_id => domain%info%config%fmodel_id) ! use domain as truth where possible
  ! ---------------------------------------------------------------------------------------
  err=0; message='setup_model_definition/'

  ! ----- define characteristics of the current model -------------------------------------

  ! Define model attributes (valid for all models)
  CALL UNIQUEMODL(NMOD)            ! get nmod unique models: stored in module model_defn; NMOD is intent(out)
  CALL GETPARMETA(ERR,CMESSAGE)    ! read parameter metadata from constraints txt file (parameter bounds etc.)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! Identify a single model: FMODEL_ID is read from the control file and used to build string for zDecisions
  CALL SELECTMODL(FMODEL_ID,ERR=ERR,MESSAGE=CMESSAGE) ! FMODEL_ID is intent(in)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif
  
  ! Define list of states and parameters for the current model
  CALL ASSIGN_STT()        ! state definitions are stored in module model_defn
  CALL ASSIGN_FLX()        ! flux definitions are stored in module model_defn
  CALL ASSIGN_PAR()        ! parameter definitions are stored in module multiparam
 
  ! save information in global data structures
  domain%info%config%nState    = NSTATE   ! NSTATE is in module model_defn
  domain%info%config%nParam    = NUMPAR   ! NSTATE is in module multiparam
  domain%info%config%listParam = LPARAM(1:NUMPAR)   ! (performs allocation) LPARAM is in module multiparam

  ! Compute derived model parameters (bucket sizes, etc.)
  CALL PAR_DERIVE(ERR,CMESSAGE)
  if (err/=0)then; message=trim(message)//trim(cmessage); err=20; return; endif

  ! ----- initialize parameters, statistics, and output -----------------------------------

  ! get number of parameter sets
  ! will be used to define the parameter set dimension of the NetCDF files
  select case(opts%runmode)
    
    ! options that run with a single parameter set
    case('def', 'idx', 'opt'); NUMPSET = 1

    ! use NUMPSET =1.2MAXN since final number of parameter sets produced by SCE is unknown
    case('sce');               NUMPSET = int(1.2_sp * real(MAXN, sp)) 
      
    ! check
    err=20; message=trim(message)//'opts%runmode is unknown: '//trim(opts%runmode)

  end select

  ! save the number of parameter sets in the global domain structure
  domain%info%config%nSets = NUMPSET

  ! define NetCDF files

  ! assign dimensions (use domain data for provenance/clarity)

  nx = domain%info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = domain%info%space%ny_local
  nt = domain%info%time%nt_window
  nb = domain%info%snow%n_bands

  nSet = domain%info%config%nSets
  nPar = domain%info%config%nParam


  CALL DEF_PARAMS(nSet)                ! define model parameters
  CALL DEF_OUTPUT(nx,ny,nb,nPar,nt)    ! define model output time series (nPar used for parameter derivatives)
  
  CALL DEF_SSTATS()                    ! define summary statistics (REDEF)
 
  ! get parameter bounds and random numbers
  ALLOCATE(APAR(NUMPAR),BL(NUMPAR),BU(NUMPAR))
 
  DO IPAR=1,NUMPAR
   CALL GETPAR_STR(LPARAM(IPAR)%PARNAME,PARAM_META)
   BL(IPAR)   = PARAM_META%PARLOW  ! lower boundary
   BU(IPAR)   = PARAM_META%PARUPP  ! upper boundary
   APAR(IPAR) = PARAM_META%PARDEF  ! using default parameter values
  END DO

  end associate

  end subroutine setup_model_definition

end module setup_model_definition_MODULE
