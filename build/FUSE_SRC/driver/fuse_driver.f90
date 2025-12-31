PROGRAM DISTRIBUTED_DRIVER
! ---------------------------------------------------------------------------------------
! Creators:
! Martyn Clark, 2011
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Nans Addor to include distributed modeling, 9/2016
! Modified by Nans Addor to re-enable catchment-scale modeling, 4/2017
! Modified by Martyn Clark to modularize and simplify CLI, 12/2025
! ---------------------------------------------------------------------------------------
! Purpose:
! Driver program to run FUSE with a snow module as either at the catchment-scale or
! at the grid-scale
! ---------------------------------------------------------------------------------------
! data types
USE nrtype                                                ! variable types, etc.
USE data_types, only: cli_options                         ! domain (includes "everything")
USE data_types, only: domain_type                         ! command line interface options
USE multistats, only: PCOUNT                              ! counter 

! data
USE globaldata, only: isPrint
USE globaldata, only: ncid_out
USE multiparam, only: NUMPAR
USE multiforce, only: NUMPSET
USE multiforce, only: ncid_forc, GRID_FLAG, SUB_PERIODS_FLAG
USE multiForce, only: AFORCE, gForce, gForce_3d, aValid
USE multiState, only: gState, gState_3d
USE multiRoute, only: aRoute, AROUTE_3d

! model setup: external subroutines/functions
USE netcdf                                                       ! NetCDF library
USE parse_command_args_MODULE, only: parse_command_args          ! parse command line arguments
USE setup_domain_module, only: setup_domain                      ! initialize the model domain
USE setup_model_definition_module, only: setup_model_definition  ! setup the FUSE model configuration

! model run: external subroutines/functions
USE get_fparam_module, only: GET_PRE_PARAM, GET_SCE_PARAM ! read parameters from netcdf file
USE sce_driver_MODULE, only: sce_driver                   ! SCE optimization

! model simulation modules
USE fuse_rmse_module                                      ! run model and compute the root mean squared error

#ifdef __MPI__
  use mpi
#endif

IMPLICIT NONE

! error control
integer(i4b)                          :: err              ! error code
character(len=1024)                   :: message          ! error message

! command line arguments
type(cli_options)                     :: cli_opts         ! command line argument options

! parameter set; parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BL      ! vector of lower parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: BU      ! vector of upper parameter bounds
REAL(SP), DIMENSION(:), ALLOCATABLE    :: APAR    ! model parameter set

! function  evaluation
REAL(SP)                               :: RMSE    ! sim-obs differences

! model output
LOGICAL(LGT)                           :: OUTPUT_FLAG     ! .TRUE. = write time series output
INTEGER(I4B)                           :: ONEMOD=1        ! just specify one model

! global domaia type
type(domain_type)                      :: domain          ! includes "everything"

! ---------------------------------------------------------------------------------------
! ----- model preliminaries (initialize) ------------------------------------------------
! ---------------------------------------------------------------------------------------

! ----- initialize MPI ------------------------------------------------------------------

#ifdef __MPI__
  domain%info%mpi%enabled = .true.
  call MPI_Init(err); call MPI_check(err, "MPI_Init")
  call MPI_Comm_rank(MPI_COMM_WORLD, domain%info%mpi%rank, err);  call MPI_check(err, "MPI_Comm_rank")
  call MPI_Comm_size(MPI_COMM_WORLD, domain%info%mpi%nproc, err); call MPI_check(err, "MPI_Comm_size")
#else
  domain%info%mpi%enabled = .false.
  domain%info%mpi%rank    = 0
  domain%info%mpi%nproc   = 1
#endif

! suppress printing for higher ranks
if(domain%info%mpi%rank > 0) isPrint=.false.

! ----- parse command line arguments ----------------------------------------------------

call parse_command_args(cli_opts,err,message)
if(err/=0) stop trim(message)

if(isPrint)then
  print*, 'Control file = ', cli_opts%control_file
  print*, 'Run mode     = ', cli_opts%runmode
endif

! ----- initialize the model domain -----------------------------------------------------

! read forcing metadata (space/time/coords), apply MPI decomposition, and allocate domain arrays
call setup_domain(cli_opts, domain, err, message) 
if(err/=0) stop trim(message)

! ----- initialize model configurations -------------------------------------------------

! choose model, load parameter metadata, derive parameters, and define NetCDF output files
call setup_model_definition(cli_opts, domain, APAR, BL, BU, err, message)
if(err/=0) stop trim(message)

! ----- set initial counters ------------------------------------------------------------

! Define output and parameter files
ONEMOD=1                 ! one file per model (i.e., model dimension = 1)
PCOUNT=0                 ! counter for parameter sets evaluated (shared in MODULE multistats)

! ----- allocate data structures --------------------------------------------------------



! ---------------------------------------------------------------------------------------
! ----- run different FUSE modes --------------------------------------------------------
! ---------------------------------------------------------------------------------------

! select fuse mode
select case(cli_opts%runmode)

  ! ----- single parameter set ----------------------------------------------------------

  case('def', 'idx', 'opt')

    OUTPUT_FLAG=.TRUE.

    ! load specific parameter set given index in vector into APAR
    if (cli_opts%runmode=='idx') then
     CALL GET_PRE_PARAM(cli_opts%sets_file, cli_opts%indx, ONEMOD, NUMPAR, APAR)
    endif

    ! load best parameter set from NetCDF file into APAR
    if (cli_opts%runmode=='opt') then
     CALL GET_SCE_PARAM(cli_opts%sets_file, ONEMOD, NUMPAR, APAR)
    endif

    ! run FUSE
    CALL FUSE_RMSE(APAR, GRID_FLAG, NCID_FORC, RMSE, OUTPUT_FLAG, NUMPSET)


  ! ----- SCE calibration run -----------------------------------------------------------

  case('sce')

    call sce_driver(APAR, BL, BU)

  case default
    stop "cannot identify FUSE mode"
  
end select ! (FUSE mode) 
  
! ----- finalize ------------------------------------------------------------------------
  
! deallocate space
DEALLOCATE(APAR, BL, BU, stat=err)
if(err/=0)then; write(*,*) 'unable to deallocate space for parameter vectors'; stop; endif

DEALLOCATE(aForce, aRoute, aValid, stat=err)
if(err/=0)then; write(*,*) 'unable to deallocate space for catchment modeling'; stop; endif

DEALLOCATE(gForce, gState, gForce_3d, gState_3d, AROUTE_3d, stat=err)
if(err/=0)then; write(*,*) 'unable to deallocate space for grid modeling'; stop; endif

! close NetCDF files
IF(GRID_FLAG)THEN
  PRINT *, 'Closing forcing file'
  err = nf90_close(ncid_forc)
  if(err/=0)then; message=trim(message)//' nf90_close failed: '//trim(nf90_strerror(err)); return; endif
ENDIF

PRINT *, 'Closing output file'
err = nf90_close(ncid_out)
if(err/=0)then; message=trim(message)//' nf90_close failed: '//trim(nf90_strerror(err)); return; endif

PRINT *, 'Done'
STOP

! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------

contains

! ----- MPI checker ---------------------------------------------------------------------

subroutine mpi_check(ierr, callee)
#ifdef __MPI__
  use mpi
#endif
  implicit none
  integer(i4b), intent(in) :: ierr
  character(len=*), intent(in) :: callee
#ifdef __MPI__
  integer(i4b) :: slen, ierr2
  character(len=256) :: errstr
  if (ierr /= MPI_SUCCESS) then
    call MPI_Error_string(ierr, errstr, slen, ierr2)
    write(*,*) "MPI error at ", trim(callee), ": ", trim(errstr(1:slen))
    call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
  end if
#else
  ! serial build: do nothing
#endif
end subroutine mpi_check

END PROGRAM DISTRIBUTED_DRIVER
