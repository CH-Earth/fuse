FUNCTION FUNCTN(NOPT,A)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2009
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Wrapper for SCE (used to compute the objective function)
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE fuse_metric_module                                ! run model and compute the metric chosen as objective function
USE multiforce, only: ncid_forc                       ! NetCDF forcing file ID
USE fuse_fileManager,only:METRIC, TRANSFO             ! metric and transformation requested in the filemanager

IMPLICIT NONE
! input
INTEGER(I4B)                           :: NOPT        ! number of parameters
REAL(MSP), DIMENSION(100), INTENT(IN)  :: A            ! model parameter set - can be bumped up to 100 elements

! internal
REAL(SP), DIMENSION(:), ALLOCATABLE    :: SCE_PAR     ! sce parameter set
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate
INTEGER(I4B)                           :: ERR         ! error code for fuse_metric
CHARACTER(LEN=256)                     :: MESSAGE     ! error message for fuse_metric
LOGICAL(LGT)                           :: OUTPUT_FLAG ! .TRUE. = write model time series
REAL(SP)                               :: METRIC_VAL  ! value of the metric chosen as objective function

! output
REAL(MSP)                              :: FUNCTN      ! objective function value

! ---------------------------------------------------------------------------------------
! get SCE parameter set
ALLOCATE(SCE_PAR(NOPT), STAT=IERR); IF (IERR.NE.0) STOP ' problem allocating space '
SCE_PAR(1:NOPT) = A(1:NOPT)  ! convert from MSP used in SCE to SP used in FUSE

OUTPUT_FLAG=.FALSE.   ! do not produce *runs.nc files only, param.nc files

CALL FUSE_METRIC(SCE_PAR,.FALSE.,NCID_FORC,METRIC_VAL,OUTPUT_FLAG,1) ! 2nd argument FALSE, always return METRIC value

! deallocate parameter set
DEALLOCATE(SCE_PAR, STAT=IERR); IF (IERR.NE.0) STOP ' problem deallocating space '
print *, 'METRIC_VAL [Metric:',METRIC,' / Transfo:',TRANSFO,'] =', METRIC_VAL

! save objective function value: SCE is a minimization algorithm
IF (METRIC=="KGE" .OR. METRIC=="KGEP" .OR. METRIC=="NSE") THEN
  FUNCTN = -METRIC_VAL
ELSE IF (METRIC=="MAE" .OR. METRIC=="RMSE" ) THEN
  FUNCTN = METRIC_VAL
ELSE 
   STOP 'The requested metric is not available in metrics module'
END IF

! ---------------------------------------------------------------------------------------
END FUNCTION FUNCTN
