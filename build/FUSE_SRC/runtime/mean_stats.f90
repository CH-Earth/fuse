SUBROUTINE MEAN_STATS()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! Modified by Nans Addor to deal with NA values in QOBS, 2016
! Modified by Cyril Thébault to allow different metrics as objective function, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes summary statistics from model simulations
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multistats -- summary statistics stored in MODULE multistats
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE metrics                                           ! available metrics and transformations
USE fuse_fileManager,only:METRIC, TRANSFO             ! metric and transformation requested in the filemanager
! FUSE modules
USE multiforce                                        ! model forcing data (obs streamflow)
USE multiroute                                        ! routed runoff
USE multi_flux                                        ! fluxes
USE multistats                                        ! summary statistics
USE model_numerix                                     ! model numerix parameters and data

IMPLICIT NONE
! internal
INTEGER(I4B)                           :: I           ! looping
INTEGER(I4B)                           :: NS          ! number of samples
INTEGER(I4B)                           :: IERR        ! error code for allocate/deallocate statements
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QOBS        ! observed runoff - whole time series
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QOBS_AVAIL  ! observed runoff - only time steps with QOBS available
INTEGER(I4B)                           :: NUM_AVAIL   ! number of time steps with QOBS available
LOGICAL(SP), DIMENSION(:), ALLOCATABLE :: QOBS_MASK   ! boolean: is QOBS available?
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QSIM        ! simulated runoff - whole time series
REAL(SP), DIMENSION(:), ALLOCATABLE    :: QSIM_AVAIL  ! simulated runoff - only time steps with QOBS available
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DOBS        ! observed runoff anomalies
REAL(SP), DIMENSION(:), ALLOCATABLE    :: DSIM        ! simulated runoff anomalies
REAL(SP), DIMENSION(:), ALLOCATABLE    :: RAWD        ! observed-simulated differences in flow
REAL(SP), DIMENSION(:), ALLOCATABLE    :: LOGD        ! observed-simulated differences in LOG flow
REAL(SP)                               :: XB_OBS      ! mean observed runoff
REAL(SP)                               :: XB_SIM      ! mean simulated runoff
REAL(SP)                               :: SS_OBS      ! sum of squared observed runoff anomalies
REAL(SP)                               :: SS_SIM      ! sum of squared simulated runoff anomalies
REAL(SP)                               :: SS_LOBS     ! sum of squared lagged differences in observed runoff
REAL(SP)                               :: SS_LSIM     ! sum of squared lagged differences in simulated runoff
REAL(SP)                               :: SS_RAW      ! sum of squared differences in observed - simulated
REAL(SP)                               :: SS_LOG      ! sum of squared differences in LOG observed - LOG simulated
REAL(SP)                               :: NO_ZERO     ! avoid divide by zero

! ---------------------------------------------------------------------------------------
! (1) PRELIMINARIES
! ---------------------------------------------------------------------------------------
! define sample size
NS =  eval_end-eval_beg+1
PRINT *, 'Number of time steps in evaluation period (EP) = ', NS

! allocate space for observed and simulated runoff
ALLOCATE(QOBS(NS),QOBS_MASK(NS),QSIM(NS),STAT=IERR)
IF (IERR.NE.0) STOP ' PROBLEM ALLOCATING SPACE IN MEAN_STATS.F90 '

! extract OBS and SIM for evaluation period, note that sim_beg, eval_beg, sim_end,
! eval_end are all with respect to julian_day_input
QSIM = AROUTE_3d(1,1,eval_beg-sim_beg+1:eval_end-sim_beg+1)%Q_ROUTED
QOBS = aValid(1,1,eval_beg-sim_beg+1:eval_end-sim_beg+1)%OBSQ

! check for missing QOBS values
QOBS_MASK = QOBS.ne.REAL(NA_VALUE, KIND(SP)) ! find the time steps for which QOBS is available
NUM_AVAIL = COUNT(QOBS_MASK)	! number of time steps for which QOBS is available

PRINT *, 'Number of time steps with observed streamflow in EP = ', NUM_AVAIL

IF (NUM_AVAIL.EQ.0) THEN

  PRINT *, 'Skiping computation of error statistics because no observed streamflow data'
  MSTATS%NASH_SUTT=-9999
  MSTATS%RAW_RMSE=-9999
  MSTATS%KGE=-9999
  MSTATS%METRIC_VAL=-9999

ELSE

  ! extract elements from QOBS and QSIM for time steps with QOBS available
  ALLOCATE(QOBS_AVAIL(NUM_AVAIL),QSIM_AVAIL(NUM_AVAIL),DOBS(NUM_AVAIL),DSIM(NUM_AVAIL),RAWD(NUM_AVAIL),LOGD(NUM_AVAIL),STAT=IERR)
  IF (IERR /= 0) STOP ' PROBLEM ALLOCATING SPACE FOR AVAILABLE DATA IN MEAN_STATS.F90 '

  
  QOBS_AVAIL=PACK(QOBS,QOBS_MASK,QOBS_AVAIL)  ! moves QOBS time steps indicated by QOBS_MASK to QOBS_AVAIL,
  											                      ! if no values is missing (i.e. NS = NUM_AVAIL) then QOBS_AVAIL
  											                      ! should be a copy of QOBS
  QSIM_AVAIL=PACK(QSIM,QOBS_MASK,QSIM_AVAIL)  ! moves QSIM time steps indicated by QOBS_MASK to QSIM_AVAIL
  											                      ! if no values is missing (i.e. NS = NUM_AVAIL) then QSIM_AVAIL
  										                      	! should be a copy of QSIM
                                                                 
  ! compute mean
  XB_OBS  = SUM(QOBS_AVAIL(:)) / INT(NUM_AVAIL, KIND(SP))
  XB_SIM  = SUM(QSIM_AVAIL(:)) / INT(NUM_AVAIL, KIND(SP))
  
  ! define NO_ZERO as 1% of the observed mean flow
  NO_ZERO = XB_OBS/100

  ! compute the sum of squares of simulated and observed vectors
  DOBS(:) = QOBS_AVAIL(:) - XB_OBS
  DSIM(:) = QSIM_AVAIL(:) - XB_SIM
  SS_OBS  = DOT_PRODUCT(DOBS,DOBS)  ! = SUM( DOBS(:)*DOBS(:) )
  SS_SIM  = DOT_PRODUCT(DSIM,DSIM)  ! = SUM( DSIM(:)*DSIM(:) )

  ! compute the sum of squares of lagged differences
  SS_LOBS = DOT_PRODUCT(DOBS(2:NUM_AVAIL),DOBS(1:NUM_AVAIL-1))
  SS_LSIM = DOT_PRODUCT(DSIM(2:NUM_AVAIL),DSIM(1:NUM_AVAIL-1))
  ! compute sum of squared differences between model and observations
  RAWD(:) = QSIM_AVAIL(:) - QOBS_AVAIL(:)
  LOGD(:) = LOG(QSIM_AVAIL(:)+NO_ZERO) - LOG(QOBS_AVAIL(:)+NO_ZERO)
  SS_RAW  = DOT_PRODUCT(RAWD,RAWD)  ! = SUM( RAWD(:)*RAWD(:) )
  SS_LOG  = DOT_PRODUCT(LOGD,LOGD)  ! = SUM( LOGD(:)*LOGD(:) )
  
  

  ! ---------------------------------------------------------------------------------------
  ! (2) COMPUTE ERROR STATISTICS
  ! ---------------------------------------------------------------------------------------
  ! compute the mean
  MSTATS%QOBS_MEAN = XB_OBS
  MSTATS%QSIM_MEAN = XB_SIM
  ! compute the coefficient of variation
  MSTATS%QOBS_CVAR = SQRT( SS_OBS / INT(NUM_AVAIL-1, KIND(SP)) ) / (XB_OBS+NO_ZERO)
  MSTATS%QSIM_CVAR = SQRT( SS_SIM / INT(NUM_AVAIL-1, KIND(SP)) ) / (XB_SIM+NO_ZERO)
  ! compute the lag-1 correlation coefficient
  MSTATS%QOBS_LAG1 = SS_LOBS / (SQRT(SS_OBS*SS_OBS)+NO_ZERO)
  MSTATS%QSIM_LAG1 = SS_LSIM / (SQRT(SS_SIM*SS_SIM)+NO_ZERO)
  
  ! Compute RMSE using the metrics module
  MSTATS%RAW_RMSE = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, '1')   ! No transformation
    
  ! Compute log RMSE using the metrics module
  MSTATS%LOG_RMSE = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, 'log') ! Log transformation

  ! Compute NSE using the metrics module
  MSTATS%NASH_SUTT = get_NSE(QOBS_AVAIL, QSIM_AVAIL, '1')   ! No transformation
  
  ! Compute KGE using the metrics module
  MSTATS%KGE = get_KGE(QOBS_AVAIL, QSIM_AVAIL, '1')         ! No transformation

  ! Compute KGEp using the metrics module
  MSTATS%KGEP = get_KGEP(QOBS_AVAIL, QSIM_AVAIL, '1')       ! No transformation
  
  ! Compute MAE using the metrics module
  MSTATS%MAE = get_MAE(QOBS_AVAIL, QSIM_AVAIL, '1')         ! No transformation
  
  ! Compute the metric chosen as objective function using the metrics module
  
  IF (METRIC == "KGE") THEN
    MSTATS%METRIC_VAL = get_KGE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO) 
  ELSE IF (METRIC == "KGEP") THEN
    MSTATS%METRIC_VAL = get_KGEP(QOBS_AVAIL, QSIM_AVAIL, TRANSFO) 
  ELSE IF (METRIC == "NSE") THEN
    MSTATS%METRIC_VAL = get_NSE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO) 
  ELSE IF (METRIC == "RMSE") THEN
    MSTATS%METRIC_VAL = get_RMSE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO) 
  ELSE IF (METRIC == "MAE") THEN
    MSTATS%METRIC_VAL = get_MAE(QOBS_AVAIL, QSIM_AVAIL, TRANSFO) 
  ELSE
    STOP 'The requested metric is not available in metrics module'
  END IF
  ! ---------------------------------------------------------------------------------------
  DEALLOCATE(QOBS,QOBS_AVAIL,QSIM,QSIM_AVAIL,DOBS,DSIM,RAWD,LOGD,STAT=IERR)
  IF (IERR.NE.0) STOP ' PROBLEM DEALLOCATING SPACE IN MEAN_STATS.F90 '

END IF

PRINT *, 'NSE = ',          MSTATS%NASH_SUTT
PRINT *, 'KGE = ',          MSTATS%KGE
PRINT *, 'KGEP = ',         MSTATS%KGEP
PRINT *, 'MAE = ',          MSTATS%MAE
PRINT *, 'RAW_RMSE = ',     MSTATS%RAW_RMSE
PRINT *, 'LOG_RMSE = ',     MSTATS%LOG_RMSE
PRINT *, 'METRIC_VAL [Metric:',METRIC,' / Transfo:',TRANSFO,'] =',   MSTATS%METRIC_VAL

! ---------------------------------------------------------------------------------------
! (3§) COMPUTE STATISTICS ON NUMERICAL ACCURACY AND EFFICIENCY
! ---------------------------------------------------------------------------------------
! compute RMSE between "more accurate" and "less accurate" solutions
! NA: note that this section uses only simulated discharge, hence uses NS, QSIM and QOBS
! instead of NUM_AVAIL, QOBS_AVAIL and QSIM_AVAIL, respectively
!QOBS = AROUTE(ISTART:NUMTIM_SIM)%Q_ACCURATE ! TODO: MISSING AT THE MOMENT
!RAWD(:) = QSIM(:) - QOBS(:)
!SS_RAW  = DOT_PRODUCT(RAWD,RAWD)    ! = SUM( RAWD(:)*RAWD(:) )
!MSTATS%NUM_RMSE = SQRT( SS_RAW / REAL(NS, KIND(SP)) )
! compute summary statistics for efficiency
MSTATS%NUM_FUNCS     = MSTATS%NUM_FUNCS     / REAL(NUMTIM_SIM, KIND(SP)) ! number of function calls
MSTATS%NUM_JACOBIAN  = MSTATS%NUM_JACOBIAN  / REAL(NUMTIM_SIM, KIND(SP)) ! number of times Jacobian is calculated
MSTATS%NUMSUB_ACCEPT = MSTATS%NUMSUB_ACCEPT / REAL(NUMTIM_SIM, KIND(SP)) ! number of sub-steps accepted (taken)
MSTATS%NUMSUB_REJECT = MSTATS%NUMSUB_REJECT / REAL(NUMTIM_SIM, KIND(SP)) ! number of sub-steps tried but rejected
MSTATS%NUMSUB_NOCONV = MSTATS%NUMSUB_NOCONV / REAL(NUMTIM_SIM, KIND(SP)) ! number of sub-steps tried that did not converge
! compute cumulative probability distributions
MSTATS%NUMSUB_PROB   = REAL(PRB_NSUBS(:), KIND(SP)) / REAL(NUMTIM_SIM, KIND(SP))

! ---------------------------------------------------------------------------------------
END SUBROUTINE MEAN_STATS
