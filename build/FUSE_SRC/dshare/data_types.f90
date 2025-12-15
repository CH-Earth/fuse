module data_types

 use nrtype
 use model_defn, only:NTDH_MAX

 ! --------------------------------------------------------------------------------------
 ! model time structure
 ! --------------------------------------------------------------------------------------
 TYPE M_TIME
  REAL(SP)                             :: STEP       ! (time interval to advance model states)
 END TYPE M_TIME

 ! --------------------------------------------------------------------------------------
 ! model forcing structures
 ! --------------------------------------------------------------------------------------
 
 ! the time data structure (will have no spatial dimension)
 TYPE TDATA
    INTEGER(I4B)                         :: IY         ! year
    INTEGER(I4B)                         :: IM         ! month
    INTEGER(I4B)                         :: ID         ! day
    INTEGER(I4B)                         :: IH         ! hour
    INTEGER(I4B)                         :: IMIN       ! minute
    REAL(SP)                             :: DSEC       ! second
    REAL(SP)                             :: DTIME      ! time in seconds since year dot
 ENDTYPE TDATA
 
 ! the response structure (will not have a spatial dimension)
 TYPE VDATA
    REAL(SP)                             :: OBSQ       ! observed runoff (mm day-1)
 END TYPE VDATA
 
 ! ancillary forcing variables used to compute ET (will have a spatial dimension)
 TYPE ADATA
    REAL(SP)                             :: AIRTEMP    ! air temperature (K)
    REAL(SP)                             :: SPECHUM    ! specific humidity (g/g)
    REAL(SP)                             :: AIRPRES    ! air pressure (Pa)
    REAL(SP)                             :: SWDOWN     ! downward sw radiation (W m-2)
    REAL(SP)                             :: NETRAD     ! net radiation (W m-2)
 END TYPE ADATA
 
 ! the forcing data structure (will have a spatial dimension)
 TYPE FDATA
    REAL(SP)                             :: PPT        ! water input: rain + melt (mm day-1)
    REAL(SP)                             :: TEMP       ! temperature for snow model (deg.C)
    REAL(SP)                             :: PET        ! energy input: potential ET (mm day-1)
 ENDTYPE FDATA

 ! --------------------------------------------------------------------------------------
 ! model state structure
 ! --------------------------------------------------------------------------------------
 TYPE STATEV
  ! snow layer
  REAL(SP)                             :: SWE_TOT    ! total storage as snow (mm)
  ! upper layer
  REAL(SP)                             :: WATR_1     ! total storage in layer1 (mm)
  REAL(SP)                             :: TENS_1     ! tension storage in layer1 (mm)
  REAL(SP)                             :: FREE_1     ! free storage in layer 1 (mm)
  REAL(SP)                             :: TENS_1A    ! storage in the recharge zone (mm)
  REAL(SP)                             :: TENS_1B    ! storage in the lower zone (mm)
  ! lower layer
  REAL(SP)                             :: WATR_2     ! total storage in layer2 (mm)
  REAL(SP)                             :: TENS_2     ! tension storage in layer2 (mm)
  REAL(SP)                             :: FREE_2     ! free storage in layer2 (mm)
  REAL(SP)                             :: FREE_2A    ! storage in the primary resvr (mm)
  REAL(SP)                             :: FREE_2B    ! storage in the secondary resvr (mm)
 END TYPE STATEV

 ! --------------------------------------------------------------------------------------
 ! model flux structure
 ! --------------------------------------------------------------------------------------
 TYPE FLUXES
  REAL(SP)                             :: EFF_PPT     ! effective precipitation (mm day-1)
  REAL(SP)                             :: SATAREA     ! saturated area (-)
  REAL(SP)                             :: QSURF       ! surface runoff (mm day-1)
  REAL(SP)                             :: EVAP_1A     ! evaporation from soil excess zone (mm day-1)
  REAL(SP)                             :: EVAP_1B     ! evaporation from soil recharge zone (mm day-1)
  REAL(SP)                             :: EVAP_1      ! evaporation from upper soil layer (mm day-1)
  REAL(SP)                             :: EVAP_2      ! evaporation from lower soil layer (mm day-1)
  REAL(SP)                             :: RCHR2EXCS   ! flow from recharge to excess (mm day-1)
  REAL(SP)                             :: TENS2FREE_1 ! flow from tension storage to free storage (mm day-1)
  REAL(SP)                             :: TENS2FREE_2 ! flow from tension storage to free storage (mm day-1)
  REAL(SP)                             :: QINTF_1     ! interflow from free water (mm day-1)
  REAL(SP)                             :: QPERC_12    ! percolation from upper to lower soil layers (mm day-1)
  REAL(SP)                             :: QBASE_2     ! baseflow (mm day-1)
  REAL(SP)                             :: QBASE_2A    ! baseflow from primary linear resvr (mm day-1)
  REAL(SP)                             :: QBASE_2B    ! baseflow from secondary linear resvr (mm day-1)
  REAL(SP)                             :: OFLOW_1     ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2     ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2A    ! bucket overflow (mm day-1)
  REAL(SP)                             :: OFLOW_2B    ! bucket overflow (mm day-1)
  REAL(SP)                             :: ERR_WATR_1  ! excessive extrapolation: total storage in layer1 (mm day-1)
  REAL(SP)                             :: ERR_TENS_1  ! excessive extrapolation: tension storage in layer1 (mm day-1)
  REAL(SP)                             :: ERR_FREE_1  ! excessive extrapolation: free storage in layer 1 (mm day-1)
  REAL(SP)                             :: ERR_TENS_1A ! excessive extrapolation: storage in the recharge zone (mm day-1)
  REAL(SP)                             :: ERR_TENS_1B ! excessive extrapolation: storage in the lower zone (mm day-1)
  REAL(SP)                             :: ERR_WATR_2  ! excessive extrapolation: total storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_TENS_2  ! excessive extrapolation: tension storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_FREE_2  ! excessive extrapolation: free storage in layer2 (mm day-1)
  REAL(SP)                             :: ERR_FREE_2A ! excessive extrapolation: storage in the primary resvr (mm day-1)
  REAL(SP)                             :: ERR_FREE_2B ! excessive extrapolation: storage in the secondary resvr (mm day-1)
  REAL(SP)                             :: CHK_TIME    ! time elapsed during time step (days)
 ENDTYPE FLUXES

 ! --------------------------------------------------------------------------------------
 ! model runoff structure
 ! --------------------------------------------------------------------------------------
 TYPE RUNOFF
  REAL(SP)                             :: Q_INSTNT   ! instantaneous runoff
  REAL(SP)                             :: Q_ROUTED   ! routed runoff
  REAL(SP)                             :: Q_ACCURATE ! "accurate" runoff estimate (mm day-1)
 END TYPE RUNOFF

 ! --------------------------------------------------------------------------------------
 ! parameter metadata
 ! --------------------------------------------------------------------------------------

 ! data structure to hold metadata for adjustable model parameters
 TYPE PARATT
  LOGICAL(LGT)                         :: PARFIT      ! flag to determine if parameter is fitted
  INTEGER(I4B)                         :: PARSTK      ! flag (0=deterministic, 1=stochastic)
  REAL(SP)                             :: PARDEF      ! default parameter set
  REAL(SP)                             :: PARLOW      ! lower limit of each parameter
  REAL(SP)                             :: PARUPP      ! upper limit of each parameter
  REAL(SP)                             :: FRSEED      ! fraction param space for "reasonable" bounds
  REAL(SP)                             :: PARSCL      ! typical scale of parameter
  INTEGER(I4B)                         :: PARVTN      ! method used for variable transformation
  INTEGER(I4B)                         :: PARDIS      ! parametric form of prob dist used for prior/hyper
  INTEGER(I4B)                         :: PARQTN      ! transformation applied before use of prob dist
  INTEGER(I4B)                         :: PARLAT      ! number of latent variables (0=onePerStep, -1=from data)
  INTEGER(I4B)                         :: PARMTH      ! imeth for all variables ???what is this???
  INTEGER(I4B)                         :: NPRIOR      ! number of prior/hyper-parameters
  CHARACTER(LEN=256)                   :: P_NAME      ! parameter name
  CHARACTER(LEN=256)                   :: CHILD1      ! name of 1st parameter child
  CHARACTER(LEN=256)                   :: CHILD2      ! name of 2nd parameter child
 END TYPE PARATT

 ! data structure to hold metadata for each parameter
 TYPE PARINFO
  ! rainfall error parameters (adjustable)
  TYPE(PARATT)                         :: RFERR_ADD   !  additive rainfall error (mm day-1)
  TYPE(PARATT)                         :: RFERR_MLT   ! multiplicative rainfall error (-)
  TYPE(PARATT)                         :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  TYPE(PARATT)                         :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  TYPE(PARATT)                         :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  TYPE(PARATT)                         :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  TYPE(PARATT)                         :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  TYPE(PARATT)                         :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  TYPE(PARATT)                         :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  TYPE(PARATT)                         :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  TYPE(PARATT)                         :: FRACTEN     ! frac total storage as tension storage (-)
  TYPE(PARATT)                         :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  TYPE(PARATT)                         :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  TYPE(PARATT)                         :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable)
  TYPE(PARATT)                         :: PERCRTE     ! percolation rate (mm day-1)
  TYPE(PARATT)                         :: PERCEXP     ! percolation exponent (-)
  TYPE(PARATT)                         :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  TYPE(PARATT)                         :: PERCFRAC    ! fraction of percolation to tension storage (-)
  TYPE(PARATT)                         :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  TYPE(PARATT)                         :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  TYPE(PARATT)                         :: BASERTE     ! baseflow rate (mm day-1)
  TYPE(PARATT)                         :: QB_POWR     ! baseflow exponent (-)
  TYPE(PARATT)                         :: QB_PRMS     ! baseflow depletion rate (day-1)
  TYPE(PARATT)                         :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  TYPE(PARATT)                         :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  TYPE(PARATT)                         :: SAREAMAX    ! maximum saturated area
  TYPE(PARATT)                         :: AXV_BEXP    ! ARNO/VIC "b" exponent
  TYPE(PARATT)                         :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  TYPE(PARATT)                         :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  TYPE(PARATT)                         :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model (adjustable)
  TYPE(PARATT)                         :: MBASE       ! base melt temperature (deg. C)
  TYPE(PARATT)                         :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  TYPE(PARATT)                         :: PXTEMP      ! rain-snow partition temperature (deg. C)
  TYPE(PARATT)                         :: OPG         ! precipitation gradient (-)
  TYPE(PARATT)                         :: LAPSE       ! temperature gradient (deg. C)
 ENDTYPE PARINFO

 ! --------------------------------------------------------------------------------------
 ! adjustable parameters
 ! --------------------------------------------------------------------------------------
 TYPE PARADJ
  ! rainfall error parameters (adjustable)
  REAL(SP)                             :: RFERR_ADD   ! additive rainfall error (mm day-1)
  REAL(SP)                             :: RFERR_MLT   ! multiplicative rainfall error (-)
  REAL(SP)                             :: RFH1_MEAN   ! hyper parameter1: mean rainfall multiplier (-)
  REAL(SP)                             :: RFH2_SDEV   ! hyper parameter2: sdev rainfall multiplier (-)
  REAL(SP)                             :: RH1P_MEAN   ! prior param1 of hyper param1: prior mean of hypermean
  REAL(SP)                             :: RH1P_SDEV   ! prior param2 of hyper param1: prior sdev of hypermean
  REAL(SP)                             :: RH2P_MEAN   ! prior param1 of hyper param2: lower bound of hypersdev
  REAL(SP)                             :: RH2P_SDEV   ! prior param2 of hyper param2: upper bound of hypersdev
  ! bucket sizes (adjustable)
  REAL(SP)                             :: MAXWATR_1   ! maximum total storage in layer1 (mm)
  REAL(SP)                             :: MAXWATR_2   ! maximum total storage in layer2 (mm)
  REAL(SP)                             :: FRACTEN     ! frac total storage as tension storage (-)
  REAL(SP)                             :: FRCHZNE     ! PRMS: frac tension storage in recharge zone (-)
  REAL(SP)                             :: FPRIMQB     ! SAC: fraction of baseflow in primary resvr (-)
  ! evaporation (adjustable)
  REAL(SP)                             :: RTFRAC1     ! fraction of roots in the upper layer (-)
  ! percolation (adjustable)
  REAL(SP)                             :: PERCRTE     ! percolation rate (mm day-1)
  REAL(SP)                             :: PERCEXP     ! percolation exponent (-)
  REAL(SP)                             :: SACPMLT     ! multiplier in the SAC model for dry lower layer (-)
  REAL(SP)                             :: SACPEXP     ! exponent in the SAC model for dry lower layer (-)
  REAL(SP)                             :: PERCFRAC    ! fraction of percolation to tension storage (-)
  REAL(SP)                             :: FRACLOWZ    ! fraction of soil excess to lower zone (-)
  ! interflow (adjustable)
  REAL(SP)                             :: IFLWRTE     ! interflow rate (mm day-1)
  ! baseflow (adjustable)
  REAL(SP)                             :: BASERTE     ! baseflow rate (mm day-1)
  REAL(SP)                             :: QB_POWR     ! baseflow exponent (-)
  REAL(SP)                             :: QB_PRMS     ! baseflow depletion rate (day-1)
  REAL(SP)                             :: QBRATE_2A   ! baseflow depletion rate for primary resvr (day-1)
  REAL(SP)                             :: QBRATE_2B   ! baseflow depletion rate for secondary resvr (day-1)
  ! surface runoff (adjustable)
  REAL(SP)                             :: SAREAMAX    ! maximum saturated area
  REAL(SP)                             :: AXV_BEXP    ! ARNO/VIC "b" exponent
  REAL(SP)                             :: LOGLAMB     ! mean value of the log-transformed topographic index (m)
  REAL(SP)                             :: TISHAPE     ! shape parameter for the topo index Gamma distribution (-)
  ! time delay in runoff
  REAL(SP)                             :: TIMEDELAY   ! time delay in runoff (days)
  ! snow model
  REAL(SP)                             :: MBASE       ! base melt temperature (deg. C)
  REAL(SP)                             :: MFMAX       ! maximum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(SP)                             :: MFMIN       ! minimum melt factor (mm melt deg C.-1 6hrs-1)
  REAL(SP)                             :: PXTEMP      ! rain-snow partition temperature (deg. C)
  REAL(SP)                             :: OPG         ! precipitation gradient (-)
  REAL(SP)                             :: LAPSE       ! temperature gradient (deg. C)
 END TYPE PARADJ

 ! --------------------------------------------------------------------------------------
 ! derived parameters
 ! --------------------------------------------------------------------------------------
 TYPE PARDVD
  ! bucket sizes (derived)
  REAL(SP)                             :: MAXTENS_1   ! maximum tension storage in layer1 (mm)
  REAL(SP)                             :: MAXTENS_2   ! maximum tension storage in layer2 (mm)
  REAL(SP)                             :: MAXFREE_1   ! maximum free storage in layer 1 (mm)
  REAL(SP)                             :: MAXFREE_2   ! maximum free storage in layer2 (mm)
  REAL(SP)                             :: MAXTENS_1A  ! maximum storage in the recharge zone (mm)
  REAL(SP)                             :: MAXTENS_1B  ! maximum storage in the lower zone (mm)
  REAL(SP)                             :: MAXFREE_2A  ! maximum storage in the primary resvr (mm)
  REAL(SP)                             :: MAXFREE_2B  ! maximum storage in the secondary resvr (mm)
  ! evaporation
  REAL(SP)                             :: RTFRAC2     ! fraction of roots in the lower layer (-)
  ! percolation/baseflow
  REAL(SP)                             :: QBSAT       ! baseflow at saturation
  ! surface runoff
  REAL(SP)                             :: POWLAMB     ! mean value of the power-transformed topographic index (m**(1/n))
  REAL(SP)                             :: MAXPOW      ! max value of the power-transformed topographic index (m**(1/n))
  ! routing
  REAL(SP), DIMENSION(NTDH_MAX)        :: FRAC_FUTURE ! fraction of runoff in future time steps
  INTEGER(I4B)                         :: NTDH_NEED   ! number of time-steps with non-zero routing contribution
 END TYPE PARDVD

 ! --------------------------------------------------------------------------------------
 ! list of parameters for a given model
 ! --------------------------------------------------------------------------------------
 TYPE PAR_ID
  CHARACTER(LEN=9)                     :: PARNAME     ! list of parameter names
 ENDTYPE PAR_ID

 ! --------------------------------------------------------------------------------------
 ! model statistics structure
 ! --------------------------------------------------------------------------------------
 TYPE SUMMARY
  ! DMSL diagnostix
  REAL(SP)                             :: VAR_RESIDUL   ! variance of the model residuals
  REAL(SP)                             :: LOGP_SIMULN   ! log density of the model simulation
  REAL(SP)                             :: JUMP_TAKEN    ! defines a jump in the MCMC production run
  ! comparisons between model output and observations
  REAL(SP)                             :: QOBS_MEAN     ! mean observed runoff (mm day-1)
  REAL(SP)                             :: QSIM_MEAN     ! mean simulated runoff (mm day-1)
  REAL(SP)                             :: QOBS_CVAR     ! coefficient of variation of observed runoff (-)
  REAL(SP)                             :: QSIM_CVAR     ! coefficient of variation of simulated runoff (-)
  REAL(SP)                             :: QOBS_LAG1     ! lag-1 correlation of observed runoff (-)
  REAL(SP)                             :: QSIM_LAG1     ! lag-1 correlation of simulated runoff (-)
  REAL(SP)                             :: RAW_RMSE      ! root-mean-squared-error of flow (mm day-1)
  REAL(SP)                             :: LOG_RMSE      ! root-mean-squared-error of LOG flow (mm day-1)
  REAL(SP)                             :: NASH_SUTT     ! Nash-Sutcliffe score
  ! attributes of model output
  REAL(SP)                             :: NUM_RMSE      ! error of the approximate solution
  REAL(SP)                             :: NUM_FUNCS     ! number of function calls
  REAL(SP)                             :: NUM_JACOBIAN  ! number of times Jacobian is calculated
  REAL(SP)                             :: NUMSUB_ACCEPT ! number of sub-steps taken
  REAL(SP)                             :: NUMSUB_REJECT ! number of sub-steps taken
  REAL(SP)                             :: NUMSUB_NOCONV ! number of sub-steps tried that did not converge
  INTEGER(I4B)                         :: MAXNUM_ITERNS ! maximum number of iterations in implicit scheme
  REAL(SP), DIMENSION(20)              :: NUMSUB_PROB   ! probability distribution for number of sub-steps
  ! error checking
  CHARACTER(LEN=1024)                  :: ERR_MESSAGE   ! error message
 ENDTYPE SUMMARY

 ! --------------------------------------------------------------------------------------
 ! parent FUSE structure
 ! --------------------------------------------------------------------------------------
 type parent
  type(tdata)                         :: time           ! time data
  type(fdata)                         :: force          ! model forcing data
  type(statev)                        :: state0         ! state variables (start of step)
  type(statev)                        :: state1         ! state variables (end of step)
  type(statev)                        :: dx_dt          ! time derivative in state variables
  type(fluxes)                        :: flux           ! fluxes
  type(fluxes), allocatable           :: df_dS(:)       ! derivative in fluxes w.r.t. states
  type(runoff)                        :: route          ! hillslope routing
  type(par_id)                        :: param_name     ! parameter names
  type(parinfo)                       :: param_meta     ! metadata on model parameters
  type(paradj)                        :: param_adjust   ! adjustable model parametrs
  type(pardvd)                        :: param_derive   ! derived model parameters
  type(summary)                       :: sim_stats      ! simulation statistics
 end type parent

end module data_types
