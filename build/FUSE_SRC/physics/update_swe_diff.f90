module update_swe_DIFF_MODULE

  USE model_defn                                        ! model definition structure
  USE model_defnames                                    ! integer model definitions
  USE globaldata, only : NA_VALUE_SP                    ! missing vale

  implicit none

  private
  public :: update_swe_diff

contains

  ! ---------------------------------------------------------------------------------------
  pure logical function is_leap_year(y)
   integer, intent(in) :: y
   is_leap_year = (mod(y,4) == 0 .and. (mod(y,100) /= 0 .or. mod(y,400) == 0))
  end function is_leap_year
  ! ---------------------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE UPDATE_SWE_DIFF(fuseStruct, DT, want_dparam)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Brian Henn, as part of FUSE snow model implementation, 6/2013
  ! Based on subroutines QSATEXCESS and UPDATSTATE, by Martyn Clark
  !
  ! Modified by Nans Addor to enable distributed modeling, 9/2016
  !
  ! Modified by Martyn Clark to extend to a differentiable model, 12/2025
  !
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the snow accumulation and melt from forcing data
  ! Then updates the SWE band states based on the fluxes
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                               ! variable types, etc. (includes PI)
  USE data_types, only: parent                             ! fuse parent data type
  use smoothers,  only: smax, dsmax                        ! max smoothers
  use smoothers,  only: smin, dsmin                        ! min smoothers (based on smax, dsmax)
  use smoothers,  only: sigmoid, dsigmoid                  ! sigmoid smoothers
  USE globaldata, only: NP => NPAR_SNOW                    ! number of snow parameters
  USE globaldata, only: iMBASE, iMFMAX, iMFMIN, iPXTEMP, iOPG, iLAPSE, &  ! indices in vectors
                        iPERR ! not a snow parameter but used in the snow model
  USE multibands, only: N_BANDS                            ! number of elevation bands
  IMPLICIT NONE
  ! input
  type(parent) , intent(inout)       :: fuseStruct         ! parent fuse data structure
  REAL(SP), INTENT(IN)               :: DT                 ! length of the time step
  logical(lgt), intent(in), optional :: want_dparam        ! if we want parameter derivatives
  ! ----- internal variables -----------------------------------------------------------------------------
  ! general
  INTEGER(I4B)                       :: ISNW               ! loop through snow model bands
  REAL(SP)                           :: DZ                 ! vert. distance from forcing
  real(sp)                           :: SWE_prev           ! SWE at start of band update (mm)
  ! melt factor
  LOGICAL(LGT)                       :: LEAP               ! leap year flag
  REAL(SP)                           :: JDAY               ! Julian day of year
  integer(i4b)                       :: days_in_year       ! number of days in year (365 or 366)
  integer(i4b)                       :: phase_shift        ! shift in sine curve in days (80 or 81)
  real(sp)                           :: season01           ! seasonal cycle scaled to [0,1]
  REAL(SP)                           :: MF                 ! melt factor (mm/deg.C-6hr) -- NOTE: check units
  ! adjusted precipitation (after precipitation multiplier)
  real(sp), parameter                :: ms_mult=1.e-4_sp   ! smoothing in smax function (additive precip error)
  real(sp)                           :: precip_adj         ! adjusted precipitation (after multiplicative/additive error)
  ! temperature lapse (simple)
  real(sp)                           :: xLapse             ! scaled temperature lapse rate
  REAL(SP)                           :: TEMP_Z             ! band temperature at timestep
  ! orographic precipitation multiplier (OPG)
  real(sp)                           :: xOPG               ! DZ * MPARAM%OPG/1000 -- scaled OPG (dimensionless)
  real(sp)                           :: gate               ! hard [0,1] gate on DZ
  real(sp)                           :: fpos               ! positive-side formula: 1 + x
  real(sp)                           :: fneg               ! megative-side formula: 1/(1-x)
  real(sp)                           :: inv                ! 1-x: demominator in negative-side formula: 1/(1-x)
  real(sp)                           :: inv_safe           ! safe denominator: max(1-x, eps_inv)
  real(sp), parameter                :: eps_inv=1.e-6_sp   ! denominator floor: dimensionless
  real(sp)                           :: OPG_mult           ! final OPG multiplier
  REAL(SP)                           :: PRECIP_Z           ! band precipitation at timestep
  ! partition rain from snow
  real(sp)                           :: fsnow              ! fraction of precip falling as snow (0–1)
  real(sp)                           :: snow               ! snowfall rate (mm/day) for this band
  real(sp)                           :: rain               ! rainfall rate (mm/day) for this band
  real(sp), parameter                :: beta_px=0.01_sp    ! sigmoid width for snow/rain partition (degC)
  ! snowmelt
  real(sp), parameter                :: ms_temp=1.e-4_sp   ! smoothing in smax function (temperature)
  real(sp)                           :: posTemp            ! positive-part temperature term used for melt (degC), smoothed
  real(sp)                           :: potMelt            ! potential melt rate before capping (mm/day)
  real(sp)                           :: meltCap            ! maximum feasible melt rate from availability (mm/day)
  real(sp)                           :: snowmelt           ! final (capped) melt rate (mm/day)
  real(sp)                           :: swe_eps=1.e-12_sp  ! small value for the derivative switch in u_swe clamp
  real(sp)                           :: u_swe              ! pre-clamp SWE update
  integer(i4b), parameter :: cumdays0(12) = [ &            ! cumulative days before the start of each month
   0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 ]
  integer(i4b)                       :: cumdays(12)        ! cumulative days adjust for leap year
  ! internal variables: paraneter derivatives
  logical(lgt)            :: comp_dparam  ! flag to compute parameter derivatives
  real(sp)                :: df_dz        ! precip partitioning
  real(sp)                :: active, dfpos_dOPG, dinv_dOPG, dfneg_dOPG, dmult_dOPG  ! OPG
  real(sp)                :: dMF(NP), dPadj(NP), dPrecZ(NP), dTempZ(NP)  ! derivative vectors
  real(sp)                :: dfsnow(NP), dsnow(NP), drain(NP)            ! derivative vectors
  real(sp)                :: g_pos, dposTemp(NP), dpotMelt(NP), dsnowmelt(NP)   ! derivative vectors
  real(sp)                :: g_u, dSWE(NP), dSWE_new(NP)   ! persist dSWE between timesteps for each band
  ! ---------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   TIMDAT => fuseStruct%time         , &  ! time information
   MFORCE => fuseStruct%force        , &  ! forcing data
   Z_FORC => fuseStruct%z_forcing    , &  ! elevation of the forcing data
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   MBANDS => fuseStruct%sbands       , &  ! elevation band variables: MBANDS(i)%var, MBANDS(i)info
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! ---------------------------------------------------------------------------------------
  ! snow accumulation and melt calculations for each band
  ! also calculates effective precipitation
  ! ---------------------------------------------------------------------------------------

  ! check the need to compute flux derivatives
  comp_dparam = .false.; if(present(want_dparam)) comp_dparam = want_dparam

  ! zero derivatives for fluxes constant over elevation bands
  if(comp_dparam)then
    dMF(:) = 0._sp; dPadj(:) = 0._sp
  endif

  ! ----- compute the melt factor ---------------------------------------------------------

  ! adjust cumulative days for leap year
  leap    = is_leap_year(timDat%IY)
  cumdays = cumdays0; if (leap) cumdays(3:12) = cumdays(3:12) + 1

  ! calculate day of year for melt factor calculation
  jday = cumdays(timDat%IM) + timDat%ID

  ! seasonal cycle scaled to [0,1]
  days_in_year = merge(366, 365, leap)
  phase_shift  = merge(81, 80, leap)   ! keeps peak timing aligned across leap/non-leap
  season01     = 0.5_sp * ( sin( (real(jday - phase_shift, sp) * 2._sp * PI) / real(days_in_year, sp) ) + 1._sp )

  ! melt factor calculations
  mf = MPARAM%MFMIN + season01*(MPARAM%MFMAX - MPARAM%MFMIN)

  ! compute derivatives
  if(comp_dparam)then

    ! NOTE: MF = (1−season01)*MFMIN + season01*MFMAX

    dMF(iMFMIN) = 1._sp - season01
    dMF(iMFMAX) = season01

  endif  ! computing derivatives

  ! ----- add error to the precipiation ---------------------------------------------------

  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e); precip_adj = smax(MFORCE%PPT + MPARAM%RFERR_ADD, 0._sp, ms_mult)   ! additive error
   CASE(iopt_multiplc_e); precip_adj = MFORCE%PPT*MPARAM%RFERR_MLT                      ! multiplicative error
   CASE DEFAULT; stop "swe_update_diff: unable to identify precip error model"
  END SELECT

  ! compute derivatives
  if(comp_dparam)then
   
     ! NOTE: parameter vector interprets theta(iPERR) as either RFERR_ADD or RFERR_MLT depending on SMODL%iRFERR

     SELECT CASE(SMODL%iRFERR)
      CASE(iopt_additive_e); dPadj(iPERR) = dsmax(MFORCE%PPT + MPARAM%RFERR_ADD, 0._sp, ms_mult)  ! additive error
      CASE(iopt_multiplc_e); dPadj(iPERR) = MFORCE%PPT                                       ! multiplicative error
      CASE DEFAULT; stop "swe_update_diff: unable to identify precip error model"
     END SELECT

  endif  ! computing derivatives

  ! ----- check OPG -----------------------------------------------------------------------
 
  if (MPARAM%OPG < 0._sp) then
    stop "swe_update_diff: OPG < 0 not allowed with hard-gate OPG scheme"
  end if

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  ! initialize effective precip
  M_FLUX%EFF_PPT = 0._sp

  ! check band rea fractions sum to 1
  if (abs(sum(MBANDS(:)%info%AF) - 1._sp) > 1.e-6_sp) stop "Band area fractions do not sum to 1"

  ! loop through model bands
  DO ISNW=1,N_BANDS
 
   ! save SWE
   SWE_prev = MBANDS(ISNW)%var%SWE
  
   ! zero derivatives for elevation band fluxes
   if(comp_dparam)then
    dPrecZ(:) = 0._sp; dTempZ(:) = 0._sp
    dfsnow(:) = 0._sp; dsnow(:) = 0._sp; drain(:) = 0._sp
    dposTemp(:)=0._sp; dpotMelt(:)=0._sp; dsnowmelt(:)=0._sp 
  endif

   ! copy the stored sensitivity of SWE from the previous timestep to propagate it forward
   if (comp_dparam) dSWE(:) = MBANDS(ISNW)%var%dSWE_dparam(:)

   ! --- use the Orographic Precipitation Gradient (OPG) to adjust precip for elevation ---

   ! dimensionless OPG
   DZ       = MBANDS(ISNW)%info%Z_MID - Z_FORC
   xOPG     = DZ * MPARAM%OPG / 1000._sp
   
   ! hard [0,1] gate by DZ sign (no smoothing): preserves original code from Henn et al.
   gate     = merge(1._sp, 0._sp, DZ >= 0._sp)   ! gate = 1 if DZ >= 0
   
   ! positive-side formula: 1 + x
   fpos     = 1._sp + xOPG
   
   ! negative-side formula: 1/(1-x), but with hard floor on denominator
   inv      = 1._sp - xOPG
   inv_safe = max(inv, eps_inv)     ! hard floor
   fneg     = 1._sp / inv_safe
   
   ! blended multiplier and band precip
   OPG_mult = gate * fpos + (1._sp - gate) * fneg
   PRECIP_Z = precip_adj * OPG_mult 

   ! compute derivatives
   if(comp_dparam)then

     ! derivative of fpos wrt OPG
     dfpos_dOPG = DZ  / 1000._sp

     ! derivative of fneg wrt OPG
     active     = merge(1._sp, 0._sp, inv >= eps_inv)  ! deriv is zero if inv is clamped at eps_inv
     dinv_dOPG  = -(DZ / 1000._sp) ! inv = 1 - xOPG,  xOPG = DZ*OPG/1000
     dfneg_dOPG = -(1._sp/(inv_safe*inv_safe)) * (active * dinv_dOPG)

     ! derivative of OPG_mult (ignore derivative of the hard gate)
     dmult_dOPG = gate*dfpos_dOPG + (1._sp-gate)*dfneg_dOPG

     ! final derivatives
     dPrecZ(:)    = dPadj(:) * OPG_mult
     dPrecZ(iOPG) = dPrecZ(iOPG) + precip_adj*dmult_dOPG

   endif  ! computing derivatives
   
   ! ----- use the temperature lapse rate to adjust temperature for elevation -------------

   xLapse = MPARAM%LAPSE/1000._sp          ! scaled temperature lapse rate
   TEMP_Z = MFORCE%TEMP + DZ*xLapse        ! adjust for elevation using lapse rate

   ! compute derivatives
   if(comp_dparam) dTempZ(iLAPSE) = DZ/1000._sp

   ! ----- calculate the (smoothed) snow accumulation -------------------------------------

   ! snowfall and rainfall fluxes
   fsnow = sigmoid(MPARAM%PXTEMP - TEMP_Z, beta_px) ! beta_px is the width, set small because originally a step function
   snow  = PRECIP_Z*fsnow
   rain  = PRECIP_Z*(1._sp - fsnow)

   MBANDS(ISNW)%var%SNOWACCMLTN = snow

   ! compute derivatives
   if(comp_dparam)then

     df_dz = dsigmoid(fsnow, beta_px)        ! d(fsnow)/d(z), z=PXTEMP - TEMP_Z
    
     dfsnow(iPXTEMP) = df_dz
     dfsnow(:) = dfsnow(:) - df_dz * dTempZ(:)   ! minus because z depends on -TEMP_Z
    
     dsnow(:) = dPrecZ(:)*fsnow + PRECIP_Z*dfsnow(:)
     drain(:) = dPrecZ(:)*(1._sp - fsnow) - PRECIP_Z*dfsnow(:)

   endif  ! computing derivatives

   ! ----- calculate the (smoothed) snow melt ---------------------------------------------

   ! potenital melt
   posTemp = smax(TEMP_Z - MPARAM%MBASE, 0._sp, ms_temp)   ! smoothed max(TEMP_Z - MPARAM%MBASE, 0)
   potMelt = MF*posTemp   !  mm day-1
 
   ! cap snowmelt
   meltCap  = SWE_prev/DT
   snowmelt = min(potMelt, meltCap) ! hard clamp: allow a kink at SWE=0 to avoid "ghost snow"
   MBANDS(ISNW)%var%SNOWMELT = snowmelt
  
   ! compute derivatives
   if(comp_dparam)then

     ! positive temperature: smoothed max(TEMP_Z - MPARAM%MBASE, 0)
     g_pos            = dsmax(TEMP_Z - MPARAM%MBASE, 0._sp, ms_temp)
     dposTemp(:)      = g_pos * dTempZ(:)
     dposTemp(iMBASE) = dposTemp(iMBASE) - g_pos

     ! potential melt
     dpotMelt(:) = dMF(:)*posTemp + MF*dposTemp(:)
     
     ! melt cap
     dsnowmelt(:) = merge(dpotMelt(:), dSWE(:)/DT, potMelt <= meltcap)

   endif  ! computing derivatives

   ! ----- update SWE ---------------------------------------------------------------------
  
   u_swe = SWE_prev + DT*(snow - snowmelt)
   MBANDS(ISNW)%var%SWE = max(u_swe, 0._sp)  ! hard clamp just removes numerical noise

   if(comp_dparam)then
     g_u = merge(1._sp, 0._sp, u_swe > swe_eps) ! sensitivities zero in snow free periods
     dSWE_new(:) = g_u * ( dSWE(:) + DT*(dsnow(:) - dsnowmelt(:)) )
     MBANDS(ISNW)%var%dSWE_dparam(:) = dSWE_new(:)
   endif

   ! ----- calculate effective precip (rain + melt)  ---------------------------------------

   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%info%AF * (rain + snowmelt)

   if(comp_dparam)then
     fuseStruct%df_dPar(1:NP)%EFF_PPT = fuseStruct%df_dPar(1:NP)%EFF_PPT + & 
                                        MBANDS(ISNW)%info%AF * (drain(:) + dsnowmelt(:))
   endif

  END DO  ! looping through elevation bands  

  end associate
  
  ! TEMPORARY: save the derivative as a "fake" loss function
  fuseStruct%dL_dPar(:)    = NA_VALUE_SP 
  fuseStruct%dL_dPar(1:NP) = fuseStruct%df_dPar(1:NP)%EFF_PPT

  END SUBROUTINE UPDATE_SWE_DIFF

end module update_swe_DIFF_MODULE
