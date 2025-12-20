module update_swe_DIFF_MODULE

  USE model_defn                                        ! model definition structure
  USE model_defnames                                    ! integer model definitions

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
  ! internal variables
  LOGICAL(LGT)                       :: LEAP               ! leap year flag
  REAL(SP)                           :: JDAY               ! Julian day of year
  integer(i4b)                       :: days_in_year       ! number of days in year (365 or 366)
  integer(i4b)                       :: phase_shift        ! shift in sine curve in days (80 or 81)
  real(sp)                           :: season01           ! seasonal cycle scaled to [0,1]
  REAL(SP)                           :: MF                 ! melt factor (mm/deg.C-6hr) -- NOTE: check units
  REAL(SP)                           :: DZ                 ! vert. distance from forcing
  real(sp)                           :: xOPG               ! scaled Orographic Precipitation Gradient (OPG)
  real(sp)                           :: xLapse             ! scaled temperature lapse rate
  real(sp)                           :: precip_adj         ! adjusted precipitation (after multiplicative/additive error)
  real(sp)                           :: xEXP               ! exponential scaling factor
  REAL(SP)                           :: PRECIP_Z           ! band precipitation at timestep
  REAL(SP)                           :: TEMP_Z             ! band temperature at timestep
  INTEGER(I4B)                       :: ISNW               ! loop through snow model bands
  real(sp)                           :: fsnow              ! fraction of precip falling as snow (0–1)
  real(sp)                           :: snow               ! snowfall rate (mm/day) for this band
  real(sp)                           :: rain               ! rainfall rate (mm/day) for this band
  real(sp), parameter                :: beta_px=0.1_sp     ! sigmoid width for snow/rain partition (degC)
  real(sp), parameter                :: ms=1.e-4_sp        ! smoothing in smax function
  real(sp)                           :: posTemp            ! positive-part temperature term used for melt (degC), smoothed
  real(sp)                           :: potMelt            ! potential melt rate before capping (mm/day)
  real(sp)                           :: meltCap            ! maximum feasible melt rate from availability (mm/day)
  real(sp)                           :: snowmelt           ! final (capped) melt rate (mm/day)
  integer(i4b), parameter :: cumdays0(12) = [ &            ! cumulative days before the start of each month
   0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 ]
  integer(i4b)                       :: cumdays(12)        ! cumulative days adjust for leap year
  ! internal variables: paraneter derivatives
  logical(lgt)            :: comp_dparam  ! flag to compute parameter derivatives
  real(sp)                :: SWE_prev     ! SWE at start of band update (mm)
  real(sp)                :: dMF(NP), dPadj(NP), dPrecZ(NP), dTempZ(NP)  ! derivative vectors
  real(sp)                :: dfsnow(NP), dsnow(NP), drain(NP)            ! derivative vectors
  real(sp)                :: df_dz
  real(sp)                :: dposTemp(NP), dpotMelt(NP), dmeltCap(NP), dsnowmelt(NP)
  real(sp)                :: dSWE(NP), dSWE_new(NP)   ! persist dSWE between timesteps for each band
  real(sp)                :: w_pot, w_cap             ! smooth-min weights
  real(sp)                :: g_pos, g_cap, g_u        ! dsmax factors
  real(sp)                :: u_swe                    ! pre-clamp SWE update
  ! ---------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   TIMDAT => fuseStruct%time         , &  ! time information
   MFORCE => fuseStruct%force        , &  ! forcing data
   Z_FORC => fuseStruct%z_forcing    , &  ! elevation of the forcing data
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   MBANDS => fuseStruct%sbands       , &  ! elevation band variables: MBANDS(i)%var%x
   DERIVS => fuseStruct%sbands       , &  ! parameter derivatives: DERIVS(i)%dx%x
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
   CASE(iopt_additive_e); precip_adj = smax(MFORCE%PPT + MPARAM%RFERR_ADD, 0._sp, ms)   ! additive error
   CASE(iopt_multiplc_e); precip_adj = MFORCE%PPT*MPARAM%RFERR_MLT                      ! multiplicative error
   CASE DEFAULT; stop "swe_update_diff: unable to identify precip error model"
  END SELECT

  ! compute derivatives
  if(comp_dparam)then
   
     ! NOTE: parameter vector interprets theta(iPERR) as either RFERR_ADD or RFERR_MLT depending on SMODL%iRFERR

     SELECT CASE(SMODL%iRFERR)
      CASE(iopt_additive_e); dPadj(iPERR) = dsmax(MFORCE%PPT + MPARAM%RFERR_ADD, 0._sp, ms)  ! additive error
      CASE(iopt_multiplc_e); dPadj(iPERR) = MFORCE%PPT                                       ! multiplicative error
      CASE DEFAULT; stop "swe_update_diff: unable to identify precip error model"
     END SELECT

  endif  ! computing derivatives

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  ! initialize effective precip
  M_FLUX%EFF_PPT = 0._sp

  ! check band rea fractions sum to 1
  if (abs(sum(MBANDS(:)%var%AF) - 1._sp) > 1.e-6_sp) stop "Band area fractions do not sum to 1"

  ! loop through model bands
  DO ISNW=1,N_BANDS
 
   ! save SWE
   SWE_prev = MBANDS(ISNW)%var%SWE
  
   ! zero derivatives for elevation band fluxes
   if(comp_dparam)then
    dPrecZ(:) = 0._sp; dTempZ(:) = 0._sp
    dfsnow(:) = 0._sp; dsnow(:) = 0._sp; drain(:) = 0._sp
    dposTemp(:)=0._sp; dpotMelt(:)=0._sp; dmeltCap(:)=0._sp; dsnowmelt(:)=0._sp 
  endif

   ! copy the stored sensitivity of SWE from the previous timestep to propagate it forward
   if (comp_dparam) dSWE(:) = DERIVS(ISNW)%dx%dSWE_dparam(:)

   ! --- use the Orographic Precipitation Gradient (OPG) to adjust precip for elevation ---

   DZ       = MBANDS(ISNW)%var%Z_MID - Z_FORC
   xOPG     = MPARAM%OPG / 1000._sp        ! scaled OPG
   xEXP     = exp(DZ * xOPG)               ! exponential scaling factor
   PRECIP_Z = precip_adj * xEXP  ! NOTE: modified from the original branch structure

   ! compute derivatives
   if(comp_dparam)then

     dPrecZ(:) = dPadj(:) * xEXP           ! chain from precip_adj
     dPrecZ(iOPG) = dPrecZ(iOPG) + PRECIP_Z * (DZ/1000._sp)

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
   posTemp = smax(TEMP_Z - MPARAM%MBASE, 0._sp, ms)   ! smoothed max(TEMP_Z - MPARAM%MBASE, 0)
   potMelt = MF*posTemp   !  mm day-1
  
   ! melt capped by availability of snow
   meltCap = smax(snow + SWE_prev/DT, 0._sp, ms)

   ! smooth snowmelt
   snowmelt = -smax(-potMelt, -meltCap, ms)   ! smooth min(potMelt, meltCap)
   MBANDS(ISNW)%var%SNOWMELT = snowmelt
  
   ! compute derivatives
   if(comp_dparam)then

     ! positive temperature: smoothed max(TEMP_Z - MPARAM%MBASE, 0)
     g_pos            = dsmax(TEMP_Z - MPARAM%MBASE, 0._sp, ms)
     dposTemp(:)      = g_pos * dTempZ(:)
     dposTemp(iMBASE) = dposTemp(iMBASE) - g_pos

     ! potential melt
     dpotMelt(:) = dMF(:)*posTemp + MF*dposTemp(:)
     
     ! melt cap
     g_cap   = dsmax(snow + SWE_prev/DT, 0._sp, ms)
     dmeltCap(:) = g_cap * (dsnow(:) + dSWE(:)/DT)

     ! cap on snowmelt: smooth min weights
     w_pot        = dsmax(-potMelt, -meltCap, ms)   ! ∂snowmelt/∂potMelt -- NOTE: minus sign cancels
     w_cap        = 1._sp - w_pot                   ! ∂snowmelt/∂meltCap
     dsnowmelt(:) = w_pot*dpotMelt(:) + w_cap*dmeltCap(:)

   endif  ! computing derivatives

   ! ----- update SWE ---------------------------------------------------------------------
   
   u_swe = SWE_prev + DT*(snow - snowmelt)
   MBANDS(ISNW)%var%SWE = smax(u_swe, 0._sp, ms)

   if(comp_dparam)then
     g_u = dsmax(u_swe, 0._sp, ms)
     dSWE_new(:) = g_u * ( dSWE(:) + DT*(dsnow(:) - dsnowmelt(:)) )
     DERIVS(ISNW)%dx%dSWE_dparam(:) = dSWE_new(:)
   endif

   ! ----- calculate effective precip (rain + melt)  ---------------------------------------

   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%var%AF * (rain + snowmelt)
 
   if(comp_dparam)then
     DERIVS(ISNW)%dx%dEffP_dParam(1:NP) = DERIVS(ISNW)%dx%dEffP_dParam(1:NP) + & 
                                          MBANDS(ISNW)%var%AF * (drain(:) + dsnowmelt(:))
   endif

  END DO  ! looping through elevation bands  
 
  end associate

  END SUBROUTINE UPDATE_SWE_DIFF

end module update_swe_DIFF_MODULE
