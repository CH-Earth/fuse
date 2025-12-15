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
  SUBROUTINE UPDATE_SWE_DIFF(fuseStruct, DT)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Brian Henn, as part of FUSE snow model implementation, 6/2013
  ! Based on subroutines QSATEXCESS and UPDATSTATE, by Martyn Clark
  !
  ! Modified by Nans Addor to enable distributed modeling, 9/2016
  !
  ! Modified by Martyn Clark to extend to a differentiable model, 9/2016
  !
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the snow accumulation and melt from forcing data
  ! Then updates the SWE band states based on the fluxes
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                               ! variable types, etc. (includes PI)
  USE data_types, only: parent                             ! fuse parent data type
  use smoothers,  only: smax, sigmoid                      ! max and sigmoid smoothers
  USE multibands   ! NOTE: include in fuseStruct           ! model basin band structure
  IMPLICIT NONE
  ! input
  type(parent) , intent(inout)       :: fuseStruct         ! parent fuse data structure
  REAL(SP), INTENT(IN)               :: DT                 ! length of the time step
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
  REAL(SP)                           :: PRECIP_Z           ! band precipitation at timestep
  REAL(SP)                           :: TEMP_Z             ! band temperature at timestep
  INTEGER(I4B)                       :: ISNW               ! loop through snow model bands
  real(sp)                           :: fsnow              ! fraction of precip falling as snow (0–1)
  real(sp)                           :: snow               ! snowfall rate (mm/day) for this band
  real(sp)                           :: rain               ! rainfall rate (mm/day) for this band
  real(sp), parameter                :: beta_px=10._sp     ! sigmoid sharpness for snow/rain partition (1/degC)
  real(sp)                           :: posTemp            ! positive-part temperature term used for melt (degC), smoothed
  real(sp)                           :: potMelt            ! potential melt rate before capping (mm/day)
  real(sp)                           :: meltCap            ! maximum feasible melt rate from availability (mm/day)
  real(sp)                           :: snowmelt           ! final (capped) melt rate (mm/day)
  integer(i4b), parameter :: cumdays0(12) = [ &            ! cumulative days before the start of each month
   0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 ]
  integer(i4b)                       :: cumdays(12)        ! cumulative days adjust for leap year
  ! ---------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   TIMDAT => fuseStruct%time         , &  ! fluxes
   MFORCE => fuseStruct%force        , &  ! fluxes
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! ---------------------------------------------------------------------------------------
  ! snow accumulation and melt calculations for each band
  ! also calculates effective precipitation
  ! ---------------------------------------------------------------------------------------

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

  ! ----- add error to the precipiation ---------------------------------------------------

  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e); precip_adj = MAX(0.0_sp, MFORCE%PPT + MPARAM%RFERR_ADD)  ! additive error
   CASE(iopt_multiplc_e); precip_adj = MFORCE%PPT*MPARAM%RFERR_MLT                 ! multiplicative error
   CASE DEFAULT; stop "swe_update_diff: unable to identify precip error model"
  END SELECT

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  ! initialize effective precip
  M_FLUX%EFF_PPT = 0._sp

  ! check band rea fractions sum to 1
  if (abs(sum(MBANDS(:)%AF) - 1._sp) > 1.e-6_sp) stop "Band area fractions do not sum to 1"

  ! loop through model bands
  DO ISNW=1,N_BANDS
  
   ! --- use the Orographic Precipitation Gradient (OPG) to adjust precip for elevation ---

   DZ       = MBANDS(ISNW)%Z_MID - Z_FORCING
   xOPG     = MPARAM%OPG / 1000._sp        ! scaled OPG
   PRECIP_Z = precip_adj * exp(DZ * xOPG)  ! NOTE: modified from the original branch structure

   ! ----- use the temperature lapse rate to adjust temperature for elevation -------------

   xLapse = MPARAM%LAPSE/1000._sp          ! scaled temperature lapse rate
   TEMP_Z = MFORCE%TEMP + DZ*xLapse        ! adjust for elevation using lapse rate

   ! ----- calculate the (smoothed) snow accumulation -------------------------------------

   ! snowfall and rainfall fluxes
   fsnow = sigmoid(MPARAM%PXTEMP - TEMP_Z, beta_px) ! beta_px is the sharpness, set large because originally a step function
   snow  = PRECIP_Z*fsnow
   rain  = PRECIP_Z*(1._sp - fsnow)

   MBANDS(ISNW)%SNOWACCMLTN = snow

   ! ----- calculate the (smoothed) snow melt ---------------------------------------------

   ! potenital melt
   posTemp = smax(TEMP_Z - MPARAM%MBASE, 0._sp)   ! smoothed max(TEMP_Z - MPARAM%MBASE, 0)
   potMelt = MF*posTemp   !  mm day-1
  
   ! melt capped by availability of snow
   meltCap  = snow + MBANDS(ISNW)%SWE / DT
   snowmelt = -smax(-potMelt, -meltCap)    ! smooth min

   MBANDS(ISNW)%SNOWMELT = snowmelt 
  
   ! ----- update SWE ---------------------------------------------------------------------
   
   MBANDS(ISNW)%DSWE_DT = MBANDS(ISNW)%SNOWACCMLTN - MBANDS(ISNW)%SNOWMELT
   MBANDS(ISNW)%SWE     = MBANDS(ISNW)%SWE + MBANDS(ISNW)%DSWE_DT*DT
   MBANDS(ISNW)%SWE     = smax(MBANDS(ISNW)%SWE, 0._sp) ! safety: clamp for small roundoff

   ! ----- calculate effective precip (rain + melt)  ---------------------------------------

   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * (rain + snowmelt)
  
  END DO  ! looping through elevation bands  
  
  end associate

  END SUBROUTINE UPDATE_SWE_DIFF

end module update_swe_DIFF_MODULE
