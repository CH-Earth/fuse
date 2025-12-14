module smoothers

  implicit none

  private
  public:: LOGISMOOTH
  public:: smoother
  public:: smax,dsmax
  public:: sfrac,dsfrac

contains

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  PURE FUNCTION sfrac(x,xmax) result(xf)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Use smoothed min function to compute smooth fraction
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  implicit none
  real(sp),   intent(in)                 :: x           ! x value
  real(sp),   intent(in)                 :: xmax        ! maximum value
  real(sp)                               :: xp          ! smooth min(x,xmax)
  real(sp)                               :: xf          ! smooth fraction x/xmax
  xp = xmax - smax(xmax - x, 0._sp)   ! smooth version of min(x, xmax)
  xf = max(0._sp, xp) / xmax          ! use max(0._sp, xp) to account for small neg values at zero
  end function sfrac

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  PURE FUNCTION dsfrac(x,xmax) result(dxf_dx)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Get derivative of the smooth fraction
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  implicit none
  real(sp),   intent(in)                 :: x           ! x value
  real(sp),   intent(in)                 :: xmax        ! maximum value
  real(sp)                               :: dxp_dx      ! derivative of the max smoother
  real(sp)                               :: dxf_dx      ! derivative of the smoothed fraction
  ! NOTE: ignore the hard clamp at zero (very small differences and not worth the extra expense)
  dxp_dx = dsmax(xmax - x, 0._sp)  ! note signs cancel out
  dxf_dx = dxp_dx / xmax
  end function dsfrac

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  PURE FUNCTION smax(x,xmin) result(xp)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Compute smoothed max function following Kavetski and Kuczera (2007)
  !
  ! Kavetski, D. and Kuczera, G., 2007. Model smoothing strategies to remove microscale
  ! discontinuities and spurious secondary optima in objective functions in hydrological
  ! calibration. Water Resources Research, 43(3).
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  implicit none
  real(sp),   intent(in)                 :: x           ! x value
  real(sp),   intent(in)                 :: xmin        ! minimum value
  real(sp),   parameter                  :: ms=1.e-4_sp ! smoothing parameter
  real(sp)                               :: srt         ! sqrt(x*x + ms)
  real(sp)                               :: xp          ! smooth max(x,xmin)
  srt = sqrt((x-xmin)**2 + ms)
  xp  = 0.5_sp*(x + xmin + srt)              ! smooth max(x,xmin)
  end function smax

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  PURE FUNCTION dsmax(x,xmin) result(dxp)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Compute derivative of smoothed max function of Kavetski and Kuczera (2007)
  !
  ! Kavetski, D. and Kuczera, G., 2007. Model smoothing strategies to remove microscale
  ! discontinuities and spurious secondary optima in objective functions in hydrological
  ! calibration. Water Resources Research, 43(3).
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  implicit none
  real(sp),   intent(in)                 :: x           ! x value
  real(sp),   intent(in)                 :: xmin        ! minimum value
  real(sp),   parameter                  :: ms=1.e-4_sp ! smoothing parameter
  real(sp)                               :: u           ! x-xmin
  real(sp)                               :: srt         ! sqrt(x*x + ms)
  real(sp)                               :: dxp         ! derivative of smooth max(x,xmin)
  u   = x-xmin
  srt = sqrt(u*u + ms)
  dxp = 0.5_sp*(1._sp + u/srt)              ! derivative of smooth max(x,xmin)
  end function dsmax

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  
  PURE FUNCTION smoother(STATE,STATE_MAX,PSMOOTH) result(w_func)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Provides the option of different smoothers
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN)                   :: STATE       ! model state
  REAL(SP), INTENT(IN)                   :: STATE_MAX   ! maximum model state
  REAL(SP), INTENT(IN)                   :: PSMOOTH     ! smoothing parameter (fraction of state)
  real(sp)                               :: w_func      ! smoothed threshold
  real(sp)                               :: delta       ! scale factor

  ! logistic smoothing (original)
  w_func = LOGISMOOTH(STATE,STATE_MAX,PSMOOTH)

  ! qintic smoother (plays better with Newton)
  !delta  = MAX(PSMOOTH*STATE_MAX, 1.0e-6_SP*STATE_MAX)
  !w_func = SMOOTHSTEP5_W(STATE,STATE_MAX,delta)

  end function smoother

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  PURE FUNCTION LOGISMOOTH(STATE,STATE_MAX,PSMOOTH)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Uses a logistic function to smooth the threshold at the top of a bucket
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN)                   :: STATE       ! model state
  REAL(SP), INTENT(IN)                   :: STATE_MAX   ! maximum model state
  REAL(SP), INTENT(IN)                   :: PSMOOTH     ! smoothing parameter (fraction of state)
  real(sp)                               :: arg         ! clamp argument
  REAL(SP)                               :: ASMOOTH     ! actual smoothing
  REAL(SP)                               :: LOGISMOOTH  ! FUNCTION name
  ! ---------------------------------------------------------------------------------------
  ASMOOTH = PSMOOTH*STATE_MAX                           ! actual smoothing
  arg     = -(STATE - (STATE_MAX - 5*ASMOOTH))/ASMOOTH  ! argument
  !arg     = max(min(arg, 50._SP), -50._SP)              ! clamp
  LOGISMOOTH = 1._SP / ( 1._SP + EXP(arg) )
  ! ---------------------------------------------------------------------------------------
  END FUNCTION LOGISMOOTH

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  
  PURE FUNCTION SMOOTHSTEP5_W(STATE, STATE_MAX, DELTA) RESULT(W)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Uses a qintic function to smooth the threshold at the top of a bucket
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: STATE, STATE_MAX, DELTA
  REAL(SP) :: W, x

  x = (STATE - (STATE_MAX - DELTA)) / DELTA
  IF (x <= 0._SP) THEN
     W = 0._SP
  ELSEIF (x >= 1._SP) THEN
     W = 1._SP
  ELSE
     W = x*x*x*(10._SP + x*(-15._SP + 6._SP*x))
  END IF
  END FUNCTION

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  
  PURE FUNCTION SMOOTHSTEP5_DWDS(STATE, STATE_MAX, DELTA) RESULT(DWDS)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Compute the derivative of the qintic function
  ! ---------------------------------------------------------------------------------------
  USE nrtype
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: STATE, STATE_MAX, DELTA
  REAL(SP) :: DWDS, x

  IF (DELTA <= 0._SP) THEN
     DWDS = 0._SP
     RETURN
  END IF

  x = (STATE - (STATE_MAX - DELTA)) / DELTA
  IF (x <= 0._SP .OR. x >= 1._SP) THEN
     DWDS = 0._SP
  ELSE
     DWDS = (30._SP * x*x * (1._SP - x)*(1._SP - x)) / DELTA
  END IF
  END FUNCTION

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

end module smoothers
