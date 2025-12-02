module smoothers

  implicit none

  private
  public:: LOGISMOOTH
  public:: smoother

contains

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
