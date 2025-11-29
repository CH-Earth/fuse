module overshoot_module

  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn, only: CSTATE,NSTATE,SMODL             ! model definition structures
  USE model_defnames
  implicit none

  private
  public :: get_bounds 
  public :: fix_ovshoot

contains

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! Numerically-stable softplus with sharpness alpha
  pure real(sp) function softplus(x, alpha) result(y)
    implicit none
    real(sp), intent(in) :: x, alpha
    real(sp) :: ax
    ax = alpha * x
    if (ax > 0.0_sp) then
      y = (ax + log(1.0_sp + exp(-ax))) / alpha
    else
      y = log(1.0_sp + exp(ax)) / alpha
    end if
  end function softplus
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE fix_ovshoot(X_TRY, lower, upper)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Apply soft constraints to model state variables
  ! ---------------------------------------------------------------------------------------
  ! input/output
  REAL(SP), DIMENSION(:), INTENT(INOUT)  :: X_TRY       ! vector of model states
  real(sp), dimension(:), intent(in)     :: lower       ! lower bound
  real(sp), dimension(:), intent(in)     :: upper       ! upper bound
  ! internal
  integer(i4b)                           :: i           ! index of model state variable
  real(sp), parameter                    :: alpha=10_sp ! controls sharpness in smoothing

  ! apply soft constraint to model states
  do i=1,NSTATE
   x_try(i) = lower(i) + softplus(x_try(i)-lower(i), alpha) - softplus(x_try(i)-upper(i), alpha)
  end do  ! looping through model state variables

  end subroutine fix_ovshoot

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  SUBROUTINE get_bounds(fuseStruct, lower, upper)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified to return lower and upper bounds by Martyn Clark, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Identify lower and upper bounds for the vector of model states
  ! ---------------------------------------------------------------------------------------
  USE model_numerix                                     ! model numerix
  IMPLICIT NONE
  ! input/output
  type(parent), intent(in)               :: fuseStruct  ! parent fuse data structure
  real(sp), dimension(:), intent(out)    :: lower       ! lower bound for states
  real(sp), dimension(:), intent(out)    :: upper       ! upper bound for states
  ! internal
  REAL(SP)                               :: XMIN        ! very small number
  INTEGER(I4B)                           :: ISTT        ! loop through model states
  ! ---------------------------------------------------------------------------------------
  associate(MPARAM => fuseStruct%param_adjust, &        ! adjuustable model parameters
            DPARAM => fuseStruct%param_derive)          ! derived model parameters
  ! ---------------------------------------------------------------------------------------
  XMIN=FRACSTATE_MIN ! used to avoid zero derivatives
  ! ---------------------------------------------------------------------------------------
  ! loop through model states
  DO ISTT=1,NSTATE
   SELECT CASE(CSTATE(ISTT)%iSNAME)
    ! upper tanks
    CASE (iopt_TENS1A)
      lower(ISTT) = XMIN*DPARAM%MAXTENS_1A
      upper(ISTT) =      DPARAM%MAXTENS_1A
    CASE (iopt_TENS1B)
      lower(ISTT) = XMIN*DPARAM%MAXTENS_1B
      upper(ISTT) =      DPARAM%MAXTENS_1B
    CASE (iopt_TENS_1)
      lower(ISTT) = XMIN*DPARAM%MAXTENS_1
      upper(ISTT) =      DPARAM%MAXTENS_1
    CASE (iopt_FREE_1)
      lower(ISTT) = XMIN*DPARAM%MAXFREE_1
      upper(ISTT) =      DPARAM%MAXFREE_1
    CASE (iopt_WATR_1)
      lower(ISTT) = XMIN*MPARAM%MAXWATR_1
      upper(ISTT) =      MPARAM%MAXWATR_1
    ! lower tanks
    CASE (iopt_TENS_2)
      lower(ISTT) = XMIN*DPARAM%MAXTENS_2
      upper(ISTT) =      DPARAM%MAXTENS_2
    CASE (iopt_FREE2A)
      lower(ISTT) = XMIN*DPARAM%MAXFREE_2A
      upper(ISTT) =      DPARAM%MAXFREE_2A
    CASE (iopt_FREE2B)
      lower(ISTT) = XMIN*DPARAM%MAXFREE_2B
      upper(ISTT) =      DPARAM%MAXFREE_2B
    CASE (iopt_WATR_2)
      ! *** SET LOWER LIMITS ***
      IF (SMODL%iARCH2.NE.iopt_topmdexp_2) THEN
       ! enforce lower limit
       lower(ISTT) = XMIN*MPARAM%MAXWATR_2
      ELSE
       ! MPARAM%MAXWATR_2 is just a scaling parameter, but don't allow stupid values
       lower(ISTT) = -MPARAM%MAXWATR_2*10._sp
      ENDIF
      ! *** SET UPPER LIMITS ***
      IF (SMODL%iARCH2.EQ.iopt_tens2pll_2 .OR. SMODL%iARCH2.EQ.iopt_fixedsiz_2) THEN
       ! cannot exceed capacity
       upper(ISTT) = MPARAM%MAXWATR_2
      ELSE
       ! unlimited storage, but make sure the values are still sensible
       upper(ISTT) = MPARAM%MAXWATR_2*1000._sp
      ENDIF
   END SELECT
  END DO ! (loop through states)
  end associate  ! end association with variables in the data structures
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE get_bounds

END MODULE overshoot_module
