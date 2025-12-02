module Q_MISSCELL_DIFF_module

  implicit none

  private
  public :: Q_MISSCELL_DIFF

contains

  SUBROUTINE Q_MISSCELL_DIFF(fuseStruct)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007 
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes miscellaneous fluxes:
  !   RCHR2EXCS   = flow from recharge to excess (mm day-1)
  !   TENS2FREE_1 = flow from tension storage to free storage in the upper layer (mm day-1)
  !   TENS2FREE_2 = flow from tension storage to free storage in the lower layer (mm day-1)
  !   OFLOW_1     = overflow from the upper soil layer (mm day-1)
  !   OFLOW_2     = overflow from the lower soil layer (mm day-1)
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames
  USE smoothers, only: smoother                         ! smoothing function
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  ! internal
  REAL(SP), PARAMETER                    :: PSMOOTH=0.05_SP ! smoothing parameter
  REAL(SP)                               :: W_FUNC      ! result from smoother
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! ---------------------------------------------------------------------------------------
 
  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH1)
   CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
    ! compute flow from recharge to excess (mm s-1)
    W_FUNC = SMOOTHER(TSTATE%TENS_1A,DPARAM%MAXTENS_1A,PSMOOTH)
    M_FLUX%RCHR2EXCS   = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
    ! compute flow from tension storage to free storage (mm s-1)
    W_FUNC = SMOOTHER(TSTATE%TENS_1B,DPARAM%MAXTENS_1B,PSMOOTH)
    M_FLUX%TENS2FREE_1 = W_FUNC * M_FLUX%RCHR2EXCS
    ! compute over-flow of free water
    W_FUNC = SMOOTHER(TSTATE%FREE_1,DPARAM%MAXFREE_1,PSMOOTH)
    M_FLUX%OFLOW_1     = W_FUNC * M_FLUX%TENS2FREE_1
   CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
    ! no separate recharge zone (flux should never be used)
    M_FLUX%RCHR2EXCS   = 0._SP
    ! compute flow from tension storage to free storage (mm s-1)
    W_FUNC = SMOOTHER(TSTATE%TENS_1,DPARAM%MAXTENS_1,PSMOOTH)
    M_FLUX%TENS2FREE_1 = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
    ! compute over-flow of free water
    W_FUNC = SMOOTHER(TSTATE%FREE_1,DPARAM%MAXFREE_1,PSMOOTH)
    M_FLUX%OFLOW_1     = W_FUNC * M_FLUX%TENS2FREE_1
   CASE(iopt_onestate_1) ! upper layer defined by a single state variable
    ! no tension stores
    M_FLUX%RCHR2EXCS   = 0._SP
    M_FLUX%TENS2FREE_1 = 0._SP
    ! compute over-flow of free water
    W_FUNC = SMOOTHER(TSTATE%WATR_1,MPARAM%MAXWATR_1,PSMOOTH)
    M_FLUX%OFLOW_1     = W_FUNC * (M_FLUX%EFF_PPT - M_FLUX%QSURF)
   CASE DEFAULT
    print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
    STOP
  END SELECT

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH2)
   CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
    ! compute flow from tension storage to free storage (mm s-1)
    W_FUNC = SMOOTHER(TSTATE%TENS_2,DPARAM%MAXTENS_2,PSMOOTH)
    M_FLUX%TENS2FREE_2 = W_FUNC * M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC)
    ! compute over-flow of free water in the primary reservoir
    W_FUNC = SMOOTHER(TSTATE%FREE_2A,DPARAM%MAXFREE_2A,PSMOOTH)
    M_FLUX%OFLOW_2A    = W_FUNC * (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP)
    ! compute over-flow of free water in the secondary reservoir
    W_FUNC = SMOOTHER(TSTATE%FREE_2B,DPARAM%MAXFREE_2B,PSMOOTH)
    M_FLUX%OFLOW_2B    = W_FUNC * (M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP)
    ! compute total overflow
    M_FLUX%OFLOW_2     = M_FLUX%OFLOW_2A + M_FLUX%OFLOW_2B
   CASE(iopt_fixedsiz_2)
    ! no tension store
    M_FLUX%TENS2FREE_2 = 0._SP
    M_FLUX%OFLOW_2A    = 0._SP
    M_FLUX%OFLOW_2B    = 0._SP
    ! compute over-flow of free water
    W_FUNC = SMOOTHER(TSTATE%WATR_2,MPARAM%MAXWATR_2,PSMOOTH)
    M_FLUX%OFLOW_2     = W_FUNC * M_FLUX%QPERC_12
   CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2) ! unlimited size
    M_FLUX%TENS2FREE_2 = 0._SP
    M_FLUX%OFLOW_2     = 0._SP
    M_FLUX%OFLOW_2A    = 0._SP
    M_FLUX%OFLOW_2B    = 0._SP
   CASE DEFAULT
    print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
    print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
    STOP
  END SELECT
  
  end associate  ! end association with variables in the data structures
  END SUBROUTINE Q_MISSCELL_DIFF

end module Q_MISSCELL_DIFF_module
