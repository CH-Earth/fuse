module Q_BASEFLOW_DIFF_module

  implicit none

  private
  public :: Q_BASEFLOW_DIFF

contains


  SUBROUTINE Q_BASEFLOW_DIFF(fuseStruct)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the baseflow from the lower soil layer
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH2)
   ! --------------------------------------------------------------------------------------
   CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
    M_FLUX%QBASE_2A = MPARAM%QBRATE_2A * TSTATE%FREE_2A    ! qbrate_2a is a fraction (T-1)
    M_FLUX%QBASE_2B = MPARAM%QBRATE_2B * TSTATE%FREE_2B    ! qbrate_2b is a fraction (T-1)
    M_FLUX%QBASE_2  = M_FLUX%QBASE_2A + M_FLUX%QBASE_2B    ! total baseflow
   ! --------------------------------------------------------------------------------------
   CASE(iopt_unlimfrc_2) ! baseflow resvr of unlimited size (0-HUGE), frac rate
    M_FLUX%QBASE_2  = MPARAM%QB_PRMS * TSTATE%WATR_2       ! qb_prms is a fraction (T-1)
   ! --------------------------------------------------------------------------------------
   CASE(iopt_unlimpow_2) ! baseflow resvr of unlimited size (0-HUGE), power recession
    M_FLUX%QBASE_2  = DPARAM%QBSAT * (TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%QB_POWR
   ! --------------------------------------------------------------------------------------
   CASE(iopt_topmdexp_2) ! topmodel exponential reservoir (-HUGE to HUGE)
    M_FLUX%QBASE_2  = DPARAM%QBSAT * EXP( -(1. - TSTATE%WATR_2/MPARAM%MAXWATR_2) )
   ! --------------------------------------------------------------------------------------
   CASE(iopt_fixedsiz_2) ! baseflow reservoir of fixed size
    M_FLUX%QBASE_2  = MPARAM%BASERTE * (TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%QB_POWR
   ! --------------------------------------------------------------------------------------
   CASE DEFAULT
    print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
    print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
    STOP
   ! --------------------------------------------------------------------------------------
  END SELECT
  ! ---------------------------------------------------------------------------------------
  
  end associate  ! end association with variables in the data structures
  END SUBROUTINE Q_BASEFLOW_DIFF

end module Q_BASEFLOW_DIFF_module
