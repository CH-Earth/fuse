module MSTATE_RHS_DIFF_module

  implicit none

  private
  public :: MSTATE_RHS_DIFF

contains

  SUBROUTINE MSTATE_RHS_DIFF(fuseStruct)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes time derivatives of all states for all model combinations
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DX_DT  => fuseStruct%dx_dt          &  ! time derivative in states
   ) ! (associate)

  ! ---------------------------------------------------------------------------------------
  ! (1) COMPUTE TIME DERIVATIVES FOR STATES IN THE UPPER LAYER
  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH1)
   CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
    DX_DT%TENS_1A = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1A - M_FLUX%RCHR2EXCS
    DX_DT%TENS_1B = M_FLUX%RCHR2EXCS - M_FLUX%EVAP_1B - M_FLUX%TENS2FREE_1
    DX_DT%FREE_1  = M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1
   CASE(iopt_tension1_1) ! upper layer broken up into tension and free storage
    DX_DT%TENS_1  = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1 - M_FLUX%TENS2FREE_1
    DX_DT%FREE_1  = M_FLUX%TENS2FREE_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 - M_FLUX%OFLOW_1
   CASE(iopt_onestate_1) ! upper layer defined by a single state variable
    DX_DT%WATR_1  = M_FLUX%EFF_PPT - M_FLUX%QSURF - M_FLUX%EVAP_1 - M_FLUX%QPERC_12 - M_FLUX%QINTF_1 &
                    - M_FLUX%OFLOW_1
   CASE DEFAULT
    print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
    STOP
  END SELECT  ! (upper layer architechure)

  ! ---------------------------------------------------------------------------------------
  ! (2) COMPUTE TIME DERIVATIVES FOR STATES IN THE LOWER LAYER
  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH2)
   CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
    DX_DT%TENS_2  = M_FLUX%QPERC_12*(1._SP-MPARAM%PERCFRAC) - M_FLUX%EVAP_2 - M_FLUX%TENS2FREE_2
    DX_DT%FREE_2A = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2A &
                    - M_FLUX%OFLOW_2A
    DX_DT%FREE_2B = M_FLUX%QPERC_12*(MPARAM%PERCFRAC/2._SP) + M_FLUX%TENS2FREE_2/2._SP - M_FLUX%QBASE_2B &
                    - M_FLUX%OFLOW_2B
   CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2,iopt_fixedsiz_2) ! single state
    ! (NOTE: M_FLUX%OFLOW_2=0 for 'unlimfrc_2','unlimpow_2','topmdexp_2') 
    DX_DT%WATR_2  = M_FLUX%QPERC_12 - M_FLUX%EVAP_2 - M_FLUX%QBASE_2 - M_FLUX%OFLOW_2
   CASE DEFAULT
    print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
    print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
    STOP
  END SELECT
  ! ---------------------------------------------------------------------------------------

  end associate  ! end association with variables in the data structures
  END SUBROUTINE MSTATE_RHS_DIFF

end module MSTATE_RHS_DIFF_module
