module QINTERFLOW_DIFF_module

  implicit none

  private
  public :: QINTERFLOW_DIFF

contains

  SUBROUTINE QINTERFLOW_DIFF(fuseStruct)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the interflow from free water in the upper soil layer
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
  SELECT CASE(SMODL%iQINTF)
   CASE(iopt_intflwsome) ! interflow
    M_FLUX%QINTF_1 = MPARAM%IFLWRTE * (TSTATE%FREE_1/DPARAM%MAXFREE_1)
   CASE(iopt_intflwnone) ! no interflow
    M_FLUX%QINTF_1 = 0.
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iQINTF must be either iopt_intflwsome or iopt_intflwnone"
    STOP
  END SELECT
  ! ---------------------------------------------------------------------------------------

  end associate  ! end association with variables in the data structures
  END SUBROUTINE QINTERFLOW_DIFF

end module QINTERFLOW_DIFF_module
