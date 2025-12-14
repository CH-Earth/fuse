module MOD_DERIVS_DIFF_module

  USE nrtype
  USE data_types, only: parent, statev
  USE qsatexcess_diff_module, only: qsatexcess_diff
  USE evap_upper_diff_module, only: evap_upper_diff
  USE evap_lower_diff_module, only: evap_lower_diff
  USE qinterflow_diff_module, only: qinterflow_diff
  USE qpercolate_diff_module, only: qpercolate_diff
  USE q_baseflow_diff_module, only: q_baseflow_diff
  USE q_misscell_diff_module, only: q_misscell_diff
  USE mstate_rhs_diff_module, only: mstate_rhs_diff

  implicit none

  private
  public :: MOD_DERIVS_DIFF

contains

  SUBROUTINE MOD_DERIVS_DIFF(fuseStruct, g_x, J_g)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified to include snow model by Brian Henn, 6/2013
  ! Modified to include analytical derivatives by Martyn Clark, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! compute the time derivative (dx/dt) of all model states (x)
  ! --------------------------------------------------------------------------------------
  implicit none
  ! input
  type(parent) , intent(inout)            :: fuseStruct  ! parent fuse data structure
  ! output
  real(sp)     , intent(out)              :: g_x(:)      ! dx/dt=g(x)
  real(sp)     , intent(out)  , optional  :: J_g(:,:)    ! flux Jacobian matrix
  ! internal
  logical(lgt)                            :: comp_dflux  ! flag to compute flux derivatives
  ! --------------------------------------------------------------------------------------

  ! check if Jacobian is desired
  comp_dflux = present(J_g)

  ! compute fluxes
  call qsatexcess_diff(fuseStruct, comp_dflux)   ! compute the saturated area and surface runoff
  call evap_upper_diff(fuseStruct, comp_dflux)   ! compute evaporation from the upper layer
  call evap_lower_diff(fuseStruct, comp_dflux)   ! compute evaporation from the lower layer
  call qinterflow_diff(fuseStruct, comp_dflux)   ! compute interflow from free water in the upper layer
  call qpercolate_diff(fuseStruct, comp_dflux)   ! compute percolation from the upper to lower soil layers
  call q_baseflow_diff(fuseStruct, comp_dflux)   ! compute baseflow from the lower soil layer
  call q_misscell_diff(fuseStruct, comp_dflux)   ! compute miscellaneous fluxes (NOTE: need sat area, evap, and perc)

  ! compute the time derivative (dx/dt) of all model states (x)
  if(comp_dflux)then
   call mstate_rhs_diff(fuseStruct, g_x, J_g)
  else
   call mstate_rhs_diff(fuseStruct, g_x)
  endif

  END SUBROUTINE MOD_DERIVS_DIFF

end module MOD_DERIVS_DIFF_module
