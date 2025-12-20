module QPERCOLATE_DIFF_module

  implicit none

  private
  public :: QPERCOLATE_DIFF

contains

  SUBROUTINE QPERCOLATE_DIFF(fuseStruct, want_dflux)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the percolation from the upper soil layer to the lower soil layer
  ! ---------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames                                    ! model definition names
  use smoothers, only : sfrac, dsfrac                   ! smoothed fraction, derivative
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  logical(lgt), intent(in), optional     :: want_dflux  ! if we want flux derivatives
  ! internal
  logical(lgt)                           :: comp_dflux  ! flag to compute flux derivatives
  integer(i4b)                           :: iState      ! state index
  real(sp)                               :: phi         ! smoothed fraction of free water
  real(sp)                               :: dphi_dx     ! derivative in smoothed fraction of free water
  real(sp)                               :: df_dpsi     ! derivative of flux w.r.t. fraction
  real(sp)                               :: dqperc_dx   ! derivative of percolation fux w.r.t. water state
  REAL(SP)                               :: LZ_PD       ! lower zone percolation demand
  real(sp), parameter                    :: ms=1.e-4_sp ! smoothing in sfrac(smax) function
  ! ---------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   dfx_dS => fuseStruct%df_dS        , &  ! deriv in fluxes w.r.t. states
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! ---------------------------------------------------------------------------------------

  ! check the need to compute flux derivatives
  comp_dflux = .false.; if(present(want_dflux)) comp_dflux = want_dflux

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iQPERC)

   ! --------------------------------------------------------------------------------------
   ! upper zone control
   ! --------------------------------------------------------------------------------------
   CASE(iopt_perc_w2sat, iopt_perc_f2sat)
  
     ! short-cuts
     associate(k=>MPARAM%PERCRTE, c=>MPARAM%PERCEXP)
 
      ! compute fractions
      select case(SMODL%iQPERC)
       case(iopt_perc_w2sat); phi = sfrac(TSTATE%WATR_1, MPARAM%MAXWATR_1, ms)
       case(iopt_perc_f2sat); phi = sfrac(TSTATE%FREE_1, DPARAM%MAXFREE_1, ms)
      end select ! no need for default since already in block
     
      ! ----- compute flux ----------------------------------------------------------------
      M_FLUX%QPERC_12 = k*phi**c

     ! ----- compute derivative ----------------------------------------------------------
     if(comp_dflux)then

      ! compute derivative in the fractions
      select case(SMODL%iQPERC)
       case(iopt_perc_w2sat); dphi_dx = dsfrac(TSTATE%WATR_1, MPARAM%MAXWATR_1, ms)
       case(iopt_perc_f2sat); dphi_dx = dsfrac(TSTATE%FREE_1, DPARAM%MAXFREE_1, ms)
      end select ! no need for default since already in block

      ! compute derivatives in the percolation flux
      df_dpsi   = k*c*phi**(c - 1._sp) ! derivative of flux w.r.t. fraction
      dqperc_dx = df_dpsi*dphi_dx

      ! populate derivative vector
      do iState=1,nState
       select case(cState(iState)%iSNAME)
        case (iopt_FREE_1); dfx_dS(iState)%QPERC_12 = dqperc_dx ! exists if separate free tank
        case (iopt_WATR_1); dfx_dS(iState)%QPERC_12 = dqperc_dx ! exists if one state in the upper layer
       end select  ! no default needed
      end do ! looping through states

     endif  ! if computing derivatives
  
     end associate

   ! --------------------------------------------------------------------------------------
   ! lower zone control
   ! --------------------------------------------------------------------------------------
   CASE(iopt_perc_lower) ! perc defined by moisture content in lower layer (SAC)
    
    ! ----- compute flux ----------------------------------------------------------------
    LZ_PD = 1._SP + MPARAM%SACPMLT*(1._SP - TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%SACPEXP
    M_FLUX%QPERC_12 = DPARAM%QBSAT*LZ_PD * (TSTATE%FREE_1/DPARAM%MAXFREE_1)

    ! ----- compute derivatives ---------------------------------------------------------------------
    if(comp_dflux) stop "qpercolate: derivatives for iopt_perc_lower not implemented yet"

   CASE DEFAULT; stop "qpercolate: SMODL%iQPERC must be iopt_perc_f2sat, iopt_perc_w2sat, or iopt_perc_lower"
  END SELECT
  ! --------------------------------------------------------------------------------------

  end associate  ! end association with variables in the data structures
  END SUBROUTINE QPERCOLATE_DIFF

end module QPERCOLATE_DIFF_module
