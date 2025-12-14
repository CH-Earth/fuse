module EVAP_UPPER_DIFF_module

  implicit none

  private
  public :: EVAP_UPPER_DIFF

contains

  SUBROUTINE EVAP_UPPER_DIFF(fuseStruct, want_dflux)
  ! -------------------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! -------------------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes evaporation from the upper soil layer
  ! -------------------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames                                    ! model definition names
  use smoothers, only : sfrac, dsfrac                   ! smoothed fraction, derivative
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  logical(lgt), intent(in), optional     :: want_dflux  ! if we want flux derivatives
  ! local variables
  logical(lgt)                           :: comp_dflux  ! flag to compute flux derivatives
  integer(i4b)                           :: iState      ! state index
  real(sp)                               :: phi         ! smoothed fraction of total tension storage (0,1]
  real(sp)                               :: phi_1a      ! smoothed fraction of primary tension storage (0,1]
  real(sp)                               :: phi_1b      ! smoothed fraction of secondary tension storage (0,1]
  real(sp)                               :: maxRate     ! maximum forcing
  real(sp)                               :: maxRate_1a  ! maximum forcing for the primary tension tank
  real(sp)                               :: maxRate_1b  ! maximum forcing for the secondary tension tank
  real(sp)                               :: dphi_dx     ! derivative in fraction w.r.t. storage
  real(sp)                               :: devap_dx    ! derivative in evaporation w.r.t. storage
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   MFORCE => fuseStruct%force        , &  ! model forcing data
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   dfx_dS => fuseStruct%df_dS        , &  ! deriv in fluxes w.r.t. states
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! -------------------------------------------------------------------------------------------------

  ! check the need to compute flux derivatives
  comp_dflux = .false.; if(present(want_dflux)) comp_dflux = want_dflux

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH1)  ! upper layer architecture

   ! --------------------------------------------------------------------------------------
   CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
   ! --------------------------------------------------------------------------------------
 
    ! calculate the smoothed fraction of tension storage (NOTE: use WATR_1)
    phi_1a = sfrac(TSTATE%TENS_1A, DPARAM%MAXTENS_1A)
    phi_1b = sfrac(TSTATE%TENS_1B, DPARAM%MAXTENS_1B)
    
    ! calculate the maximum evap rate for the storage
    SELECT CASE(SMODL%iESOIL)
     CASE(iopt_sequential)
      maxrate_1a = MFORCE%PET
      maxrate_1b = MFORCE%PET - MFORCE%PET*phi_1a
     CASE(iopt_rootweight)
      maxrate_1a = MFORCE%PET * MPARAM%RTFRAC1
      maxrate_1b = MFORCE%PET * DPARAM%RTFRAC2
     CASE DEFAULT; stop "evap_upper: SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
    END SELECT

    ! ----- compute flux ----------------------------------------------------------------  
    M_FLUX%EVAP_1A = maxrate_1a*phi_1a
    M_FLUX%EVAP_1B = maxrate_1b*phi_1b
    M_FLUX%EVAP_1  = M_FLUX%EVAP_1A + M_FLUX%EVAP_1B

    ! ----- compute derivatives ---------------------------------------------------------------------
    if(comp_dflux) stop "evap_upper: derivatives for iopt_tension2_1 not implemented yet"

   ! --------------------------------------------------------------------------------------
   CASE(iopt_tension1_1,iopt_onestate_1)   ! single tension store or single state
   ! --------------------------------------------------------------------------------------
 
    ! zero fluxes not used
    M_FLUX%EVAP_1A = 0._sp
    M_FLUX%EVAP_1B = 0._sp

    select case(SMODL%iARCH1)
     case(iopt_tension1_1); phi = sfrac(TSTATE%TENS_1, DPARAM%MAXTENS_1)
     case(iopt_onestate_1); phi = sfrac(TSTATE%WATR_1, DPARAM%MAXTENS_1) ! NOTE: use WATR_1
    end select ! no need for default because checked above

    ! calculate the maximum evap rate for the upper layer
    SELECT CASE(SMODL%iESOIL)
     CASE(iopt_sequential); maxRate = MFORCE%PET
     CASE(iopt_rootweight); maxRate = MFORCE%PET*MPARAM%RTFRAC1
     CASE DEFAULT; stop "evap_upper: SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
    END SELECT  ! (evaporation schemes)

    ! ----- compute flux ----------------------------------------------------------------
    M_FLUX%EVAP_1  = maxRate*phi 

    ! ----- compute derivatives ---------------------------------------------------------
    if(comp_dflux)then
     
     ! calculate the derivative in the smoothed fraction of tension storage
     select case(SMODL%iARCH1)
      case(iopt_tension1_1); dphi_dx = dsfrac(TSTATE%TENS_1, DPARAM%MAXTENS_1)
      case(iopt_onestate_1); dphi_dx = dsfrac(TSTATE%WATR_1, DPARAM%MAXTENS_1) ! NOTE: use WATR_1
     end select ! no need for default because checked above

     ! calculate the derivative in the maximum rate
     devap_dx = maxRate*dphi_dx
    
     ! populate derivative vector 
     do iState=1,nState
      select case(cState(iState)%iSNAME)
       case (iopt_TENS_1); dfx_dS(iState)%EVAP_1 = devap_dx ! exists if one tension tank
       case (iopt_WATR_1); dfx_dS(iState)%EVAP_1 = devap_dx ! exists if one state in the upper layer
      end select  ! no default needed
     end do ! looping through states

    endif  ! if computing derivatives

   CASE DEFAULT; stop "evap_upper: SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
  END SELECT  ! (upper-layer architecture)


  end associate  ! end association with variables in the data structures
  END SUBROUTINE EVAP_UPPER_DIFF

end module EVAP_UPPER_DIFF_module
