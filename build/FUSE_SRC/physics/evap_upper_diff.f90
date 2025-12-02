module EVAP_UPPER_DIFF_module

  implicit none

  private
  public :: EVAP_UPPER_DIFF

contains

  SUBROUTINE EVAP_UPPER_DIFF(fuseStruct)
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
  USE model_defnames
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   MFORCE => fuseStruct%force        , &  ! model forcing data
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! -------------------------------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH1)  ! upper layer architecture

   ! --------------------------------------------------------------------------------------
   CASE(iopt_tension2_1) ! tension storage sub-divided into recharge and excess
   ! --------------------------------------------------------------------------------------
 
    ! use different evaporation schemes for the upper layer
    ! -----------------------------------------------------
    SELECT CASE(SMODL%iESOIL)
     CASE(iopt_sequential)
      M_FLUX%EVAP_1A = MFORCE%PET * TSTATE%TENS_1A/DPARAM%MAXTENS_1A
      M_FLUX%EVAP_1B = (MFORCE%PET - M_FLUX%EVAP_1A) * TSTATE%TENS_1B/DPARAM%MAXTENS_1B
      M_FLUX%EVAP_1  = M_FLUX%EVAP_1A + M_FLUX%EVAP_1B
     CASE(iopt_rootweight)
      M_FLUX%EVAP_1A = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1A/DPARAM%MAXTENS_1A
      M_FLUX%EVAP_1B = MFORCE%PET * DPARAM%RTFRAC2 * TSTATE%TENS_1B/DPARAM%MAXTENS_1B
      M_FLUX%EVAP_1  = M_FLUX%EVAP_1A + M_FLUX%EVAP_1B
     CASE DEFAULT
      print *, "SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
      STOP
    END SELECT
   ! --------------------------------------------------------------------------------------
   CASE(iopt_tension1_1,iopt_onestate_1)   ! single tension store or single state
   ! --------------------------------------------------------------------------------------
 
    ! use different evaporation schemes for the upper layer
    ! -----------------------------------------------------
    SELECT CASE(SMODL%iESOIL)
     CASE(iopt_sequential)
      M_FLUX%EVAP_1A = 0._sp
      M_FLUX%EVAP_1B = 0._sp
      M_FLUX%EVAP_1  = MFORCE%PET * TSTATE%TENS_1/DPARAM%MAXTENS_1
     CASE(iopt_rootweight)
      M_FLUX%EVAP_1A = 0._sp
      M_FLUX%EVAP_1B = 0._sp
      M_FLUX%EVAP_1  = MFORCE%PET * MPARAM%RTFRAC1 * TSTATE%TENS_1/DPARAM%MAXTENS_1
     CASE DEFAULT
      print *, "SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
    END SELECT  ! (evaporation schemes)
 
   ! --------------------------------------------------------------------------------------
   CASE DEFAULT
    print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
    STOP
   ! --------------------------------------------------------------------------------------
  
  END SELECT  ! (upper-layer architechure)

  end associate  ! end association with variables in the data structures
  END SUBROUTINE EVAP_UPPER_DIFF

end module EVAP_UPPER_DIFF_module
