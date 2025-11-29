module EVAP_LOWER_DIFF_MODULE

  implicit none

  private
  public :: EVAP_LOWER_DIFF

contains

  SUBROUTINE EVAP_LOWER_DIFF(fuseStruct)
  ! -------------------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! -------------------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes evaporation from the lower soil layer
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
  SELECT CASE(SMODL%iARCH2)  ! lower layer architecture
   CASE(iopt_tens2pll_2,iopt_fixedsiz_2)

    ! -------------------------------------------------------------------------------------
    SELECT CASE(SMODL%iARCH1)
     ! ------------------------------------------------------------------------------------
     CASE(iopt_tension1_1,iopt_onestate_1) ! lower-layer evap is valid
 
     ! ------------------------------------------------------------------------------------
     ! use different evaporation schemes for the lower layer
     ! -----------------------------------------------------
     SELECT CASE(SMODL%iESOIL)
      CASE(iopt_sequential)
       M_FLUX%EVAP_2 = (MFORCE%PET-M_FLUX%EVAP_1) * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
      CASE(iopt_rootweight)
       M_FLUX%EVAP_2 = MFORCE%PET * DPARAM%RTFRAC2 * (TSTATE%TENS_2/DPARAM%MAXTENS_2)
      CASE DEFAULT
       print *, "SMODL%iESOIL must be either iopt_sequential or iopt_rootweight"
     END SELECT  ! (evaporation schemes)
     
     ! ------------------------------------------------------------------------------------
     CASE(iopt_tension2_1)               ! lower-layer evap is zero
      M_FLUX%EVAP_2 = 0._sp
     
     ! ------------------------------------------------------------------------------------
     CASE DEFAULT
      print *, "SMODL%iARCH1 must be iopt_tension2_1, iopt_tension1_1, or iopt_onestate_1"
      STOP
     
    ! ------------------------------------------------------------------------------------
    END SELECT  ! (upper-layer architechure)

   ! --------------------------------------------------------------------------------------
   CASE(iopt_unlimfrc_2,iopt_unlimpow_2,iopt_topmdexp_2)
    M_FLUX%EVAP_2 = 0._sp
   
   ! --------------------------------------------------------------------------------------
   CASE DEFAULT
    print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
    print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
    STOP
  
  END SELECT
  ! ---------------------------------------------------------------------------------------
    
  end associate  ! end association with variables in the data structures
  END SUBROUTINE EVAP_LOWER_DIFF

end module EVAP_LOWER_DIFF_module
