module QSATEXCESS_DIFF_MODULE

  implicit none

  private
  public :: QSATEXCESS_DIFF

contains 

  SUBROUTINE QSATEXCESS_DIFF(fuseStruct, want_dflux)
  ! -------------------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to create a differentiable model, 12/25
  ! -------------------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Computes the saturated area and surface runoff
  ! -------------------------------------------------------------------------------------------------
  USE nrtype                                            ! variable types, etc.
  USE data_types, only: parent                          ! fuse parent data type
  USE model_defn                                        ! model definition structure
  USE model_defnames
  USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
  USE smoothers, only : smax,dsmax                      ! smoothed max function, derivative
  IMPLICIT NONE
  ! input-output
  type(parent), intent(inout)            :: fuseStruct  ! parent fuse data structure
  logical(lgt), intent(in), optional     :: want_dflux  ! if we want flux derivatives
  ! internal variables -- vic
  real(sp)                               :: u,xp        ! temporary variables
  real(sp)                               :: ds_dx       ! derivative of saturated area w.r.t. x
  real(sp)                               :: dx_du       ! derivative of smooth max(u,0) w.r.t. u
  real(sp)                               :: du_dw       ! derivative of u w.r.t. w
  real(sp)                               :: ds_dw       ! derivative of saturated area w.r.t. w
  ! internal variables -- topmodel
  REAL(SP)                               :: TI_SAT      ! topographic index where saturated
  REAL(SP)                               :: TI_LOG      ! critical value of topo index in log space
  REAL(SP)                               :: TI_OFF      ! offset in the Gamma distribution
  REAL(SP)                               :: TI_SHP      ! shape of the Gamma distribution
  REAL(SP)                               :: TI_CHI      ! CHI, see Sivapalan et al., 1987
  REAL(SP)                               :: TI_ARG      ! argument of the Gamma function
  REAL(SP)                               :: NO_ZERO=1.E-8  ! avoid divide by zero
  ! derivatives
  logical(lgt)                           :: comp_dflux  ! flag to compute flux derivatives
  integer(i4b)                           :: iState      ! state index  
  real(sp), parameter                    :: ms=1.e-4_sp ! smoothing in smax function 
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   dfx_dS => fuseStruct%df_dS        , &  ! deriv in fluxes w.r.t. states
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)
  ! -------------------------------------------------------------------------------------------------
 
  ! check the need to compute flux derivatives
  comp_dflux = .false.; if(present(want_dflux)) comp_dflux = want_dflux

  ! saturated area method
  SELECT CASE(SMODL%iQSURF)
  
   ! ------------------------------------------------------------------------------------------------
   ! ----- ARNO/Xzang/VIC parameterization (upper zone control) -------------------------------------
   ! ------------------------------------------------------------------------------------------------
   CASE(iopt_arno_x_vic)
  
    ! define variables
    associate(w=>TSTATE%WATR_1, wmax=>MPARAM%MAXWATR_1, b=>MPARAM%AXV_BEXP)
    
    ! ----- compute flux ----------------------------------------------------------------------------
    u  = 1._sp - w/wmax
    xp = smax(u, 0._sp, ms)   ! smooth version of max(u,0)
    M_FLUX%SATAREA = 1._sp - xp**b

    ! ----- compute derivatives ---------------------------------------------------------------------
    if(comp_dflux)then

      ! compute derivative w.r.t. saturated area
      ds_dx = -b*xp**(b - 1._sp)  ! derivative of saturated area w.r.t. xp
      dx_du = dsmax(u, 0._sp, ms) ! derivative of smooth max(u,0) w.r.t. u
      du_dw = -1._sp/wmax         ! derivative of u w.r.t. w
      ds_dw = du_dw*dx_du*ds_dx   ! derivative of saturated area w.r.t. w

      ! since WATR_1 is the sum of individual state variables (e.g., WATR_1=TENS_1+FREE_1) simply copy derivative
      do iState=1,nState
       select case(cState(iState)%iSNAME)
         case (iopt_TENS1A); dfx_dS(iState)%SATAREA = ds_dw ! exists if two tension tanks
         case (iopt_TENS1B); dfx_dS(iState)%SATAREA = ds_dw ! exists if two tension tanks
         case (iopt_TENS_1); dfx_dS(iState)%SATAREA = ds_dw ! exists if one tension tank
         case (iopt_FREE_1); dfx_dS(iState)%SATAREA = ds_dw ! exists if separate free storage
         case (iopt_WATR_1); dfx_dS(iState)%SATAREA = ds_dw ! exists if one state in the upper layer
       end select  ! no default needed
      end do ! looping through states

    endif  ! if want derivatives

    end associate

   ! ------------------------------------------------------------------------------------------------
   ! ----- PRMS variant (fraction of upper tension storage) ----------------------------------------- 
   ! ------------------------------------------------------------------------------------------------
   CASE(iopt_prms_varnt)
    
    ! ----- compute flux ----------------------------------------------------------------------------
    M_FLUX%SATAREA = MIN(TSTATE%TENS_1/DPARAM%MAXTENS_1, 1._sp) * MPARAM%SAREAMAX
  
    ! ----- compute derivatives ---------------------------------------------------------------------
    if(comp_dflux) stop "qsatexcess: derivatives for iopt_prms_varnt not implemented yet"

   ! ------------------------------------------------------------------------------------------------
   ! ----- TOPMODEL parameterization (only valid for TOPMODEL qb) ----------------------------------- 
   ! ------------------------------------------------------------------------------------------------
   CASE(iopt_tmdl_param)
  
    ! ----- compute flux ----------------------------------------------------------------------------
  
    ! compute the minimum value of the topographic index where the basin is saturated
    ! (this is correct, as MPARAM%MAXWATR_2 is m*n -- units are meters**(1/n)
    TI_SAT = DPARAM%POWLAMB / (TSTATE%WATR_2/MPARAM%MAXWATR_2 + NO_ZERO)
    ! compute the saturated area
    IF (TI_SAT.GT.DPARAM%MAXPOW) THEN
     M_FLUX%SATAREA = 0.
    ELSE
     ! convert the topographic index to log space
     TI_LOG = LOG( TI_SAT**MPARAM%QB_POWR )   
     ! compute the saturated area (NOTE: critical value of the topographic index is in log space)
     TI_OFF = 3._sp           ! offset in the Gamma distribution (the "3rd" parameter)
     TI_SHP = MPARAM%TISHAPE  ! shape of the Gamma distribution (the "2nd" parameter)
     TI_CHI = (MPARAM%LOGLAMB - TI_OFF) / MPARAM%TISHAPE ! Chi -- loglamb is the first parameter (mean)
     TI_ARG = MAX(0._sp, TI_LOG - TI_OFF) / TI_CHI       ! argument to the incomplete Gamma function
     M_FLUX%SATAREA = 1._sp - GAMMP(TI_SHP, TI_ARG)      ! GAMMP is the incomplete Gamma function
    ENDIF
  
    ! ----- compute derivatives ---------------------------------------------------------------------
    if(comp_dflux) stop "qsatexcess: derivatives for iopt_tmdl_param not implemented yet"

   ! ------------------------------------------------------------------------------------------------
   ! ------------------------------------------------------------------------------------------------
   ! check processed surface runoff selection 
   CASE DEFAULT
    print *, "SMODL%iQSURF must be iopt_arno_x_vic, iopt_prms_varnt, or iopt_tmdl_param"
    STOP
  
  END SELECT  ! (different surface runoff options)

  ! ...and, compute surface runoff
  ! ------------------------------
  M_FLUX%QSURF = M_FLUX%EFF_PPT * M_FLUX%SATAREA

  end associate  ! end association with variables in the data structures
  END SUBROUTINE QSATEXCESS_DIFF

end module QSATEXCESS_DIFF_MODULE
