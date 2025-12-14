module Q_BASEFLOW_DIFF_module

  implicit none

  private
  public :: Q_BASEFLOW_DIFF

contains


  SUBROUTINE Q_BASEFLOW_DIFF(fuseStruct, want_dflux)
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
  logical(lgt), intent(in), optional     :: want_dflux  ! if we want flux derivatives
  ! derivatives
  logical(lgt)                           :: comp_dflux  ! flag to compute flux derivatives
  integer(i4b)                           :: iState      ! state index
  real(sp)                               :: phi         ! scaled water storage, phi=w/ws
  real(sp)                               :: dqb_dw      ! derivative in baseflow flux w.r.t. water store
  ! -------------------------------------------------------------------------------------------------
  ! associate variables with elements of data structure
  associate(&
   M_FLUX => fuseStruct%flux         , &  ! fluxes
   dfx_dS => fuseStruct%df_dS        , &  ! deriv in fluxes w.r.t. states
   TSTATE => fuseStruct%state1       , &  ! trial state variables (end of step)
   MPARAM => fuseStruct%param_adjust , &  ! adjustable model parameters
   DPARAM => fuseStruct%param_derive   &  ! derived model parameters
   ) ! (associate)

  ! check the need to compute flux derivatives
  comp_dflux = .false.; if(present(want_dflux)) comp_dflux = want_dflux

  ! ---------------------------------------------------------------------------------------
  SELECT CASE(SMODL%iARCH2)
   
   ! --------------------------------------------------------------------------------------
   CASE(iopt_tens2pll_2) ! tension reservoir plus two parallel tanks
    M_FLUX%QBASE_2A = MPARAM%QBRATE_2A * TSTATE%FREE_2A    ! qbrate_2a is a fraction (T-1)
    M_FLUX%QBASE_2B = MPARAM%QBRATE_2B * TSTATE%FREE_2B    ! qbrate_2b is a fraction (T-1)
    M_FLUX%QBASE_2  = M_FLUX%QBASE_2A + M_FLUX%QBASE_2B    ! total baseflow
    if(comp_dflux) stop "q_baseflow: derivative not implemented yet for iopt_tens2pll_2"

   ! --------------------------------------------------------------------------------------
   CASE(iopt_unlimfrc_2) ! baseflow resvr of unlimited size (0-HUGE), frac rate
    M_FLUX%QBASE_2  = MPARAM%QB_PRMS * TSTATE%WATR_2       ! qb_prms is a fraction (T-1)
    if(comp_dflux) stop "q_baseflow: derivative not implemented yet for iopt_unlimfrc_2"
   
   ! --------------------------------------------------------------------------------------
   CASE(iopt_unlimpow_2) ! baseflow resvr of unlimited size (0-HUGE), power recession

    associate(qbsat=>DPARAM%QBSAT, w=>TSTATE%WATR_2, ws=>MPARAM%MAXWATR_2, p=>MPARAM%QB_POWR)
 
    ! ----- compute flux ------------------------------------------------------------------
    phi             = w/ws
    M_FLUX%QBASE_2  = qbsat*phi**p

    ! ----- compute derivative ------------------------------------------------------------
    if(comp_dflux) dqb_dw = (qbsat*p/ws)*phi**(p - 1._sp)

    end associate
   
   ! --------------------------------------------------------------------------------------
   CASE(iopt_topmdexp_2) ! topmodel exponential reservoir (-HUGE to HUGE)
    M_FLUX%QBASE_2  = DPARAM%QBSAT * EXP( -(1. - TSTATE%WATR_2/MPARAM%MAXWATR_2) )
    if(comp_dflux) stop "q_baseflow: derivative not implemented yet for iopt_topmdexp_2"

   ! --------------------------------------------------------------------------------------
   CASE(iopt_fixedsiz_2) ! baseflow reservoir of fixed size
    M_FLUX%QBASE_2  = MPARAM%BASERTE * (TSTATE%WATR_2/MPARAM%MAXWATR_2)**MPARAM%QB_POWR
    if(comp_dflux) stop "q_baseflow: derivative not implemented yet for iopt_fixedsiz_2"
   
   ! --------------------------------------------------------------------------------------
   CASE DEFAULT
    print *, "SMODL%iARCH2 must be iopt_tens2pll_2, iopt_unlimfrc_2, iopt_unlimpow_2"
    print *, "  iopt_topmdexp_2, or iopt_fixedsiz_2"
    STOP
   ! --------------------------------------------------------------------------------------
  
  END SELECT
  ! ---------------------------------------------------------------------------------------
  
       ! populate derivative vector
  if(comp_dflux)then
    do iState=1,nState
      select case(cState(iState)%iSNAME)
       case (iopt_WATR_2); dfx_dS(iState)%QBASE_2 = dqb_dw ! exists if one state in the upper layer
      end select  ! no default needed
    end do ! looping through states
  endif

  end associate  ! end association with variables in the data structures
  END SUBROUTINE Q_BASEFLOW_DIFF

end module Q_BASEFLOW_DIFF_module
