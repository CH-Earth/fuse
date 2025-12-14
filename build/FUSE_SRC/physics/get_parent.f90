module get_parent_module
  use nrtype
  use data_types, only: parent
  USE model_defn, ONLY:NSTATE
  implicit none

contains

  subroutine get_parent(fuseStruct)
  use multiforce, only: mForce
  use multistate, only: mState
  use multi_flux, only: m_flux
  use multiparam, only: parMeta,mParam,dParam
  implicit none
  type(parent), intent(inout) :: fuseStruct
  integer(i4b)                :: iState

  ! populate parent fuse structures
  fuseStruct%force        = mForce
  fuseStruct%state0       = mState
  fuseStruct%state1       = mState
  fuseStruct%flux         = m_flux  ! initialized at zero
  fuseStruct%param_meta   = parMeta
  fuseStruct%param_adjust = mParam
  fuseStruct%param_derive = dParam

  ! initialize derivatives
  do iState=1,nState
   fuseStruct%df_dS(iState) = m_flux ! initialized at zero
  end do

  end subroutine get_parent


end module get_parent_module
