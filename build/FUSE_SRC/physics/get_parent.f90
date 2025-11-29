module get_parent_module
  use data_types, only: parent
  implicit none

contains

  subroutine get_parent(fuseStruct)
  use multiforce, only: mForce
  use multistate, only: mState
  use multi_flux, only: m_flux
  use multiparam, only: parMeta,mParam,dParam
  implicit none
  type(parent), intent(out)    :: fuseStruct
  ! populate parent fuse structures
  fuseStruct%force        = mForce
  fuseStruct%state0       = mState
  fuseStruct%state1       = mState
  fuseStruct%flux         = m_flux
  fuseStruct%param_meta   = parMeta
  fuseStruct%param_adjust = mParam
  fuseStruct%param_derive = dParam

  end subroutine get_parent


end module get_parent_module
