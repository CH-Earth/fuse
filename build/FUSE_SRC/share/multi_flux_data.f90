MODULE multi_flux

 USE nrtype
 
 USE multi_flux_types, only: FLUXES

 implicit none
 private

 public :: M_FLUX, FLUX_0, FLUX_1, FDFLUX, W_FLUX, W_FLUX_3d
 public :: CURRENT_DT

 TYPE(FLUXES)                          :: M_FLUX      ! model fluxes
 TYPE(FLUXES)                          :: FLUX_0      ! model fluxes at start of step
 TYPE(FLUXES)                          :: FLUX_1      ! model fluxes at end of step
 TYPE(FLUXES), DIMENSION(:), POINTER   :: FDFLUX=>NULL() ! finite difference fluxes
 TYPE(FLUXES)                          :: W_FLUX      ! weighted sum of model fluxes over a time step
 TYPE(FLUXES), dimension(:,:,:), allocatable  :: W_FLUX_3d   ! weighted sum of model fluxes over a time step for several time steps
 
 REAL(SP)                              :: CURRENT_DT  ! current time step (days)

END MODULE multi_flux
