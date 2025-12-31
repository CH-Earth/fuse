MODULE multi_flux
 USE nrtype
 use data_types, only: fluxes
 TYPE(FLUXES)                                 :: M_FLUX      ! model fluxes
 TYPE(FLUXES)                                 :: FLUX_0      ! model fluxes at start of step
 TYPE(FLUXES)                                 :: FLUX_1      ! model fluxes at end of step
 TYPE(FLUXES)                                 :: W_FLUX      ! weighted sum of model fluxes over a time step
 TYPE(FLUXES), dimension(:,:,:), allocatable  :: W_FLUX_3d   ! weighted sum of model fluxes over a time step for several time steps
 TYPE(FLUXES), dimension(:)    , pointer      :: FDFLUX => null() ! finite difference fluxes
END MODULE multi_flux
