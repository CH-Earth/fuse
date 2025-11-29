MODULE multistate
 USE nrtype
 use data_types, only: statev, m_time   ! <— import canonical types
 
 ! variable definitions
 type(statev),dimension(:,:),pointer   :: gState     ! (grid of model states)
 type(statev),dimension(:,:,:),pointer :: gState_3d  ! (grid of model states with a time dimension)
 TYPE(STATEV)                          :: ASTATE     ! (model states at the start of full timestep)
 TYPE(STATEV)                          :: FSTATE     ! (model states at start of sub-timestep)
 TYPE(STATEV)                          :: MSTATE     ! (model states at start/middle of sub-timestep)
 TYPE(STATEV)                          :: TSTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: BSTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: ESTATE     ! (temporary copy of model states)
 TYPE(STATEV)                          :: DSTATE     ! (default model states)
 TYPE(STATEV)                          :: DYDT_0     ! (derivative of model states at start of sub-step)
 TYPE(STATEV)                          :: DYDT_1     ! (derivative of model states at end of sub-step)
 TYPE(STATEV)                          :: DY_DT      ! (derivative of model states)
 TYPE(STATEV)                          :: DYDT_OLD   ! (derivative of model states for final solution)
 TYPE(M_TIME)                          :: HSTATE     ! (time interval to advance model states)

 ! NetCDF
 integer(i4b)                          :: ncid_out=-1              ! NetCDF output file ID

 ! initial store fraction (initialization)
 real(sp),parameter::fracState0=0.25_sp

END MODULE multistate
