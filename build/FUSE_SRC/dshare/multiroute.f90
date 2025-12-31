MODULE multiroute
 USE nrtype
 USE model_defn,ONLY:NTDH_MAX
 Use data_types,only:runoff
 REAL(SP)     , DIMENSION(NTDH_MAX)              :: FUTURE     ! runoff placed in future time steps
 TYPE(RUNOFF) , DIMENSION(:)     , allocatable   :: AROUTE     ! runoff for all time steps
 TYPE(RUNOFF)  ,dimension(:,:,:) , allocatable   :: AROUTE_3d  ! runoff for all time steps on a grid
 TYPE(RUNOFF)                                    :: MROUTE     ! runoff for one time step
END MODULE multiroute
