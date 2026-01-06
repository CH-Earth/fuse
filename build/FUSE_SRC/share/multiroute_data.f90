MODULE multiroute

 USE nrtype
 USE model_defn,ONLY:NTDH_MAX
 USE multiroute_types, only: RUNOFF

 implicit none
 private

 public :: FUTURE
 public :: AROUTE, AROUTE_3d
 public :: MROUTE

 REAL(SP), DIMENSION(NTDH_MAX)                :: FUTURE     ! runoff placed in future time steps

 TYPE(RUNOFF), DIMENSION(:), POINTER          :: AROUTE     ! runoff for all time steps
 TYPE(RUNOFF),dimension(:,:,:), allocatable   :: AROUTE_3d  ! runoff for all time steps on a grid
 TYPE(RUNOFF)                                 :: MROUTE     ! runoff for one time step

END MODULE multiroute
