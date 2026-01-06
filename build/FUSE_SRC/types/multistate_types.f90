MODULE multistate_types

 USE nrtype

 implicit none
 private

 public :: STATEV, M_TIME

 ! --------------------------------------------------------------------------------------
 ! model state structure
 ! --------------------------------------------------------------------------------------
 TYPE STATEV
  ! snow layer
  REAL(SP)                             :: SWE_TOT    ! total storage as snow (mm)
  ! upper layer
  REAL(SP)                             :: WATR_1     ! total storage in layer1 (mm)
  REAL(SP)                             :: TENS_1     ! tension storage in layer1 (mm)
  REAL(SP)                             :: FREE_1     ! free storage in layer 1 (mm)
  REAL(SP)                             :: TENS_1A    ! storage in the recharge zone (mm)
  REAL(SP)                             :: TENS_1B    ! storage in the lower zone (mm)
  ! lower layer
  REAL(SP)                             :: WATR_2     ! total storage in layer2 (mm)
  REAL(SP)                             :: TENS_2     ! tension storage in layer2 (mm)
  REAL(SP)                             :: FREE_2     ! free storage in layer2 (mm)
  REAL(SP)                             :: FREE_2A    ! storage in the primary resvr (mm)
  REAL(SP)                             :: FREE_2B    ! storage in the secondary resvr (mm)
 END TYPE STATEV
 
 ! --------------------------------------------------------------------------------------
 ! model time structure
 ! --------------------------------------------------------------------------------------
 TYPE M_TIME
  REAL(SP)                             :: STEP       ! (time interval to advance model states)
 END TYPE M_TIME

END MODULE multistate_types
