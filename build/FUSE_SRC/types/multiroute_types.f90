MODULE multiroute_types

 USE nrtype

 implicit none
 private

 public :: RUNOFF
 
 TYPE RUNOFF
  REAL(SP)                             :: Q_INSTNT   ! instantaneous runoff
  REAL(SP)                             :: Q_ROUTED   ! routed runoff
  REAL(SP)                             :: Q_ACCURATE ! "accurate" runoff estimate (mm day-1)
 END TYPE RUNOFF
 
END MODULE multiroute_types
