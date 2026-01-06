MODULE model_defn_types

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Martyn Clark to separate data tyoes from data store, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype
 
 implicit none
 private

 public :: DESC, UMODEL, SNAMES, FNAMES

 ! description of model component
 TYPE DESC
  CHARACTER(LEN=16)                    :: MCOMPONENT      ! description of model component
 END TYPE DESC
 
 ! structure that holds (x) unique combinations
 TYPE UMODEL
  INTEGER(I4B)                         :: MODIX           ! model index
  CHARACTER(LEN=256)                   :: MNAME           ! model name
  INTEGER(I4B)                         :: iRFERR
  INTEGER(I4B)                         :: iARCH1
  INTEGER(I4B)                         :: iARCH2
  INTEGER(I4B)                         :: iQSURF
  INTEGER(I4B)                         :: iQPERC
  INTEGER(I4B)                         :: iESOIL
  INTEGER(I4B)                         :: iQINTF
  INTEGER(I4B)                         :: iQ_TDH
  INTEGER(I4B)                         :: iSNOWM           ! snow
 END TYPE UMODEL

 ! structure to hold model state names
 TYPE SNAMES
  INTEGER(I4B)                         :: iSNAME          ! integer value of state name
 END TYPE SNAMES
 
 ! structure to hold model flux names
 TYPE FNAMES
  CHARACTER(LEN=16)                    :: FNAME           ! state name
 END TYPE FNAMES

END MODULE model_defn_types
