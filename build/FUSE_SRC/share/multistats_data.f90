MODULE multistats

 USE nrtype
 USE multistats_types, only: SUMMARY

 implicit none
 private

 public :: MSTATS, MOD_IX, PCOUNT, FCOUNT
 
 TYPE(SUMMARY)                         :: MSTATS        ! (model summary statistics)
 INTEGER(I4B)                          :: MOD_IX = 1    ! (model index)
 INTEGER(I4B)                          :: PCOUNT        ! (number of parameter sets in model output files)
 INTEGER(I4B)                          :: FCOUNT        ! (number of model simulations)

END MODULE multistats
