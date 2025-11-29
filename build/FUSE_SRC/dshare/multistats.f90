MODULE multistats
 USE nrtype
 Use data_types, only: summary
 ! final data structures
 TYPE(SUMMARY)                         :: MSTATS        ! (model summary statistics)
 INTEGER(I4B)                          :: MOD_IX=1      ! (model index)
 INTEGER(I4B)                          :: PCOUNT        ! (number of parameter sets in model output files)
 INTEGER(I4B)                          :: FCOUNT        ! (number of model simulations)
END MODULE multistats
