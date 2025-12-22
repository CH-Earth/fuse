! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark
! Modified by Brian Henn to include snow model, 6/2013
! Modified by Martyn Clark to separate derived types from shard data, 12/2025
! ---------------------------------------------------------------------------------------
MODULE multiparam
 USE nrtype
 USE data_types,ONLY:par_id,parinfo,paradj,pardvd
 ! --------------------------------------------------------------------------------------
 INTEGER(I4B), PARAMETER               :: MAXPAR=50   ! maximum number of parameters for a single model
 TYPE(PARADJ), DIMENSION(:), POINTER   :: APARAM=>null()  ! all model parameter sets; DK/2008/10/21: explicit null
 TYPE(PARADJ)                          :: MPARAM      ! single model parameter set
 TYPE(PARDVD)                          :: DPARAM      ! derived model parameters
 TYPE(PARINFO)                         :: PARMETA     ! parameter metadata (all parameters)
 TYPE(PAR_ID), DIMENSION(MAXPAR)       :: LPARAM      ! list of model parameter names (need to modify to 16 for SCE)
 INTEGER(I4B)                          :: NUMPAR      ! number of model parameters for current model
 INTEGER(I4B)                          :: SOBOL_INDX  ! code to re-assemble Sobol parameters
 integer(i4b)                          :: MAXN        ! maximum number of trials before optimization is terminated
 integer(i4b)                          :: KSTOP       ! number of shuffling loops the value must change by PCENTO
 REAL(MSP)                             :: PCENTO      ! the percentage
 ! --------------------------------------------------------------------------------------
END MODULE multiparam
