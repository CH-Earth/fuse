MODULE multiparam

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------

 USE nrtype
 USE multiparam_types, only: PARATT ! included for legacy for routines that USE multiparam
 USE multiparam_types, only: PARADJ, PARDVD, PARINFO, PAR_ID

 implicit none
 private

 public :: PARATT, PARADJ, PARDVD, PARINFO, PAR_ID

 public :: MAXPAR, NUMPAR
 public :: APARAM, MPARAM, DPARAM
 public :: PARMETA, LPARAM
 public :: SOBOL_INDX

 INTEGER(I4B), PARAMETER               :: MAXPAR=50   ! maximum number of parameters for a single model
 INTEGER(I4B)                          :: NUMPAR      ! number of model parameters for current model
 
 TYPE(PARADJ), DIMENSION(:), POINTER   :: APARAM=>null()  ! all model parameter sets; DK/2008/10/21: explicit null
 TYPE(PARADJ)                          :: MPARAM      ! single model parameter set
 TYPE(PARDVD)                          :: DPARAM      ! derived model parameters
 
 TYPE(PARINFO)                         :: PARMETA     ! parameter metadata (all parameters)
 TYPE(PAR_ID), DIMENSION(MAXPAR)       :: LPARAM      ! list of model parameter names (need to modify to 16 for SCE)
 
 INTEGER(I4B)                          :: SOBOL_INDX  ! code to re-assemble Sobol parameters

END MODULE multiparam
