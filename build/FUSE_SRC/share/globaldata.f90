MODULE globaldata

 USE nrtype
 
 implicit none
 include "fuseversion.inc"

 ! time step
 REAL(SP), save          :: CURRENT_DT            ! current time step (days)

 ! missing values
 INTEGER(I4B),PARAMETER  :: NA_VALUE=-9999        ! integer designating missing values - TODO: retrieve from NetCDF file
 REAL(SP),PARAMETER      :: NA_VALUE_SP=-9999_sp  ! integer designating missing values - TODO: retrieve from NetCDF file

 ! NetCDF
 integer(i4b), save      :: ncid_out=-1           ! NetCDF output file ID

 ! initial store fraction (initialization)
 real(sp), parameter     :: fracState0=0.25_sp

 ! original code
 logical(lgt), save      :: isOriginal=.true.

 ! print flag
 logical(lgt), save      :: isPrint=.true.
 logical(lgt), save      :: isDebug=.false.

 ! snow parameters
 integer(i4b), parameter :: NPAR_SNOW=7
 integer(i4b), parameter :: iMBASE=1, iMFMAX=2, iMFMIN=3, iPXTEMP=4, iOPG=5, iLAPSE=6  ! indices in vectors
 integer(i4b), parameter :: iPERR=7   ! not a snow parameter, but used here

 ! number of fuse evaluations
 integer(i4b), save      :: nFUSE_eval    

end MODULE globaldata
