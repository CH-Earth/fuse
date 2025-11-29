MODULE globaldata

 USE nrtype

 ! time step
 REAL(SP), save         :: CURRENT_DT            ! current time step (days)

 ! missing values
 INTEGER(I4B),PARAMETER :: NA_VALUE=-9999        ! integer designating missing values - TODO: retrieve from NetCDF file
 REAL(SP),PARAMETER     :: NA_VALUE_SP=-9999_sp  ! integer designating missing values - TODO: retrieve from NetCDF file

 ! NetCDF
 integer(i4b), save     :: ncid_out=-1           ! NetCDF output file ID

 ! initial store fraction (initialization)
 real(sp), parameter    :: fracState0=0.25_sp

 ! print flag
 logical(lgt)           :: isPrint=.true.

end MODULE globaldata
