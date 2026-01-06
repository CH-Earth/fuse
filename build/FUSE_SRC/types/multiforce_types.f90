MODULE multiforce_types

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Nans Addor to enable distributed modeling, 9/2016
 ! Modified by Cyril Thébault to allow different metrics as objective function, 2024
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype
 
 implicit none
 private

 public :: TDATA, VDATA, ADATA, FDATA

 ! the time data structure (will have no spatial dimension)
 TYPE TDATA
    INTEGER(I4B)                         :: IY         ! year
    INTEGER(I4B)                         :: IM         ! month
    INTEGER(I4B)                         :: ID         ! day
    INTEGER(I4B)                         :: IH         ! hour
    INTEGER(I4B)                         :: IMIN       ! minute
    REAL(SP)                             :: DSEC       ! second
    REAL(SP)                             :: DTIME      ! time in seconds since year dot
 ENDTYPE TDATA
 
 ! the response structure (will not have a spatial dimension)
 TYPE VDATA
    REAL(SP)                             :: OBSQ       ! observed runoff (mm day-1)
 END TYPE VDATA
 
 ! ancillary forcing variables used to compute ET (will have a spatial dimension)
 TYPE ADATA
    REAL(SP)                             :: AIRTEMP    ! air temperature (K)
    REAL(SP)                             :: SPECHUM    ! specific humidity (g/g)
    REAL(SP)                             :: AIRPRES    ! air pressure (Pa)
    REAL(SP)                             :: SWDOWN     ! downward sw radiation (W m-2)
    REAL(SP)                             :: NETRAD     ! net radiation (W m-2)
 END TYPE ADATA
 
 ! the forcing data structure (will have a spatial dimension)
 TYPE FDATA
    REAL(SP)                             :: PPT        ! water input: rain + melt (mm day-1)
    REAL(SP)                             :: TEMP       ! temperature for snow model (deg.C)
    REAL(SP)                             :: PET        ! energy input: potential ET (mm day-1)
 ENDTYPE FDATA

END MODULE multiforce_types
