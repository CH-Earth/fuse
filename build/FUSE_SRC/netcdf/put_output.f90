MODULE PUT_OUTPUT_MODULE

  USE nrtype  ! variable types, etc.

  implicit none

  private
  public :: PUT_GOUTPUT_3D

  contains

  SUBROUTINE PUT_GOUTPUT_3D(istart_sim,istart_in,numtim)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's 2007 PUT_OUTPUT
  ! Modified by Marytn Clark to use the elevation band dimension, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! write a 3D (or 4D) data structure to the NetCDF output file
  ! ---------------------------------------------------------------------------------------

  ! subroutines
  USE varextract_module, only: VAREXTRACT_3d      ! interface for the function to extract variables

  ! data
  USE model_defn, only: FNAME_NETCDF_RUNS         ! model definition (includes filename)
  USE metaoutput, only: NOUTVAR                   ! number of output variables
  USE metaoutput, only: VNAME, LNAME, VUNIT       ! metadata for all model variables
  USE metaoutput, only: isBand                    ! logical flag to define vars with elevation dimension
  USE multibands, only: MBANDS_VAR_4d, N_BANDS    ! variables for elevation bands
  USE multiforce, only: timDat,time_steps         ! time data
  USE multiforce, only: nspat1,nspat2,startSpat2  ! spatial dimensions
  USE multiforce, only: gForce_3d                 ! test only
  USE multiforce, only: GRID_FLAG                 ! .true. if distributed
  USE globaldata, only: ncid_out                  ! NetCDF output file ID
  USE fuse_fileManager, only: Q_ONLY              ! only write streamflow to output file?

  IMPLICIT NONE

  ! input
  INTEGER(I4B), INTENT(IN)                    :: istart_sim     ! index start time step relative to numtim_sim
  INTEGER(I4B), INTENT(IN)                    :: istart_in      ! index start time step relative to numtim_in - for time dimension
  INTEGER(I4B), INTENT(IN)                    :: numtim         ! number of time steps to write

  ! internal
  LOGICAL(LGT)                                :: WRITE_VAR      ! used to denote if the variable is written
  INTEGER(I4B)                                :: IERR           ! error code
  integer(i4b), dimension(3)                  :: start3         ! start indices: exclude elevation bands
  integer(i4b), dimension(3)                  :: count3         ! count indices: exclude elevation bands
  integer(i4b), dimension(4)                  :: start4         ! start indices: include elevation bands
  integer(i4b), dimension(4)                  :: count4         ! count indices: include elevation bands
  INTEGER(I4B)                                :: IVAR           ! loop through variables
  REAL(SP)                                    :: XVAR           ! desired variable (SP NOT NECESSARILY SP)
  REAL(MSP)                                   :: AVAR           ! desired variable (SINGLE PRECISION)
  REAL(SP),  DIMENSION(nspat1,nspat2,numtim)  :: XVAR_3d        ! desired 3-d variable (SINGLE PRECISION)
  REAL(MSP), DIMENSION(nspat1,nspat2,numtim)  :: AVAR_3d        ! desired 3-d variable (SINGLE PRECISION)
  REAL(SP),  DIMENSION(nspat1,nspat2,n_bands,numtim)  :: XVAR_4d        ! desired 4-d variable (SINGLE PRECISION)
  REAL(MSP), DIMENSION(nspat1,nspat2,n_bands,numtim)  :: AVAR_4d        ! desired 4-d variable (SINGLE PRECISION)
  REAL(MSP), DIMENSION(numtim)                :: tDat           ! time data
  REAL(SP),  DIMENSION(numtim)                :: time_steps_sub ! time data
  INTEGER(I4B)                                :: IVAR_ID        ! variable ID

  INCLUDE 'netcdf.inc'  ! use netCDF libraries


  ! define dimension list (exclude elevation bands)
  ! NOTE: if enabling parallel output you need 1,startSpat2 instead of 1,1 below
  start3 = (/1,1,istart_sim/)
  count3 = (/nspat1,nspat2,numtim/)
  
  ! define dimension list (exclude elevation bands)
  start4 = (/1,1,1,istart_sim/)
  count4 = (/nspat1,nspat2,n_bands,numtim/)

  ! open file
  IERR = NF_OPEN(TRIM(FNAME_NETCDF_RUNS),NF_WRITE,ncid_out)
  CALL HANDLE_ERR(IERR)

  ! loop through variables with time-varying model output
  DO IVAR=1,NOUTVAR

    ! check if there is a need to write the variable - see also def_output
    IF (Q_ONLY) THEN
       WRITE_VAR=.FALSE.
       IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
       IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
       IF (.NOT.WRITE_VAR) CYCLE ! start new iteration of do loop, i.e. skip writting variable
    ENDIF

    ! get variable ID
    IERR = NF_INQ_VARID(ncid_out,TRIM(VNAME(IVAR)),IVAR_ID)
    CALL HANDLE_ERR(IERR)
    
    ! 3-d variables
    if(.not.isBand(iVar))then

      ! write 3-d matrix
      XVAR_3d = VAREXTRACT_3d(VNAME(IVAR), nspat1, nspat2, numtim); AVAR_3d = XVAR_3d ! get variable and convert format
      IERR    = NF_PUT_VARA_REAL(ncid_out, IVAR_ID, start3, count3, AVAR_3d)        ! write data
      CALL HANDLE_ERR(IERR)

    ! 4-d variables
    else

      ! extract variable from 4-D elevation band matrix
      select case (trim(VNAME(IVAR)))
       case ('swe_z'    ); XVAR_4d = MBANDS_VAR_4d(:,:,:,1:numtim)%SWE
       case ('snwacml_z'); XVAR_4d = MBANDS_VAR_4d(:,:,:,1:numtim)%SNOWACCMLTN
       case ('snwmelt_z'); XVAR_4d = MBANDS_VAR_4d(:,:,:,1:numtim)%SNOWMELT
       case default; stop "put_output.f90: cannot identify elevation band variable: "//trim(VNAME(IVAR))
      end select
      aVar_4d = xVar_4d ! use MSP to write single precision

      ! write 4-d matrix
      IERR = NF_PUT_VARA_REAL(ncid_out, IVAR_ID, START4, COUNT4, AVAR_4d)
      call HANDLE_ERR(IERR)

    endif  ! (switch between 3-d and 4-d variables)

  END DO  ! (ivar)

  ! write the time
  time_steps_sub = time_steps(istart_in:(istart_in+numtim-1)) ! extract time for subperiod
  tDat = time_steps_sub ! convert to actual single precision
  ierr = nf_inq_varid(ncid_out,'time',ivar_id); CALL handle_err(ierr)             ! get variable ID for time
  ierr = nf_put_vara_real(ncid_out,ivar_id,(/istart_sim/),(/numtim/),tDat); CALL handle_err(ierr)  ! write time variable

  ! close NetCDF file
  IERR = NF_CLOSE(ncid_out)

  END SUBROUTINE PUT_GOUTPUT_3D

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  SUBROUTINE PUT_OUTPUT(iSpat1, iSpat2, ITIM)

  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! write NetCDF output files
  ! ---------------------------------------------------------------------------------------

  ! subroutines
  USE varextract_module, only: VAREXTRACT         ! interface for the function to extract variables

  ! data
  USE model_defn, only: FNAME_NETCDF_RUNS         ! model definition (includes filename)
  USE metaoutput, only: NOUTVAR                   ! number of output variables
  USE metaoutput, only: VNAME, LNAME, VUNIT       ! metadata for all model variables
  USE metaoutput, only: isBand                    ! logical flag to define vars with elevation dimension
  USE multibands, only: MBANDS, N_BANDS           ! variables for elevation bands
  USE multiforce, only: timDat,time_steps         ! time data
  USE multiforce, only: nspat1,nspat2,startSpat2  ! spatial dimensions
  USE multiforce, only: gForce_3d                 ! test only
  USE multiforce, only: GRID_FLAG                 ! .true. if distributed
  USE globaldata, only: ncid_out                  ! NetCDF output file ID
  USE fuse_fileManager, only: Q_ONLY              ! only write streamflow to output file?

  IMPLICIT NONE
  ! input
  INTEGER(I4B), INTENT(IN)               :: iSpat1      ! index of 1st spatial dimension
  INTEGER(I4B), INTENT(IN)               :: iSpat2      ! index of 2nd spatial dimension
  INTEGER(I4B), INTENT(IN)               :: ITIM        ! time step index
  ! internal
  LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
  INTEGER(I4B)                           :: IERR        ! error code
  !INTEGER(I4B), DIMENSION(5)             :: INDX        ! indices for time series write
  INTEGER(I4B), DIMENSION(3)             :: INDX        ! indices for time series write
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  REAL(SP)                               :: XVAR        ! desired variable (SP NOT NECESSARILY SP)
  REAL(MSP)                              :: AVAR        ! desired variable (SINGLE PRECISION)
  REAL(MSP)                              :: tDat        ! time data
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID
  INCLUDE 'netcdf.inc'                                  ! use netCDF libraries
  ! ---------------------------------------------------------------------------------------
  
  ! open file
  IERR = NF_OPEN(TRIM(FNAME_NETCDF_RUNS),NF_WRITE,ncid_out); CALL HANDLE_ERR(IERR)

  ! define indices for model output
  INDX = (/iSpat1, iSpat2, ITIM/)

  ! loop through time-varying model output
  DO IVAR=1,NOUTVAR

     ! check if there is a need to write the variable - see also def_output
     IF (Q_ONLY) THEN
        WRITE_VAR=.FALSE.
        IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
        IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
        IF (.NOT.WRITE_VAR) CYCLE
     ENDIF

     ! write the variable
     XVAR = VAREXTRACT(VNAME(IVAR)); AVAR=XVAR                                  ! get variable ivar
     IERR = NF_INQ_VARID(ncid_out,TRIM(VNAME(IVAR)),IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
     IERR = NF_PUT_VAR1_REAL(ncid_out,IVAR_ID,INDX,AVAR); CALL HANDLE_ERR(IERR)     ! write data

  END DO  ! (ivar)

  ! write the time
  tDat = timDat%dtime ! convert to actual single precision
  ierr = nf_inq_varid(ncid_out,'time',ivar_id); CALL handle_err(ierr)        ! get variable ID for time
  ierr = nf_put_var1_real(ncid_out,ivar_id,(/itim/),tDat); CALL handle_err(ierr) ! write time variable

  ! close NetCDF file
  IERR = NF_CLOSE(ncid_out)

  END SUBROUTINE PUT_OUTPUT
  
END MODULE PUT_OUTPUT_MODULE
