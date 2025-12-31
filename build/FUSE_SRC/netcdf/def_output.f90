MODULE DEF_OUTPUT_MODULE

  USE nrtype                                            ! variable types, etc.
  
  implicit none

  private
  public :: DEF_OUTPUT

  contains

  SUBROUTINE DEF_OUTPUT(nSpat1,nSpat2,n_bands,NUMPAR,NTIM)

  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Martyn Clark to include elevation bands, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Define NetCDF output files -- time-varying model output
  ! ---------------------------------------------------------------------------------------

  ! subroutines
  USE metaoutput, only: VARDESCRIBE                     ! define metadata for model variables
  
  ! data modules
  USE globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
  USE metaoutput, only: NOUTVAR                         ! number of output variables
  USE metaoutput, only: VNAME, LNAME, VUNIT             ! metadata for all model variables
  USE metaoutput, only: isBand, isFlux                  ! logical flag to define vars with band/flux dimension
  USE model_defn, only: FNAME_NETCDF_RUNS               ! model definition (includes filename)
  USE fuse_fileManager, only: Q_ONLY                    ! only write streamflow to output file?
  USE multiforce, only: GRID_FLAG                       ! .true. if distributed
  USE multiforce, only: latitude,longitude              ! dimension arrays
  USE multiforce, only: name_psets,time_steps           ! dimension arrays
  USE multiforce, only: latUnits,lonUnits               ! lat/lon units string
  USE multiforce, only: timeUnits                       ! time units string
  USE globaldata, only: ncid_out                        ! NetCDF output file ID

  IMPLICIT NONE

  ! input
  INTEGER(I4B), INTENT(IN)               :: NTIM           ! number of time steps
  INTEGER(I4B), INTENT(IN)               :: nSpat1,nSpat2  ! length of spatial dimensions
  INTEGER(I4B), INTENT(IN)               :: n_bands        ! number of elevation bands
  INTEGER(I4B), INTENT(IN)               :: NUMPAR         ! number of model parameters

  ! internal
  integer(i4b), dimension(n_bands)       :: band_i               ! coordinate variable
  integer(i4b), dimension(NUMPAR)        :: param_i              ! coordinate variable
  REAL(MSP),DIMENSION(nspat1)            :: longitude_msp        ! coordinate variable (SINGLE PRECISION)
  REAL(MSP),DIMENSION(nspat2)            :: latitude_msp         ! coordinate variable (SINGLE PRECISION)
  REAL(SP),parameter                     :: NA_VALUE_OUT= -9999. ! NA_VALUE for output file
  REAL(MSP)                              :: NA_VALUE_OUT_MSP     ! NA_VALUE for output file

  LOGICAL(LGT)                           :: WRITE_VAR   ! used to denote if the variable is written
  INTEGER(I4B)                           :: IERR        ! error code
  INTEGER(I4B)                           :: NTIM_DIM    ! time
  INTEGER(I4B)                           :: lon_dim     ! 1st spatial dimension
  INTEGER(I4B)                           :: lat_dim     ! 2nd spatial dimension
  INTEGER(I4B)                           :: par_dim     ! parameter dimension
  INTEGER(I4B)                           :: band_dim    ! band dimension
  INTEGER(I4B), DIMENSION(3)             :: TVAR        ! dimension list: exclude band, param
  INTEGER(I4B), DIMENSION(4)             :: EVAR        ! dimension list: include band
  INTEGER(I4B), DIMENSION(4)             :: PVAR        ! dimension list: include param
  integer(i4b)                           :: ib          ! loop through bands
  integer(i4b)                           :: ip          ! loop through parameters
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID

  include 'netcdf.inc'                                  ! use netCDF libraries

  ! ---------------------------------------------------------------------------------------
  CALL VARDESCRIBE()  ! get list of variable descriptions
  ! ---------------------------------------------------------------------------------------
  
  ! put file in define mode
  print *, 'Create NetCDF file for runs:'
  PRINT *, FNAME_NETCDF_RUNS

  IERR = NF_CREATE(TRIM(FNAME_NETCDF_RUNS),NF_CLOBBER,ncid_out); CALL HANDLE_ERR(IERR)

  ! define dimensions
  IERR = NF_DEF_DIM(ncid_out, 'time', NF_UNLIMITED, NTIM_DIM);   CALL HANDLE_ERR(IERR) !record dimension (unlimited length)
  IERR = NF_DEF_DIM(ncid_out, 'band',       n_bands, band_dim);  CALL HANDLE_ERR(IERR)
  IERR = NF_DEF_DIM(ncid_out, 'param',      NUMPAR,  par_dim);   CALL HANDLE_ERR(IERR)
  IERR = NF_DEF_DIM(ncid_out, 'longitude',  nSpat1,  lon_dim);   CALL HANDLE_ERR(IERR)
  IERR = NF_DEF_DIM(ncid_out, 'latitude',   nSpat2,  lat_dim);   CALL HANDLE_ERR(IERR)

  ! define dimension vector
  TVAR = (/lon_dim, lat_dim, NTIM_DIM/) 
  PVAR = (/lon_dim, lat_dim, par_dim,  NTIM_DIM/) 
  EVAR = (/lon_dim, lat_dim, band_dim, NTIM_DIM/) 

  ! define time-varying output variables
  DO IVAR=1,NOUTVAR

    ! check if there is a need to write the variable - see also put_output
    ! uncomment variables that should be written to output file
    IF (Q_ONLY) THEN
      WRITE_VAR=.FALSE.
      IF (TRIM(VNAME(IVAR)).EQ.'q_instnt') WRITE_VAR=.TRUE.
      IF (TRIM(VNAME(IVAR)).EQ.'q_routed') WRITE_VAR=.TRUE.
      IF (.NOT.WRITE_VAR) CYCLE ! start new iteration of do loop, i.e. skip writting variable
    ENDIF

    ! write the variable
    if(isBand(iVar))then
      IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)),NF_REAL,4,EVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
    ELSE
      IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)),NF_REAL,3,TVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
    ENDIF

    ! define missing value
    NA_VALUE_OUT_MSP=NA_VALUE_OUT
    
    ! write metadata
    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(LNAME(IVAR)),TRIM(LNAME(IVAR)));  CALL HANDLE_ERR(IERR)
    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(VUNIT(IVAR)),TRIM(VUNIT(IVAR)));      CALL HANDLE_ERR(IERR)
    IERR = NF_PUT_ATT_REAL(ncid_out,IVAR_ID,'_FillValue',NF_FLOAT,1,NA_VALUE_OUT_MSP);             CALL HANDLE_ERR(IERR)
    
    ! define the parameter sensitivity for each flux: extra variable
    if(isFlux(iVar))then
      IERR = NF_DEF_VAR(ncid_out,TRIM(VNAME(IVAR)//'__dFlux_dParam'),NF_REAL,4,PVAR,IVAR_ID); CALL HANDLE_ERR(IERR)
      IERR = NF_PUT_ATT_REAL(ncid_out,IVAR_ID,'_FillValue',NF_FLOAT,1,NA_VALUE_OUT_MSP);      CALL HANDLE_ERR(IERR)
    endif

  END DO  ! ivar

  ! define the time variable
  ierr = nf_def_var(ncid_out,'time',nf_real,1,(/ntim_dim/),ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',len_trim(timeUnits),trim(timeUnits))
  call handle_err(ierr)

  ! define the latitude variable
  ierr = nf_def_var(ncid_out,'latitude',nf_real,1,(/lat_dim/),ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',8,'degreesN'); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'axis',1,'Y'); call handle_err(ierr)

  ! define the longitude variable
  ierr = nf_def_var(ncid_out,'longitude',nf_real,1,(/lon_dim/),ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',8,'degreesE'); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'axis',1,'X'); call handle_err(ierr)

  ! define the parameter set variable
  ierr = nf_def_var(ncid_out,'param',nf_int,1,(/band_dim/),ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',1,'-'); call handle_err(ierr)

  ! define the band variable
  ierr = nf_def_var(ncid_out,'band',nf_int,1,(/band_dim/),ivar_id); call handle_err(ierr)
  ierr = nf_put_att_text(ncid_out,ivar_id,'units',1,'-'); call handle_err(ierr)

  ! add global attributes
  ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "software",        len("FUSE"),              "FUSE");               call HANDLE_ERR(ierr)
  ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_version",    len_trim(FUSE_VERSION),   trim(FUSE_VERSION));   call HANDLE_ERR(ierr)
  ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_build_time", len_trim(FUSE_BUILDTIME), trim(FUSE_BUILDTIME)); call HANDLE_ERR(ierr)
  ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_git_branch", len_trim(FUSE_GITBRANCH), trim(FUSE_GITBRANCH)); call HANDLE_ERR(ierr)
  ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_git_hash",   len_trim(FUSE_GITHASH),   trim(FUSE_GITHASH));   call HANDLE_ERR(ierr)

  ! end definitions
  IERR = NF_ENDDEF(ncid_out); call handle_err(ierr)

  latitude_msp=latitude ! convert to actual single precision
  IERR = NF_INQ_VARID(ncid_out,'latitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat2,latitude_msp); CALL HANDLE_ERR(IERR) ! write data

  longitude_msp=longitude ! convert to actual single precision
  IERR = NF_INQ_VARID(ncid_out,'longitude',IVAR_ID); CALL HANDLE_ERR(IERR) ! get variable ID
  IERR = NF_PUT_VARA_REAL(ncid_out,IVAR_ID,1,nspat1,longitude_msp); CALL HANDLE_ERR(IERR) ! write data

  band_i = [(ib, ib=1,n_bands)]   ! 1..n_bands
  ierr = NF_INQ_VARID(ncid_out, 'band', ivar_id); call HANDLE_ERR(ierr)
  ierr = NF_PUT_VARA_INT(ncid_out, ivar_id, (/1/), (/n_bands/), band_i); call HANDLE_ERR(ierr)

  param_i = [(ip, ip=1,NUMPAR)]   ! 1..NUMPAR
  ierr = NF_INQ_VARID(ncid_out, 'param', ivar_id); call HANDLE_ERR(ierr)
  ierr = NF_PUT_VARA_INT(ncid_out, ivar_id, (/1/), (/NUMPAR/), param_i); call HANDLE_ERR(ierr)

  PRINT *, 'NetCDF file for model runs defined with dimensions', n_bands, nSpat1 , nSpat2, NUMPAR, NTIM

  ! close output file
  IERR = NF_CLOSE(ncid_out)

  stop "DEF_OUTPUT"

  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE DEF_OUTPUT

END MODULE DEF_OUTPUT_MODULE
