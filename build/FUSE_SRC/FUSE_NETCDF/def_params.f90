MODULE DEF_PARAMS_MODULE

  USE nrtype                                            ! variable types, etc.

  implicit none

  private
  public :: DEF_PARAMS

  contains

  SUBROUTINE DEF_PARAMS(NPAR)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Nans Addor to include snow routine
  ! Modified by Matyn Clark to include band dimension, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Define NetCDF output files -- parameter variables
  ! ---------------------------------------------------------------------------------------
  
  ! subroutines
  USE metaparams, only: PARDESCRIBE                     ! define metadata for model parameters
  
  ! data modules
  USE metaparams, only: NOUTPAR                         ! number of model parameters
  USE metaparams, only: PNAME, PDESC, PUNIT             ! metadata for all model parameters
  USE metaparams, only: isBand                          ! logical flag to define vars with elevation dimension
  USE model_defn, only: FNAME_NETCDF_PARA               ! model definition (includes filename)
  USE multistats, ONLY: MSTATS                          ! model statistics structure
  USE multibands, ONLY: N_BANDS                         ! number of elevation bands
  USE globaldata, only: ncid_out                        ! NetCDF output file ID
  USE globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH

  IMPLICIT NONE
  
  ! input
  INTEGER(I4B), INTENT(IN)               :: NPAR        ! number of parameter sets
  
  ! internal
  INTEGER(I4B)                           :: IERR        ! error code
  INTEGER(I4B)                           :: PAR_DIM     ! parameter set dimension
  INTEGER(I4B)                           :: BAND_DIM    ! elevation band dimension
  INTEGER(I4B), DIMENSION(1)             :: DIMS1       ! 1-d parameter vector
  INTEGER(I4B), DIMENSION(2)             :: DIMS2       ! 2-d parameter-bands matrix
  INTEGER(I4B)                           :: IVAR        ! loop through variables
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID
  
  include 'netcdf.inc'                                  ! use netCDF libraries
  
  ! ---------------------------------------------------------------------------------------
  CALL PARDESCRIBE()               ! get list of parameter descriptions
  ! ---------------------------------------------------------------------------------------
  
  PRINT *, 'Define NetCDF output files - parameter variables = ', TRIM(FNAME_NETCDF_PARA)
  
  ! Create file
  IERR = NF_CREATE(TRIM(FNAME_NETCDF_PARA),NF_CLOBBER,ncid_out); CALL HANDLE_ERR(IERR)
   
   ! define dimensions
   IERR = NF_DEF_DIM(ncid_out, 'par',  NPAR,    PAR_DIM);  CALL HANDLE_ERR(IERR)
   IERR = NF_DEF_DIM(ncid_out, 'band', N_BANDS, BAND_DIM); CALL HANDLE_ERR(IERR)
   
   ! assign dimensions to indices
   DIMS1 = (/PAR_DIM/)              ! 1-d parameter vector 
   DIMS2 = (/PAR_DIM, BAND_DIM/)    ! 2-d parameter-bands matrix
   
   ! define fixed output variables
   DO IVAR=1,NOUTPAR
  
    ! define variables
    if(isBand(iVar))then
     IERR = NF_DEF_VAR(ncid_out, TRIM(PNAME(IVAR)), NF_REAL, 2, DIMS2, IVAR_ID); CALL HANDLE_ERR(IERR)
    else
     IERR = NF_DEF_VAR(ncid_out, TRIM(PNAME(IVAR)), NF_REAL, 1, DIMS1, IVAR_ID); CALL HANDLE_ERR(IERR)
    endif
   
    ! define metadata
    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'long_name',LEN_TRIM(PDESC(IVAR)),TRIM(PDESC(IVAR))); CALL HANDLE_ERR(IERR)
    IERR = NF_PUT_ATT_TEXT(ncid_out,IVAR_ID,'units',LEN_TRIM(PUNIT(IVAR)),TRIM(PUNIT(IVAR)));     CALL HANDLE_ERR(IERR)
    IERR = NF_PUT_ATT_REAL(ncid_out,IVAR_ID,'_FillValue',NF_REAL,1,-9999.);                       CALL HANDLE_ERR(IERR)
  
    END DO  ! ivar
  
    ! add global attributes
    ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "software",        len("FUSE"),              "FUSE");               call HANDLE_ERR(ierr)
    ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_version",    len_trim(FUSE_VERSION),   trim(FUSE_VERSION));   call HANDLE_ERR(ierr)
    ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_build_time", len_trim(FUSE_BUILDTIME), trim(FUSE_BUILDTIME)); call HANDLE_ERR(ierr)
    ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_git_branch", len_trim(FUSE_GITBRANCH), trim(FUSE_GITBRANCH)); call HANDLE_ERR(ierr)
    ierr = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, "fuse_git_hash",   len_trim(FUSE_GITHASH),   trim(FUSE_GITHASH));   call HANDLE_ERR(ierr)
  
  ! end definitions and close file
  IERR = NF_ENDDEF(ncid_out)
  IERR = NF_CLOSE(ncid_out)
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE DEF_PARAMS

END MODULE DEF_PARAMS_MODULE
