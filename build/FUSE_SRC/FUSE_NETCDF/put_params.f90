MODULE PUT_PARAMS_MODULE

  USE nrtype                                            ! variable types, etc.

  implicit none

  private
  public :: PUT_PARAMS

  contains

  SUBROUTINE PUT_PARAMS(IPAR)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Nans Addor to include snow module
  ! Modified by Martyn Clark to write snow bands as a vector, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! write NetCDF output files  -- model parameters
  ! ---------------------------------------------------------------------------------------
  USE model_defn, only: FNAME_NETCDF_PARA               ! model definition structures (includes filename)
  USE metaparams, only: NOUTPAR                         ! number of model parameters
  USE metaparams, only: PNAME, PDESC, PUNIT             ! metadata for all model parameters
  USE metaparams, only: isBand                          ! logical flag to define vars with elevation dimension
  USE multibands, only: MBANDS, N_BANDS                 ! information for elevation bands
  USE parextract_module                                 ! extract parameters
  IMPLICIT NONE
  ! input
  INTEGER(I4B), INTENT(IN)               :: IPAR        ! parameter set index
  ! internal
  INTEGER(I4B)                           :: IERR,NCID   ! error code; NetCDF ID
  INTEGER(I4B), DIMENSION(1)             :: INDX        ! indices for parameter write
  integer(i4b), dimension(2)             :: start2      ! 2-d start vector
  integer(i4b), dimension(2)             :: count2      ! 2-d count vector
  INTEGER(I4B)                           :: IVAR        ! loop through parameters
  REAL(SP)                               :: XPAR        ! desired parameter
  REAL(MSP)                              :: APAR        ! convert to SP (need for SP write)
  integer(i4b)                           :: ib          ! index of elevation bands
  REAL(SP)     , DIMENSION(N_BANDS)      :: XVEC        ! desired vector
  REAL(MSP)    , DIMENSION(N_BANDS)      :: AVEC        ! convert to SP (need for SP write)
  INTEGER(I4B)                           :: IVAR_ID     ! variable ID
  include 'netcdf.inc'                                  ! use netCDF libraries
  ! ---------------------------------------------------------------------------------------
  
  ! open file
  IERR = NF_OPEN(TRIM(FNAME_NETCDF_PARA),NF_WRITE,NCID)
  CALL HANDLE_ERR(IERR)
  
   ! define indices for model output
   INDX = (/IPAR/)
  
   ! loop through model parameters
   DO IVAR=1,NOUTPAR  ! NOUTPAR is stored in module metaparams
  
    ! get variable ID
    IERR = NF_INQ_VARID(NCID,TRIM(PNAME(IVAR)),IVAR_ID)
    CALL HANDLE_ERR(IERR)
    
    ! standard scalar parameters
    if(.not.isBand(iVar))then
  
      ! extract parameter and write data
      XPAR = PAREXTRACT(PNAME(IVAR)); APAR=XPAR ! get parameter PNAME(IVAR)
      IERR = NF_PUT_VAR1_REAL(NCID, IVAR_ID, INDX, APAR); CALL HANDLE_ERR(IERR)     ! write data
  
    ! elevation band parameters
    else
  
      ! extract vector
      select case (trim(PNAME(IVAR)))
        case ('AF')   ; xVec(1:n_bands) = [ (MBANDS(ib)%AF,    ib=1,n_bands) ]
        case ('Z_MID'); xVec(1:n_bands) = [ (MBANDS(ib)%Z_MID, ib=1,n_bands) ]
        case default; stop "put_params.f90: cannot identify elevation band variable"
      end select
      aVec = xVec ! use MSP to write single precision
  
      ! write row at par=IPAR
      start2 = (/ IPAR, 1 /)
      count2 = (/ 1, n_bands /)
      IERR = NF_PUT_VARA_REAL(NCID, IVAR_ID, start2, count2, aVec(1:n_bands))
      CALL HANDLE_ERR(IERR)
  
    endif  ! elevation band switch
  
   END DO  ! (ivar)
  
  ! close NetCDF file
  IERR = NF_CLOSE(NCID)
  ! ---------------------------------------------------------------------------------------
  END SUBROUTINE PUT_PARAMS

END MODULE PUT_PARAMS_MODULE
