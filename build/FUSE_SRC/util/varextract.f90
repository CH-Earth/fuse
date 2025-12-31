MODULE VAREXTRACT_MODULE

  USE nrtype               ! variable types, etc.

  IMPLICIT NONE

  private
  public :: VAREXTRACT_3d
  public :: VAREXTRACT

  CONTAINS
   
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  PURE FUNCTION VAREXTRACT_3d(VARNAME,nspat1,nspat2,numtim)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's 2007 VAREXTRACT
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Extracts variable "VARNAME" from relevant data structures
  ! ---------------------------------------------------------------------------------------
  USE model_numerix                                           ! model numerix parameters
  USE globaldata, only: NA_VALUE_SP                           ! missing value
  USE multiforce, only: gForce_3d, aValid                     ! model forcing data
  USE multistate, only: gState_3d                             ! model states
  USE multi_flux, only: w_flux_3d                             ! model fluxes
  USE multiroute, only: aroute_3d                             ! routed runoff
  IMPLICIT NONE
  ! input
  CHARACTER(*), INTENT(IN)                  :: VARNAME        ! variable name
  INTEGER(i4b), INTENT(IN)                  :: nspat1,nspat2  ! number of elements in spat1, spat2 (lon, lat)
  INTEGER(i4b), INTENT(IN)                  :: numtim         ! number of time steps
  ! internal
  real(sp), DIMENSION(nspat1,nspat2,numtim) :: XVAR_3d        ! variable
  integer(i4b)                              :: ierr           ! error code
  CHARACTER(LEN=1024)                       :: MESSAGE        ! error message
  ! output
  real(sp), DIMENSION(nspat1,nspat2,numtim) :: VAREXTRACT_3d  ! FUNCTION name
  
  ! ---------------------------------------------------------------------------------------
  ! the length of the temporal dimension of the state variables (gState_3d and MBANDS_VAR_4d)
  ! is greater by one time step, so only keeping first numtim time steps, i.e. not writing
  ! last value the output file
  
  SELECT CASE (TRIM(VARNAME))
  
   ! extract forcing data
   CASE ('ppt')        ; XVAR_3d = gForce_3d%PPT
   CASE ('temp')       ; XVAR_3d = gForce_3d%TEMP
   CASE ('pet')        ; XVAR_3d = gForce_3d%PET
   
   ! extract response data
   CASE ('obsq')       ; XVAR_3d = aValid%OBSQ
   
   ! extract model states
   CASE ('tens_1')     ; XVAR_3d = gState_3d(:,:,1:numtim)%TENS_1
   CASE ('tens_1a')    ; XVAR_3d = gState_3d(:,:,1:numtim)%TENS_1A
   CASE ('tens_1b')    ; XVAR_3d = gState_3d(:,:,1:numtim)%TENS_1B
   CASE ('free_1')     ; XVAR_3d = gState_3d(:,:,1:numtim)%FREE_1
   CASE ('watr_1')     ; XVAR_3d = gState_3d(:,:,1:numtim)%WATR_1
   CASE ('tens_2')     ; XVAR_3d = gState_3d(:,:,1:numtim)%TENS_2
   CASE ('free_2')     ; XVAR_3d = gState_3d(:,:,1:numtim)%FREE_2
   CASE ('free_2a')    ; XVAR_3d = gState_3d(:,:,1:numtim)%FREE_2A
   CASE ('free_2b')    ; XVAR_3d = gState_3d(:,:,1:numtim)%FREE_2B
   CASE ('watr_2')     ; XVAR_3d = gState_3d(:,:,1:numtim)%WATR_2
   CASE ('swe_tot')    ; XVAR_3d = gState_3d(:,:,1:numtim)%swe_tot
   
   ! extract model fluxes
   CASE ('eff_ppt')    ; XVAR_3d = W_FLUX_3d%EFF_PPT
   CASE ('satarea')    ; XVAR_3d = W_FLUX_3d%SATAREA
   CASE ('qsurf')      ; XVAR_3d = W_FLUX_3d%QSURF
   CASE ('evap_1a')    ; XVAR_3d = W_FLUX_3d%EVAP_1A
   CASE ('evap_1b')    ; XVAR_3d = W_FLUX_3d%EVAP_1B
   CASE ('evap_1')     ; XVAR_3d = W_FLUX_3d%EVAP_1
   CASE ('evap_2')     ; XVAR_3d = W_FLUX_3d%EVAP_2
   CASE ('rchr2excs')  ; XVAR_3d = W_FLUX_3d%RCHR2EXCS
   CASE ('tens2free_1'); XVAR_3d = W_FLUX_3d%TENS2FREE_1
   CASE ('oflow_1')    ; XVAR_3d = W_FLUX_3d%OFLOW_1
   CASE ('tens2free_2'); XVAR_3d = W_FLUX_3d%TENS2FREE_2
   CASE ('qintf_1')    ; XVAR_3d = W_FLUX_3d%QINTF_1
   CASE ('qperc_12')   ; XVAR_3d = W_FLUX_3d%QPERC_12
   CASE ('qbase_2')    ; XVAR_3d = W_FLUX_3d%QBASE_2
   CASE ('qbase_2a')   ; XVAR_3d = W_FLUX_3d%QBASE_2A
   CASE ('qbase_2b')   ; XVAR_3d = W_FLUX_3d%QBASE_2B
   CASE ('oflow_2')    ; XVAR_3d = W_FLUX_3d%OFLOW_2
   CASE ('oflow_2a')   ; XVAR_3d = W_FLUX_3d%OFLOW_2A
   CASE ('oflow_2b')   ; XVAR_3d = W_FLUX_3d%OFLOW_2B
   
   ! extract extrapolation errors
   CASE ('err_tens_1') ; XVAR_3d = W_FLUX_3d%ERR_TENS_1
   CASE ('err_tens_1a'); XVAR_3d = W_FLUX_3d%ERR_TENS_1A
   CASE ('err_tens_1b'); XVAR_3d = W_FLUX_3d%ERR_TENS_1B
   CASE ('err_free_1') ; XVAR_3d = W_FLUX_3d%ERR_FREE_1
   CASE ('err_watr_1') ; XVAR_3d = W_FLUX_3d%ERR_WATR_1
   CASE ('err_tens_2') ; XVAR_3d = W_FLUX_3d%ERR_TENS_2
   CASE ('err_free_2') ; XVAR_3d = W_FLUX_3d%ERR_FREE_2
   CASE ('err_free_2a'); XVAR_3d = W_FLUX_3d%ERR_FREE_2A
   CASE ('err_free_2b'); XVAR_3d = W_FLUX_3d%ERR_FREE_2B
   CASE ('err_watr_2') ; XVAR_3d = W_FLUX_3d%ERR_WATR_2
   
   ! time check
   CASE ('chk_time')   ; XVAR_3d = W_FLUX_3d%CHK_TIME
   
   ! extract model runoff
   CASE ('q_instnt')   ; XVAR_3d = AROUTE_3d%Q_INSTNT
   CASE ('q_routed')   ; XVAR_3d = AROUTE_3d%Q_ROUTED
   
   ! extract information on numerical solution (shared in MODULE model_numerix)
   CASE ('num_funcs')  ; XVAR_3d = NUM_FUNCS
   CASE ('numjacobian'); XVAR_3d = NUM_JACOBIAN
   CASE ('sub_accept') ; XVAR_3d = NUMSUB_ACCEPT
   CASE ('sub_reject') ; XVAR_3d = NUMSUB_REJECT
   CASE ('sub_noconv') ; XVAR_3d = NUMSUB_NOCONV
   CASE ('max_iterns') ; XVAR_3d = MAXNUM_ITERNS

   ! default
   case default;         XVAR_3d = NA_VALUE_SP 

  END SELECT
  
  ! save the output
  VAREXTRACT_3d = XVAR_3d
  
  ! ---------------------------------------------------------------------------------------
  END FUNCTION VAREXTRACT_3d

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------------------
  PURE FUNCTION VAREXTRACT(VARNAME)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Brian Henn to include snow model, 6/2013
  ! Modified by Nans Addor to enable distributed modeling, 9/2016
  ! Modified by Martyn Clark to use dimension for elevation bands, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Extracts variable "VARNAME" from relevant data structures
  ! ---------------------------------------------------------------------------------------
  USE model_numerix                                     ! model numerix parameters
  USE globaldata, only: NA_VALUE_SP                     ! missing value
  USE multiforce, only: MFORCE, valDat                  ! model forcing data
  USE multistate, only: FSTATE                          ! model states
  USE multi_flux, only: W_FLUX                          ! model fluxes
  USE multiroute, only: MROUTE                          ! routed runoff
  IMPLICIT NONE
  ! input
  CHARACTER(*), INTENT(IN)               :: VARNAME     ! variable name
  ! internal
  REAL(SP)                               :: XVAR        ! variable
  ! output
  REAL(SP)                               :: VAREXTRACT  ! FUNCTION name
  ! ---------------------------------------------------------------------------------------
  SELECT CASE (TRIM(VARNAME))
   
   ! extract forcing data
   CASE ('ppt')        ; XVAR = MFORCE%PPT
   CASE ('temp')       ; XVAR = MFORCE%TEMP
   CASE ('pet')        ; XVAR = MFORCE%PET
   
   ! extract response data
   CASE ('obsq')       ; XVAR = valDat%OBSQ
   
   ! extract model states
   CASE ('tens_1')     ; XVAR = FSTATE%TENS_1
   CASE ('tens_1a')    ; XVAR = FSTATE%TENS_1A
   CASE ('tens_1b')    ; XVAR = FSTATE%TENS_1B
   CASE ('free_1')     ; XVAR = FSTATE%FREE_1
   CASE ('watr_1')     ; XVAR = FSTATE%WATR_1
   CASE ('tens_2')     ; XVAR = FSTATE%TENS_2
   CASE ('free_2')     ; XVAR = FSTATE%FREE_2
   CASE ('free_2a')    ; XVAR = FSTATE%FREE_2A
   CASE ('free_2b')    ; XVAR = FSTATE%FREE_2B
   CASE ('watr_2')     ; XVAR = FSTATE%WATR_2
   CASE ('swe_tot')    ; XVAR = FSTATE%swe_tot
   
   ! extract model fluxes
   CASE ('eff_ppt')    ; XVAR = W_FLUX%EFF_PPT
   CASE ('satarea')    ; XVAR = W_FLUX%SATAREA
   CASE ('qsurf')      ; XVAR = W_FLUX%QSURF
   CASE ('evap_1a')    ; XVAR = W_FLUX%EVAP_1A
   CASE ('evap_1b')    ; XVAR = W_FLUX%EVAP_1B
   CASE ('evap_1')     ; XVAR = W_FLUX%EVAP_1
   CASE ('evap_2')     ; XVAR = W_FLUX%EVAP_2
   CASE ('rchr2excs')  ; XVAR = W_FLUX%RCHR2EXCS
   CASE ('tens2free_1'); XVAR = W_FLUX%TENS2FREE_1
   CASE ('oflow_1')    ; XVAR = W_FLUX%OFLOW_1
   CASE ('tens2free_2'); XVAR = W_FLUX%TENS2FREE_2
   CASE ('qintf_1')    ; XVAR = W_FLUX%QINTF_1
   CASE ('qperc_12')   ; XVAR = W_FLUX%QPERC_12
   CASE ('qbase_2')    ; XVAR = W_FLUX%QBASE_2
   CASE ('qbase_2a')   ; XVAR = W_FLUX%QBASE_2A
   CASE ('qbase_2b')   ; XVAR = W_FLUX%QBASE_2B
   CASE ('oflow_2')    ; XVAR = W_FLUX%OFLOW_2
   CASE ('oflow_2a')   ; XVAR = W_FLUX%OFLOW_2A
   CASE ('oflow_2b')   ; XVAR = W_FLUX%OFLOW_2B
   
   ! extract extrapolation errors
   CASE ('err_tens_1') ; XVAR = W_FLUX%ERR_TENS_1
   CASE ('err_tens_1a'); XVAR = W_FLUX%ERR_TENS_1A
   CASE ('err_tens_1b'); XVAR = W_FLUX%ERR_TENS_1B
   CASE ('err_free_1') ; XVAR = W_FLUX%ERR_FREE_1
   CASE ('err_watr_1') ; XVAR = W_FLUX%ERR_WATR_1
   CASE ('err_tens_2') ; XVAR = W_FLUX%ERR_TENS_2
   CASE ('err_free_2') ; XVAR = W_FLUX%ERR_FREE_2
   CASE ('err_free_2a'); XVAR = W_FLUX%ERR_FREE_2A
   CASE ('err_free_2b'); XVAR = W_FLUX%ERR_FREE_2B
   CASE ('err_watr_2') ; XVAR = W_FLUX%ERR_WATR_2
   
   ! time check
   CASE ('chk_time')   ; XVAR = W_FLUX%CHK_TIME
   
   ! extract model runoff
   CASE ('q_instnt')   ; XVAR = MROUTE%Q_INSTNT
   CASE ('q_routed')   ; XVAR = MROUTE%Q_ROUTED
   
   ! extract information on numerical solution (shared in MODULE model_numerix)
   CASE ('num_funcs')  ; XVAR = NUM_FUNCS
   CASE ('numjacobian'); XVAR = NUM_JACOBIAN
   CASE ('sub_accept') ; XVAR = NUMSUB_ACCEPT
   CASE ('sub_reject') ; XVAR = NUMSUB_REJECT
   CASE ('sub_noconv') ; XVAR = NUMSUB_NOCONV
   CASE ('max_iterns') ; XVAR = MAXNUM_ITERNS

   ! default
   case default;         XVAR = NA_VALUE_SP 

   END SELECT
  
  ! and, save the output
  VAREXTRACT = XVAR
  ! ---------------------------------------------------------------------------------------
  END FUNCTION VAREXTRACT



END MODULE VAREXTRACT_MODULE
