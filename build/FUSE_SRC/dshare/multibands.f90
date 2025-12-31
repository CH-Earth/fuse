! Created by Brian Henn to allow multi-band snow modeling, 6/2013
! Based on module MULTIFORCE by Martyn Clark
MODULE multibands
 USE nrtype
 USE data_types, only: BANDS, BANDS_INFO, BANDS_VAR
 ! --------------------------------------------------------------------------------------
 TYPE(BANDS)      , DIMENSION(:)       , allocatable   :: MBANDS          ! basin band information
 type(BANDS_INFO) , dimension(:,:,:)   , allocatable   :: MBANDS_INFO_3d  ! basin band information in space
 type(BANDS_VAR)  , dimension(:,:,:,:) , allocatable   :: MBANDS_VAR_4d   ! basin band information in space plus time
 LOGICAL(LGT)     , DIMENSION(:,:)     , allocatable   :: elev_mask       ! mask domain - TRUE means the cell must be masked, i.e. not run
 REAL(SP)         , DIMENSION(:,:)     , allocatable   :: Z_FORCING_grid  ! elevation of forcing data (m) for the 2D domain
 INTEGER(I4B)                                          :: N_BANDS=0       ! number of bands, initialize to zero
 REAL(SP)                                              :: Z_FORCING       ! elevation of forcing data (m)
 ! --------------------------------------------------------------------------------------
END MODULE multibands
