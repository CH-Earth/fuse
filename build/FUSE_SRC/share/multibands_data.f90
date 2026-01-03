MODULE multibands
 
 ! Created by Brian Henn to allow multi-band snow modeling, 6/2013
 ! Based on module MULTIFORCE by Martyn Clark
 
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026

 USE nrtype

 USE multibands_types, only: BANDS, BANDS_INFO, BANDS_VAR

 implicit none
 private

 public :: N_BANDS
 public :: MBANDS, MBANDS_INFO_3d, MBANDS_VAR_4d
 public :: Z_FORCING, Z_FORCING_grid, elev_mask

 ! --------------------------------------------------------------------------------------
 TYPE(BANDS),DIMENSION(:),ALLOCATABLE  :: MBANDS          ! basin band information
 type(BANDS_INFO),dimension(:,:,:),ALLOCATABLE :: MBANDS_INFO_3d    ! basin band information in space
 type(BANDS_VAR),dimension(:,:,:,:),ALLOCATABLE :: MBANDS_VAR_4d    ! basin band information in space plus time

 INTEGER(I4B)                          :: N_BANDS=0       ! number of bands, initialize to zero
 REAL(SP)                              :: Z_FORCING       ! elevation of forcing data (m)
 REAL(SP),DIMENSION(:,:),ALLOCATABLE   :: Z_FORCING_grid  ! elevation of forcing data (m) for the 2D domain
 LOGICAL(LGT),DIMENSION(:,:),ALLOCATABLE   :: elev_mask   ! mask domain - TRUE means the cell must be masked, i.e. not run
 ! --------------------------------------------------------------------------------------

END MODULE multibands
