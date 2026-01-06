MODULE model_defn

 ! ---------------------------------------------------------------------------------------
 ! Creator:
 ! --------
 ! Martyn Clark
 ! Modified by Brian Henn to include snow model, 6/2013
 ! Modified by Martyn Clark to separate type definitions from data storage, 01/2026
 ! ---------------------------------------------------------------------------------------
 
 USE nrtype
 USE model_defn_types, only: DESC, UMODEL, SNAMES, FNAMES

 USE globaldata, only: FUSE_VERSION 
 
 implicit none
 private

 public :: NDEC, NTDH_MAX, NSTATE, N_FLUX
 public :: LIST_RFERR, LIST_ARCH1, LIST_ARCH2, LIST_QSURF, LIST_QPERC, LIST_ESOIL, LIST_QINTF, LIST_Q_TDH, LIST_SNOWM
 public :: FNAME_PREFIX, FNAME_TEMPRY, FNAME_ASCII
 public :: FNAME_NETCDF_RUNS, FNAME_NETCDF_PARA, FNAME_NETCDF_PARA_SCE, FNAME_NETCDF_PARA_PRE
 public :: AMODL, SMODL, CSTATE, C_FLUX

 ! list of combinations in each model component
 INTEGER, PARAMETER :: NDEC = 9                           ! number of model decisions
 TYPE(DESC), DIMENSION(2)              :: LIST_RFERR      ! rainfall error
 TYPE(DESC), DIMENSION(3)              :: LIST_ARCH1      ! upper-layer architecture
 TYPE(DESC), DIMENSION(4)              :: LIST_ARCH2      ! lower-layer architecture
 TYPE(DESC), DIMENSION(3)              :: LIST_QSURF      ! surface runoff
 TYPE(DESC), DIMENSION(3)              :: LIST_QPERC      ! percolation
 TYPE(DESC), DIMENSION(2)              :: LIST_ESOIL      ! evaporation
 TYPE(DESC), DIMENSION(2)              :: LIST_QINTF      ! interflow
 TYPE(DESC), DIMENSION(2)              :: LIST_Q_TDH      ! time delay in runoff
 TYPE(DESC), DIMENSION(2)              :: LIST_SNOWM      ! snow model
 
 ! max steps in routing function
 INTEGER(I4B),PARAMETER::NTDH_MAX=500
 
 ! model definitions
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_RUNS    ! NETCDF output filename for model runs
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA    ! NETCDF output filename for model parameters
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA_SCE   ! NETCDF output filename for model parameters produced by SCE
 CHARACTER(LEN=256)                    :: FNAME_NETCDF_PARA_PRE   ! NETCDF filename for pre-defined model parameters set
 CHARACTER(LEN=256)                    :: FNAME_PREFIX    ! prefix for desired output files
 CHARACTER(LEN=256)                    :: FNAME_TEMPRY    ! prefix for temporary output files
 CHARACTER(LEN=256)                    :: FNAME_ASCII     ! ASCII output filename
 TYPE(UMODEL),DIMENSION(5000)          :: AMODL           ! (model definition -- all)
 TYPE(UMODEL)                          :: SMODL           ! (model definition -- single model)
 TYPE(SNAMES),DIMENSION(7)             :: CSTATE          ! (list of model states for SMODL)
 TYPE(FNAMES),DIMENSION(50)            :: C_FLUX          ! (list of model fluxes for SMODL)
 INTEGER(I4B)                          :: NSTATE=0        ! number of model states
 INTEGER(I4B)                          :: N_FLUX=0        ! number of model fluxes
 ! --------------------------------------------------------------------------------------

END MODULE model_defn
