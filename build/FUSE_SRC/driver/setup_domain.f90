module setup_domain_module

  USE nrtype
  USE data_types, only: cli_options
  USE data_types, only: domain_type
  USE globaldata, only: isPrint

  implicit none

  private
  public :: setup_domain

contains

  subroutine setup_domain(opts, domain, ierr, message)

  ! access subroutines
  use netcdf, only: nf90_open, nf90_nowrite, nf90_strerror     ! NetCDF functions
  USE fuse_fileManager,     only: fuse_SetDirsUndPhiles        ! sets directories and filenames
  USE fuse_fileManager,     only: export_filemanager_to_domain ! populates domain%info structure
  USE fuse_fileManager,     only: finalize_domain_config       ! compute additional filenames/variables 
  USE fuse_fileManager,     only: export_domain_to_legacy      ! populates legacy modules

  USE force_info_module,    only: force_info                   ! get forcing info for NetCDF files

  USE get_gforce_module,    only: read_ginfo                   ! get dimension lengths from the NetCDF file
  USE get_gforce_module,    only: get_varID                    ! list of var ids
  
  USE get_mbands_module,    only: GET_MBANDS_INFO              ! get elevation bands for snow modeling 
  
  USE domain_decomp_module, only: read_forcing_dimensions      ! get forcing dimensions for MPI domain decomposition
  USE domain_decomp_module, only: get_domain_decomp_indices    ! get MPI domain decomposition indices
  
  USE time_windows_module,  only: get_time_windows             ! get info on the rolling time windows
  USE time_windows_module,  only: export_time_to_multiforce    ! populate legacy multiforce modules

  USE alloc_domain_module,  only: allocate_domain_data         ! allocate space for data arrays in the domain structure
  USE alloc_domain_module,  only: set_legacy_arrays            ! copy arrays in the domain%data structure to legacy arrays 

  ! shared data: TODO move into domain structure
  USE multiforce, only: ncid_forc
  
  implicit none
  
  ! input
  type(cli_options)   , intent(in)                  :: opts            ! command line interface options
  type(domain_type)   , intent(inout)               :: domain          ! the structure that holds "everything"
  
  ! output
  integer(i4b)        , intent(out)                 :: ierr            ! error code
  character(len=1024) , intent(out)                 :: message         ! error message
  
  ! ----- internal -----------------------------------------------------------------------
  CHARACTER(LEN=1024)                               :: CMESSAGE        ! error message
  ! ---------------------------------------------------------------------------------------
  associate(INPUT_PATH => domain%info%files%input_path, &
            forcefile  => domain%info%files%forcing_file)
  ! ---------------------------------------------------------------------------------------
  ierr=0; message='setup_domain/'

  ! ----- set paths and file names --------------------------------------------------------
  
  ! set directories and filenames for control files
  call fuse_SetDirsUndPhiles(fuseFileManagerIn=opts%control_file, err=ierr, message=cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
 
  ! copy global file information to the domain structure
  call export_filemanager_to_domain(opts, domain)

  ! derive filenames + parse config strings 
  call finalize_domain_config(opts, domain, ierr, message)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! populate legacy modules
  call export_domain_to_legacy(domain)

  ! ----- read information on numerical decisions and forcing files -----------------------
  
  ! defines method/parameters used for numerical solution based on numerix file
  ! NOTE: This routine supports the legacy FUSE v1 numerics experiments
  CALL GETNUMERIX(IERR,CMESSAGE)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
  
  ! get forcing info from the txt file
  !   -- forcing info is text that describes the forcing NetCDF files (TODO: needs improvement)
  call force_info(ierr,cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
  if(isPrint) print *, 'Open forcing file:', trim(INPUT_PATH)//trim(forcefile)
  
  ! ----- read grid info and define indices for MPI domain decomposition ------------------
 
  ! open NetCDF forcing file
  ierr = nf90_open(trim(INPUT_PATH)//trim(forcefile), nf90_nowrite, ncid_forc)
  if (ierr/=0)then; message=trim(message)//' nf90_open failed: '//trim(nf90_strerror(ierr)); return; endif
  if(isPrint) PRINT *, 'NCID_FORC is', ncid_forc
 
  ! get NetCDF ID for each variable of the forcing file
  ! NOTE: populates data structures in multiforce
  call get_varID(ncid_forc, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
 
  ! populate domain structure with (x,y,t) dimension lengths
  !   -- nx_global, ny_global, nt_global
  call read_forcing_dimensions(ncid_forc, domain%info, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
 
  ! get indices for MPI decomposition of the spatial domain: y_start_global, ny_local 
  ! NOTE: These indices will be used later to read different subsets of forcing data for different ranks
  call get_domain_decomp_indices(domain)
  
  ! ----- Compute time indices for sim/eval windows and subperiod chunk size --------------
  
  call get_time_windows(ncid_forc, domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! export domain%info%time -> multiforce to keep legacy code working
  call export_time_to_multiforce(domain)

  ! ----- Get information on elevation bands ----------------------------------------------
 
  ! get elevation band info, in particular N_BANDS, from NetCDF file
  CALL GET_MBANDS_INFO(domain%info, ierr, cmessage) 
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif
  
  ! ----- Allocate space for domain data --------------------------------------------------

  ! allocate space for the arrays in the domain%data structure
  call allocate_domain_data(domain, ierr, cmessage)
  if (ierr/=0)then; message=trim(message)//trim(cmessage); ierr=20; return; endif

  ! copy arrays in the domain%data structure to legacy arrays
  call set_legacy_arrays(domain)

  end associate
  end subroutine setup_domain

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

end module setup_domain_module
