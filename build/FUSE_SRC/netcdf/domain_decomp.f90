module domain_decomp_module

  use nrtype
  use data_types, only: domain_info

  implicit none

  private
  public :: read_forcing_dimensions
  public :: get_domain_decomp_indices

contains

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- read forcing dimensions (used for domain decomposition) -----------------------

  subroutine read_forcing_dimensions(ncid, info, ierr, message)
    use netcdf
    USE multiforce,only:vname_aprecip ! name of precip variable
    implicit none
    integer(i4b),      intent(in)    :: ncid
    type(domain_info), intent(inout) :: info
    integer(i4b), intent(out)   :: ierr
    character(*), intent(out)   :: message
   
    integer(i4b) :: ivarid
    integer(i4b), parameter :: ndims=3
    integer(i4b) :: dimids(ndims), dimLen
    associate(nx_global => info%space%nx_global, &
              ny_global => info%space%ny_global, &
              nt_global => info%time%nt_global   )

    ierr=0; message="read_forcing_dimensions/"
  
    ! pick one required variable to identify shape (in this case precip)
    ierr = nf90_inq_varid(ncid, trim(vname_aprecip), ivarid)

    ! get dimension IDs (x,y,t)
    ierr = nf90_inquire_variable(ncid, ivarid, dimids=dimids)
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr)); return; endif
   
    ! get dimsension lengths (nx,ny,nt)
    ierr = nf90_inquire_dimension(ncid, dimids(1), len=nx_global); if(ierr/=0) return
    ierr = nf90_inquire_dimension(ncid, dimids(2), len=ny_global); if(ierr/=0) return
    ierr = nf90_inquire_dimension(ncid, dimids(3), len=nt_global); if(ierr/=0) return

    end associate
  end subroutine read_forcing_dimensions

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- get indices to decompose the spatial domain -----------------------------------
  ! 1) Determine global run mode (grid vs catchment)
  ! 2) Apply MPI decomposition (y dimension) and store local dims + offsets

  subroutine get_domain_decomp_indices(domain)
    use data_types, only: domain_type
    implicit none
    type(domain_type), intent(inout) :: domain
   
    associate(&
        grid_flag      => domain%info%space%grid_flag,      &
        nx_global      => domain%info%space%nx_global,      &
        ny_global      => domain%info%space%ny_global,      &
        nx_local       => domain%info%space%nx_local,       &
        ny_local       => domain%info%space%ny_local,       &
        y_start_global => domain%info%space%y_start_global, &
        y_end_global   => domain%info%space%y_end_global,   &
        mpi_enabled    => domain%info%mpi%enabled,          &
        nproc          => domain%info%mpi%nproc,            &
        rank           => domain%info%mpi%rank   )
   
    ! Set flag to toggle between grid and lumped catchment modes
    grid_flag = (nx_global>1 .or. ny_global>1) 
   
    ! Copy globals
    nx_local = nx_global
    ny_local = ny_global
    y_start_global = 1
   
    ! Get indices for split dimensions
    if(grid_flag .and. mpi_enabled .and. nproc>1) then
      call split_1d(ny_global, rank, nproc, &  ! input
                    y_start_global, ny_local)  ! output
    endif
    y_end_global = y_start_global + ny_local - 1
   
    end associate
  end subroutine get_domain_decomp_indices

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- split the dimensions for each MPI rank ----------------------------------------
  ! Purpose: Split domain to allow for MPI.
  !          Given rank, nproc, and n_global, provide start and n_local indices
  ! Creator: Ethan Gutmann, 2020
  ! Modified by Martyn Clark to simplify code and input/output, 12/2025

  subroutine split_1d(n_global, rank, nproc, start, n_local, verbose)
  use nrtype
  implicit none
  integer(i4b), intent(in)  :: n_global, rank, nproc
  logical(lgt), intent(in), optional :: verbose
  integer(i4b), intent(out) :: start, n_local

  integer(i4b) :: base, extra
  logical(lgt) :: talk

  talk = .false.; if(present(verbose)) talk = verbose

  ! --- sanity checks ---
  if(nproc   <= 0)                    stop "split_1d: nproc must be > 0"
  if(rank     < 0 .or. rank >= nproc) stop "split_1d: rank out of range"
  if(n_global < 1)                    stop "split_1d: n_global must be >= 1"

  base  = n_global / nproc                      ! floor(n_global / nproc) rows per rank
  extra = mod(n_global, nproc)                  ! remainder; first 'extra' ranks get +1 row

  n_local = base + merge(1, 0, rank < extra)    ! add 1 row for ranks 0..extra-1
  start   = rank*base + min(rank, extra) + 1    ! shift start by #extra rows assigned before this rank

  if(talk) then
    write(*,'(a,i0,a,i0)') "split_1d: nproc=", nproc, " rank  =", rank
    write(*,'(a,i0,a,i0)') "split_1d: base =", base,  " extra =", extra
    write(*,'(a,i0,a,i0)') "split_1d: start=", start, " nLocal=", n_local
    write(*,'(a,i0,a,i0)') "split_1d: global rows=", start, ":", start+n_local-1
  endif
  end subroutine split_1d

end module domain_decomp_module
