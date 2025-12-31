module alloc_domain_module

  USE nrtype
  USE data_types, only: domain_type

  implicit none
  private
  public :: allocate_domain_data
  public :: set_legacy_arrays

CONTAINS

  subroutine allocate_domain_data(domain, ierr, message)

  implicit none

  type(domain_type), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b) :: nx, ny, nt, nb

  ierr=0; message="allocate_domain_data/"

  ! define dimensions
  nx = domain%info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = domain%info%space%ny_local
  nt = domain%info%time%nt_window
  nb = domain%info%snow%n_bands

  ! allocate validity mask
  allocate(domain%data%valid(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate valid"; return; endif

  ! allocate ancillary forcing
  allocate(domain%data%ancil(nx,ny), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate ancil"; return; endif

  ! allocate forcing window
  allocate(domain%data%force(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate force"; return; endif

  ! allocate state window
  allocate(domain%data%state(nx,ny,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate state"; return; endif

  ! allocate flux window
  allocate(domain%data%flux(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate flux"; return; endif

  ! allocate basin averages
  allocate(domain%data%aForce(nt), domain%data%aRoute(nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate aForce/aRoute"; return; endif

  ! allocate routing if needed
  allocate(domain%data%route(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate route"; return; endif

  ! allocate bands
  allocate(domain%data%bands(nx,ny,nb,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate bands"; return; endif

  end subroutine allocate_domain_data

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- copy arrays in the domain%data structure to legacy arrays ---------------------

  subroutine set_legacy_arrays(domain)

    ! legacy modules
    use multiforce, only: nSpat1, nSpat2, numtim_sub
    use multiForce, only: gForce_3d, ancilF, aValid
    use multiState, only: gState_3d
    use multiRoute, only: AROUTE_3d
    use multiBands, only: MBANDS_VAR_4d, N_BANDS
    implicit none

    type(domain_type), intent(inout)  :: domain

    ! ensure the spatial dimensions match what is in domain%info
    nSpat1        = domain%info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
    nSpat2        = domain%info%space%ny_local
    numtim_sub    = domain%info%time%nt_window
    n_bands       = domain%info%snow%n_bands

    ! copy arrays in the domain%data structure to legacy arrays
    aValid        = domain%data%valid   ! validity mask
    ancilF        = domain%data%ancil   ! ancillary forcing
    gForce_3d     = domain%data%force   ! forcing window
    gState_3d     = domain%data%state   ! state window
    AROUTE_3d     = domain%data%route   ! routing window
    MBANDS_VAR_4d = domain%data%bands   ! elevation band window

  end subroutine set_legacy_arrays

end module alloc_domain_module
