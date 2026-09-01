module alloc_domain_module

  USE nrtype
  USE info_types, only: fuse_info
  USE domain_types, only: domain_data
 
  use fuse_globaldata, only: NVAR_HYDROMET

  implicit none
  private
  public :: allocate_domain_data
  public :: set_legacy_arrays

CONTAINS

  subroutine allocate_domain_data(info, domain, ierr, message)

  implicit none

  type(fuse_info),   intent(inout) :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b) :: nx, ny, nt, nb

  ierr=0; message="allocate_domain_data/"

  ! define dimensions
  nx = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = info%space%ny_local
  nt = info%time%nt_window
  nb = info%snow%n_bands

  ! allocate forcing window
  allocate(domain%force(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate force"; return; endif

  ! allocate state window
  allocate(domain%state(nx,ny,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate state"; return; endif

  ! allocate flux window
  allocate(domain%flux(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate flux"; return; endif

  ! allocate routing if needed
  allocate(domain%route(nx,ny,nt), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate route"; return; endif

  ! allocate elevation arrays
  allocate(domain%z_forcing(nx,ny), domain%elev_mask(nx,ny), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elevation arrays"; return; endif

  ! allocate overlap area
  allocate(domain%olap_area(nx,ny), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate overlap area"; return; endif

  ! allocate elevation bands (info)
  allocate(domain%bands_info(nx,ny,nb), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands (info)"; return; endif

  ! allocate elevation bands (var)
  allocate(domain%bands_var(nx,ny,nb,nt+1), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands (var)"; return; endif

  ! allocate hydromet lookup table
  allocate(info%files%hydromet%name(NVAR_HYDROMET),  &
           info%files%hydromet%units(NVAR_HYDROMET), &
           info%files%hydromet%varid(NVAR_HYDROMET), &
           info%files%hydromet%ndims(NVAR_HYDROMET), &
           info%files%hydromet%multiplier(NVAR_HYDROMET), &
           info%files%hydromet%fill_value(NVAR_HYDROMET), &
           info%files%hydromet%has_fill_value(NVAR_HYDROMET), &
           stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate hydromet lookup table"; return; endif
  info%files%hydromet%multiplier(:) = -1._wp
  info%files%hydromet%fill_value(:)     = 0._wp
  info%files%hydromet%has_fill_value(:) = .false.

  end subroutine allocate_domain_data

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- copy arrays in the domain structure to legacy arrays ---------------------

  subroutine set_legacy_arrays(info, domain, ierr, message)

  ! legacy modules
  use multiforce, only: nSpat1, nSpat2, numtim_sub
  USE multiforce, only: startSpat2                         ! starting y index for data read
  USE multiforce, only: ncid_hydromet
  use multiforce, only: timeUnits
  use multiBands, only: N_BANDS, MBANDS
  implicit none

  type(fuse_info),   intent(in)    :: info
  type(domain_data), intent(inout) :: domain
  integer(i4b),      intent(out)   :: ierr
  character(*),      intent(out)   :: message

  integer(i4b) :: nx, ny, nt, nb

  ierr = 0
  message = 'set_legacy_arrays/'

  ! define dimensions
  nx = info%space%nx_local  ! NOTE: local to rank (MPI parallelization)
  ny = info%space%ny_local
  nt = info%time%nt_window
  nb = info%snow%n_bands

  ! set variables in multiforce
  nSpat1     = nx
  nSpat2     = ny
  numtim_sub = nt
  startSpat2 = info%space%y_start_global

  ncid_hydromet  = info%files%ncid_hydromet

  ! set bands
  N_BANDS = nb

  ! allocate elevation bands
  allocate(MBANDS(nb), stat=ierr)
  if(ierr/=0)then; message=trim(message)//"cannot allocate elev bands"; return; endif

  end subroutine set_legacy_arrays

end module alloc_domain_module
