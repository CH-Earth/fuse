module read_coords_module

  use nrtype
  use netcdf

  implicit none

  private
  public :: read_latlon_1d

contains

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- read 1-D coordinate vectors (lat(y), lon(x)) for a rectilinear 2-D grid -------

  subroutine read_latlon_1d(ncid, domain, ierr, message)
    use data_types, only: domain_type
    implicit none

    integer(i4b),      intent(in)    :: ncid
    type(domain_type), intent(inout) :: domain
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    integer(i4b) :: lat_id, lon_id
    integer(i4b) :: y0, ny, nx
    character(len=32) :: lat_name, lon_name

    ierr = 0
    message = "read_latlon_1d/"

    associate( &
      y0     => domain%info%space%y_start_global, &  ! 1-based index into global forcing file
      ny     => domain%info%space%ny_local,       &  ! number of local y rows
      nx     => domain%info%space%nx_local,       &  ! should equal nx_global (split only along y)
      lat_1d => domain%data%coords%lat_1d,        &
      lon_1d => domain%data%coords%lon_1d )

      ! ----- allocate storage if needed ------------------------------------------------

      if(.not. allocated(domain%data%coords%lat_1d)) print*, 'hello'

      if(.not. allocated(lat_1d)) then
        allocate(lat_1d(ny), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'allocate(lat_1d) failed'; return; endif
      endif

      if(.not. allocated(lon_1d)) then
        allocate(lon_1d(nx), stat=ierr)
        if(ierr/=0)then; message=trim(message)//'allocate(lon_1d) failed'; return; endif
      endif

      ! ----- find latitude and longitude variables -------------------------------------

      lat_name = "latitude"
      ierr = nf90_inq_varid(ncid, trim(lat_name), lat_id)
      if(ierr /= nf90_noerr) then
        lat_name = "lat"
        ierr = nf90_inq_varid(ncid, trim(lat_name), lat_id)
        if(ierr /= nf90_noerr) then
          message = trim(message)//"cannot find latitude variable (latitude or lat)"
          return
        endif
      endif

      lon_name = "longitude"
      ierr = nf90_inq_varid(ncid, trim(lon_name), lon_id)
      if(ierr /= nf90_noerr) then
        lon_name = "lon"
        ierr = nf90_inq_varid(ncid, trim(lon_name), lon_id)
        if(ierr /= nf90_noerr) then
          message = trim(message)//"cannot find longitude variable (longitude or lon)"
          return
        endif
      endif

      ! ----- read latitude slice (local y slab) ----------------------------------------

      ierr = nf90_get_var(ncid, lat_id, lat_1d, start=(/y0/), count=(/ny/))
      if(ierr /= nf90_noerr) then
        message = trim(message)//"nf90_get_var("//trim(lat_name)//") failed: " // &
                  trim(nf90_strerror(ierr))
        return
      endif

      ! ----- read all longitudes (same on every rank) ----------------------------------

      ierr = nf90_get_var(ncid, lon_id, lon_1d, start=(/1/), count=(/nx/))
      if(ierr /= nf90_noerr) then
        message = trim(message)//"nf90_get_var("//trim(lon_name)//") failed: " // &
                  trim(nf90_strerror(ierr))
        return
      endif

    end associate

  end subroutine read_latlon_1d

end module read_coords_module
