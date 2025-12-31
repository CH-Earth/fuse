module get_mbands_module

  USE nrtype

  implicit none

  private
  public :: GET_MBANDS_INFO

contains

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- get the number of elevation bands (used for allocate statements later) --------

  subroutine GET_MBANDS_INFO(info, ierr, message)

    use data_types, only: domain_info
    use netcdf,     only: nf90_open, nf90_nowrite, nf90_close, nf90_inq_dimid, &
                          nf90_inquire_dimension, nf90_strerror
   
    implicit none
   
    type(domain_info), intent(inout) :: info
    integer(i4b)     , intent(out)   :: ierr
    character(*)     , intent(out)   :: message
   
    integer(i4b) :: ncid_eb, dimid_eb, dimLen
    character(len=1024) :: cfile
   
    ierr=0; message="GET_MBANDS_INFO/"
   
    cfile = trim(info%files%input_path)//trim(info%files%elevbands_file)
   
    ierr = nf90_open(cfile, nf90_nowrite, ncid_eb)
    if(ierr/=0) then
      message=trim(message)//"nf90_open failed: "//trim(nf90_strerror(ierr))
      return
    endif
   
    ierr = nf90_inq_dimid(ncid_eb, "elevation_band", dimid_eb)
    if(ierr/=0) then
      message=trim(message)//"nf90_inq_dimid failed: "//trim(nf90_strerror(ierr))
      return
    endif
   
    ierr = nf90_inquire_dimension(ncid_eb, dimid_eb, len=dimLen)
    if(ierr/=0) then
      message=trim(message)//"nf90_inquire_dimension failed: "//trim(nf90_strerror(ierr))
      return
    endif
   
    ierr = nf90_close(ncid_eb)
    if(ierr/=0) then
      message=trim(message)//"nf90_close failed: "//trim(nf90_strerror(ierr))
      return
    endif
   
    info%snow%n_bands = dimLen

  end subroutine GET_MBANDS_INFO

end module get_mbands_module
