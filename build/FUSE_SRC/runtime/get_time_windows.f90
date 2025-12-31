module time_windows_module

  use nrtype
  use netcdf
  use data_types, only: domain_type
  use fuse_fileManager, only: date_start_sim, date_end_sim, date_start_eval, date_end_eval, numtim_sub_str
  use time_io, only: date_extractor, juldayss

  implicit none

  private
  public :: get_time_windows
  public :: export_time_to_multiforce

  contains

  subroutine get_time_windows(ncid, domain, ierr, message)
    
    integer(i4b),      intent(in)    :: ncid
    type(domain_type), intent(inout) :: domain
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    integer(i4b) :: nt
    real(sp), allocatable :: time_steps(:)
    character(len=1024) :: units_local
    integer(i4b) :: ios
    character(len=1024) :: cmessage
    
    ierr=0; message="get_time_windows/"

    ! ----- read forcing time axis ------------------------------------------------------

    call read_time_axis(ncid, time_steps, units_local, nt, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    domain%info%time%nt_global = nt
    domain%info%time%units     = trim(units_local)

    ! ----- build julian-day axis -------------------------------------------------------
    
    call build_julian_axis(time_steps, trim(units_local), domain%info%time%jdate_ref, domain%info%time%jdate, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- compute indices for sim/eval windows ----------------------------------------

    ! simulation indices
    call map_dates_to_indices(domain%info%time%jdate, date_start_sim, date_end_sim, &
                              domain%info%time%sim_beg, domain%info%time%sim_end, ierr, cmessage)
    if (ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! evaluation indices
    call map_dates_to_indices(domain%info%time%jdate, date_start_eval, date_end_eval, &
                              domain%info%time%eval_beg, domain%info%time%eval_end, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- validate window consistency -------------------------------------------------
    
    call validate_windows(domain%info%time, ierr, cmessage)
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- derive simulation length ----------------------------------------------------
    
    domain%info%time%nt_sim = domain%info%time%sim_end - domain%info%time%sim_beg + 1

    ! ----- configure sub-period windowing ----------------------------------------------
    
    ! convert sub-period string to integer
    read(numtim_sub_str,*,iostat=ios) domain%info%time%nt_window
    if(ios/=0) then
      ierr=1; message=trim(message)//"cannot parse numtim_sub_str"; return
    endif

    ! handle cases where sub-periods are undefined
    if(domain%info%time%nt_window == -9999) then
      domain%info%time%use_subperiods = .false.
      domain%info%time%nt_window      = domain%info%time%nt_sim
    else
      domain%info%time%use_subperiods = .true.
      ! keep nt_window as user-chosen chunk size
    endif

    ! ----- validate time-window configuration (subperiods allowed only in grid mode) ---
    if( (.not. domain%info%space%grid_flag) .and. domain%info%time%use_subperiods ) then
      ierr = 1
      message = trim(message)// &
                "catchment mode requires running the full time series in one chunk; " // &
                "set numtim_sub = -9999 in the filemanager."
      return
    endif

    ! ----- finalize --------------------------------------------------------------------

    if(allocated(time_steps)) deallocate(time_steps)

  end subroutine get_time_windows

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- backwards compatibility: export to multiforce globals -------------------------
 
  ! - New code stores all time-window metadata in domain%info%time (source of truth).
  ! - Legacy routines still read multiforce globals (sim_beg, sim_end, numtim_sub, ...).

  subroutine export_time_to_multiforce(domain)
    use multiforce, only: sim_beg, sim_end, eval_beg, eval_end, numtim_sim, numtim_sub, &
                          SUB_PERIODS_FLAG, istart
    implicit none
    type(domain_type), intent(in) :: domain
   
    sim_beg = domain%info%time%sim_beg
    sim_end = domain%info%time%sim_end
    eval_beg = domain%info%time%eval_beg
    eval_end = domain%info%time%eval_end
   
    numtim_sim = domain%info%time%nt_sim
    numtim_sub = domain%info%time%nt_window
    SUB_PERIODS_FLAG = domain%info%time%use_subperiods
   
    istart = sim_beg
  end subroutine

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! ----- helper routines ---------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: read time axis from NetCDF --------------------------------------------

  subroutine read_time_axis(ncid, time_steps, units, nt, ierr, message)
    
    integer(i4b), intent(in) :: ncid
    real(sp), allocatable, intent(out) :: time_steps(:)
    character(len=*), intent(out) :: units
    integer(i4b), intent(out) :: nt, ierr
    character(*), intent(out) :: message

    integer(i4b) :: varid, dimids(1)

    ierr=0; message="read_time_axis/"

    ierr = nf90_inq_varid(ncid, "time", varid)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot find time variable"; return
    endif

    ierr = nf90_inquire_variable(ncid, varid, dimids=dimids)
    if(ierr/=nf90_noerr) then
      message=trim(message)//trim(nf90_strerror(ierr)); return
    endif

    ierr = nf90_inquire_dimension(ncid, dimids(1), len=nt)
    if(ierr/=nf90_noerr) then
      message=trim(message)//trim(nf90_strerror(ierr)); return
    endif

    allocate(time_steps(nt), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"allocate(time_steps) failed"; return
    endif

    ierr = nf90_get_var(ncid, varid, time_steps)
    if(ierr/=nf90_noerr) then
      message=trim(message)//trim(nf90_strerror(ierr)); return
    endif

    ierr = nf90_get_att(ncid, varid, "units", units)
    if(ierr/=nf90_noerr) then
      message=trim(message)//"cannot read time units attribute"; return
    endif

  end subroutine read_time_axis

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: build julian axis -----------------------------------------------------

  subroutine build_julian_axis(time_steps, units, jref, jdate, ierr, message)
    
    real(sp), intent(in) :: time_steps(:)
    character(len=*), intent(in) :: units
    real(sp), intent(out) :: jref
    real(sp), allocatable, intent(out) :: jdate(:)
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: iy,im,id,ih
    character(len=1024) :: cmessage
    real(sp) :: scale_to_days

    ierr=0; message="build_julian_axis/"

    ! extract reference date from the units string
    call date_extractor(trim(units), iy, im, id, ih)
    call juldayss(iy,im,id,ih, jref, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! determine scaling factor to convert time_steps into days
    scale_to_days = time_units_to_days(units, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! build julian axis
    allocate(jdate(size(time_steps)), stat=ierr)
    if(ierr/=0) then; message=trim(message)//"allocate(jdate) failed"; return; endif
    jdate = jref + time_steps * scale_to_days

  end subroutine build_julian_axis

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: determine scaling factor to convert time_steps into days --------------

  real(sp) function time_units_to_days(units, ierr, message)
    implicit none
    character(len=*), intent(in) :: units
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
   
    character(len=:), allocatable :: u
    integer(i4b) :: p
   
    ierr=0; message="time_units_to_days/"
   
    ! lower-case copy (simple approach)
    u = tolower_str( trim(adjustl(units)) )
   
    ! Look at the first token before a space
    p = index(u, " ")
    if(p <= 1) then
      ierr=1; message=trim(message)//"cannot parse units string: "//trim(units)
      time_units_to_days = 0._sp
      return
    endif
   
    select case (trim(u(1:p-1)))
      case ("days", "day")
        time_units_to_days = 1._sp
      case ("hours", "hour")
        time_units_to_days = 1._sp / 24._sp
      case ("minutes", "minute", "mins", "min")
        time_units_to_days = 1._sp / 1440._sp
      case ("seconds", "second", "secs", "sec")
        time_units_to_days = 1._sp / 86400._sp
      case default
        ierr=1
        message=trim(message)//"unsupported time unit: "//trim(u(1:p-1))
        time_units_to_days = 0._sp
    end select

  end function time_units_to_days

  pure function tolower_str(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i
    do i=1,len(s)
      select case(s(i:i))
        case("A":"Z"); out(i:i) = achar(iachar(s(i:i)) + 32)
        case default;  out(i:i) = s(i:i)
      end select
    end do
  end function tolower_str

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: map start/end date strings to indices ---------------------------------

  subroutine map_dates_to_indices(jdate, date_start, date_end, i_beg, i_end, ierr, message)

    real(sp), intent(in) :: jdate(:)
    character(len=*), intent(in) :: date_start, date_end
    integer(i4b), intent(out) :: i_beg, i_end
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: iy,im,id,ih
    real(sp) :: j_start, j_end
    character(len=1024) :: cmessage

    ierr=0; message="map_dates_to_indices/"

    ! start date
    call date_extractor(trim(date_start), iy,im,id,ih)
    call juldayss(iy,im,id,ih, j_start, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! end date
    call date_extractor(trim(date_end), iy,im,id,ih)
    call juldayss(iy,im,id,ih, j_end, ierr, cmessage)
    if(ierr/=0) then;  message=trim(message)//trim(cmessage); return; endif

    ! validate

    if(j_start > j_end) then
      ierr=1; message=trim(message)//"start date > end date"; return
    endif

    if(j_start < minval(jdate) .or. j_end > maxval(jdate)) then
      ierr=1; message=trim(message)//"requested window outside forcing range"; return
    endif

    ! get indices in jdate vector
    i_beg = minloc(abs(jdate - j_start), 1)
    i_end = minloc(abs(jdate - j_end  ), 1)

  end subroutine map_dates_to_indices


  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  ! ----- helper: validate sim/eval logic -----------------------------------------------

  subroutine validate_windows(ti, ierr, message)
 
    use data_types, only: time_info
    type(time_info), intent(in) :: ti
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    ierr=0; message="validate_windows/"

    if(ti%eval_beg < ti%sim_beg) then
      ierr=1; message=trim(message)//"eval start < sim start"; return
    endif
    if(ti%eval_end > ti%sim_end) then
      ierr=1; message=trim(message)//"eval end > sim end"; return
    endif

  end subroutine validate_windows

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

end module time_windows_module
