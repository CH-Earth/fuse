module read_hydromet_module

use nrtype
use info_types, only: fuse_info
use domain_types, only: domain_data

use netcdf
use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

use fuse_globaldata, only: NA_VALUE
use fuse_globaldata, only: NVAR_HYDROMET
use fuse_globaldata, only: iPRECIP, iTEMP, iPET, iQOBS

implicit none

private

public :: read_met_data
public :: read_qobs_data
public :: read_hydromet_metadata

contains

  subroutine read_hydromet_metadata(ncid, info, ierr, message)
  implicit none
  integer(i4b), intent(in)       :: ncid
  type(fuse_info), intent(inout) :: info
  integer(i4b), intent(out)      :: ierr
  character(*), intent(out)      :: message

  integer(i4b) :: ivar
  integer(i4b) :: ierr_att
  character(len=256) :: cmessage

  ierr = 0
  message = "read_hydromet_metadata/"

  ! get table of name/varid pairs (names set in TOML read)
  info%files%hydromet%name(iPRECIP) = info%files%precip_name
  info%files%hydromet%name(iTEMP)   = info%files%temp_name
  info%files%hydromet%name(iPET)    = info%files%pet_name
  info%files%hydromet%name(iQOBS)   = info%files%qobs_name

  info%files%hydromet%varid(:) = -1

  ! get varid for each forcing variable
  do ivar = 1, NVAR_HYDROMET

    ! get the variable id
    ierr = nf90_inq_varid(ncid, trim(info%files%hydromet%name(ivar)), info%files%hydromet%varid(ivar))
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//' ['//trim(info%files%hydromet%name(ivar))//']'; return; endif

    ! get rank of variable in NetCDF file
    ierr = nf90_inquire_variable(ncid, info%files%hydromet%varid(ivar), ndims=info%files%hydromet%ndims(ivar) )
    if(ierr/=0)then; message=trim(message)//trim(nf90_strerror(ierr))//' ['//trim(info%files%hydromet%name(ivar))//']'; return; endif

    ! get units
    ierr = nf90_get_att(ncid, info%files%hydromet%varid(ivar), "units", &
                              info%files%hydromet%units(ivar) )
    
    if (ierr /= nf90_noerr) then
      message = trim(message)// "cannot read units for variable '"// &
                trim(info%files%hydromet%name(ivar))//"': "//            &
                trim(nf90_strerror(ierr))
      return
    end if

    ! get missing value metadata
    ierr_att = nf90_get_att(ncid, info%files%hydromet%varid(ivar), "_FillValue", &
                            info%files%hydromet%fill_value(ivar))

    if (ierr_att == nf90_noerr) then

      info%files%hydromet%has_fill_value(ivar) = .true.

    else if (ierr_att == nf90_enotatt) then

      ! Backward compatibility with NetCDF files using missing_value
      ierr_att = nf90_get_att(ncid, info%files%hydromet%varid(ivar), "missing_value", &
                              info%files%hydromet%fill_value(ivar))

      if (ierr_att == nf90_noerr) then
        info%files%hydromet%has_fill_value(ivar) = .true.

      else if (ierr_att /= nf90_enotatt) then
        ierr = ierr_att
        message = trim(message)//"cannot read missing_value for variable '"// &
                  trim(info%files%hydromet%name(ivar))//"': "//             &
                  trim(nf90_strerror(ierr))
        return
      end if

    else

      ierr = ierr_att
      message = trim(message)//"cannot read _FillValue for variable '"// &
                trim(info%files%hydromet%name(ivar))//"': "//           &
                trim(nf90_strerror(ierr))
      return

    end if

    ! get unit conversion
    select case (ivar)
    
      case (iTEMP) ! do nothing
    
      case (iPRECIP, iPET, iQOBS)
     
        call get_input_flux_multiplier(info%files%hydromet%units(ivar),      &
                                       info%files%hydromet%multiplier(ivar), &
                                       ierr, cmessage)
        if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

      case default
        message = trim(message)//"unexpected variable while reading units"
        ierr = 20; return
    
    end select

  end do  ! ivar

  end subroutine read_hydromet_metadata

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  
  SUBROUTINE read_met_data(info, itim_start, numtim, &
                           domain, ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Nans Addor, based on Martyn Clark's get_gforce
  ! Modified by Martyn Clark to simplify and update to new data structures, 02/2026
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Read NetCDF gridded forcing data for a range of time steps
  ! ---------------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input
  type(fuse_info),   intent(in)           :: info        ! info data structure that holds spatial indices
  integer(i4b),      intent(in)           :: itim_start  ! index of model time step - start of the period to extract
  integer(i4b),      intent(in)           :: numtim      ! number of model time steps to extract

  ! output
  type(domain_data), intent(inout)        :: domain      ! domain data structure that holds 3-d arrays
  integer(i4b),      intent(out)          :: ierr        ! error code
  character(*),      intent(out)          :: message     ! error message

  ! internal
  integer(i4b)                            :: iVar        ! loop through hydromet data
  integer(i4b), parameter                 :: ndim_2d=2   ! named variable for 2 dimensions
  integer(i4b), parameter                 :: ndim_3d=3   ! named variable for 3 dimensions
  real(wp),dimension(:,:,:),allocatable   :: gTemp       ! temporary 3d grid
  integer(i4b)                            :: nx, ny      ! grid dimensions
  integer(i4b)                            :: ystart      ! start index iin input file (for MPI)
  integer(i4b), allocatable               :: nc_start(:) ! start indices in NetCDF file
  integer(i4b), allocatable               :: nc_count(:) ! count in NetCDF file

  ! initialize error control
  ierr=0; message='read_met_data/'
  ! ---------------------------------------------------------------------------------------

  ! get indices in the input file
  nx     = info%space%nx_local
  ny     = info%space%ny_local
  ystart = info%space%y_start_global
  
  ! allocate space for gridded forcing buffer
  allocate(gtemp(nx,ny,numtim), stat=ierr)
  if (ierr /= 0) then
    message = trim(message)//'unable to allocate hydromet buffer'
    return
  end if

  ! loop through hydromet variables
  do ivar = 1, NVAR_HYDROMET

    ! just meteorological forcing data here
    if (ivar == iQOBS) cycle

    ! separate read for 2d and 3d input format
    select case ( info%files%hydromet%ndims(ivar) )

      ! 2-d catchments (hru, time) -> (1,nSpat2,time)
      case (ndim_2d)

        nc_start = (/ystart, itim_start/)
        nc_count = (/    ny,     numtim/)

        ierr = nf90_get_var(info%files%ncid_hydromet,        &
                            info%files%hydromet%varid(ivar), &
                            gtemp(1,:,1:numtim),             &
                            start=nc_start, count=nc_count)
     
      ! 3-d grid (x,y,t)                    
      case (ndim_3d)
        
        nc_start = (/ 1, ystart, itim_start/)
        nc_count = (/nx,     ny,     numtim/)

        ierr = nf90_get_var(info%files%ncid_hydromet,        &
                            info%files%hydromet%varid(ivar), &
                            gtemp(:,:,1:numtim),             &
                            start=nc_start, count=nc_count)
    
      case default
        message=trim(message)//'unknown dimensions for variable '//trim( info%files%hydromet%name(ivar) )
        ierr=20; return

    end select

    if (ierr /= nf90_noerr) then
      message = trim(message)//trim(nf90_strerror(ierr))
      return
    end if

    ! save the data in the structure -- and convert fluxes to mm/day
    select case(ivar)
  
      case (iPRECIP);  domain%force(:,:,1:numtim)%ppt  = gTemp(:,:,:) * info%files%hydromet%multiplier(iPRECIP)
      case (iTEMP);    domain%force(:,:,1:numtim)%temp = gTemp(:,:,:)
      case (iPET);     domain%force(:,:,1:numtim)%pet  = gTemp(:,:,:) * info%files%hydromet%multiplier(iPET)
      case default
        message=trim(message)//'unable to identify forcing variable'
        ierr=10; return
  
    end select  ! identify forcing variable

  end do  ! (loop thru forcing variables)

  end subroutine read_met_data

  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------

  SUBROUTINE read_qobs_data(info, itim_start, numtim, &
                            qobs, ierr, message)
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Read NetCDF streamflow for a range of time steps
  ! ---------------------------------------------------------------------------------------

  IMPLICIT NONE

  ! input
  type(fuse_info),   intent(in)           :: info        ! info data structure that holds spatial indices
  integer(i4b),      intent(in)           :: itim_start  ! index of model time step - start of the period to extract
  integer(i4b),      intent(in)           :: numtim      ! number of model time steps to extract

  ! output
  real(wp),          intent(out)          :: qobs(:,:)   ! data array that holds qobs (nq, nt)
  integer(i4b),      intent(out)          :: ierr        ! error code
  character(*),      intent(out)          :: message     ! error message

  ! internal
  integer(i4b)                            :: nc_start(2) ! start indices in NetCDF file
  integer(i4b)                            :: nc_count(2) ! count in NetCDF file

  ! initialize error control
  ierr=0; message='read_qobs_data/'
  ! ---------------------------------------------------------------------------------------

  nc_start = (/              1, itim_start/)
  nc_count = (/info%space%nobs,     numtim/)

  ierr = nf90_get_var(info%files%ncid_hydromet,         &
                      info%files%hydromet%varid(iQOBS), &
                      qobs(:,1:numtim),                 &
                      start=nc_start, count=nc_count)

  if (ierr /= nf90_noerr) then
    message = trim(message)//trim(nf90_strerror(ierr))
    return
  end if

  ! Normalize missing observations before unit conversion
  where (ieee_is_nan(qobs(:,1:numtim)))
    qobs(:,1:numtim) = real(NA_VALUE, wp)
  end where

  if (info%files%hydromet%has_fill_value(iQOBS)) then
    where (qobs(:,1:numtim) == info%files%hydromet%fill_value(iQOBS))
      qobs(:,1:numtim) = real(NA_VALUE, wp)
    end where
  end if

  ! Convert valid streamflow observations to FUSE internal units
  where (qobs(:,1:numtim) /= real(NA_VALUE, wp))
    qobs(:,1:numtim) = qobs(:,1:numtim) * &
                       info%files%hydromet%multiplier(iQOBS)
  end where

  end subroutine read_qobs_data

  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------

  subroutine get_input_flux_multiplier(cunits, amult, ierr, message)

   character(len=*), intent(in)    :: cunits
   real(wp),         intent(out)   :: amult
   integer(i4b),     intent(inout) :: ierr
   character(*),     intent(inout) :: message

   character(len=128) :: units_lc
   character(len=32)  :: length_unit
   character(len=32)  :: time_unit
   real(wp)           :: length_to_mm
   real(wp)           :: units_per_day
   integer(i4b)       :: ipos
   integer(i4b)       :: i
   integer(i4b)       :: ichar_value

   integer(i4b), parameter :: syntax_unknown = 0
   integer(i4b), parameter :: syntax_slash   = 1
   integer(i4b), parameter :: syntax_per     = 2
   integer(i4b), parameter :: syntax_inverse = 3

   integer(i4b) :: syntax_type
   integer(i4b) :: delimiter_len
   integer(i4b) :: ispace

   character(len=:), allocatable :: delimiter

   ierr  = 0
   amult = -1._wp

   units_lc = adjustl(trim(cunits))

   ! Convert ASCII upper case to lower case
   do i = 1, len_trim(units_lc)
     ichar_value = iachar(units_lc(i:i))

     if (ichar_value >= iachar("A") .and. &
         ichar_value <= iachar("Z")) then
       units_lc(i:i) = achar(ichar_value + 32)
     end if
   end do

   ! Current supported syntax: length/time, length per time and length time-1
   syntax_type  = syntax_unknown
   delimiter_len = 0
   ipos = 0

   ! Identify the unit syntax.
   if (index(units_lc, " per ") > 0) then
     syntax_type = syntax_per

   else if (index(units_lc, "/") > 0) then
     syntax_type = syntax_slash

   else if (index(units_lc, "-1") > 0) then
     syntax_type = syntax_inverse
   end if

   select case (syntax_type)

     case (syntax_per)

       delimiter = " per "
       delimiter_len = len(delimiter)
       ipos = index(units_lc, delimiter)

       if (ipos <= 1 .or. &
           ipos + delimiter_len > len_trim(units_lc)) then
         ierr = 20
         message = trim(message)// &
           "unsupported flux units '"//trim(cunits)//"': empty length or time unit"
         return
       end if

       length_unit = adjustl(trim(units_lc(:ipos-1)))
       time_unit = adjustl(trim( &
         units_lc(ipos+delimiter_len:len_trim(units_lc)) ))

     case (syntax_slash)

       delimiter = "/"
       delimiter_len = len(delimiter)
       ipos = index(units_lc, delimiter)

       if (ipos <= 1 .or. &
           ipos + delimiter_len > len_trim(units_lc)) then
         ierr = 20
         message = trim(message)//"unsupported flux units '"//trim(cunits)// "': empty length or time unit"
         return
       end if

       length_unit = adjustl(trim(units_lc(:ipos-1)))
       time_unit = adjustl(trim( &
         units_lc(ipos+delimiter_len:len_trim(units_lc)) ))

     case (syntax_inverse)

       ! Expected syntax: length time-1
       ispace = index(trim(units_lc), " ")

       if (ispace <= 1) then
         ierr = 20
         message = trim(message)// &
           "unsupported flux units '"//trim(cunits)//"': expected length time-1"
         return
       end if

       length_unit = adjustl(trim(units_lc(:ispace-1)))
       time_unit = adjustl(trim(units_lc(ispace+1:)))

       if (len_trim(time_unit) <= 2) then
         ierr = 20
         message = trim(message)// "unsupported inverse-time unit '"//trim(time_unit)//"' in '"//trim(cunits)//"'"
         return
       end if

       if (time_unit(len_trim(time_unit)-1:len_trim(time_unit)) /= "-1") then
         ierr = 20
         message = trim(message)// &
           "unsupported inverse-time unit '"//trim(time_unit)// "' in '"//trim(cunits)//"'"
         return
       end if

       ! Remove the trailing "-1".
       time_unit = trim(time_unit(:len_trim(time_unit)-2))

     case default

       ierr = 20
       message = trim(message)// "unsupported flux units '"//trim(cunits)// "': expected length/time, length per time, or length time-1"
       return

   end select

   ! Convert the numerator to millimetres
   select case (trim(length_unit))

     case ("mm", "millimeter", "millimeters", "millimetre", "millimetres")
       length_to_mm = 1._wp

     case ("m", "meter", "meters", "metre", "metres")
       length_to_mm = 1000._wp

     case default
       ierr = 20
       message = trim(message)// "unsupported length unit '"//trim(length_unit)// "' in '"//trim(cunits)//"'. Supported length units are: mm or m."
       return

   end select

   ! Convert the denominator to units per day
   select case (trim(time_unit))

     case ("s", "sec", "secs", "second", "seconds")
       units_per_day = 86400._wp

     case ("min", "mins", "minute", "minutes")
       units_per_day = 1440._wp

     case ("h", "hr", "hrs", "hour", "hours")
       units_per_day = 24._wp

     case ("d", "day", "days")
       units_per_day = 1._wp

     case default
       ierr = 20
       message = trim(message)// "unsupported time unit '"//trim(time_unit)// "' in '"//trim(cunits)//"'. Supported time units are: s, min, h, d."
       return

   end select

   amult = length_to_mm * units_per_day

  end subroutine get_input_flux_multiplier

end module read_hydromet_module
