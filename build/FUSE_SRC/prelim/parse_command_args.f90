module parse_command_args_MODULE

  USE nrtype
  USE data_types, only: cli_options

  implicit none

  private
  public :: parse_command_args

contains

  subroutine parse_command_args(opts, err, message)
  implicit none
  ! dummies
  type(cli_options)   , intent(out)   :: opts          ! command line interface options
  integer(i4b)        , intent(out)   :: err           ! error code
  character(len=1024) , intent(out)   :: message       ! error message
  ! internal
  integer(i4b)                        :: i             ! index of command line argument
  character(len=:)    , allocatable   :: a, v          ! command line arguments
  character(len=:)    , allocatable   :: cIndex        ! character index
  integer(i4b)                        :: nArg          ! number of command line arguments
  character(len=:)    , allocatable   :: cmessage
  ! initialize error control
  err=0; message='parse_command_args/'
  
  ! ----- parse command line arguments ------------------------------------------------------
  
  ! -----------------------------------------------------------------------------------------
  !   CLI parsing for FUSE run modes
  !    -c/--control   <file>     (required unless --version)
  !    -m/--runmode   <def|idx|opt|sce> (required unless --version)
  !    -d/--domid     <string>   (required unless --version)
  !    -s/--sets      <file>     (required for idx,opt)
  !    -i/--index     <int>      (required for idx)
  !    -r/--restart   <y|m|d|e|never>   (optional)
  !    -p/--progress  <m|d|h|never>     (optional)
  !    -t/--tag       <string>          (optional)
  !    -v/--version              (prints version info and exits)
  !    -h/--help                 (prints help and exits)
  ! -----------------------------------------------------------------------------------------

  nArg = command_argument_count()
  if (nArg < 1) call printCommandHelp()

  i = 1
  do while (i <= narg)
    call get_arg(i,a)

    select case (trim(a))

    case ('-h','--help')
      opts%show_help = .true.
      i = i + 1

    case ('-v','--version')
      opts%show_version = .true.
      i = i + 1

     case ('-t','--tag')
      call require_next(i, narg, a, v, err, cmessage)
      opts%tag = trim(v)
      i = i + 2

    case ('-c','--control')
      call require_next(i, narg, a, v, err, cmessage)
      opts%control_file = trim(v)
      i = i + 2

    case ('-m','--runmode')
      call require_next(i, narg, a, v, err, cmessage)
      opts%runmode = to_lower(trim(v))
      i = i + 2

    case ('-d','--domid')
      call require_next(i, narg, a, v, err, cmessage)
      opts%domain_id = trim(v)
      i = i + 2

    case ('-s','--sets','--param-sets')
      call require_next(i, narg, a, v, err, cmessage)
      opts%sets_file = trim(v)
      i = i + 2

    case ('-i','--index')
      call require_next(i, narg, a, cIndex, err, cmessage)
      i = i + 2

    case ('-r','--restart')
      call require_next(i, narg, a, v, err, cmessage)
      opts%restart_freq = to_lower(trim(v))
      i = i + 2

    case ('-p','--progress')
      call require_next(i, narg, a, v, err, cmessage)
      opts%progress_freq = to_lower(trim(v))
      i = i + 2

    case default
      if (len_trim(a) > 0 .and. a(1:1) == '-') then
        err = 1
        cmessage = "unknown option: "//trim(a)//"; type 'fuse.exe --help' for usage"
      else
        err = 1
        cmessage = "unexpected positional argument: "//trim(a)//"; type 'fuse.exe --help' for usage"
      end if
    end select

    ! process error code
    if(err/=0)then
     message=trim(message)//trim(cmessage)
     err=20; return
    endif

  end do  ! looping through arguments

  ! Early exits
  if (opts%show_help) then
    call printCommandHelp()
    stop 0
  end if
  if (opts%show_version) then
    call printVersionInfo()
    stop 0
  end if

  ! Parse parameter index
  if(allocated(cIndex))then
   call parse_int(cIndex, opts%indx, err, cmessage)
   if(err/=0)then
     message=trim(message)//trim(cmessage)
     err=20; return
    endif
  endif

  ! Validate required args
  if (.not. allocated(opts%control_file)) then
    err = 1; message = trim(message)//"missing required --control; type 'fuse.exe --help' for usage"; return
  end if
  if (.not. allocated(opts%domain_id)) then
    err = 1; message = trim(message)//"missing required --domid;   type 'fuse.exe --help' for usage"; return
  end if
  if (.not. allocated(opts%runmode)) then
    err = 1; message = trim(message)//"missing required --runmode; type 'fuse.exe --help' for usage"; return
  end if

  if (.not. is_valid_mode(opts%runmode)) then
    err = 1; message = trim(message)//"invalid --runmode: "//trim(opts%runmode)//" (expect def|idx|opt|sce)"; return
  end if

  ! Mode-dependent requirements
  select case (trim(opts%runmode))
  case ('idx')
    if (.not. allocated(opts%sets_file)) then
      err = 1; message = trim(message)//"runmode idx requires --sets <file>"; return
    end if
    if (opts%indx < 0) then
      err = 1; message = trim(message)//"runmode idx requires --index <int>"; return
    end if
  case ('opt')
    if (.not. allocated(opts%sets_file)) then
      err = 1; message = trim(message)//"runmode opt requires --sets <file>"; return
    end if
  case ('def','sce')
    ! no extra requirements
  end select

  ! Validate frequencies if provided (optional)
  if (allocated(opts%restart_freq)) then
    if (.not. is_valid_restart(opts%restart_freq)) then
      err = 1; message = trim(message)//"invalid --restart: "//trim(opts%restart_freq)//" (expect y|m|d|e|never)"; return
    end if
  end if
  if (allocated(opts%progress_freq)) then
    if (.not. is_valid_progress(opts%progress_freq)) then
      err = 1; message = trim(message)//"invalid --progress: "//trim(opts%progress_freq)//" (expect m|d|h|never)"; return
    end if
  end if

  end subroutine parse_command_args

  ! ----- list version ----------------------------------------------------------------------
  
  subroutine printVersionInfo()
    ! Assumes these are available, e.g. from:
    !   include "fuseversion.inc"
    ! somewhere in a used module (e.g., globaldata) OR add that include here.
    use globaldata, only: FUSE_VERSION, FUSE_BUILDTIME, FUSE_GITBRANCH, FUSE_GITHASH
    implicit none
    print '(A)', repeat('-', 70)
    print '(A)', 'FUSE'
    print '("  ",A12," : ",A)', 'Version',    trim(FUSE_VERSION)
    print '("  ",A12," : ",A)', 'Build time', trim(FUSE_BUILDTIME)
    print '("  ",A12," : ",A)', 'Git branch', trim(FUSE_GITBRANCH)
    print '("  ",A12," : ",A)', 'Git hash',   trim(FUSE_GITHASH)
    print '(A)', repeat('-', 70)
  end subroutine printVersionInfo

  ! ----- list command usage ----------------------------------------------------------------

  subroutine printCommandHelp()
    implicit none
    print "(A)", ""
    print "(A)", "Usage:"
    print "(A)", "  fuse.exe -d domain_id -c control_file -m {def|idx|opt|sce} [options]"
    print "(A)", ""
    
    print "(A)", "Run modes:"
    print "(A)", "  def : run with default parameter sets"
    print "(A)", "  idx : run using a given index from a parameter sets file"
    print "(A)", "  opt : run using best simulation from a parameter sets file"
    print "(A)", "  sce : optimize (SCE)"
    print "(A)", ""
    
    print "(A)", "Required:"
    print "(A)", "  -d, --domid        <string>   Domain ID"
    print "(A)", "  -c, --control      <file>     Control file"
    print "(A)", "  -m, --runmode      <mode>     def|idx|opt|sce"
    print "(A)", ""
    
    print "(A)", "Conditional:"
    print "(A)", "  -s, --sets         <file>   Parameter sets file (required for idx,opt)"
    print "(A)", "  -i, --index        <int>    Index (required for idx)"
    print "(A)", ""
    
    print "(A)", "Optional:"
    print "(A)", "  -r, --restart      <freq>   y|m|d|e|never"
    print "(A)", "  -p, --progress     <freq>   m|d|h|never"
    print "(A)", "  -t, --tag          <string> Add tag to output filename"
    print "(A)", "  -v, --version               Print version info and exit"
    print "(A)", "  -h, --help                  Print this help and exit"
    print "(A)", ""
    
    print "(A)", "Examples:"
    print "(A)", "  Default run (no parameter-set file):"
    print "(A)", "  fuse.exe -d camels-12345 -c ./control/FUSE_control.txt -m def"
    print "(A)", ""
    
    print "(A)", "  Default run with restart and progress output:"
    print "(A)", "  fuse.exe -d camels-12345 -c ./control/FUSE_control.txt -m def -r d -p h"
    print "(A)", ""
    
    print "(A)", "  Run using parameter set index 17 from a sets file:"
    print "(A)", "  fuse.exe -d camels-12345 -c ./control/FUSE_control.txt -m idx -s ./params/sets.nc -i 17"
    print "(A)", ""
    
    print "(A)", "  Run using the best simulation from a sets file:"
    print "(A)", "  fuse.exe -d camels-12345 -c ./control/FUSE_control.txt -m opt -s ./params/sets.nc"
    print "(A)", ""
    
    print "(A)", "  Optimize using SCE:"
    print "(A)", "  fuse.exe -d camels-12345 -c ./control/FUSE_control.txt -m sce"
    print "(A)", ""
    
    print "(A)", "  Print version information:"
    print "(A)", "  fuse.exe --version"
    print "(A)", ""
  end subroutine printCommandHelp

  ! -----------------------------------------------------------------------------------------
  ! Helpers
  ! -----------------------------------------------------------------------------------------

  subroutine get_arg(i, out)
    integer, intent(in) :: i
    character(len=:), allocatable, intent(out) :: out
    integer :: L
    call get_command_argument(i, length=L)
    allocate(character(len=L) :: out)
    call get_command_argument(i, out)
  end subroutine get_arg

  subroutine require_next(i, narg, opt, val, err, message)
    integer, intent(in) :: i, narg
    character(len=*), intent(in) :: opt
    character(len=:), allocatable, intent(out) :: val
    integer, intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    err = 0
    message = ""
    if (i+1 > narg) then
      err = 1
      message = "missing value after "//trim(opt)//"; type 'fuse.exe --help' for usage"
      return
    end if
    call get_arg(i+1, val)
  end subroutine require_next

  subroutine parse_int(s, x, err, message)
    character(len=*), intent(in) :: s
    integer, intent(out) :: x
    integer, intent(out) :: err
    character(len=:), allocatable, intent(out) :: message
    integer :: ios
    err = 0
    message = ""
    read(s, *, iostat=ios) x
    if (ios /= 0) then
      err = 1
      message = "invalid integer: "//trim(s)
    end if
  end subroutine parse_int

  pure function to_lower(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: k, c
    t = s
    do k = 1, len(s)
      c = iachar(t(k:k))
      if (c >= iachar('A') .and. c <= iachar('Z')) then
        t(k:k) = achar(c + (iachar('a') - iachar('A')))
      end if
    end do
  end function to_lower

  pure logical function is_valid_mode(m)
    character(len=*), intent(in) :: m
    is_valid_mode = (trim(m) == 'def' .or. trim(m) == 'idx' .or. trim(m) == 'opt' .or. trim(m) == 'sce')
  end function is_valid_mode

  pure logical function is_valid_restart(f)
    character(len=*), intent(in) :: f
    is_valid_restart = (trim(f) == 'y' .or. trim(f) == 'm' .or. trim(f) == 'd' .or. trim(f) == 'e' .or. trim(f) == 'never')
  end function is_valid_restart

  pure logical function is_valid_progress(f)
    character(len=*), intent(in) :: f
    is_valid_progress = (trim(f) == 'm' .or. trim(f) == 'd' .or. trim(f) == 'h' .or. trim(f) == 'never')
  end function is_valid_progress

end module parse_command_args_MODULE



