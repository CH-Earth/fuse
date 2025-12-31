module alloc_scratch_module


  USE nrtype
  use data_types, only: domain_info, fuse_work

  implicit none
  private
  public :: init_fuse_work

CONTAINS

  subroutine init_fuse_work(info, work, ierr, message)

    use globaldata, only: NPAR_SNOW
    implicit none

    type(domain_info), intent(in)    :: info
    type(fuse_work),   intent(inout) :: work
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    integer(i4b) :: ib
    integer(i4b) :: nBands, nState, nPar

    ierr=0; message="init_fuse_work/"

    ! identify dimensions
    nBands = info%snow%n_bands
    nState = info%config%nState
    nPar   = info%config%nParam

    ! If already initialized, don't reallocate unless sizes mismatch
    if (work%is_initialized) then
      if (size(work%state0)==nState .and. size(work%state1)==nState) return
      call free_fuse_work(work, ierr, message)
      if(ierr/=0) return
    endif

    ! ---- allocate core state vectors ----
    allocate(work%state0(nState), work%state1(nState), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate state0/state1"
      return
    endif

    ! optional debug scratch
    ! allocate(work%dSdt(nState), work%J(nState,nState), stat=ierr)

    ! ---- allocate differentiable parent derivatives ----
    allocate(work%fuseStruct%df_dS(nState), &
             work%fuseStruct%df_dPar(nPar), &
             work%fuseStruct%dL_dPar(nPar), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate fuseStruct derivatives"
      return
    endif

    ! ---- allocate elevation band containers ----
    allocate(work%fuseStruct%sbands(nBands), stat=ierr)
    if(ierr/=0) then
      message=trim(message)//"cannot allocate fuseStruct sbands"
      return
    endif

    ! ---- allocate per-band parameter derivative vectors ----
    do ib=1,nBands
      allocate(work%fuseStruct%sbands(ib)%var%dSWE_dParam(nPar_snow), stat=ierr)
      if(ierr/=0) then
        message=trim(message)//"cannot allocate dSWE_dParam for band"
        return
      endif
      work%fuseStruct%sbands(ib)%var%dSWE_dParam(:) = 0._sp
    enddo

    ! ---- initialize the band snow vars once ----
    work%fuseStruct%sbands(:)%var%SWE         = 0._sp
    work%fuseStruct%sbands(:)%var%SNOWACCMLTN = 0._sp
    work%fuseStruct%sbands(:)%var%SNOWMELT    = 0._sp
    work%fuseStruct%sbands(:)%var%DSWE_DT     = 0._sp

    work%is_initialized = .true.

  end subroutine init_fuse_work

  ! -------------------------------------------------------------------------------------

  subroutine free_fuse_work(work, ierr, message)

    implicit none
    type(fuse_work), intent(inout) :: work
    integer(i4b),    intent(out)   :: ierr
    character(*),    intent(out)   :: message

    integer(i4b) :: ib, istat

    ierr    = 0
    message = "free_fuse_work/"

    ! ---- state vectors ----
    if (allocated(work%state0)) then
      deallocate(work%state0, stat=istat)
      call note_fail("state0", istat)
    endif

    if (allocated(work%state1)) then
      deallocate(work%state1, stat=istat)
      call note_fail("state1", istat)
    endif

    ! ---- derivative arrays ----
    if (allocated(work%fuseStruct%df_dS)) then
      deallocate(work%fuseStruct%df_dS, stat=istat)
      call note_fail("fuseStruct%df_dS", istat)
    endif

    if (allocated(work%fuseStruct%df_dPar)) then
      deallocate(work%fuseStruct%df_dPar, stat=istat)
      call note_fail("fuseStruct%df_dPar", istat)
    endif

    if (allocated(work%fuseStruct%dL_dPar)) then
      deallocate(work%fuseStruct%dL_dPar, stat=istat)
      call note_fail("fuseStruct%dL_dPar", istat)
    endif

    ! ---- elevation band structures ----
    if (allocated(work%fuseStruct%sbands)) then

      do ib = 1, size(work%fuseStruct%sbands)
        if (allocated(work%fuseStruct%sbands(ib)%var%dSWE_dParam)) then
          deallocate(work%fuseStruct%sbands(ib)%var%dSWE_dParam, stat=istat)
          call note_fail("sbands%var%dSWE_dParam", istat)
        endif
      enddo

      deallocate(work%fuseStruct%sbands, stat=istat)
      call note_fail("fuseStruct%sbands", istat)

    endif

    work%is_initialized = .false.

    contains

      subroutine note_fail(where, istat)
        character(*), intent(in) :: where
        integer(i4b), intent(in) :: istat

        if (istat /= 0) then
          ! preserve the first nonzero stat as ierr
          if (ierr == 0) ierr = istat

          ! append context (do not overwrite)
          message = trim(message)//" dealloc_fail("//trim(where)//")"
        endif
      end subroutine note_fail

  end subroutine free_fuse_work

end module alloc_scratch_module
