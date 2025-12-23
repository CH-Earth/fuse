MODULE GET_TIME_INDICES_MODULE

  USE nrtype
  use time_io
  implicit none

  contains
  SUBROUTINE GET_TIME_INDICES

    ! EXTRACT DATES AND DETERMINE ASSOCIATED INDICES

    ! convert start and end date of the NetCDF input file to julian day (Julian day is the continuous
    ! count of days since the beginning of the Julian Period around 4700 BC)

    USE multiforce, ONLY: timeUnits,time_steps,julian_day_input    ! time data
    USE multiforce, only: numtim_in, itim_in, istart               ! length of input time series and associated index
    USE multiforce, only: numtim_sim, itim_sim                ! length of simulated time series and associated index
    USE multiforce, only: numtim_sub, itim_sub                ! length of subperiod time series and associated index
    USE multiforce, only: sim_beg,sim_end                     ! timestep indices
    USE multiforce, only: eval_beg,eval_end                   ! timestep indices
    USE multiforce, only: SUB_PERIODS_FLAG                    ! .true. if subperiods are used to run FUSE

    USE fuse_fileManager,only:date_start_sim,date_end_sim,&
              date_start_eval,date_end_eval,&
              numtim_sub_str

    real(sp)                               :: jdate_ref_netcdf
    INTEGER(I4B)                           :: ERR             ! error code
    CHARACTER(LEN=1024)                    :: MESSAGE         ! error message

    ! dummies
    integer(i4b)                           :: iy,im,id,ih,imin  ! to temporarily store year, month, day, hour, min
    real(sp)                               :: isec              ! to temporarily store sec
    real(sp)                               :: jdate             ! to temporarily store a julian date
    real(sp)                               :: jdate_start_sim    ! date start simulation
    real(sp)                               :: jdate_end_sim      ! date end simulation
    real(sp)                               :: jdate_start_eval   ! date start evaluation period
    real(sp)                               :: jdate_end_eval     ! date end evaluation period

    ! ---------------------------------------------------------------------------------------
    ! process time data from foring file
    call date_extractor(trim(timeUnits),iy,im,id,ih) ! break down reference date of NetCDF file
    call juldayss(iy,im,id,ih,            &          ! convert it to julian day
                    jdate_ref_netcdf,err,message)

    julian_day_input=jdate_ref_netcdf+time_steps ! julian day of each time step of the input file

    call caldatss(julian_day_input(1),iy,im,id,ih,imin,isec)
    print *, 'Start date input file=',iy,im,id

    call caldatss(julian_day_input(numtim_in),iy,im,id,ih,imin,isec)
    print *, 'End date input file=  ',iy,im,id

    ! convert dates for simulation into julian day
    call date_extractor(trim(date_start_sim),iy,im,id,ih)        ! break down date
    call juldayss(iy,im,id,ih,jdate_start_sim,err,message)       ! convert it to julian day
    if(jdate_start_sim.lt.minval(julian_day_input))then          ! check forcing available
      call caldatss(jdate_start_sim,iy,im,id,ih,imin,isec)
       print *, 'Error: hydrologic simulation cannot start on ',iy,im,id,' because atmospheric forcing starts later (see above)';stop;
    endif
    sim_beg= minloc(abs(julian_day_input-jdate_start_sim),1)     ! find correponding index
    call caldatss(julian_day_input(sim_beg),iy,im,id,ih,imin,isec)
    print *, 'Date start sim=       ',iy,im,id

    call date_extractor(trim(date_end_sim),iy,im,id,ih)          ! break down date
    call juldayss(iy,im,id,ih,jdate_end_sim,err,message)         ! convert it to julian day
    if(jdate_end_sim.gt.maxval(julian_day_input))then         ! check forcing available
      call caldatss(jdate_end_sim,iy,im,id,ih,imin,isec)
       print *, 'Error: hydrologic simulation cannot end on ',iy,im,id,' because atmospheric forcing ends earlier (see above)';stop;
    endif
    sim_end= minloc(abs(julian_day_input-jdate_end_sim),1)      ! find correponding index
    call caldatss(julian_day_input(sim_end),iy,im,id,ih,imin,isec)
    print *, 'Date end sim=         ',iy,im,id

    call date_extractor(trim(date_start_eval),iy,im,id,ih)       ! break down date
    call juldayss(iy,im,id,ih,jdate_start_eval,err,message)      ! convert it to julian day
    eval_beg= minloc(abs(julian_day_input-jdate_start_eval),1)  ! find correponding index
    call caldatss(julian_day_input(eval_beg),iy,im,id,ih,imin,isec)
    print *, 'Date start eval=      ',iy,im,id

    call date_extractor(trim(date_end_eval),iy,im,id,ih)         ! break down date
    call juldayss(iy,im,id,ih,jdate_end_eval,err,message)        ! convert it to julian day
    eval_end= minloc(abs(julian_day_input-jdate_end_eval),1)    ! find correponding index
    call caldatss(julian_day_input(eval_end),iy,im,id,ih,imin,isec)
    print *, 'Date end eval=        ',iy,im,id

    ! check start before end
    if(jdate_start_sim.gt.jdate_end_sim)then; print *, 'Error: date_start_sim > date_end_sim '; stop; endif
    if(jdate_start_eval.gt.jdate_end_eval)then; print *, 'Error: date_start_eval > date_end_eval '; stop; endif

    ! check input data available for desired runs
    if(jdate_start_sim.lt.julian_day_input(1))then; print *, 'Error: date_start_sim is before the start if the input data'; stop; endif
    if(jdate_end_sim.gt.julian_day_input(numtim_in))then; print *, 'Error: the date_stop_sim is after the end of the input data'; stop; endif

    ! check input data available for desired runs
    if(jdate_start_eval.lt.jdate_start_sim)then; print *, 'Error: date_start_eval < date_start_sim'; stop; endif
    if(jdate_end_eval.gt.jdate_end_sim)then; print *, 'Error: date_end_eval > date_end_sim'; stop; endif

    ! determine length of simulations
    numtim_sim=sim_end-sim_beg+1
    istart=sim_beg

    ! determine length of subperiods
    read(numtim_sub_str,*,iostat=err) numtim_sub ! convert string to integer

    if(numtim_sub.eq.-9999)then

      print *, 'numtim_sub = -9999, FUSE will be run in 1 chunk of ',numtim_sim, 'time steps'
      SUB_PERIODS_FLAG=.FALSE.

      numtim_sub=numtim_sim ! no subperiods, run the whole time series

    else

      print *, 'FUSE will be run in chunks of ',numtim_sub, 'time steps'
      SUB_PERIODS_FLAG=.TRUE.

    end if

  END SUBROUTINE GET_TIME_INDICES
END MODULE
