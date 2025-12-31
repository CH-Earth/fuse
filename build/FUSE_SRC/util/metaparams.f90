MODULE metaparams
  
  ! ---------------------------------------------------------------------------------------
  ! Creator:
  ! --------
  ! Martyn Clark, 2007
  ! Modified by Brian Henn to include snow model, 6/2013
  ! Modified by Martyn Clark to avoid per-band parameters, 12/2025
  ! ---------------------------------------------------------------------------------------
  ! Purpose:
  ! --------
  ! Describe all parameters used in the model (used to define NetCDF output files, etc.)
  ! ---------------------------------------------------------------------------------------
  
  ! variable definitions
  USE nrtype
  
  IMPLICIT NONE

  private
  public :: PARDESCRIBE                   ! make subroutine public
  public :: PNAME, PDESC, PUNIT, isBand   ! make metadata variables public
  public :: NOUTPAR                       ! make number of output parameters public

  CHARACTER(LEN=11), DIMENSION(200)      :: PNAME       ! parameter names
  CHARACTER(LEN=52), DIMENSION(200)      :: PDESC       ! parameter long names (description of variable)
  CHARACTER(LEN= 8), DIMENSION(200)      :: PUNIT       ! parameter units
  logical(lgt)     , DIMENSION(200)      :: isBand      ! flag for the parameter dimension
  INTEGER(I4B)                           :: NOUTPAR     ! number of model parameters for output
  
  CONTAINS
  ! ---------------------------------------------------------------------------------------
  
  SUBROUTINE PARDESCRIBE()
  implicit none
  INTEGER(I4B)                           :: I           ! loop through parameter sets
  
  I=0  ! initialize counter
  
  ! adjustable model parameters
  I=I+1; PNAME(I)='RFERR_ADD  '; PDESC(I)='additive rainfall error                            '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='RFERR_MLT  '; PDESC(I)='multiplicative rainfall error                      '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXWATR_1  '; PDESC(I)='maximum total storage in the upper layer           '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXWATR_2  '; PDESC(I)='maximum total storage in the lower layer           '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='FRACTEN    '; PDESC(I)='fraction total storage as tension storage          '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='FRCHZNE    '; PDESC(I)='fraction tension storage in recharge zone          '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='FPRIMQB    '; PDESC(I)='fraction of baseflow in primary reservoir          '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='RTFRAC1    '; PDESC(I)='fraction of roots in the upper layer               '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='PERCRTE    '; PDESC(I)='percolation rate                                   '; PUNIT(I)='mm day-1'; isBand(i)=.false.
  I=I+1; PNAME(I)='PERCEXP    '; PDESC(I)='percolation exponent                               '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='SACPMLT    '; PDESC(I)='percolation multiplier in the SAC model            '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='SACPEXP    '; PDESC(I)='percolation exponent in the SAC model              '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='PERCFRAC   '; PDESC(I)='fraction of percolation to tension storage         '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='FRACLOWZ   '; PDESC(I)='fraction of soil excess to lower zone              '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='IFLWRTE    '; PDESC(I)='interflow rate                                     '; PUNIT(I)='mm day-1'; isBand(i)=.false.
  I=I+1; PNAME(I)='BASERTE    '; PDESC(I)='baseflow rate                                      '; PUNIT(I)='mm day-1'; isBand(i)=.false.
  I=I+1; PNAME(I)='QB_POWR    '; PDESC(I)='baseflow exponent                                  '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='QB_PRMS    '; PDESC(I)='baseflow depletion rate                            '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='QBRATE_2A  '; PDESC(I)='baseflow depletion rate for primary reservoir      '; PUNIT(I)='day-1   '; isBand(i)=.false.
  I=I+1; PNAME(I)='QBRATE_2B  '; PDESC(I)='baseflow depletion rate for secondary reservoir    '; PUNIT(I)='day-1   '; isBand(i)=.false.
  I=I+1; PNAME(I)='SAREAMAX   '; PDESC(I)='maximum saturated area                             '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='AXV_BEXP   '; PDESC(I)='ARNO/VIC b exponent                                '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='LOGLAMB    '; PDESC(I)='mean value of the log-transformed topographic index'; PUNIT(I)='log m   '; isBand(i)=.false.
  I=I+1; PNAME(I)='TISHAPE    '; PDESC(I)='shape parameter for the topo index Gamma distribtn '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='TIMEDELAY  '; PDESC(I)='time delay in runoff (routing)                     '; PUNIT(I)='day     '; isBand(i)=.false.
  I=I+1; PNAME(I)='MBASE      '; PDESC(I)='snow model base melt temperature                   '; PUNIT(I)='deg.C   '; isBand(i)=.false.
  I=I+1; PNAME(I)='MFMAX      '; PDESC(I)='snow model maximum melt factor                     '; PUNIT(I)='mm/(C-d)'; isBand(i)=.false.
  I=I+1; PNAME(I)='MFMIN      '; PDESC(I)='snow model minimum melt factor                     '; PUNIT(I)='mm/(C-d)'; isBand(i)=.false.
  I=I+1; PNAME(I)='PXTEMP     '; PDESC(I)='rain-snow partition temperature                    '; PUNIT(I)='deg.C   '; isBand(i)=.false.
  I=I+1; PNAME(I)='OPG        '; PDESC(I)='maximum relative precip difference across the bands'; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='LAPSE      '; PDESC(I)='maximum temperature difference across the bands    '; PUNIT(I)='deg.C   '; isBand(i)=.false.
  
  ! derived model parameters
  I=I+1; PNAME(I)='MAXTENS_1  '; PDESC(I)='maximum tension storage in the upper layer         '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXTENS_1A '; PDESC(I)='maximum storage in the recharge zone               '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXTENS_1B '; PDESC(I)='maximum storage in the lower zone                  '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXFREE_1  '; PDESC(I)='maximum free storage in the upper layer            '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXTENS_2  '; PDESC(I)='maximum tension storage in the lower layer         '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXFREE_2  '; PDESC(I)='maximum free storage in the lower layer            '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXFREE_2A '; PDESC(I)='maximum storage in the primary baseflow reservoir  '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXFREE_2B '; PDESC(I)='maximum storage in the secondary baseflow reservoir'; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='RTFRAC2    '; PDESC(I)='fraction of roots in the lower layer               '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='QBSAT      '; PDESC(I)='baseflow at saturation (derived parameter)         '; PUNIT(I)='mm day-1'; isBand(i)=.false.
  I=I+1; PNAME(I)='POWLAMB    '; PDESC(I)='mean value of power-transformed topographic index  '; PUNIT(I)='m**(1/n)'; isBand(i)=.false.
  I=I+1; PNAME(I)='MAXPOW     '; PDESC(I)='max value of power-transformed topographic index   '; PUNIT(I)='m**(1/n)'; isBand(i)=.false.
  
  ! model bands parameters 
  I=I+1; PNAME(I)='N_BANDS    '; PDESC(I)='number of basin bands in model                     '; PUNIT(I)='=       '; isBand(i)=.false.
  I=I+1; PNAME(I)='Z_FORCING  '; PDESC(I)='elevation of model forcing data                    '; PUNIT(I)='m       '; isBand(i)=.false.
  I=I+1; PNAME(I)='Z_MID      '; PDESC(I)='basin band mid-point elevation   (bands)           '; PUNIT(I)='m       '; isBand(i)=.true.
  I=I+1; PNAME(I)='AF         '; PDESC(I)='basin band area fraction         (bands)           '; PUNIT(I)='-       '; isBand(i)=.true.
  
  ! numerical solution parameters
  I=I+1; PNAME(I)='SOLUTION   '; PDESC(I)='0=explicit euler; 1=implicit euler                 '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='TIMSTEP_TYP'; PDESC(I)='0=fixed time steps; 1=adaptive time steps          '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='INITL_GUESS'; PDESC(I)='0=old state; 1=explicit half-step; 2=expl full-step'; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='JAC_RECOMPT'; PDESC(I)='0=variable; 1=constant sub-step; 2=const full step '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='CK_OVRSHOOT'; PDESC(I)='0=always take full newton step; 1=line search      '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='SMALL_ESTEP'; PDESC(I)='0=step truncation; 1=look-ahead; 2=step absorption '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='ERRTRUNCABS'; PDESC(I)='absolute temporal truncation error tolerance       '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='ERRTRUNCREL'; PDESC(I)='relative temporal truncation error tolerance       '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='ERRITERFUNC'; PDESC(I)='iteration convergence tolerance for function values'; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='ERR_ITER_DX'; PDESC(I)='iteration convergence tolerance for dx             '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='THRESH_FRZE'; PDESC(I)='threshold for freezing the Jacobian                '; PUNIT(I)='mm      '; isBand(i)=.false.
  I=I+1; PNAME(I)='FSTATE_MIN '; PDESC(I)='fractional minimum value of state                  '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='STEP_SAFETY'; PDESC(I)='safety factor in step-size equation                '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='RMIN       '; PDESC(I)='minimum step size multiplier                       '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='RMAX       '; PDESC(I)='maximum step size multiplier                       '; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='NITER_TOTAL'; PDESC(I)='maximum number of iterations in the implicit scheme'; PUNIT(I)='-       '; isBand(i)=.false.
  I=I+1; PNAME(I)='MIN_TSTEP  '; PDESC(I)='minimum time step length                           '; PUNIT(I)='day     '; isBand(i)=.false.
  I=I+1; PNAME(I)='MAX_TSTEP  '; PDESC(I)='maximum time step length                           '; PUNIT(I)='day     '; isBand(i)=.false.
  
  ! parameter identifier
  I=I+1; PNAME(I)='SOBOL_INDX '; PDESC(I)='indentifier for Sobol parameter set                '; PUNIT(I)='-       '; isBand(i)=.false.
  
  NOUTPAR=I
  
  END SUBROUTINE PARDESCRIBE
END MODULE metaparams
