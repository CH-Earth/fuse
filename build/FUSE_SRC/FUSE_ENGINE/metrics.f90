module metrics
  use nrtype
  implicit none
  private
  public :: get_KGE, get_KGEp, get_NSE, get_MAE, get_RMSE

contains

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Allow different transformation on streamflow before using metrics
! ---------------------------------------------------------------------------------------
  
  ! Compute transformation on streamflow
  subroutine apply_transformation(obs, sim, transfo)
  
    real(SP), intent(inout)      :: obs(:), sim(:)     ! observed and simulated streamflow time series
    character(len=*), intent(in) :: transfo            ! transformation applied to streamflow time series
    real(SP)                     :: eps, transfo_val   ! epsilon value to avoid infinite number, power transformation as float
    integer                      :: i                  ! loop index
  
    ! log transformation
    if (transfo == 'log') then
  
      eps = sum(obs) / (100 * size(obs)) ! use to avoid log(0)
    
      do i = 1, size(obs)
          obs(i) = log(eps + obs(i))
          sim(i) = log(eps + sim(i))
      end do
    
    ! box-cox transformation  
    else if (transfo == 'boxcox') then
      do i = 1, size(obs)
        obs(i) = (obs(i)**0.25 - 1) / 0.25
        sim(i) = (sim(i)**0.25 - 1) / 0.25
      end do
      
    ! power transformation
    else
      transfo_val = char_to_float(transfo)
    
      if (transfo_val < 0) then
    
        eps = sum(obs) / (100 * size(obs)) ! use to avoid division by 0
      
        do i = 1, size(obs)
          obs(i) = (eps + obs(i)) ** transfo_val
          sim(i) = (eps+ sim(i)) ** transfo_val
        end do
      
      else
    
        do i = 1, size(obs)
          obs(i) = obs(i) ** transfo_val
          sim(i) = sim(i) ** transfo_val
        end do
      
      end if
    
    end if
  
  end subroutine apply_transformation

! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Allow different metrics as objective function
! ---------------------------------------------------------------------------------------

  ! Compute KGE
  function get_KGE(obs, sim, transfo) result(kge)
  
    real(SP), intent(in)                   :: obs(:), sim(:)                                    ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo                                           ! transformation applied to streamflow time series
    real(SP)                               :: kge                                               ! KGE
    real(SP), allocatable                  :: obs_t(:), sim_t(:)                                      ! transformed time series without NA values
    real(SP)                               :: sd_sim, sd_obs, m_sim, m_obs, r, alpha, beta      ! KGE variables
    integer                                :: n  
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    
    ! compute mean
    m_obs  = sum(obs_t) / n
    m_sim  = sum(sim_t) / n
      
    ! compute standard deviation
    sd_obs  = standard_deviation(obs_t)
    sd_sim  = standard_deviation(sim_t)
              
    ! compute correlation
    r = correlation(sim_t, obs_t)  ! correlation
    alpha = sd_sim / sd_obs        ! variation
    beta = m_sim / m_obs           ! bias
    
    kge = 1.0 - sqrt((r-1)**2 + (alpha-1)**2 + (beta-1)**2)
    
    if (isnan(kge)) then
        kge = -1.0e6                ! bad value
    end if
    
  end function get_KGE


  ! Compute KGEp
  function get_KGEp(obs, sim, transfo) result(kgep)
  
    real(SP), intent(in)                   :: obs(:), sim(:)                                    ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo                                           ! transformation applied to streamflow time series
    real(SP)                               :: kgep                                              ! Kling-Gupta Efficiency' (Kling et al. in 2012)
    real(SP), allocatable                  :: obs_t(:), sim_t(:)                                      ! transformed time series without NA values
    real(SP)                               :: sd_sim, sd_obs, m_sim, m_obs, r, alphap, beta     ! KGE' variables
    integer                                :: n                                                 ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    
    ! compute mean
    m_obs  = sum(obs_t) / n
    m_sim  = sum(sim_t) / n
      
    ! compute standard deviation
    sd_obs  = standard_deviation(obs_t)
    sd_sim  = standard_deviation(sim_t)
              
    ! KGEp component
    r = correlation(sim_t, obs_t)             ! correlation
    alphap = (sd_sim/m_sim) / (sd_obs/m_obs)  ! relative variation
    beta = m_sim / m_obs                      ! bias
    
    kgep = 1.0 - sqrt((r-1)**2 + (alphap-1)**2 + (beta-1)**2)
    
    if (isnan(kgep)) then
        kgep = -1.0e6                          ! bad value
    end if
        
  end function get_KGEp


  ! Compute NSE
  function get_NSE(obs, sim, transfo) result(nse)
  
    real(SP), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(SP)                               :: nse                 ! Nash-Sutcliff score
    real(SP), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    nse = 1 - (sum((obs_t - sim_t)**2) / sum((obs_t - (sum(obs_t)/n))**2))
    
    if (isnan(nse)) then
        nse = -1.0e6                          ! bad value
    end if
    
  end function get_NSE


  ! Compute MAE
  function get_MAE(obs, sim, transfo) result(mae)
  
    real(SP), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(SP)                               :: mae                 ! mean absolute error
    real(SP), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values

        
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    mae = sum(abs(obs_t - sim_t)) / n
    
    if (isnan(mae)) then
        mae = 1.0e6                          ! bad value
    end if
    
  end function get_MAE


  ! Compute RMSE
  function get_RMSE(obs, sim, transfo) result(rmse)
  
    real(SP), intent(in)                   :: obs(:), sim(:)      ! observed and simulated streamflow time series
    character(len=*), intent(in), optional :: transfo             ! transformation applied to streamflow time series
    real(SP)                               :: rmse                ! root mean square error
    real(SP), allocatable                  :: obs_t(:), sim_t(:)  ! transformed time series without NA values
    integer                                :: n                   ! length of the time series without NA values
    
    n = count(.not. (isnan(obs) .or. isnan(sim)))
    
    allocate(obs_t(n), sim_t(n))
    
    obs_t = pack(obs, .not. (isnan(obs) .or. isnan(sim)))
    sim_t = pack(sim, .not. (isnan(obs) .or. isnan(sim)))
    
    call apply_transformation(obs_t, sim_t, transfo)
    
    rmse = sqrt(sum((obs_t - sim_t)**2) / n)
    
    if (isnan(rmse)) then
        rmse = 1.0e6                          ! bad value
    end if
    
  end function get_RMSE


! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Cyril Thébault, 2024
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Additional functions
! ---------------------------------------------------------------------------------------

  ! Function to convert character to float
  function char_to_float(char_val) result(float_val)
  
    character(len=*), intent(in) :: char_val
    real(SP)                     :: float_val
    integer                      :: io_stat

    read(char_val, *, iostat=io_stat) float_val
    if (io_stat /= 0) then
      print *, "Error: Unable to convert '", char_val, "' to float"
      float_val = 1.0  ! set to 1 for no transformation
    end if
  end function char_to_float


  ! Function to compute standard deviation
  function standard_deviation(x) result(sd)
    real(SP), intent(in) :: x(:)
    real(SP)             :: sd
    real(SP)             :: mean_x
    integer              :: n
    
    n = size(x)
    mean_x = sum(x) / n
    sd = sqrt(sum((x - mean_x)**2) / (n - 1))
  end function standard_deviation

  ! Function to compute correlation
  function correlation(x, y) result(r)
    real(SP), intent(in) :: x(:), y(:)
    real(SP)             :: r
    real(SP)             :: mean_x, mean_y, sd_x, sd_y
    integer              :: n, i
    
    n = size(x)
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    sd_x = standard_deviation(x)
    sd_y = standard_deviation(y)
    
    r = sum((x - mean_x) * (y - mean_y)) / (n * sd_x * sd_y)
  end function correlation

  elemental function isnan(x) result(is_nan)
    real(SP), intent(in) :: x
    logical              :: is_nan
  
    is_nan = (x /= x)
  end function isnan

end module metrics