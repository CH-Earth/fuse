module implicit_solve_module

 ! data types
 use nrtype                              ! variable types, etc.
 use data_types, only: parent            ! parent fuse data structure

 ! modules
 use xtry_2_str_module                   ! puts state vector into FUSE state structure
 use str_2_xtry_module                   ! puts FUSE state structure into state vector

 ! global data
 use model_defn, only: nState            ! number of state variables
 use multiforce, only: dt => deltim      ! time step
 use globaldata, only: isDebug           ! print flag

 use model_numerix, only: NUM_FUNCS      ! number of function calls
 use model_numerix, only: NUM_JACOBIAN   ! number of times Jacobian is calculated

 implicit none

 private
 public :: implicit_solve

 contains

 ! ----- calculate dx/dt=g(x) -----------------------------------------------------------
 subroutine dx_dt(fuseStruct, x_try, g_x, J_g)
 use MOD_DERIVS_DIFF_module, only: MOD_DERIVS_DIFF      ! compute dx/dt
 implicit none
 ! input
 type(parent) , intent(inout)            :: fuseStruct  ! parent fuse data structure
 real(sp)     , intent(in)               :: x_try(:)    ! trial state vector
 ! output
 real(sp)     , intent(out)              :: g_x(:)      ! dx/dt=g(x)
 real(sp)     , intent(out)  , optional  :: J_g(:,:)    ! flux Jacobian matrix
 ! internal
 logical(lgt)                            :: comp_dflux  ! flag to compute flux derivatives
 ! --------------------------------------------------------------------------------------

 comp_dflux = present(J_g)

 ! put data in structure
 call XTRY_2_STR(x_try, fuseStruct%state1)

 ! run the fuse physics
 if (present(J_g)) then
  call mod_derivs_diff(fuseStruct, g_x, J_g)
 else
  call mod_derivs_diff(fuseStruct, g_x)
 end if
 
 ! track the total number of function calls
 NUM_FUNCS = NUM_FUNCS + 1

 end subroutine dx_dt

 ! ----- calculate the Jacobian of g(x) -------------------------------------------------
 SUBROUTINE jac_flux(fuseStruct, x_try, g_x, lower, upper, Jac)
 IMPLICIT NONE
 ! input-output
 type(parent) , intent(in) :: fuseStruct            ! parent fuse data structure
 REAL(SP), DIMENSION(:), INTENT(IN) :: g_x, lower, upper
 REAL(SP), DIMENSION(:), INTENT(IN) :: x_try
 REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Jac
 ! locals
 type(parent) :: fuseStruct_local
 real(sp), parameter :: eps_rel = 1e-4_sp
 real(sp), parameter :: eps_abs = 1e-6_sp   ! or smaller, but NOT 1e-9 scale
 real(sp), parameter :: h_min = 1e-8_sp 
 INTEGER(I4B) :: j,n
 REAL(SP), DIMENSION(size(x_try)) :: x, xsav, g_ph
 real(sp) :: h_try, h_act

 ! preliminaries
 n = size(x)
 fuseStruct_local = fuseStruct
 x    = x_try
 xsav = x

 ! loop through columns
 do j=1,n

  ! propose one-sided step (NOTE: negative)
  h_try = -max(eps_rel*abs(xsav(j)), eps_abs)

  ! flip sign if necessary
  if(xsav(j) + h_try < lower(j)) h_try = -h_try

  ! compute function from the perturbed vector
  x(j)  = xsav(j) + h_try
  call dx_dt(fuseStruct_local, x, g_ph)
  h_act = x(j) - xsav(j)

  ! compute column in the Jacobian
  Jac(:,j) = (g_ph - g_x) / h_act

  ! safety: save full vector and data structure
  fuseStruct_local = fuseStruct ! restores consistency after finite differencing
  x = xsav

 end do ! looping through Jacobian columns

 NUM_JACOBIAN = NUM_JACOBIAN + 1   ! keep track of the number of iterations
 end SUBROUTINE jac_flux

 ! ----- simple implicit solve for differentiable model  --------------------------

 subroutine implicit_solve(fuseStruct, x0, x1, nx, ierr, message, isVerbose)
 USE nr, ONLY : lubksb,ludcmp
 USE overshoot_module, only : get_bounds             ! get state bounds
 USE overshoot_module, only : fix_ovshoot            ! fix overshoot (soft clamp)
 USE conserve_clamp_module, only: conserve_clamp     ! fix overshoot and disaggregate fluxes to conserve mass
 USE model_numerix, only:  ERR_ITER_FUNC             ! Iteration convergence tolerance for function values
 USE model_numerix, only:  ERR_ITER_DX               ! Iteration convergence tolerance for dx
 implicit none
 ! input-output
 type(parent), intent(inout)        :: fuseStruct    ! parent fuse data structure
 real(sp)    , intent(in)           :: x0(:)         ! state vector at start of step
 real(sp)    , intent(out)          :: x1(:)         ! state vector at end of step
 integer(i4b), intent(in)           :: nx            ! number of state variables
 ! error control
 integer(i4b), intent(out)          :: ierr          ! error code
 character(*), intent(out)          :: message       ! error message
 logical(lgt), intent(in), optional :: isVerbose     ! flag for printing (subroutine argument)
 logical(lgt)                       :: isPrint       ! flag for printing (local flag)
 ! internal: newton iterations       
 real(sp)                           :: x_old(nx)     ! old trial state vector
 real(sp)                           :: x_try(nx)     ! trial state vector
 real(sp)                           :: g_x(nx)       ! dx/dt=g(x)
 real(sp)                           :: res(nx)       ! residual vector
 real(sp)                           :: Ja(nx,nx)     ! Jacobian matrix (flux)
 real(sp)                           :: Jg(nx,nx)     ! Jacobian matrix (flux)
 real(sp)                           :: Jac(nx,nx)    ! Jacobian matrix (full)
 real(sp)                           :: dx(nx)        ! state update
 real(sp)                           :: phi           ! half squared residual norm
 real(sp)                           :: d             ! determinant sign tracker
 integer(i4b)                       :: indx(nx)      ! LU pivot indices (row-swap bookkeeping)
 integer(i4b)                       :: i             ! index of state
 integer(i4b)                       :: it            ! index of newton iteration
 integer(i4b), parameter            :: maxit=100     ! maximum number of iterations
 logical(lgt)                       :: converged     ! flag for convergence
 ! internal: backtracking line search w/ overshoot reject 
 real(sp)                           :: xnorm         ! norm used in maximum step
 real(sp)                           :: dxnorm        ! norm used to evaluate step size
 real(sp)                           :: stpmax        ! the maximum step
 real(sp)                           :: dxScale       ! used to scale dx if dxnorm > stpmax
 real(sp)                           :: gpsi(nx)      ! function gradient: func = 0.5*sum(res*res)
 real(sp)                           :: slope         ! direction of decrease
 real(sp)                           :: lambda        ! backtrack length multiplier (lambda*dx)
 real(sp)                           :: alamin        ! minimum lambda
 real(sp)                           :: lam_i         ! maximum lambda for the i-th state
 real(sp)                           :: lam_max       ! maximum lambda
 real(sp)                           :: lower(nx)     ! lower bound
 real(sp)                           :: upper(nx)     ! lower bound
 real(sp)                           :: dclamp(nx)    ! derivative in the clamp
 real(sp)                           :: x_trial(nx)   ! state vector for backtrack
 real(sp)                           :: g_trial(nx)   ! dx/dt=g(x) for backtrack
 real(sp)                           :: res_trial(nx) ! residual for backtrack
 real(sp)                           :: phi_new       ! half squared residual norm
 integer(i4b)                       :: ls_it         ! index of line search iteration
 logical(lgt)                       :: ovshoot       ! flag for overshoot
 logical(lgt)                       :: accepted      ! flag for accepting newton step
 real(sp)                           :: phi_best      ! best function evaluation
 real(sp)                           :: x_best(nx)    ! best state vector
 real(sp)                           :: g_best(nx)    ! dx/dt = g(x_best)
 logical(lgt)                       :: have_best     ! check if found a state vector
 logical(lgt)                       :: isClamped     ! check if fallback is clamped
 ! algorithmic control parameters (most passed through MODULE model_numerix)
 REAL(SP), PARAMETER                :: TOLMIN=1.0e-10_sp ! check for spurious minima
 REAL(SP), PARAMETER                :: STPMX=100.0_sp    ! maximum step in lnsrch
 real(sp), parameter                :: shrink   = 0.5_sp
 real(sp), parameter                :: dampen   = 0.1_sp
 real(sp), parameter                :: phi_rel_tol = 1e-5_sp  ! 0.001%
 real(sp), parameter                :: phi_abs_tol = 1e-6_sp
 real(sp), parameter                :: epsb = 1.e-10_sp        ! small safety margin
 integer(i4b), parameter            :: ls_max = 5
 ! ----- procedure starts here --------------------------------------------------------------------
 ! initialize error control
 ierr=0; message='implicit_solve/'

 ! check dimension size
 if (nx /= nState) stop "implicit_solve: nx /= nState"

 ! initialize check for best function evaluation
 phi_best = huge(1._sp); have_best=.false. 

 ! initialize number of calls
 NUM_FUNCS    = 0  ! number of function calls
 NUM_JACOBIAN = 0  ! number of times Jacobian is calculated

 ! get the flag for printing
 isPrint = .false.; if (present(isVerbose)) isPrint = isVerbose

 ! get the bounds for the state variables
 ! NOTE: This can be done outside of the time and iteration loops (keeping here for now)
 call get_bounds(fuseStruct, lower, upper)

 ! put state vector into the fuse data structure
 call XTRY_2_STR(x0, fuseStruct%state0)
 
 ! intialize state vector (and soft clamp)
 x_try  = x0
 x_old  = x_try
 dclamp = 1._sp
 
 ! fix overshoot (only if necessary)
 if(any(x_try < lower) .or. any(x_try > upper)) &
 call fix_ovshoot(x_try, lower, upper, dclamp) 

 ! define maximum step
 xnorm  = sqrt( sum(x_try*x_try) )
 stpmax = STPMX * max( xnorm, real(nx, sp) )

 ! initialize flags
 accepted  = .false.
 converged = .false.

 ! --- F(x), J(x), and objective phi
 call dx_dt(fuseStruct, x_try, g_x, Jg)  ! compute analytical Jacobian
 res = x_try - (x0 + g_x*dt)
 phi = 0.5_sp * dot_product(res, res)

 ! iterate
 do it = 1, maxit

   ! save x
   x_old = x_try

   ! check convergence
   if (phi < ERR_ITER_FUNC) then
     converged = .true.
     exit ! exit iteration loop
   end if

   ! --- compute residual Jacobian J(x) from flux Jacobian Jg(x) ----
   !call jac_flux(fuseStruct, x_try, g_x, lower, upper, Jg)
   do i=1,nx
     Jac(:,i) = -dt*Jg(:,i) 
     Jac(i,i) = Jac(i,i) + 1.0_sp
   end do

   ! --- function gradient: before Jac is modified in ludcmp
   gpsi  = matmul(transpose(Jac), res) ! assumes func =  0.5_sp * sum(res*res)

   ! --- Solve J dx = -F (Newton step)
   dx = -res
   call ludcmp(Jac, indx, d)     ! J overwritten with LU
   call lubksb(Jac, indx, dx)    ! dx becomes solution
     
   ! --- Modify dx

   ! modify dx if norm > stpmax
   dxnorm = sqrt( sum(dx*dx) )
   if (dxnorm > stpmax) then
    dxScale = stpmax / dxnorm
    dx = dxScale * dx
   end if

   ! modify dx if Newton step not descending for psi
   slope = dot_product(gpsi, dx)
   if (slope >= 0._sp) dx = -gpsi ! fallback

   ! implement active-set methods
   do i=1,nx
     if (x_try(i) <= lower(i)+epsb .and. dx(i) < 0._sp) dx(i)=0._sp
     if (x_try(i) >= upper(i)-epsb .and. dx(i) > 0._sp) dx(i)=0._sp
   end do

   ! ---- backtracking line search  --------------

   ! line search control
   accepted = .false. ! flag to check if line search is accepted
   alamin   = ERR_ITER_DX / maxval( abs(dx) / max(abs(x_try), 1.0_sp) )

   lambda = 1.0_sp
   do ls_it = 1, ls_max

     ! update x 
     x_trial   = x_try + lambda*dx

     ! shrink lambda until find a value in the feasible space
     if(any(x_trial < lower) .or. any(x_trial > upper))then
      lambda = lambda * shrink
      cycle
     endif

     ! compute function and function eval -- no need for the Jacobian here
     call dx_dt(fuseStruct, x_trial, g_trial)
     res_trial  = x_trial - (x0 + dt*g_trial)
     phi_new    = 0.5_sp * dot_product(res_trial, res_trial)

     ! save best function evaluation   
     if (phi_new < phi_best) then
       phi_best  = phi_new
       x_best    = x_trial
       g_best    = g_trial
       have_best = .true.
     endif

     if (phi_new <= phi + phi_abs_tol) then
       accepted = .true.; exit
     endif 
     
     ! update lambda
     lambda = lambda * shrink
     if (lambda < alamin) exit   ! give up shrinking

   end do  ! line search

   ! ----- fallback: try a small step  -----
   if(.not. accepted)then
     x_trial = x_try + dampen*dx
     if(any(x_trial < lower) .or. any(x_trial > upper)) &
     call fix_ovshoot(x_trial, lower, upper, dclamp) 
   end if ! (if accepted)

   ! recompute dx_dt because we need the Jacobian
   x_try  = x_trial
   call dx_dt(fuseStruct, x_try, g_x, Jg)  ! compute analytical Jacobian
   res = x_try - (x0 + g_x*dt)
   phi = 0.5_sp * dot_product(res, res)

   ! save best function evaluation   
   if (phi < phi_best) then
     phi_best  = phi
     x_best    = x_try
     g_best    = g_x
     have_best = .true.
   endif

   ! tiny-step convergence
   if (maxval( abs(x_try - x_old) / max(abs(x_try), 1._sp)  ) < ERR_ITER_DX) then
     converged = .true.
     exit ! exit iteration loop
   end if

 end do  ! loop through iterations

 ! ----- handle the extremely rare case of non-convergence -----
 if( .not. converged)then

   ! use explicit Euler if did not find anything
   if( .not. have_best) call dx_dt(fuseStruct, x0, g_best)

   ! use dx/dt = g(x_best)
   x_try = x0 + dt*g_best

   ! test bounds violations: if bounds exceeded, then clamp and disaggregate fluxes (conserve mass)
   call XTRY_2_STR(x_try, fuseStruct%state1)
   call conserve_clamp(fuseStruct, dt, isClamped)
   print*, 'WARNING: '//trim(message)//"failed to converge: use best function evaluation. Clamp = ", isClamped

 endif  ! if not converged

 ! save final state
 x1 = x_try

 end subroutine implicit_solve

end module implicit_solve_module
