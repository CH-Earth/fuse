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
 function dx_dt(fuseStruct, x_try) result(g_x)
 use MOD_DERIVS_DIFF_module, only: MOD_DERIVS_DIFF     ! compute dx/dt
 implicit none
 ! input
 type(parent) , intent(inout) :: fuseStruct            ! parent fuse data structure
 real(sp)     , intent(in)    :: x_try(:)              ! trial state vector
 ! output
 real(sp)                     :: g_x(size(x_try))      ! dx/dt=g(x)

 ! put data in structure
 call XTRY_2_STR(x_try, fuseStruct%state1)

 ! run the fuse physics
 call mod_derivs_diff(fuseStruct)
 
 ! extract dx_dt from fuse structure
 call STR_2_XTRY(fuseStruct%dx_dt, g_x)

 ! track the total number of function calls
 NUM_FUNCS = NUM_FUNCS + 1

 end function dx_dt

 ! ----- calculate the Jacobian of g(x) -------------------------------------------------
 SUBROUTINE jac_flux(fuseStruct, x_try, g_x, lower, upper, Jac)
 IMPLICIT NONE
 ! input-output
 type(parent) , intent(inout) :: fuseStruct            ! parent fuse data structure
 REAL(SP), DIMENSION(:), INTENT(IN) :: g_x, lower, upper
 REAL(SP), DIMENSION(:), INTENT(IN) :: x_try
 REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Jac
 ! locals
 type(parent) :: ctx_sav
 real(sp), parameter :: eps_rel = 1e-4_sp
 real(sp), parameter :: eps_abs = 1e-6_sp   ! or smaller, but NOT 1e-9 scale
 real(sp), parameter :: h_min = 1e-8_sp 
 INTEGER(I4B) :: j,n
 REAL(SP), DIMENSION(size(x_try)) :: x, xsav, g_ph
 real(sp) :: h_try, h_act

 ! preliminaries
 n = size(x)
 ctx_sav = fuseStruct
 x    = x_try
 xsav = x

 ! loop through columns
 do j=1,n

  ! safety: save full vector and data structure
  fuseStruct = ctx_sav
  x=xsav

  ! propose one-sided step
  h_try = -max(eps_rel*abs(xsav(j)), eps_abs)

  ! flip sign if necessary
  if(xsav(j) + h_try < lower(j)) h_try = -h_try

  ! compute function from the perturbed vector
  x(j)  = xsav(j) + h_try
  g_ph  = dx_dt(fuseStruct, x) 
  h_act = x(j) - xsav(j)

  ! compute column in the Jacobian
  Jac(:,j) = (g_ph - g_x) / h_act

 end do ! looping through Jacobian columns

 NUM_JACOBIAN = NUM_JACOBIAN + 1   ! keep track of the number of iterations
 fuseStruct = ctx_sav ! restores consistency after finite differencing
 end SUBROUTINE jac_flux

 ! ----- simple implicit solve for differentiable model  --------------------------

 subroutine implicit_solve(fuseStruct, x0, x1, nx, ierr, message, isVerbose)
 USE nr, ONLY : lubksb,ludcmp
 USE overshoot_module, only : get_bounds             ! get state bounds
 USE overshoot_module, only : fix_ovshoot            ! fix overshoot (soft clamp)
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
 type(parent)                       :: ctx           ! save the fuse structure
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

 ! initialize number of calls
 NUM_FUNCS    = 0  ! number of function calls
 NUM_JACOBIAN = 0  ! number of times Jacobian is calculated

 ! get the flag for printing
 isPrint = .false.; if (present(isVerbose)) isPrint = isVerbose

 ! save the fuse structure
 ctx = fuseStruct

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

 if(isPrint) isDebug = .true.
 
 ! --- F(x) and objective phi
 g_x = dx_dt(fuseStruct, x_try)
 res = x_try - (x0 + g_x*dt)
 phi = 0.5_sp * dot_product(res, res)

 if(isPrint) isDebug = .false.

 ! iterate
 do it = 1, maxit

   ! save x
   x_old = x_try

   if(isPrint) print*, '***** start of iteration *****'

   ! check
   if(isPrint)then
     print*, 'x_try = ', x_try
     print*, 'g_x = ', g_x
     print*, 'res = ', res
     print*, 'phi = ', phi
     print*, 'dclamp = ', dclamp
     if(it > 10) stop 1
   endif

   if (phi < ERR_ITER_FUNC) then
     converged = .true.
     exit ! exit iteration loop
   end if

   if(isPrint) print*, 'x_try 0 = ', x_try

   ! --- J(x)
   call jac_flux(fuseStruct, x_try, g_x, lower, upper, Jg)
   do i=1,nx
     Jac(:,i) = -dt*Jg(:,i) !* dclamp(i)     ! multiply dt and clamp derivative
     Jac(i,i) = Jac(i,i) + 1.0_sp
   end do

   if(isPrint) print*, 'x_try 1 = ', x_try

   ! --- function gradient: before Jac is modified in ludcmp
   gpsi  = matmul(transpose(Jac), res) ! assumes func =  0.5_sp * sum(res*res)

   ! --- Solve J dx = -F (Newton step)
   dx = -res
   call ludcmp(Jac, indx, d)     ! J overwritten with LU
   call lubksb(Jac, indx, dx)    ! dx becomes solution
     
   if(isPrint) print*, 'x_try 2 = ', x_try

   if(isPrint)then
     print*, 'dx     = ', dx
     print*, 'Jg     = ', Jg
   endif

   ! --- Modify dx

   ! modify dx if norm > stpmax
   dxnorm = sqrt( sum(dx*dx) )
   if (dxnorm > stpmax) then
    dxScale = stpmax / dxnorm
    dx = dxScale * dx
   end if

   ! implement active-set methods
   do i=1,nx
     if (x_try(i) <= lower(i)+epsb .and. dx(i) < 0._sp) dx(i)=0._sp
     if (x_try(i) >= upper(i)-epsb .and. dx(i) > 0._sp) dx(i)=0._sp
   end do

   ! modify dx if Newton step not descending for psi
   slope = dot_product(gpsi, dx)
   if (slope >= 0._sp) dx = -gpsi ! fallback

     if(isPrint) print*, 'x_try 3 = ', x_try
   ! ---- backtracking line search  --------------

   ! save the fuse structure to re-use in subsequent linesearch calls
   ctx = fuseStruct      ! <--- snapshot *now*, for this Newton step

   ! line search control
   accepted = .false. ! flag to check if line search is accepted
   alamin   = ERR_ITER_DX / maxval( abs(dx) / max(abs(x_try), 1.0_sp) )

   ! compute maximum lambda
   lam_max = 1.0_sp
   ! do i=1,nx
   !   if (dx(i) > 0._sp) then
   !     lam_i = (upper(i) - x_try(i) - epsb) / dx(i)
   !     if (lam_i < lam_max) lam_max = max(0._sp, lam_i)
   !   else if (dx(i) < 0._sp) then
   !     lam_i = (lower(i) - x_try(i) + epsb) / dx(i)   ! dx<0 so this is positive
   !     if (lam_i < lam_max) lam_max = max(0._sp, lam_i)
   !   end if
   ! end do

   ! check
   if(isPrint)then
     print*, 'alamin = ', alamin
     print*, 'slope  = ', slope
     print*, 'gpsi   = ', gpsi
   endif

   if(isPrint) isDebug = .true.

   lambda = lam_max
   do ls_it = 1, ls_max

     if(isPrint)then
       print*, '***** new linesearch *****', ls_it
       print*, 'dx     = ', dx
     endif

     ! update x 
     x_trial   = x_try + lambda*dx

     if(isPrint)then
       print*, 'x_try = ', x_try
       print*, 'x_trial = ', x_trial
       print *, "delta = ", x_trial - x_try
       print *, "lambda*dx = ", lambda*dx
     endif

     ! exit if violate bounds: line search direction is not valid
     if(any(x_trial < lower) .or. any(x_trial > upper))then
      accepted = .false.; exit
      !lambda = lambda * shrink
      !cycle
     endif

     ! compute function and function eval
     fuseStruct = ctx
     g_trial    = dx_dt(fuseStruct, x_trial)
     res_trial  = x_trial - (x0 + dt*g_trial)
     phi_new    = 0.5_sp * dot_product(res_trial, res_trial)

     if(isPrint)then
       print*, 'ls_it, lambda, phi, phi_new', ls_it, lambda, phi, phi_new
       print*, 'phi, phi_new, slope=', phi, phi_new, slope
       print*, 'x_trial = ', x_trial
       print*, 'g_trial = ', g_trial
       print*, 'res _trial= ', res_trial
     endif
    
     if (phi_new <= phi + phi_abs_tol) then
       accepted = .true.; exit
     endif 
     
     ! update lambda
     lambda = lambda * shrink
     if (lambda < alamin) exit   ! give up shrinking

   end do  ! line search

   if(isPrint) isDebug = .false.
   
   if (accepted) then
     x_try   = x_trial
     g_x     = g_trial
     res     = res_trial
     phi     = phi_new
   else
     ! ----- fallback: try a small step along the direction of steepest descent -----
     !dx = -gpsi ! use steepest descent
     x_trial = x_try + dampen*dx
     if(any(x_trial < lower) .or. any(x_trial > upper)) &
     call fix_ovshoot(x_trial, lower, upper, dclamp) 
     ! get new function evaluation
     x_try      = x_trial
     fuseStruct = ctx
     g_x        = dx_dt(fuseStruct, x_try)
     res        = x_try - (x0 + g_x*dt)
     phi        = 0.5_sp * dot_product(res, res)
   end if

   ! tiny-step convergence
   if (maxval( abs(x_try - x_old) / max(abs(x_try), 1._sp)  ) < ERR_ITER_DX) then
     converged = .true.
     exit ! exit iteration loop
   end if

 end do  ! loop through iterations

 ! save final state
 x1 = x_try

 ! check convergence
 if( .not. converged)then
  message=trim(message)//"failed to converge"
  ierr=10; return
 endif

 end subroutine implicit_solve

end module implicit_solve_module
