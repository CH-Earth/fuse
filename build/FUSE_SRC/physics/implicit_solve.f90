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

 use model_numerix, only: NUM_FUNCS      ! number of function calls
 use model_numerix, only: NUM_JACOBIAN   ! number of times Jacobian is calculated

 implicit none

 ! provide access to the fuse parent structure
 type(parent), pointer, save :: ctx => null()

 private
 public :: implicit_solve

 contains

 ! ----- point to the fuse parent structure ---------------------------------------------

  subroutine set_dxdt_context(fuseStruct)
    type(parent), target, intent(inout) :: fuseStruct
    ctx => fuseStruct
  end subroutine set_dxdt_context

  subroutine clear_dxdt_context()
    nullify(ctx)
  end subroutine clear_dxdt_context

 ! --------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------
 
 ! ----- calculate dx/dt=g(x) -----------------------------------------------------------
 function dx_dt(x_try) result(g_x)
 use MOD_DERIVS_DIFF_module, only: MOD_DERIVS_DIFF     ! compute dx/dt
 implicit none
 ! input
 real(sp)     , intent(in)    :: x_try(:)              ! trial state vector
 ! output
 real(sp)                     :: g_x(size(x_try))      ! dx/dt=g(x)

 ! check made the association to ctx (ctx=>fuseStruct)
 if (.not. associated(ctx)) stop "dx_dt: context not set"

 ! put data in structure
 call XTRY_2_STR(x_try, ctx%state1)

 ! run the fuse physics
 call mod_derivs_diff(ctx)
 
 ! extract dx_dt from fuse structure
 call STR_2_XTRY(ctx%dx_dt, g_x)

 ! track the total number of function calls
 NUM_FUNCS = NUM_FUNCS + 1

 end function dx_dt

 ! ----- calculate the Jacobian of g(x) -------------------------------------------------
 SUBROUTINE jac_flux(x,g_x,Jac)
 IMPLICIT NONE
 REAL(SP), DIMENSION(:), INTENT(IN) :: g_x
 REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
 REAL(SP), DIMENSION(:,:), INTENT(OUT) :: Jac
 REAL(SP), PARAMETER :: EPS=-1.0e-4_sp  ! NOTE force h to be negative
 INTEGER(I4B) :: j,n
 REAL(SP), DIMENSION(size(x)) :: xsav,xph,h
 xsav=x
 n=size(x)
 h=EPS*abs(xsav)
 where (h == 0.0) h=EPS
 xph=xsav+h
 h=xph-xsav
 do j=1,n
  x(j)=xph(j)
  Jac(:,j)=(dx_dt(x)-g_x(:))/h(j)
  x(j)=xsav(j)
 end do
 NUM_JACOBIAN = NUM_JACOBIAN + 1   ! keep track of the number of iterations
 call XTRY_2_STR(xsav, ctx%state1) ! restores consistency after finite differencing
 end SUBROUTINE jac_flux

 ! ----- simple implicit solve for differentiable model  --------------------------

 subroutine implicit_solve(fuseStruct, x0, x1, nx)
 USE nr, ONLY : lubksb,ludcmp
 USE overshoot_module, only : get_bounds       ! get state bounds
 USE overshoot_module, only : fix_ovshoot      ! fix overshoot (soft clamp)
 USE model_numerix, only:  ERR_ITER_FUNC       ! Iteration convergence tolerance for function values
 USE model_numerix, only:  ERR_ITER_DX         ! Iteration convergence tolerance for dx
 implicit none
 ! input-output
 type(parent), intent(inout)  :: fuseStruct    ! parent fuse data structure
 real(sp)    , intent(in)     :: x0(:)         ! state vector at start of step
 real(sp)    , intent(out)    :: x1(:)         ! state vector at end of step
 integer(i4b), intent(in)     :: nx            ! number of state variables
 ! internal: newton iterations
 real(sp)                     :: x_try(nx)     ! trial state vector
 real(sp)                     :: g_x(nx)       ! dx/dt=g(x)
 real(sp)                     :: res(nx)       ! residual vector
 real(sp)                     :: Jg(nx,nx)     ! Jacobian matrix (flux)
 real(sp)                     :: Jac(nx,nx)    ! Jacobian matrix (full)
 real(sp)                     :: dx(nx)        ! state update
 real(sp)                     :: phi           ! half squared residual norm
 real(sp)                     :: d             ! determinant sign tracker
 integer(i4b)                 :: indx(nx)      ! LU pivot indices (row-swap bookkeeping)
 integer(i4b)                 :: i             ! index of state
 integer(i4b)                 :: it            ! index of newton iteration
 integer(i4b), parameter      :: maxit=100     ! maximum number of iterations
 logical(lgt)                 :: converged     ! flag for convergence
 ! internal: backtracking line search w/ overshoot reject 
 real(sp)                     :: lambda        ! backtrack length multiplier (lambda*dx)
 real(sp)                     :: lower(nx)     ! lower bound
 real(sp)                     :: upper(nx)     ! lower bound
 real(sp)                     :: x_trial(nx)   ! state vectorfor backtrack
 real(sp)                     :: g_trial(nx)   ! dx/dt=g(x) for backtrack
 real(sp)                     :: res_trial(nx) ! residual for backtrack
 real(sp)                     :: phi_new       ! half squared residual norm
 integer(i4b)                 :: ls_it         ! index of line search iteration
 logical(lgt)                 :: ovshoot       ! flag for overshoot
 logical(lgt)                 :: accepted      ! flag for accepting newton step
 ! line search params
 real(sp), parameter          :: shrink   = 0.5_sp
 real(sp), parameter          :: c_armijo = 1e-4_sp
 integer(i4b), parameter      :: ls_max = 5

 ! check dimension size
 if (nx /= nState) stop "implicit_solve: nx /= nState"

 ! initialize number of calls
 NUM_FUNCS    = 0  ! number of function calls
 NUM_JACOBIAN = 0  ! number of times Jacobian is calculated

 ! get the bounds for the state variables
 ! NOTE: This can be done outside of the time and iteration loops (keeping here for now)
 call get_bounds(fuseStruct, lower, upper)

 ! point to the fuse parent structure so that it is available in other routines
 call set_dxdt_context(fuseStruct)

 ! put state vector into the fuse data structure
 call XTRY_2_STR(x0, fuseStruct%state0)
 
 ! intialize state vector and convergence flag
 x_try     = x0
 accepted  = .false.
 converged = .false.

 ! --- F(x) and objective phi = 0.5*||F||^2
 g_x = dx_dt(x_try)
 res = x_try - (x0 + g_x*dt)
 phi = 0.5_sp * sum(res*res)

 ! iterate
 do it = 1, maxit

   if (sqrt(2.0_sp*phi) < ERR_ITER_FUNC) then
     converged = .true.
     exit ! exit iteration loop
   end if

   ! --- J(x)
   call jac_flux(x_try, g_x, Jg)
   Jac = -dt*Jg                                     ! multiply dt
   do i=1,nx; Jac(i,i) = Jac(i,i) + 1.0_sp; end do  ! add identity matrix

   ! --- Solve J dx = -F (Newton step)
   dx = -res
   call ludcmp(Jac, indx, d)     ! J overwritten with LU
   call lubksb(Jac, indx, dx)    ! dx becomes solution

   ! initialize flag to check if line search is accepted
   accepted = .false.

   ! ---- backtracking line search w/ overshoot reject ----
   lambda = 1.0_sp
   do ls_it = 1, ls_max
     x_trial = x_try + lambda*dx

     ! check overshoot
     ovshoot = any(x_trial < lower) .or. any(x_trial > upper)
     if (.not. ovshoot) then
       ! new function and residual
       g_trial   = dx_dt(x_trial)
       res_trial = x_trial - (x0 + dt*g_trial)
       phi_new   = 0.5_sp * sum(res_trial*res_trial)
       ! check for sufficient decrease (Armijo-lite)
       if (phi_new <= (1.0_sp - c_armijo*lambda) * phi)then
         accepted = .true.
         exit
       endif
     end if
     lambda = lambda * shrink
   end do  ! line search

   if (accepted) then
     x_try   = x_trial
     g_x     = g_trial
     res     = res_trial
     phi     = phi_new
   else
     ! ----- fallback: soft clamp a very small Newton step -----
     x_trial = x_try + lambda*dx
     call fix_ovshoot(x_trial, lower, upper) 
     ! get new function evaluation
     x_try = x_trial
     g_x   = dx_dt(x_try)
     res   = x_try - (x0 + g_x*dt)
     phi   = 0.5_sp * sum(res*res)
   end if

   ! re-populate fuse data structure
   call XTRY_2_STR(x_try, fuseStruct%state1)

   ! tiny-step convergence
   if (maxval(abs(lambda*dx)) < ERR_ITER_DX) then
     converged = .true.
     exit ! exit iteration loop
   end if

 end do  ! loop through iterations

 ! save final state
 x1 = x_try

 ! nullify pointer to the fuse structure
 call clear_dxdt_context()

 ! check convergence
 if( .not. converged) STOP "failed to converge in implicit_solve"

 end subroutine implicit_solve

end module implicit_solve_module
