!=======================================================================
!
!  odeab90 module
!
!  Version 1.0
!
!  Main interface subroutine:
!
!     odeab (t, tout)
!                      
!  This subroutine generates strong solutions to ordinary
!  differential equations of the form
!    
!     dy = f(y, t) dt
!
!  where y is a vector, and f is a vector-valued functions.
!  This routine implements Adams-Bashforth formulae with adaptive
!  steps, and is modelled after the DE/STEP/INTERP integrator by 
!  Shampine and Gordon and the LSODE integrator by Hindmarsh.
!  This integrator is appropriate for high-accuracy calculations,
!  and should be very efficient when the derivative function call
!  is expensive.  This integrator is also appropriate for nonstiff
!  or mildly stiff problems.
!          
!
!  References:
!    
!  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
!     Solvers, in Scientific Computing, R. S. Stepleman et al., Eds.
!     (North-Holland, Amsterdam, 1983) pp. 55-64.
!
!  L. F. Shampine and M. K. Gordon, Computer Solution of Ordinary
!     Differential Equations: The Initial Value Problem
!     (W. H. Freeman, San Francisco, 1975).
!
!                      
!  Interface:
!
!  The interface seems to be somewhat complicated, but allows the
!    integrator to handle a solution array of any real or complex
!    kind, as well as any shape or size.  This is because the user
!    defines and allocates the solution and storage arrays, and
!    ODEAB only uses array operations internally.  This also has
!    the advantage that this routine can be easily adapted to
!    parallel form (e.g., in HPF) merely by adding directives.
!    Furthermore, configuration should not be too difficult
!    by following the sample module included with this integrator.
!
!  Most of the communication with this integrator is through a
!    user-defined module called "odeab_support". A sample module
!    is included as a separate file with this software.  This module
!    should define several parameters, variables, and functions:
!                    
!    odeab_prec -> (default integer parameter) the real/complex kind
!                  for the precision type to use for all computations,
!                  e.g., this should be the output of selected_real_kind.
!
!    odeab_atol -> (real(odeab_prec)) absolute tolerance for 
!                  integration, combined with odeab_rtol to
!                  control error management.
!
!    odeab_rtol -> (real(odeab_prec)) relative tolerance for 
!                  integration, combined with odeab_atol to
!                  control error management.
!
!    odeab_stop -> (logical) controls whether or not the integrator
!                  can go past tout.  Set this to true if there is
!                  a singularity or other problem with the equations
!                  past tout.  Otherwise a false setting will 
!                  improve efficiency.
!
!    odeab_maxstp -> (default integer) maximum number of steps
!                    to be taken before reaching tout.  If this
!                    number is exceeded, the routine will return
!                    with the solution computed when exceeding
!                    this value.  A value of 500 is reasonable.
!
!    odeab_istate -> (default integer) is both input and output.
!                    This should be set to 1 on the first call
!                    to initialize the integrator.  On return,
!                    a value of 2 indicates success (t is set to
!                    tout, and odeab_y is set to the solution at
!                    tout), while a negative value indicates a
!                    problem:
!                      -1 -> too many steps attempted before
!                            reaching tout (more than odeab_maxstp);
!                            in this case the user can reset
!                            odeab_istate = 2 and continue
!                      -2 -> error tolerances too small 
!                      -3 -> equivalent to (-1), but the problem
!                            is probably stiff
!                    see notes below for input values on subsequent
!                    calls, but normally this should not be
!                    modified after the first call.
!                    
!    odeab_y -> (real/complex array) the solution array, where
!               the type/dimension are set according to the
!               user's requirements.  Because the initial
!               value is set by the user, the user must control
!               the allocation for this array.  May be declared
!               as a pointer to a more sensible array name.
!    
!    odeab_yy   ->  \
!    odeab_p    ->   |  
!    odeab_yp   ->   |
!    odeab_wt   ->  /
!                   auxilliary storage arrays, of the same type and
!                   dimensions as odeab_y.  These should be   
!                   declared and allocated by the user.  
!                   These should not be modified between calls,
!                   although it is possible to directly put
!                   error tolerance weights into odeab_wt.
!
!    odeab_type -> (named type) declare as, e.g.,
!                    type odeab_type
!                      real(odeab_prec), dimension(:), allocatable :: phi
!                    end type
!                  where phi is to have the same type, size, and 
!                  shape as odeab_y (and need not be allocatable).
!
!    odeab_idx -> (type odeab_type array) declare as:
!                    type(odeab_type), dimension(16), save :: odeab_idx
!                 if phi was declared to be allocatable in odeab_type,
!                 then the subroutine will need to allocate it, e.g.,
!                    do j = 1, 16
!                      allocate odeab_idx(j)%phi(n)
!                    end do
!                 if odeab_y is a 1-D array of length n.
!
!    odeab_func -> subroutine of the form
!                    subroutine odeab_func(t, y, ydot)
!                      real(sderk_prec), intent(in) :: t
!                      <solution type>, intent(in)  :: y
!                      <solution type>, intent(out) :: ydot
!                    which returns the derivative array f(y, t). Here,
!                    <solution type> is the array type/dimension
!                    defined above for the solution arrays, e.g., this
!                    could be "real(odeab_prec), dimension(100,3)".
!
!
!  Calling interface:
!
!  The integrator is called to advance odeab_y in time via
!      call odeab (t, tout)
!    where
!
!    t -> (real(odeab_prec)) initial value of the time parameter;
!       on successful integration, this is updated to tout,
!       and on failed integration, this is the value of t
!       to which the solution was advanced.
!
!    tout -> (real(odeab_prec)) the desired output time.
!
!  On output, the value of odeab_istate should be checked for
!    a negative value, in case a problem occurred.  If odeab_istate
!    is 2 on output, everything is ok; to advance the solution
!    again, leave odeab_istate unchanged and call odeab again.
!    To restart the integrator (e.g., because of a discontinuity
!    at the current time), set odeab_istate to 1 and call odeab
!    again.
!
!  The subroutine handle_odeab_error() is available for checking
!    odeab_istate, which simply returns on a positive value and
!    terminates execution with a diagnostic message on a negative
!    value.
!
!  An additional subroutine, print_odeab_stats, is available
!    to print performance statistics to standard error after
!    at least one call to ODEAB.  Set optional input 
!    since_last = .true. to print cumulative statistics since
!    the last initialization of ODEAB, and set the optional input
!    cumulative = .true. to print cumulative statistics since the
!    beginning of execution.
!
!
!=======================================================================


module odeab90

  use odeab_support, only :                                            &
           odeab_func,                                                 &
           odeab_prec, odeab_atol, odeab_rtol, odeab_maxstp,           &
           odeab_y, odeab_yy, odeab_p, odeab_yp, odeab_wt,             &
           odeab_idx, odeab_istate, odeab_stop


  private
  public :: odeab, print_odeab_stats, handle_odeab_error


  !!! local storage for integrator stuff

  ! aliases
  integer, parameter  :: wp = odeab_prec

  ! other storage
  real(wp), save :: x, eps, h
  real(wp), save :: xold, told, hold
  integer, save :: ord, oldord, ordle4, old_istate, ns
  logical, save  :: start, stiff, crash, no_rnd_control
  real(wp), dimension(12), save :: psi

  ! statistics
  integer, save :: function_calls, function_calls_last, function_calls_all = 0
  integer, save :: success_steps, success_steps_last, success_steps_all = 0
  integer, save :: failed_steps, failed_steps_last, failed_steps_all = 0


contains


!-----------------------------------------------------------------------
!
!  main subroutine: odeab
!
!-----------------------------------------------------------------------

subroutine odeab(t, tout)

  implicit none

  real(wp), intent(inout) :: t
  real(wp), intent(in)    :: tout

  real(wp) :: delt, tend, aeps, reps
  integer  :: stepctr


  !!! Check inputs
  if ( t .ge. tout ) then
    write(0,*) 'Error (ODEAB): t must be before tout'
    stop
  end if

  eps = max(odeab_atol, odeab_rtol)
  if ( (odeab_atol .lt. 0 .or. odeab_rtol .lt. 0) .or. eps .le. 0 ) then
    write(0,*) 'Error (ODEAB): bad tolerance settings'
    stop
  end if

  if ( odeab_istate .le. 0 ) then
    write(0,*) 'Error (ODEAB): called with bad istate'
    stop
  end if

  !!! Set tolerances and counters to be set on each step
  delt = tout - t
  tend = t + 10*delt  ! don't integrate past this
  if ( odeab_stop ) tend = min(tend, tout)
  stepctr = 0
  ordle4 = 0
  stiff = .false.
  reps = odeab_rtol / eps
  aeps = odeab_atol / eps

  !!! reset statistics
  if ( odeab_istate .eq. 1 ) then
    function_calls = 0
    success_steps = 0
    failed_steps = 0
  end if
  function_calls_last = function_calls
  success_steps_last = success_steps
  failed_steps_last = failed_steps

  !!! Set work storage on start
  if ( odeab_istate .eq. 1 .or. old_istate .lt. 0 ) then
    start = .true.
    x = t
    odeab_yy = odeab_y
    h = max(tout-x, 4*epsilon(1._wp)*abs(x))
  end if

  !!!!!! main loop starts here
  do 
 
  !!! if already past tout, interpolate and return
  if ( x-t .ge. delt ) then
    call odeab_interp(tout)
    odeab_istate = 2
    t = tout
    told = t
    old_istate = odeab_istate
    return
  end if

  !!! if impossible to go past tout but close, extrapolate and return
  if ( odeab_stop .and. abs(tout-x) .lt. 4*epsilon(1._wp)*abs(x) ) then
    h = tout - x
    call increment_function_calls()
    call odeab_func(x, odeab_yy, odeab_yp)
    odeab_y = odeab_yy + h*odeab_yp
    t = tout
    told = t
    old_istate = odeab_istate
    return
  end if

  !!! test for too much work
  if ( stepctr .ge. odeab_maxstp ) then
    odeab_istate = -1
    if (stiff) odeab_istate = -3
    odeab_y = odeab_yy
    t = x
    told = t
    old_istate = 2
    return
  end if

  !!! limit step size, set weights, and take a step
  h = min( abs(h), abs(tend-x) )
  odeab_wt = reps*abs(odeab_yy) + aeps
  call odeab_step(h)

  !!! check tolerances
  if ( crash ) then
    odeab_istate = -2
    odeab_y = odeab_yy
    t = x
    told = t
    old_istate = -2
    return
  end if

  !!! increment counter and test for stiffness
  stepctr = stepctr + 1
  ordle4 = ordle4 + 1
  if ( oldord .gt. 4 ) ordle4 = 0
  if ( ordle4 .ge. 50 ) stiff = .true.

  end do

end subroutine odeab


!-----------------------------------------------------------------------
!
!  subroutine odeab_step(h)
!
!  Internal stepping routine for odeab.
!
!-----------------------------------------------------------------------

subroutine odeab_step(h)

  implicit none

  real(wp), intent(inout) :: h

  logical, save :: phase1
  integer :: ifail, i, j, neword
  real(wp) :: round, ysum, ah, tmp1=0, tmp2=0, hnew
  real(wp) :: errord, errordm1, errordm2, errordp1, err
  real(wp), dimension(13) :: gstr = (/ 0.50000_wp, 0.08330_wp, 0.04170_wp, &
                                       0.02640_wp, 0.01880_wp, 0.01430_wp, &
                                       0.01140_wp, 0.00936_wp, 0.00789_wp, &
                                       0.00679_wp, 0.00592_wp, 0.00524_wp, &
                                       0.00468_wp /)
  real(wp), dimension(12), save :: alpha, beta, w, v
  real(wp), dimension(13) :: sig = (/ 1,0,0,0,0,0,0,0,0,0,0,0,0 /)
  real(wp), dimension(13) :: g = (/ 1._wp,0.5_wp,0._wp,0._wp,0._wp,0._wp,&
                                    0._wp,0._wp,0._wp,0._wp,0._wp,0._wp,0._wp /)
  real(wp), dimension(13) :: two = &
                           (/ 2,4,8,16,32,64,128,256,512,1024,2048,4096,8192 /)


  !!! check step size
  if ( abs(h) .lt. 4*epsilon(1._wp)*abs(x) ) then
    h = 4*epsilon(1._wp)*abs(x)
    crash = .true.
    return
  end if

  !!! check tolerance
  round = sum( abs(odeab_yy / odeab_wt)**2 )
  round = 2*epsilon(1._wp)*sqrt(round)
  if ( 0.5_wp*eps .lt. round ) then
    crash = .true.
    return
  end if

  !!! block 0
  crash = .false.
  g(1:2) = (/ 1._wp, 0.5_wp /)
  sig(1) = 1
  if ( start ) then
    call increment_function_calls()
    call odeab_func(x, odeab_yy, odeab_yp)
    odeab_idx(1)%phi = odeab_yp
    odeab_idx(2)%phi = 0
    ysum = sqrt( sum( abs(odeab_yp / odeab_wt)**2 ) )
    ah = h
    if ( eps .lt. 16*ysum*h*h )  ah = 0.25_wp*sqrt(eps/ysum)
    h = max(ah, 4*epsilon(1._wp))
    hold = 0
    ord = 1
    oldord = 0
    start = .false.
    phase1 = .true.
    no_rnd_control = .true.

    if ( 0.5_wp*eps .le. 100*round ) then
      no_rnd_control = .false.
      odeab_idx(15)%phi = 0
    end if
  end if
  ifail = 0

  !!! block 1
  !!! compute coefficients
  do 
  if ( h .ne. hold )  ns = 0
  if ( ns .le. oldord ) ns = ns + 1
  if ( ord .ge. ns ) then
    !!! compute components of alpha, beta, psi, sig, that are changed
    beta(ns) = 1
    alpha(ns) = 1.0_wp/ns
    tmp1 = h*ns
    sig(ns+1) = 1
    do i = ns+1, ord
      tmp2 = psi(i-1)
      psi(i-1) = tmp1
      beta(i) = beta(i-1)*psi(i-1)/tmp2
      tmp1 = tmp2 + h
      alpha(i) = h/tmp1
      sig(i+1) = i*alpha(i)*sig(i)
    end do
    psi(ord) = tmp1

    !!! compute coefficients g; initialize v and set w; g(2) is already set
    if ( ns .le. 1 ) then
      do i = 1, ord
        v(i) = 1.0_wp / (i*(i+1))
        w(i) = v(i)
      end do
    else
      if ( ord .gt. oldord ) then
        v(ord) = 1.0_wp / (ord*(ord+1))
        do i = 1, ns-2
          v(ord-i) = v(ord-i) - alpha(i+1)*v(ord-i+1)
        end do
      end if
      do i = 1, ord+1-ns
        v(i) = v(i) - alpha(ns) * v(i+1)
        w(i) = v(i)
      end do
      g(ns+1) = w(1)
    end if

    !!! compute the g in the work vector w
    if ( ord .ge. ns+1 ) then
      do i = ns+2, ord+1
        do j = 1, ord+2-i
          w(j) = w(j) - alpha(i-1)*w(j+1)
        end  do
        g(i) = w(1)
      end do
    end if
  end if


  !!! block 2
  !!! predict solution odeab_p, estimate errors 
  do i = ns+1, ord
    odeab_idx(i)%phi = beta(i) * odeab_idx(i)%phi
  end do

  !!! predict solution and differences
  odeab_idx(ord+2)%phi = odeab_idx(ord+1)%phi
  odeab_idx(ord+1)%phi = 0
  odeab_p = 0
  do j = 1, ord
    odeab_p = odeab_p + g(ord+1-j) * odeab_idx(ord+1-j)%phi
    odeab_idx(ord+1-j)%phi = odeab_idx(ord+1-j)%phi + odeab_idx(ord+2-j)%phi 
  end do
  if ( no_rnd_control ) then
    odeab_p = odeab_yy + h * odeab_p
  else
    odeab_idx(16)%phi = h*odeab_p - odeab_idx(15)%phi
    odeab_p = odeab_yy + odeab_idx(16)%phi
    odeab_idx(16)%phi = (odeab_p - odeab_yy) - odeab_idx(16)%phi
  end if
  xold = x
  x = x + h
  ah = abs(h)
  call increment_function_calls()
  call odeab_func(x, odeab_p, odeab_yp) 

  !!! estimate errors
  if ( ord .gt. 2 ) then
     errordm2 = &
        sum( abs((odeab_idx(ord-1)%phi+odeab_yp-odeab_idx(1)%phi)/odeab_wt)**2 )
     errordm2 = ah*sig(ord-1)*gstr(ord-2)*sqrt(errordm2)
  end if
  if ( ord .gt. 1 ) then
    errordm1 = &
          sum( abs((odeab_idx(ord)%phi+odeab_yp-odeab_idx(1)%phi)/odeab_wt)**2 )
    errordm1 = ah*sig(ord)*gstr(ord-1)*sqrt(errordm1)
  end if
  errord = sum( abs((odeab_yp-odeab_idx(1)%phi)/odeab_wt)**2 )
  err = ah*sqrt(errord)*(g(ord)-g(ord+1))
  errord = ah*sqrt(errord)*sig(ord+1)*gstr(ord)
  neword = ord

  !!! test if order should be lowered
  if ( ord .gt. 2 ) then
    if ( max(errordm1, errordm2) .le. errord ) neword = ord-1
  else if ( ord .eq. 2 ) then
    if ( errordm1 .le. 0.5_wp*errord ) neword = ord-1
  end if

  !!! test if step successful
  if ( err .gt. eps ) then
  
    !!! block 3 -- unsuccessful
    call increment_failed_steps()
    phase1 = .false.
    x = xold
    do i = 1, ord
      odeab_idx(i)%phi = (odeab_idx(i)%phi - odeab_idx(i+1)%phi)/beta(i)
    end do
    if ( ord .ge. 2 ) then
      do i = 2, ord
        psi(i-1) = psi(i) - h
      end do
    end if
    ifail = ifail + 1
    tmp2 = 0.5_wp
    if ( ifail .gt. 3 .and. eps .lt. errord/2 ) tmp2 = sqrt(0.5_wp*eps/errord)
    if ( ifail .gt. 2 )  neword = 1
    h = tmp2*h
    ord = neword
    if ( h .lt. 4*epsilon(1._wp)*abs(x) ) then
      crash = .true.
      h = 4*epsilon(1._wp)*abs(x)
      eps = eps + eps
      return
    end if
    cycle
  end if

  !!! block 4 -- successful
  call increment_success_steps()
  oldord = ord
  hold = h

  !!! correct and evaluate
  if ( no_rnd_control ) then
    odeab_yy = odeab_p + h*g(ord+1)*(odeab_yp - odeab_idx(1)%phi)
  else
    odeab_idx(15)%phi = &
                    h*g(ord+1)*(odeab_yp - odeab_idx(1)%phi) - odeab_idx(16)%phi
    odeab_yy = odeab_p + odeab_idx(15)%phi
    odeab_idx(15)%phi = (odeab_yy - odeab_p) - odeab_idx(15)%phi
  end if
  call increment_function_calls()
  call odeab_func(x, odeab_yy, odeab_yp)

  !!! update differences for next step
  odeab_idx(ord+1)%phi = odeab_yp - odeab_idx(1)%phi
  odeab_idx(ord+2)%phi = odeab_idx(ord+1)%phi - odeab_idx(ord+2)%phi
  do i = 1, ord
    odeab_idx(i)%phi = odeab_idx(i)%phi + odeab_idx(ord+1)%phi
  end do

  !!! estimate error at order ord+1 under appropriate conditions
  errordp1 = 0
  if ( neword .eq. ord-1 .or. ord .eq. 12 )  phase1 = .false.
  if ( phase1 ) then
    ! raise order
    ord = ord + 1
    errord = errordp1
  else if ( neword .eq. ord-1 ) then
    ! lower order
    ord = ord - 1
    errord = errordm1
  else if ( ord+1 .le. ns ) then
    errordp1 = sum( abs((odeab_idx(ord+2)%phi)/odeab_wt)**2 )
    errordp1 = ah*gstr(ord+1)*sqrt(errordp1)

    ! determine appropriate order for next step
    if ( ord .le. 1 ) then
      if ( errordp1 .lt. 0.5_wp*errord ) then
        ! raise order
        ord = ord + 1
        errord = errordp1
      end if
    else if ( errordm1 .le. min(errord, errordp1) ) then
      ! lower order
      ord = ord - 1
      errord = errordm1
    else if ( errordp1 .lt. errord .and. ord .ne. 12 ) then
      ! raise order
      ord = ord + 1
      errord = errordp1
    end if
  end if

  !!! with new order determine next step size
  hnew = h + h
  if ( phase1 .or. 0.5_wp*eps .ge. errord*two(ord+1) ) then
    h = hnew
    return
  end if

  hnew = h
  if ( 0.5_wp*eps .ge. errord ) return
  hnew = ah*max(0.5_wp, min(0.9_wp, (0.5_wp*eps/errord)**(1._wp/(ord+1))))
  hnew = max(hnew, 4*epsilon(1._wp)*abs(x), h)
  h = hnew
  return

  end do

end subroutine odeab_step


!-----------------------------------------------------------------------
!
!  subroutine odeab_interp(tout)
!
!  Internal interpolation routine for odeab, computes solution 
!  odeab_y at time tout.
!
!-----------------------------------------------------------------------

subroutine odeab_interp(tout)

  implicit none

  real(wp), intent(in) :: tout

  real(wp), dimension(13) :: g = (/ 1,0,0,0,0,0,0,0,0,0,0,0,0 /)
  real(wp), dimension(13) :: rho = (/ 1,0,0,0,0,0,0,0,0,0,0,0,0 /)
  real(wp), dimension(13) :: w
  real(wp) :: hi, term
  integer :: ki, i, j

  hi = tout - x
  ki = oldord + 1

  ! initialize w for computing g
  do i = 1, ki
    w(i) = 1._wp / i
  end do
  term = 0

  ! compute g
  do j = 2, ki
    do i = 1, ki+1-j
      w(i) = ((hi+term)*w(i) - hi*w(i+1))/psi(j-1)
    end do
    g(j) = w(1)
    rho(j) = (hi+term)/psi(j-1)
    term = psi(j-1)
  end do

  ! interpolate
  odeab_y = 0
  do j = 1, ki
    odeab_y = odeab_y + g(ki+1-j)*odeab_idx(ki+1-j)%phi
  end do
  odeab_y = odeab_yy + hi*odeab_y

  return

end subroutine odeab_interp


!-----------------------------------------------------------------------
!
!  subroutine increment_function_calls()
!
!  Internal routine to increment function_call counters (module 
!  globals)
!
!-----------------------------------------------------------------------

subroutine increment_function_calls()
  implicit none
  function_calls = function_calls + 1
  function_calls_all = function_calls_all + 1
end subroutine increment_function_calls 


!-----------------------------------------------------------------------
!
!  subroutine increment_success_steps()
!
!  Internal routine to increment success_steps counters (module 
!  globals)
!
!-----------------------------------------------------------------------

subroutine increment_success_steps()
  implicit none
  success_steps = success_steps + 1
  success_steps_all = success_steps_all + 1
end subroutine increment_success_steps 


!-----------------------------------------------------------------------
!
!  subroutine increment_failed_steps()
!
!  Internal routine to increment failed_steps counters (module 
!  globals)
!
!-----------------------------------------------------------------------

subroutine increment_failed_steps()
  implicit none
  failed_steps = failed_steps + 1
  failed_steps_all = failed_steps_all + 1
end subroutine increment_failed_steps 


!-----------------------------------------------------------------------
!
!  subroutine print_odeab_stats(optional since_last, optional cumulative)
!
!  Subroutine to print various statistics about the performance
!  of the odeab integrator.  Optional inputs are since_last (to
!  display cumulative statistics since the last initialization of
!  the integrator) and cumulative (to display cumulative statistics
!  for the integrator since the start of execution).
!
!-----------------------------------------------------------------------

subroutine print_odeab_stats(since_last, cumulative)

  logical, optional, intent(in) :: since_last, cumulative

  logical :: since_last_in, cumulative_in

  since_last_in = .false.
  cumulative_in = .false.
  if ( present(since_last) )  since_last_in = since_last
  if ( present(cumulative) )  cumulative_in = cumulative

  write(0,*) 'ODEAB statistics:'
  write(0,*) '  Last step size: ', hold
  write(0,*) '  Last order: ', oldord
  write(0,*) '  Using rounding-error abatement: ', .not. no_rnd_control
  write(0,*) '  Successful steps (since last ODEAB call): ', &
                                              success_steps - success_steps_last
  write(0,*) '  Failed steps (since last ODEAB call): ', &
                                                failed_steps - failed_steps_last
  write(0,*) '  Function calls (since last ODEAB call): ', &
                                            function_calls - function_calls_last
  if ( since_last_in ) then
    write(0,*) '  Successful steps (since last initialization): ', success_steps
    write(0,*) '  Failed steps (since last initialization): ', failed_steps
    write(0,*) '  Function calls (since last initialization): ', function_calls
  end if
  if ( cumulative_in ) then
    write(0,*) '  Successful steps (cumulative): ', success_steps_all
    write(0,*) '  Failed steps (cumulative): ', failed_steps_all
    write(0,*) '  Function calls (cumulative): ', function_calls_all
  end if

end subroutine print_odeab_stats


!-----------------------------------------------------------------------
!
!  subroutine handle_odeab_error()
!
!  Subroutine to check value of odeab_istate after a call to odeab,
!  and terminate execution with an error message is something is
!  wrong.
!
!-----------------------------------------------------------------------

subroutine handle_odeab_error()
  implicit none
  if ( odeab_istate .eq. 2 )  return
  if ( odeab_istate .eq. -1 ) then
     write(0,*) &
      'Error(HANDLE_ODEAB_ERROR): too many steps attempted before reaching tout'
  else if ( odeab_istate .eq. -2 ) then
    write(0,*) &
      'Error(HANDLE_ODEAB_ERROR): error tolerances too small'
  else if ( odeab_istate .eq. -3 ) then
    write(0,*) &
    'Error(HANDLE_ODEAB_ERROR): too many steps attempted; problem appears stiff'
  else
    write(0,*) &
    'Error(HANDLE_ODEAB_ERROR): impossible ODEAB_ISTATE after ODEAB call: ', &
      odeab_istate
  end if
  stop
end subroutine handle_odeab_error

end module odeab90
