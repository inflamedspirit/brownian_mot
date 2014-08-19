!-----------------------------------------------------------------------
!
!  RK4 module.  Includes Runge-Kutta integrator, uses odeab_support
!               to store the equations of motion.
!
!-----------------------------------------------------------------------

module rk4

  ! get global variables
  use globals
  use odeab_support
  ! storage for coefficients for integrator
  real(wp) :: c2, c3
  real(wp) :: b1, b2, b3, b4
  real(wp) :: a21, a31, a32, a41, a42, a43


contains

  ! Subroutine to initialize coefficients; call before calling integrator

  subroutine rk4init()
    implicit none
    c2 = 2._wp/5
    c3 = 7._wp/8._wp - 3._wp*sqrt(5._wp)/16
    b1 = (263._wp+24._wp*sqrt(5._wp))/1812
    b2 = 125._wp*(1._wp-8._wp*sqrt(5._wp))/3828
    b3 = 1024._wp/(4869._wp*sqrt(5._wp)-10038)
    b4 = 1._wp/(9._wp/2+3._wp/sqrt(5._wp))
    a21 = 2._wp/5
    a31 = 3._wp*(476._wp*sqrt(5._wp)-963)/1024
    a32 = 5._wp*(757._wp-324._wp*sqrt(5._wp))/1024
    a41 = (2094._wp*sqrt(5._wp)-3365)/6040
    a42 = -(975._wp+3046._wp*sqrt(5._wp))/2552
    a43 = 32._wp*(14595._wp+6374._wp*sqrt(5._wp))/240845
  end subroutine rk4init


  ! Subroutine to perform one Runge-Kutta step from t to t+dt

  subroutine eulerstep(y, t, dt)
    implicit none
    real(wp), dimension(3,2), intent(inout) :: y
    real(wp), dimension(3,2) :: ydot
    real(wp), intent(in) :: t, dt
  

    call odeab_func(t, y, ydot)
    y = y + dt * ydot

  end subroutine eulerstep

  subroutine rk4step(y, t, dt)
    implicit none
    real(wp), dimension(3,2), intent(inout) :: y
    real(wp), dimension(3,2) :: ydot
    real(wp), intent(in) :: t, dt
  
    real(wp), dimension(3,2) :: k1, k2, k3, k4

! 1D version
!    k1 = dt * ydot(t, y)
!    k2 = dt * ydot(t + c2*dt, y + a21*k1)
!    k3 = dt * ydot(t + c3*dt, y + a31*k1 + a32*k2)
!    k4 = dt * ydot(t + dt, y + a41*k1 + a42*k2 + a43*k3)
!    y = y + b1*k1 + b2*k2 + b3*k3 + b4*k4

    call odeab_func(t, y, ydot)
    k1 = dt * ydot
 
    call odeab_func(t + c2*dt, y + a21*k1, ydot)
    k2 = dt * ydot

    call odeab_func(t + c3*dt, y + a31*k1 + a32*k2, ydot)
    k3 = dt * ydot

    call odeab_func(t + dt, y + a41*k1 + a42*k2 + a43*k3, ydot)
    k4 = dt * ydot

    y = y + b1*k1 + b2*k2 + b3*k3 + b4*k4



  end subroutine rk4step

  

end module rk4
