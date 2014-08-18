!-----------------------------------------------------------------------
!
!  Support module for hosc.f90 sample program to demonstrate the 
!    odeab integrator
!
!  This module contains certain setup stuff for the integrator,
!    as well as the equations of motion to solve.
!
!-----------------------------------------------------------------------

module odeab_support

  use globals

  ! declare precision, max # of steps, and error tolerances for integration
  integer, parameter  :: odeab_prec = wp
  integer, parameter  :: odeab_maxstp = 500     !500 default
  real(wp), parameter :: odeab_atol = 1.e-13_wp !e-13 default
  real(wp), parameter :: odeab_rtol = 1.e-13_wp !e-13 default

  ! the following lines declare internal storage for the odeab
  ! integrator; just leave it, except:
  !   1. change each dimension to match the number of integration variables
  !      (this could be dimensions like '(2,2)' for a 2x2 matrix)
  !   2. if your variables are all real, change "complex"es to "real"
  type odeab_type
    real(wp), dimension(3,2) :: phi
  end type
  type(odeab_type), dimension(16), save :: odeab_idx
  real(wp), dimension(3,2), target :: odeab_y
  real(wp), dimension(3,2), save   :: odeab_yy, odeab_p, odeab_yp
  real(wp), dimension(3,2), save      :: odeab_wt
 
  ! declare basic management stuff for integrator (just leave it)
  integer :: odeab_istate = 1
  logical, parameter :: odeab_stop = .false.


contains


! Subroutine that implements equations of motion, to be called
!   by odeab; implements simple complex rotation

subroutine odeab_func(t, rp, rpdot)
  implicit none
  real(wp), intent(in) :: t
  real(wp), dimension(3,2), intent(in) :: rp
  real(wp), dimension(3,2), intent(out) :: rpdot
  real(wp), dimension(3) :: r,p,rdot,pdot
  complex(wp) :: sum1, sum2, sum3, sum4
  integer j, jj

  real(wp), dimension(num_beams) :: Delta
  real(wp), dimension(3) :: dHdr = 0.0_wp
  real(wp), dimension(3) :: dHdp = 0.0_wp

  r = rp(:,1)
  p = rp(:,2)

  ! Detuning Matrix Construction
  do j=1,num_beams
     Delta(j) = Delta_base(j) - dot_product(k(:,j),p)/m - (hbar * dot_product(k(:,j), k(:,j)) / (2.0_wp*m))
  end do

  ! The following are the construction of H, dH/dr, and dH/dp. H isn't explicitely needed,
  ! but it should be a good indicator of any cooling due to the spontaneous emission.


  ! Hamiltonian construction
  sum1 = 0.0_wp
  sum2 = 0.0_wp
  do j=1,num_beams
     sum1 = sum1 + (hbar/2.0_wp)*Rabi(j)*exp(i*dot_product(k(:,j),r))/(2.0_wp*Delta(j)+i*Gamma)
     sum2 = sum2 + conjg(Rabi(j))*exp(-i*dot_product(k(:,j),r))
  end do
  H = dot_product(p,p)/(2.0_wp*m) + 2.0_wp*real(sum1*sum2)

  ! dH/dr construction
  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp
     sum3 = 0.0_wp
     sum4 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*i*Rabi(jj)*k(j,jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
        sum3 = sum3 + (hbar/2.0_wp)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
        sum4 = sum4 + (-i)*k(j,jj)*conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdr(j) = 2.0_wp*real(sum1*sum2)+2.0_wp*real(sum3*sum4)
  end do

  ! dH/dp construction
  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(j,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(j) = p(j)/m + 2.0_wp*real(sum1*sum2)
  end do

  ! Classical Hamiltonian Evolution
  rdot = dHdp
!  rdot = p/m ! WHOAH!? Why was I doing this?
  pdot = -dHdr

  ! Load back into rp matrix.
  rpdot(:,1) = rdot
  rpdot(:,2) = pdot

!  write(0,*) 'r:', r 
!  write(0,*) 'p:', p
!  write(0,*) 'H:', H
!  write(0,*) 'k:', k
!  write(0,*) 'Delta:', Delta
!  write(0,*) 'dH/dx', dHdr
!  write(0,*) 'dH/dp:', dHdp

  return
end subroutine odeab_func

end module odeab_support
