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



  if ( hamiltonian_version .eq. 0 ) then
  ! Detuning Matrix Construction

  ! The following are the construction of H, dH/dr, and dH/dp. H isn't explicitely needed,
  ! but it should be a good indicator of any cooling due to the spontaneous emission.
  call calc_Delta( Delta, p)

  ! Hamiltonian construction ! Don't need to calculate this here, can do it only for datapoints
  ! call calc_H( Delta, r, p )

  ! dH/dr construction
!  call calc_dHdr( dHdr, Delta, r )
  call calc_dHdr2( dHdr, Delta, r )
 ! call calc_dHdr3( dHdr, Delta, r )
 ! call calc_dHdr4( dHdr, Delta, r )

  ! dH/dp construction
!  call calc_dHdp( dHdp, Delta, r, p )
  call calc_dHdp2( dHdp, Delta, r, p )
!  call calc_dHdp3( dHdp, Delta, r, p )
!  call calc_dHdp4( dHdp, Delta, r, p )
  else if ( hamiltonian_version .eq. 1 ) then
     call calc_Delta_dipole( Delta )
     call calc_dHdr_dipole( dHdr, Delta, r )
     call calc_dHdp_dipole( dHdp, p )
  end if

  ! Classical Hamiltonian Evolution
  rdot = dHdp
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

! dHdr calculation (most slowdown!)
subroutine calc_dHdr(dHdr, Delta, r)
  implicit none

  complex(wp) :: sum1, sum2, sum3, sum4
  real(wp), dimension(3), intent(out) :: dHdr
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  integer j, jj

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
end subroutine calc_dHdr

! dHdr calculation trial method
subroutine calc_dHdr2(dHdr, Delta, r)
  implicit none

  complex(wp) :: sum1, sum2, sum3, sum4
  real(wp), dimension(3), intent(out) :: dHdr
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*i*Rabi(jj)*k(j,jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
     end do
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     sum3 = 0.0_wp
     do jj=1,num_beams
        sum3 = sum3 + (hbar/2.0_wp)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
     end do
     sum4 = 0.0_wp
     do jj=1,num_beams
        sum4 = sum4 + (-i)*k(j,jj)*conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdr(j) = 2.0_wp*real(sum1*sum2)+2.0_wp*real(sum3*sum4)
  end do
end subroutine calc_dHdr2

subroutine calc_dHdr3(dHdr, Delta, r)
  implicit none

  complex(wp) :: sum1, sum2, sum3, sum4
  real(wp), dimension(3), intent(out) :: dHdr
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*i*Rabi(jj)*k(j,jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
     end do
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdr(j) = 2.0_wp*real(sum1*sum2)
     sum3 = 0.0_wp
     do jj=1,num_beams
        sum3 = sum3 + (hbar/2.0_wp)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)
     end do
     sum4 = 0.0_wp
     do jj=1,num_beams
        sum4 = sum4 + (-i)*k(j,jj)*conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdr(j) = dHdr(j) + 2.0_wp*real(sum3*sum4)

  end do
end subroutine calc_dHdr3

subroutine calc_dHdr4(dHdr, Delta, r)
  implicit none

  complex(wp) :: sum1, sum2, sum3, sum4
  real(wp), dimension(3), intent(out) :: dHdr
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  integer j, jj
  complex(wp), dimension(6) :: expdot, expdot2



     do jj=1,num_beams
        expdot(jj) = (hbar/2.0_wp)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)                     
     end do

     do jj=1,num_beams
        expdot2(jj) = (-i)*conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do

     sum2 = 0.0_wp
     do jj=1,num_beams
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do

     sum3 = 0.0_wp
     do jj=1,num_beams
        sum3 = sum3 + expdot(jj)
     end do

  do j=1,3
     sum1 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + i*k(j,jj)*expdot(jj)
     end do

     sum4 = 0.0_wp
     do jj=1,num_beams
        sum4 = sum4 + k(j,jj)*expdot2(jj)
     end do
     dHdr(j) = 2.0_wp*real(sum1*sum2)+2.0_wp*real(sum3*sum4)
  end do
end subroutine calc_dHdr4




! dHdp calculation
subroutine calc_dHdp(dHdp, Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(3), intent(out) :: dHdp
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(j,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(j) = p(j)/m + 2.0_wp*real(sum1*sum2)
  end do

end subroutine calc_dHdp

subroutine calc_dHdp2(dHdp, Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(3), intent(out) :: dHdp
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(j,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
     end do
     do jj=1,num_beams
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(j) = p(j)/m + 2.0_wp*real(sum1*sum2)
  end do

end subroutine calc_dHdp2

subroutine calc_dHdp3(dHdp, Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(3), intent(out) :: dHdp
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer jj

     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(1,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(1) = p(1)/m + 2.0_wp*real(sum1*sum2)

     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(2,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(2) = p(2)/m + 2.0_wp*real(sum1*sum2)

     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        sum1 = sum1 + (hbar/2.0_wp)*(2.0_wp*k(3,jj)/m)*Rabi(jj)*exp(i*dot_product(k(:,jj),r))/(2.0_wp*Delta(jj)+i*Gamma)**(2.0_wp)
        sum2 = sum2 + conjg(Rabi(jj))*exp(-i*dot_product(k(:,jj),r))
     end do
     dHdp(3) = p(3)/m + 2.0_wp*real(sum1*sum2)

end subroutine calc_dHdp3

subroutine calc_dHdp4(dHdp, Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(3), intent(out) :: dHdp
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp

     sum1 = (hbar/2.0_wp)*(2.0_wp*k(j,1)/m)*Rabi(1)*exp(i*dot_product(k(:,1),r))/(2.0_wp*Delta(1)+i*Gamma)**(2.0_wp) &
     + (hbar/2.0_wp)*(2.0_wp*k(j,2)/m)*Rabi(2)*exp(i*dot_product(k(:,2),r))/(2.0_wp*Delta(2)+i*Gamma)**(2.0_wp) &
     + (hbar/2.0_wp)*(2.0_wp*k(j,3)/m)*Rabi(3)*exp(i*dot_product(k(:,3),r))/(2.0_wp*Delta(3)+i*Gamma)**(2.0_wp) &
     + (hbar/2.0_wp)*(2.0_wp*k(j,4)/m)*Rabi(4)*exp(i*dot_product(k(:,4),r))/(2.0_wp*Delta(4)+i*Gamma)**(2.0_wp) &
     + (hbar/2.0_wp)*(2.0_wp*k(j,5)/m)*Rabi(5)*exp(i*dot_product(k(:,5),r))/(2.0_wp*Delta(5)+i*Gamma)**(2.0_wp) &
     + (hbar/2.0_wp)*(2.0_wp*k(j,6)/m)*Rabi(6)*exp(i*dot_product(k(:,6),r))/(2.0_wp*Delta(6)+i*Gamma)**(2.0_wp) 

     sum2 = conjg(Rabi(1))*exp(-i*dot_product(k(:,1),r)) &
     + conjg(Rabi(2))*exp(-i*dot_product(k(:,2),r)) &
     + conjg(Rabi(3))*exp(-i*dot_product(k(:,3),r)) &
     + conjg(Rabi(4))*exp(-i*dot_product(k(:,4),r)) &
     + conjg(Rabi(5))*exp(-i*dot_product(k(:,5),r)) &
     + conjg(Rabi(6))*exp(-i*dot_product(k(:,6),r)) 

     dHdp(j) = p(j)/m + 2.0_wp*real(sum1*sum2)
  end do

end subroutine calc_dHdp4


! H calculation
subroutine calc_H(Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer j

  sum1 = 0.0_wp
  sum2 = 0.0_wp
  do j=1,num_beams
     sum1 = sum1 + (hbar/2.0_wp)*Rabi(j)*exp(i*dot_product(k(:,j),r))/(2.0_wp*Delta(j)+i*Gamma)
     sum2 = sum2 + conjg(Rabi(j))*exp(-i*dot_product(k(:,j),r))
  end do
  H = dot_product(p,p)/(2.0_wp*m) + 2.0_wp*real(sum1*sum2)

end subroutine calc_H


! hasn't been changed yet
subroutine calc_H_dipole(Delta, r, p)
  implicit none

  complex(wp) :: sum1, sum2
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  real(wp), dimension(3), intent(in) :: p
  integer j

  sum1 = 0.0_wp

  do j=1,num_beams
     sum1 = sum1 + Rabi(j)*exp(i*dot_product(k(:,j),r))
  end do
  
  H = dot_product(p,p)/(2.0_wp*m) + (2.0_wp*hbar*Delta(1))/(4.0_wp*Delta(1)*Delta(1)+Gamma*Gamma)*sum1*conjg(sum1)

end subroutine calc_H_dipole

subroutine calc_dHdr_dipole(dHdr, Delta, r)
  implicit none

  complex(wp) :: sum1, sum2, temp1
  real(wp), dimension(3), intent(out) :: dHdr
  real(wp), dimension(num_beams), intent(in) :: Delta
  real(wp), dimension(3), intent(in) :: r
  integer j, jj

  do j=1,3
     sum1 = 0.0_wp
     sum2 = 0.0_wp
     do jj=1,num_beams
        temp1 = Rabi(jj)*exp(i*dot_product(k(:,jj),r))
        sum1 = sum1 + temp1
        sum2 = sum2 + i*k(j,jj)*temp1
     end do
     dHdr(j) = (2.0_wp*hbar*Delta(1))/(4.0_wp*Delta(1)*Delta(1)+Gamma*Gamma)&
     *(sum1*conjg(sum2)+sum2*conjg(sum1))
!     write(0,*) sum1,sum2,dHdr
  end do
end subroutine calc_dHdr_dipole

subroutine calc_dHdp_dipole(dHdp, p)
  implicit none

  real(wp), dimension(3), intent(out) :: dHdp
  real(wp), dimension(3), intent(in) :: p
  integer j

  do j=1,3
     dHdp(j) = p(j)/m
  end do

end subroutine calc_dHdp_dipole



! Delta calculation
subroutine calc_Delta(Delta, p)
  implicit none

  real(wp), dimension(num_beams), intent(out) :: Delta
  real(wp), dimension(3), intent(in) :: p
  integer j

  do j=1,num_beams
     Delta(j) = Delta_base(j) - dot_product(k(:,j),p)/m - (hbar * dot_product(k(:,j), k(:,j)) / (2.0_wp*m))
  end do
end subroutine calc_Delta


subroutine calc_Delta_dipole(Delta)
  implicit none

  real(wp), dimension(num_beams), intent(out) :: Delta
  integer j

  do j=1,num_beams
     Delta(j) = Delta_base(j) - (hbar * dot_product(k(:,j), k(:,j)) / (2.0_wp*m))
  end do
end subroutine calc_Delta_dipole




end module odeab_support
