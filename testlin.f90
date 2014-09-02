!-----------------------------------------------------------------------
!
!  testlin.f90
!
!  Simple test program for the linsol module for solutions of equations.
!  
!  The 'globals' module defines the real/complex kind.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! This program caclulates the steady state populations and coherences of an
! two level atom excited by n laser fields. The Hamiltonian and the
! unconditioned master equation give us the equations of motion of the density
! operator.
! 
! For this program, we write the equations of motion in the form
! (d/dt)rho=rhodot
! where rho is a vector of all rho_ij (where i,j vary from 1 to n+1), rhodot is
! the time derivative vector, and (d/dt) is the linear derivative operator, a
! two dimensional (n+1)x(n+1) matrix. Each row is how drho_ij is calculated in
! terms of a linear combination of all components of the rho vector.
!
! Since rho has two indices, but we are writing it as a vector, we must be
! clear how we enumerate it. We can go from the double index to the single
! index by the convention  rho_ij <-> rho_ind where ind=i+j*(dim-1)
! 
! To make it easier to code, the contributions from each part of the
! Hamiltonian is split up into HAF, HA, DECAY, and NORMALIZE parts. The parts
! of the (d/dt) matrix are calculated by the calc_* functions, where the ir,ic
! inputs determine which rho_ij we are calculating the time derivative of.
! That is, calc_*(i,j,*) calculates the row in (d/dt) corresponding to
! (d/dt)rho_ij.
!
! A population constraint popconstraint replaces one of the population
! derivatives. This row should be populated such that the ijth element
! is equal to delta(ij).
!
!-----------------------------------------------------------------------


module testlin

use globals
use linsol

implicit none

integer, parameter :: n = num_beams
integer, parameter :: dim = n + 1
integer, parameter :: dim2 = dim*dim


contains

subroutine calc_populations(rho, Delta)

  complex(wp), dimension(dim2,dim2) :: ddt
  complex(wp), dimension(dim2)   :: rho
  complex(wp), dimension(dim2)   :: ddtHAF, ddtHA, ddtDECAY, ddtNORMALIZE, popconstraint, rhodot
  real(wp), dimension(num_beams) :: Delta
  integer :: sindex

  real(wp) :: rval = 0
  integer :: icount, count, count_rate
  integer :: ir, ic

  ! So, I'm trying to find the steady state values of rho. I currently have
  ! the position and momentum of the atom. I can write down the master equation
  ! which is essentially the matrix (d/dt) of the equation (d/dt)rho = rhodot. I remove one of the
  ! populations, replace it with a constraint on the populations, set rhodot=0
  ! except for the constraint on the populations. Then, I find the 


  do ir=1,dim
     do ic=1,dim

        ddtHAF = 0.0_wp
        call calc_ddtHAF(ir,ic,ddtHAF)

        ddtHA = 0.0_wp
        call calc_ddtHA(ir,ic,Delta,ddtHA)

        ddtDECAY = 0.0_wp
        call calc_ddtDECAY(ir,ic,ddtDECAY)
        
        ddtNORMALIZE = 0.0_wp
        call calc_ddtNORMALIZE(ir,ic,ddtNORMALIZE)
        
        ddt(ind(ir,ic),:) = ddtHAF + ddtHA + ddtDECAY + ddtNORMALIZE
        
     end do
  end do
  

call calc_constraint(popconstraint)

ddt(1,:) = popconstraint

rhodot = 0.0_wp
rhodot(1) = 1.0_wp

! test vector linsolve
call system_clock(icount, count_rate)
rho = linsolve(ddt, rhodot, resid=rval, tol=0._wp)
call system_clock(count, count_rate)

print *, 'Time to solve system of equations: ', real(count-icount)/real(count_rate)

! compute residual
write(*,*) 'max computed residual: ', rval
write(*,*) 'max check residual: ', maxval(abs(matmul(ddt,rho)-rhodot))
write(*,*) '(machine epsilon * N = ', epsilon(1.0_wp) * N, ')'
write(*,*)

end subroutine calc_populations


!---------------

function ind(ir,ic)
  integer :: ind
  integer, intent(in) :: ir, ic
  ind = ir+(ic-1)*dim
  return
end function ind

!---------------

subroutine calc_ddtHAF(ir,ic,array)

  complex(wp), dimension(dim2), intent(out) :: array
  integer, intent(in) :: ir, ic
  integer :: ibeam

  if ( ic .ne. 1 ) then
     array(ind(ir,1)) = i*0.5_wp*Rabi(ic-1)
  end if

  if ( ir .ne. 1 ) then
     array(ind(1,ic)) = -i*0.5_wp*conjg(Rabi(ir-1))
  end if

  do ibeam=1,n

     if ( ir .eq. 1 ) then
        array(ind(ibeam+1,ic)) = -i*0.5_wp*Rabi(ibeam)
     end if

     if ( ic .eq. 1 ) then
        array(ind(ir,ibeam+1)) = i*0.5_wp*conjg(Rabi(ibeam))
     end if

  end do

end subroutine calc_ddtHAF

!---------------

subroutine calc_ddtHA(ir,ic,Delta,array)
  complex(wp), dimension(dim), intent(out) :: array
  real(wp), dimension(num_beams), intent(in) :: Delta
  integer, intent(in) :: ir, ic
  complex(wp) :: sum 
  integer :: irho, isum

  if ( ir .gt. 1 .and. ic .eq. 1 ) then
     array(ind(ir,ic)) = -i*Delta(ir-1)
  end if

  if ( ir .eq. 1 .and. ic .gt. 1 ) then
     array(ind(ir,ic)) = i*Delta(ic-1)
  end if

  if ( ir .gt. 1 .and. ic .gt. 1 ) then
     array(ind(ir,ic)) = i*(Delta(ic-1)-Delta(ir-1))
  end if

end subroutine calc_ddtHA

!---------------

subroutine calc_ddtDECAY(ir,ic,array)

  complex(wp), dimension(dim2), intent(out) :: array
  integer, intent(in) :: ir, ic


  ! is the folloing really correcy? Seems suspect to me.
  if ( ir .gt. 1 .and. ic .gt. 1 ) then
     array(ind(ir,ic)) = -Gamma
  end if

  if ( ir .gt. 1 .and. ic .eq. 1 ) then
     array(ind(ir,ic)) = -Gamma*0.5_wp
  end if

  if ( ir .eq. 1 .and. ic .gt. 1 ) then
     array(ind(ir,ic)) = -Gamma*0.5_wp
  end if

end subroutine calc_ddtDECAY

!---------------

subroutine calc_ddtNORMALIZE(ir,ic,array)

  complex(wp), dimension(dim2), intent(out) :: array
  integer, intent(in) :: ir, ic
  integer :: ibeam

  if ( ir .eq. 1 .and. ic .eq. 1 ) then
     do ibeam=1,n
        array(ind(ibeam+1,ibeam+1))=Gamma
     end do
  end if

end subroutine calc_ddtNORMALIZE

!---------------

subroutine calc_constraint(array)

  complex(wp), dimension(dim2), intent(out) :: array
  complex(wp) :: sum 
  integer :: iindex

  ! just an identity matrix
  do iindex=1,dim
     array(ind(iindex,iindex)) = 1
  end do

end subroutine calc_constraint

!---------------

end module testlin
