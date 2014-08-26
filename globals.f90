!-----------------------------------------------------------------------
!
!  Globals module.  Put all variables to be accessible everywhere
!    here. 
!
!-----------------------------------------------------------------------

module globals

  integer, parameter     :: wp = selected_real_kind(p=14)   ! this sets the floating-point precision to double
  complex(wp), parameter :: i = (0.0_wp, 1.0_wp)            ! sqrt (-1)
  integer, parameter     :: num_beams = 6                   ! Number of trapping beams
  real(wp), parameter    :: pi = 3.14159265358979_wp        ! pi! REAL: VERY WEIRD BUGS HAPPEN IF YOU SET PI=3

  ! flag to choose hamiltonian
  integer, parameter :: hamiltonian_version = 1   ! 0: use the full effective hamiltonian
                                                  ! 1: use the large delta hamiltonian (real only)
  ! flag to choose integration method
  integer, parameter :: step_method = 0           ! 0: use the adaptive stepping integrator odeab90 (default)
                                                  ! 1: use the rk4 integrator
                                                  ! 2: use the euler integrator
                                                  ! 3: use the verlet step integrator

  ! Constants in alternative unit system for Rb:
  ! - mass measured in atomic mass of Rb (~10^-25 kg)
  ! - length measured in wavelength of D2 Transition (~780 nm)
  ! - time measured in the natural decay rate (~10^6 s^-1)
  real(wp) :: m          = 1.0_wp               ! Atomic mass of Rb
  real(wp) :: lambda0    = 1.0_wp               ! D2 transition wavelength
  real(wp) :: Gamma      = 1.0_wp               ! Spontaneous emission rate
  real(wp) :: hbar       = 0.0000314909_wp      ! Reduced Plank's constant in m*Gamma*lambda0^2
  real(wp) :: photon_detection_efficiency = 0.0007510 ! Detection efficiency due to solid angle and detector efficiency.
!  real(wp) :: photon_detection_efficiency = 0.7510 ! Detection efficiency due to solid angle and detector efficiency.  1000 times better than what we really have

  real(wp) :: kval                              ! k-vector magnitude
  real(wp), dimension(3)   :: r0                ! position at t=0
  real(wp), dimension(3)   :: p0                ! momentum at t=0
  real(wp), dimension(3,2) :: rp0               ! initial position and momentum matrix at t=0

  real(wp) :: focalradius                       ! Radius of where photons are collected from trapped atoms
  real(wp) :: focalthreshhold                   ! Atoms must travel at LEAST this far before they can count as a return
  real(wp) :: motradius                         ! Radius which atoms can be declared as escaped from the MOT (simulation stops)

  real(wp) :: H                                       ! Hamiltonian
  complex(wp), dimension(num_beams)   :: Rabi         ! Rabi frequency for each beam
  real(wp), dimension(num_beams)      :: Delta_base   ! laser detuning from resonance (before momentum considerations).
  real(wp), dimension(3,num_beams)    :: k            ! k variable is a matrix of column k vectors, one for each beam.

  ! settings for the portable random number generator rand_pl
  integer, parameter :: rand_pl_wp = wp    ! get double-precision random #s
  integer, parameter :: rand_pl_mf = 100   ! use simple, fast generator

end module globals
