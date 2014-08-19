!-----------------------------------------------------------------------
! 
!  mot.f90
!
!  Code for solving an atom in a basic optical lattice. Can
!  do basic cooling!
! 
!
!  Updates:
!
!  8/14/2014 - Decided it's about time to add a changelog to this thing.
!              Currently, it is a simple optical molasses simulation
!              that does classical hamiltonian evolution with added
!              momentum kicks for photon absorbtion and emission.
!
!  8/14/2014 - Added the function for monitoring diffusion in and out
!              of a particular radius. See escape_state.
!
!  8/15/2014 - Changed the method for picking a random emission direction
!              in the hopes to make things run faster. I think it has a
!              marginal improvement. I think the main slowdown is when
!              the odeab integrator is reset after each jump.
!
!  8/15/2014 - Added an rk4 solver for stepping (instead of the adaptive
!              stepping using odeab90). This gives qualitatively similar
!              evolution. Unfortunately, it doesn't seem to be faster, as
!              I had been hoping.
!
!-----------------------------------------------------------------------


program mot

  ! we are using variables and subroutines in these other "modules,"
  ! defined in the files of the same name, ending in .f90
  ! (where the modules are declared)
  ! so we notify the compiler to load that information
  ! (this is one of the nicest features of Fortran 90)
  use globals
  use random_pl
  use odeab_support
  use odeab90
  use rk4
  use utilities

  implicit none

  ! declare storage for solution vector; this is specific to the odeab
  ! integrator-- only modify if you want to use real variables instead
  ! of complex (change complex -> real)
  real(wp), dimension(:,:), pointer :: y

  ! declare variables to use below
  real(wp)       :: t = 0
  character(256) :: format1, buff
  integer        :: l, nout, j

  ! declare command-line arguments:
  ! Rabi                (in globals.f90)  -> lattice Rabi frequency
  ! Delta_base          (in globals.f90)  -> base detuning omega-omega0
  ! kval                (in globals.f90)  -> lattice wavelength
  ! Gamma               (in globals.f90)  -> spontaneous emission rate
  ! v0                  (in globals.f90)  -> velocity at t=0
  ! focalradius         (in globals.f90)  -> Radius of where photons are collected from trapped atoms
  ! focalthreshhold     (in globals.f90)  -> Atoms must travel at LEAST this far before they can count as a return
  ! motradius           (in globals.f90)  -> Radius which atoms can be declared as escaped from the MOT (simulation stops)
  real(wp) :: tstep, tfinal        !      -> time steps and final time for integrator
  real(wp), dimension(3,2) :: y0   !      -> initial conditions storage
  integer :: save_interval         !      -> Save every N data points
  integer(rpk) :: seed             !      -> Seed for RNG


  ! declare variables for spontaneous emission:
  real(wp), dimension(num_beams) :: Delta
  real(wp), dimension(num_beams) :: sds_array
  real(wp), dimension(3)         :: pkick
  real(wp), dimension(3)         :: pkick_avg
  real(wp), dimension(3)         :: pabs_avg
  real(wp), dimension(6)         :: pabs_count
  real(wp)                       :: sds_sum
  real(wp)                       :: rand_mode
  real(wp)                       :: pmag
  real(wp)                       :: phi
  real(wp)                       :: theta
  integer                        :: num_emissions
  integer                        :: spontaneous_emission_flag = 1

  ! declare variables for return/escape probabilities
  integer :: escape_state = 0  ! 0: initial state, atom still in focalradius
                               ! 1: atom has passed beyond theshhold and is diffusing
                               ! 2: atom has returned to focalradius after being in state 1 (end simulation)
                               ! 3: atom has escaped beyond the motradius                   (end simulation)




  integer :: min_arguments = 13 !

  pkick_avg(1) = 0.0_wp
  pkick_avg(2) = 0.0_wp
  pkick_avg(3) = 0.0_wp
  pabs_avg(1) = 0.0_wp
  pabs_avg(2) = 0.0_wp
  pabs_avg(3) = 0.0_wp


  !!!!!! Get arguments from the command line
  if ( command_argument_count() .ne. min_arguments .and. command_argument_count() .ne. (min_arguments + 1) )  call usage()
  call get_command_argument(1, buff)
  Rabi(1) = s2r(buff)
  do l=1,num_beams
     Rabi(l)=Rabi(1)
  end do

  call get_command_argument(2, buff)
  Delta_base(1) = s2r(buff)
  do l=2,num_beams
     Delta_base(l)=Delta_base(1)
  end do

  call get_command_argument(3, buff)
  kval = s2r(buff)
  !!
  k(1,1) = kval
  k(2,1) = 0
  k(3,1) = 0
  !!
  k(1,2) = 0
  k(2,2) = kval
  k(3,2) = 0
  !!
  k(1,3) = 0
  k(2,3) = 0
  k(3,3) = kval
  !!
  k(1,4) = -kval
  k(2,4) = 0
  k(3,4) = 0
  !!
  k(1,5) = 0
  k(2,5) = -kval
  k(3,5) = 0
  !!
  k(1,6) = 0
  k(2,6) = 0
  k(3,6) = -kval

  call get_command_argument(4, buff)
  Gamma = s2r(buff)

  call get_command_argument(5, buff)
  p0(1) = s2r(buff)

  call get_command_argument(6, buff)
  p0(2) = s2r(buff)

  call get_command_argument(7, buff)
  p0(3) = s2r(buff)

  call get_command_argument(8, buff)
  focalradius = s2r(buff)

  call get_command_argument(9, buff)
  focalthreshhold = s2r(buff)

  call get_command_argument(10, buff)
  motradius = s2r(buff)

  call get_command_argument(11, buff)
  tstep = s2r(buff)
  if ( tstep .le. 0.0_wp ) then 
    write(0,*)  "Error: inappropriate tstep"
    call usage()
  end if

  call get_command_argument(12, buff)
  tfinal = s2r(buff)
  if ( tfinal .le. 0.0_wp .or. tfinal .lt. tstep ) then 
    write(0,*)  "Error: inappropriate tfinal"
    call usage()
  end if

  nout = nint(tfinal/tstep)
  if ( nout .gt. 50000000 ) then
    write(0,*)  "Error: too many steps, are you crrrrazy??!"
    call usage()
  end if

  call get_command_argument(13, buff)
  save_interval = s2r(buff)

  ! get random number seed, if specified
  if ( command_argument_count() .eq. (min_arguments + 1) ) then
    call get_command_argument(14, buff)
    seed = s2i(buff)
    call init_rand_pl(seed1=seed)
  end if


  ! the integrator needs you to do this (just do it)
  ! note: this storage is used for both odeab90 AND rk4 methods 
  y => odeab_y

  ! initialize the rk4 if we are using it
  if( step_method .eq. 1 ) then
     call rk4init()
  end if


  ! Set initial values
  r0(1) = 0.0_wp ! 128165.0_wp
  r0(2) = 0.0_wp
  r0(3) = 0.0_wp
  y = 0.0_wp
  y(:,1) = r0
  y(:,2) = p0

  ! Output initial conditions
  format1 = "(A,d10.3,A,d10.3,A,d10.3,A,d10.3,A,d10.3,A,d10.3,A,d10.3)"
  !write(*,format1) '#Rb=', real(Rabi(1)), 'D=', Delta_base(1), 'k=', kval, 'G=', Gamma, 'px=', p0(1), 'py=', p0(2), 'pz=', p0(3)
  call get_command(buff)
  write(*,*) "#", buff

  

  write(*,*) t, real(y(1,1)), real(y(2,1)), real(y(3,1)), real(y(1,2)), real(y(2,2)), real(y(3,2)), H, escape_state


  !!!!!! Main integration loop
  do l = 1, nout
    ! write(0,*) "l=", l, "out of nout=", nout
    ! this advances the solution array 'y' from t to t+tstep
    ! note that the variable 't' gets automatically updated to t+tstep

     if ( step_method .eq. 0 ) then
        call odeab(t, t+tstep)
     else if ( step_method .eq. 1) then
        call rk4step(y, t, tstep)
        t = t+tstep
     else if ( step_method .eq. 2) then
        call eulerstep(y, t, tstep)
        t = t+tstep
     end if

    ! Spontaneous Emission 
    ! First calculate <sigma^dagger sigma> (<SDS>), and while were at it, we store the running
    ! sum into sds_array. This will be used only if the decay occurs, but might as well.
    if ( spontaneous_emission_flag .eq. 1 ) then
    sds_sum = 0.0_wp
    do j=1,num_beams
       ! Detuning Matrix Construction
       ! (we could make Delta global, as it is computed in odeab_support )
       Delta(j) = Delta_base(j) - dot_product(k(:,j),y(:,2))/m - (hbar * dot_product(k(:,j), k(:,j)) / (2.0_wp*m))

       !Sum and running sum of SDS:
       sds_sum = sds_sum + Rabi(j)*conjg(Rabi(j))/(Gamma**2+4*Delta(j)**2)
       sds_array(j) = sds_sum
    end do
  
    ! Reset trajectory and give momentum kicks if spontaneous emission occurs:
    if ( rand_pl() .lt. Gamma*sds_sum*tstep ) then           
       num_emissions = num_emissions + 1
       !write(0,*) "SPONTANEOUS EMISSION"
       ! Atom decayed, so it must have absorbed a photon from one of the
       ! trapping beams. The probability is proportional to that beam's
       ! contribution to the <SDS> sum. So we normalize the running sum
       ! and pick a random number to bin the most likely momentum state.

       sds_array = sds_array/sds_sum
       rand_mode = rand_pl()

       !write(0,*) "rand_mode = ", rand_mode
       do j=1,num_beams
          if ( rand_mode .lt. sds_array(j) ) then
             y(:,2) = y(:,2) + hbar*k(:,j)
             !write(0,*) "ABSORBTION DIRECTION KICK = ", (hbar*k(:,j))
             !write(0,*) "sds_array(j) = ", j, sds_array(j)
             pabs_avg = pabs_avg + k(:,j)
             pabs_count(j) = pabs_count(j) + 1
             exit
          end if
       end do


       ! And the atom also emitted a photon, so we pick a random direction,
       ! give atom a single photon momentum kick. Maybe use dipole emission
       ! pattern in the future, but start with uniform direction:


!       phi = 2.0_wp*pi*rand_pl()
!       theta = acos(2.0_wp*rand_pl()-1.0_wp)
!       pkick(1) = cos(phi)*sin(theta)
!       pkick(2) = sin(phi)*sin(theta)
!       pkick(3) = cos(theta)

       do
          pkick(1) = (rand_pl()*2.0_wp)-1.0_wp
          pkick(2) = (rand_pl()*2.0_wp)-1.0_wp
          pkick(3) = (rand_pl()*2.0_wp)-1.0_wp
          pmag = dot_product(pkick,pkick)

          if ( pmag < 1.0_wp ) then
             pmag = sqrt(pmag)
             pkick(1) = pkick(1)/pmag
             pkick(2) = pkick(2)/pmag
             pkick(3) = pkick(3)/pmag
             exit
          end if

       end do

       pkick_avg = pkick_avg + pkick
       y(:,2) = y(:,2) + hbar*kval*pkick
       !write(0,*) "EMISSION DIRECTION KICK = ", (hbar*kval*pkick)

       ! Reset the integrator since it doesn't appriciate us adding discontinuities.
       ! If we don't do this, it basically ignores whatever changes we make! (This
       ! essentially will make it forget the history of all previous steps).
       odeab_istate = 1
    end if
    end if

    ! Do escape/return probability checks
    if ( escape_state .eq. 0 .and. dot_product(y(:,1),y(:,1)) .gt. focalthreshhold**2  ) then
       write(0,*) "escape_state = 1; atom has exited threshhold and is diffusing! "     
       escape_state = 1;
    else if ( escape_state .eq. 1 .and. dot_product(y(:,1),y(:,1)) .lt. focalradius**2  ) then
       write(0,*) "escape_state = 2; atom has returned to the focal volume! (stopping...) "     
       escape_state = 2;
    else if ( escape_state .eq. 1 .and. dot_product(y(:,1),y(:,1)) .gt. motradius**2  ) then
       write(0,*) "escape_state = 3; atom has escaped from the mot volume! (stopping...) "     
       escape_state = 3;
    end if

    ! Checking For Errors while integrating.
    ! We shouldn't do this, or we should only do it BEFORE we modify odeab_istate
    ! call handle_odeab_error()

    ! Calculate exact solution and output results; needs reformatting?
    call calc_H( Delta, y(:,1), y(:,2) )
    if ( mod( l, save_interval ) .eq. 0 .or. escape_state .eq. 2 .or. escape_state .eq. 3 ) then
       write(*,*) t, real(y(1,1)), real(y(2,1)), real(y(3,1)), real(y(1,2)), real(y(2,2)), real(y(3,2)), H, escape_state
    end if

    ! end the simulation if the atom has returned or escaped
    if ( escape_state .eq. 2 .or. escape_state .eq. 3 ) then
       write(0,*) "Atom has returned or escaped, ending simulation"
       exit
    end if

  end do

  ! print performance statistics, if you're curious, to standard error
  call print_odeab_stats(cumulative = .true.)
  write(0,*) 'MOT Stats:'
  write(0,*) 'num_emissions: ', num_emissions
  write(0,*) 'avg emissions per step (if too big, decrease step size!): ', num_emissions*1.0_wp/(tfinal/tstep)
  write(0,*) 'avg emission direction: ', pkick_avg/(1.0_wp*num_emissions)
  write(0,*) 'avg absorbtion direction: ', pabs_avg/(1.0_wp*num_emissions)
  write(0,*) 'absorbtion direction counts: ', pabs_count/(1.0_wp*num_emissions)

  contains

  subroutine usage()
    ! write (0,*) means write to standard error (0), using default
    !   formatting (*)
    write(0,*) ''
    write(0,*) ''
    write(0,*) 'Usage: mot <Rabi> <Delta_base> <kval> <Gamma> <px> <py> <pz> <tstep> <tfinal> [seed]'
    write(0,*) ' Rabi        ->  scaled Rabi frequency'
    write(0,*) ' Delta_base  ->  base detuning = omega-omega0'
    write(0,*) ' kval        ->  magnitude of k vectors'
    write(0,*) ' Gamma       ->  spontaneous emission rate'
    write(0,*) ' px          ->  px at t=0'
    write(0,*) ' py          ->  py at t=0'
    write(0,*) ' pz          ->  pz at t=0'
    write(0,*) ' focalradius ->  radius of where photons are collected from trapped atoms'
    write(0,*) ' focalthreshhold -> radius beyond atoms must travel before being counted'
    write(0,*) ' motradius   ->  radius where atoms escape and simulation ends'
    write(0,*) ' tstep       ->  time step for integrator output (important for spontaneous emission)'
    write(0,*) ' tfinal      ->  final time of integrator output'
    write(0,*) ' seed        ->  (optional) set the random number generator seed'
    write(0,*) ''
    write(0,*) ' Units are in terms of the Rb characteristic constants:'
    write(0,*) ' -Atomic mass, Spontaneous Emission Rate, and D2 transition wavelength.'
    write(0,*) ' Example invocation: ./mot 0.5 -1 1 1 0.0259 0.0 0.0 .10 50000 > out'
    write(0,*) ''
    write(0,*) 'Output is 7 columns to standard output, except the first'
    write(0,*) '  scaled time, x, y, z, px, py, pz, H'
    write(0,*) ''
    write(0,*) ''
    stop
  end subroutine usage

end program mot

