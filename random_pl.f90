!-----------------------------------------------------------------------
!
!  random_pl.f90
!
!  Version 2.0.6
!
!  Module of routines to generate uniform random numbers, implementing
!    3 different 32-bit portable algorithms, in order of increasing
!    period and complexity:
!   
!    1. The subtractive, lagged-Fibonacci generator discussed by 
!      Donald Knuth, _The Art of Computer Programming, v. 2: 
!      Seminumerical Algorithms_, 3rd ed. (Addison-Wesley, Reading, 
!      1998). This generator has a period of 2**29 * (2**100 - 1), or 
!      10**38.  Different seeds give different outputs for at least the 
!      first 2**70 numbers.  (These numbers go down by a factor of 10 
!      when skipping numbers to improve randomness.)  Initialization
!      is according to Knuth's revised recommendations in the 9th
!      printing of his book, making the number robust to restarting
!      often with a new seed.  Also, the implementation here uses
!      the improvement by M. Luescher of generating 1009 numbers
!      at a time but only using the first 100.  The range of the 
!      integer output of this routine is [0, 2**30 - 1].  
!
!    2. Pierre L'Ecuyer's combined multiple recursive generator
!      (Operations Research _44_, 816 (1996)), which has a period
!      of about 2**185.  Distinct seeds will cover adjacent sections
!      of the full sequence (by 'jumping ahead' in the sequence in
!      the initialization routine); thus sequences from distinct 
!      seeds will not overlap for at least about 2**152 deviates
!      (with minimum length for adjacent seeds).  The corresponding
!      integer output range of this routine is [0, 2**31 - 2].
!
!    3. The 'Mersenne Twister' of Matsumoto and Nishimura 
!      (ACM Trans. on Modeling and Comp. Sim. _8_, 3 (1998),
!      with an incredible period of 2**19937 - 1.  The routine
!      is initialized using the seed initialization in the updated
!      version 'mt19937ar.c', with no special effort to prevent 
!      overlapping sequences for distinct seeds (of which there 
!      are about 2**32); however, the period is so long that 
!      overlap is very unlikely.  The corresponding integer
!      output range of this routine is [0, 2**31 - 1] (the sign
!      bit is deleted in this implementation).
!     
!
!  The outputs of the two generators can be combined in several ways;
!    This must be selected via an integer parameter "rand_pl_mf", which 
!    has three digits, with the following meaning:
!
!    digit 1: select which generator is "primary output", i.e., which
!             one nominally generates the random deviates
!           allowed values: 1 -> lagged Fibonacci generator
!                           2 -> L'Ecuyer multiple recursive generator
!                           3 -> Mersenne Twister
!
!    digit 2: select which generator to subtractively mix with primary
!             output, if any, before entering shuffle table
!           allowed values: 0 -> none (primary generator is unmodified)
!                           1 -> lagged Fibonacci generator
!                           2 -> L'Ecuyer multiple recursive generator
!                           3 -> Mersenne Twister
!               constraint: this digit must be distinct from the first 
!                           digit
!
!    digit 3: select which generator manages shuffle table
!           allowed values: 0 -> none (no shuffling implemented)
!                           1 -> lagged Fibonacci generator
!                           2 -> L'Ecuyer multiple recursive generator
!                           3 -> Mersenne Twister
!                           4 -> Bays-Durham shuffle (shuffles output
!                                without an extra random sequence)
!               constraint: if nonzero, this digit must be distinct from 
!                           the first two digits
!
!  Some examples of selected combinations:
!    1. "economical generators:"
!    rand_pl_mf = 100 -> simple lagged Fibonacci generator
!    rand_pl_mf = 200 -> simple L'Ecuyer CMRG generator
!    rand_pl_mf = 300 -> simple Mersenne twister generator
!
!    2. "economical shuffled generators:"
!    rand_pl_mf = 104 -> lagged Fibonacci generator w/ Bays-Durham shuffle
!    rand_pl_mf = 204 -> L'Ecuyer CMRG generator w/ Bays-Durham shuffle
!    rand_pl_mf = 304 -> Mersenne twister generator w/ Bays-Durham shuffle
!
!    3. "subtractively combined generators:"
!    rand_pl_mf = 120 -> mix lagged Fibonacci w/ L'Ecuyer outputs
!    rand_pl_mf = 130 -> mix lagged Fibonacci w/ Mersenne Twister outputs
!    rand_pl_mf = 230 -> mix L'Ecuyer w/ Mersenne Twister outputs
!
!    4. "shuffle-combined generators:"
!    rand_pl_mf = 102 -> shuffle lagged Fibonacci w/ L'Ecuyer output
!    rand_pl_mf = 103 -> shuffle lagged Fibonacci w/ Mersenne Twister output
!    rand_pl_mf = 203 -> shuffle L'Ecuyer w/ Mersenne Twister output
!    rand_pl_mf = 302 -> shuffle Mersenne Twister w/ L'Ecuyer output
!
!    5. "complicated generator:"
!    rand_pl_mf = 231 -> mix L'Ecuyer output with Mersenne Twister output,
!                        shuffle result with Fibonacci generator
!   
!    There are 3*(2*3 + 1*4) = 30 different combinations in all.  The
!    three primary generators are of a sufficiently different character
!    that switching among possibilities in an application should give
!    confidence that things are working well.  A reasonable method flag
!    for most applications would be 301.
!
!  Following is a representative list of combinations; all have 
!    passed Marsaglia's 'Diehard' test battery 
!    (http://stat.fsu.edu/~geo/diehard.html), and are given 
!    with the time taken to generate 100 million deviates on a 2.2 
!    GHz Pentium IV Xeon, code compiled with ifc, as well as the
!    final deviate (both as reported by the 'diehard_data' example), 
!    to validate the operation of the generator:
!
!      rand_pl_mf     time  deviate  comment
!      ----------    -----  -------  -------
!             100    4.3 s  0.85311  (raw generators)
!             200   12.0 s  0.13756
!             300    4.8 s  0.75958
!
!             104   12.8 s  0.84113  (self-shuffled generators)
!             204   23.4 s  0.42919
!             304   14.1 s  0.57834
!
!             102   23.4 s  0.91295  (one generator shuffles another)
!             203   25.9 s  0.75207
!             301   16.8 s  0.49370
!
!             120   14.3 s  0.71556  (subtractive combination of 2)
!             230   16.2 s  0.37797
!             310    7.5 s  0.90647
!
!             123   28.4 s  0.36292  (combinations of all 3)
!             231   28.6 s  0.62976
!             312   26.7 s  0.36917
!
!
!  The interface routines are as follows:
!
!    rand_pl -> function returning a pseudorandom real(rand_pl_wp) 
!                   number in the range (0, 1), endpoints excluded;
!                   call as result = rand_pl().
!
!    get_rand_pl_state -> subroutine to get the state of the generator;
!                   call as get_rand_pl_state(ivec), where ivec is an
!                   integer vector of type rpk and length 
!                   rand_pl_state_sz (defined by this module).
!                   Useful, e.g., for checkpointing.
!
!    set_rand_pl_state -> subroutine to restore the state of the generator;
!                   call as set_rand_pl_state(ivec), where ivec is an
!                   integer vector of type rpk and length 
!                   rand_pl_state_sz (defined by this module), and
!                   was obtained by a previous call to get_rand_pl_state.
!
!    init_rand_pl -> subroutine to initialize random number generator; 
!                   call as init_rand_pl(seed1=n1, seed2=n2, seed3=n3), 
!                   where the three arguments of type integer(rpk) are 
!                   seeds to start out the sequence.  The first seed
!                   initializes the Fibonacci generator, the second seed
!                   initializes the L'Ecuyer generator, and the last seed
!                   initializes the Mersenne Twister.  All three
!                   arguments are optional; missing arguments are set
!                   equal to the first argument, so that all three 
!                   generators can be initialized from a single number.
!                   If no arguments are present, the generators are
!                   initialized using default seeds.
!
!  Additionally, there are two interface routines for generating 
!    vectors of random numbers:
!
!    rand_pl_vec -> function returning a vector of random numbers
!                     in the same way as rand_pl; call as
!                     out_vector = rand_pl_vec(n),
!                     where n is the desired quantity of random numbers.
!    nrand_pl_vec -> returns a vector of standard normal deviates;
!                      call as out_vector = rand_pl_vec(n),
!                      where n is an *even* integer.
!
!  In addition to this single random-number generator, this module
!    has facilities for multiple, simultaneous, identical generators in
!    addition to this first one.  To initialize the multiple
!    generators, call the allocation/initialization subroutine:
!
!    allocate_rand_pl_multi -> subroutine to allocate space for
!                   and initialize an extra set of random number
!                   generators.  Call as allocate_rand_pl_multi(m, 
!                   seed1, seed2, seed3, no_init_single), where the
!                   last 4 arguments are optional.  Here m is the
!                   number of generators to initialize, the seed
!                   arguments are as described above for init_rand_pl,
!                   and no_init_single is a logical flag for whether
!                   or not to reinitialize the single generator.
!                   The m+1 generators (including the single 
!                   generator) are initialized as follows: the
!                   single generator is initialized with seeds
!                   seed1, seed2, seed3; to initialize the multiple
!                   generators, the sequential seeds are used
!                   (i.e., generator i uses seeds seed1+i, seed2+i,
!                   and seed3+i).  The generators are designed so
!                   that they will be independent when initialized
!                   in this way.  To initialize them with other
!                   seeds, use init_rand_pl_multi described below.
!
!  Then these generators are accessed by routines as in the single
!    generators above, with _multi appended to the routine names,
!    and the first argument is an integer (from 1 to m) that 
!    indicates which generator to use.
!
!    rand_pl_multi -> function returning a pseudorandom real(rand_pl_wp)
!                   number in the range (0, 1), endpoints excluded,
!                   call as result = rand_pl_multi(j) for integer j.
!                   
!    get_rand_pl_state_multi -> subroutine to get the state of one of 
!                   the multiple generators; call as
!                   get_rand_pl_state_multi(j, ivec), where ivec is an
!                   integer vector of type rpk and length  
!                   rand_pl_state_sz (defined by this module).
!                   Useful, e.g., for checkpointing.
!                   
!    set_rand_pl_state_multi -> subroutine to restore the state of one
!                   of the multiple generators; call as
!                   set_rand_pl_state_multi(j, ivec), where ivec is an
!                   integer vector of type rpk and length  
!                   rand_pl_state_sz (defined by this module), and
!                   was obtained by a previous call to 
!                   get_rand_pl_state_multi.
!                   
!    init_rand_pl_multi -> subroutine to initialize one of the multiple
!                   random number generators; call as
!                   init_rand_pl_multi(j, seed1=n1, seed2=n2, seed3=n3),
!                   where the three arguments of type integer(rpk) are 
!                   seeds to start out the sequence.  The first seed
!                   initializes the Fibonacci generator, the second seed
!                   initializes the L'Ecuyer generator, and the last seed
!                   initializes the Mersenne Twister.  All three
!                   arguments are optional; missing arguments are set
!                   equal to the first argument, so that all three 
!                   generators can be initialized from a single number.
!                   If no arguments are present, the generators are
!                   initialized using default seeds.
!
!  Of course, there are again two interface routines for generating 
!    vectors of random numbers from the multiple generators:
!                   
!    rand_pl_vec_multi -> function returning a vector of random numbers
!                     in the same way as rand_pl_multi; call as
!                     out_vector = rand_pl_vec_multi(j, n),
!                     where n is the desired quantity of random numbers.
!    nrand_pl_vec_multi -> returns a vector of standard normal deviates; 
!                      call as out_vector = rand_pl_vec(j, n),
!                      where n is an *even* integer.
!
!  This module also defines the integer type "rpk", which it expects
!    for all integers 
!    (e.g., call init_rand_pl(1873334_rpk, 2983_rpk, 58954_rpk)).
!
!  This module requires a "globals" module, which defines the 
!    following parameters (mentioned above); actual values here are 
!    examples only:
!
!    module globals
!
!    ! rand_pl_wp is the real number kind
!    ! rand_pl_mf selects which random number strategy to use
!    integer, parameter :: rand_pl_wp = selected_real_kind(p=14) 
!    integer, parameter :: rand_pl_mf = 11
!
!    end module globals
!
!-----------------------------------------------------------------------


module random_pl

  use globals, only : rand_pl_wp, rand_pl_mf

  private
  public :: rpk, rand_pl_state_sz
  public :: init_rand_pl, rand_pl
  public :: get_rand_pl_state, set_rand_pl_state
  public :: rand_pl_vec, nrand_pl_vec
  public :: allocate_rand_pl_multi
  public :: init_rand_pl_multi, rand_pl_multi
  public :: get_rand_pl_state_multi, set_rand_pl_state_multi
  public :: rand_pl_vec_multi, nrand_pl_vec_multi

  ! kind is big enough for 32-bit integer arithmetic
  integer, parameter :: rpk = selected_int_kind(9)

  ! rename real kind
  integer, parameter :: wp = rand_pl_wp

  ! default seeds
  integer(rpk), parameter :: dseed1 = 310952
  integer(rpk), parameter :: dseed2 = 314159
  integer(rpk), parameter :: dseed3 = 4357

  ! parse out digits of rand_pl_mf
  integer, parameter :: rand_pl_mf_pri = rand_pl_mf/100
  integer, parameter :: rand_pl_mf_mix = rand_pl_mf/10 - rand_pl_mf_pri*10
  integer, parameter :: rand_pl_mf_scr = rand_pl_mf - (rand_pl_mf/10)*10

  ! temporary storage for fibonacci generator
  integer(rpk), dimension(1009) :: aa

  ! define storage size parameters
  integer, parameter :: scramble_stor_sz = 64
  integer, parameter :: fib_stor_sz = 100
  integer, parameter :: lecuyer_stor_sz = 6
  integer, parameter :: mersenne_stor_sz = 624
  integer, parameter :: rand_pl_state_sz = &
                scramble_stor_sz+fib_stor_sz+lecuyer_stor_sz+mersenne_stor_sz+5

  ! define storage for the single generator
  integer(rpk), dimension(lecuyer_stor_sz), target :: lecuyer_stor_s
  integer(rpk), dimension(scramble_stor_sz), target :: scramble_stor_s
  integer(rpk), dimension(fib_stor_sz), target :: fib_stor_s
  integer(rpk), dimension(mersenne_stor_sz), target :: mersenne_stor_s
  integer(rpk), target :: initialized_s = 0
  integer(rpk), target :: fib_ctr_s, bd_lastout_s, mersenne_ptr_s

  ! define storage for multiple generators
  integer(rpk), allocatable, dimension(:,:), target :: lecuyer_stor_m
  integer(rpk), allocatable, dimension(:,:), target :: scramble_stor_m
  integer(rpk), allocatable, dimension(:,:), target :: fib_stor_m
  integer(rpk), allocatable, dimension(:,:), target :: mersenne_stor_m
  integer(rpk), allocatable, dimension(:), target :: initialized_m
  integer(rpk), allocatable, dimension(:), target :: fib_ctr_m, bd_lastout_m, &
                                                     mersenne_ptr_m

  ! define storage pointers
  integer(rpk), dimension(:), pointer :: lecuyer_stor
  integer(rpk), dimension(:), pointer :: scramble_stor
  integer(rpk), dimension(:), pointer :: fib_stor
  integer(rpk), dimension(:), pointer :: mersenne_stor
  integer(rpk), pointer :: initialized
  integer(rpk), pointer :: fib_ctr, bd_lastout, mersenne_ptr

  ! to keep track of where pointers are associated
  integer :: curr_ptr = -1


  contains


!=======================================================================
!
!  subroutine init_rand_pl(optional seed1, optional seed2,
!                          optional seed3) 
!                               
!  I: seed1,2,3  int(rpk)     optional seeds for initializing generator
!  
!  O: globals only
!  
!  Initializes all 3 random number generators, with the
!  specified seeds (1 for Fibonacci, 2 for L'Ecuyer, 3 for Mersenne).
!  If seeds are not specified, default values are used, except if
!  only some seeds are specified, the rest are set to the first 
!  specified.
!  
!=======================================================================

subroutine init_rand_pl(seed1, seed2, seed3)

  implicit none
  
  integer(rpk), optional, intent(in) :: seed1, seed2, seed3
  integer(rpk) :: lseed1, lseed2, lseed3, firstseed

  ! set seeds
  if ( present(seed1) .or. present(seed2) .or. present(seed3) ) then
    if ( present(seed3) ) firstseed = seed3
    if ( present(seed2) ) firstseed = seed2
    if ( present(seed1) ) firstseed = seed1
    lseed1 = firstseed; lseed2 = firstseed; lseed3 = firstseed
    if ( present(seed1) ) lseed1 = seed1
    if ( present(seed2) ) lseed2 = seed2
    if ( present(seed3) ) lseed3 = seed3
  else
    lseed1 = dseed1; lseed2 = dseed2; lseed3 = dseed3
  end if

  call associate_single()
  call init_rand_pl_main(lseed1, lseed2, lseed3)

end subroutine init_rand_pl


!=======================================================================
!
!  subroutine init_rand_pl_multi(m, optional seed1, optional seed2,
!                                optional seed3)
!
!  I: m          integer      selects which generator to initialize
!     seed1,2,3  int(rpk)     optional seeds for initializing generator
!
!  O: globals only
!
!  Analogous routine to init_rand_pl, to initialize one of the set
!  of multiple random number generators, specified by m.
!
!=======================================================================

subroutine init_rand_pl_multi(m, seed1, seed2, seed3)

  implicit none

  integer, intent(in) :: m
  integer(rpk), optional, intent(in) :: seed1, seed2, seed3
  integer(rpk) :: lseed1, lseed2, lseed3, firstseed

  ! set seeds
  if ( present(seed1) .or. present(seed2) .or. present(seed3) ) then
    if ( present(seed3) ) firstseed = seed3
    if ( present(seed2) ) firstseed = seed2
    if ( present(seed1) ) firstseed = seed1
    lseed1 = firstseed; lseed2 = firstseed; lseed3 = firstseed
    if ( present(seed1) ) lseed1 = seed1
    if ( present(seed2) ) lseed2 = seed2
    if ( present(seed3) ) lseed3 = seed3
  else
    lseed1 = dseed1+m; lseed2 = dseed2+m; lseed3 = dseed3+m
  end if

  call associate_multi(m)
  call init_rand_pl_main(lseed1, lseed2, lseed3)

end subroutine init_rand_pl_multi


!=======================================================================
!
!  subroutine allocate_rand_pl_multi(m, optional seed1, optional seed2,
!                              optional seed3, optional no_init_single)
!
!  I: m          integer      selects how many generators to allocate
!     seed1,2,3  int(rpk)     optional seeds for initializing generator
!     no_init_single  logical selects whether or not to initialize
!                             the single generator
!
!  O: globals only
!
!  Allocate storage for m multiple random number generators.  Also
!  initializes all generators, including the single generator
!  (unless no_init_single = .true.), as follows: the single
!  generator is set using the specified/default seeds, and
!  the multi generators 1 through m are set using the same seeds
!  incremented by the generator number.
!
!  Any existing set of multiple generators will be obliterated.
!  
!=======================================================================

subroutine allocate_rand_pl_multi(m, seed1, seed2, seed3, no_init_single)

  implicit none

  integer, intent(in) :: m
  integer(rpk), optional, intent(in) :: seed1, seed2, seed3
  logical, optional, intent(in) :: no_init_single
  integer(rpk) :: lseed1, lseed2, lseed3, firstseed
  logical :: init_single
  integer :: j
  integer(rpk) :: jrpk
  
  ! set seeds
  if ( present(seed1) .or. present(seed2) .or. present(seed3) ) then
    if ( present(seed3) ) firstseed = seed3
    if ( present(seed2) ) firstseed = seed2
    if ( present(seed1) ) firstseed = seed1
    lseed1 = firstseed; lseed2 = firstseed; lseed3 = firstseed
    if ( present(seed1) ) lseed1 = seed1
    if ( present(seed2) ) lseed2 = seed2
    if ( present(seed3) ) lseed3 = seed3
  else 
    lseed1 = dseed1; lseed2 = dseed2; lseed3 = dseed3
  end if 
  
  if ( allocated( initialized_m ) ) then
    deallocate( lecuyer_stor_m, scramble_stor_m, fib_stor_m, mersenne_stor_m)
    deallocate( initialized_m, fib_ctr_m, bd_lastout_m, mersenne_ptr_m )
  end if
  allocate( lecuyer_stor_m(lecuyer_stor_sz, m) )
  allocate( scramble_stor_m(scramble_stor_sz, m) )
  allocate( fib_stor_m(fib_stor_sz, m), mersenne_stor_m(mersenne_stor_sz, m) )
  allocate( initialized_m(m), fib_ctr_m(m) )
  allocate( bd_lastout_m(m), mersenne_ptr_m(m) )

  init_single = present( no_init_single )
  if ( init_single ) init_single = .not. no_init_single
  if ( init_single ) then
    call associate_single()
    call init_rand_pl_main(lseed1, lseed2, lseed3)
  end if

  do j = 1, m
    jrpk = j
    call associate_multi(j)
    call init_rand_pl_main(lseed1+jrpk, lseed2+jrpk, lseed3+jrpk)
    initialized = 1
  end do

end subroutine allocate_rand_pl_multi


!=======================================================================
!
!  real(kind=wp) function rand_pl()
!
!  I: none
!
!  O: rand_pl     real(kind=wp)         random number in (0,1)
!
!  This is a wrapper for rand_pl_main that associates pointers for
!  the single generator before calling the main routine.
!
!=======================================================================

function rand_pl() 
  implicit none
  real(wp) :: rand_pl
  if ( curr_ptr .ne. 0 ) call associate_single()
  rand_pl = rand_pl_main()
end function rand_pl


!=======================================================================
!
!  real(kind=wp) function rand_pl_multi(m)
!
!  I: m              integer          selects which generator to use
!  
!  O: rand_pl_multi  real(kind=wp)    random number in (0,1)
!
!  This is a wrapper for rand_pl_main that associates pointers for
!  the multiple generators before calling the main routine. 
!
!=======================================================================

function rand_pl_multi(m) 
  implicit none
  integer, intent(in) :: m
  real(wp) :: rand_pl_multi
  if ( curr_ptr .ne. m ) call associate_multi(m)
  if ( initialized .ne. 1 ) call init_rand_pl_multi(m)
  rand_pl_multi = rand_pl_main()
end function rand_pl_multi


!=======================================================================
!
!  subroutine associate_single()
!
!  I: globals only
!
!  O: globals only
!
!  Associates pointers so that operations are on the single generator.
! 
!=======================================================================

subroutine associate_single()
  implicit none
  curr_ptr = 0
  lecuyer_stor  => lecuyer_stor_s
  scramble_stor => scramble_stor_s
  fib_stor      => fib_stor_s
  mersenne_stor => mersenne_stor_s
  initialized   => initialized_s
  fib_ctr       => fib_ctr_s
  bd_lastout    => bd_lastout_s
  mersenne_ptr  => mersenne_ptr_s
end subroutine associate_single


!=======================================================================
!
!  subroutine associate_multi(m)
!
!  I: m        integer        selects which of the set of multiple 
!                             generators to use
!
!  O: globals only
!
!  Associates pointers so that operations are on one of the multiple
!  generators, specified by the integer n.
!  
!=======================================================================
  
subroutine associate_multi(m)

  implicit none

  integer, intent(in) :: m

  if ( .not. allocated( initialized_m ) ) then
    write(0,*) 'Error (ASSOCIATE_MULTI): attempt to use multiple ', &
               'generators before initialization'
    stop
  end if

  if ( m .lt. 1 .or. m .gt. size(initialized_m) ) then
    write(0,*) 'Error (ASSOCIATE_MULTI): selector is out of allocated range'
    stop
  end if
  
  curr_ptr = m
  lecuyer_stor  => lecuyer_stor_m(:,m)
  scramble_stor => scramble_stor_m(:,m)
  fib_stor      => fib_stor_m(:,m)
  mersenne_stor => mersenne_stor_m(:,m)
  initialized   => initialized_m(m)
  fib_ctr       => fib_ctr_m(m)
  bd_lastout    => bd_lastout_m(m)
  mersenne_ptr  => mersenne_ptr_m(m)

end subroutine associate_multi


!=======================================================================
!
!  integer(kind=rpk) function mersenne()
!
!  I: only through module globals
!
!  O: mersenne  int(rpk)    random integer, in the range
!                           [0, 2**31 - 1] inclusive
!
!  Implements the Mersenne twister generator, returning one random
!  integer.
!
!=======================================================================

function mersenne()

  implicit none

  integer(rpk) :: mersenne

  integer(rpk), parameter :: n = 624, m = 397
  integer(rpk), parameter :: unity      =           1_rpk  ! 00000001
  integer(rpk), parameter :: matrix_a   = -1727483681_rpk  ! 9908b0df
  integer(rpk), parameter :: upp_p1     = -2147483647_rpk  ! 80000001
  integer(rpk), parameter :: upp_mask   =    upp_p1-unity  ! 80000000
  integer(rpk), parameter :: low_mask   =  2147483647_rpk  ! 7fffffff
  integer(rpk), parameter :: tmp_mask_b = -1658038656_rpk  ! 9d2c5680
  integer(rpk), parameter :: tmp_mask_c =  -272236544_rpk  ! efc60000
  integer(rpk), dimension(0:1), parameter :: mag01 = (/ 0_rpk, matrix_a /)

  integer, parameter :: shift_u = -11
  integer, parameter :: shift_s =   7
  integer, parameter :: shift_t =  15
  integer, parameter :: shift_l = -18

  integer(rpk) :: y
  integer :: k

  if ( mersenne_ptr .ge. n+1 ) then
    do k = 1, n-m
      y = ior( iand(mersenne_stor(k),   upp_mask), &
               iand(mersenne_stor(k+1), low_mask) )
      mersenne_stor(k) = ieor( ieor( mersenne_stor(k+m), ishft(y,-1) ), &
                               mag01( iand(y, unity) ) )
    end do

    do k = n-m+1, n-1
      y = ior( iand(mersenne_stor(k),   upp_mask), &
               iand(mersenne_stor(k+1), low_mask) )
      mersenne_stor(k) = ieor( ieor( mersenne_stor(k+m-n), ishft(y,-1) ), &
                               mag01( iand(y, unity) ) )
    end do

    y = ior( iand(mersenne_stor(n), upp_mask), &
             iand(mersenne_stor(1), low_mask) )
    mersenne_stor(n) = ieor( ieor( mersenne_stor(m), ishft(y,-1) ), &
                             mag01( iand(y, unity) ) )

    mersenne_ptr = 1
  end if

  y = mersenne_stor(mersenne_ptr)
  mersenne_ptr = mersenne_ptr + 1
  y = ieor( y, ishft(y, shift_u) )
  y = ieor( y, iand(ishft(y, shift_s), tmp_mask_b) )
  y = ieor( y, iand(ishft(y, shift_t), tmp_mask_c) )
  y = ieor( y, ishft(y, shift_l) )

  ! shift off 1 bit, so we get a 31-bit integer
  mersenne = ishft(y, -1)

end function mersenne


!=======================================================================
!
!  subroutine init_mersenne(seed)
!
!  I: seed     int(rpk)       initialization seed for generator,
!                              any nonzero integer is ok
!
!  O: none
!
!  Initializes Mersenne twister generator as in the sample code by
!  Matsumoto and Nishimura, using Knuth's LCG.
!
!=======================================================================
  
subroutine init_mersenne(seed)
  
  implicit none 
  
  integer(rpk), intent(in) :: seed

  integer(rpk), parameter :: mask32 = -1_rpk ! ffffffff

  ! check seed
  if ( iand(seed, mask32) .eq. 0_rpk ) then
    write(0,*) 'Error (INIT_MERSENNE): nonzero seed required'
    stop
  end if

  ! just in case the storage uses more than 32 bits
  mersenne_stor(1) = iand( seed, mask32 )
  
  do mersenne_ptr = 2, mersenne_stor_sz
    mersenne_stor(mersenne_ptr) = 1812433253_rpk * &
      ieor( mersenne_stor(mersenne_ptr-1), &
            ishft( mersenne_stor(mersenne_ptr-1), -30 ) ) + mersenne_ptr - 1
    mersenne_stor(mersenne_ptr) = iand( mersenne_stor(mersenne_ptr), mask32 )
  end do

  mersenne_ptr = mersenne_stor_sz + 1

end subroutine init_mersenne
  

!=======================================================================
!  
!  integer(kind=rpk) function lecuyer()
!
!  I: only through module globals
!
!  O: lecuyer   int(rpk)    random integer, in the range
!                           [0, 2**31 - 2] inclusive
!
!  Implements one iteration of the L'Ecuyer generator.
!
!=======================================================================

function lecuyer()

  implicit none

  integer(rpk) :: lecuyer

  integer(rpk), parameter :: m1  = 2147483647_rpk ! 2**31 - 1
  integer(rpk), parameter :: m2  = 2145483479_rpk

  integer(rpk), parameter :: a12 =      63308_rpk
  integer(rpk), parameter :: a13 =    -183326_rpk
  integer(rpk), parameter :: a21 =      86098_rpk
  integer(rpk), parameter :: a23 =    -539608_rpk

  integer(rpk), parameter :: q12 =      33921_rpk
  integer(rpk), parameter :: q13 =      11714_rpk
  integer(rpk), parameter :: q21 =      24919_rpk
  integer(rpk), parameter :: q23 =       3976_rpk

  integer(rpk), parameter :: r12 =      12979_rpk
  integer(rpk), parameter :: r13 =       2883_rpk
  integer(rpk), parameter :: r21 =       7417_rpk
  integer(rpk), parameter :: r23 =       2071_rpk

  integer(rpk) :: h, p12, p13, p21, p23
  integer(rpk), pointer :: x10, x11, x12, x20, x21, x22
  x10 => lecuyer_stor(1); x11 => lecuyer_stor(2); x12 => lecuyer_stor(3)
  x20 => lecuyer_stor(4); x21 => lecuyer_stor(5); x22 => lecuyer_stor(6)

  ! Component 1
  h = x10 / q13
  p13 = -a13 * (x10 - h * q13) - h * r13
  h = x11 / q12
  p12 =  a12 * (x11 - h * q12) - h * r12
  if ( p13 .lt. 0 ) p13 = p13 + m1
  if ( p12 .lt. 0 ) p12 = p12 + m1
  x10 = x11
  x11 = x12
  x12 = p12 - p13
  if ( x12 .lt. 0 ) x12 = x12 + m1

  ! Component 2
  h = x20 / q23
  p23 = -a23 * (x20 - h * q23) - h * r23
  h = x22 / q21
  p21 =  a21 * (x22 - h * q21) - h * r21
  if ( p23 .lt. 0 ) p23 = p23 + m2
  if ( p21 .lt. 0 ) p21 = p21 + m2
  x20 = x21
  x21 = x22
  x22 = p21 - p23
  if ( x22 .lt. 0 ) x22 = x22 + m2

  ! Combination
  if ( x12 .le. x22 ) then
    lecuyer = x12 - x22 + m1 - 1
  else
    lecuyer = x12 - x22 - 1
  end if

end function lecuyer


!=======================================================================
!
!  subroutine init_lecuyer(seed)
!
!  I: seed     int(rpk)       initialization seed for generator,
!                              different seeds give distinct sequences
!                              by partitioning the full sequence;
!                              any number in [0, 2**31 - 2] is ok
!
!  O: none
!
!  Initializes lecuyer generator by setting the seeds from a single
!  integer, by partitioning the full sequence and jumping ahead
!  by an amount proportional to the seed.
!
!=======================================================================

subroutine init_lecuyer(seed)

  implicit none
  
  integer(rpk), intent(in) :: seed

  integer(rpk) :: lseed

  integer(rpk), parameter :: m1  = 2147483647_rpk ! 2**31 - 1
  integer(rpk), parameter :: m2  = 2145483479_rpk

  integer(rpk), parameter :: a12 =      63308_rpk
  integer(rpk), parameter :: a13 =    -183326_rpk + m1
  integer(rpk), parameter :: a21 =      86098_rpk
  integer(rpk), parameter :: a23 =    -539608_rpk + m2

  integer(rpk), parameter :: logk = 184 - 32
  integer(rpk), dimension(3,3) :: A1, A2, A1k, A2k, A1ks, A2ks
  integer(rpk), dimension(3) :: tmparr

  integer :: j

  
  ! check seed
  if ( seed .lt. 0_rpk .or. seed .gt. 2147483646_rpk ) then
    write(0,*) 'Error (INIT_LECUYER): seed outside acceptable range'
    stop
  end if
  lseed = seed + 1

  ! single-iteration matrices
  A1(1,:) = (/   0_rpk,   1_rpk,   0_rpk /)
  A1(2,:) = (/   0_rpk,   0_rpk,   1_rpk /)
  A1(3,:) = (/   0_rpk,     a12,     a13 /)
  A2(1,:) = (/   0_rpk,   1_rpk,   0_rpk /)
  A2(2,:) = (/   0_rpk,   0_rpk,   1_rpk /)
  A2(3,:) = (/     a21,   0_rpk,     a23 /)

  ! compute k-iteration matrices, to jump between different partitions
  A1k = A1
  A2k = A2
  do j = 1, logk
    A1k = modmatsq(A1k, m1)  ! A1k = A1k * A1k (mod m1)
    A2k = modmatsq(A2k, m2)  ! A2k = A2k * A2k (mod m2)
  end do

  ! compute (A**k)**lseed matrices, to jump ahead to the location 
  !   marked by lseed
  A1ks = divconq(A1k, lseed, m1)
  A2ks = divconq(A2k, lseed, m2)

  ! set default seeds
  lecuyer_stor = (/ 1852689663, 1962642687, 580869375, &
                    2039711750, 1671394257, 879888250  /)

  ! advance by k*seed iterations (tmparr used to work around ifc bug)
  tmparr = lecuyer_stor(1:3)
  lecuyer_stor(1:3) = modmvmul(A1ks, tmparr, m1)
  tmparr = lecuyer_stor(4:6)
  lecuyer_stor(4:6) = modmvmul(A2ks, tmparr, m2)

end subroutine init_lecuyer


  !!!!
  !!!! routines for integer modular matrix/vector/scalar arithmetic
  !!!!

  ! divide-and-conquer to calculate A**j mod m
  recursive function divconq(A, j, m) result (Aj)
    implicit none
    integer(rpk), dimension(:,:), intent(in) :: A
    integer(rpk), intent(in) :: j, m
    integer(rpk), dimension(size(A,1),size(A,2)) :: Aj
    if ( j .eq. 1 ) then
      Aj = A
    else if ( isodd(j) ) then
      Aj = modmatmul(A, divconq(A, j-1_rpk, m), m)
    else
      Aj = modmatsq(divconq(A, j/2_rpk, m), m)
    end if
  end function divconq
    
  ! check if integer scalar is odd
  function isodd(a)
    implicit none
    integer(rpk), intent(in) :: a
    logical :: isodd
    isodd = (iand(a, 1_rpk) .eq. 1_rpk)
  end function isodd

  ! wrapper for modmatmul to square a matrix
  function modmatsq(A, m)
    implicit none
    integer(rpk), dimension(:,:), intent(in) :: A
    integer(rpk), intent(in) :: m
    integer(rpk), dimension(size(A,1),size(A,2)) :: modmatsq
    modmatsq = modmatmul(A, A, m)
  end function modmatsq

  ! calculates the matrix product A * B mod m, assuming nonnegative
  ! entries; A and B are assumed square
  function modmatmul(A, B, m)
    implicit none
    integer(rpk), dimension(:,:), intent(in) :: A, B
    integer(rpk), intent(in) :: m
    integer(rpk), dimension(size(A,1),size(A,2)) :: modmatmul
    integer(rpk) :: n, j
    n = size(A,1)
    if ( n.ne.size(A,2) .or. n.ne.size(B,1) .or. n.ne.size(B,2) ) then
      write(0,*) 'Error (MODMATMUL): nonsquare or incompatible matrix sizes'
      stop
    end if
    do j = 1, n
      modmatmul(:,j) = modmvmul(A, B(:,j), m)
    end do
  end function modmatmul

  ! calculates the matrix-vector product A * x mod m, assuming 
  ! nonnegative entries; A is assumed square
  function modmvmul(A, x, m)
    implicit none
    integer(rpk), dimension(:,:), intent(in) :: A
    integer(rpk), dimension(:),   intent(in) :: x
    integer(rpk), intent(in) :: m
    integer(rpk), dimension(size(A,1)) :: modmvmul
    integer(rpk) :: n, j
    n = size(A,1)
    if ( n.ne.size(A,2) .or. n.ne.size(x) ) then
      write(0,*) 'Error (MODMVMUL): nonsquare matrix or incompatible vector'
      stop
    end if
    do j = 1, n
      modmvmul(j) = moddotprod(A(j,:), x, m)
    end do
  end function modmvmul

  ! calculates the vector dot product x * y mod m, assuming 
  ! nonnegative entries
  function moddotprod(x, y, m)
    implicit none
    integer(rpk), dimension(:),   intent(in) :: x, y
    integer(rpk), intent(in) :: m
    integer(rpk) :: moddotprod
    integer(rpk) :: n, j
    n = size(x)
    if ( n.ne.size(y) ) then
      write(0,*) 'Error (MODDOTPROD): incompatible vector sizes'
      stop
    end if
    moddotprod = modmult(x(1), y(1), m)
    do j = 2, n
      moddotprod = modadd( moddotprod, modmult(x(j), y(j), m), m )
    end do
  end function moddotprod

  ! calculates x+y mod m without overflow, assuming x, y >= 0
  function modadd(x, y, m)
    implicit none
    integer(rpk), intent(in) :: x, y, m
    integer(rpk) :: modadd
    modadd = (x - m) + y; if ( modadd .lt. 0 ) modadd = modadd + m
  end function modadd
    
  ! calculates x*y mod m without overflow via the decomposition
  ! algorithm; see L'Ecuyer and Cote, ACM Trans. Math. Soft. 17, 98 
  ! (1991). Assumes x, y >= 0.  Note the redundant p = p+m statements
  ! that were necessary to work around an intel compiler bug.
  function modmult(x, y, m)
    implicit none
    integer(rpk), intent(in) :: x, y, m
    integer(rpk) :: modmult
    integer(rpk), parameter :: h = 32768 ! for 32-bit integers
    integer(rpk) :: a, s, a0, a1, q, qh, rh, k, p
    if ( x .lt. 0 .or. y .lt. 0 ) then
      write(0,*) 'Error (MODMULT): inputs must be positive'
      stop
    end if
    if ( x .eq. 0 .or. y .eq. 0 ) then
      modmult = 0; return
    end if
    a = x; s = y
    if ( a .lt. h ) then
      a0 = a; p = 0
    else
      a1 = a / h; a0 = a - h*a1
      qh = m / h; rh = m - h*qh
      if ( a1 .ge. h ) then
        a1 = a1 - h; k = s / qh
        p = h * (s - k*qh) - k * rh
        if ( p .lt. 0 ) p = p + m
        do while ( p .lt. 0 ); p = p + m; end do
      else
        p = 0
      end if
      if ( a1 .ne. 0 ) then
        q = m / a1; k = s / q
        p = p - k * (m - a1*q)
        if ( p .gt. 0 ) p = p - m
        p = p + a1 * (s - k*q)
        if ( p .lt. 0 ) p = p + m
        do while ( p .lt. 0 ); p = p + m; end do
      end if
      k = p / qh; p = h * (p - k*qh) - k * rh
      if ( p .lt. 0 ) p = p + m
      do while ( p .lt. 0 ); p = p + m; end do
    end if
    if ( a0 .ne. 0 ) then
      q = m / a0; k = s / q
      p = p - k * (m - a0*q)
      if ( p .gt. 0 ) p = p - m
      p = p + a0 * (s - k*q)
      if ( p .lt. 0 ) p = p + m
      do while ( p .lt. 0 ); p = p + m; end do
    end if
    modmult = p
  end function modmult


!=======================================================================
!  
!  integer(kind=rpk) function fibonacci(optional force_init)
!
!  I: force_init  logical    (optional) force a new set of 1009 
!                            numbers to be generated
!
!  O: fibonacci   int(rpk)    random integer, in the range
!                           [0, 2**31 - 2] (note that the lsb is 
!                                           always 0)
!
!  Implements one iteration of Knuth's subtractive lagged Fibonacci
!    generator, with skipping to improve randomness.
!
!=======================================================================

function fibonacci(force_init)

  implicit none

  logical, optional :: force_init

  integer(rpk) :: fibonacci

  integer(rpk), parameter :: lag = 37, sbst = 100, total = 1009
  integer(rpk), parameter :: mdls = 1073741824_rpk ! 2**30
  logical :: reinit
  integer(rpk) :: j

  ! check for forcing reinitialization
  reinit = .false.
  if ( present( force_init ) ) reinit = force_init

  ! generate 1009 numbers, if needed
  if ( fib_ctr .gt. sbst .or. reinit ) then
    aa(1:sbst) = fib_stor(1:sbst)
    do j = sbst+1, total
      aa(j) = iand(aa(j-sbst) - aa(j-lag), mdls-1_rpk)
    end do
    do j = 1, lag
      fib_stor(j) = iand(aa(total+j-sbst) - aa(total+j-lag), mdls-1_rpk)
    end do
    do j = lag+1, sbst
      fib_stor(j) = iand(aa(total+j-sbst) - fib_stor(j-lag), mdls-1_rpk)
    end do
    fib_ctr = 1
  end if

  fibonacci = 2 * fib_stor(fib_ctr)
  fib_ctr = fib_ctr + 1

end function fibonacci
  

!=======================================================================
!  
!  subroutine init_fibonacci(seed)
!
!  I: seed    int(rpk)   initialization seed, in the range
!                        [0, 2**30 - 3 = 1 073 741 821]
!
!  O: only through module globals
!
!  Initializes Knuth's Fibonacci generator, according to his
!    recommendations.
!
!=======================================================================

subroutine init_fibonacci(seed)

  implicit none

  integer(rpk), intent(in) :: seed

  integer(rpk), parameter :: tt = 70, lag = 37 
  integer(rpk), parameter :: mdls = 1073741824_rpk ! 2**30
  integer(rpk), parameter :: maxseed = mdls-3_rpk
  integer, parameter :: kk = fib_stor_sz, kkk = kk+kk-1
  integer(rpk), dimension(kkk) :: xtmp
  integer :: t, j, ss
  integer(rpk) :: dummy


  ! check seed
  if ( seed .lt. 1 .or. seed .gt. maxseed ) then
    write(0,*) 'Error (INIT_FIBONACCI): inappropriate seed'
  end if

  ! set pointer
  fib_ctr = 1

  ! do knuthian stuff
  ss = seed - mod(seed,2_rpk) + 2
  do j = 1, kk
    xtmp(j) = ss
    ss = ss + ss
    if ( ss .ge. mdls ) ss = ss - mdls + 2
  end do
  xtmp(2) = xtmp(2) + 1
  ss = seed
  t = tt - 1
10 continue
  do j = kk, 2, -1
    xtmp(j+j-1) = xtmp(j)
    xtmp(j+j-2) = 0
  end do
  do j = kkk, kk+1, -1
    xtmp(j-(kk-lag)) = xtmp(j-(kk-lag)) - xtmp(j)
    if ( xtmp(j-(kk-lag) ) .lt. 0)  xtmp(j-(kk-lag)) = xtmp(j-(kk-lag)) + mdls
    xtmp(j-kk) = xtmp(j-kk) - xtmp(j)
    if ( xtmp(j-kk) .lt. 0 )  xtmp(j-kk) = xtmp(j-kk) + mdls
  end do
  if ( mod(ss, 2) .eq. 1 ) then
    do j = kk, 1, -1
      xtmp(j+1) = xtmp(j)
    end do
    xtmp(1) = xtmp(kk+1)
    xtmp(lag+1) = xtmp(lag+1) - xtmp(kk+1)
    if ( xtmp(lag+1) .lt. 0 )  xtmp(lag+1) = xtmp(lag+1) + mdls
  end if
  if ( ss .ne. 0 ) then
    ss = ss / 2
  else
    t = t - 1
  end if
  if ( t .gt. 0 ) go to 10
  do j = 1, lag
    fib_stor(j+kk-lag) = xtmp(j)
  end do
  do j = lag+1, kk
    fib_stor(j-lag) = xtmp(j)
  end do

  ! discard the first 2020 or so numbers
  do j = 1, 2
    dummy = fibonacci(force_init=.true.)
  end do

end subroutine init_fibonacci


!=======================================================================
!
!  subroutine init_rand_pl_main(seed1, seed2, seed3)
!
!  I: seed1,2,3  int(rpk)     seeds for initializing generator
!
!  O: through module globals
!
!  Initializes all 3 random number generators, with the
!  specified seeds (1 for Fibonacci, 2 for L'Ecuyer, 3 for Mersenne).  
!
!=======================================================================

subroutine init_rand_pl_main(seed1, seed2, seed3)

  implicit none

  integer(rpk), intent(in) :: seed1, seed2, seed3
  integer :: j


  ! check method flag
  if ( rand_pl_mf_pri .lt. 1 .or. rand_pl_mf_pri .gt. 3 .or. &
       rand_pl_mf_mix .lt. 0 .or. rand_pl_mf_mix .gt. 3 .or. &
       rand_pl_mf_scr .lt. 0 .or. rand_pl_mf_scr .gt. 4 .or. &
       rand_pl_mf_pri .eq. rand_pl_mf_mix .or. &
       rand_pl_mf_pri .eq. rand_pl_mf_scr .or. &
       (rand_pl_mf_mix .ne. 0 .and. rand_pl_mf_mix .eq. rand_pl_mf_scr) ) then
    write(0,*) 'Error (INIT_RAND_PL): illegal rand_pl_mf'
    stop
  end if

  ! initialize all generators
  call init_fibonacci(seed1)
  call init_lecuyer(seed2)
  call init_mersenne(seed3)

  ! fill scramble table
  bd_lastout = 0
  scramble_stor = 0
  if ( rand_pl_mf_scr .ne. 0 ) then
    do j = 1, scramble_stor_sz
      scramble_stor(j) = rand_pl_raw()
    end do
    if ( rand_pl_mf_scr .eq. 4 ) bd_lastout = rand_pl_raw()
  end if

  initialized = 1

end subroutine init_rand_pl_main


!=======================================================================
!
!  real(kind=wp) function rand_pl_raw()
!
!  I: none
!
!  O: rand_pl     integer(rpk)         random integer 
!                              
!  Constructs a random integer from the proper generator(s) and
!  subtractively combines them, if appropriate.  Basically this
!  does everything short of shuffling the output numbers.
!
!=======================================================================

function rand_pl_raw() result(r)

  implicit none

  integer(rpk) :: r

  integer(rpk), parameter :: modmask = 2147483647_rpk  ! 7fffffff = 2**31-1


  if ( rand_pl_mf_pri .eq. 1 ) then
    r = fibonacci()
  else if ( rand_pl_mf_pri .eq. 2 ) then
    r = lecuyer()
  else if ( rand_pl_mf_pri .eq. 3 ) then
    r = mersenne()
  end if

  if ( rand_pl_mf_mix .ne. 0 ) then
    if ( rand_pl_mf_mix .eq. 1 ) then
      r = iand( r - fibonacci(), modmask )
    else if ( rand_pl_mf_mix .eq. 2 ) then
      r = iand( r - lecuyer(), modmask )
    else if ( rand_pl_mf_mix .eq. 3 ) then
      r = iand( r - mersenne(), modmask )
    end if
  end if

end function rand_pl_raw


!=======================================================================
!
!  real(kind=wp) function rand_pl_main()
!
!  I: none
!
!  O: rand_pl     real(kind=wp)         random number in (0,1)
!
!=======================================================================

function rand_pl_main()

  implicit none

  real(wp) :: rand_pl_main

  integer(rpk), parameter :: fibonacci_range = 1073741823_rpk ! 2**31 - 2
  real(wp), parameter :: fibonacci_mult = 4.65661287307739e-10_wp ! 1/(2**31)
  integer(rpk), parameter :: lecuyer_range   = 2147483646_rpk ! 2**31 - 2
  real(wp), parameter :: lecuyer_mult = 4.65661287307739e-10_wp ! 1/(2**31)
  integer(rpk), parameter :: mersenne_range  = 2147483647_rpk ! 2**31 - 1
  real(wp), parameter :: mersenne_mult = 4.65661287090899e-10_wp ! 1/(2**31 + 1)

  real(wp), parameter :: gen_mult = 4.65661287090899e-10_wp ! 1/(2**31 + 1)

  integer :: indx

 
  ! this can only be true if this is the single generator
  if ( initialized .ne. 1 ) call init_rand_pl()

  if ( rand_pl_mf_scr .eq. 0 ) then
    rand_pl_main = (real(rand_pl_raw(), kind=wp) + 1) * gen_mult

  else if ( rand_pl_mf_scr .eq. 1 ) then
    ! scramble output using Fibonacci output
    indx = ceiling(fibonacci() * gen_mult * scramble_stor_sz)
    rand_pl_main = (real(scramble_stor(indx), kind=wp) + 1) * gen_mult
    scramble_stor(indx) = rand_pl_raw()

  else if ( rand_pl_mf_scr .eq. 2 ) then
    ! scramble output using L'Ecuyer output
    indx = ceiling(lecuyer() * gen_mult * scramble_stor_sz)
    rand_pl_main = (real(scramble_stor(indx), kind=wp) + 1) * gen_mult
    scramble_stor(indx) = rand_pl_raw()

  else if ( rand_pl_mf_scr .eq. 3 ) then
    ! scramble output using Mersenne output
    indx = ceiling(mersenne() * gen_mult * scramble_stor_sz)
    rand_pl_main = (real(scramble_stor(indx), kind=wp) + 1) * gen_mult
    scramble_stor(indx) = rand_pl_raw()

  else if ( rand_pl_mf_scr .eq. 4 ) then
    ! scramble output using Bays-Durham scrambling
    indx = ceiling(bd_lastout * gen_mult * scramble_stor_sz)
    bd_lastout = scramble_stor(indx)
    rand_pl_main = (real(bd_lastout, kind=wp) + 1) * gen_mult
    scramble_stor(indx) = rand_pl_raw()

  end if

end function rand_pl_main


!=======================================================================
!
!  real(kind=wp) function rand_pl_vec(n)
!
!  I: n            default integer      length of output vector
!
!  O: rand_pl_vec  real(kind=wp) array  array of random numbers in (0,1)
!
!  Single wrapper for rand_pl_vec_main()
!
!=======================================================================

function rand_pl_vec(n)
  implicit none
  integer, intent(in) :: n
  real(wp), dimension(n) :: rand_pl_vec
  if ( curr_ptr .ne. 0 ) call associate_single()
  rand_pl_vec = rand_pl_vec_main(n)
end function rand_pl_vec


!=======================================================================
!
!  real(kind=wp) function rand_pl_vec_multi(m, n)
!
!  I: m            default integer      selects which generator to use
!     n            default integer      length of output vector
!
!  O: rand_pl_vec  real(kind=wp) array  array of random numbers in (0,1)
!
!  Single wrapper for rand_pl_vec_main()
!
!=======================================================================

function rand_pl_vec_multi(m, n)
  implicit none
  integer, intent(in) :: m, n
  real(wp), dimension(n) :: rand_pl_vec_multi
  if ( curr_ptr .ne. m ) call associate_multi(m)
  rand_pl_vec_multi = rand_pl_vec_main(n)
end function rand_pl_vec_multi


!=======================================================================
!
!  real(kind=wp) function rand_pl_vec_main(n)
!
!  I: n            default integer      length of output vector
!
!  O: rand_pl_vec_main  real(kind=wp) array  array of random numbers 
!                                            in (0,1)
!
!=======================================================================

function rand_pl_vec_main(n)
  implicit none
  integer, intent(in) :: n
  real(wp), dimension(n) :: rand_pl_vec_main
  integer :: j
  do j = 1, n
    rand_pl_vec_main(j) = rand_pl_main()
  end do
end function rand_pl_vec_main


!=======================================================================
!
!  real(kind=wp) function nrand_pl_vec(n)
!
!  I: n            default integer      length of output vector, must
!                                         be even
!
!  O: nrand_pl_vec real(kind=wp) array  array of standard normal
!                                         deviates
!
!  Single wrapper for nrand_pl_vec_main
!
!=======================================================================

function nrand_pl_vec(n)
  implicit none
  integer, intent(in) :: n
  real(wp), dimension(n) :: nrand_pl_vec
  if ( curr_ptr .ne. 0 ) call associate_single()
  nrand_pl_vec = nrand_pl_vec_main(n)
end function nrand_pl_vec


!=======================================================================
!
!  real(kind=wp) function nrand_pl_vec_multi(m, n)
!
!  I: m            default integer      selects which generator to use
!     n            default integer      length of output vector, must
!                                         be even
!
!  O: nrand_pl_vec real(kind=wp) array  array of standard normal
!                                         deviates
!  
!  Multiple wrapper for nrand_pl_vec_main
!  
!=======================================================================

function nrand_pl_vec_multi(m, n)
  implicit none
  integer, intent(in) :: m, n
  real(wp), dimension(n) :: nrand_pl_vec_multi
  if ( curr_ptr .ne. m ) call associate_multi(m)
  nrand_pl_vec_multi = nrand_pl_vec_main(n)
end function nrand_pl_vec_multi


!=======================================================================
!
!  real(kind=wp) function nrand_pl_vec_main(n)
!
!  I: n            default integer      length of output vector, must
!                                         be even
!
!  O: nrand_pl_vec_main real(kind=wp) array  array of standard normal
!                                            deviates
!
!=======================================================================

function nrand_pl_vec_main(n) result(nd)

  implicit none

  integer, intent(in) :: n
  real(wp), dimension(n) :: nd

  real(wp), parameter :: pi = 3.141592653589793_wp, twopi = 2*pi
  integer :: j
  real(wp) :: tmp

  
  if ( (n/2)*2 .ne. n ) then
    write(0,*) 'Error (NRAND_PL_VEC): requested length not even'
    stop
  end if

  nd = rand_pl_vec_main(n)

  do j = 1, n, 2
    tmp = sqrt((-2.0_wp) * log(nd(j)))
    nd(j)   = cos(twopi * nd(j+1)) * tmp
    nd(j+1) = sin(twopi * nd(j+1)) * tmp
  end do

end function nrand_pl_vec_main


!=======================================================================
!
!  subroutine get_rand_pl_state(ivec)
!
!  I: only through globals
!
!  O: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator
!
!  Single wrapper for get_rand_pl_state_main
!
!=======================================================================

subroutine get_rand_pl_state(ivec)
  implicit none
  integer(rpk), dimension(rand_pl_state_sz), intent(out) :: ivec
  call associate_single()
  call get_rand_pl_state_main(ivec)
end subroutine get_rand_pl_state


!=======================================================================
!
!  subroutine get_rand_pl_state_multi(m, ivec)
!
!  I: m      integer            selects which generator to use
!
!  O: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator
!
!  Single wrapper for get_rand_pl_state_main
!
!=======================================================================
    
subroutine get_rand_pl_state_multi(m, ivec)
  implicit none
  integer(rpk), dimension(rand_pl_state_sz), intent(out) :: ivec
  integer, intent(in) :: m
  if ( curr_ptr .ne. m ) call associate_multi(m)
  call get_rand_pl_state_main(ivec)
end subroutine get_rand_pl_state_multi


!=======================================================================
!
!  subroutine get_rand_pl_state_main(ivec)
!
!  I: only through globals
!
!  O: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator
!
!=======================================================================

subroutine get_rand_pl_state_main(ivec)

  implicit none

  integer(rpk), dimension(rand_pl_state_sz), intent(out) :: ivec
  integer :: ptr1, ptr2

  if ( initialized .eq. 0 ) then
    ivec = 0
    ivec(1) = rand_pl_mf
  else
    ivec(1) = rand_pl_mf
    ivec(2) = initialized
    ivec(3) = fib_ctr
    ivec(4) = bd_lastout
    ivec(5) = mersenne_ptr
    ptr1 = 6
    ptr2 = ptr1 + fib_stor_sz - 1
    ivec(ptr1:ptr2) = fib_stor
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + lecuyer_stor_sz - 1
    ivec(ptr1:ptr2) = lecuyer_stor
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + scramble_stor_sz - 1
    ivec(ptr1:ptr2) = scramble_stor
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + mersenne_stor_sz - 1
    ivec(ptr1:ptr2) = mersenne_stor
  end if

end subroutine get_rand_pl_state_main


!=======================================================================
!
!  subroutine set_rand_pl_state(ivec)
!
!  I: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator,
!                               from previous call to get_rand_pl_state
!
!  O: only through globals
!
!  Single wrapper for set_rand_pl_state_main.
!
!=======================================================================

subroutine set_rand_pl_state(ivec)
  implicit none
  integer(rpk), dimension(rand_pl_state_sz), intent(in) :: ivec
  call associate_single()
  call set_rand_pl_state_main(ivec)
end subroutine set_rand_pl_state


!=======================================================================
!
!  subroutine set_rand_pl_state_multi(m, ivec)
!  
!  I: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator,
!                               from previous call to get_rand_pl_state
!     m      integer            selects which generator to use
!                               !  O: only through globals
!  
!  Multiple generator wrapper for set_rand_pl_state_main.
!  
!=======================================================================

subroutine set_rand_pl_state_multi(m, ivec)
  implicit none
  integer(rpk), dimension(rand_pl_state_sz), intent(in) :: ivec
  integer, intent(in) :: m
  if ( curr_ptr .ne. m ) call associate_multi(m)
  call set_rand_pl_state_main(ivec)
end subroutine set_rand_pl_state_multi


!=======================================================================
!
!  subroutine set_rand_pl_state_main(ivec)
!
!  I: ivec   int(rpk) array     contains sufficient information to
!                               restore the state of the generator,
!                               from previous call to get_rand_pl_state
!
!  O: only through globals
!
!=======================================================================

subroutine set_rand_pl_state_main(ivec)

  implicit none

  integer(rpk), dimension(rand_pl_state_sz), intent(in) :: ivec
  integer :: ptr1, ptr2

  if ( ivec(1) .ne. rand_pl_mf ) then
    write(0,*) 'Error (SET_RAND_PL_STATE): rand_pl_mf changed since ' &
               // 'state vector was taken'
    stop
  else
    initialized   = ivec(2)
    fib_ctr       = ivec(3)
    bd_lastout    = ivec(4)
    mersenne_ptr  = ivec(5)
    ptr1 = 6
    ptr2 = ptr1 + fib_stor_sz - 1
    fib_stor      = ivec(ptr1:ptr2)
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + lecuyer_stor_sz - 1
    lecuyer_stor  = ivec(ptr1:ptr2)
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + scramble_stor_sz - 1
    scramble_stor = ivec(ptr1:ptr2)
    ptr1 = ptr2 + 1
    ptr2 = ptr1 + mersenne_stor_sz - 1
    mersenne_stor = ivec(ptr1:ptr2)
  end if

end subroutine set_rand_pl_state_main

end module random_pl
