!-----------------------------------------------------------------------
!
!  module linsol
!
!  Module of utilities for finding solutions to linear systems of
!  equations described by square matrices.
!
!  Precision kind is obtained through the 'globals' module, which
!  must define the integer parameter 'wp'.
!
!  Routines are based on the ones in the Bowdler et al. article in 
!  Wilkinson and Reinsch, "Linear Algebra" 
!  (Springer-Verlag, Berlin, 1971), p. 93;
!  also based on some of the routines in the LINPACK library for
!  improved performance.
!
!  Usage: 
!
!  Call any of several routines to perform desired computations
!  based on an LU decomposition.  Here, the optional output
!  'success' in any routine is .false. if problems were encountered;
!  otherwise a warning message goes to standard error.  The
!  input matrix A is never modified.
!
!    determinant(A, optional success) 
!       -> function to return determinant of A.
!
!    linsolve(A, b, optional resid, optional success,
!                                   optional tol, optional niter)
!       -> function to take input column vector and replace with 
!          solution to Ax=b by backsubstitution of the LU decomposition 
!          and iteration to improve accuracy.  b can also be a 2D arrray
!          of column vectors, all of which will be solved.  The maximum
!          amplitude of the residuals is available through 'resid'.
!          Note that for a square matrix B, linsolve(A, B) returns the 
!          matrix product A^(-1) B (more efficiently and accurately 
!          than computing the inverse of A and then multiplying by B).
!          This routine uses iteration to improve the accuracy of the
!          solutions.  Iteration is controlled by the inputs 'tol'
!          and 'niter', which have these effects:
!            neither present -> no iteration, just simple backsub
!            tol present -> perform automatic iteration until relative
!                           max residual is less than tol; set tol=0 
!                           for auto tolerance (max resid/max soln = 
!                           2*N*eps).
!            niter present -> iterate this fixed number of times to 
!                           attempt accuracy improvement.  If niter 
!                           is 0, this is equivalent to simple 
!                           backsubstitution.
!            both present -> not allowed
!
!    matrix_inverse(A, b, optional resid, optional success)
!       -> function interface to linsolve, uses iterated back-
!          substitution to find the matrix inverse.
!
!  Description of input/output arguments:
!
!    A -> (input only) input matrix to be LU factorized; should
!         be square and either real or complex (of type wp).
!
!    b -> (input only) input column vector or 2D array of column
!         vectors to be solved; should be of same type as A.
!
!    resid -> (optional, output only) max residual of solution vector;
!         should be a real(wp) scalar.  Note that evalutating this
!         requires extra computation (except when iterating).
!
!    tol -> (optional, input only) iteration tolerance, as 
!         described above; should be a real(wp) scalar.
!
!    niter -> (optional, input only) fixed number of iterations
!         to perform; should be an integer scalar.
!
!    successs -> (optional, output only) logical output, indicating
!         whether or not problems were encountered that would
!         affect backsubstituted solutions.
!
!  Note that all of these calls make automatic calls to an LU
!  decomposition routine.  Each new call makes a fresh call to
!  the LU routine; there is currently no facility for reusing
!  the factorization other than using linsolve to solve many
!  vectors at once.
!
!-----------------------------------------------------------------------


module linsol

  use globals, only : wp

  private
  public determinant, linsolve, matrix_inverse

  ! generic procedure names
  interface determinant
    module procedure rdeterminant
    module procedure zdeterminant
  end interface determinant

  interface linsolve
    module procedure rlinsolve, rlinsolve_multi
    module procedure zlinsolve, zlinsolve_multi
  end interface linsolve

  interface matrix_inverse
    module procedure rmatrix_inverse
    module procedure zmatrix_inverse
  end interface matrix_inverse

  ! module global storage
  real(wp),    dimension(:,:), allocatable :: rLU
  complex(wp), dimension(:,:), allocatable :: zLU
  integer, dimension(:), allocatable :: indx

contains

subroutine initialize_module_linsol()
  character(127), parameter :: RCS_ID = &
    "@(#)$Id:$"
  write(0,*) trim(RCS_ID)
end subroutine initialize_module_linsol


!-----------------------------------------------------------------------
!
!  subroutine rLU_decomp(A, optional success)
!
!  Routine for LU decomposition of a real matrix.  Prerequisite for
!  calls to other real-matrix routines in this module.
!
!  The variable success is an optional output; when LU decomposition
!  fails, the presence of this variable ensures that execution will
!  continue (with the output value set to false).
!
!  Based on the dgefa routine in the LINPACK library.
! 
!-----------------------------------------------------------------------

subroutine rLU_decomp(A, success)

  implicit none

  real(wp), dimension(:,:), intent(in) :: A
  logical, optional, intent(out) :: success

  integer :: n, j, k, l
  real(wp) :: x
  logical :: locsuccess

  n = size(A, dim=1)
  if ( n .ne. size(A, dim=2 ) ) then
    write(0,*) 'Error (LU_DECOMP): input matrix not square'
    stop
  end if

  allocate( rLU(n, n), indx(n) )
  rLU = A
  if ( n .eq. 1 )  return

  locsuccess = .true.; if ( present(success) )  success = .true.
  do k = 1, n-1
    l = maxloc(abs(rLU(k:n,k)), dim=1) + k - 1 ! pivot row
    indx(k) = l
    if ( rLU(l,k) .ne. 0 ) then 
      !if ( l .ne. k ) then
        x = rLU(l,k); rLU(l,k) = rLU(k,k); rLU(k,k) = x
      !end if
      rLU((k+1):n,k) = rLU((k+1):n,k) * (-1/rLU(k,k))
      do j = k+1, n
        x = rLU(l,j)
        !if ( l .ne. k ) then
          rLU(l,j) = rLU(k,j)
          rLU(k,j) = x
        !end if
        if (x.ne.0) rLU((k+1):n,j) = rLU((k+1):n,j) + x * rLU((k+1):n,k)
      end do !j
    else
      locsuccess = .false.
    end if
  end do !k
  indx(n) = n

  if ( .not. locsuccess .or. rLU(n,n) .eq. 0 ) then
    if ( present(success) ) then
      success = .false.
    else
      write(0,*) 'Warning (LU_DECOMP): input matrix singular'
    end if
  end if

end subroutine rLU_decomp


!-----------------------------------------------------------------------
!
!  real(wp) function rdeterminant(A, optional success)
!
!  Returns determinant of matrix A by calling LU_decomp. 
!  Returns as a 2-element array, with the result being
!  out(1) * 10**out(2).
!
!  Based on LINPACK routine dgedi.f
!
!-----------------------------------------------------------------------

function rdeterminant(A, success)

  implicit none

  real(wp), dimension(:,:), intent(in) :: A
  logical, intent(out), optional :: success
  real(wp), dimension(2) :: rdeterminant
  real(wp) :: rd1, rd2
  integer :: j, n

  n = size(A, dim=1)
  if ( present(success) ) then
    call rLU_decomp(A, success=success)
  else
    call rLU_decomp(A)
  end if
  rd1 = 1; rd2 = 0
  do j = 1, n
    if ( indx(j) .ne. j )  rd1 = -rd1
    rd1 = rLU(j,j)*rd1
    if ( rd1 .eq. 0 )  exit
    do while ( abs(rd1) .lt. 1 ); rd1 = 10*rd1; rd2 = rd2-1; end do
    do while ( abs(rd1) .gt. 10 ); rd1 = 0.1_wp*rd1; rd2 = rd2+1; end do
  end do
  rdeterminant = (/ rd1, rd2 /) 
  deallocate( rLU, indx )

end function rdeterminant
 
        
!-----------------------------------------------------------------------
!
!  function rlinsolve(A, b, optional resid, optional success, 
!                     optional tol, optional niter) 
!
!  Produces solutions to Ax=b for column vector b, by making an
!  LU decomposition and backsubstituting with iteration.
!  Here, b is not modified.
!
!  The infinity norm (i.e., max abs value) of the residuals is 
!  available in the optional output resid.
!
!  The variable success is an optional output; when the system is
!  ill-conditioned, the presence of this variable ensures that 
!  execution will continue (with the output value set to false).
!
!  Note that the precision of this routine is limited by the
!  'bb = b - matmul(A, x)' statement, where the residuals are 
!  calculated.  This can be improved by performing this in
!  a higher precision mode than everything else, to get
!  residuals that are nearly machine epsilon.  On some architectures
!  there is a more accurate accumulator (e.g., Pentium), in which
!  case this is automatic.
!
!  To handle this, specify optionals 'tol' (real) and 'niter' (int) 
!  as follows:
!
!    neither present -> no iteration, just simple backsub
!    tol present -> perform automatic iteration until relative
!                   max residual is less than tol; set tol=0 for 
!                   auto tolerance (max resid/max soln = 2*N*eps).
!    niter present -> iterate this fixed number of times to attempt
!                   accuracy improvement.  If niter is 0, this
!                   is equivalent to simple backsubstitution.
!    both present -> not allowed
!
!  The output success indicates success of iteration and LU
!  decomposition together, as appropriate.
!
!-----------------------------------------------------------------------

function rlinsolve(A, b, resid, success, tol, niter) result (x)

  implicit none

  real(wp), dimension(:,:), intent(in   ) :: A
  real(wp), dimension(:),   intent(in   ) :: b
  real(wp),                 intent(  out), optional :: resid
  logical,                  intent(  out), optional :: success
  real(wp),                 intent(in   ), optional :: tol
  integer,                  intent(in   ), optional :: niter

  real(wp), dimension(size(b)) :: x
  logical :: rptflag, LUsuccess, locsuccess, usetol
  integer :: n, l, maxiter
  real(wp) :: xmax, bbmax, dd0, dd1, tolfac
  real(wp), dimension(size(b)) :: bb

  n = size(A, dim=1)
  if ( n .ne. size(b, dim=1) ) then
    write(0,*) 'Error (LINSOLVE): input vector of wrong length'
    stop
  end if

  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (LINSOLVE): tol and niter cannot both be present'
    stop
  end if
  usetol = .false.
  tolfac = 0
  maxiter = 100
  if ( present(tol) ) then
    usetol = .true.
    if ( tol .gt. 0 ) then
      tolfac = tol
    else if ( tol .eq. 0 ) then
      tolfac = 2*n*epsilon(1.0_wp)
    else
      write(0,*) 'Error (LINSOLVE): illegal value for tol'
      stop
    end if
  else if ( present(niter) ) then
    if ( niter .lt. 0 ) then
      write(0,*) 'Error (LINSOLVE): illegal value for niter'
      stop
    end if
    maxiter = niter
  else
    maxiter = 0
  end if
  
  locsuccess = .true.
  call rLU_decomp(A, success=LUsuccess)
  if ( .not. present(success) .and. .not. LUsuccess ) then
    write(0,*) 'Warning (LINSOLVE): problems with call to LUdecomp'
  end if

  l = 0
  dd0 = 0
  rptflag = .true.
  bb = b
  x = 0

  do while ( rptflag )
    call rlinsolve_backsub(bb)
    l = l + 1
    rptflag = .false.
    x = x + bb
    xmax = maxval(abs(x)); bbmax = maxval(abs(bb))
    bb = b - matmul(A, x)
    dd1 = bbmax / xmax
    if ( bbmax .gt. tolfac*xmax )  rptflag = .true.

    if ( dd1 .gt. dd0*0.5_wp .and. l .ne. 1 .and. usetol .and. rptflag )  then
      if ( .not. present(success) ) then
        write(0,*) 'Warning (LINSOLVE): slow convergence to tolerance'
      end if
      locsuccess = .false.
      rptflag = .false.
    end if

    if ( l .gt. maxiter )  rptflag = .false.

    dd0 = dd1

  end do

  if ( present(resid) )  resid = maxval(abs(bb))
  if ( present(success) ) success = locsuccess .and. LUsuccess

  deallocate( rLU, indx )

end function rlinsolve


!-----------------------------------------------------------------------
!
!  subroutine rlinsolve_backsub(b)
!
!  Subroutine to solve Ax=b, based on back-substitution of
!  LU decomposition.
!  Intended to be called indirectly by linsolve_iter for additional 
!  iterative refinements.
!  Based on dgesl routine in the LINPACK library.
!
!-----------------------------------------------------------------------

subroutine rlinsolve_backsub(b)

  implicit none

  real(wp), dimension(:), intent(inout) :: b

  integer :: n, k
  real(wp) :: x

  n = size(b, dim=1)
  if ( n .eq. 1 ) then; b(1) = b(1)/rLU(1,1); return; end if

  ! undo pivot permutations
  ! solve L*y=b
  do k = 1, n-1
    x = b(indx(k))
    if ( indx(k) .ne. k ) then
      b(indx(k)) = b(k); b(k) = x
    end if
    if (x.ne.0) b((k+1):n) = b((k+1):n) + x*rLU((k+1):n,k)
  end do !k
  ! solve Ux=y
  do k = n, 2, -1
    b(k) = b(k)/rLU(k,k)
    if (b(k).ne.0) b(1:(k-1)) = b(1:(k-1)) - b(k)*rLU(1:(k-1),k)
  end do
  b(1) = b(1)/rLU(1,1)

end subroutine rlinsolve_backsub


!-----------------------------------------------------------------------
!       
!  function rlinsolve_multi(A, b, optional resid, optional success,
!                           optional tol, optional niter) 
!       
!  Produces solutions to Ax=b for set of column vectors b, by making an
!  LU decomposition and backsubstituting with iteration.
!  Here, b (the 2D array of column vectors) is unmodifed.
!  This is a Fortran 90 implementation of the 'unsym' routine in
!  Wilkinson and Reinsch.
! 
!  The variable success is an optional output; when the system is
!  ill-conditioned, the presence of this variable ensures that
!  execution will continue (with the output value set to false).
!
!  Other dummy variable act as in rlinsolve.
!
!-----------------------------------------------------------------------

function rlinsolve_multi(A, b, resid, success, tol, niter) result (x)

  implicit none

  real(wp), dimension(:,:), intent(in   ) :: A
  real(wp), dimension(:,:), intent(in   ) :: b 
  real(wp),                 intent(  out), optional :: resid
  logical,                  intent(  out), optional :: success
  real(wp),                 intent(in   ), optional :: tol
  integer,                  intent(in   ), optional :: niter

  real(wp), dimension(size(b, dim=1), size(b, dim=2)) :: x
  logical rptflag, LUsuccess, locsuccess, usetol
  integer :: n, r, j, l, maxiter
  real(wp) :: xmax, bbmax, dd0, dd1, tolfac
  real(wp), dimension(size(b, dim=1), size(b, dim=2)) :: bb

  n = size(A, dim=1)
  r = size(b, dim=2)
  if ( n .ne. size(b, dim=1) ) then
    write(0,*) 'Error (LINSOLVE_MULTI): input array of wrong length'
    stop
  end if
  if ( present(resid) )  resid = 0
  
  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (LINSOLVE_MULTI): tol and niter cannot both be present'
    stop
  end if
  usetol = .false.
  tolfac = 0
  maxiter = 100
  if ( present(tol) ) then
    usetol = .true.
    if ( tol .gt. 0 ) then
      tolfac = tol
    else if ( tol .eq. 0 ) then
      tolfac = 2*n*epsilon(1.0_wp)
    else
      tolfac = 2*epsilon(1.0_wp)
    end if
  else if ( present(niter) ) then
    if ( niter .lt. 0 ) then
      write(0,*) 'Error (LINSOLVE): illegal value for niter'
      stop
    end if
    maxiter = niter
  else
    maxiter = 0
  end if

  locsuccess = .true.
  call rLU_decomp(A, success=LUsuccess)
  if ( .not. present(success) .and. .not. LUsuccess ) then
    write(0,*) 'Warning (LINSOLVE): problems with call to LUdecomp'
  end if

  ! if no iteration, skip all the residuals overhead
  if ( maxiter .eq. 0 .and. .not. present(resid) ) then
    success=.true.;  x = b
    do j = 1, r; call rlinsolve_backsub(x(:,j)); end do
    deallocate( rLU, indx )
    return
  end if

  l = 0
  dd0 = 0
  rptflag = .true.
  bb = b
  x = 0

  do while ( rptflag )
    do j = 1, r; call rlinsolve_backsub(bb(:,j)); end do
    l = l + 1
    rptflag = .false.
    dd1 = 0
    x = x + bb
    xmax = maxval(abs(x)); bbmax = maxval(abs(bb))
    bb = b - matmul(A, x)
    if ( bbmax .gt. dd1*xmax )  dd1 = bbmax / xmax
    if ( bbmax .gt. tolfac*xmax ) rptflag = .true.

    if ( dd1 .gt. dd0*0.5_wp .and. l .ne. 1 .and. usetol .and. rptflag )  then
      if ( .not. present(success) .and. locsuccess ) then
        write(0,*) 'Warning (LINSOLVE_MULTI): slow convergence to tolerance'
      end if
      locsuccess = .false.
      rptflag = .false.
    end if

    if ( l .gt. maxiter )  rptflag = .false.

    dd0 = dd1

  end do

  if ( present(resid) )  resid = maxval(abs(bb))
  if ( present(success) ) success = locsuccess .and. LUsuccess
  deallocate( rLU, indx )

end function rlinsolve_multi


!-----------------------------------------------------------------------
!
!  function rmatrix_inverse(A, optional resid, optional success,
!                           optional tol, optional niter)
!
!  Interface to rlinsolve to compute the inverse of the matrix A.
!
!  Other dummy variables as in rlinsolve
!
!-----------------------------------------------------------------------

function rmatrix_inverse(A, resid, success, tol, niter) result (inv)

  implicit none

  real(wp), dimension(:,:), intent(in   ) :: A
  real(wp),                 intent(  out), optional :: resid
  logical,                  intent(  out), optional :: success
  real(wp),                 intent(in   ), optional :: tol
  integer,                  intent(in   ), optional :: niter

  real(wp), dimension(size(A, dim=1), size(A, dim=1)) :: identity, inv
  integer :: n, j
  logical :: locsuccess

  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (MATRIX_INVERSE): tol and niter cannot both be present'
    stop
  end if

  n = size(A, dim=1)
  identity = 0
  forall( j = 1:n ); identity(j,j) = 1; end forall
  if ( present(tol) .and. present(resid) ) then
    inv = rlinsolve_multi(A, identity, resid=resid, success=locsuccess, tol=tol)
  else if ( present(niter) .and. present(resid) ) then
    inv = rlinsolve_multi(A, identity, resid=resid, success=locsuccess, &
                          niter=niter)
  else if ( present(tol) .and. present(resid) ) then
    inv = rlinsolve_multi(A, identity, resid=resid, success=locsuccess, tol=tol)
  else if ( present(tol) ) then
    inv = rlinsolve_multi(A, identity, success=locsuccess, tol=tol)
  else if ( present(niter) ) then
    inv = rlinsolve_multi(A, identity, success=locsuccess, niter=niter)
  else if ( present(resid) ) then
    inv = rlinsolve_multi(A, identity, resid=resid, success=locsuccess)
  else
    inv = rlinsolve_multi(A, identity, success=locsuccess)
  end if

  if ( present(success) ) then
    success = locsuccess
  else if ( .not. locsuccess ) then
    write(0,*) 'Warning (MATRIX_INVERSE): problems in underlying routines'
  end if

end function rmatrix_inverse


!-----------------------------------------------------------------------
!
!  subroutine zLU_decomp(A, optional success)
!
!  Complex version of rLU_decomp.
!
!-----------------------------------------------------------------------

subroutine zLU_decomp(A, success)

  implicit none

  complex(wp), dimension(:,:), intent(in) :: A
  logical, optional, intent(out) :: success

  integer :: n, j, k, l
  complex(wp) :: x
  logical :: locsuccess

  n = size(A, dim=1)
  if ( n .ne. size(A, dim=2 ) ) then
    write(0,*) 'Error (LU_DECOMP): input matrix not square'
    stop
  end if

  allocate( zLU(n, n), indx(n) )
  zLU = A
  if ( n .eq. 1 )  return

  locsuccess = .true.; if ( present(success) )  success = .true.
  do k = 1, n-1
    l = maxloc(abs(zLU(k:n,k)), dim=1) + k - 1 ! pivot row
    indx(k) = l
    if ( zLU(l,k) .ne. 0 ) then
      !if ( l .ne. k ) then
        x = zLU(l,k); zLU(l,k) = zLU(k,k); zLU(k,k) = x
      !end if
      zLU((k+1):n,k) = zLU((k+1):n,k) * (-1/zLU(k,k))
      do j = k+1, n
        x = zLU(l,j)
        !if ( l .ne. k ) then
          zLU(l,j) = zLU(k,j)
          zLU(k,j) = x
        !end if
        if (x.ne.0) zLU((k+1):n,j) = zLU((k+1):n,j) + x * zLU((k+1):n,k)
      end do !j
    else
      locsuccess = .false.
    end if
  end do !k
  indx(n) = n

  if ( .not. locsuccess .or. zLU(n,n) .eq. 0 ) then
    if ( present(success) ) then
      success = .false.
    else
      write(0,*) 'Warning (LU_DECOMP): input matrix singular'
    end if
  end if

end subroutine zLU_decomp


!-----------------------------------------------------------------------
!
!  complex(wp) function zdeterminant(A, optional success)
!
!  Complex version of rdeterminant.
!
!-----------------------------------------------------------------------

function zdeterminant(A, success)

  implicit none

  complex(wp), dimension(:,:), intent(in) :: A
  logical, intent(out), optional :: success
  complex(wp), dimension(2) :: zdeterminant
  complex(wp) :: zd1, zd2
  integer :: j, n

  n = size(A, dim=1)
  if ( present(success) ) then
    call zLU_decomp(A, success=success)
  else
    call zLU_decomp(A)
  end if
  zd1 = 1; zd2 = 0
  do j = 1, n
    if ( indx(j) .ne. j )  zd1 = -zd1
    zd1 = zLU(j,j)*zd1
    if ( zd1 .eq. 0 )  exit
    do while ( abs(zd1) .lt. 1 ); zd1 = 10*zd1; zd2 = zd2-1; end do
    do while ( abs(zd1) .gt. 10 ); zd1 = 0.1_wp*zd1; zd2 = zd2+1; end do
  end do
  zdeterminant = (/ zd1, zd2 /)
  deallocate( zLU, indx )

end function zdeterminant


!-----------------------------------------------------------------------
!
!  function zlinsolve(A, b, optional resid, optional success,
!                     optional tol, optional niter)
!
!  Complex version of rlinsolve.
!
!-----------------------------------------------------------------------

function zlinsolve(A, b, resid, success, tol, niter) result (x)
  
  implicit none
  
  complex(wp), dimension(:,:), intent(in   ) :: A
  complex(wp), dimension(:),   intent(in   ) :: b
  real(wp),                    intent(  out), optional :: resid
  logical,                     intent(  out), optional :: success
  real(wp),                    intent(in   ), optional :: tol
  integer,                     intent(in   ), optional :: niter
  
  complex(wp), dimension(size(b)) :: x
  logical :: rptflag, LUsuccess, locsuccess, usetol
  integer :: n, l, maxiter
  real(wp) :: xmax, bbmax, dd0, dd1, tolfac
  complex(wp), dimension(size(b)) :: bb

  n = size(A, dim=1)
  if ( n .ne. size(b, dim=1) ) then
    write(0,*) 'Error (LINSOLVE): input vector of wrong length'
    stop
  end if

  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (LINSOLVE): tol and niter cannot both be present'
    stop
  end if
  usetol = .false.
  tolfac = 0
  maxiter = 100
  if ( present(tol) ) then
    usetol = .true.
    if ( tol .gt. 0 ) then
      tolfac = tol
    else if ( tol .eq. 0 ) then
      tolfac = 2*n*epsilon(1.0_wp)
    else
      tolfac = 2*epsilon(1.0_wp)
    end if
  else if ( present(niter) ) then
    if ( niter .lt. 0 ) then
      write(0,*) 'Error (LINSOLVE): illegal value for niter'
      stop
    end if
    maxiter = niter
  else
    maxiter = 0
  end if

  locsuccess = .true.
  call zLU_decomp(A, success=LUsuccess)
  if ( .not. present(success) .and. .not. LUsuccess ) then
    write(0,*) 'Warning (LINSOLVE): problems with call to LUdecomp'
  end if

  l = 0
  dd0 = 0
  rptflag = .true.
  bb = b
  x = 0

  do while ( rptflag )
    call zlinsolve_backsub(bb)
    l = l + 1
    rptflag = .false.
    x = x + bb
    xmax = maxval(abs(x)); bbmax = maxval(abs(bb))
    bb = b - matmul(A, x)
    dd1 = bbmax / xmax
    if ( bbmax .gt. tolfac*xmax )  rptflag = .true.

    if ( dd1 .gt. dd0*0.5_wp .and. l .ne. 1 .and. usetol .and. rptflag )  then
      if ( .not. present(success) ) then
        write(0,*) 'Warning (LINSOLVE): slow convergence to tolerance'
      end if
      locsuccess = .false.
      rptflag = .false.
    end if

    if ( l .gt. maxiter )  rptflag = .false.

    dd0 = dd1

  end do

  if ( present(resid) )  resid = maxval(abs(bb))
  if ( present(success) ) success = locsuccess .and. LUsuccess

  deallocate( zLU, indx )

end function zlinsolve


!-----------------------------------------------------------------------
!
!  subroutine zlinsolve_backsub(b)
!
!  Complex version of rlinsolve_backsub.
!
!-----------------------------------------------------------------------

subroutine zlinsolve_backsub(b)

  implicit none

  complex(wp), dimension(:), intent(inout) :: b

  integer :: n, k
  complex(wp) :: x

  n = size(b, dim=1)
  if ( n .eq. 1 ) then; b(1) = b(1)/zLU(1,1); return; end if

  ! undo pivot permutations
  ! solve L*y=b
  do k = 1, n-1
    x = b(indx(k))
    if ( indx(k) .ne. k ) then
      b(indx(k)) = b(k); b(k) = x
    end if
    if (x.ne.0) b((k+1):n) = b((k+1):n) + x*zLU((k+1):n,k)
  end do !k
  ! solve Ux=y 
  do k = n, 2, -1
    b(k) = b(k)/zLU(k,k)
    if (b(k).ne.0) b(1:(k-1)) = b(1:(k-1)) - b(k)*zLU(1:(k-1),k)
  end do 
  b(1) = b(1)/zLU(1,1)

end subroutine zlinsolve_backsub



!-----------------------------------------------------------------------
!
!  function zlinsolve_multi(A, b, optional resid, optional success,
!                           optional tol, optional niter)
!
!  Complex version of rlinsolve_multi.
!
!-----------------------------------------------------------------------

function zlinsolve_multi(A, b, resid, success, tol, niter) result (x)

  implicit none

  complex(wp), dimension(:,:), intent(in   ) :: A
  complex(wp), dimension(:,:), intent(in   ) :: b 
  real(wp),                    intent(  out), optional :: resid
  logical,                     intent(  out), optional :: success
  real(wp),                    intent(in   ), optional :: tol
  integer,                     intent(in   ), optional :: niter

  complex(wp), dimension(size(b, dim=1), size(b, dim=2)) :: x
  logical rptflag, LUsuccess, locsuccess, usetol
  integer :: n, r, j, l, maxiter
  real(wp) :: xmax, bbmax, dd0, dd1, tolfac
  complex(wp), dimension(size(b, dim=1), size(b, dim=2)) :: bb
  
  n = size(A, dim=1)
  r = size(b, dim=2)
  if ( n .ne. size(b, dim=1) ) then
    write(0,*) 'Error (LINSOLVE_MULTI): input array of wrong length'
    stop 
  end if
  if ( present(resid) )  resid = 0
      
  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (LINSOLVE_MULTI): tol and niter cannot both be present'
    stop 
  end if
  usetol = .false.
  tolfac = 0
  maxiter = 100
  if ( present(tol) ) then
    usetol = .true.
    if ( tol .gt. 0 ) then
      tolfac = tol
    else if ( tol .eq. 0 ) then
      tolfac = 2*n*epsilon(1.0_wp)
    else
      tolfac = 2*epsilon(1.0_wp)
    end if
  else if ( present(niter) ) then
    if ( niter .lt. 0 ) then
      write(0,*) 'Error (LINSOLVE): illegal value for niter'
      stop
    end if
    maxiter = niter
  else
    maxiter = 0
  end if

  locsuccess = .true.
  call zLU_decomp(A, success=LUsuccess)
  if ( .not. present(success) .and. .not. LUsuccess ) then
    write(0,*) 'Warning (LINSOLVE): problems with call to LUdecomp'
  end if

  ! if no iteration, skip all the residuals overhead
  if ( maxiter .eq. 0 .and. .not. present(resid) ) then
    success=.true.;  x = b
    do j = 1, r; call zlinsolve_backsub(x(:,j)); end do
    deallocate( zLU, indx )
    return
  end if

  l = 0
  dd0 = 0
  rptflag = .true.
  bb = b
  x = 0

  do while ( rptflag )
    do j = 1, r; call zlinsolve_backsub(bb(:,j)); end do
    l = l + 1
    rptflag = .false.
    dd1 = 0
    x = x + bb
    xmax = maxval(abs(x)); bbmax = maxval(abs(bb))
    bb = b - matmul(A, x)
    if ( bbmax .gt. dd1*xmax )  dd1 = bbmax / xmax
    if ( bbmax .gt. tolfac*xmax ) rptflag = .true.

    if ( dd1 .gt. dd0*0.5_wp .and. l .ne. 1 .and. usetol .and. rptflag )  then
      if ( .not. present(success) .and. locsuccess ) then
        write(0,*) 'Warning (LINSOLVE_MULTI): slow convergence to tolerance'
      end if
      locsuccess = .false.
      rptflag = .false.
    end if

    if ( l .gt. maxiter )  rptflag = .false.

    dd0 = dd1

  end do

  if ( present(resid) )  resid = maxval(abs(bb))
  if ( present(success) ) success = locsuccess .and. LUsuccess
  deallocate( zLU, indx )

end function zlinsolve_multi



!-----------------------------------------------------------------------
!
!  function zmatrix_inverse(A, optional resid, optional success,
!                           optional tol, optional niter)
!
!  Complex version of rmatrix_inverse.
!
!-----------------------------------------------------------------------

function zmatrix_inverse(A, resid, success, tol, niter) result (inv)

  implicit none

  complex(wp), dimension(:,:), intent(in   ) :: A
  real(wp),                    intent(  out), optional :: resid
  logical,                     intent(  out), optional :: success
  real(wp),                    intent(in   ), optional :: tol
  integer,                     intent(in   ), optional :: niter

  complex(wp), dimension(size(A, dim=1), size(A, dim=1)) :: identity, inv
  integer :: n, j
  logical :: locsuccess

  if ( present(tol) .and. present(niter) ) then
    write(0,*) 'Error (MATRIX_INVERSE): tol and niter cannot both be present'
    stop
  end if

  n = size(A, dim=1)
  identity = 0
  forall( j = 1:n ); identity(j,j) = 1; end forall
  if ( present(tol) .and. present(resid) ) then
    inv = zlinsolve_multi(A, identity, resid=resid, success=locsuccess, tol=tol)
  else if ( present(niter) .and. present(resid) ) then
    inv = zlinsolve_multi(A, identity, resid=resid, success=locsuccess, &
                          niter=niter)
  else if ( present(tol) .and. present(resid) ) then
    inv = zlinsolve_multi(A, identity, resid=resid, success=locsuccess, tol=tol)
  else if ( present(tol) ) then
    inv = zlinsolve_multi(A, identity, success=locsuccess, tol=tol)
  else if ( present(niter) ) then
    inv = zlinsolve_multi(A, identity, success=locsuccess, niter=niter)
  else if ( present(resid) ) then
    inv = zlinsolve_multi(A, identity, resid=resid, success=locsuccess)
  else
    inv = zlinsolve_multi(A, identity, success=locsuccess)
  end if

  if ( present(success) ) then
    success = locsuccess
  else if ( .not. locsuccess ) then
    write(0,*) 'Warning (MATRIX_INVERSE): problems in underlying routines'
  end if

end function zmatrix_inverse


end module linsol
