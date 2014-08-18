!-----------------------------------------------------------------------
!
!     utilities.f90
!
!     $Id: utilities.f90,v 2.1 2003/02/14 03:02:25 dsteck Exp $
!
!     Utilities for handling strings and input parameters
!    
!     Requires 'wp' to define the real kind, to be read from the
!       'globals' module.
!
!-----------------------------------------------------------------------


module utilities


use globals, only : wp


contains


subroutine initialize_module_utilities()
  character(127), parameter :: &
    RCS_ID = "@(#)$Id: utilities.f90,v 2.1 2003/02/14 03:02:25 dsteck Exp $"
  write(0,*) trim(RCS_ID)
end subroutine initialize_module_utilities


!=======================================================================
!
!  integer function  indlnb(s)
!
!  I: s                 character string
!
!  O: indlnb            integer
!
!     Returns index of last non-blank character in S.
!     (This routine written by Matthew Choptuik.)
!
!=======================================================================

function indlnb(s)
  character(*) :: s
  integer :: indlnb
  do indlnb = len(s) , 1 , -1
    if( s(indlnb:indlnb) .ne. ' ' ) return
  end do
  indlnb = 0
end function indlnb


!=======================================================================
!
!  integer function  indfnb(s)
!
!  I: s                 character string
!
!  O: indfnb            integer
!
!     Returns index of first non-blank character in S.
!     (Based on routine written by Matthew Choptuik.)
!
!=======================================================================

function indfnb(s)
  character(*) :: s
  integer :: indfnb
  do indfnb = 1, len(s) , 1
    if( s(indfnb:indfnb) .ne. ' ' ) return
  end do
  indfnb = 0
end function indfnb


!=======================================================================
!
!  integer function  getu()
!
!  I: none
!
!  O: getu          integer
!
!     Returns first unit number .ge. umin not attached to
!     a file.
!     (This routine written by Matthew Choptuik, translated to f90
!        by D. Steck)
!
!=======================================================================

integer function getu()   

   implicit none

   integer, parameter :: umin = 10
   integer, parameter :: umax = 99
   integer u
   logical opened

   getu = -1
   do u = umin, umax
      inquire(unit=u,opened=opened)
      if( .not. opened ) then
         getu = u
         return
      end if
   end do
   write(0,*) 'Error (GETU): no available unit number' 
   stop

end function getu


!=======================================================================
!
!  integer function  s2i(s)
!
!  I: s             character string
!
!  O: s2i           integer
!
!     Converts string to integer.
!     (This routine written by Matthew Choptuik.)
!
!=======================================================================

function s2i(s)

   implicit none

   character(*) :: s     
   integer :: s2i

   integer, parameter :: i4_never  = -2000000000
   character(32) ::  buffer
   integer       :: lens

   integer, parameter :: default = 0

   lens = indlnb(s)
   if( lens .gt. 99 ) then
      write(0,*) 'Warning (S2I): String too long for conversion.'
      s2i = default
   else 
     if( lens .le. 9 ) then
       write(buffer,100) lens
     else 
       write(buffer,101) lens
     end if
100  format('(I',i1,')')
101  format('(I',i2,')')
     read(s,fmt=buffer,end=900,err=900) s2i
   end if

   return

900 s2i = i4_never
    return

end function s2i


!=======================================================================
!
!  real(wp) function  s2r(s) 
!       
!  I: s             character string
!  
!  O: s2r           real(wp)
!
!     Converts string to real(kind=wp), where wp is defined in
!     the module globals.
!     (This routine written by Matthew Choptuik, ported to F90 by
!     D. Steck.)
!
!=======================================================================


function s2r(s)

  implicit none

  real(wp) :: s2r
  character(*) :: s

  real(wp), parameter :: r_never = -1.0e-20_wp
  character(32) :: buffer    
  integer       :: lens

  real(wp), parameter :: default = 0.0_wp

  lens = indlnb(s)
  if( lens .gt. 99 ) then
     write(0,*) 'Warning (S2R): String too long for conversion.'
     s2r = default
  else
     if( lens .le. 9 ) then
        write(buffer,100) lens     
     else
        write(buffer,101) lens
     end if    
100 format('(G',i1,'.0)')
101 format('(G',i2,'.0)')
     read(s,fmt=buffer,end=900,err=900) s2r
  end if

  return

900 s2r = r_never
    return

end function s2r


!=======================================================================
!
!  subroutine  get_real_param(p_file, p_name, param, trace)
!     
!  I: p_file            character string          parameter file name
!     p_name            character string          parameter name
!     trace             logical                   sets tracing mode
!
!  O: param             real(wp)
!
!     Parses through parameter file (p_file), with lines of form
!               p_name := param
!       and returns the value param corresponding to the
!       parameter name in p_name.
!
!=======================================================================

subroutine  get_real_param(p_file, p_name, param, trace)

  implicit none     

  character(*)      p_file, p_name    
  real(wp)          param
  logical           trace
  integer           fp   
  character(128)    buffer
  integer           length, err, locat
  integer           frst, lst

  fp = getu()  

  open( unit=fp, file=p_file(1:indlnb(p_file)), status='old', err=87 )

85 read( fp, '(a)', iostat=err, end=88 ) buffer

  if ( err .eq. 0 ) then

!   find length of line
    length = indlnb(buffer)

!   check to see if this is correct line in file
    if ( length .le. 4 ) go to 85
    if ( index(buffer(:length), p_name) .eq. 0 )  go to 85

!   locate position of ':=' assignment string
    locat = index(buffer, ':=')
    if ( locat .eq. 0 ) go to 85

! check to make sure this is correct line in file
    frst = indfnb(buffer)
    lst = indlnb(buffer(1:(locat-1)))
    if ( len(p_name) .ne. len(buffer(frst:lst)) ) go to 85

!   find position of parameter (skip white space after :=)      
    locat = locat + 2
86  if ( buffer(locat:locat) .eq. ' ' ) then
      locat = locat + 1
      go to 86
    end if

!   convert character data to real
    param = s2r(buffer(locat:length))

!   if tracing is enabled, write results to stdout
    if ( trace .eqv. .true. ) then    
      write(0,*) 'parameter ',p_name,' set to ',param
    end if

    close(fp)
    return

  end if

87 continue
   write(0,*)  'Error while reading parameter file "', p_file,'"'
   stop        

88 continue
   write(0,*) 'Parameter "',p_name,'" not found in file "', &
               p_file(1:indlnb(p_file)),'"'    
   stop

end subroutine get_real_param


!=======================================================================
!
!  subroutine  get_int_param(p_file, p_name, param, trace)
!
!  I: p_file            character string          parameter file name
!     p_name            character string          parameter name
!     trace             logical                   sets tracing mode
!
!  O: param             integer
! 
!     Parses through parameter file (p_file), with lines of form
!               p_name := param      
!       and returns the value param corresponding to the  
!       parameter name in p_name.
!   
!=======================================================================

subroutine  get_int_param(p_file, p_name, param, trace)

  implicit none

  character(*)      p_file, p_name  
  integer           param
  logical           trace
  integer           fp
  character(128)    buffer
  integer           length, err, locat
  integer           frst, lst

  fp = getu()

  open( unit=fp, file=p_file(1:indlnb(p_file)), status='old', err=87 )

85 read( fp, '(a)', iostat=err, end=88 ) buffer

  if ( err .eq. 0 ) then

!   find length of line
    length = indlnb(buffer)

!   check to see if this is correct line in file
    if ( length .le. 4 ) go to 85
    if ( index(buffer(:length), p_name) .eq. 0 )  go to 85

!   locate position of ':=' assignment string
    locat = index(buffer, ':=')
    if ( locat .eq. 0 ) go to 85

!   check to make sure this is correct line in file
    frst = indfnb(buffer)
    lst = indlnb(buffer(1:(locat-1)))
    if ( len(p_name) .ne. len(buffer(frst:lst)) ) go to 85

!   find position of parameter (skip white space after :=)      
    locat = locat + 2
86  if ( buffer(locat:locat) .eq. ' ' ) then
      locat = locat + 1
      go to 86
    end if

!   convert character data to integer
    param = s2i(buffer(locat:length))

!   if tracing is enabled, write results to stdout
    if ( trace .eqv. .true. ) then  
      write(0,*) 'parameter ',p_name,' set to ',param
    end if

    close(fp)
    return

  end if

87 continue
   write(0,*)  'Error while reading parameter file "', p_file,'"'
   stop      

88 continue
   write(0,*) 'Parameter "',p_name,'" not found in file "', &
               p_file(1:indlnb(p_file)),'"'    
   stop

end subroutine get_int_param


!=======================================================================
!
!  subroutine  get_str_param(p_file, p_name, param, trace)
!
!  I: p_file            character string          parameter file name
!     p_name            character string          parameter name
!     trace             logical                   sets tracing mode
!
!  O: param             character string
!   
!     Parses through parameter file (p_file), with lines of form
!               p_name := param      
!       and returns the value param corresponding to the  
!       parameter name in p_name.
!   
!=======================================================================

subroutine  get_str_param(p_file, p_name, param, trace)

  implicit none

  character(*)      p_file
  character(*)      p_name
  character(*)      param
  logical           trace
  integer           fp
  character(256)    buffer
  integer           length, err, locat
  integer           frst, lst

  fp = getu()

  open( unit=fp, file=p_file, status='old', err=87 )

85 read( fp, '(a)', iostat=err, end=88 ) buffer
   if ( err .eq. 0 ) then

!    find length of line
     length = indlnb(buffer)

!    check to see if this is correct line in file
     if ( length .le. 4 ) go to 85
     if ( index(buffer(1:length), p_name) .eq. 0 )  go to 85

!    locate position of ':=' assignment string
     locat = index(buffer, ':=')
     if ( locat .eq. 0 ) go to 85

!    check to make sure this is correct line in file
     frst = indfnb(buffer)
     lst = indlnb(buffer(1:(locat-1)))
     if ( len(p_name) .ne. len(buffer(frst:lst)) ) go to 85

!    find position of parameter (skip white space after :=)
     locat = locat + 2
86   if ( buffer(locat:locat) .eq. ' ' ) then
       locat = locat + 1
       go to 86
     end if

     param = buffer(locat:length)

!    if tracing is enabled, write results to stdout
     if ( trace .eqv. .true. ) then
       write(0,*) 'parameter ',p_name,' set to ', param(1:indlnb(param))
     end if

     close(fp)
     return

   end if

87 continue
   write(0,*)  'Error while reading parameter file "', p_file,'"'
   stop

88 continue  
   write(0,*) 'Parameter "',p_name,'" not found in file "', &
               p_file(1:indlnb(p_file)),'"' 
   stop

end subroutine get_str_param


!=======================================================================
!
!  subroutine  get_log_param(p_file, p_name, param, trace)
!
!  I: p_file            character string          parameter file name
!     p_name            character string          parameter name
!     trace             logical                   sets tracing mode
!
!  O: param             logical  
!
!     Parses through parameter file (p_file), with lines of form
!               p_name := param
!       and returns the value param corresponding to the
!       parameter name in p_name.
!
!     This routine differs from the real and integer routines in
!       that it does not terminate the program if it doesn't
!       find the parameter;  rather, it defaults to false if
!       the parameter isn't found, so that unwanted options
!       need not be mentioned.
!
!=======================================================================

subroutine  get_log_param(p_file, p_name, param, trace)

  implicit none

  character(*)      p_file
  character(*)      p_name
  logical           param
  logical           trace
  integer           fp
  character(128)    buffer
  integer           length, err, locat
  integer           frst, lst

  fp = getu()

  open( unit=fp, file=p_file, status='old', err=87 )

85 read( fp, '(a)', iostat=err, end=88 ) buffer
   if ( err .eq. 0 ) then

!  find length of line
   length = indlnb(buffer)

!  check to see if this is correct line in file
   if ( length .le. 4 ) go to 85
   if ( index(buffer(:length), p_name) .eq. 0 )  go to 85       

!  locate position of ':=' assignment string
   locat = index(buffer, ':=')   

   if ( locat .eq. 0 ) go to 85

!  check to make sure this is correct line in file      
   frst = indfnb(buffer)
   lst = indlnb(buffer(1:(locat-1)))
   if ( len(p_name) .ne. len(buffer(frst:lst)) ) go to 85       

!  find position of parameter (skip white space after :=)   
   locat = locat + 2
86 if ( buffer(locat:locat) .eq. ' ' ) then
     locat = locat + 1
     go to 86
   end if

   if ( buffer(locat:length) .eq. 'true' ) then
     param = .true.
   else if ( buffer(locat:length) .eq. 'false' ) then
     param = .false.      
   else
     write(0,*) 'Illegal value for parameter ', buffer(locat:length)
     stop
   end if

!  if tracing is enabled, write results to stdout
   if ( trace .eqv. .true. ) then
     write(0,*) 'parameter ',p_name,' set to ', param
   end if

   close(fp)
   return    

  end if

87 continue
   write(0,*)  'Error while reading parameter file "', p_file,'"'
   stop

88 continue
   if ( trace .eqv. .true. ) then
     write(0,*) 'Parameter "',p_name,'" not found in file "', &
                 p_file(1:indlnb(p_file)),'"; defaulting to FALSE'
   end if
   close(fp)
   param = .false.
   return

end subroutine get_log_param


end module utilities
