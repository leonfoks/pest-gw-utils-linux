
!*****************************************************************************
! subroutines comprising the generic subroutine CHAR2NUM ------->
!*****************************************************************************


! -- The subroutines comprising char2num convert a string to either an integer,
!    a real number, or a double precision number.

! -- Arguments are as follows:-
!      ifail:   indicates failure if returned as non-zero
!      string:  a character string containing a number
!      num:     an integer (for a2i), real (for a2r), or double precision (for
!               a2d) number extracted from the string.

! -- Revision history:-
!       June-November, 1995: version 1.


subroutine a2i(ifail,string,num)

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: string
	integer, intent(out)            :: num
	character (len=10)              :: afmt

	ifail=0
	afmt='(i    )'
	write(afmt(3:6),'(i4)')len(string)
	read(string,afmt,err=10) num
	return

10      ifail=1
	return

end subroutine a2i


subroutine a2r(ifail,string,num)

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: string
	real, intent(out)               :: num
	character (len=10)              :: afmt

	ifail=0
	afmt='(f    .0)'
	write(afmt(3:6),'(i4)')len(string)
	read(string,afmt,err=10) num
	return

10      ifail=1
	return

end subroutine a2r


subroutine a2d(ifail,string,num)

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: string
	double precision, intent(out)   :: num
	character (len=10)              :: afmt

	ifail=0
	afmt='(f    .0)'
	write(afmt(3:6),'(i4)')len(string)
	read(string,afmt,err=10) num
	return

10      ifail=1
	return

end subroutine a2d

 