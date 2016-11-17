subroutine casetrans(string,hi_or_lo)

! -- Subroutine casetrans converts a string to upper or lower case.

! -- Arguments are as follows:-
!      string:	  contains the string whose case must be changed
!      hi_or_lo:  must be either 'lo' or 'hi' to indicate
!                 change of case direction.

! -- Revision history:-
!       June-November, 1995: version 1.

	use inter

	character (len=*), intent(inout)        :: string
	character (len=*), intent(in)           :: hi_or_lo
	character                               :: alo, ahi
	integer                                 :: inc,i

	if(hi_or_lo.eq.'lo') then
	  alo='A'; ahi='Z'; inc=iachar('a')-iachar('A')
	else if(hi_or_lo.eq.'hi') then
	  alo='a'; ahi='z'; inc=iachar('A')-iachar('a')
	else
	  call sub_error('CASETRANS')
	endif

	do i=1,len_trim(string)
	  if((string(i:i).ge.alo).and.(string(i:i).le.ahi)) &
	  string(i:i)=achar(iachar(string(i:i))+inc)
	end do

	return

end subroutine casetrans
 