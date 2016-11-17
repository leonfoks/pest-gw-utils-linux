subroutine sub_error(subname)

! -- Subroutine sub_error names the subroutine causing a run-time error.

! -- Arguments are as follows:-
!       subname:  name of offending subroutine

! -- Revision history:-
!       June-November, 1995: version 1.

	character (len=*)               ::subname

	write(6,10) trim(subname)
10      format(/,' *** PROGRAMMING ERROR CALLING SUBROUTINE ',a,' ***')
	stop

end subroutine sub_error
 