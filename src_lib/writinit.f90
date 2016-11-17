subroutine write_initial_message(leadspace,endspace)

! -- Subroutine write_initial_message writes the first part of a two-part
!    error message.

! -- Arguments are as follows:-
!       leadspace:  if "yes" write a leading space to message
!       endspace:   if "no" write a trailing space to message

! -- Revision history:-
!       June-November, 1995: version1. 

	use defn
	use inter

	character (len=*), intent(in), optional :: leadspace,endspace
	character (len=500)     :: atemp

	atemp=amessage
	amessage=initial_message
	call write_message(leadspace=leadspace,endspace=endspace)
	amessage=atemp

end subroutine write_initial_message
 