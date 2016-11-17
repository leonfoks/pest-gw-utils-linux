subroutine char_add(astring,achar)

! -- Subroutine CHAR_ADD concatenates a string to the end of another string.

! -- Arguments are as follows:-
!       astring:  a string to the end of which the second string is to be 
!                 concatenated
!       achar:    a string to be joined to the end of the first string.

! -- Revision history:-
!       June-November 1995: version 1.

	character (len=*), intent(inout)        :: astring
	character (len=*), intent(in)           :: achar
	integer                                 :: k

	k=len_trim(astring)+1
	astring(k:)=achar
	return

end subroutine char_add
 