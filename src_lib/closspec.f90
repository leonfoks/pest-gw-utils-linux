subroutine close_spec_file(gridspec,ok)

! -- Subroutine close_spec_file closes the grid specification file.

! -- Arguments are as follows:-
!       gridspec:     defined type holding grid specification information 
!                     (including unit number)
!       ok:           if supplied as "yes" a message is written to the screen

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	type (modelgrid)                :: gridspec
	logical                         :: lopen
	character (len=*),optional      :: ok

	inquire(unit=gridspec%specunit, opened=lopen)
	if(lopen) close(unit=gridspec%specunit)
	if(present(ok))then
	  if(ok.eq.'yes')then
	    write(amessage,10) trim(gridspec%specfile)
10          format('  - grid specifications read from file ',a)
	    call write_message
	  endif
	endif
	return

end subroutine close_spec_file
 