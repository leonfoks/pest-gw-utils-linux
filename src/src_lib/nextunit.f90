integer function nextunit()

! -- Function nextunit determines the lowest unit number available for
! -- opening.

! -- Revision history:-
!       June-November, 1995: version 1.

	logical::lopen

	do nextunit=10,100
	  inquire(unit=nextunit,opened=lopen)
	  if(.not.lopen) return
	end do
	write(6,10)
10      format(' *** No more unit numbers to open files ***')
	stop

end function nextunit
 