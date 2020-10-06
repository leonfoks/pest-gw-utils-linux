subroutine close_files

! -- Subroutine close_files closes all open files.

! -- Revision history:-
!       June-November, 1995: version 1.

	integer         :: i,ierr

	do i=10,100
	  close(unit=i,iostat=ierr)
	end do
	return

end subroutine close_files
 