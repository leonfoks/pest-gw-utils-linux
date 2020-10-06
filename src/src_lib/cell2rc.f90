subroutine cell2rc(icellno,irow,icol,gridspec)

! -- Subroutine cell2rc converts a cell number to a row and column number.

! -- Arguments are as follows:-
!       icellno:     cell number supplied to cell2rc
!       irow,icol:   row and column numbers calculated by cell2rc
!       gridspec:    a derived type variable holding the grid specifications.

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	implicit none

	integer, intent(in)             :: icellno
	integer, intent(out)            :: irow,icol
	type (modelgrid), intent(in)    :: gridspec

	irow=(icellno-1)/gridspec%ncol+1
	icol=icellno-((irow-1)*gridspec%ncol)
	return

end subroutine cell2rc
 