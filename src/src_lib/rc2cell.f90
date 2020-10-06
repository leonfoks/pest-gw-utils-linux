subroutine rc2cell(icellno,irow,icol,gridspec)

! -- Subroutine rc2cell converts a row and column number to a cell number.

! -- Arguments are as follows:-
!       icellno:      cell number calculated by rc2cell
!       irow, icol:   the cell row and column number supplied to rc2cell
!       gridspec:     a derived type variable holding the grid specifications

	use defn
	implicit none

	integer, intent(out)            :: icellno
	integer, intent(in)             :: irow,icol
	type(modelgrid), intent(in)     :: gridspec

	icellno=(irow-1)*gridspec%ncol+icol

	return

end subroutine rc2cell


subroutine rc2cell3d(icellno,irow,icol,ilay,gridspec)

! -- Subroutine rc2cell3d converts a row, column and layer number to a cell number.

! -- Arguments are as follows:-
!       icellno:            cell number calculated by rc2cell3d
!       irow, icol, ilay:   the cell row,column and layer number supplied to rc2cell3d
!       gridspec:           a derived type variable holding the grid specifications

	use defn
	implicit none

	integer, intent(out)            :: icellno
	integer, intent(in)             :: irow,icol,ilay
	type(modelgrid), intent(in)     :: gridspec

	icellno=(irow-1)*gridspec%ncol+icol
	if(ilay.gt.1)then
	  icellno=icellno+gridspec%ncol*gridspec%nrow*(ilay-1)
	end if

	return

end subroutine rc2cell3d
