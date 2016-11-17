subroutine corner(ind,e,n,icol,irow,gridspec,ecg,ncg)

! -- Subroutine corner calculates the coordinates of one of the four corners
!    of a grid cell. The coordinate system has its origin at the top left of
!    the finite difference grid; its x-direction is oriented easterly.

! -- Arguments are as follows:-
!       ind:       indicates corner for which coordinates are required [1-4]
!       e,n:       computed x and y coordinates of corresponding corner
!       icol,irow: row and column numbers of cell for which corner coordinates
!                  are required
!       gridspec:  defined type holding grid specifications
!       ecg, ncg:  array holding grid cell centre x and y coordinates
!                  expressed in coordinate system for which the x-direction is
!                  oriented easterly and the origin is at the top left corner
!                  of the grid.

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(in)                     :: ind
	real, intent(out)                       :: e,n
	integer, intent(in)                     :: icol,irow
	type (modelgrid), intent(in)            :: gridspec
	real, dimension(:,:), intent(in)        :: ecg,ncg
	real                                    :: dr,dc,eg,ng

	dr=gridspec%delr(icol)/2.0
	dc=gridspec%delc(irow)/2.0
	eg=ecg(icol,irow)
	ng=ncg(icol,irow)
	if(ind.eq.1) Then
	  e=eg-dr*gridspec%cosang-dc*gridspec%sinang
	  n=ng-dr*gridspec%sinang+dc*gridspec%cosang
	else if(ind.eq.2)Then
	  e=eg+dr*gridspec%cosang-dc*gridspec%sinang
	  n=ng+dr*gridspec%sinang+dc*gridspec%cosang
	else if(ind.eq.3)Then
	  e=eg-dr*gridspec%cosang+dc*gridspec%sinang
	  n=ng-dr*gridspec%sinang-dc*gridspec%cosang
	else
	  e=eg+dr*gridspec%cosang+dc*gridspec%sinang
	  n=ng+dr*gridspec%sinang-dc*gridspec%cosang
	end if

	return
end subroutine corner
 