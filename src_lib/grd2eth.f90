subroutine grid2earth(east,north,gridspec)

! -- Subroutine grid2earth calls subroutine transform_to_earth for each
!    grid cell centre in order to express cell centre coordinates in a system
!    whose origin is at the top left corner of the grid and whose x-direction
!    is oriented easterly.

! -- Arguments are as follows:-
!       east, north:  arrays holding grid cell centre x and y coordinates
!       gridspec:     defined type holding grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	real, intent(inout),dimension(:,:)      :: east,north
	type(modelgrid), intent(in)             :: gridspec
	integer                                 :: icol,irow
	real                                    :: earth_east,earth_north


	do irow=1,size(gridspec%delc)
	  do icol=1,size(gridspec%delr)
	    call transform_to_earth(earth_east,earth_north, &
	    east(icol,irow),north(icol,irow),gridspec)
	    east(icol,irow)=earth_east
	    north(icol,irow)=earth_north
	  end do
	end do

	return

end subroutine grid2earth
 