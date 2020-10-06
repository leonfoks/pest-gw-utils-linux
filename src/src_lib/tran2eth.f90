subroutine transform_to_earth(earth_east,earth_north,grid_east,grid_north,&
			      gridspec)

! -- Subroutine transform_to_earth transforms coordinates from those 
!    pertaining to a coordinate system oriented with its x-direction
!    coincident with the grid row direction and with its origin at the top
!    left corner of the grid, to a coordinate system with its x-direction
!    oriented easterly and with its origin at the top left corner of the grid.

! -- Arguments are as follows:-
!       earth_east, earth_north: coordinates calculated for coordinate
!                                system whose x direction is coincident
!                                with the grid row direction.
!       grid_east, grid_north:   coordinates supplied in local grid
!                                coordinate system
!       gridspec                 defined type holding grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	real, intent(out)               :: earth_east, earth_north
	real, intent(in)                :: grid_east, grid_north
	type(modelgrid), intent(in)     :: gridspec

	earth_east= grid_east*gridspec%cosang-grid_north*gridspec%sinang
	earth_north=grid_east*gridspec%sinang+grid_north*gridspec%cosang

	return

end subroutine transform_to_earth
 