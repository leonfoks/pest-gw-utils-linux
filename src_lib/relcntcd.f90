subroutine rel_centre_coords(east,north,gridspec)

! -- Subroutine rel_centre_coords evaluates the x and y coordinates of the 
!    centre of each cell of the finite difference grid relative to a
!    coordinate system where the top left corner of the grid is at (0,0) and 
!    the grid row direction is oriented in an easterly direction.

! -- Arguments are as follows:-
!       east:     array holding cell east coordinates (output only)
!       north:    array holding cell north coordinates (output only)
!       gridspec: defined type holding grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	real, intent(out), dimension(:,:)       :: east,north
	type(modelgrid), intent(in)             :: gridspec
	real                                    :: temp
	integer                                 :: i

	temp=gridspec%delr(1)/2.0
	east(1,:)=temp
	do i=2,size(gridspec%delr)
	  temp=temp+(gridspec%delr(i)+gridspec%delr(i-1))/2.0
	  east(i,:)=temp
	end do

	temp=-gridspec%delc(1)/2.0
	north(:,1)=temp
	do i=2,size(gridspec%delc)
	  temp=temp-(gridspec%delc(i)+gridspec%delc(i-1))/2.0
	  north(:,i)=temp
	end do

	return

end subroutine rel_centre_coords
 