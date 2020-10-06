subroutine free_grid_mem(gridspec)

! -- Subroutine free_grid_mem frees memory and deallocates pointers pertinent
!    to the storage of grid data.

! -- Arguments are as follows:-
!       gridspec:  defined variable holding grid specifications

! -- Revision history:-
!       June-November, 1995: version1.

	use defn
	use inter

	type (modelgrid)        :: gridspec
	integer                 ::ierr

	deallocate(gridspec%delr,gridspec%delc,stat=ierr)
	nullify (gridspec%delr,gridspec%delc)

end subroutine free_grid_mem
 