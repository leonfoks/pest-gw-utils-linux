subroutine free_bore_mem

! -- Subroutine free_bore_mem frees memory and deallocates pointers pertinent
!    to the storage of bore data.

! -- Revision history:-
!      June-November, 1995: version 1.

	use defn
	use inter

	integer                 :: ierr

	deallocate(bore_coord_id,bore_coord_east,bore_coord_north, &
	bore_list_id,bore_coord_layer,bore_diff_id,stat=ierr)
	nullify (bore_coord_id,bore_coord_east,bore_coord_north,bore_list_id,&
	bore_coord_layer,bore_diff_id)

end subroutine free_bore_mem
