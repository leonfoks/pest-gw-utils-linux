!     Last change:  JD   21 Dec 2000    5:30 pm

subroutine free_point_mem

! -- Subroutine free_point_mem frees memory and deallocates pointers pertinent
!    to the storage of pilot point data.

	use defn
	use inter

	integer                 :: ierr

	deallocate(pilot_point_id,pilot_point_east,pilot_point_north, &
	pilot_point_zone,pilot_point_val,stat=ierr)
	nullify   (pilot_point_id,pilot_point_east,pilot_point_north, &
        pilot_point_zone,pilot_point_val)
        if(associated(pilot_point_val1))then
          deallocate(pilot_point_val1,stat=ierr)
          nullify   (pilot_point_val1)
        end if

end subroutine free_point_mem

