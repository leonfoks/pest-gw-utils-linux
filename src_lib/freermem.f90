!     Last change:  JD    8 Feb 2001    1:19 pm

subroutine free_replace_mem

! -- Subroutine free_replace_mem frees memory and deallocates pointers pertinent
!    to MODFLOW 2000 replacement point data.

	use defn
	use inter

	integer                 :: ierr,irep


!        if(numparamreplace.ne.0)then
!          if(associated(replace)) then
!            do irep=1,numparamreplace
!              if(associated(replace(irep)%maxiarray))then
!                deallocate(replace(irep)%maxiarray,stat=ierr)
!                nullify   (replace(irep)%maxiarray)
!              end if
!              if(associated(replace(irep)%miniarray))then
!                deallocate(replace(irep)%miniarray,stat=ierr)
!                nullify   (replace(irep)%miniarray)
!              end if
!            end do
!          end if
!        end if

        if(associated(removezone)) then
          deallocate(removezone)
 !         nullify(removezone)
        end if
        if(associated(removemult)) then
          deallocate(removemult)
 !         nullify(removemult)
        end if
        if(associated(replace)) then
          deallocate(replace)
 !         nullify(replace)
        end if

end subroutine free_replace_mem

