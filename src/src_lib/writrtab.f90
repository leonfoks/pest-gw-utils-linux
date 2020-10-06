!     Last change:  JD   13 Jun 2002   11:53 am
subroutine write_real_table_file(ifail,outunit,outfile,intarray,realarray)

! -- Subroutine write_real_table_file writes a real array table file containing
!    the "active" cell values of a real array.

! -- Arguments are as follows:-
!       ifail:        returned as non-zero if error condition encountered
!       outunit:      unit number of output file
!       outfile:      name of output file
!       intarray:     window integer array
!       realarray:    data integer array


	use defn
	use inter
	implicit none

	integer, intent(out)                    :: ifail
	integer, intent(in)                     :: outunit
	character (len=*), intent(in)           :: outfile
	integer, dimension(:,:), intent(in)     :: intarray
	real, dimension(:,:), intent(in)        :: realarray

	integer irow,icol

	ifail=0
	do irow=1,size(intarray,2)
	  do icol=1,size(intarray,1)
	    if(intarray(icol,irow).eq.0) cycle
	    write(outunit,50,err=900) irow,icol,realarray(icol,irow)
50          format(1x,i5,3x,i5,3x,1pe14.6)
	  end do
	end do
	return

900     ifail=1
	write(amessage,910) trim(outfile)
910     format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes',endspace='no')
	return

end subroutine write_real_table_file
 