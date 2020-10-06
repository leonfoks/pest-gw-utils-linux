subroutine write_integer_mid_file(ifail,outunit,outfile,intarraywin,intarraydat)

! -- Subroutine write_integer_mid_file writes a MAPINFO "MID" file containing
!    the "active" cell values of an integer array.

! -- Arguments are as follows:-
!       ifail:        returned as non-zero if error condition encountered
!       outunit:      unit number of output file
!       outfile:      name of output file
!       intarraywin:  window integer array
!       intarraydat:  data integer array

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter
	implicit none

	integer, intent(out)                    :: ifail
	integer, intent(in)                     :: outunit
	character (len=*), intent(in)           :: outfile
	integer, dimension(:,:), intent(in)     :: intarraywin,intarraydat

	integer irow,icol

	ifail=0
	do irow=1,size(intarraywin,2)
	  do icol=1,size(intarraywin,1)
	    if(intarraywin(icol,irow).eq.0) cycle
	    write(outunit,50,err=900) irow,icol,intarraydat(icol,irow)
50          format(1x,i5,',',1x,i5,',',1x,i10)
	  end do
	end do
	return

900     ifail=1
	write(amessage,910) trim(outfile)
910     format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes',endspace='no')
	return

end subroutine write_integer_mid_file
 