subroutine write_dxf_polyhead(iunit)

! -- Subroutine write_dxf_polyhead writes the header to a polyline on a DXF
!    file.

! -- Arguments are as follows:-
!       iunit:    unit number of output file

! Revision history:-
!       June-November, 1995: version 1.

	integer, intent(in)     :: iunit

10      format(i3)
20      format(a)

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'POLYLINE'
	write(iunit,10,err=100) 8
	write(iunit,20,err=100) '0'
	write(iunit,10,err=100) 6
	write(iunit,20,err=100) 'CONTINUOUS'
	write(iunit,10,err=100) 62
	write(iunit,20,err=100) '7'
	write(iunit,10,err=100) 66
	write(iunit,20,err=100) '1'

	return

100     write(6,110)
110     format(' Cannot write DXF file: file inaccessible or disk full.',/)
	stop

end subroutine write_dxf_polyhead
 