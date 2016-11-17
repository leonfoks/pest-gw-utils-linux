subroutine write_dxf_polyfin(iunit)

! -- Subroutine write_dxf_polyfin writes the footer to a polyline on a DXF
!    output file.

! -- Arguments are as follows:-
!       iunit:  unit number of DXF output file

! -- Revision history:-
!       June-November, 1995: version 1.

	integer, intent(in)     :: iunit

10      format(i3)
20      format(a)

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'SEQEND'
	write(iunit,10,err=100) 8
	write(iunit,20,err=100) '0'

	return

100     write(6,110)
110     format(' Cannot write DXF file: file inaccessible or disk full.',/)
	stop

end subroutine write_dxf_polyfin
 