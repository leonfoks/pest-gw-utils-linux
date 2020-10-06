subroutine write_dxf_vertex(iunit,x,y)

! -- Subroutine write_dxf_vertex writes a vertex within a polyline to a DXF
!    file.

! -- Subroutine arguments are as follows:-
!       iunit:   unit number of output file
!       x,y:     x and y coordinates of vertex

! -- Revision history:-
!       June-November, 1995:  version 1.

	integer, intent(in)             :: iunit
	double precision, intent(in)    :: x,y

10      format(i3)
20      format(a)
30      format(1pg20.12)

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'VERTEX'
	write(iunit,10,err=100) 8
	write(iunit,20,err=100) '0'
	write(iunit,10,err=100) 10
	write(iunit,30,err=100) x
	write(iunit,10,err=100) 20
	write(iunit,30,err=100) y
	write(iunit,10,err=100) 30
	write(iunit,30,err=100) 0.0d0

	return

100     write(6,110)
110     format(' Cannot write DXF file: file inaccessible or disk full.',/)
	stop

end subroutine write_dxf_vertex
 