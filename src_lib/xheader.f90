subroutine write_dxf_header(iunit,xlo,xhi,ylo,yhi)

! -- Subroutine write_dxf_header writes the header to a DXF file.

! -- Arguments are as follows:-
!       iunit:     the unit number representing the file (input only)
!       xlo, ylo:  minimum x and y coordinates of drawing
!       xhi, yhi:  maximum x and y coordinates of drawing

! -- Revision history:-
!       June-November, 1995: version 1.

	integer, intent(in)             :: iunit
	double precision, intent(in)    :: xlo,xhi,ylo,yhi

10      format(i3)
20      format(a)
30      format(1pg20.12)

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'SECTION'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'HEADER'
	write(iunit,10,err=100) 9
	write(iunit,20,err=100) '$EXTMIN'
	write(iunit,10,err=100) 10
	write(iunit,30,err=100) xlo
	write(iunit,10,err=100) 20
	write(iunit,30,err=100) ylo
	write(iunit,10,err=100) 9
	write(iunit,20,err=100) '$EXTMAX'
	write(iunit,10,err=100) 10
	write(iunit,30,err=100) xhi
	write(iunit,10,err=100) 20
	write(iunit,30,err=100) yhi
	write(iunit,10,err=100) 9
	write(iunit,20,err=100) '$LIMMIN'
	write(iunit,10,err=100) 10
	write(iunit,30,err=100) xlo
	write(iunit,10,err=100) 20
	write(iunit,30,err=100) ylo
	write(iunit,10,err=100) 9
	write(iunit,20,err=100) '$LIMMAX'
	write(iunit,10,err=100) 10
	write(iunit,30,err=100) xhi
	write(iunit,10,err=100) 20
	write(iunit,30,err=100) yhi
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'ENDSEC'

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'SECTION'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'TABLES'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'TABLE'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'LTYPE'
	write(iunit,10,err=100) 70
	write(iunit,20,err=100) '2'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'LTYPE'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'CONTINUOUS'
	write(iunit,10,err=100) 70
	write(iunit,20,err=100) '0'
	write(iunit,10,err=100) 3
	write(iunit,20,err=100) 'Solid Line'
	write(iunit,10,err=100) 72
	write(iunit,20,err=100) '65'
	write(iunit,10,err=100) 73
	write(iunit,20,err=100) '0'
	write(iunit,10,err=100) 40
	write(iunit,20,err=100) '0.0'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'LTYPE'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'DASHED'
	write(iunit,10,err=100) 70
	write(iunit,20,err=100) '0'
	write(iunit,10,err=100) 3
	write(iunit,20,err=100) '_ _ _ _ _'
	write(iunit,10,err=100) 72
	write(iunit,20,err=100) '65'
	write(iunit,10,err=100) 73
	write(iunit,20,err=100) '2'
	write(iunit,10,err=100) 40
	write(iunit,20,err=100) '0.75'
	write(iunit,10,err=100) 49
	write(iunit,20,err=100) '0.5'
	write(iunit,10,err=100) 49
	write(iunit,20,err=100) '-0.25'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'ENDTAB'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'ENDSEC'

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'SECTION'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'BLOCKS'
	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'ENDSEC'

	write(iunit,10,err=100) 0
	write(iunit,20,err=100) 'SECTION'
	write(iunit,10,err=100) 2
	write(iunit,20,err=100) 'ENTITIES'
	
	return

100     write(6,110)
110     format(' Cannot write DXF file: file inaccessible or disk full.',/)
	stop

end subroutine write_dxf_header
 