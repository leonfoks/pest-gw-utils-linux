!     Last change:  JD   16 Dec 2000    6:42 pm
program rotdxf

! -- Program ROTDXF rotates and displaces the coordinates found in a DXF file
!    about the coordinates of the upper left corner of a finite difference
!    grid.

	use defn
	use inter

	implicit none

	type (modelgrid) 			:: gridspec
	integer					:: ifail,inunit,outunit,iline,&
						   icoord,iwork,i11
	real					:: rotation,delta
	double precision			:: easting,northing,rot_easting,&
						   rot_northing
	character (len=80)			:: aprompt,infile,outfile
	character (len=20)			:: anum
	character (len=300)			:: cline1,cline2,cline3


	write(amessage,5)
5	format(' Program ROTDXF rotates and displaces the coordinates found ',&
	'in a DXF file about the coordinates of the upper left corner of a ',&
	'finite difference grid.')

	call readfig(gridspec%specfile)
10	call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	rotation=gridspec%rotation
	delta=spacing(rotation)
	if((gridspec%rotation-delta.le.0.0d0).and. &
	   (gridspec%rotation+delta.ge.0.0d0)) then
	   write(amessage,20) trim(gridspec%specfile)
20	   format(' The grid defined in grid specification file ',a, &
	   ' is oriented east-west. There is no need to rotate and displace ',&
	   'a DXF file.')
	   call write_message(leadspace='yes')
	   go to 9900
	end if

30	write(6,*)
	aprompt=' Enter name of DXF input file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  call free_grid_mem(gridspec)
	  escset=0
	  go to 10
	end if

	aprompt=' Enter name for rotated DXF output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  close(unit=inunit)
	  escset=0
	  go to 30
	end if

	iwork=0
	iline=0
	icoord=0
	do
	  iline=iline+1
	  if((iwork.eq.0).and.(iline.gt.2000))then
	    write(6,40)
40	    format(/,' Working .....',/)
	    iwork=1
	  endif
	  read(inunit,'(a)',err=9000,end=500) cline
	  call linesplit(ifail,1)
	  if(ifail.lt.0) then
            write(outunit,'(a)') '  '
            cycle
	  end if
	  if(ifail.ne.0) go to 9000
	  if((cline(left_word(1):right_word(1)).eq.'10').or.&
	     (cline(left_word(1):right_word(1)).eq.'11')) then
	    i11=0
	    if(cline(left_word(1):right_word(1)).eq.'11') i11=1
	    cline1=cline
	    icoord=icoord+1
	    iline=iline+1
	    read(inunit,'(a)',err=9000,end=9500) cline
	    call linesplit(ifail,1)
	    if(ifail.ne.0) go to 50
	    easting=char2double(ifail,1)
	    if(ifail.ne.0) go to 50
	    go to 60
50	    write(outunit,'(a)',err=9300) trim(cline1)
	    write(outunit,'(a)',err=9300) trim(cline)
	    go to 112
60	    cline3=cline
	    iline=iline+1
	    read(inunit,'(a)',err=9000,end=9500) cline
	    call linesplit(ifail,1)
	    if(ifail.ne.0) go to 70
	    if(i11.eq.0)then
	      if(cline(left_word(1):right_word(1)).ne.'20')go to 70
	    else
	      if(cline(left_word(1):right_word(1)).ne.'21')go to 70
	    end if
	    cline2=cline
	    go to 80
70	    write(outunit,'(a)',err=9300) trim(cline1)
	    write(outunit,'(a)',err=9300) trim(cline3)
	    write(outunit,'(a)',err=9300) trim(cline)
	    go to 112
80	    iline=iline+1
	    read(inunit,'(a)',err=9000,end=9500) cline
	    call linesplit(ifail,1)
	    if(ifail.ne.0) go to 90
	    northing=char2double(ifail,1)
	    if(ifail.ne.0) go to 90
	    go to 95
90	    write(outunit,'(a)',err=9300) trim(cline1)
	    write(outunit,'(a)',err=9300) trim(cline3)
	    write(outunit,'(a)',err=9300) trim(cline2)
	    write(outunit,'(a)',err=9300) trim(cline)
	    go to 112
95	    easting=easting-gridspec%east_corner
	    northing=northing-gridspec%north_corner
	    rot_easting=easting*gridspec%cosang+northing*gridspec%sinang
	    rot_northing=northing*gridspec%cosang-easting*gridspec%sinang
	    write(outunit,'(a)',err=9300) trim(cline1)
	    write(outunit,100,err=9300) rot_easting
100	    format(1pg20.12)
	    write(outunit,'(a)',err=9300) trim(cline2)
	    write(outunit,100,err=9300) rot_northing
	  else
	    write(outunit,'(a)',err=9300) trim(cline)
	  end if
112	  continue
	end do

500	write(6,*)
	call num2char(iline-1,anum)
	write(amessage,115) trim(anum),trim(infile)
115	format('  - ',a,' lines read from DXF input file ',a)
	call write_message
	write(amessage,120) trim(outfile)
120	format('  - rotated DXF output file ',a,' written ok.')
	call write_message
	if(icoord.eq.0)then
	   write(amessage,140) trim(infile)
140	   format(' Warning: no eastings and northings found in DXF input file ',a)
	   call write_message
	end if
	    
	go to 9900

9000	call num2char(iline,anum)
	write(amessage,9010) trim(anum),trim(infile)
9010	format(' Error reading line ',a,' of DXF file ',a)
	call write_message(leadspace='yes')
	go to 9900
9300	write(amessage,9310) trim(outfile)
9310	format(' Error writing to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 9900
9500	write(amessage,9510)trim(infile)
9510	format(' Unexpected end encountered to DXF file ',a)
	call write_message(leadspace='yes')
	go to 9900

9900	call close_files
	call free_grid_mem(gridspec)

end program rotdxf

 