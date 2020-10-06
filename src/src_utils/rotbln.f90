!     Last change:  JD   16 Dec 2000    6:39 pm
program rotbln

! -- Program ROTBLN rotates and displaces the coordinates found in a SURFER
!    blanking file about the coordinates of the upper left corner of a finite-
!    difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) 			:: gridspec
	integer					:: ifail,inunit,outunit,iline,&
						   iseg,numpts,i,iwork
	real					:: rotation,delta
	double precision			:: easting,northing,rot_easting,&
						   rot_northing
	character (len=80)			:: aprompt,infile,outfile
	character (len=20)			:: anum


	write(amessage,5)
5	format(' Program ROTBLN rotates and displaces the coordinates found ',&
	'in a SURFER blanking file about the coordinates of the upper left ',&
 	'corner of a finite difference grid.')

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
	   'a blanking file.')
	   call write_message(leadspace='yes')
	   go to 9900
	end if

30	write(6,*)
	aprompt=' Enter name of SURFER blanking input file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  call free_grid_mem(gridspec)
	  escset=0
	  go to 10
	end if

	aprompt=' Enter name for rotated SURFER blanking output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  close(unit=inunit)
	  escset=0
	  go to 30
	end if

	iseg=0
	iwork=0
	iline=0
40	iline=iline+1
	read(inunit,'(a)',err=9000,end=500) cline
	call linesplit(ifail,1)
	if(ifail.ne.0) go to 40
	numpts=char2int(ifail,1)
	if(ifail.ne.0) go to 9400
	iseg=iseg+1
	write(outunit,'(a)',err=9300) trim(cline)
	do i=1,numpts
50	  iline=iline+1
	  if((iwork.eq.0).and.(iline.gt.2000))then
	    write(6,60)
60	    format(/,' Working.....')
	    iwork=1
	  end if
	  read(inunit,'(a)',err=9000,end=9600) cline
	  call linesplit(ifail,2)
	  if(ifail.lt.0) go to 50
	  if(ifail.gt.0) then
	    call num2char(iline,anum)
	    write(amessage,70) trim(anum),trim(infile)
70	    format(' Two entries expected on line ',a,' of file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  easting=char2double(ifail,1)
	  if(ifail.ne.0) then
	    call num2char(iline,anum)
	    write(amessage,90) trim(anum),trim(infile)
90	    format(' Error reading east coordinate from line ',a,' of file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  northing=char2double(ifail,2)
	  if(ifail.ne.0) then
	    call num2char(iline,anum)
	    write(amessage,110) trim(anum),trim(infile)
110	    format(' Error reading north coordinate from line ',a,' of file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  easting=easting-gridspec%east_corner
	  northing=northing-gridspec%north_corner
	  rot_easting=easting*gridspec%cosang+northing*gridspec%sinang
	  rot_northing=northing*gridspec%cosang-easting*gridspec%sinang
	  write(outunit,130,err=9300) rot_easting,rot_northing
130 	  format(1x,1pg20.12,1x,1pg20.12)
	end do
	go to 40

500	if(iseg.eq.0)then
	  write(amessage,510) trim(infile)
510	  format(' No line segments found in SURFER blanking input file ',a)
	  call write_message(leadspace='yes')
	  go to 9900
	end if
	write(6,*)
	call num2char(iline-1,anum)
	write(amessage,115) trim(anum),trim(infile)
115	format('  - ',a,' lines read from SURFER blanking input file ',a)
	call write_message
	write(amessage,120) trim(outfile)
120	format('  - rotated SURFER blanking output file ',a,' written ok.')
	call write_message

	go to 9900

9000	call num2char(iline,anum)
	write(amessage,9010) trim(anum),trim(infile)
9010	format(' Error reading line ',a,' of SURFER blanking file ',a)
	call write_message(leadspace='yes')
	go to 9900
9300	write(amessage,9310) trim(outfile)
9310	format(' Error writing to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 9900
9400	call num2char(iline,anum)
	write(amessage,9410) trim(anum),trim(infile)
9410	format(' Integer expected at line ',a,' of SURFER blanking file ',a)
	call write_message(leadspace='yes')
	go to 9900
9600	write(amessage,9510) trim(infile)
9510	format(' Unexpected end to SURFER blanking input file ',a)
	call write_message(leadspace='yes')
	go to 9900

9900	call close_files
	call free_grid_mem(gridspec)

end program rotbln

 