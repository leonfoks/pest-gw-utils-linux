!     Last change:  JD   16 Dec 2000    6:31 pm
program rotdat

! -- Program ROTDAT rotates and displaces the coordinates found in a 
!    columnar data file about the coordinates of the upper left corner of a
!    finite difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) 			:: gridspec
	integer					:: ifail,inunit,outunit,iline,&
						   east_col,north_col,max_col,&
						   numcol,j
	real					:: rotation,delta
	double precision			:: easting,northing,rot_easting,&
						   rot_northing
	character (len=80)			:: aprompt,infile,outfile
	character (len=20)			:: anum,anum1,anum2


	write(amessage,5)
5	format(' Program ROTDAT rotates and displaces the coordinates found ',&
	'in a columnar data file about the coordinates of the upper left ',&
 	'corner of a finite difference grid.')
	call write_message(leadspace='yes',endspace='yes')

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
	   'a data file.')
	   call write_message(leadspace='yes')
	   go to 9900
	end if

30	write(6,*)
	aprompt=' Enter name of input data file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  call free_grid_mem(gridspec)
	  escset=0
	  go to 10
	end if

210	write(6,220,advance='no')
220	format(' In which column are the east coordinates? ')
	if(key_read(east_col).ne.0) go to 210
	if(escset.ne.0)then
	  escset=0
	  close(unit=inunit)
	  go to 30
	end if
	if(pos_test(east_col,'column number').ne.0) go to 210
	if(east_col.ge.NUM_WORD_DIM)then
	  call num2char(NUM_WORD_DIM,anum)
	  write(amessage,230) trim(anum)
230	  format(' Column number must be less than ',a,'  - try again.')
	  call write_message
	  go to 210
	end if
250	write(6,260,advance='no')
260	format(' In which column are the north coordinates? ')
	if(key_read(north_col).ne.0) go to 250
	if(escset.eq.1)then
	  escset=0
	  write(6,*)
	  go to 210
	end if
	if(pos_test(north_col,'column number').ne.0) go to 250
	if(east_col.eq.north_col)then
	  write(amessage,270)
270	  format(' Cannot be same as east coordinate column  - try again.')
	  call write_message
	  go to 250
	end if
	if(north_col.ge.NUM_WORD_DIM)then
	  call num2char(NUM_WORD_DIM,anum)
	  write(amessage,230) trim(anum)
	  call write_message
	  go to 250
	end if
	max_col=max(east_col,north_col)
	numcol=max_col

	write(6,*)
	aprompt=' Enter name for rotated output data file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
	  go to 250
	end if

	iline=0
	do
	  iline=iline+1
	  if(iline.eq.1000)then
	    write(6,280)
280	    format(/,' Working.....')
	  end if
	  read(inunit,'(a)',err=9000,end=500) cline 
	  call linesplit(ifail,max_col)
	  if(ifail.lt.0) then
	    write(outunit,*,err=9300)
	    cycle
	  end if
	  if(ifail.gt.0) then
	    call num2char(iline,anum)
	    write(amessage,300) trim(anum),trim(infile)
300	    format(' Insufficient data columns in line ',a,' of file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  easting=char2double(ifail,east_col)
	  if(ifail.ne.0)then
	    call num2char(east_col,anum1)
	    call num2char(iline,anum2)
	    write(amessage,320) trim(anum1),trim(anum2),trim(infile)
320	    format(' Cannot read easting from column ',a,' of line ',a,&
	    ' of input data file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  northing=char2double(ifail,north_col)
	  if(ifail.ne.0)then
	    call num2char(north_col,anum1)
	    call num2char(iline,anum2)
	    write(amessage,325) trim(anum1),trim(anum2),trim(infile)
325	    format(' Cannot read northing from column ',a,' of line ',a,&
	    ' of input data file ',a)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  easting=easting-gridspec%east_corner
	  northing=northing-gridspec%north_corner
	  rot_easting=easting*gridspec%cosang+northing*gridspec%sinang
	  rot_northing=northing*gridspec%cosang-easting*gridspec%sinang
	  call num2char(rot_easting,anum1,12)
	  call num2char(rot_northing,anum2,12)
	  do j=numcol,NUM_WORD_DIM
	    call linesplit(ifail,j)
	    if(ifail.ne.0) go to 350
	  end do	   
	  call num2char(iline,anum)
	  write(amessage,330) trim(anum),trim(infile)
330	  format(' Too many data columns in line ',a,' of file ',a)
	  call write_message(leadspace='yes')
	  go to 9900
350	  if(j.eq.numcol) then
	    do j=numcol-1,max_col+1,-1
	      call linesplit(ifail,j)
	      if(ifail.eq.0) then
	        numcol=j
	        go to 360
	      end if 
	    end do
	    numcol=max_col
360	    continue
	  else
	    numcol=j-1
	    call linesplit(ifail,numcol)
	  end if
	  do j=1,numcol-1
	    if(j.eq.east_col)then
	      write(outunit,370,advance='no',err=9300) trim(anum1)
370	      format(2x,a,',')
	    else if(j.eq.north_col)then
	      write(outunit,370,advance='no',err=9300) trim(anum2)
	    else
	      write(outunit,370,advance='no',err=9300) &
	      cline(left_word(j):right_word(j))
	    end if
	  end do
	  j=numcol
	  if(j.eq.east_col)then
	    write(outunit,380,advance='no',err=9300) trim(anum1)
380	    format(2x,a)
	  else if(j.eq.north_col)then
	    write(outunit,380,advance='no',err=9300) trim(anum2)
	  else
	    write(outunit,380,advance='no',err=9300) &
	    cline(left_word(j):right_word(j))
	  end if
	  write(outunit,*,err=9300)
	end do

500	write(6,*)
	call num2char(iline-1,anum)
	write(amessage,115) trim(anum),trim(infile)
115	format('  - ',a,' lines read from input data file ',a)
	call write_message
	write(amessage,120) trim(outfile)
120	format('  - rotated output data file ',a,' written ok.')
	call write_message

	go to 9900

9000	call num2char(iline,anum)
	write(amessage,9010) trim(anum),trim(infile)
9010	format(' Error reading line ',a,' of input data file ',a)
	call write_message(leadspace='yes')
	go to 9900
9300	write(amessage,9310) trim(outfile)
9310	format(' Error writing to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 9900

9900	call close_files
	call free_grid_mem(gridspec)

end program rotdat

 