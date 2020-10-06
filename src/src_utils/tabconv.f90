!     Last change:  JD   16 Dec 2000    2:14 pm
program tabconv

! -- Program TABCONV translates from row/column number to cell number format
!    and from cell number to row/column number format for real and integer 
!    array table files.


	use defn
	use inter

	implicit none

	integer	:: ichoice, itemp, inunit, outunit, ifail,ncol,nrow,&
		   iline,icellnum,maxcell,irow,icol
	character (len=80)			:: infile,outfile,aprompt
	type (modelgrid) 			:: gridspec
	character (len=15)			:: aline


	imessage=0
	write(amessage,5)
5	format(' Program TABCONV translates from row/column number ',&
	'to cell number format and from cell number to row/column number ',&
        'format for real and integer array table files.')
	call write_message(leadspace='yes',endspace='yes')

	call readfig(gridspec%specfile)
10	call spec_open(ifail,gridspec)
	if((ifail.ne.0).or.(escset.eq.1)) go to 9950
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	maxcell=ncol*nrow
	write(6,*)

50	write(6,70)
70	format(' For cell identifier translation:')
	write(6,80)
80	format('   if from row/column number to cell number  -  enter 1')
	write(6,90)
90	format('   if from cell number to row/column number  -  enter 2')
100	write(6,110,advance='no')
110	format(' Enter your choice: ')
	itemp=key_read(ichoice)
	if(itemp.ne.0) go to 100
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 10
	end if
	if((ichoice.ne.1).and.(ichoice.ne.2)) go to 100
	write(6,*)

180	aprompt=' Enter name of array table input file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9950
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 50
	end if
	aprompt=' Enter name for array table output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9950
	if(escset.eq.1) then
	  escset=0
	  close(unit=inunit)
	  write(6,*)
	  go to 180
 	end if
	write(6,*)

	write(initial_message,320) trim(infile)
320	format(' Errors in array table input file ',a,':-')
	iline=0
	imessage=0
330	iline=iline+1
	if((iline.eq.1000).and.(imessage.eq.0)) write(6,350)
350	format(' Working.....',/)
	read(inunit,'(a)',err=9000,end=1000) cline
	call linesplit(ifail,4)
	if(ifail.lt.0) go to 330
	if(ifail.eq.0) then
	  if(imessage.eq.0) call write_initial_message
	  call num2char(iline,aline)	    
	  write(amessage,370) trim(aline)
370	  format('   Line ',a,': too many entries.')
	  call write_message(increment=1)
	  if(imessage.gt.25) go to 9900
	  go to 330
	end if
	call linesplit(ifail,3)
	if(ichoice.eq.1) then
	  if(ifail.gt.0) then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)	    
	    write(amessage,380) trim(aline)
380	    format('   Line ',a,': insufficient entries.')
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	    go to 330
	  end if
	else
	  if(ifail.eq.0)then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)	    
	    write(amessage,370) trim(aline)
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	    go to 330
	  end if
	  call linesplit(ifail,2)
	  if(ifail.gt.0) then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)	    
	    write(amessage,380) trim(aline)
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	    go to 330
	  end if
	end if

	if(ichoice.eq.2) then
	  icellnum=char2int(ifail,1)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)
	    write(amessage,400) trim(aline)
400	    format('   Line ',a,': cannot read cell number.')
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	  else
	    if(icellnum.gt.maxcell) then
	      if(imessage.eq.0) call write_initial_message
	      call num2char(iline,aline)
	      write(amessage,420) trim(aline)
420	      format('   Line ',a,': cell number exceeds maximum number of ',&
              'cells in grid.')
	      call write_message(increment=1)
	      if(imessage.gt.25) go to 9900
	    else if (icellnum.lt.1) then
	      if(imessage.eq.0) call write_initial_message
	      call num2char(iline,aline)
	      write(amessage,440) trim(aline)
440	      format('   Line ',a,': cell number cannot be less than 1.')
	      call write_message(increment=1)
	      if(imessage.gt.25) go to 9900
	    end if
	    if(imessage.eq.0) then
	      call cell2rc(icellnum,irow,icol,gridspec)
	      write(outunit,460) irow,icol,cline(left_word(2):right_word(2))
460	      format(1x,i6,1x,i6,2x,a)
	    end if
	  end if
	else
	  irow=char2int(ifail,1)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)
	    write(amessage,480) trim(aline)
480	    format('   Line ',a,': cannot read row number.')
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	  else
	    if((irow.lt.1).or.(irow.gt.nrow)) then
	      if(imessage.eq.0) call write_initial_message
	      call num2char(iline,aline)
	      write(amessage,500) trim(aline)
500	      format('   Line ',a,': row number out of range.')
	      call write_message(increment=1)
	      if(imessage.gt.25) go to 9900
	    end if
	  end if
	  icol=char2int(ifail,2)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) call write_initial_message
	    call num2char(iline,aline)
	    write(amessage,520) trim(aline)
520	    format('   Line ',a,': cannot read column number.')
	    call write_message(increment=1)
	    if(imessage.gt.25) go to 9900
	  else
	    if((icol.lt.1).or.(icol.gt.ncol)) then
	      if(imessage.eq.0) call write_initial_message
	      call num2char(iline,aline)
	      write(amessage,540) trim(aline)
540	      format('   Line ',a,': column number out of range.')
	      call write_message(increment=1)
	      if(imessage.gt.25) go to 9900
	    end if
	  end if
	  if(imessage.eq.0) then
	    call rc2cell(icellnum,irow,icol,gridspec)
	    write(outunit,560) icellnum,cline(left_word(3):right_word(3))
560	    format(1x,i6,2x,a)
	  end if
	end if
	go to 330

1000	if(imessage.gt.0) go to 9900
	call num2char(iline-1,aline)
	write(amessage,1040) trim(aline),trim(infile)
1040	format('  - ',a,' lines read from file ',a)
	call write_message
	write(amessage,1060) trim(outfile)
1060	format('  - file ',a,' written ok.')
	call write_message
	go to 9950

9000	call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(infile)
9010	format(' Cannot read line ',a,' of file ',a)
	call write_message(leadspace='yes', increment=1)
9900	if(imessage.ne.0) then
	  write(amessage,9910) trim(outfile),trim(infile)
9910	  format(' TABCONV output file ',a,' only partially written as ',&
	  'errors were encountered in input table file ',a)
 	  call write_message(leadspace='yes')
	end if
9950	write(6,*)
	call close_files

end program tabconv

 