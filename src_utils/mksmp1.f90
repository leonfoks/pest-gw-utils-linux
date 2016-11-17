!     Last change:  JD   17 Dec 2000   10:08 pm
program mksmp1

! -- Program MKSMP1 builds a bore sample file from a Queensland Department
!    of Natural Resources borehole water elevation file.

	use defn
	use inter

	implicit none

	integer	:: ifail,inunit,iline,cols,outunit,dd,mm,yy,hhh,mmm,sss,idate, &
		   idatespec,iheader
	double precision :: value
	character (len=80) :: aprompt,infile,outfile
	character (len=20) :: aword,atime
	character (len=15) :: aline
	character (len=4)  :: atemp


3	format(/,' Errors in water elevation file ',a,' ------->')

	write(amessage,5)
5       format(' Program MKSMP1 builds a bore sample file from a Queensland ',&
	'Department of Natural Resources borehole water elevation file.')
	call write_message(leadspace='yes',endspace='yes')

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          call write_message(endspace='yes')
          go to 9900
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          call write_message(endspace='yes')
          go to 9900
        endif
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
          'settings.fig')
          call write_message(endspace='yes')
          go to 9900
        end if
	idatespec=datespec
	datespec=1

10	aprompt=' Enter name of water elevation file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

15	aprompt=' Enter name for bore sample output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=inunit)
	  write(6,*)
	  go to 10
	end if	

	write(6,*)
20	write(6,30,advance='no')
30	format(' Enter notional time for all water level readings ',&
	'[hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 20
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) then
	  close(unit=outunit)
	  write(6,*)
	  go to 15
	end if
	call char2time(ifail,atime,hhh,mmm,sss)
	if(ifail.ne.0) then
	  write(6,40)
40	  format(' Illegal time  - try again.')
	  go to 20
	end if
	call time2char(ifail,hhh,mmm,sss,atime)

	imessage=0
	iline=0
	read_a_line: do
	  if(imessage.gt.40) then
	    write(6,*)
	    go to 1025
	  end if
	  iline=iline+1
	  if((iline.eq.2000).and.(imessage.eq.0)) then
	    write(6,50)
50	    format(/,' Working.....')
	  end if
	  cols=5
	  read(inunit,'(a)',err=9000,end=1000) cline
	  call linesplit(ifail,5)
	  if(ifail.eq.-1) cycle read_a_line
	  if(ifail.ne.0) then
	    cols=4
	    call linesplit(ifail,4)
	    if(ifail.ne.0) then
	      if(imessage.eq.0) write(6,3) trim(infile)
	      call num2char(iline,aline)
	      write(amessage,70) trim(aline)
70	      format('   Line ',a,': insufficient items (minimum of 4 ',&
	      'expected).')
	      call write_message(increment=1)
	      cycle read_a_line
	    end if
	  end if
	  if(right_word(1)-left_word(1).gt.8)then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,80) trim(aline)
80	    format('   Line ',a,': bore number greater than 9 characters long.')
	    call write_message(increment=1)
	  end if
	  if(right_word(2).ne.left_word(2))then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,90) trim(aline)
90	    format('   Line ',a,': second data entry (ie. pipe id) ',&
	    'should be only one character long.')
	    call write_message(increment=1)
	  end if
	  call char2date(ifail,cline(left_word(3):right_word(3)),dd,mm,yy)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,110) trim(aline)
110	    format('   Line ',a,': cannot read date.')
	    call write_message(increment=1)
	  end if
	  value=char2double(ifail,4)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,120) trim(aline)
120	    format('   Line ',a,': cannot read sample value.')
	    call write_message(increment=1)
	  end if

	  if(cols.eq.5)then
	    atemp=cline(left_word(5):left_word(5)+3)
	    if((index(atemp,'P').eq.0).and.(index(atemp,'p').eq.0).and.&
	       (index(atemp,'D').eq.0).and.(index(atemp,'d').eq.0)) then
	       if(imessage.eq.0) write(6,3) trim(infile)
	       call num2char(iline,aline)
	       write(amessage,130) trim(aline)
130	       format('   Line ',a,': unrecognized character(s) in optional ',&
	       'fifth column.')
	       call write_message(increment=1)
	    end if
	  end if

	  if(imessage.ne.0) cycle read_a_line
	  aword=cline(left_word(1):right_word(1))// &
	    cline(left_word(2):right_word(2))

	  datespec=idatespec
	  if(cols.eq.4)then
	    if(datespec.eq.1) then
	      write(outunit,140,err=9200) trim(aword),&
	      cline(left_word(3):right_word(3)),trim(atime),&
	      cline(left_word(4):right_word(4))
140	      format(1x,a10,t16,a,t29,a,t39,a)
	    else
	      write(outunit,145,err=9200) trim(aword),mm,dd,yy, &
              trim(atime),cline(left_word(4):right_word(4))
145	      format(1x,a10,t16,i2.2,'/',i2.2,'/',i4.4,t29,a,t39,a)
	    end if
	  else
	    if(datespec.eq.1) then
	      write(outunit,150,err=9200) trim(aword),&
	      cline(left_word(3):right_word(3)),trim(atime),&
	      cline(left_word(4):right_word(4)),'x'
150	      format(1x,a10,t16,a,t29,a,t39,a,2x,a)
	    else
	      write(outunit,155,err=9200) trim(aword),mm,dd,yy, &
              trim(atime),cline(left_word(4):right_word(4)),'x'
155	      format(1x,a10,t16,i2.2,'/',i2.2,'/',i4.4,t29,a,t39,a,2x,a)
	    end if
	  end if
	  datespec=1

	end do read_a_line  

1000	write(6,*)
	call num2char (iline-1,aline)
	write(amessage,1010) trim(aline), trim(infile)
1010	format('  - ',a,' lines read from water elevation file ',a)
	call write_message

1025    if(imessage.eq.0) then
          write(amessage,1020) trim(outfile)
1020      format('  - bore sample file ',a,' written ok.')
          call write_message
          write(amessage,1030) trim(outfile)
1030      format(' It is adviseable to check file ',a,' with program SMPCHEK ',&
          'before using it, as program MKSMP1 undertakes only limited error ',&
          'checking in translating a water elevation file to a bore pumping file.')
          call write_message(leadspace='yes',endspace='yes')
          go to 9900
        else
          write(amessage,1050) trim(outfile),trim(infile)
1050      format('  - bore sample file ',a,' only partially written as ',&
          'errors were encountered in water elevation file ',a)
          call write_message(endspace='yes')
        end if
        go to 9900
 
9000	call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(infile)
9010	format(' Unable to read line ',a,' of water elevation file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 9900
9200	write(amessage,9210) trim(outfile)
9210	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9900

9900    call close_files

end program mksmp1

 