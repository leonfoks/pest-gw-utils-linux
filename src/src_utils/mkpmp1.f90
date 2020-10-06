!     Last change:  JD   17 Dec 2000    9:55 pm
program mkpmp1

! --Program MKPMP1 builds a bore pumping file from a Queensland Department
!   of Natural Resources metered use file.


	use defn
	use inter

	implicit none

	integer	:: ifail,inunit,iline,outunit,dd1,mm1,yy1,dd2,mm2,yy2, &
	hhh,mmm,sss,idate,iheader
	double precision :: pumpage
	character (len=80) :: aprompt,infile,outfile
	character (len=20) :: atime
	character(len=10)  :: abore,aline
	character(len=11)  :: anum


3	format(/,' Errors in metered use file ',a,' ------->')

	write(amessage,5)
5       format(' Program MKPMP1 builds a bore pumping file from a Queensland ',&
	'Department of Natural Resources metered use file.')
	call write_message(leadspace='yes',endspace='yes')

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
	  call write_message(endspace='yes')
	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message(endspace='yes')
	  go to 9900
	endif
	if((idate.ne.0).or.(datespec.eq.0)) then
	  write(amessage,9)
9	  format(' Cannot read date format from settings file ', &
	  'settings.fig')
	  call write_message(endspace='yes')
	  go to 9900
	end if

10	aprompt=' Enter name of metered use file: '
	call open_input_file(ifail,aprompt,infile,inunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

15	aprompt=' Enter name for bore pumping output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=inunit)
	  write(6,*)
	  go to 10
	end if	

        write(6,*)
20      write(6,30,advance='no')
30      format(' Enter notional time for all pumpage readings ',&
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
40        format(' Illegal time  - try again.')
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
	  read(inunit,'(a)',err=9000,end=1000) cline
	  call linesplit(ifail,8)
	  if(ifail.eq.-1) cycle read_a_line
	  if(ifail.gt.0) then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,70) trim(aline)
70	    format('   Line ',a,': insufficient entries.')
	    call write_message(increment=1)
	    cycle read_a_line
	  end if
	  if(right_word(1)-left_word(1).gt.9) then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,80) trim(aline)
80	    format('   Line ',a,': bore identifier greater than 10 ',&
	    'characters in length.')
	    call write_message(increment=1)
	  else
	    abore=cline(left_word(1):right_word(1))
	  end if
	  dd2=char2int(ifail,2)
	  if(ifail.eq.0) mm2=char2int(ifail,3)
	  if(ifail.eq.0) yy2=char2int(ifail,4)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,90) trim(aline)
90	    format('   Line ',a,': cannot read first date.')
	    call write_message(increment=1)
	    cycle read_a_line
	  end if
	  dd1=char2int(ifail,5)
	  if(ifail.eq.0)mm1=char2int(ifail,6)
	  if(ifail.eq.0)yy1=char2int(ifail,7)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,110) trim(aline)
110	    format('   Line ',a,': cannot read second date.')
	    call write_message(increment=1)
	    cycle read_a_line
	  end if
	  pumpage=char2double(ifail,8)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,3) trim(infile)
	    call num2char(iline,aline)
	    write(amessage,130) trim(aline)
130	    format('   Line ',a,': cannot read pumping value.')
	    call write_message(increment=1)
	    cycle read_a_line
	  end if
	  if(imessage.ne.0) cycle
	  call num2char(pumpage,anum)
	  if(datespec.eq.1) then
	    write(outunit,150,err=9200) trim(abore),dd1,mm1,yy1, &
	    trim(atime),dd2,mm2,yy2,trim(atime),trim(anum)
150	    format(1x,a10,t16,i2.2,'/',i2.2,'/',i4,t27,a8,t38, &
            i2.2,'/',i2.2,'/',i4,t49,a,t59,a)
	  else
	    write(outunit,150,err=9200) trim(abore),mm1,dd1,yy1, &
	    trim(atime),mm2,dd2,yy2,trim(atime),trim(anum)
	  end if
	end do read_a_line

1000	write(6,*)
	call num2char (iline-1,anum)
	write(amessage,1010) trim(anum), trim(infile)
1010	format('  - ',a,' lines read from metered use file ',a)
	call write_message
1025	if(imessage.eq.0) then
	  write(amessage,1020) trim(outfile)
1020	  format('  - bore pumping file ',a,' written ok.')
	  call write_message
	  write(amessage,1030) trim(outfile)
1030	  format(' It is adviseable to check file ',a,' with program PMPCHEK ',&
	  'before using it, as program MKPMP1 undertakes only limited error ',&
	  'checking in translating a metered use file to a bore pumping file.')
	  call write_message(leadspace='yes',endspace='yes') 
	  go to 9900
	else
	  write(amessage,1050) trim(outfile),trim(infile)
1050	  format('  - bore pumping file ',a,' only partially written as ',&
	  'errors were encountered in metered use file ',a)
	  call write_message(endspace='yes')
	end if
	go to 9900
 
9000	call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(infile)
9010	format(' Unable to read line ',a,' of metered use file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 9900
9200	write(amessage,9210) trim(outfile)
9210	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9900

9900    call close_files

end program mkpmp1

 