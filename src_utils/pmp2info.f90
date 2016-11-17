!     Last change:  JD   17 Dec 2000   10:01 pm
program pmp2info

! -- Program PMP2INFO evaluates the water pumped from a list of bores between
!    two user-specifed times, writing the results to a bore information file.

	use defn
	use inter

	implicit none

	integer	:: ifail,pumpunit,outunit,idate, iheader, &
	 	   iline,icollect,ibore,i,ierr,ilist,j,icoord,&
		   dd,mm,yy,intday1,intday2,hhh,mmm,sss,inttime1,inttime2
	character (len=80)	:: specfile,sampfile,aprompt,outfile,pumpfile
	character(len=15)	:: aline,startdate,findate,starttime,fintime
	character(len=11)	:: aeast,anorth
	character(len=10)	:: abore,aoldbore
	character(len=30)	:: aval
	character(len=1)	:: aout
	double precision	:: pumped1,pumped2
	real, allocatable,dimension(:)	:: pmpinterp

	integer, parameter				:: NUM_PMP_BORE=3000
	integer, dimension(NUM_PMP_BORE)		:: ndays,nsecs
	double precision, dimension(NUM_PMP_BORE)	:: pumped


	write(amessage,5)
5       format(' Program PMP2INFO evaluates the water pumped from a list ',&
	'of bores between two user-specifed times, writing the results to ',&
	'a bore information file.')
	call write_message(leadspace='yes',endspace='yes')

	call readfig(specfile,bore_coord_file,sampfile,pumpfile)

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message
	  go to 9900
	endif
	if((idate.ne.0).or.(datespec.eq.0)) then
	  write(amessage,9)
9	  format(' Cannot read date format from settings file ', &
	  'settings.fig')
	  call write_message
 	  go to 9900
	end if

30      call read_bore_coord_file(ifail, &
       ' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

50     call read_bore_list_file(ifail, &
       ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 30
	end if

	allocate(pmpinterp(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,60)
60	  format(' Cannot allocate sufficient memory to run PMP2INFO.')
	  call write_message(leadspace='yes')
	  go to 9900
	end if
	pmpinterp=-6.1e37

70	call open_named_input_file(ifail, &
	' Enter name of bore pumping file: ',pumpfile,pumpunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  deallocate(pmpinterp,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,80)
80	    format(' Memory management error: cannot continue execution.')
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  go to 50
	end if

300	write(6,*)
310	if(datespec.eq.1) then
	  write(6,320,advance='no')
320	  format(' Enter time interval starting date  [dd/mm/yyyy]: ')
	else
	  write(6,321,advance='no')
321	  format(' Enter time interval starting date  [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') startdate
	if(startdate.eq.' ') go to 310
	startdate=adjustl(startdate)
	if(index(eschar,startdate(1:2)).ne.0) then
	  write(6,*)
	  close(unit=pumpunit,err=9100)
	  go to 70
	end if
	call char2date(ifail,startdate,dd,mm,yy)
	if(ifail.ne.0) then
	  write(amessage,325)
325	  format(' Illegal date  - try again.')
	  call write_message
	  go to 310
	end if
	intday1=numdays(1,1,1970,dd,mm,yy)

510	write(6,520,advance='no')
520	format(' Enter time interval starting time  [hh:mm:ss]:   ')
	read(5,'(a)') starttime
	if(starttime.eq.' ') go to 510
	starttime=adjustl(starttime)
	if(index(eschar,starttime(1:2)).ne.0) go to 300
	call char2time(ifail,starttime,hhh,mmm,sss)
	if(ifail.ne.0) then
	  write(amessage,525)
525	  format(' Illegal time  - try again.')
	  call write_message
	  go to 510
	end if
	inttime1=numsecs(0,0,0,hhh,mmm,sss)

410	if(datespec.eq.1) then
	  write(6,420,advance='no')
420	  format(' Enter time interval finishing date [dd/mm/yyyy]: ')
	else
	  write(6,421,advance='no')
421	  format(' Enter time interval finishing date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') findate
	if(findate.eq.' ') go to 410
	findate=adjustl(findate)
	if(index(eschar,findate(1:2)).ne.0) then
	  write(6,*)
	  go to 510
	end if
	call char2date(ifail,findate,dd,mm,yy)
	if(ifail.ne.0) then
	  write(amessage,325)
	  call write_message
	  go to 410
	end if
	intday2=numdays(1,1,1970,dd,mm,yy)
	if(intday2.lt.intday1) then
	  write(amessage,425)
425	  format(' Finishing date must equal or postdate starting ', &
	  'date  - try again.')
	  call write_message
	  go to 410
	end if

610	write(6,620,advance='no')
620	format(' Enter time interval finishing time [hh:mm:ss]:   ')
	read(5,'(a)') fintime
	if(fintime.eq.' ') go to 610
	fintime=adjustl(fintime)
	if(index(eschar,fintime(1:2)).ne.0) then
	  write(6,*)
	  go to 410
	end if
	call char2time(ifail,fintime,hhh,mmm,sss)
	if(ifail.ne.0) then
	  write(amessage,525)
	  call write_message
	  go to 610
	end if
	inttime2=numsecs(0,0,0,hhh,mmm,sss)

	if((intday2.lt.intday1).or. &
	  ((intday2.eq.intday1).and.(inttime2.le.inttime1))) then
	  write(amessage,335)
335	  format(' Finishing date/time must postdate starting ',&
	  'date/time  - try again.')
	  call write_message
	  go to 410
	end if

330	write(6,*)
	aprompt=' Enter name for output bore information file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 610
	end if	

340	write(6,350,advance='no')
350	format(' Record any uninterpolated bores to output file?  [y/n]: ')
	read(5,'(a)') aout
	if(aout.eq.' ') go to 340
	call casetrans(aout,'lo')
	if(aout.eq.'e')then
	  close(unit=outunit)
	  go to 330
	end if
	if((aout.ne.'y').and.(aout.ne.'n')) go to 340

	icollect=0
	iline=0
	abore=' '
	read_line: do
	  aoldbore=abore
	  iline=iline+1
	  if(iline.eq.1000) write(6,95)
95	  format(/,' Working.....')
	  read(pumpunit,'(a)',err=9000,end=1000) cline
	  call linesplit(ifail,6)
	  if(ifail.lt.0) cycle read_line
	  if(ifail.gt.0)then
	    call num2char(iline,aline)
	    write(amessage,100) trim(aline),trim(pumpfile)
100	    format('insufficient items on line ',a,' of bore pumping file ',a)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if

	  if(right_word(1)-left_word(1).gt.9) then
	    call num2char(iline,aline)
	    write(amessage,130) trim(aline),trim(pumpfile)
130	    format('bore identifier greater than 10 characters in length ',&
	    'on line ',a,' of bore pumping file ',a)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  abore=cline(left_word(1):right_word(1))
	  call casetrans(abore,'hi')
	  if(abore.eq.aoldbore)then
	    if(icollect.eq.1)then
	      ibore=ibore+1
	      if(ibore.gt.NUM_PMP_BORE)then
		call num2char(NUM_PMP_BORE,aline)
		write(amessage,140) trim(pumpfile),trim(abore),trim(aline)
140		format(' There are more entries in file ',a,' for bore ',&
		a,' than the maximum allowed number of ',a,&
		': edit source code file and assign a greater value to ',&
		'parameter NUM_PMP_BORE.  Then recompile program.')
		call write_message(leadspace='yes')
		go to 9900
	      end if
	      call read_rest_of_pump_line(ifail,ibore,ndays,nsecs,&
	      pumped,iline,pumpfile)
	      if(ifail.ne.0) go to 9900	      
	    else
	      cycle read_line
	    end if
	  else
	    if(icollect.eq.1)then
	      if((pmpinterp(ilist).lt.-6.2e37).or.&
		 (pmpinterp(ilist).gt.-6.05e37))then
		 write(amessage,143) trim(bore_list_id(ilist)),&
		 trim(pumpfile)
143		 format('entries for bore ',a,' in bore pumping file ',&
		 a,' are not all in sequence.')
		 call write_message(error='yes',leadspace='yes')
	    	 go to 9900
	      end if
	      call time_interp(ifail,ibore,ndays,nsecs,pumped,intday1, &
	      inttime1,huge(1.0),0.0,pumped1)
	      if(ifail.eq.1)then
		write(amessage,145) trim(bore_list_id(ilist)),trim(pumpfile)
145		format('non-increasing pumping times for bore ',a,&
		' in bore pumping file ',a)
		call write_message(error='yes',leadspace='yes')
	 	go to 9900
	      end if
	      if(pumped1.gt.-8.05e37) then
		call time_interp(ifail,ibore,ndays,nsecs,pumped,intday2, &
	        inttime2,huge(1.0),0.0,pumped2)
		if(pumped2.gt.-9.05e37) then
		  pmpinterp(ilist)=pumped2-pumped1
		else
		  pmpinterp(ilist)=-9.1e37
		end if
	      else
		pmpinterp(ilist)=pumped1
	      end if
	      icollect=0
	      do i=1,num_bore_list
		if(i.eq.ilist) cycle
		if(bore_list_id(i).eq.bore_list_id(ilist)) pmpinterp(i)=pmpinterp(ilist)
	      end do
	    end if
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) then
		ilist=i
		icollect=1
		ibore=1
		call read_rest_of_pump_line(ifail,ibore,ndays,nsecs,pumped,&
		iline,pumpfile)
		if(ifail.ne.0) go to 9900
		go to 150
	      end if
	    end do
	    icollect=0
150	    continue
	  end if
	end do read_line

1000	if(icollect.eq.1)then
	  if((pmpinterp(ilist).lt.-6.2e37).or.&
	  (pmpinterp(ilist).gt.-6.05e37))then
	    write(amessage,143) trim(bore_list_id(ilist)),&
	    trim(pumpfile)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  call time_interp(ifail,ibore,ndays,nsecs,pumped,intday1, &
	  inttime1,huge(1.0),0.0,pumped1)
	  if(ifail.eq.1)then
	    write(amessage,145) trim(bore_list_id(ilist)),trim(pumpfile)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  if(pumped1.gt.-8.05e37) then
	    call time_interp(ifail,ibore,ndays,nsecs,pumped,intday2, &
	    inttime2,huge(1.0),0.0,pumped2)
	    if(pumped2.gt.-9.05e37) then
	      pmpinterp(ilist)=pumped2-pumped1
	    else
	      pmpinterp(ilist)=-9.1e37
	    end if
	  else
	    pmpinterp(ilist)=pumped1
	  end if
	  do i=1,num_bore_list
	    if(i.eq.ilist) cycle
	    if(bore_list_id(i).eq.bore_list_id(ilist)) pmpinterp(i)=pmpinterp(ilist)
	  end do
	end if

	ibore=0
	do i=1,num_bore_list
	  if((aout.eq.'n').and.(pmpinterp(i).lt.-5.0e37)) cycle
	  ibore=ibore+1
	  icoord=0
	  do j=1,num_bore_coord
	    if(bore_list_id(i).eq.bore_coord_id(j)) then
	      icoord=1
	      exit
	    end if
	  end do
	  if(icoord.eq.0)then
	    aeast='no_coords'
	    anorth='no_coords'
	  else
	    call num2char(bore_coord_east(j),aeast)
	    call num2char(bore_coord_north(j),anorth)
	  end if
	  if(pmpinterp(i).lt.-9.0e37) then
	    aval='after_last_sample'
	  else if(pmpinterp(i).lt.-8.0e37)then
	    aval='before_first_sample'
	  else if(pmpinterp(i).lt.-6.0e37)then
	    aval='not_in_pumping_file'
	  else
	    write(aval,'(1pg12.5)') pmpinterp(i)
!	    call num2char(pmpinterp(i),aval,9)
	  end if
	  write(outunit,1050) bore_list_id(i),aeast,anorth,trim(aval)
1050	  format(1x,a10,2x,a11,2x,a11,2x,a)
	end do

	call num2char(ibore,aline)
	write(amessage,1080) trim(aline),trim(outfile)
1080	format('  - data for ',a,' bores written to file ',a)
	call write_message(leadspace='yes')
	go to 9900

9000	call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(pumpfile)
9010	format(' Error reading line ',a,' of bore pumping file ',a)
	call write_message(leadspace='yes')
	go to 9900
9100	write(amessage,9110)
9110	format(' File management error: cannot continue operation.')
	call write_message(leadspace='yes')
	go to 9900

9900    call close_files
	call free_bore_mem
	deallocate(pmpinterp,stat=ierr)
	write(6,*)

end program pmp2info

 