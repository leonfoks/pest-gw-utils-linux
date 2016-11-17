!     Last change:  JD   16 Dec 2000   11:06 pm
program smp2info

! -- Program SMP2INFO interpolates data contained in a bore sample file to a
!    user-specified date and time.

	use defn
	use inter

	implicit none

	integer	:: ifail,sampunit,itemp,dd,mm,yy,intday,outunit,&
	 	   iline,cols,icollect,ibore,i,ierr,ilist,j,icoord,intsec, &
		   hhh,mmm,sss,idate,iheader
	real	:: rnear,rrtemp,rconst
	character (len=80)	:: specfile,sampfile,aprompt,outfile
	character(len=15)	:: adate,aline,atime
	character(len=12)	:: aeast,anorth
	character(len=10)	:: abore,aoldbore
	character(len=30)	:: aval
	character(len=1)	:: aout
	double precision, allocatable,  dimension(:)   :: valinterp

	integer, parameter			:: NUM_SAMP_BORE=10000
	integer, dimension(NUM_SAMP_BORE)	:: ndays,nsecs
	double precision, dimension(NUM_SAMP_BORE)  :: value


	write(amessage,5)
5       format(' Program SMP2INFO interpolates data contained in a bore ', &
	'sample file to a user-specified date and time.')
	call write_message(leadspace='yes',endspace='yes')

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          call write_message
          go to 9900
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          call write_message
          go to 9900
        endif
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
          'settings.fig')
          call write_message
          go to 9900
        end if

	call readfig(specfile,bore_coord_file,sampfile)

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

	allocate(valinterp(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,60)
60	  format(' Cannot allocate sufficient memory to run SMP2INFO.')
	  call write_message(leadspace='yes')
	  go to 9900
	end if
	valinterp=-6.1e37

70	call open_named_input_file(ifail, &
	' Enter name of bore sample file: ',sampfile,sampunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  deallocate(valinterp,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,80)
80	    format(' Memory management error: cannot continue execution.')
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	  go to 50
	end if

190	write(6,*)
200	write(6,210)
210	format(' Enter the following time-interpolation parameters --->')
215	write(6,230,advance='no')
230	format('   maximum days to nearest sample  [90.000]: ')
	itemp=key_read(rnear)
	if(escset.ne.0)then
	  escset=0
	  close(unit=sampunit)
	  write(6,*)
	  go to 70
	end if
	if(itemp.lt.0)then
	  rnear=90.0
	else if(itemp.gt.0) then
	  go to 215
	endif
	if(nneg_test(rnear,'number of days').ne.0) go to 215

250	write(6,260)
260	format('   days over which sample can be assumed constant if linear ',&
	'interpolation')
	rrtemp=min(rnear,3.0000)
	call num2char(rrtemp,aline,6)
	write(6,261,advance='no') trim(aline)
261	format('   cannot take place  [',a,']: ')
	itemp=key_read(rconst)
	if(escset.ne.0)then
	  escset=0
	  go to 190
	end if
	if(itemp.lt.0)then
	  rconst=rrtemp
	else if(itemp.gt.0) then
	  go to 250
	endif
	if(nneg_test(rconst,'number of days').ne.0) go to 250
	if(rconst.gt.rnear) then
	  write(amessage,270)
270	  format('this number cannot exceed days to nearest sample - try again.')
	  call write_message(error='yes')
	  go to 250
	end if

300	write(6,*)
310	if(datespec.eq.1) then
	  write(6,320,advance='no')
320	  format(' Enter interpolation date [dd/mm/yyyy]: ')
	else
	  write(6,321,advance='no')
321	  format(' Enter interpolation date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') adate
	if(adate.eq.' ') go to 310
	adate=adjustl(adate)
	if(index(eschar,adate(1:2)).ne.0) then
	  write(6,*)
	  go to 250
	end if
	call char2date(ifail,adate,dd,mm,yy)
	if(ifail.ne.0) then
	  write(amessage,325)
325	  format(' Illegal date  - try again.')
	  call write_message
	  go to 310
	end if
	intday=numdays(1,1,1970,dd,mm,yy)

410	write(6,420,advance='no')
420	format(' Enter interpolation time   [hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 410
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) go to 300
	call char2time(ifail,atime,hhh,mmm,sss)
	if(ifail.ne.0) then
	  write(amessage,425)
425	  format(' Illegal time  - try again.')
	  call write_message
	  go to 410
	end if
	intsec=numsecs(0,0,0,hhh,mmm,sss)

330	write(6,*)
	aprompt=' Enter name for output bore information file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 410
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
	  if(iline.eq.2000) write(6,95)
95	  format(/,' Working.....')
	  read(sampunit,'(a)',err=9000,end=1000) cline
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_line
	  if(ifail.gt.0)then
	    cols=4
	    call linesplit(ifail,4)
	    if(ifail.ne.0)then
	      call num2char(iline,aline)
	      write(amessage,100) trim(aline),trim(sampfile)
100	      format('insufficient items on line ',a,' of bore sample file ',a)
	      call write_message(error='yes',leadspace='yes')
	      go to 9900
	    end if
	  end if

	  if(right_word(1)-left_word(1).gt.9) then
	    call num2char(iline,aline)
	    write(amessage,130) trim(aline),trim(sampfile)
130	    format('bore identifier greater than 10 characters in length ',&
	    'on line ',a,' of bore sample file ',a)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  abore=cline(left_word(1):right_word(1))
	  call casetrans(abore,'hi')
	  if(abore.eq.aoldbore)then
	    if(icollect.eq.1)then
	      ibore=ibore+1
	      if(ibore.gt.NUM_SAMP_BORE)then
		call num2char(NUM_SAMP_BORE,aline)
		write(amessage,140) trim(sampfile),trim(abore),trim(aline)
140		format(' There are more samples in file ',a,' for bore ',&
		a,' than the maximum allowed number of ',a,&
		': edit source code file and assign a greater value to ',&
		'parameter NUM_SAMP_BORE.  Then recompile program.')
		call write_message(leadspace='yes')
		go to 9900
	      end if
	      call read_rest_of_sample_line(ifail,cols,ndays(ibore), &
	      nsecs(ibore),value(ibore),iline,sampfile)
	      if(ifail.ne.0) go to 9900	      
	    else
	      cycle read_line
	    end if
	  else
	    if(icollect.eq.1)then
	      if((valinterp(ilist).lt.-6.2e37).or.&
		 (valinterp(ilist).gt.-6.05e37))then
		 write(amessage,143) trim(bore_list_id(ilist)),&
		 trim(sampfile)
143		 format('entries for bore ',a,' in bore sample file ',&
		 a,' are not all in sequence.')
		 call write_message(error='yes',leadspace='yes')
	    	 go to 9900
	      end if
	      call time_interp(ifail,ibore,ndays,nsecs,value,intday,intsec, &
	      rnear,rconst,valinterp(ilist))
	      if(ifail.eq.1)then
		write(amessage,145) trim(bore_list_id(ilist)),trim(sampfile)
145		format('non-increasing sample times for bore ',a,&
		' in bore sample file ',a)
		call write_message(error='yes',leadspace='yes')
	 	go to 9900
	      end if
	      icollect=0
	      do i=1,num_bore_list
		if(i.eq.ilist) cycle
		if(bore_list_id(i).eq.bore_list_id(ilist)) valinterp(i)=valinterp(ilist)
	      end do
	    end if
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) then
		ilist=i
		icollect=1
		ibore=1
		call read_rest_of_sample_line(ifail,cols,ndays(1),nsecs(1), &
		value(1),iline,sampfile)
		if(ifail.ne.0) go to 9900
		go to 150
	      end if
	    end do
	    icollect=0
150	    continue
	  end if
	end do read_line

1000	if(icollect.eq.1)then
	  if((valinterp(ilist).lt.-6.2e37).or.&
	  (valinterp(ilist).gt.-6.05e37))then
	    write(amessage,143) trim(bore_list_id(ilist)),&
	    trim(sampfile)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  call time_interp(ifail,ibore,ndays,nsecs,value,intday,intsec, &
	  rnear,rconst,valinterp(ilist))
	  if(ifail.eq.1)then
	    write(amessage,145) trim(bore_list_id(ilist)),trim(sampfile)
	    call write_message(error='yes',leadspace='yes')
	    go to 9900
	  end if
	  do i=1,num_bore_list
	    if(i.eq.ilist) cycle
	    if(bore_list_id(i).eq.bore_list_id(ilist)) &
	    valinterp(i)=valinterp(ilist)
	  end do
	end if

	ibore=0
	do i=1,num_bore_list
	  if((aout.eq.'n').and.(valinterp(i).lt.-5.0e37)) cycle
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
	  if(valinterp(i).lt.-1.0e38)then
	    aval='x-affected'
	  else if(valinterp(i).lt.-9.0e37) then
	    aval='after_last_sample'
	  else if(valinterp(i).lt.-8.0e37)then
	    aval='before_first_sample'
	  else if(valinterp(i).lt.-7.0e37)then
	    aval='too_long_to_nearest_sample'
	  else if(valinterp(i).lt.-6.0e37)then
	    aval='not_in_sample_file'
	  else
	    write(aval,'(1pg14.7)')valinterp(i)
!	    call num2char(valinterp(i),aval,9)
	  end if
	  write(outunit,1050) bore_list_id(i),aeast,anorth,trim(aval)
1050	  format(1x,a10,2x,a12,2x,a12,2x,a)
	end do

	call num2char(ibore,aline)
	write(amessage,1080) trim(aline),trim(outfile)
1080	format('  - data for ',a,' bores written to file ',a)
	call write_message(leadspace='yes')
	go to 9900

9000	call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(sampfile)
9010	format(' Error reading line ',a,' of bore sample file ',a)
	call write_message(leadspace='yes')
	go to 9900
9900    call close_files
	call free_bore_mem
	deallocate(valinterp,stat=ierr)
	write(6,*)

end program smp2info

 