!     Last change:  JD   13 Feb 2001    4:33 pm
program smp2hyd

! -- Program SMP2HYD reads a bore sample file. It then writes a series of
!    output files from which bore hydrographs can be plotted for
!    user-specified bores.

	use defn
	use inter

	implicit none

	logical  :: active,lexist
	integer  :: ifail,sampunit,ibore,i,nbore,day0,mon0,year0,hour0,min0, &
	            sec0,daystart,monstart,yearstart,hourstart,minstart, &
	            secstart,dayfin,monfin,yearfin,hourfin,minfin,secfin, &
		    iline,cols,iactive,ndays,nsecs,itime,ierr,ndaystart,idate, &
		    nsecstart,ndayfin,nsecfin,ndayzero,nseczero,outunit,iwrite, &
                    iheader,ifail1,nbb
	real	 :: dayfac,time,secfac
	double precision value,rtemp,radd
	character (len=1)  :: atime,awindow,x,aa
	character (len=20) :: anum,atemp,prevbore,timestring
	character (len=80) :: specfile,coordfile,sampfile,afile,bfile

	integer, parameter			        :: TOTBORE = 500 !!was 30
	logical, allocatable, dimension(:)		:: recorded
	character (len=10), allocatable, dimension(:)	:: abore
	character (len=80), allocatable, dimension(:)	:: outfile


	write(amessage,5)
5       format(' Program SMP2HYD reads a bore sample file. It then writes a ',&
	'series of output files from which bore hydrographs can be plotted ',&
	'for user-specified bores.')
	call write_message(leadspace='yes',endspace='yes')

	allocate(recorded(TOTBORE),abore(TOTBORE),outfile(TOTBORE),stat=ierr)
	if(ierr.ne.0) go to 9300
	recorded=.false.
	secfac=1.0/86400.0
        prevbore=' '

	call readfig(specfile,coordfile,sampfile)

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

30      call open_named_input_file(ifail, &
	' Enter name of bore sample file: ',sampfile,sampunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

	write(6,*)
50	write(6,70)
70	format(' Enter bores for which hydrographs are ',&
	'required (Press <ENTER> if no more):-')
	iwrite=1
	ibore=0
100	ibore=ibore+1
	if(iwrite.eq.0) write(6,*)
110	call num2char(ibore,anum)
120	write(6,130,advance='no') trim(anum)
130	format('   Enter bore for hydrograph number ',a,': ')
	iwrite=0
	read(5,'(a)') atemp
	if(atemp.eq.' ')then
	  if(ibore.eq.1) go to 120
	  go to 225
	end if
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  if(ibore.eq.1) then
	    close(unit=sampunit,iostat=ierr)
	    if(ierr.ne.0) then
	      write(amessage,140)
140	      format(' File management error: cannot continue execution.')
	      go to 9890
	    end if
	    go to 30
	  else if(ibore.eq.2) then
	    go to 50
	  else
	    ibore=ibore-1
	    go to 110
	  end if
	end if
	if(ibore.gt.TOTBORE) then
	  call num2char(TOTBORE,anum)
	  write(amessage,105) trim(anum)
105	  format(' SMP2HYD can write hydrograph files for only ',a, &
	  ' bores at a time. To increase this number, edit the SMP2HYD ',&
	  'source code, assign a higher value to parameter TOTBORE ',&
	  'and recompile.')
	  call write_message(leadspace='yes')
	  go to 225
	end if
	if(len_trim(atemp).gt.10) then
	  write(amessage,160)
160	  format('   Bore identifier must be 10 characters or less in ',&
	  'length  - try again.')
	  call write_message
	  go to 120
	end if
	call casetrans(atemp,'hi')
	abore(ibore)=atemp
	if(ibore.ne.1) then
	  do i=1,ibore-1
	    if(abore(ibore).eq.abore(i))then
	      write(amessage,180)
180	      format('   This bore identifier has already been supplied  -',&
	      ' try again.')
	      call write_message
	      go to 120
	    end if
	  end do
	end if
200	write(6,220,advance='no') trim(anum)
220	format('   Enter output file for hydrograph number ',a,': ')
	read(5,'(a)') afile
	if(afile.eq.' ') go to 200
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0) then
	  write(6,*)
	  go to 120
	end if
        bfile=afile
        nbb=len_trim(bfile)
        call getfile(ifail1,bfile,afile,1,nbb)
        if(ifail1.ne.0) go to 200
	if(ibore.ne.1) then
	  do i=1,ibore-1
	    if(afile.eq.outfile(i)) then
	      write(amessage,223)
223	      format('   This filename has already been supplied  - try again.')
	      call write_message
	      go to 200
	    end if
	  end do
	end if
!	inquire(file=afile,exist=lexist)
!	if(lexist) then
!183	  write(6,185,advance='no')
!185	  format('   File already exists; overwrite it?  [y/n] ')
!	  read(5,'(a)') aa
!	  if(aa.eq.' ') go to 183
!	  call casetrans(aa,'lo')
!	  if(index(eschar,aa).ne.0) then
!	    write(6,*)
!	    go to 200
!	  end if
!	  if(aa.eq.'n') go to 200
!	  if(aa.ne.'y') go to 183
!	end if
	outfile(ibore)=afile
	go to 100

225	nbore=ibore-1
	call num2char(nbore,anum)
	write(amessage,227) trim(sampfile),trim(anum)
227	format(' - data will be extracted from bore sample file ',a, &
	' for ',a,' bores.')
	call write_message(leadspace='yes',endspace='yes')
230	write(6,240,advance='no')
240	format(' Use all samples for nominated bores, or specify sample ',&
	'window?  [a/w]: ')
	read(5,'(a)') awindow
	if(awindow.eq.' ') go to 230
	call casetrans(awindow,'lo')
	if(index(eschar,awindow).ne.0) then
	  write(6,*)
	  if(nbore.eq.1) then
	    go to 50
	  else
	    ibore=nbore
	    go to 110
	  end if
	end if
	if(awindow.eq.'a')then
	  daystart=-999
	  go to 300
	else if(awindow.ne.'w') then
	  go to 230
	end if
250	if(datespec.eq.1) then
	  write(6,260,advance='no')
260	  format('   Enter sample window start date  [dd/mm/yyyy]: ')
	else
	  write(6,261,advance='no')
261	  format('   Enter sample window start date  [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 250
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 230
	end if
	call char2date(ifail,atemp,daystart,monstart,yearstart)
	if(ifail.ne.0) then
	  write(6,420)
	  go to 250
	end if
	ndaystart=numdays(1,1,1970,daystart,monstart,yearstart)
265	write(6,266,advance='no')
266	format('   Enter sample window start time    [hh:mm:ss]: ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 265
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 250
	end if
	call char2time(ifail,atemp,hourstart,minstart,secstart)
	if(ifail.ne.0) then
	  write(6,520)
	  go to 265
	end if
	nsecstart=numsecs(0,0,0,hourstart,minstart,secstart)
270	if(datespec.eq.1) then
	  write(6,280,advance='no')
280	  format('   Enter sample window finish date [dd/mm/yyyy]: ')
	else
	  write(6,281,advance='no')
281	  format('   Enter sample window finish date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 270
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 265
	end if
	call char2date(ifail,atemp,dayfin,monfin,yearfin)
	if(ifail.ne.0) then
	  write(6,420)
	  go to 270
	end if
	ndayfin=numdays(1,1,1970,dayfin,monfin,yearfin)
	if(ndayfin.lt.ndaystart)then
	  write(amessage,282)
282	  format('   Sample window finish date must not precede starting ', &
	  'date  - try again.')
	  call write_message
	  go to 270
	end if
285	write(6,290,advance='no')
290	format('   Enter sample window finish time   [hh:mm:ss]: ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 285
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 270
	end if
	call char2time(ifail,atemp,hourfin,minfin,secfin)
	if(ifail.ne.0) then
	  write(6,520)
	  go to 285
	end if
	nsecfin=numsecs(0,0,0,hourfin,minfin,secfin)
	if((ndaystart.eq.ndayfin).and.(nsecfin.le.nsecstart))then
	  write(amessage,295)
295	  format('   Sample window finish time must follow starting ', &
	  'time  - try again.')
	  call write_message
	  go to 285
	end if
	
300	write(6,*)
330	write(6,350)
350	format(' When is zero time?')
370	if(datespec.eq.1) then
	  write(6,390,advance='no')
390	  format('   Enter reference date [dd/mm/yyyy]: ')
	else
	  write(6,391,advance='no')
391	  format('   Enter reference date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 370
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  if(awindow.eq.'w') then
	    go to 285
	  else
	    go to 230
	  end if
	end if
	call char2date(ifail,atemp,day0,mon0,year0)
	if(ifail.ne.0) then
	  write(6,420)
420	  format('   Illegal date  - try again.')
	  go to 370
	end if
	ndayzero=numdays(1,1,1970,day0,mon0,year0)
470	write(6,490,advance='no')
490	format('   Enter reference time   [hh:mm:ss]: ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 470
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 330
	end if
	call char2time(ifail,atemp,hour0,min0,sec0)
	if(ifail.ne.0) then
	  write(6,520)
520	  format('   Illegal time  - try again.')
	  go to 470
	end if
	nseczero=numsecs(0,0,0,hour0,min0,sec0)

	write(6,*)
550	write(6,560,advance='no')
560	format(' Enter output time units (yr/day/hr/min/sec) [y/d/h/m/s]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 550
	if(index(eschar,atime).ne.0)then
	  write(6,*)
	  go to 470
	end if
	call casetrans(atime,'lo')
	if(atime.eq.'y')then
	  dayfac=1.0/365.25
	  timestring='YEARS'
	else if(atime.eq.'d')then
	  dayfac=1.0
	  timestring='DAYS'
	else if(atime.eq.'h') then
	  dayfac=24.0
	  timestring='HOURS'
	else if(atime.eq.'m')then
	  dayfac=1440.0
	  timestring='MINUTES'
	else if(atime.eq.'s') then
	  dayfac=86400.0
	  timestring='SECONDS'
	else
	  go to 550
	end if

561     write(6,562,advance='no')
562     format(' Enter number to add to elapsed time column (<Enter> if none): ')
        read(5,'(a)') atemp
        if(atemp.eq.' ')then
          radd=0.0
        else
          call casetrans(atemp,'lo')
          if(atemp(1:2).eq.'e ')then
            write(6,*)
            go to 550
          end if
          call char2num(ifail,atemp,radd)
          if(ifail.ne.0) go to 561
        end if

	active=.false.
	iwrite=0
	itime=0
	iline=0
	imessage=0
	read_sample: do
	  iline=iline+1
	  if((iwrite.eq.0).and.(iline.eq.2000)) write(6,590)
590	  format(/,' Working.....')
	  read(sampunit,'(a)',err=9100,end=1000) cline
	  if(active) then
	    cols=5
	    call linesplit(ifail,5)
	    if(ifail.lt.0) cycle read_sample
	    if(ifail.gt.0) then
	      call linesplit(ifail,4)
	      if(ifail.ne.0) go to 9150
	      cols=4
	    end if
	    if(right_word(1)-left_word(1).gt.9) go to 9200
	    atemp=cline(left_word(1):right_word(1))
	    call casetrans(atemp,'hi')
	    if(atemp.eq.abore(iactive))then
	      call read_rest_of_sample_line(ifail,cols,ndays,nsecs,rtemp, &
	      iline,sampfile)
	      if(ifail.ne.0) go to 9900
	      if(daystart.ne.-999) then
		if((ndays.lt.ndaystart).or. &
		  ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) &
		  cycle read_sample
		if((ndays.gt.ndayfin).or. &
		  ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin)))then
	   	  close(unit=outunit,iostat=ierr)
		  if(imessage.eq.0) write(6,*)
		  call num2char(itime,anum)
		  write(amessage,640) trim(anum),trim(outfile(iactive)), &
		  trim(abore(iactive))
640		  format(' - ',a,' lines of data written to file ',a, &
		  ' for bore ',a)
		  call write_message(increment=1)
		  iwrite=1
		  if(itime.eq.0) then
	    	    outunit=nextunit()
	    	    open(unit=outunit,file=outfile(iactive),iostat=ierr)
	    	    if(ierr.ne.0) then
	      	      write(amessage,570) trim(outfile(iactive))
	      	      go to 9890
	    	    end if
	    	    write(outunit,620,err=9400) trim(timestring),'DATE', &
		    'TIME',trim(abore(iactive))
	    	    close(unit=outunit,iostat=ierr)
	   	  end if
		  itime=0
		  active=.false.
		  prevbore=abore(iactive)
		  cycle read_sample
		end if
	      end if
	      itime=itime+1
	      time=((ndays-ndayzero)+(nsecs-nseczero)*secfac)*dayfac
	      if(rtemp.lt.-1.0e38)then
	        value=char2double(ifail,4)
	        x='x'
	      else
	        value=rtemp
	        x=' '
	      end if
	      if(itime.eq.1) then
	        outunit=nextunit()
	        open(unit=outunit,file=outfile(iactive),iostat=ierr)
	        if(ierr.ne.0) then
	          write(amessage,570) trim(outfile(iactive))
570	          format(' Unable to open file ',a,' for output.')
	          go to 9890
	        end if
	        write(outunit,620,err=9400) trim(timestring),'DATE','TIME', &
	        trim(abore(iactive))
620	        format(1x,'TIME_IN_',a,t20,a,t35,a,t50,'BORE_',a)
	      end if
	      write(atemp,'(1pg14.6)') value
	      write(outunit,650,err=9400) time+radd, &
	      cline(left_word(2):right_word(2)), &
	      cline(left_word(3):right_word(3)), trim(atemp),x
650	      format(2x,1pg14.6,t20,a,t35,a,t48,a,a)
	      cycle read_sample
	    else
	      close(unit=outunit,iostat=ierr)
	      if(imessage.eq.0) write(6,*)
	      call num2char(itime,anum)
	      write(amessage,640) trim(anum),trim(outfile(iactive)), &
	      trim(abore(iactive))
	      call write_message(increment=1)
	      iwrite=1
	      if(itime.eq.0) then
	    	outunit=nextunit()
	    	open(unit=outunit,file=outfile(iactive),iostat=ierr)
	    	if(ierr.ne.0) then
	      	  write(amessage,570) trim(outfile(iactive))
	      	  go to 9890
	    	end if
	        write(outunit,620,err=9400) trim(timestring),'DATE', &
	        'TIME',trim(abore(iactive))
	        close(unit=outunit,iostat=ierr)
	      end if
	      active=.false.
	      itime=0
	      go to 580
	    end if
	  end if
	  call linesplit(ifail,1)
	  if(ifail.ne.0) cycle read_sample
	  if(right_word(1)-left_word(1).gt.9) go to 9200
	  atemp=cline(left_word(1):right_word(1))
	  call casetrans(atemp,'hi')
	  if(atemp.eq.prevbore) cycle read_sample
580	  do ibore=1,nbore
	    if(atemp.eq.abore(ibore)) go to 600
	  end do
	  prevbore=atemp
	  cycle read_sample
600	  active=.true.
	  iactive=ibore
	  recorded(iactive)=.true.
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_sample
	    if(ifail.gt.0) then
	    call linesplit(ifail,4)
	    if(ifail.ne.0) go to 9150
	    cols=4
	  end if
	  call read_rest_of_sample_line(ifail,cols,ndays,nsecs,rtemp,iline, &
	  sampfile)
	  if(ifail.ne.0) go to 9900
	  if(daystart.ne.-999) then
	    if((ndays.lt.ndaystart).or. &
	      ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) &
	      cycle read_sample
	    if((ndays.gt.ndayfin).or. &
	      ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin))) then
		active=.false.
		prevbore=abore(iactive)
		itime=0
		call num2char(itime,anum)
		if(imessage.eq.0) write(6,*)
		write(amessage,640) trim(anum),trim(outfile(iactive)), &
		trim(abore(iactive))
		call write_message(increment=1)
		iwrite=1
	        if(itime.eq.0) then
	    	  outunit=nextunit()
	    	  open(unit=outunit,file=outfile(iactive),iostat=ierr)
	    	  if(ierr.ne.0) then
	      	    write(amessage,570) trim(outfile(iactive))
	      	    go to 9890
	    	  end if
	          write(outunit,620,err=9400) trim(timestring),'DATE', &
	          'TIME',trim(abore(iactive))
	          close(unit=outunit,iostat=ierr)
	        end if
		cycle read_sample
	    end if
	  end if
	  itime=1
	  outunit=nextunit()
	  open(unit=outunit,file=outfile(iactive),iostat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,570) trim(outfile(iactive))
	    go to 9890
	  end if
	  write(outunit,620,err=9400) trim(timestring),'DATE','TIME', &
	  trim(abore(iactive))
	  recorded(iactive)=.true.
	  time=((ndays-ndayzero)+(nsecs-nseczero)*secfac)*dayfac
	  if(rtemp.lt.-1.0e38)then
	    value=char2double(ifail,4)
	    x='x'
	  else
	    value=rtemp
	    x=' '
	  end if
	  write(atemp,'(1pg14.6)') value
	  write(outunit,650,err=9400) time+radd, &
	  cline(left_word(2):right_word(2)), &
	  cline(left_word(3):right_word(3)),trim(atemp),x
	end do read_sample

1000	if(active) then
	  close(unit=outunit,iostat=ierr)
	  if(imessage.eq.0) write(6,*)
	  call num2char(itime,anum)
	  write(amessage,640) trim(anum),trim(outfile(iactive)), &
	  trim(abore(iactive))
	  call write_message(increment=1)
	  if(itime.eq.0) then
	    outunit=nextunit()
	    open(unit=outunit,file=outfile(iactive),iostat=ierr)
	    if(ierr.ne.0) then
	      write(amessage,570) trim(outfile(iactive))
	      go to 9890
	    end if
	    write(outunit,620,err=9400) trim(timestring),'DATE', &
	    'TIME',trim(abore(iactive))
	    close(unit=outunit,iostat=ierr)
	  end if
	end if

	call num2char(iline,anum)
	write(amessage,1020) trim(anum),trim(sampfile)
1020	format(' - ',a,' lines read from bore sample file ',a)
	call write_message(leadspace='yes')
	imessage=0
	do ibore=1,nbore
	  if(.not.recorded(ibore)) then
	    if(imessage.eq.0) write(6,*)
	    write(amessage,1030) trim(abore(ibore)),trim(sampfile)
1030	    format(' Warning: no entries found for bore ',a,' in bore sample ',&
	    'file ',a)
	    call write_message(increment=1)
	    outunit=nextunit()
	    open(unit=outunit,file=outfile(ibore),iostat=ierr)
	    if(ierr.ne.0) then
	      write(amessage,570) trim(outfile(ibore))
	      go to 9890
	    end if
	    write(outunit,620,err=9400) trim(timestring),'DATE','TIME', &
	    trim(abore(ibore))
	    close(unit=outunit,iostat=ierr)
	  end if
	end do
	go to 9900

9100	call num2char(iline,anum)
	write(amessage,9110) trim(anum),trim(sampfile)
9110	format(' Error reading line ',a,' of bore sample file ',a)
	go to 9890
9150	call num2char(iline,anum)
	write(amessage,9160) trim(anum),trim(sampfile)
9160	format(' Error on line ',a,' of bore sample file ',a,': ',&
	'insufficient entries.')
	go to 9890
9200	call num2char(iline,anum)
	write(amessage,9210) trim(anum),trim(sampfile)
9210	format(' Error on line ',a,' of bore sample file ',a,': bore ',&
	'identifier greater than 10 characters in length.')
	go to 9890
9300	write(amessage,9310)
9310	format(' Cannot allocate sufficient memory to run program SMP2HYD.')
	go to 9890
9400	write(amessage,9410) trim(outfile(iactive))
9410	format(' Error writing to file ',a,': file inaccessible or disk full.')
	go to 9890

9890	call write_message(leadspace='yes')
9900	call close_files
	deallocate(recorded,abore,outfile,stat=ierr)
	write(6,*)

end program smp2hyd

 