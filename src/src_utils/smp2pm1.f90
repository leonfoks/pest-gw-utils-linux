!     Last change:  JD   18 Dec 2000    0:24 am
program smp2pm1

! -- Program SMP2PM1 reads a bore sample file. It re-writes the data
!    contained in this file in PMWIN-compatible form as borehole observation
!    data.


	use defn
	use inter

	implicit none

	logical  :: active
	integer  :: ifail,sampunit,ibore,i,day0,mon0,year0,hour0,min0,j, &
		    sec0,daystart,monstart,yearstart,hourstart,minstart, &
		    secstart,dayfin,monfin,yearfin,hourfin,minfin,secfin, &
		    iline,cols,iactive,ndays,nsecs,itime,ierr,ndaystart, &
		    nsecstart,ndayfin,nsecfin,ndayzero,nseczero, &
		    outunit1,outunit2,outunit3,nobs,numbore,iibore,iobstart,&
	            iobsfin,idate,iheader
        integer, dimension(13)  		:: icolour
	integer, allocatable, dimension(:)	:: iobs_out
	real     				:: dayfac,time,secfac
	double precision 			:: value
	character (len=1)  :: awindow,atype,atime
	character (len=20) :: anum,atemp,prevbore
	character (len=80) :: specfile,sampfile,outfile1, &
			      outfile2,outfile3,aprompt

	integer, parameter		:: MAX_PM_BOR=1000
	integer, parameter		:: MAX_PM_OBS=6000
	real, dimension(MAX_PM_OBS)	:: time_out,value_out

        data icolour /255,16711680,65535,16776960,65280,8388736,12632256, &
        16711935,32896,8388608,4227072,128,0/


	write(amessage,5)
5       format(' Program SMP2PM1 reads a bore sample file. It re-writes the ',&
	'data contained in this file in PMWIN-compatible form as ', &
	'borehole observation data.')
	call write_message(leadspace='yes',endspace='yes')

	secfac=1.0/86400.0

	call readfig(specfile,bore_coord_file,sampfile)

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
32      write(6,33,advance='no')
33      format(' Does this file contain head, drawdown, concentration, ',&
	'compaction, ',/,'   preconsolidated head, or subsidence data?  ', &
	'[h/d/c/m/p/s]: ')
	read(5,'(a)') atype
	if(atype.eq.' ') go to 32
	if(index(eschar,atype).ne.0) then
	  write(6,*)
	  close(unit=sampunit,err=9500)
	  go to 30
	end if
	call casetrans(atype,'lo')
	if((atype.ne.'h').and.(atype.ne.'d').and.(atype.ne.'c').and. &
	(atype.ne.'m').and.(atype.ne.'p').and.(atype.ne.'s')) go to 32

	write(6,*)
35      call read_bore_coord_file(ifail, &
	' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 32
	end if

	do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) go to 36
	end do
	go to 49
36      write(amessage,38) trim(bore_coord_file)
38      format(' The following bores were not provided with layer numbers in ',&
	'bore coordinates file ',a,' ----->')
	call write_message(leadspace='yes')
	j=2
	imessage=0
	write_bores: do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) then
	    write(amessage(j:),'(a)') trim(bore_coord_id(i))
	    j=j+11
	    if(j.ge.71) then
	      call write_message(increment=1)
	      if(imessage.gt.12) goto 9900
	      j=2
	    end if
	  end if
	end do write_bores
	if(j.ne.2) call write_message
	go to 9900

49      write(6,*)
50      call read_bore_list_file(ifail, &
	' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 35
	end if

	if(num_bore_list.gt.MAX_PM_BOR) then
	  call num2char(MAX_PM_BOR,anum)
	  write(amessage,55) trim(anum),trim(bore_list_file)
55        format(' PMWIN can read a maximum of ',a,' bores: there are more ', &
	  'than this cited in bore listing file ',a)
	  go to 9890
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
		write(amessage,60) trim(bore_list_file)
60              format(' Execution of program SMP2PM1 cannot proceed as ',&
		'there are multiple occurrences of the same bore in bore ',&
		'listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	allocate(iobs_out(num_bore_list),stat=ierr)
	if(ierr.ne.0) go to 9300

	write(6,*)
230     write(6,240,advance='no')
240     format(' Use all samples for listed bores, or specify sample ',&
	'window?  [a/w]: ')
	read(5,'(a)') awindow
	if(awindow.eq.' ') go to 230
	call casetrans(awindow,'lo')
	if(index(eschar,awindow).ne.0) then
	  deallocate(iobs_out,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,245)
245	    format(' Memory allocation error: cannot continue execution.')
	    go to 9890
	  end if
	  go to 49
	end if
	if(awindow.eq.'a')then
	  daystart=-99999
	  go to 300
	else if(awindow.ne.'w') then
	  go to 230
	end if
250	if(datespec.eq.1) then
	  write(6,260,advance='no')
260       format('   Enter sample window start date  [dd/mm/yyyy]: ')
	else
	  write(6,261,advance='no')
261       format('   Enter sample window start date  [mm/dd/yyyy]: ')
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
265     write(6,266,advance='no')
266     format('   Enter sample window start time    [hh:mm:ss]: ')
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
280       format('   Enter sample window finish date [dd/mm/yyyy]: ')
	else
	  write(6,281,advance='no')
281       format('   Enter sample window finish date [mm/dd/yyyy]: ')
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
282       format('   Sample window finish date must not precede starting ', &
	  'date  - try again.')
	  call write_message
	  go to 270
	end if
285     write(6,290,advance='no')
290     format('   Enter sample window finish time   [hh:mm:ss]: ')
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
295       format('   Sample window finish time must follow starting ', &
	  'time  - try again.')
	  call write_message
	  go to 285
	end if
	
300     write(6,*)
330     write(6,350)
350     format(' When is zero time?')
370	if(datespec.eq.1) then
	  write(6,390,advance='no')
390       format('   Enter reference date [dd/mm/yyyy]: ')
	else
	  write(6,391,advance='no')
391       format('   Enter reference date [mm/dd/yyyy]: ')
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
420       format('   Illegal date  - try again.')
	  go to 370
	end if
	ndayzero=numdays(1,1,1970,day0,mon0,year0)
470     write(6,490,advance='no')
490     format('   Enter reference time   [hh:mm:ss]: ')
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
520       format('   Illegal time  - try again.')
	  go to 470
	end if
	nseczero=numsecs(0,0,0,hour0,min0,sec0)

	write(6,*)
550     write(6,560,advance='no')
560     format(' Enter output time units (yr/day/hr/min/sec) [y/d/h/m/s]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 550
	if(index(eschar,atime).ne.0)then
	  write(6,*)
	  go to 470
	end if
	call casetrans(atime,'lo')
	if(atime.eq.'y')then
	  dayfac=1.0/365.25
	else if(atime.eq.'d')then
	  dayfac=1.0
	else if(atime.eq.'h') then
	  dayfac=24.0
	else if(atime.eq.'m')then
	  dayfac=1440.0
	else if(atime.eq.'s') then
	  dayfac=86400.0
	else
	  go to 550
	end if

1500    write(6,*)
	write(6,1505)
1505    format(' Provide the name for the following SMP2PM1 output files:-')
	aprompt=' Enter name for output PMWIN observation file: '
1510    call open_output_file(ifail,aprompt,outfile1,outunit1)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 550
	end if  

1530    aprompt=' Enter name for output PMWIN bore coordinates file: '
	call open_output_file(ifail,aprompt,outfile2,outunit2)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=outunit1,err=9350)
	  go to 1500
	end if

1550    aprompt=' Enter name for output bore_id to PMWIN_id file: '
	call open_output_file(ifail,aprompt,outfile3,outunit3)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=outunit2,err=9350)
	  write(6,*)
	  go to 1530
	end if

	iobs_out=-999
	nobs=0
	active=.false.
	itime=0
	iline=0
	prevbore=' '
	read_sample: do
	  iline=iline+1
	  if(iline.eq.2000) write(6,590)
590       format(/,' Working.....')
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
	    if(atemp.eq.bore_list_id(iactive))then
	      call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value, &
	      iline,sampfile)
	      if(ifail.ne.0) go to 9900
	      if(value.lt.-1.0e38) cycle read_sample
	      if(daystart.ne.-99999) then
		if((ndays.lt.ndaystart).or. &
		  ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) &
		  cycle read_sample
		if((ndays.gt.ndayfin).or. &
		  ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin)))then
		  itime=0
		  active=.false.
		  prevbore=bore_list_id(iactive)
		  cycle read_sample
		end if
	      end if
	      time=((ndays-ndayzero)+(nsecs-nseczero)*secfac)*dayfac
	      itime=itime+1
	      nobs=nobs+1
	      if(nobs.gt.MAX_PM_OBS) go to 9700
	      time_out(nobs)=time
	      value_out(nobs)=value
	      if(itime.eq.1) iobs_out(iactive)=nobs
	      cycle read_sample
	    else
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
580       do ibore=1,num_bore_list
	    if(atemp.eq.bore_list_id(ibore)) go to 600
	  end do
	  prevbore=atemp
	  cycle read_sample
600       active=.true.
	  iactive=ibore
	  itime=0
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_sample
	    if(ifail.gt.0) then
	    call linesplit(ifail,4)
	    if(ifail.ne.0) go to 9150
	    cols=4
	  end if
	  call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value,iline, &
	  sampfile)
	  if(ifail.ne.0) go to 9900
	  if(value.lt.-1.0e38) cycle read_sample
	  if(daystart.ne.-99999) then
	    if((ndays.lt.ndaystart).or. &
	      ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) &
	      cycle read_sample
	    if((ndays.gt.ndayfin).or. &
	      ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin))) then
		active=.false.
		prevbore=bore_list_id(iactive)
		itime=0
		cycle read_sample
	    end if
	  end if
	  time=((ndays-ndayzero)+(nsecs-nseczero)*secfac)*dayfac
	  itime=1
	  nobs=nobs+1
	  if(nobs.gt.MAX_PM_OBS) go to 9700
	  time_out(nobs)=time
	  value_out(nobs)=value
	  iobs_out(iactive)=nobs
	end do read_sample

1000	call num2char(iline,anum)
	write(amessage,1020) trim(anum),trim(sampfile)
1020    format('  - ',a,' lines read from bore sample file ',a)
	call write_message(leadspace='yes')

	if(nobs.eq.0) then
	  write(amessage,1100) trim(sampfile)
1100	  format(' No observations could be recorded for listed ',&
	  'bores using information contained in bore sample file ',a)
	  go to 9890
	end if

	write(outunit1,1150,err=9550)
1150	format('PMWIN4000_OBS_FILE')
	i=1
	if(atype.eq.'d') i=-1
	write(outunit1,1170,err=9550) nobs,i,0,0,0
1170	format(5(i6))
	numbore=count(iobs_out.ne.-999)
	write(outunit2,1200,err=9600)
1200	format('PMWIN4000_BOR_FILE')
	write(outunit2,1220,err=9600) numbore,0,0,0,0
1220	format(5(i6))
	write(outunit3,1250,err=9650)
1250	format(' bore_identifier   PMWIN_bore_identifier')

	iibore=0
	bore_list: do ibore=1,num_bore_list
	  iobstart=iobs_out(ibore)
	  if(iobstart.eq.-999)cycle bore_list
	  iibore=iibore+1
	  iobsfin=nobs+1
	  inner_bore_list: do j=1,num_bore_list
	    i=iobs_out(j)
	    if(i.le.iobstart)then
	      cycle inner_bore_list
	    else
	      if(i.lt.iobsfin) iobsfin=i
	    end if
	  end do inner_bore_list
	  iobsfin=iobsfin-1
	  do i=iobstart,iobsfin
	    if(atype.eq.'h') then
	      write(outunit1,1300,err=9550) iibore,time_out(i),1.0, &
	      value_out(i),0,0,0,0,0
1300	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,1pg13.6,5(i6))
	    else if(atype.eq.'d') then
	      write(outunit1,1310,err=9550) iibore,time_out(i),1.0,0, &
	      value_out(i),0,0,0,0
1310	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,i6,1x,1pg13.6,4(i6))
	    else if(atype.eq.'c') then
	      write(outunit1,1311,err=9550) iibore,time_out(i),1.0,0, &
	      0,value_out(i),0,0,0
1311	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,2(i6),1x,1pg13.6,3(i6))
	    else if(atype.eq.'m') then
	      write(outunit1,1312,err=9550) iibore,time_out(i),1.0,0, &
	      0,0,value_out(i),0,0
1312	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,3(i6),1x,1pg13.6,2(i6))
	    else if(atype.eq.'p') then
	      write(outunit1,1313,err=9550) iibore,time_out(i),1.0,0, &
	      0,0,0,value_out(i),0
1313	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,4(i6),1x,1pg13.6,1(i6))
	    else if(atype.eq.'s') then
	      write(outunit1,1314,err=9550) iibore,time_out(i),1.0,0, &
	      0,0,0,0,value_out(i)
1314	      format(i6,1x,1pg13.6,1x,0p,f4.1,1x,5(i6),1x,1pg13.6)
	    end if
	  end do
	  do j=1,num_bore_coord
	    if(bore_list_id(ibore).eq.bore_coord_id(j)) exit
	  end do
	  write(outunit2,1350,err=9600) -1,bore_coord_east(j), &
	  bore_coord_north(j),bore_coord_layer(j),-1, &
	  icolour(iibore-((iibore-1)/13)*13)
1350	  format(i6,1x,f15.3,1x,f15.3,2(1x,i5),1x,i12)
	  write(outunit3,1370,err=9650) bore_list_id(ibore),iibore
1370	  format(1x,a10,t20,i6)
	end do bore_list

	write(6,*)
	imessage=0
	do ibore=1,num_bore_list
	  if(iobs_out(ibore).eq.-999) then
	    write(amessage,1610) trim(bore_list_id(ibore))
1610	    format(' Warning: no observation information could be written ',&
	    'for bore ',a,'.')
	    call write_message(increment=1)
	  end if
	end do
	if(imessage.ne.0) write(6,*)

	call num2char(nobs,anum)
        write(amessage,1420) trim(anum),trim(outfile1)
1420    format('  - ',a,' observations written to PMWIN observation file ',a)
        call write_message
        call num2char(iibore,anum)
        write(amessage,1410) trim(anum),trim(outfile2)
1410    format('  - ',a,' bore coordinates written to PMWIN bore cordinates ',&
        'file ',a)
        call write_message
        write(amessage,1450) trim(outfile3)
1450    format('  - User bore_id''s and PMWIN bore_id''s written to file ',a)
        call write_message

	go to 9900

9100    call num2char(iline,anum)
	write(amessage,9110) trim(anum),trim(sampfile)
9110    format(' Error reading line ',a,' of bore sample file ',a)
	go to 9890
9150    call num2char(iline,anum)
	write(amessage,9160) trim(anum),trim(sampfile)
9160    format(' Error on line ',a,' of bore sample file ',a,': ',&
	'insufficient entries.')
	go to 9890
9200    call num2char(iline,anum)
	write(amessage,9210) trim(anum),trim(sampfile)
9210    format(' Error on line ',a,' of bore sample file ',a,': bore ',&
	'identifier greater than 10 characters in length.')
	go to 9890
9300    write(amessage,9310)
9310    format(' Cannot allocate sufficient memory to run program SMP2PM1.')
	go to 9890
9350	write(amessage,9360)
9360	format(' File handling error: cannot continue execution.')
	go to 9890
9500	write(amessage,9510)
9510	format(' File handling error: cannot continue execution.')
	go to 9890
9550	write(amessage,9560) trim(outfile1)
9560	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	go to 9890
9600	write(amessage,9560) trim(outfile2)
	go to 9890
9650	write(amessage,9560) trim(outfile3)
	go to 9890
9700	call num2char(MAX_PM_OBS,anum)
	write(amessage,9710) trim(anum)
9710	format(' There are more than ',a,' observations pertaining to the ', &
	'bores listed in the bore listing file. This exceeds the ', &
 	'limit for a PMWIN observation file.')
	go to 9890

9890    call write_message(leadspace='yes')
9900    call close_files
	deallocate(iobs_out,stat=ierr)
	write(6,*)

end program smp2pm1

 