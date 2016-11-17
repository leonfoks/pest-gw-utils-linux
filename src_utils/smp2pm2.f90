!     Last change:  JD   18 Dec 2000    0:31 am
program smp2pm2

! -- Program SMP2PM2 interpolates data contained in a bore sample file to a
!    set of times specified in a MODBORE/MT3D output file, creating a PMWIN
!    observation file.

	use defn
	use inter

	implicit none

	integer :: sampunit,modunit,daystart,monstart,yearstart,hourstart, &
		   minstart,secstart,ntime,i,j,itemp,intdaystart, &
		   intsecstart,itime,icollect,iline,cols,ibore,ilist,iobs, &
		   ierr,ifail,outunit1,outunit2,outunit3,iibore,newbore,idate,iheader
	integer, dimension(13)  :: icolour
	real	:: timefac,rnear,rconst,rtemp
	double precision   :: dtemp
	character (len=80) :: specfile,coordfile,sampfile,modfile,aprompt, &
			      outfile1,outfile2,outfile3
	character (len=15) :: atemp,abore,aoldbore,aline
	character (len=4)  :: aobs
	character (len=1)  :: aunit,atype

	integer, parameter               	       :: NUMTIME = 1000
	integer, allocatable, dimension(:)	       :: intday,intsec
	character (len=15), allocatable, dimension(:)  :: atime
	real, allocatable, dimension(:)                :: time

	integer, parameter			       :: NUM_SAMP_BORE=3000
	integer, allocatable, dimension(:)	       :: ndays,nsecs
	double precision, allocatable, dimension(:)    :: value

	real, allocatable, dimension(:,:)              :: valinterp

	integer, parameter		:: MAX_PM_BORE=1000
	integer, parameter		:: MAX_PM_OBS=6000

	data icolour /255,16711680,65535,16776960,65280,8388736,12632256, &
        16711935,32896,8388608,4227072,128,0/ 


	write(amessage,5)
5       format(' Program SMP2PM2 interpolates data contained in a bore ', &
	'sample file to a set of times specified in a MODBORE/MT3D output ',&
	'file, creating a PMWIN observation file.')
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

	allocate(atime(NUMTIME), time(NUMTIME), ndays(NUM_SAMP_BORE), &
	nsecs(NUM_SAMP_BORE), value(NUM_SAMP_BORE), stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,10)
10	  format(' Insufficient memory to continue SMP2PM2 execution. ', &
	  'Edit SMP2PM2 source code, reduce value of variables NUMTIME ', &
	  'and NUM_SAMP_BORE; then recompile.')
	  go to 9890
	end if

	call readfig(specfile,bore_coord_file,sampfile)

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
35	call read_bore_coord_file(ifail, &
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

49	write(6,*)
50      call read_bore_list_file(ifail, &
        ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 35
	end if

	if(num_bore_list.gt.MAX_PM_BORE) then
	  call num2char(MAX_PM_BORE,aline)
	  write(amessage,55) trim(aline),trim(bore_list_file)
55	  format(' PMWIN can read a maximum of ',a,' bores: there are more ', &
	  'than this in bore listing file ',a)
	  go to 9890
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
		write(amessage,60) trim(bore_list_file)
60              format(' Execution of program SMP2PM2 cannot proceed as ',&
		'there are multiple occurrences of the same bore in bore ',&
		'listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	write(6,*)
70	call open_input_file(ifail,&
	' Enter name of MODBORE/MT3BORE-generated file: ',modfile,modunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 50
	end if

	call get_modbore_times(ifail,NUMTIME,ntime,time,atime,modunit,modfile)
	if(ifail.ne.0) go to 9900
	close(unit=modunit,err=9300)

	write(6,*)
100     write(6,110)
110     format(' Provide the following model timing data ---->')
130     write(6,150,advance='no')
150     format('   Enter time units used by model (yr/day/hr/min/sec) ', &
	'[y/d/h/m/s]: ')
	read(5,'(a)') aunit
	if(aunit.eq.' ') go to 130
	call casetrans(aunit,'lo')
	if(index(eschar,aunit).ne.0) then
	  write(6,*)
	  go to 70
	end if
	if(aunit.eq.'y')then
	  timefac=365.25
	else if(aunit.eq.'d')then
	  timefac=1.0
	else if(aunit.eq.'h')then
	  timefac=1.0/24.0
	else if(aunit.eq.'m')then
	  timefac=1.0/1440.0
	else if(aunit.eq.'s')then
	  timefac=1.0/86400.0
	else
	  go to 130
	end if

200     if(datespec.eq.1) then
	  write(6,210,advance='no')
210       format('   Enter simulation starting date [dd/mm/yyyy]: ')
	else
	  write(6,211,advance='no')
211       format('   Enter simulation starting date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 200
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 130
	end if
	call char2date(ifail,atemp,daystart,monstart,yearstart)
	if(ifail.ne.0) then
	  write(6,240)
240       format('   Improper date  - try again.')
	  go to 200
	end if

250     write(6,260,advance='no')
260     format('   Enter simulation starting time   [hh:mm:ss]: ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 250
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 200
	end if
	call char2time(ifail,atemp,hourstart,minstart,secstart)
	if(ifail.ne.0) then
	  write(6,280)
280       format('   Improper time  - try again.')
	  go to 250
	end if

300     write(6,*)
310     write(6,320)
320     format(' Provide the following time-interpolation parameters ---->')
350     write(6,360,advance='no')
360     format('   Enter maximum days to nearest sample (fractional if ', &
	'necessary): ')
	itemp=key_read(rnear)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 250
	end if
	if(itemp.lt.0)then
	  go to 350
	else if(itemp.gt.0) then
	  write(6,370)
370       format('   Data input error  - try again.')
	  go to 350
	end if
	if(rnear.le.0.0) then
	  write(6,375)
375	  format('   Error: number of days must be greater than zero',&
	  '  - try again.')
	  go to 350
	end if

380     write(6,400)
400     format('   Enter days over which sample can be assumed constant if ',&
	'linear interpolation')
	write(6,410,advance='no')
410     format('   cannot take place (fractional if necessary): ')
	itemp=key_read(rconst)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 350
	end if
	if(itemp.lt.0)then
	  go to 380
	else if(itemp.gt.0) then
	  write(6,370)
	  go to 380
	end if
	if(rconst.lt.0) then
	  write(6,420)
420	  format('   Error: number of days must not be negative  - try again.')
	  go to 380
	end if
	if(rconst.gt.rnear) then
	  write(amessage,430)
430       format('   Error: this number cannot exceed days to nearest sample', &
	  '  - try again.')
	  call write_message
	  go to 380
	end if

500     write(6,*)
	write(6,505)
505	format(' Provide the name for the following SMP2PM2 output files:-')
	aprompt=' Enter name for output PMWIN observation file: '
510	call open_output_file(ifail,aprompt,outfile1,outunit1)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 380
	end if  

530	aprompt=' Enter name for output PMWIN bore coordinates file: '
	call open_output_file(ifail,aprompt,outfile2,outunit2)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=outunit1,err=9300)
	  go to 500
	end if

550	aprompt=' Enter name for output bore_id to PMWIN_id file: '
	call open_output_file(ifail,aprompt,outfile3,outunit3)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=outunit2,err=9300)
	  write(6,*)
	  go to 530
	end if

	allocate(intday(ntime),intsec(ntime),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,10)
	  go to 9890
	end if
	intdaystart=numdays(1,1,1970,daystart,monstart,yearstart)
	intsecstart=numsecs(0,0,0,hourstart,minstart,secstart)
	do itime=1,ntime
	  time(itime)=time(itime)*timefac
	  itemp=floor(time(itime))
	  intday(itime)=intdaystart+itemp
	  rtemp=time(itime)-itemp
	  rtemp=rtemp*86400.0
	  itemp=nint(rtemp)
	  intsec(itime)=intsecstart+itemp
	  if(intsec(itime).ge.86400) then
	    intsec(itime)=intsec(itime)-86400
	    intday(itime)=intday(itime)+1
	  end if
	end do

	deallocate(time,stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,525)
525	  format(' Memory management error: cannot continue SMP2PM2 execution.')
	  go to 9890
	end if

	allocate(valinterp(num_bore_list,ntime), stat=ierr)
	if(ierr.ne.0) then
	   write(amessage,520)
520        format(' Insufficient memory to continue SMP2PM2 execution. ',&
	   'Use fewer bores and/or model output times. It may also help ', &
	   'to edit SMP2PM2 source code, lower parameters NUMTIME and ', &
	   'NUM_SAMP_BORE as appropriate, and recompile.')
	   go to 9890
	end if
	valinterp=-6.1e37

	write(6,*)
	icollect=0
	iline=0
	abore=' '
	read_line: do
	  aoldbore=abore
	  iline=iline+1
	  if(iline.eq.2000) write(6,650)
650       format(' Working.....')
	  read(sampunit,'(a)',err=9000,end=1000) cline
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_line
	  if(ifail.gt.0)then
	    cols=4
	    call linesplit(ifail,4)
	    if(ifail.ne.0)then
	      call num2char(iline,aline)
	      write(amessage,680) trim(aline),trim(sampfile)
680           format('insufficient items on line ',a,' of bore sample file ',a)
	      call write_message(error='yes')
	      go to 9900
	    end if
	  end if

	  if(right_word(1)-left_word(1).gt.9) then
	    call num2char(iline,aline)
	    write(amessage,700) trim(aline),trim(sampfile)
700         format('bore identifier greater than 10 characters in length ',&
	    'on line ',a,' of bore sample file ',a)
	    call write_message(error='yes')
	    go to 9900
	  end if
	  abore=cline(left_word(1):right_word(1))
	  call casetrans(abore,'hi')
	  if(abore.eq.aoldbore)then
	    if(icollect.eq.1)then
	      ibore=ibore+1
	      if(ibore.gt.NUM_SAMP_BORE)then
		call num2char(NUM_SAMP_BORE,aline)
		write(amessage,720) trim(sampfile),trim(abore),trim(aline)
720             format(' There are more samples in file ',a,' for bore ',&
		a,' than the maximum allowed number of ',a,&
		': edit source code file and assign a greater value to ',&
		'parameter NUM_SAMP_BORE.  Then recompile program.')
	        call write_message
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
	      if((valinterp(ilist,1).lt.-6.2e37).or.&
		 (valinterp(ilist,1).gt.-6.05e37))then
		 write(amessage,740) trim(bore_list_id(ilist)),&
		 trim(sampfile)
740              format('entries for bore ',a,' in bore sample file ',&
		 a,' are not all in sequence.')
		 call write_message(error='yes')
		 go to 9900
	      end if
	      do itime=1,ntime
	        call time_interp(ifail,ibore,ndays,nsecs,value,intday(itime), &
	        intsec(itime),rnear,rconst,dtemp)
		valinterp(ilist,itime)=dtemp
	        if(ifail.eq.1)then
		  write(amessage,750) trim(bore_list_id(ilist)),trim(sampfile)
750               format('non-increasing sample times for bore ',a,&
		  ' in bore sample file ',a)
		 call write_message(error='yes')
		  go to 9900
	        end if
	      end do
	      icollect=0
	    end if
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) then
		ilist=i
		icollect=1
		ibore=1
		call read_rest_of_sample_line(ifail,cols,ndays(1),nsecs(1), &
		value(1),iline,sampfile)
		if(ifail.ne.0) go to 9900
		go to 800
	      end if
	    end do
	    icollect=0
800         continue
	  end if
	end do read_line

1000    if(icollect.eq.1)then
	  if((valinterp(ilist,1).lt.-6.2e37).or.&
	  (valinterp(ilist,1).gt.-6.05e37))then
	    write(amessage,740) trim(bore_list_id(ilist)),&
	    trim(sampfile)
	    call write_message(error='yes')
	    go to 9900
	  end if
	  do itime=1,ntime
	    call time_interp(ifail,ibore,ndays,nsecs,value,intday(itime), &
	    intsec(itime),rnear,rconst,dtemp)
	    valinterp(ilist,itime)=dtemp
	    if(ifail.eq.1)then
	      write(amessage,750) trim(bore_list_id(ilist)),trim(sampfile)
	      call write_message(error='yes')
	      go to 9900
	    end if
	  end do
	end if

! -- The number of observations is next determined for the file header.

	iobs=0
	do ibore=1,num_bore_list
	  do itime=1,ntime
	    if(valinterp(ibore,itime).lt.-1.0e37) cycle
	    iobs=iobs+1
	  end do
	end do
	if(iobs.gt.MAX_PM_OBS) then
	  call num2char(MAX_PM_OBS,aline)
	  write(amessage,1010) trim(aline)
1010	  format(' There are more than ',a,' observations. This exceeds ',&
	  'the limit for a PMWIN observation file.')
	  call write_message
	  go to 9900
	end if
	write(outunit1,1020,err=9200)
1020	format('PMWIN4000_OBS_FILE')
	write(outunit1,1025,err=9200) iobs,-1,0,0,0
1025	format(5(i6))

	iibore=0
	iobs=0
	do ibore=1,num_bore_list
	  newbore=1
	  do itime=1,ntime
	    if(valinterp(ibore,itime).lt.-1.0e37) cycle
	    if(newbore.eq.1) then
	      iibore=iibore+1
	      newbore=0
	    end if
	    iobs=iobs+1
	    if(atype.eq.'h') then
	      write(outunit1,1100,err=9200) iibore,atime(itime),1.0, &
	      valinterp(ibore,itime),0,0,0,0,0
1100	      format(i6,1x,a15,1x,f4.1,1x,1pg13.6,5(i6))
	    else if(atype.eq.'d') then
	      write(outunit1,1101,err=9200) iibore,atime(itime),1.0,0, &
	      valinterp(ibore,itime),0,0,0,0
1101	      format(i6,1x,a15,1x,f4.1,1x,i6,1x,1pg13.6,4(i6))
	    else if(atype.eq.'c') then
	      write(outunit1,1102,err=9200) iibore,atime(itime),1.0,0, &
	      0,valinterp(ibore,itime),0,0,0
1102	      format(i6,1x,a15,1x,f4.1,1x,2(i6),1x,1pg13.6,3(i6))
	    else if(atype.eq.'m') then
	      write(outunit1,1103,err=9200) iibore,atime(itime),1.0,0, &
	      0,0,valinterp(ibore,itime),0,0
1103	      format(i6,1x,a15,1x,f4.1,1x,3(i6),1x,1pg13.6,2(i6))
	    else if(atype.eq.'p') then
	      write(outunit1,1104,err=9200) iibore,atime(itime),1.0,0, &
	      0,0,0,valinterp(ibore,itime),0
1104	      format(i6,1x,a15,1x,f4.1,1x,4(i6),1x,1pg13.6,1(i6))
	    else if(atype.eq.'s') then
	      write(outunit1,1105,err=9200) iibore,atime(itime),1.0,0, &
	      0,0,0,0,valinterp(ibore,itime)
1105	      format(i6,1x,a15,1x,f4.1,1x,5(i6),1x,1pg13.6)
	    end if
	  end do
	end do

	write(6,*)
	imessage=0
	write(outunit2,1150,err=9400)
1150	format('PMWIN4000_BOR_FILE')
	write(outunit2,1160,err=9400) iibore,0,0,0,0
1160	format(5(i6))
	write(outunit3,1165,err=9500)
1165	format(' bore_identifier   PMWIN_bore_identifier')
	iibore=0
	borelist: do ibore=1,num_bore_list
	  do itime=1,ntime
	    if(valinterp(ibore,itime).gt.-1.0e37) go to 1180
	  end do
	  write(amessage,1240) trim(bore_list_id(ibore))
1240	  format(' Warning: no data for bore ',a,' could be interpolated ',&
	  'to model times.')
	  call write_message(increment=1)
	  cycle borelist
1180	  iibore=iibore+1
	  do j=1,num_bore_coord
	    if(bore_list_id(ibore).eq.bore_coord_id(j)) exit
	  end do
	  write(outunit2,1170,err=9400) -1,bore_coord_east(j), bore_coord_north(j), &
	  bore_coord_layer(j),-1,icolour(iibore-((iibore-1)/13)*13)
1170	  format(i6,1x,f15.3,1x,f15.3,2(1x,i5),1x,i12)
	  write(outunit3,1190,err=9500) bore_list_id(ibore),iibore
1190	  format(1x,a10,t20,i6)
	end do borelist

	if(imessage.ne.0) write(6,*)
1200	call num2char(iobs,aline)
	write(amessage,1220) trim(aline),trim(outfile1)
1220	format('  - ',a,' observations written to PMWIN observation file ',a)
	call write_message
	call num2char(iibore,aline)
	write(amessage,1210) trim(aline),trim(outfile2)
1210	format('  - ',a,' bore coordinates written to PMWIN bore cordinates ',&
        'file ',a)
	call write_message
	write(amessage,1250) trim(outfile3)
1250	format('  - User bore_id''s and PMWIN bore_id''s written to file ',a)
	call write_message

	go to 9900

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(sampfile)
9010    format(' Error reading line ',a,' of bore sample file ',a)
	go to 9890
9200	write(amessage,9210) trim(outfile1)
9210	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	go to 9890
9300	write(amessage,9310)
9310	format(' File management error: cannot continue execution.')
	go to 9890
9400	write(amessage,9210) trim(outfile2)
	go to 9890
9500	write(amessage,9210) trim(outfile3)
	go to 9890

9890	call write_message(leadspace='yes')
9900    call close_files
	call free_bore_mem
	deallocate(valinterp,intday,intsec,atime,time,ndays,nsecs,value, &
	stat=ierr)
	write(6,*)

end program smp2pm2


subroutine get_modbore_times(ifail,NUMTIME,ntime,time,atime,modunit,modfile)

! -- Subroutine get_modbore_times reads a MODBORE/MT3BORE output file to
!    ascertain MODFLOW/MT3D unformated output times.

! -- Arguments are as follows:-
!       ifail:    returned as non-zero if error condition encountered
!       numtime:  maximum allowed number of output times
!       ntime:    number of MODFLOW/MT3D output times found on MODBORE/MT3BORE
!                 output file
!       time:     MODFLOW/MT3D output times
!       atime:    character string array recording times on MODBORE/MT3BORE
!                 output file
!       modunit:  unit number of MODBORE/MT3BORE output file
!       modfile:  name of MODBORE/MT3BORE output file

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter
	implicit none

	integer, intent(out)                             :: ifail
	integer, intent(in)                              :: NUMTIME
	integer, intent(out)            	         :: ntime
	real, intent(out)		                 :: time(NUMTIME)
	character (len=*), intent(out)                   :: atime(NUMTIME)
	integer, intent(in)                              :: modunit
	character (len=*), intent(in)                    :: modfile

	integer                         :: iline,imark
	character (len=15)              :: anum1,anum2
	character (len=20)		:: atemp
	
	iline=0
	ntime=0
	read_a_line: do
	  iline=iline+1
	  read(modunit,'(a)',err=9000,end=1000) cline
	  if(index(cline,'SIMUL').eq.0) cycle read_a_line
	  ntime=ntime+1
	  if(ntime.gt.NUMTIME) cycle read_a_line
	  imark=index(cline,'=')
	  atemp=cline(imark+1:len_trim(cline))
	  atime(ntime)=adjustl(atemp)
	  call char2num(ifail,atime(ntime),time(ntime))
	  if(ifail.ne.0) go to 9200
	end do read_a_line

1000    if(ntime.gt.NUMTIME) go to 9050
	if(ntime.eq.0) then
	  write(amessage,1050) trim(modfile)
1050	  format(' No model output data has been detected in MODBORE/MT3D-', &
	  'generated file ',a)
	  go to 9890
	end if
	call num2char(iline-1,anum1)
	write(amessage,1010) trim(anum1), trim(modfile)
1010    format(' - ',a,' lines read from MODBORE/MT3BORE output file ',a)
	call write_message
	call num2char(ntime,anum2)
	write(amessage,1020) trim(anum2)
1020	format(' - model output data detected for ',a,' simulation times.')
	call write_message
	go to 9900

9000    call num2char(iline,anum1)
	write(amessage,9010) trim(anum1),trim(modfile)
9010    format(' Error reading line ',a,' of MODBORE/MT3BORE-generated ', &
	'file ',a)
	go to 9890
9050    call num2char(ntime,anum1)
	call num2char(NUMTIME,anum2)
	write(amessage,9060) trim(anum1),trim(modfile),trim(anum2)
9060    format(' A total of ',a,' model output times have been detected ', &
	'in MODBORE/MT3D output file ',a,'. However program SMP2PM2 can ',&
	'handle only ',a,' model output times. Edit SMP2PM2 source code, ', &
	'alter parameter NUMTIME, and recompile.')
	go to 9890
9200    call num2char(iline,anum1)
	write(amessage,9210) trim(anum1),trim(modfile)
9210    format(' Error reading elapsed model simulation time from line ',a, &
	' of MODBORE/MT3BORE output file ',a)
	go to 9890

9890    call write_message(leadspace='yes')
	ifail=1
9900    return

end subroutine get_modbore_times

 