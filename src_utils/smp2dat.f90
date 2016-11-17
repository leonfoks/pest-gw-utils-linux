!     Last change:  JD   17 Dec 2000   11:13 pm
program smp2dat

! -- Program SMP2DAT interpolates data contained in a bore sample file to a
!    set of times specified in a MODBORE/MT3D output file, creating a bore
!    data file for use with the PEST MODFLOW/MT3D Utilities.

	use defn
	use inter

	implicit none

	integer :: sampunit,modunit,daystart,monstart,yearstart,hourstart, &
		   minstart,secstart,ntime,i,j,itemp,outunit,intdaystart, &
		   intsecstart,itime,icollect,iline,cols,ibore,ilist,iobs, &
		   ierr,ifail,idate,nn,ii1,ii2,jj1,jj2,nnchar,nnnchar,iheader
	real	:: timefac,rnear,rconst,rtemp
	double precision   :: dtemp
	character (len=80) :: specfile,coordfile,sampfile,modfile,aprompt, &
			      outfile
	character (len=15) :: atemp,abore,aoldbore,aline,afmt
	character (len=8)  :: aachar
	character (len=8)  :: aobs
	character (len=1)  :: aunit,ause,aid

	integer, parameter               	       :: NUMTIME = 1000
	integer, allocatable, dimension(:)	       :: intday,intsec
	character (len=15), allocatable, dimension(:)  :: atime
	real, allocatable, dimension(:)                :: time

	integer, parameter			       :: NUM_SAMP_BORE=3000
	integer, allocatable, dimension(:)	       :: ndays,nsecs
	double precision, allocatable, dimension(:)    :: value

	real, allocatable, dimension(:,:)              :: valinterp


	write(amessage,5)
5       format(' Program SMP2DAT interpolates data contained in a bore ', &
	'sample file to a set of times specified in a MODBORE/MT3D output ',&
	'file, creating a bore data file for use with the PEST MODFLOW/MT3D', &
	' Utilities.')
	call write_message(leadspace='yes',endspace='yes')

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

	allocate(atime(NUMTIME), time(NUMTIME), ndays(NUM_SAMP_BORE), &
	nsecs(NUM_SAMP_BORE), value(NUM_SAMP_BORE), stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,10)
10	  format(' Insufficient memory to continue SMP2DAT execution. ', &
	  'Edit SMP2DAT source code, reduce value of variables NUMTIME ', &
	  'and NUM_SAMP_BORE; then recompile.')
	  go to 9890
	end if

	call readfig(specfile,bore_coord_file,sampfile)

30      call open_named_input_file(ifail, &
	' Enter name of bore sample file: ',sampfile,sampunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

	write(6,*)
50      call read_bore_list_file(ifail, &
        ' Enter name of bore listing file: ',coord_check='no')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  close(unit=sampunit,err=9300)
	  go to 30
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
		write(amessage,60) trim(bore_list_file)
60              format(' Execution of program SMP2DAT cannot proceed as ',&
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
	aprompt=' Enter name for output bore data file: '
510	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 380
	end if  

550	write(6,560,advance='no')
560	format(' Use numbers or bore identifiers for observation ', &
	'names?  [n/b]: ')
	read(5,'(a)') ause
	if(ause.eq.' ') go to 550
	if(index(eschar,ause).ne.0) then
	  close(unit=outunit,err=9300)
	  go to 500
	end if
	call casetrans(ause,'lo')
	if((ause.ne.'n').and.(ause.ne.'b')) go to 550
	if(ause.eq.'b')then
	  if(ntime.gt.99999)then
	    write(amessage,563)
563	    format(' There are too many model output times for ',   &
            'this option.')
	    call write_message
	    go to 550
	  end if
	  if(ntime.gt.10000)then
	    nnnchar=5
	  else if(ntime.ge.1000)then
	    nnnchar=4
	  else if(ntime.ge.100)then
	    nnnchar=3
	  else if(ntime.ge.10)then
	    nnnchar=2
	  else
	    nnnchar=1
	  end if
	  nnchar=8-nnnchar-1
	  write(aachar,'(i8)') nnchar
	  aachar=adjustl(aachar)
564	  write(6,565,advance='no') trim(aachar),trim(aachar)
565	  format(' Use first ',a,' or last ',a,' characters ',      &
          'of bore identifier?  [f/l]: ')
	  read(5,'(a)') aid
	  if(aid.eq.' ') go to 564
	  if(index(eschar,aid).ne.0) go to 550
	  call casetrans(aid,'lo')
	  if((aid.ne.'f').and.(aid.ne.'l')) go to 564
	  if(aid.eq.'f') then
	    do i=1,num_bore_list-1
	      do j=i+1,num_bore_list
	        if(bore_list_id(i)(1:nnchar).eq.                    &
                   bore_list_id(j)(1:nnchar)) then
	           write(amessage,570)
570	           format(' This will not result in unique observation ', &
	           'names.')
	           call write_message(leadspace='yes',endspace='yes')
	           go to 550
	        end if
	      end do
	    end do
	  else
	    do i=1,num_bore_list-1
	      do j=i+1,num_bore_list
	        jj1=len_trim(bore_list_id(i))
	        ii1=max(1,jj1-nnchar+1)
	        jj2=len_trim(bore_list_id(j))
	        ii2=max(1,jj2-nnchar+1)
	        if(bore_list_id(i)(ii1:jj1).eq.                     &
                   bore_list_id(j)(ii2:jj2))then
	           write(amessage,570)
	           call write_message(leadspace='yes',endspace='yes')
	           go to 550
	        end if
	      end do
	    end do
	  end if
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
	  write(amessage,530)
530	  format(' Memory management error: cannot continue SMP2DAT execution.')
	  go to 9890
	end if

	allocate(valinterp(num_bore_list,ntime), stat=ierr)
	if(ierr.ne.0) then
	   write(amessage,520)
520        format(' Insufficient memory to continue SMP2DAT execution. ',&
	   'Use fewer bores and/or model output times. It may also help ', &
	   'to edit SMP2DAT source code, lower parameters NUMTIME and ', &
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

	iobs=0
	do ibore=1,num_bore_list
	  do itime=1,ntime
	    if(valinterp(ibore,itime).lt.-1.0e37) cycle
	    iobs=iobs+1
	    if(ause.eq.'n') then
	      if(iobs.gt.99999999) then
	        write(amessage,1050) trim(outfile)
1050	        format(' Cannot write entire bore data file ',a, &
	        ': maximum of 9999999 observations can be written to this file.')
	        call write_message(endspace='yes')
	        iobs=99999999
	        go to 1200
	      end if
	      call num2char(iobs,aobs)
	    else
	      if(aid.eq.'f')then
	        ii1=1
	        jj1=min(nnchar,len_trim(bore_list_id(ibore)))
	      else
	        jj1=len_trim(bore_list_id(ibore))
	        ii1=max(1,jj1-nnchar+1)
	      end if
	      write(afmt,1060) nnnchar,nnnchar
1060	      format('(i',i5,'.',i5,')')
	      write(atemp,afmt) itime
	      atemp=adjustl(atemp)
	      aobs=bore_list_id(ibore)(ii1:jj1)//'_'//trim(atemp)
	    end if
	    write(outunit,1100,err=9200) trim(bore_list_id(ibore)), &
	    trim(atime(itime)),valinterp(ibore,itime),'1.0',trim(aobs)
1100	    format(1x,a,t15,a,t31,1pg14.6,t47,a,t53,a)
	  end do
	end do

1200	call num2char(iobs,aline)
	write(amessage,1220) trim(aline),trim(outfile)
1220	format(' - ',a,' observations written to bore data file ',a)
	call write_message
	imessage=0
	do ibore=1,num_bore_list
	  if(imessage.gt.40) go to 9900
	  if((valinterp(ibore,1).gt.-6.15e37).and. &
	     (valinterp(ibore,1).lt.-6.05e37)) then
	     if(imessage.eq.0) write(6,*)
	     write(amessage,1240) trim(bore_list_id(ibore)),trim(sampfile)
1240	     format(' Warning: no data for bore ',a, &
	     ' found in bore sample file ',a)
	     call write_message(increment=1)
	  end if
	end do
	go to 9900

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(sampfile)
9010    format(' Error reading line ',a,' of bore sample file ',a)
	go to 9890
9200	write(amessage,9210) trim(outfile)
9210	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	go to 9890
9300	write(amessage,9310)
9310	format(' File management error: cannot continue execution.')
	go to 9890

9890	call write_message(leadspace='yes')
9900    call close_files
	call free_bore_mem
	deallocate(valinterp,intday,intsec,atime,time,ndays,nsecs,value, &
	stat=ierr)
	write(6,*)

end program smp2dat


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
	'in MODBORE/MT3D output file ',a,'. However program SMP2DAT can ',&
	'handle only ',a,' model output times. Edit SMP2DAT source code, ', &
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

 