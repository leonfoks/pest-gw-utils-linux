!     Last change:  JD   13 Feb 2001    5:50 pm
program smpcal

! -- Program SMPCAL calibrates the contents of one bore sample file against the
!    contents of another.

	use defn
	use inter

	implicit none

	logical :: lexist
	integer :: ifail,idate,sunit1,sunit2,ierr,nid1,maxdatsmp,i,calunit, &
                   repunit,j,cols,icollect,iline,ndays,nsecs,istand,nstand, &
                   esecs,nsamp,iid,iend,iday,imon,iyear,ihour,imin,isec,    &
                   iloc1,iloc2,iloc,iday1,imon1,iyear1,ihour1,imin1,    &
                   isec1,iday2,imon2,iyear2,ihour2,imin2,isec2,begday,      &
                   begsec,finday,finsec,k,iheader
	real                               :: ehours
	double precision                   :: value,bbval1,bbval2,ssval1,   &
                                              ssval2,m,c
	integer, allocatable               :: ndays1(:),nsecs1(:),ndays2(:), &
                                              nsecs2(:),loc(:)
	integer, allocatable               :: sdays(:),ssecs(:),bdays(:), &
                                              bsecs(:)
	double precision, allocatable      :: sval(:),bval(:)
	character (len=10)                 :: abore,aboreold,aabore,aaboreold
	character (len=15)                 :: aline,atemp
	character (len=10), allocatable    :: aid1(:)
	character (len=80)                 :: specfile,sfile1,sfile2,calfile, &
                                              repfile1,aprompt


! -- Preparitory stuff.

	calfile=' '
	repfile1=' '

        write(amessage,5)
5       format(' Program SMPCAL calibrates the contents of one bore sample ', &
        'file against those of another.')
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

	call readfig(specfile,sampfile=sfile2)

! -- Obtain user input.

50	call open_input_file(ifail, &
	' Enter name of bore sample file requiring calibration: ',  &
	sfile1,sunit1)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

150	call open_named_input_file(ifail, &
        ' Enter name of standard bore sample file: ',sfile2,sunit2)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
	  close(unit=sunit1,iostat=ierr)
	  write(6,*)
          go to 50
        end if

	write(6,*)
180	write(6,190,advance='no')
190	format(' Enter maximum extrapolation time in hours: ')
	ifail=key_read(ehours)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  close(unit=sunit2)
	  go to 150
	end if
        if(ifail.lt.0)then
          go to 180
        else if(ifail.gt.0) then
          write(6,195)
195       format(' Data input error  - try again.')
          go to 180
        end if
        if(ehours.lt.0.0) then
          write(6,200)
200       format(' Error: number of hours must not be negative',&
          '  - try again.')
          go to 180
        end if
	esecs=ehours*3600

	write(6,*)
220	aprompt=' Enter name for calibrated bore sample file: '
	call open_output_file(ifail,aprompt,calfile,calunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 180
	end if

240	aprompt=' Enter name for report file: '
	call open_output_file(ifail,aprompt,repfile1,repunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  close(unit=calunit)
	  write(6,*)
	  go to 220
	end if


! -- The logger sample file is read to ascertain how many different boreid's
!    are cited in it. The maximum number of samples for any bore is also read.

	write(6,241)
241	format(/,' Working....')

	call get_num_ids(ifail,sunit1,sfile1,nid1,maxdatsmp,ignore_x='yes')
	if(ifail.ne.0) go to 9900
	if(nid1.eq.0) then
	  write(amessage,270) trim(sfile1)
270	  format('no measurements were found in bore sample file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9900
	end if

! -- The boreids are now read from the logger data file.

	allocate(aid1(nid1),ndays1(nid1),nsecs1(nid1),ndays2(nid1), &
                 nsecs2(nid1),loc(nid1),stat=ierr)
	loc=0

	if(ierr.ne.0) then
	  write(amessage,280)
280	  format(' Cannot allocate sufficient memory to run program.')
	  call write_message(leadspace='yes')
          go to 9900
        end if
	call get_ids_and_interval(ifail,sunit1,sfile1,nid1,aid1,ndays1,nsecs1, &
                      ndays2,nsecs2,ignore_x='yes')
	if(ifail.ne.0) go to 9900
	if(nid1.ne.1) then
	  do i=1,nid1-1
	    do j=i+1,nid1
	      if(aid1(i).eq.aid1(j)) then
	        write(amessage,300) trim(aid1(i)),trim(sfile1)
300	        format('all samples for bore ',a,' are not together ', &
                       'in file ',a)
	        call write_message(error='yes',leadspace='yes')
	        go to 9900
	      end if
	    end do
	  end do
	end if

	do i=1,nid1
	  nsecs1(i)=nsecs1(i)-esecs
	  do while(nsecs1(i).lt.0)
	    ndays1(i)=ndays1(i)-1
	    nsecs1(i)=nsecs1(i)+86400
	  end do
	  nsecs2(i)=nsecs2(i)+esecs
	  do while(nsecs2(i).ge.86400)
	    ndays2(i)=ndays2(i)+1
	    nsecs2(i)=nsecs2(i)-86400
	  end do
	end do

! -- Next the standard bore sample file is read into memory, but only for 
!    those bores and those times for which there is data to be calibrated.

	istand=0
	icollect=0
	iline=0
	aboreold=' '
	read_samp: do
	  iline=iline+1
          read(sunit2,'(a)',err=9000,end=500) cline
          cols=5
          call linesplit(ifail,5)
          if(ifail.lt.0) cycle read_samp
	  if(ifail.eq.0) then
	    call casetrans(cline(left_word(5):right_word(5)),'lo')
	    if(cline(left_word(5):right_word(5)).eq.'x') cycle read_samp
	  else
            cols=4
            call linesplit(ifail,4)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,320) trim(aline),trim(sfile2)
320           format('insufficient items on line ',a,' of bore sample file ',a)
              call write_message(error='yes',leadspace='yes')
              go to 9900
            end if
          end if
	  if(right_word(1)-left_word(1).gt.9) then
            call num2char(iline,aline)
            write(amessage,330) trim(aline),trim(sfile2)
330         format('bore identifier greater than 10 characters in length ',&
            'on line ',a,' of bore sample file ',a)
            call write_message(error='yes',leadspace='yes')
            go to 9900
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if((icollect.eq.0).and.(abore.eq.aboreold)) cycle read_samp
	  if(abore.ne.aboreold)then
	    do i=1,nid1
	      if(abore.eq.aid1(i)) go to 360
	    end do
	    aboreold=abore
	    icollect=0
	    cycle read_samp
360	    aboreold=abore
	    icollect=i
	  end if
	  if(icollect.eq.0) cycle read_samp
	  call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value, &
          iline,sfile2)
	  if(ifail.ne.0) go to 9900
	  if(((ndays.gt.ndays1(icollect)).or.                               &
             ((ndays.eq.ndays1(icollect)).and.(nsecs.ge.nsecs1(icollect)))) &
             .and.                                                          &
             ((ndays.lt.ndays2(icollect)).or.                               &
             ((ndays.eq.ndays2(icollect)).and.(nsecs.le.nsecs2(icollect)))))&
             istand=istand+1
	end do read_samp	

500	continue
	rewind(unit=sunit2)
	nstand=istand
	if(nstand.eq.0) then
	  write(amessage,510) trim(sfile1), trim(sfile2)
510	  format('either there are no bores from file ',a,' cited in ', &
          'standard sample file ',a,' or no samples within the latter file ', &
          'are within the calibration time window for the former file.')
	  call write_message(error='yes', leadspace='yes')
	  go to 9900
	end if

	allocate(sdays(nstand),ssecs(nstand),sval(nstand),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,280)
	  call write_message(leadspace='yes')
          go to 9900
        end if
	icollect=0
	istand=0
	iline=0
	aboreold=' '
	aaboreold=' '
	read_samp_1: do
	  iline=iline+1
          read(sunit2,'(a)',err=9000,end=600) cline
          cols=5
          call linesplit(ifail,5)
          if(ifail.lt.0) cycle read_samp_1
	  if(ifail.eq.0) then
	    call casetrans(cline(left_word(5):right_word(5)),'lo')
	    if(cline(left_word(5):right_word(5)).eq.'x') cycle read_samp_1
	  else
            cols=4
            call linesplit(ifail,4)
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if((icollect.eq.0).and.(abore.eq.aboreold)) cycle read_samp_1
	  if(abore.ne.aboreold)then
	    do i=1,nid1
	      if(abore.eq.aid1(i)) go to 520
	    end do
	    aboreold=abore
	    icollect=0
	    cycle read_samp_1
520	    aboreold=abore
	    icollect=i
	  end if
	  if(icollect.eq.0) cycle read_samp_1
	  call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value, &
          iline,sfile2)
	  if(ifail.ne.0) go to 9900
	  if(((ndays.gt.ndays1(icollect)).or.                               &
             ((ndays.eq.ndays1(icollect)).and.(nsecs.ge.nsecs1(icollect)))) &
             .and.                                                          &
             ((ndays.lt.ndays2(icollect)).or.                               &
             ((ndays.eq.ndays2(icollect)).and.(nsecs.le.nsecs2(icollect)))))&
	    then
            istand=istand+1
	    aabore=abore
	    if(aabore.ne.aaboreold)then
	      loc(icollect)=istand
	      aaboreold=aabore
	    end if
	    aaboreold=aabore
	    sdays(istand)=ndays
	    ssecs(istand)=nsecs
	    sval(istand)=value     
	  end if
	end do read_samp_1

600	continue
	close (unit=sunit2,iostat=ierr)

! -- First we test for if there are no standard measurement for any bore.

	ierr=0
	do i=1,nid1
	  if(loc(i).eq.0) then		!zero loc
	    write(amessage,620) trim(aid1(i)),trim(sfile1),trim(sfile2)
620	    format(' Bore ',a,' cited in file ',a,' is either not cited in ', &
            'the standard bore sample file ',a,' or there are no standard ', &
            'samples in time interval required for calibration of data for ', &
            'this bore.')
	    ierr=ierr+1
	    if(ierr.eq.1) then
	      call write_message(leadspace='yes')
	    else
	      call write_message
	    end if
	  end if
	end do
	if(ierr.ne.0) go to 9900

! -- Now we check whether there are any bores for which there is only one 
!    standard reading.

	if(nid1.eq.1)then
	  if(nstand.eq.1)then
	    write(amessage,640) trim(aid1(1)),trim(sfile1),trim(sfile2)
	    call write_message(leadspace='yes')
	    go to 9900
	  end if
	else
	  do i=1,nid1-1
	    do j=i+1,nid1
	      if(abs(loc(i)-loc(j)).eq.1)then
	        k=i
	        if(loc(j).lt.loc(i)) k=j
	        write(amessage, 640) trim(aid1(k)), trim(sfile1),trim(sfile2)
640	        format(' Data for bore ',a,' cited in file ',a,' cannot be ', &
                'calibrated as there is only one measurement for this bore ', &
                'within calibration time window in file ',a)
	        ierr=ierr+1
	        if(ierr.eq.1) then
	          call write_message(leadspace='yes')
	        else
	          call write_message
	        end if
	      end if
	    end do
	  end do
	  do i=1,nid1
	    if(loc(i).eq.nstand)then
	      write(amessage,640) trim(aid1(i)),trim(sfile1),trim(sfile2)
	      ierr=ierr+1
	      if(ierr.eq.1) then
	        call write_message(leadspace='yes')
	      else
	        call write_message
	      end if
	    end if
	  end do
	  if(ierr.ne.0) go to 9900
	end if

! -- Now we commence the calibration process.

	allocate(bdays(maxdatsmp),bsecs(maxdatsmp),bval(maxdatsmp),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,650)
650	  format(' Cannot allocate sufficient memory to continue execution.')
	  call write_message(leadspace='yes')
	  go to 9900
	end if

	write(repunit,660)
660	format(t25,'DETAILS OF DATA ADJUSTMENT')

	write(6,*)
	aboreold=' '
	iline=0
	nsamp=0
	iid=0
	iend=0
	read_samp_2: do
	  iline=iline+1
	  read(sunit1,'(a)',err=9100,end=720) cline
	  cols=5
          call linesplit(ifail,5)
          if(ifail.lt.0) cycle read_samp_2
	  if(ifail.eq.0) then
	    call casetrans(cline(left_word(5):right_word(5)),'lo')
	    if(cline(left_word(5):right_word(5)).eq.'x') cycle read_samp_2
	  else
            cols=4
            call linesplit(ifail,4)
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if(abore.eq.aboreold) then
	    nsamp=nsamp+1
            call read_rest_of_sample_line(ifail,cols,bdays(nsamp), &
            bsecs(nsamp),bval(nsamp),iline,sfile1)
	    if(ifail.ne.0) go to 9900
	    cycle read_samp_2
	  end if
	  go to 750
720	  iend=1
750	  if(nsamp.eq.0) go to 2000

! --  We've collected the readings; now lets interpolate.

	  iid=iid+1
	  write(repunit,770) trim(aid1(iid))
770	  format(/,/,/,' Data adjustment for bore ',a,':-')
	  write(repunit,780)
780	  format(/,'   Raw data ---->')
	  call newdate(bdays(1),1,1,1970,iday,imon,iyear)
	  call sectime(bsecs(1),isec,imin,ihour)
	  if(datespec.eq.1)then
	    write(repunit,800) iday,imon,iyear,ihour,imin,isec
	  else
	    write(repunit,800) imon,iday,iyear,ihour,imin,isec
	  end if
800	  format('       First sample of raw data at:      ',i2.2,'/',i2.2,'/', &
          i4,'  ',i2.2,':',i2.2,':',i2.2)
	  call newdate(bdays(nsamp),1,1,1970,iday,imon,iyear)
	  call sectime(bsecs(nsamp),isec,imin,ihour)
	  if(datespec.eq.1)then
	    write(repunit,820) iday,imon,iyear,ihour,imin,isec
	  else
	    write(repunit,820) imon,iday,iyear,ihour,imin,isec
	  end if
820	  format('       Last  sample of raw data at:      ',i2.2,'/',i2.2,'/', &
          i4,'  ',i2.2,':',i2.2,':',i2.2)
	  call num2char(nsamp,atemp)
	  write(repunit,840) trim(atemp)
840	  format('       Total number of samples: ',a)
	  write(repunit,860)
860	  format(/,'   Standard data within calibration time frame ---->')
	  iloc1=loc(iid)
	  iloc2=nstand+1
	  do k=1,nid1
	    if(loc(k).le.iloc1) cycle
	    if(loc(k).lt.iloc2) iloc2=loc(k)
	  end do
	  iloc2=iloc2-1
	  call newdate(sdays(iloc1),1,1,1970,iday,imon,iyear)
	  call sectime(ssecs(iloc1),isec,imin,ihour)
	  if(datespec.eq.1)then
	    write(repunit,880) iday,imon,iyear,ihour,imin,isec
	  else
	    write(repunit,880) imon,iday,iyear,ihour,imin,isec
	  end if
880	  format('       First sample of standard data at: ',i2.2,'/',i2.2,'/', &
          i4,'  ',i2.2,':',i2.2,':',i2.2)
	  call newdate(sdays(iloc2),1,1,1970,iday,imon,iyear)
	  call sectime(ssecs(iloc2),isec,imin,ihour)
	  if(datespec.eq.1)then
	    write(repunit,900) iday,imon,iyear,ihour,imin,isec
	  else
	    write(repunit,900) imon,iday,iyear,ihour,imin,isec
	  end if
900	  format('       Last  sample of standard data at: ',i2.2,'/',i2.2,'/', &
          i4,'  ',i2.2,':',i2.2,':',i2.2)
	  call num2char(iloc2-iloc1+1,atemp)
	  write(repunit,840) trim(atemp)
	  write(repunit,920)
920	  format(/'   Calibration equation  -   Y = M*X + C  ---->')
	  write(repunit,940)
940	  format(/,t22,' Interval', t56,'M',t68,'C')
	  write(repunit,945)
945	  format(3x,71('-'))


	  do iloc=iloc1,iloc2-1
	    if(((sdays(iloc).lt.bdays(1)).or.                            &
               ((sdays(iloc).eq.bdays(1)).and.                           &
                                    (ssecs(iloc).le.bsecs(1)))).and.     &
	       ((sdays(iloc+1).lt.bdays(1)).or.                          &
               ((sdays(iloc+1).eq.bdays(1)).and.                         &
                                    (ssecs(iloc+1).le.bsecs(1)))))       &
            then
	         m=1;c=1.0d30
	         go to 959
	    end if
	    if(((sdays(iloc).gt.bdays(nsamp)).or.                        &
               ((sdays(iloc).eq.bdays(nsamp)).and.                       &
                                    (ssecs(iloc).ge.bsecs(nsamp)))).and. &
	       ((sdays(iloc+1).gt.bdays(nsamp)).or.                      &
               ((sdays(iloc+1).eq.bdays(nsamp)).and.                     &
                                    (ssecs(iloc+1).ge.bsecs(nsamp)))))   &
            then
	         m=1;c=1.0d30
	         go to 959
	    end if
	    call time_interp(ifail,nsamp,bdays,bsecs,bval,sdays(iloc), &
            ssecs(iloc),1.0e30,1.0e30,bbval1,extrap='yes',direction='hi')
	    if(ifail.ne.0) go to 9200
	    call time_interp(ifail,nsamp,bdays,bsecs,bval,sdays(iloc+1), &
            ssecs(iloc+1),1.0e30,1.0e30,bbval2,extrap='yes',direction='lo')
	    if(ifail.ne.0) go to 9200
	    ssval1=sval(iloc)
	    ssval2=sval(iloc+1)
	    if(bbval2.eq.bbval1)then
	      write(amessage,950) trim(aid1(iid)),trim(sfile1),trim(sfile1),  &
              trim(sfile2)
950	      format(' Cannot adjust data for bore ',a,' in file ',a,      &
              ' because zero change was recorded in file ',a,' between ',  &
              'times at which standard samples were taken in file ',a)
	      call write_message(increment=1,leadspace='yes')
	      call newdate(sdays(iloc),1,1,1970,iday1,imon1,iyear1)
	      call sectime(ssecs(iloc),isec1,imin1,ihour1)
	      call newdate(sdays(iloc+1),1,1,1970,iday2,imon2,iyear2)
	      call sectime(ssecs(iloc+1),isec2,imin2,ihour2)
	      if(datespec.eq.1)then
	        write(amessage,975) trim(sfile2),                          &
                                  iday1,imon1,iyear1,ihour1,imin1,isec1,   &
                                  iday2,imon2,iyear2,ihour2,imin2,isec2
	      else
	        write(amessage,975) trim(sfile2),                          &
                                  imon1,iday1,iyear1,ihour1,imin1,isec1,   &
                                  imon2,iday2,iyear2,ihour2,imin2,isec2
	      end if
975	      format('   Sample times (',a,'): ',                          &
              i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2,' and ',     &
              i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2)
	      call write_message
	      m=1.0d0; c=1.0d30
	    else if(ssval2.eq.ssval1)then
	      write(amessage,955) trim(aid1(iid)),trim(sfile1),trim(sfile2)
955	      format(' Cannot adjust data for bore ',a,' in file ',a,      &
              ' because there is zero change recorded between two ',       &
              'neighbouring sample points for this bore in file ',a)
	      call write_message(increment=1,leadspace='yes')
	      call newdate(sdays(iloc),1,1,1970,iday1,imon1,iyear1)
	      call sectime(ssecs(iloc),isec1,imin1,ihour1)
	      call newdate(sdays(iloc+1),1,1,1970,iday2,imon2,iyear2)
	      call sectime(ssecs(iloc+1),isec2,imin2,ihour2)
	      if(datespec.eq.1)then
	        write(amessage,975) trim(sfile2),                          &
                                  iday1,imon1,iyear1,ihour1,imin1,isec1,   &
                                  iday2,imon2,iyear2,ihour2,imin2,isec2
	      else
	        write(amessage,975) trim(sfile2),                          &
                                  imon1,iday1,iyear1,ihour1,imin1,isec1,   &
                                  imon2,iday2,iyear2,ihour2,imin2,isec2
	      end if
	      call write_message
	      m=1.0d0; c=1.0d30
	    else
	      m=(ssval2-ssval1)/(bbval2-bbval1)
	      c=ssval1-m*bbval1
	    end if
959	    continue
	    if(iloc2-iloc1.eq.1) then
	      begday=-1000000; begsec=0
	      finday=1000000; finsec=0
	    else if(iloc.eq.iloc1) then
	      begday=-1000000; begsec=0
	      finday=sdays(iloc+1); finsec=ssecs(iloc+1)
	    else if(iloc.eq.iloc2-1)then
	      begday=sdays(iloc); begsec=ssecs(iloc)
	      finday=1000000; finsec=0
	    else
	      begday=sdays(iloc); begsec=ssecs(iloc)
	      finday=sdays(iloc+1); finsec=ssecs(iloc+1)
	    end if
	    k=0
	    do i=1,nsamp
	      if(((bdays(i).gt.begday).or.                                  &
                 ((bdays(i).eq.begday).and.(bsecs(i).ge.begsec))).and.     &
                 ((bdays(i).lt.finday).or.                                  &
                 ((bdays(i).eq.finday).and.(bsecs(i).lt.finsec))))then
	         k=k+1
       	         call newdate(bdays(i),1,1,1970,iday,imon,iyear)
	         call sectime(bsecs(i),isec,imin,ihour)
	         if(datespec.eq.1) then
	           write(calunit,980) trim(aid1(iid)),iday,imon,iyear,ihour,  &
                   imin,isec,m*bval(i)+c
	         else
	           write(calunit,980) trim(aid1(iid)),imon,iday,iyear,ihour,  &
                   imin,isec,m*bval(i)+c
	         end if
980	         format(1x,a,t13,i2.2,'/',i2.2,'/',i4,2x,i2.2,':',i2.2,':',  &
                 i2.2,2x,1pg14.7)
	      end if
	    end do
	    if(k.eq.0) c=1.0d30

	    call newdate(sdays(iloc),1,1,1970,iday1,imon1,iyear1)
	    call sectime(ssecs(iloc),isec1,imin1,ihour1)
	    call newdate(sdays(iloc+1),1,1,1970,iday2,imon2,iyear2)
	    call sectime(ssecs(iloc+1),isec2,imin2,ihour2)
	    if(c.lt.0.9d30)then
	      if(datespec.eq.1)then
	        write(repunit,960) iday1,imon1,iyear1,ihour1,imin1,isec1,  &
	        iday2,imon2,iyear2,ihour2,imin2,isec2,m,c
	      else
	        write(repunit,960) imon1,iday1,iyear1,ihour1,imin1,isec1,  &
	        imon2,iday2,iyear2,ihour2,imin2,isec2,m,c
	      end if
960	      format(4x,i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2,  &
              ' to ',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2,   &
              t52,1pg11.4,1x,1pg11.4)
	    else
	      if(datespec.eq.1)then
	        write(repunit,965) iday1,imon1,iyear1,ihour1,imin1,isec1,  &
	        iday2,imon2,iyear2,ihour2,imin2,isec2,'   not used'
	      else
	        write(repunit,965) imon1,iday1,iyear1,ihour1,imin1,isec1,  &
	        imon2,iday2,iyear2,ihour2,imin2,isec2,'   not used'
	      end if
965	      format(4x,i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2,  &
              ' to ',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2,   &
              t57,a)
	    end if

	  end do


2000	  if(iend.eq.1) go to 8000
	  aboreold=abore
	  nsamp=1
          call read_rest_of_sample_line(ifail,cols,bdays(nsamp), &
          bsecs(nsamp),bval(nsamp),iline,sfile1)
	  if(ifail.ne.0) go to 9900
          cycle read_samp_2
	end do read_samp_2


8000	continue
	if(imessage.eq.0)then
	  write(6,8010) trim(calfile)
8010	  format(/,' - new bore sample file ',a,' written ok.')
	  write(6,8040) trim(repfile1)
8040	  format(' - see file ',a,' for calibration details.')
	end if
	go to 9900


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(sfile2)
9010    format(' Error reading line ',a,' of bore sample file ',a)
        call write_message(leadspace='yes')
        go to 9900
9100	write(amessage, 9110) trim(sfile1)
9110	format( ' Error re-reading bore sample file ',a)
	call write_message(leadspace='yes')
	go to 9900
9200	write(amessage, 9210) trim(sfile1)
9210	format(' Error using data from bore sample file ',a,'; check ', &
        'this file with SMPCHEK.')
	call write_message(leadspace='yes')
	go to 9900


9900	call close_files
	deallocate(aid1,ndays1,nsecs1,ndays2,nsecs2,loc,stat=ierr)
	deallocate(sdays,ssecs,sval,stat=ierr)
	deallocate(bdays,bsecs,bval,stat=ierr)
	if(imessage.ne.0) then
	  write(6,*)
	  if(calfile.ne.' ') then
            inquire(file=trim(calfile),exist=lexist)
	    if(lexist) call system('del "'//trim(calfile)//'"')
	  end if
	  if(repfile1.ne.' ')then
	    inquire(file=trim(repfile1),exist=lexist)
	    if(lexist) call system('del "'//trim(repfile1)//'"')
	  end if
	end if


end program smpcal
 