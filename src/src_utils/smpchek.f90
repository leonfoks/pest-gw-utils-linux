!     Last change:  JD   14 Feb 2001    1:33 pm
program smpchek

! -- Program SMPCHEK checks the integrity of a bore sample file.

	use defn
	use inter

	implicit none

	integer	:: ifail,inunit,iline,col5,ibore,j,oldbore,no_date_test,&
		   dd,mm,yy,ddold,mmold,yyold,ndays,igt10,ierr,hhh,mmm,sss,&
		   hhhold,mmmold,sssold,nsecs,idate,iheader,nbb
	character (len=80) :: aprompt,infile,atemp1,atemp2
	character (len=20) :: adate
	character (len=10) :: aline,abore,acode
	real	:: sample

	integer, parameter			:: NUM_BORE=10000
	integer, dimension(NUM_BORE)		:: nbore
	character (len=10), dimension(NUM_BORE)	:: boreid


	write(amessage,5)
5       format(' Program SMPCHEK checks the integrity of a bore sample file.')
	call write_message(leadspace='yes',endspace='yes')

	call readfig(atemp1,atemp2,infile)

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

        yy=0
        dd=0
        mm=0
        hhh=0
        mmm=0
        sss=0

4	imessage=0
6	if(imessage.eq.5) go to 9900
13	if(infile.eq.' ')then
11	  write(6,15,advance='no')
15	  format(' Enter name of bore sample file: ')
	  read(5,'(a)') atemp1
	  if(atemp1.eq.' ') go to 11
	  atemp1=adjustl(atemp1)
	  if(index(eschar,atemp1(1:2)).ne.0) go to 9900
          nbb=len_trim(atemp1)
          call getfile(ifail,atemp1,atemp2,1,nbb)
          if(ifail.ne.0) go to 11
          atemp1=atemp2
	else
	  write(6,16,advance='no') trim(infile)
16	  format(' Enter name of bore sample file [',a,']: ')
	  read(5,'(a)') atemp1
      	  if(index(eschar,atemp1(1:2)).ne.0) go to 9900
	  if(atemp1.eq.' ') then
            atemp1=infile
          else
            nbb=len_trim(atemp1)
            call getfile(ifail,atemp1,atemp2,1,nbb)
            if(ifail.ne.0) go to 13
            atemp1=atemp2
          end if
	end if
	inunit=nextunit()
	open(unit=inunit,file=atemp1,status='old',iostat=ierr)
	if(ierr.ne.0) then
	  write(amessage,20) trim(atemp1)
20	  format(' Cannot open file ',a,'  - try again.')
	  call write_message(increment=1)
	  if(atemp1.eq.infile) infile=' '
	  go to 6
	end if
	infile=atemp1

10	format(/,' Errors in bore sample file ',a,' ----->')

	imessage=0
	iline=0
	ibore=0
	nbore=0
	read_a_line: do

	  if(imessage.gt.40) go to 9900
	  igt10=0

	  iline=iline+1
	  if((imessage.eq.0).and.(iline.eq.2000))then
	    write(6,25)
25	    format(/,' Working.....')
	  end if
	  call num2char(iline,aline)
	  col5=0
	  read(inunit,'(a)',err=9000,end=1000) cline
	  call linesplit(ifail,5)
	  if(ifail.eq.-1) cycle read_a_line
	  if(ifail.eq.0) then
	    col5=1
	  else
	    call linesplit(ifail,4)
	    if(ifail.ne.0) then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,30) trim(aline)
30	      format('   Line ',a,': insufficient data items (minimum of ',&
	      '4 items required).')
	      call write_message(increment=1)
	      ddold=-9999
	      cycle read_a_line
	    end if
	  end if

	  if(right_word(1)-left_word(1).gt.9) then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,50) trim(aline)
50	    format('   Line ',a,': bore identifier greater than 10 characters ',&
	    'in length.')
	    call write_message(increment=1)
	    igt10=1
	    go to 110
	  end if
	  abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if(ibore.eq.0)then
	    ibore=1
	    oldbore=1
	    if(ibore.gt.NUM_BORE) go to 9050
	    boreid(ibore)=abore
	    nbore(ibore)=1
	    no_date_test=1
	  else
	    if(abore.ne.boreid(oldbore)) then
	      do j=1,ibore
	        if(abore.eq.boreid(j))then
	          if(imessage.eq.0) write(6,10) trim(infile)
		  write(amessage,60) trim(aline),trim(abore), trim(abore)
60	          format('   Line ',a,': bore ',a,' cited previously - ', &
		  'entries for bore ',a,' are not consecutive.')
		  call write_message(increment=1)
	          nbore(j)=nbore(j)+1
	  	  no_date_test=1
	 	  oldbore=j
		  go to 110
	 	end if
	      end do
	      ibore=ibore+1
	      if(ibore.gt.NUM_BORE) go to 9050
	      boreid(ibore)=abore
	      nbore(ibore)=1
	      no_date_test=1
	      oldbore=ibore
	    else
	      nbore(ibore)=nbore(ibore)+1
	      no_date_test=0
	    end if
	  end if
 
110	  adate=cline(left_word(2):right_word(2))
	  call char2date(ifail,adate,dd,mm,yy)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,120) trim(aline)
120	    format('   Line ',a,': illegal date.')
	    call write_message(increment=1)
	    dd=-9999
	  end if
	  if(igt10.eq.1) dd=-9999
	  call char2time(ifail,cline(left_word(3):right_word(3)),hhh,mmm,sss)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,130) trim(aline)
130	    format('   Line ',a,': illegal time.')
	    call write_message(increment=1)
	    dd=-9999
	  end if
	  if((no_date_test.ne.1).and.(dd.ne.-9999).and.(ddold.ne.-9999))then
	    ndays=numdays(ddold,mmold,yyold,dd,mm,yy)
	    nsecs=numsecs(hhhold,mmmold,sssold,hhh,mmm,sss)
	    if((ndays.eq.0).and.(nsecs.eq.0))then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,140) trim(aline)
140	      format('   Line ',a,': sample occurs at same time as previous sample.')
	      call write_message(increment=1)
	    else if((ndays.lt.0).or.((ndays.eq.0).and.(nsecs.lt.0))) then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,160) trim(aline)
160	      format('   Line ',a,': sample occurs earlier than previous sample.')
	      call write_message(increment=1)
	    end if
	  end if
	  ddold=dd
	  mmold=mm
	  yyold=yy
	  hhhold=hhh
	  mmmold=mmm
	  sssold=sss

	  sample=char2double(ifail,4)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,180) trim(aline)
180	    format('   Line ',a,': cannot read sample value.')
	    call write_message(increment=1)
	  end if

	  if(col5.eq.1)then
	    acode=cline(left_word(5):right_word(5))
	    call casetrans(acode,'lo')
	    if((len_trim(acode).ne.1).or.(acode.ne.'x'))then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,200) trim(aline)
200	      format('   Line ',a,': optional fifth column can only contain "x".')
	      call write_message(increment=1)
	    end if
	  end if

	end do read_a_line

1000	call num2char (iline-1,aline)
	write(amessage,1010) trim(aline), trim(infile)
1010	format(' - ',a,' lines read from bore sample file ',a)
	call write_message(leadspace='yes')
	if(imessage.eq.0)then
	  write(amessage,1020)
1020	  format(' - no errors found.')
	  call write_message
	  call num2char(ibore,aline)
	  write(amessage,1030) trim(aline)
1030	  format(' - sample values found for ',a,' different bores.')
	  call write_message
	end if
	go to 9900
 
9000	if(imessage.eq.0) write(6,10) trim(infile)
	write(amessage,9010) trim(aline)
9010	format('   Line ',a,': cannot read line.')
	call write_message
	go to 9900
9050	call num2char(NUM_BORE,aline)
	write(amessage,9060) trim(infile),trim(aline)
9060	format(' Cannot continue execution. There are too many different ',&
	'bores referenced in file ',a,'. Redimension NUM_BORE greater than ',&
	a,' in SMPCHEK source code and recompile.')
	call write_message(leadspace='yes')
	go to 9900

9900    call close_files
	write(6,*)

end program smpchek

 