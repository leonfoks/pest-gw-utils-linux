!     Last change:  JD   13 Feb 2001    9:26 pm
program pmpchek

! -- Program PMPCHEK checks the integrity of a bore pumping file.


	use defn
	use inter

	implicit none

	integer	:: ifail,inunit,iline,ibore,j,oldbore,no_date_test,&
		   dd1,mm1,yy1,dd2,mm2,yy2,dd2old,mm2old,yy2old,&
		   igt10,ierr,nd,hhh1,mmm1,sss1,hhh2,mmm2,sss2,hhh2old, &
		   mmm2old,sss2old,idate,iheader,nbb,ifail1
	character (len=80) :: infile,atemp1,atemp2,atemp3,bfile,cfile
	character (len=20) :: adate,atime
	character (len=10) :: aline,abore
	real	:: sample

	integer, parameter			:: NUM_BORE=5000
	integer, dimension(NUM_BORE)		:: nbore
	character (len=10), dimension(NUM_BORE)	:: boreid


	write(amessage,5)
5       format(' Program PMPCHEK checks the integrity of a bore pumping file.')
	call write_message(leadspace='yes',endspace='yes')

	call readfig(atemp1,atemp2,atemp3,infile)

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

4	imessage=0
6	if(imessage.eq.5) go to 9900
	if(infile.eq.' ')then
11	  write(6,15,advance='no')
15	  format(' Enter name of bore pumping file: ')
	  read(5,'(a)') atemp1
	  if(atemp1.eq.' ') go to 11
	  atemp1=adjustl(atemp1)
	  if(index(eschar,atemp1(1:2)).ne.0) go to 9900
	else
	  write(6,16,advance='no') trim(infile)
16	  format(' Enter name of bore pumping file [',a,']: ')
	  read(5,'(a)') atemp1
	  if(atemp1.eq.' ') then
            call addquote(infile,atemp1)
          end if
	  atemp1=adjustl(atemp1)
	  if(index(eschar,atemp1(1:2)).ne.0) go to 9900
	end if
        bfile=atemp1
        nbb=len_trim(bfile)
        call getfile(ifail1,bfile,cfile,1,nbb)
        if(ifail1.ne.0) then
          go to 6
        else
          atemp1=cfile
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
	  if((imessage.eq.0).and.(iline.eq.1000))then
	    write(6,25)
25	    format(/,' Working.....')
	  end if
	  call num2char(iline,aline)
	  read(inunit,'(a)',err=9000,end=1000) cline
	  call linesplit(ifail,6)
	  if(ifail.eq.-1) cycle read_a_line
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,30) trim(aline)
30	    format('   Line ',a,': insufficient data items (6 items required).')
	    call write_message(increment=1)
	    dd2old=-9999
	    cycle read_a_line
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
	  call char2date(ifail,adate,dd1,mm1,yy1)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,120) trim(aline)
120	    format('   Line ',a,': illegal first date.')
	    call write_message(increment=1)
	    dd1=-9999
	  end if
	  atime=cline(left_word(3):right_word(3))
	  call char2time(ifail,atime,hhh1,mmm1,sss1)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,130) trim(aline)
130	    format('   Line ',a,': illegal first time.')
	    call write_message(increment=1)
	    dd1=-9999
	  end if
	  if(igt10.eq.1) dd1=-9999
	  if((no_date_test.ne.1).and.(dd1.ne.-9999).and.(dd2old.ne.-9999))then
	    if((dd1.ne.dd2old).or.(mm1.ne.mm2old).or.(yy1.ne.yy2old).or. &
               (hhh1.ne.hhh2old).or.(mmm1.ne.mmm2old).or.(sss1.ne.sss2old))then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,140) trim(aline)
140	      format('   Line ',a,': first date/time does not agree with ',&
	      'second date/time on previous line.')
	      call write_message(increment=1)
	    end if
	  end if
	  adate=cline(left_word(4):right_word(4))
	  call char2date(ifail,adate,dd2,mm2,yy2)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,150) trim(aline)
150	    format('   Line ',a,': illegal second date.')
	    call write_message(increment=1)
	    dd2=-9999
	  end if
	  atime=cline(left_word(5):right_word(5))
	  call char2time(ifail,atime,hhh2,mmm2,sss2)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,155) trim(aline)
155	    format('   Line ',a,': illegal second time.')
	    call write_message(increment=1)
	    dd2=-9999
	  endif
	  if((dd1.ne.-9999).and.(dd2.ne.-9999))then
	    nd=numdays(dd1,mm1,yy1,dd2,mm2,yy2)
	    if((nd.lt.0).or.((nd.eq.0).and. &
	      (numsecs(hhh1,mmm1,sss1,hhh2,mmm2,sss2).le.0))) then
	      if(imessage.eq.0) write(6,10) trim(infile)
	      write(amessage,160) trim(aline)
160	      format('   Line ',a,': second date/time must follow first date/time.')
	      call write_message(increment=1)
	    end if
	  end if 
	  dd2old=dd2
	  mm2old=mm2
	  yy2old=yy2
	  hhh2old=hhh2
	  mmm2old=mmm2
	  sss2old=sss2

	  sample=char2double(ifail,6)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) write(6,10) trim(infile)
	    write(amessage,180) trim(aline)
180	    format('   Line ',a,': cannot read pumping value.')
	    call write_message(increment=1)
	  end if

	end do read_a_line

1000	call num2char (iline-1,aline)
	write(amessage,1010) trim(aline), trim(infile)
1010	format(' - ',a,' lines read from bore pumping file ',a)
	call write_message(leadspace='yes')
	if(imessage.eq.0)then
	  write(amessage,1020)
1020	  format(' - no errors found.')
	  call write_message
	  call num2char(ibore,aline)
	  write(amessage,1030) trim(aline)
1030	  format(' - pumping values found for ',a,' different bores.')
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
	a,' in PMPCHEK source code and recompile.')
	call write_message(leadspace='yes')
	go to 9900

9900    call close_files
	write(6,*)

end program pmpchek

 