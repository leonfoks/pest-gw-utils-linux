subroutine read_rest_of_pump_line(ifail,ibore,ndays,nsecs,pumped,iline,pmpfile)

! -- Subroutine read_rest_of_pump_line reads the dates, times and pumped volumes
! -- from a line of a bore pumping file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition is encountered
!       ibore:     internal bore index
!       ndays:     number of days since 1/1/1970 for each sample
!       nsecs:     number of seconds elapsed since midnight for each sample
!       pumped:    amount pumped since previous sample for each sample
!       iline:     current line number of bore pumping file
!       pmpfile:   name of bore pumping file

!    Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                    :: ifail
	integer, intent(inout)                  :: ibore
	integer, intent(out), dimension (:)     :: ndays,nsecs
	double precision, intent(inout), dimension(:)   :: pumped
	integer, intent(in)                     :: iline
	character (len=*), intent(in)           :: pmpfile
	integer                                 :: dd,mm,yy,iday1,iday2
	integer                                 :: hhh,mmm,sss,itime1,itime2
	character (len=15)                      :: aline
	real                                    :: pumpage

	ifail=0
	call char2date(ifail,cline(left_word(2):right_word(2)),dd,mm,yy)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,150) trim(aline),trim(pmpfile)
150       format('illegal first date at line ',a,' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	iday1=numdays(1,1,1970,dd,mm,yy)

	call char2time(ifail,cline(left_word(3):right_word(3)),hhh,mmm,sss)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,110) trim(aline),trim(pmpfile)
110       format('illegal first time at line ',a,' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	itime1=numsecs(0,0,0,hhh,mmm,sss)

	if(ibore.ne.1) then
	  if((iday1.ne.ndays(ibore-1)).or.(itime1.ne.nsecs(ibore-1))) then
	    call num2char(iline,aline)
	    write(amessage,100) trim(aline),trim(pmpfile)
100         format('first date/time on line ',a,' of bore pumping file ',a,&
	    ' does not match second date/time on previous line. ')
	    call write_message(error='yes',leadspace='yes')
	    go to 9800
	  end if
	end if
	call char2date(ifail,cline(left_word(4):right_word(4)),dd,mm,yy)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,160) trim(aline),trim(pmpfile)
160       format('illegal second date at line ',a,' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	iday2=numdays(1,1,1970,dd,mm,yy)

	call char2time(ifail,cline(left_word(5):right_word(5)),hhh,mmm,sss)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,130) trim(aline),trim(pmpfile)
130       format('illegal second time at line ',a,' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	itime2=numsecs(0,0,0,hhh,mmm,sss)

	if((iday2.lt.iday1).or. &
	  ((iday2.eq.iday1).and.(itime2.le.itime1))) then
	  call num2char(iline,aline)
	  write(amessage,170) trim(aline),trim(pmpfile)
170       format('second date/time does not follow first date/time at line ',a,&
	  ' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	pumpage=char2real(ifail,6)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,180) trim(aline),trim(pmpfile)
180       format('cannot read pumping value at line ',a,' of bore pumping file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	if(pumpage.lt.-1.0e37) then
	  call num2char(iline,aline)
	  write(amessage,190) trim(aline),trim(pmpfile)
190       format('illegal pumping value at line ',a,' of bore pumping file ',a, &
	  '  - lower limit is -1.0E37.')
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	if(ibore.eq.1)then
	  ndays(1)=iday1
	  nsecs(1)=itime1
	  pumped(1)=0.0d0
	  ibore=ibore+1
	end if
	ndays(ibore)=iday2
	nsecs(ibore)=itime2
	pumped(ibore)=pumped(ibore-1)+pumpage
	return

9800    ifail=1
	return

end subroutine read_rest_of_pump_line
 