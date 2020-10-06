subroutine read_rest_of_sample_line(ifail,cols,ndays,nsecs,value,iline,sampfile)

! -- Subroutine read_rest_of_sample_line reads the date, time, value and
!    optional fifth column from a line of a bore sample file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition is encountered
!       cols:      number of data columns in the line
!       ndays:     number of days from 1/1/1970 until sample date
!       nsecs:     number of seconds from midnight until sample time
!       value:     sample value
!       iline:     current line number of bore sample file
!       sampfile:  name of bore sample file

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)            :: ifail
	integer, intent(in)             :: cols
	integer, intent(out)            :: ndays,nsecs
	double precision, intent(out)   :: value
	integer, intent(in)             :: iline
	character (len=*), intent(in)   :: sampfile
	integer                         :: dd,mm,yy,hhh,mmm,sss
	character (len=15)              :: aline
	character (len=2)               :: aa

	ifail=0
	call char2date(ifail,cline(left_word(2):right_word(2)),dd,mm,yy)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,150) trim(aline),trim(sampfile)
150       format('illegal date at line ',a,' of bore sample file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	ndays=numdays(1,1,1970,dd,mm,yy)

	call char2time(ifail,cline(left_word(3):right_word(3)),hhh,mmm,sss)
	if(ifail.ne.0) then
	  call num2char(iline,aline)
	  write(amessage,160) trim(aline),trim(sampfile)
160       format('illegal time at line ',a,' of bore sample file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	nsecs=numsecs(0,0,0,hhh,mmm,sss)

	value=char2double(ifail,4)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,180) trim(aline),trim(sampfile)
180       format('cannot read sample value at line ',a,' of bore sample file ',a)
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	if(value.lt.-1.0e37) then
	  call num2char(iline,aline)
	  write(amessage,190) trim(aline),trim(sampfile)
190       format('illegal sample value at line ',a,' of bore sample file ',a, &
	  '; lower limit is -1.0E37.')
	  call write_message(error='yes',leadspace='yes')
	  go to 9800
	end if
	if(cols.eq.5)then
	  aa=cline(left_word(5):right_word(5))
	  call casetrans(aa,'lo')
	  if(aa.eq.'x ') then
	    value=-1.1e38
	  else
	    call num2char(iline,aline)
	    write(amessage,210) trim(aline),trim(sampfile)
210         format('illegal optional fifth item on line ',a,' of bore sample ',&
	    'file ',a,'; item must be "x" if present.')
	    call write_message(error='yes',leadspace='yes')
	    go to 9800
	  end if
	end if
	return

9800    ifail=1
	return

end subroutine read_rest_of_sample_line
 