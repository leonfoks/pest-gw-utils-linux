subroutine datestring(dd,mm,yy,hhh,mmm,sss,time,at,adate,atime)

! -- Subroutine DATESTRING produces an ascii date string and time string
!    detailing the date and time after a certain amount of time has elapsed
!    since a user-supplied starting time.

! -- Arguments are as follows:-
!	dd,mm,yy:		initial day, month and year
!	hhh,mmm,sss:		initial hours, minutes and seconds
!	time:			elapsed time since initial time
!	at:			code indicating units of elapsed time
!	adate:			output date string
!	atime:			output time string

! -- Revision history:-
!	December, 1998: version 1.

	use defn
	use inter
	implicit none

	integer, intent(in)		:: dd,mm,yy,mmm,sss,hhh
	real, intent(in)		:: time
	character (len=1), intent(in)	:: at
	character (len=*), intent(out)	:: adate, atime

	integer				:: idays,isecs,dd1,mm1,yy1,isec1,  &
	                                   hhh1,mmm1,sss1
	real				:: rdays,rhours,rmins


	if(time.lt.0.0) call sub_error('DATESTRING')



! --	First the number of elapsed days and seconds are ascertained.



	if((at.eq.'y').or.(at.eq.'Y'))then
	  rdays=time*365.25
	  idays=rdays
	  isecs=nint((rdays-idays)*86400.0)
	else if((at.eq.'d').or.(at.eq.'D'))then	
	  idays=time
	  isecs=nint((time-idays)*86400.0)
	else if((at.eq.'h').or.(at.eq.'H'))then
	  rdays=time/24.0
	  idays=rdays
	  rhours=time-(idays*24)
	  isecs=nint(rhours*3600.0)
	else if((at.eq.'m').or.(at.eq.'M'))then
	  rdays=time/1440.0
	  idays=rdays
	  rmins=time-(idays*1440)
	  isecs=nint(rmins*60.0)
	else
	  rdays=time/86400.0
	  idays=rdays
	  isecs=nint(time-86400*idays)
	end if

! -- The new date is now evaluated.

	call newdate(idays,dd,mm,yy,dd1,mm1,yy1)
	isec1=hhh*3600+mmm*60+sss
	isec1=isec1+isecs
	if(isec1.ge.86400)then
	  call newdate(idays+1,dd,mm,yy,dd1,mm1,yy1)
	  isec1=isec1-86400
	end if
	hhh1=isec1/3600
	mmm1=(isec1-hhh1*3600)/60
	sss1=isec1-hhh1*3600-mmm1*60

! -- The output time and date strings are now written.
        if((yy1.ge.-999).and.(yy1.le.9999))then
          if(datespec.eq.1)then
            write(adate,50) dd1,mm1,yy1
50	    format(i2.2,'/',i2.2,'/',i4)
	  else if(datespec.eq.2) then
            write(adate,50) mm1,dd1,yy1
          endif
        else
          if(datespec.eq.1)then
            write(adate,51) dd1,mm1,yy1
51	    format(i2.2,'/',i2.2,'/',i5)
	  else if(datespec.eq.2) then
            write(adate,51) mm1,dd1,yy1
          endif
        end if
	adate=adjustl(adate)
	write(atime,60)hhh1,mmm1,sss1
60	format(i2.2,':',i2.2,':',i2.2)
	atime=adjustl(atime)

	return

end subroutine datestring
 