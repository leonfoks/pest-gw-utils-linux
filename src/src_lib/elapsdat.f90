subroutine elapsdate(eltime,dayfactor,day1,mon1,year1,hour1,min1,sec1, &
	day2,mon2,year2,hour2,min2,sec2)

! -- Subroutine elapsdate calculates a new date and time, given an original
!    date and time and an elapsed time.

! -- Arguments are as follows:-
!       eltime:              elapsed time
!       dayfactor:           factor to convert eltime to days
!       day1, mon1, year1:   day, month and year of first date
!       hour1, min1, sec1:   hours, minutes and seconds of first time
!       day2, mon2, year2:   day, month and year of second date
!       hour2, min2, sec2:   hours, minutes and seconds of second time

! -- Revision history:-
!       June-November, 1995: version 1.

	use inter
	implicit none

	real, intent(in)        :: eltime,dayfactor
	integer, intent(in)     :: day1,mon1,year1,hour1,min1,sec1
	integer, intent(out)    :: day2,mon2,year2,hour2,min2,sec2

	integer iday,isec,day3,mon3,year3
	real daytime,rday

	if(eltime.lt.0.0) call sub_error('ELAPSDATE')

	daytime=eltime*dayfactor
	iday=floor(daytime)
	if(iday.gt.0) then
	  call newdate(iday,day1,mon1,year1,day2,mon2,year2)
	else
	  day2=day1
	  mon2=mon1
	  year2=year1
	end if

	rday=iday
	if(abs(daytime-rday).lt.abs(2.0*spacing(daytime)))then
	  hour2=hour1
	  min2=min1
	  sec2=sec1
	  return
	end if

	rday=daytime-rday
	isec=nint(rday*86400.0)
	isec=isec+numsecs(0,0,0,hour1,min1,sec1)
	iday=isec/86400
	isec=isec-iday*86400
	hour2=isec/3600
	isec=isec-hour2*3600
	min2=isec/60
	sec2=isec-min2*60
	if((hour2.gt.24).or.(min2.gt.60).or.(sec2.gt.60).or.(iday.gt.1).or. &
	   (hour2.lt.0).or.(min2.lt.0).or.(sec2.lt.0)) &
	   call sub_error('ELAPSDATE')

	if(iday.eq.1) then
	  call newdate(1,day2,mon2,year2,day3,mon3,year3)
	  day2=day3
	  mon2=mon3
	  year2=year3
	end if

	return

end subroutine elapsdate
 