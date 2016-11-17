subroutine newdate(ndays,day1,mon1,year1,day2,mon2,year2)

! -- Subroutine NEWDATE evaluates the date after NDAYS days have elapsed from
!    a provided date. NDAYS may be negative.

! -- Arguments are as follows:-
!       ndays:            elapsed number of days
!       day1,mon1,year1:  days, month and year of first date
!       day2,mon2,year2:  days, month and year of second date

! -- Revision history:-
!       June-November, 1995: version 1.

	use inter
	implicit none

	integer, intent(in)     :: ndays,day1,mon1,year1
	integer, intent(out)    :: day2,mon2,year2

	integer  :: yearref,newdays,idays,iyear,jdays,i
	integer, dimension(12) :: monthdays

	data monthdays /31,28,31,30,31,30,31,31,30,31,30,31/

! -- First a reference date is chosen. This is the beginning of the first
! -- year. Alternatively the reference date is the beginning of a year prior
! -- to the likely calculated date if NDAYS is negative.

	if(ndays.ge.0) then
	  yearref=year1
	else
	  yearref=year1-abs(ndays)/365-1
	end if
	newdays=numdays(31,12,yearref-1,day1,mon1,year1)
	newdays=ndays+newdays
	if(newdays.lt.0) call sub_error('NEWDATE')

! -- Next days are counted, starting at the new reference date.

	idays=0
	iyear=yearref
	do
	  jdays=idays+365
	  if(leap(iyear)) jdays=jdays+1
	  if(jdays.ge.newdays) go to 20
	  iyear=iyear+1
	  idays=jdays
	end do
	call sub_error('NEWDATE')
20      year2=iyear

	do i=1,12
	  jdays=idays+monthdays(i)
	  if((i.eq.2).and.(leap(year2))) jdays=jdays+1
	  if(jdays.ge.newdays) go to 40
	  idays=jdays
	end do
	call sub_error('NEWDATE')
40      mon2=i
	day2=newdays-idays
	if((day2.le.0).or.(mon2.le.0).or.(year2.le.0)) call sub_error('NEWDATE')

	return

end subroutine newdate  
 