logical function leap(year)

! -- Function LEAP returns .true. if a year is a leap year.

! -- Revision history:-
!       June-November, 1995: version 1.

	integer, intent(in)     :: year

        leap = ( mod(year,4).eq.0 .and. mod(year,100).ne.0 ) .or. &
               ( mod(year,400).eq.0 .and. year.ne.0 )

	return
end function leap
 