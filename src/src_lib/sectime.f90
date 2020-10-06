subroutine sectime(nsecs,sec,min,hour)

! -- Subroutine SECTIME expresses elapsed time since midnight as hours,
!    minutes and seconds.
!    NSECS must be positive.
!    If the resulting hours are greater than 24, no correction is made.

! -- Arguments are as follows:-
!	nsecs:                 elapsed number of seconds
!       sec,min,hour:          seconds, minutes, hours of time of day

! -- Revision history:-
!	April, 1997: version 1.

	use inter
	implicit none

	integer, intent(in)   :: nsecs
	integer, intent(out)  :: sec,min,hour

	integer               :: tsecs

	if(nsecs.lt.0) call sub_error('NEWTIME')

	hour=nsecs/3600
	tsecs=nsecs-hour*3600
	min=tsecs/60
	sec=tsecs-min*60

	return

end subroutine sectime
 