integer function numsecs(h1,m1,s1,h2,m2,s2)

! -- Subroutine NUMSECS calculates the number of seconds between two times.

! -- Arguments are as follows:-
!       h1,m1,s1:   hours, minutes seconds of first time
!       h2,m2,y2:   hours, minutes seconds of second time

! -- Revision history:-
!       June-November 1995: version 1.

	integer, intent(in)             :: h1,m1,s1,h2,m2,s2

	numsecs=(h2-h1)*3600+(m2-m1)*60+s2-s1

end function numsecs
 