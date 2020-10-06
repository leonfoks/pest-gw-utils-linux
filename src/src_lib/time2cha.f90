subroutine time2char(ifail,hh,mm,ss,atime)

! -- Subroutine TIME2CHAR converts hours, minutes, seconds to a date string.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       hh,mm,ss:  hours minutes seconds
!       atime:     the date string produced by time2char

! -- Revision history:-
!       June-November, 1995: version 1.

	integer, intent(out)            :: ifail
	integer, intent(in)             :: hh,mm,ss
	character (len=*), intent(out)  :: atime

	ifail=0
	if((hh.lt.0).or.(hh.gt.23).or.(mm.lt.0).or.(mm.gt.59).or.&
	   (ss.lt.0).or.(ss.gt.59)) go to 9000
	if(len(atime).lt.8) go to 9000
	write(atime(1:2),'(i2.2)',err=9000) hh
	write(atime(4:5),'(i2.2)',err=9000) mm
	write(atime(7:8),'(i2.2)',err=9000) ss
	atime(3:3)=':'
	atime(6:6)=':'
	return

9000    ifail=1
	return

end subroutine time2char
 