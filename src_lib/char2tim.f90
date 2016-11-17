subroutine char2time(ifail,atime,hh,mm,ss)

! -- Subroutine CHAR2TIME extracts the time from a string.

! -- Arguments are as follows:-
!       ifail:     indicates failure if returned as non-zero
!       atime:     a string containing the time in ASCII format
!       hh,mm,ss   hours, minutes and seconds extracted from the atime string.

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: atime
	integer, intent(out)            :: hh,mm,ss
	integer                         :: lentime,i,j
	character (len=2)               :: asep
	character (len=20)              :: btime

	ifail=0
	asep=':.'
	if(atime.eq.' ') go to 9000
	btime=adjustl(atime)
	lentime=len_trim(btime)
	if(lentime.lt.5) go to 9000

	do i=1,lentime
	  if(index(asep,btime(i:i)).ne.0) go to 20
	end do
	go to 9000

! -- The first integer is extracted from the string. This represents hours.

20      if(i.eq.1) go to 9000
	call char2num(ifail,btime(1:i-1),hh)
	if(ifail.ne.0) go to 9000
	if((hh.lt.0).or.(hh.gt.23)) go to 9000

	i=i+1
	if(lentime-i.lt.2) go to 9000
	do j=i,lentime
	  if(index(asep,btime(j:j)).ne.0) go to 40
	end do
	go to 9000

! -- The second integer (representing minutes) is extracted from the string.

40      if(j.eq.i) go to 9000
	call char2num(ifail,btime(i:j-1),mm)
	if(ifail.ne.0) go to 9000
	if((mm.lt.0).or.(mm.gt.59)) go to 9000

! -- The third integer (representing seconds) is extracted from the string.

	j=j+1
	if(lentime-j.lt.0) go to 9000
	call char2num(ifail,btime(j:lentime),ss)
	if(ifail.ne.0) go to 9000
	if((ss.lt.0).or.(ss.gt.59)) go to 9000
	ifail=0
	return

9000    ifail=1
	return

end subroutine char2time
 