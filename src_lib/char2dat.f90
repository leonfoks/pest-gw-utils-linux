subroutine char2date(ifail,adate,dd,mm,yy)

! -- Subroutine CHAR2DATE extracts the date from a string.


! -- Arguments are as follows:-
!      ifail:      returns a non-zero value if an error condition is encountered
!      adate:      the string containing the date
!      dd,mm,yy    the day, month and year read from the date string

! --  Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)    :: ifail
	character (len=*), intent(in)   :: adate
	integer, intent(out) :: dd,mm,yy
	integer :: lendate,i,j
	character (len=2)       :: asep
	character (len=20)      :: bdate

	ifail=0
	asep=':/'
	if(adate.eq.' ') go to 9000
	bdate=adjustl(adate)
	lendate=len_trim(bdate)
	if(lendate.lt.8) go to 9000

	do i=1,lendate
	  if(index(asep,bdate(i:i)).ne.0) go to 20
	end do
	go to 9000

! -- The first integer is extracted from the date string. This is either days
!    or months depending on the contents of file settings.fig.

20      if(i.eq.1) go to 9000
	if(datespec.ne.1) then
	   call char2num(ifail,bdate(1:i-1),mm)
	else
	   call char2num(ifail,bdate(1:i-1),dd)
	end if
	if(ifail.ne.0) go to 9000

	i=i+1
	if(lendate-i.lt.5) go to 9000
	do j=i,lendate
	  if(index(asep,bdate(j:j)).ne.0) go to 40
	end do
	go to 9000

! -- The second integer is extracted from the date string. This is either months
!    or days depending on the contents of file settings.fig.

40      if(j.eq.i) go to 9000
	if(datespec.ne.1) then
	  call char2num(ifail,bdate(i:j-1),dd)
	else
	  call char2num(ifail,bdate(i:j-1),mm)
	end if
	if(ifail.ne.0) go to 9000
	if((dd.le.0).or.(dd.gt.31)) go to 9000
	if((mm.le.0).or.(mm.gt.12)) go to 9000
	if(dd.eq.31)then
	  if((mm.eq.2).or.(mm.eq.4).or.(mm.eq.6).or.(mm.eq.9).or.&
	  (mm.eq.11)) go to 9000
	end if
	if((mm.eq.2).and.(dd.eq.30)) go to 9000

! -- The third integer is extracted from the date string. This is years.

	j=j+1
!	if(lendate-j.ne.3) go to 9000
	call char2num(ifail,bdate(j:lendate),yy)
	if(ifail.ne.0) go to 9000
	if(.not.leap(yy))then
	  if((mm.eq.2).and.(dd.eq.29)) go to 9000
	end if
	ifail=0
	return

9000    ifail=1
	return

end subroutine char2date
 