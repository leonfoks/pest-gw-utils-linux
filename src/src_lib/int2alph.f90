subroutine int2alph(inum,alph,nsig)

! -- Subroutine int2alph converts a number to a string of letters in counting
!    sequence with base 26.

! -- Arguments are as follows:-
!       inum:  the number to be converted
!       alph:  the alphanumeric output string
!       nsig:  an optional externally-supplied field width

! -- Revision history:-
!       June-November, 1995: version 1.

	implicit none

	integer, intent(in)             :: inum
	character (len=*), intent(out)  :: alph
	integer, optional, intent(in)   :: nsig

	integer                         :: llen,nnsig,itemp1,itemp2,i,iinum

	alph=' '
	iinum=inum
	llen=len(alph)
	if(present(nsig)) then
	  if(nsig.gt.llen) call sub_error('INT2ALPH')
	  nnsig=nsig
	else
	  nnsig=llen
	end if

	do i=1,nnsig
	  itemp1=(iinum-1)/26
	  itemp2=iinum-itemp1*26
	  iinum=itemp1
	  alph(nnsig-i+1:nnsig-i+1)=achar(itemp2+96)
	  if(iinum.eq.0) go to 100
	end do
	do i=1,nnsig
	  alph(nnsig-i+1:nnsig-i+1)='*'
	end do

100     alph=adjustl(alph)
	return

end subroutine int2alph
 