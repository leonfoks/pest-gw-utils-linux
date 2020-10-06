subroutine write_message(increment,iunit,error,leadspace,endspace)

! -- Subroutine write_message formats and writes a message.

! -- Arguments are as follows:-
!       increment:  the increment to the message counter
!       iunit:      the unit number to which the message is written
!       error:      if "yes" precede message with "Error"
!       leadspace   if "yes" precede message with blank line
!       endspace    if "yes" follow message by blank line

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(in), optional           ::increment,iunit
	integer                                 ::jend,i,nblc,junit,leadblank
	integer                                 ::itake,j
	character (len=*), intent(in), optional ::error,leadspace,endspace
	character (len=20) ablank

	ablank=' '
	itake=0
	j=0
	if(present(increment)) imessage=imessage+increment
	if(present(iunit))then
	  junit=iunit
	else
	  junit=6
	end if
	if(present(leadspace))then
	  if(leadspace.eq.'yes') write(junit,*)
	endif
	if(present(error))then
	  if(index(error,'yes').ne.0)then
	    nblc=len_trim(amessage)
	    amessage=adjustr(amessage(1:nblc+8))
	    if(nblc+8.lt.len(amessage)) amessage(nblc+9:)=' '
	    amessage(1:8)=' Error: '
	  end if
	end if

	do i=1,20
	  if(amessage(i:i).ne.' ')exit
20      end do
	leadblank=i-1
	nblc=len_trim(amessage)
5       jend=j+78-itake
	if(jend.ge.nblc) go to 100
	do i=jend,j+1,-1
	if(amessage(i:i).eq.' ') then
	  if(itake.eq.0) then
	     write(junit,'(a)',err=200) amessage(j+1:i)
	     itake=2+leadblank
	  else
	     write(junit,'(a)',err=200) ablank(1:leadblank+2)//amessage(j+1:i)
	  end if
	  j=i
	  go to 5
	end if
	end do
	if(itake.eq.0)then
	  write(junit,'(a)',err=200) amessage(j+1:jend)
	  itake=2+leadblank
	else
	  write(junit,'(a)',err=200) ablank(1:leadblank+2)//amessage(j+1:jend)
	end if
	j=jend
	go to 5
100     jend=nblc
	if(itake.eq.0)then
	write(junit,'(a)',err=200) amessage(j+1:jend)
	  else
	write(junit,'(a)',err=200) ablank(1:leadblank+2)//amessage(j+1:jend)
	  end if
	if(present(endspace))then
	  if(endspace.eq.'yes') write(junit,*)
	end if
	return

200     call exit(100)

end subroutine write_message
 