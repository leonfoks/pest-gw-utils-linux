!     Last change:  JD   16 Dec 2000    0:53 am
subroutine read_settings(ifail,idate,iheader)

! -- Subroutine read_settings reads a settings file located in the current
! -- directory.

! -- Arguments are as follows:-
!      ifail:    returned as zero unless settings file cannot be read
!      idate:    returned as zero unless date format is incorrect in settings
!                file
!      iheader:  returned as zero unless header specifier is incorrect


	use defn
	use inter

	integer, intent(out)			:: ifail,idate,iheader
	integer                                 :: iunit,ierr,iequals,i
	character (len=40)                      :: aline

	ifail=0
	idate=0
        iheader=0
	datespec=0
        headerspec=' '

	iunit=nextunit()
	open(unit=iunit,file='settings.fig',status='old',iostat=ierr)
	if(ierr.ne.0) then
	  ifail=1
	  return
	end if

	do
	  read(iunit,'(a)',err=100,end=200) cline
	  call casetrans(cline,'lo')
	  iequals=index(cline,'=')
	  if(iequals.le.1) cycle
	  aline=cline(1:iequals-1)
	  aline=adjustl(aline)
	  if(aline(1:4).eq.'date') then
	    aline=cline(iequals+1:)
	    aline=adjustl(aline)
	    if((aline(1:2).eq.'dd').and.(aline(4:5).eq.'mm')) then
	      datespec=1
	    else if((aline(1:2).eq.'mm').and.(aline(4:5).eq.'dd')) then
	      datespec=2
	    else
	      idate=1
	    end if
	  else if(aline(1:6).eq.'colrow') then
	    aline=cline(iequals+1:)
            do i=1,len_trim(aline)
              if(aline(i:i).eq.'''')aline(i:i)=' '
            end do
	    aline=adjustl(aline)
            call casetrans(aline,'lo')
            if(aline(1:3).eq.'yes')then
              headerspec='yes'
            else if(aline(1:2).eq.'no')then
              headerspec='no '
            else
              iheader=1
            end if
          end if
	end do
	go to 200

100	ifail=2
200     close(unit=iunit,iostat=ierr)
	return

end subroutine read_settings

 