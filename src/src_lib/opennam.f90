!     Last change:  JD   23 Dec 2000    8:32 pm
subroutine open_named_input_file(ifail,aprompt,infile,inunit)

! -- Subroutine open_named_input_file opens an input file for which a default
!    name exists.

! -- Arguments are as follows:-
!       ifail:    returned as non-zero in case of failure
!       aprompt:  user prompt for filename
!       infile:   name of input file (carries default name on entry)
!       inunit:   unit number of input file

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(in)           :: aprompt
	character (len=*), intent(inout)        :: infile
	integer, intent(out)                    :: inunit
	integer                                 :: ierr,nbb,ifail1
	logical                                 :: lopened
	character (len=200)                     :: tempfile,atempf

! -- The user is prompted for the name of the file to open.

	imessage=0
	ifail=0
5       if(infile.eq.' ')then
10        write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	  read(5,'(a)') tempfile
	  if(tempfile.eq.' ') go to 10
	  tempfile=adjustl(tempfile)
	  if(index(eschar,tempfile(1:2)).ne.0) then
	    escset=1
	    return
	  end if
          nbb=len_trim(tempfile)
          call getfile(ifail1,tempfile,atempf,1,nbb)
          if(ifail1.ne.0) go to 10
          tempfile=atempf
	else
	  write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)-1)//&
	  ' ['//trim(infile)//']: '
	  read(5,'(a)') tempfile
	  if(tempfile.eq.' ')then
	    tempfile=infile
	  else
	    tempfile=adjustl(tempfile)
	    if(index(eschar,tempfile(1:2)).ne.0) then
	      escset=1
	      return
	    end if
            nbb=len_trim(tempfile)
            call getfile(ifail1,tempfile,atempf,1,nbb)
            if(ifail1.ne.0) go to 5
            tempfile=atempf
	  end if
	end if

! -- Is the file already open?

	inquire(file=tempfile,opened=lopened)
	if(lopened)then
	  write(amessage,30) trim(tempfile)
30        format(' File ',a,' is already open  - try again.')
	  call write_message(increment=1)
	  go to 5
	end if

! -- The file is opened.

	inunit=nextunit()
	open(unit=inunit,file=tempfile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.gt.5) then
	    ifail=1
	    return
	  end if
	  write(amessage,50) trim(tempfile)
50        format(' Cannot open file ',a,'  - try again.')
	  call write_message(increment=1)
	  if(infile.eq.tempfile) infile=' '
	  go to 5
	end if
	infile=tempfile

	return

end subroutine open_named_input_file
 