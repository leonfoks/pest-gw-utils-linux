!     Last change:  JD   11 Apr 2003    8:38 am
subroutine open_input_file(ifail,aprompt,infile,inunit,file_format,form_prompt,fformat)

! -- Subroutine open_input_file opens an input file.

! -- Arguments are as follows:-
!       ifail:        returned as non-zero if failure
!       aprompt:      user prompt for filename
!       infile:       name of input file
!       inunit:       unit number to communicate with file
!       file_format:  supplied as either "formatted" or "unformatted" to open
!                     formatted or unformatted file
!       form_prompt:  supplied as either "yes" or "no" to prompt for file format
!       fformat:      returns file format to main program if form_prompt is supplied

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                       :: ifail
	character (len=*), intent(in)              :: aprompt
	character (len=*), intent(out)             :: infile
	integer, intent(out)                       :: inunit
	character (len=*), intent(in), optional    :: file_format
        character (len=*), intent(in), optional    :: form_prompt
        character (len=*), intent(out),optional    :: fformat
	integer                                    :: ierr,nbb,ifail1
	character (len=1)                          :: aformat,aa
	character (len=10)                         :: aaformat
        character (len=120)                        :: bfile
	logical                                    :: lopened
        character (len=20)                         :: bin_unform

! -- The correct unformatted word is read.

        include 'bin_unform.inc'

! -- The user is prompted for the name of the file to open.

	imessage=0
	ifail=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)') infile
	if(infile.eq.' ') go to 10
	infile=adjustl(infile)
	if(index(eschar,infile(1:2)).ne.0) then
	  escset=1
	  return
	end if
        nbb=len_trim(infile)
        call getfile(ifail1,infile,bfile,1,nbb)
        if(ifail1.ne.0) go to 10
        infile=bfile


! -- Is the file already open?

	inquire(file=infile,opened=lopened)
	if(lopened)then
	  write(amessage,30) trim(infile)
30        format(' File ',a,' is already open  - try again.')
	  call write_message(increment=1)
	  go to 10
	end if

! -- The file is opened (after the format is established).

	aformat='f'
        if(present(form_prompt))then
          if(form_prompt.eq.'yes')then
            if(.not.(present(fformat))) call sub_error('OPEN_INPUT_FILE')
31          write(6,32,advance='no')
32          format(' Is this a formatted or unformatted file? [f/u]: ')
            read(5,'(a)') aa
            call casetrans(aa,'lo')
            if(aa.eq.'e') go to 10
            if(aa.eq.'f')then
              aformat='f'
              fformat='f'
            else if(aa.eq.'u')then
              aformat='u'
              fformat='u'
            else
              go to 31
            end if
          end if
        else if(present(file_format)) then
	  if(file_format.eq.'formatted') then
	    aformat='f'
	  else if(file_format.eq.'unformatted')then
	    aformat='u'
	  else
	    call sub_error('OPEN_INPUT_FILE')
	  end if
	end if
	inunit=nextunit()
	if(aformat.eq.'f')then
	  open(unit=inunit,file=infile,status='old',iostat=ierr)
	else
	  open(unit=inunit,file=infile,status='old',form=bin_unform,  &
          iostat=ierr)
	end if

! -- If the file could not be opened a test is made as to whether the file
!    is formatted or unformatted when the opposite was expected.

	if(ierr.ne.0)then
	  inquire(file=infile,unformatted=aaformat,iostat=ierr)
	  if(ierr.ne.0)go to 40
	  call casetrans(aaformat,'lo')
	  if(aaformat.eq.' ') go to 40
	  if(aformat.eq.'f')then
	    if(aaformat.eq.'yes') then
	      if(imessage.gt.5) then
		ifail=1
		return
	      end if
	      write(amessage,35) trim(infile)
35            format(' File ',a,' is an unformatted file  - try again.')
	      call write_message(increment=1)
	      go to 10
	    else
	      go to 40
	    end if
	  else
	    if(aaformat.eq.'no') then
	      if(imessage.gt.5) then
		ifail=1
		return
	      end if
	      write(amessage,36) trim(infile)
36            format(' File ',a,' is a formatted file  - try again.')
	      call write_message(increment=1)
	      go to 10
	    else
	      go to 40
	    end if
	  end if
40        if(imessage.gt.5) then
	    ifail=1
	    return
	  end if
	  write(amessage,50) trim(infile)
50        format(' Cannot open file ',a,'  - try again.')
	  call write_message(increment=1)
	  go to 10
	end if

	return

end subroutine open_input_file
 