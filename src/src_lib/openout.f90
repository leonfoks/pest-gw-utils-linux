!     Last change:  JD   23 Dec 2000    8:30 pm
subroutine open_output_file(ifail,aprompt,outfile,outunit,file_format)

! -- Subroutine open_output_file opens a file for data output.

! -- Subroutine arguments are as follows:-
!       ifail:    returned as non-zero in case of failure
!       aprompt:  user prompt for filename
!       outfile:  name of output file
!       outunit:  unit number of output file
!       file_format:  supplied as either "formatted" or "unformatted" to open
!                     formatted or unformatted file

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                       :: ifail
	character (len=*)                          :: aprompt,outfile
	integer, intent(out)                       :: outunit
	character (len=*), intent(in), optional    :: file_format
	integer                                    :: ierr,nbb,ifail1
	logical                                    :: lexist,lopened
	character (len=1)                          :: aa,aformat
        character (len=200)                        :: atempf
        character (len=20)                         :: bin_unform

! -- The correct unformatted word is read.

        include 'bin_unform.inc'

! -- The user is prompted for the name of the output file.

        aformat='f'
        if(present(file_format)) then
	  if(file_format.eq.'formatted') then
	    aformat='f'
	  else if(file_format.eq.'unformatted')then
	    aformat='u'
	  else
	    call sub_error('OPEN_OUTPUT_FILE')
	  end if
	end if

	imessage=0
	ifail=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)')outfile
	if(outfile.eq.' ') go to 10
	outfile=adjustl(outfile)
	if(index(eschar,outfile(1:2)).ne.0) then
	  escset=1
	  return
	end if
        nbb=len_trim(outfile)
        call getfile(ifail1,outfile,atempf,1,nbb)
        if(ifail1.ne.0) go to 10
        outfile=atempf
	inquire(file=outfile,opened=lopened)
	if(lopened)then
	  write(amessage,30) trim(outfile)
30        format(' File ',a,' is already open  - try again.')
	  call write_message(increment=1)
	  go to 10
	end if
	inquire(file=outfile,exist=lexist)
!       if(lexist) then
!40        write(6,50,advance='no')
!50        format(' File already exists: overwrite it?  [y/n] ')
!         read(5,'(a)') aa
!         call casetrans(aa,'lo')
!         if((aa.eq.'e').or.(aa.eq.'n'))then
!           write(6,*)
!           go to 10
!         end if
!         if(aa.ne.'y') go to 40
!       end if

! -- The file is opened.

	outunit=nextunit()
        if(aformat.eq.'f')then
          open(unit=outunit,file=outfile,status='new',iostat=ierr)
          if(ierr.ne.0) then
            open(unit=outunit,file=outfile,status='replace',iostat=ierr)
          end if
        else
          open(unit=outunit,file=outfile,status='new',form=bin_unform,iostat=ierr)
          if(ierr.ne.0) then
            open(unit=outunit,file=outfile,status='replace',form=bin_unform,iostat=ierr)
          end if
        end if
	if(ierr.ne.0) then
	  if(imessage.gt.5) then
	    ifail=1
	    return
	  end if
	  write(amessage,100) trim(outfile)
100       format(' Unable to open file ',a,' - try again.')
	  call write_message(increment=1)
	  go to 10
	end if

	return

end subroutine open_output_file
 