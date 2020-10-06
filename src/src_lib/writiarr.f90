!     Last change:  JD    9 May 2003   11:55 am
subroutine write_integer_array(ifail,aprompt,array,pm_header,rows,columns)

! -- Subroutine write_integer_array writes an integer array to an output file.

! -- Arguments are as follows:-
!       ifail:         returned as non-zero if error condition arises
!       aprompt:       prompt requesting name of output file (contains output
!                      filename on exit)
!       array:         real array for output
!       pm_header:     if "yes" put number of columns and rows at start of file
!       rows,columns:  number of rows and columns in finite-difference grid


	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(inout)        :: aprompt
	integer, intent(in),dimension(:,:)      :: array
	character (len=*), intent(in), optional :: pm_header
	integer, intent(in), optional           :: rows,columns
	character (len=120)                     :: afile,bfile
	character (len=4)                       :: extension
	logical                                 :: pm,lexist,lopened
	character (len=1)                       :: aformat,aa
	integer                                 :: i,iunit,ifail1,nbb

	ifail=0
	pm=.false.
	if(present(pm_header)) then
	  if(pm_header.ne.'yes') go to 5
	  if((.not.(present(rows))).or.(.not.(present(columns)))) &
	  call sub_error('WRITE_INTEGER_ARRAY')
	  pm=.true.
	end if

5       imessage=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)')afile
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0)then
	  escset=1
	  return
	end if
	if(afile.eq.' ') go to 10
        nbb=len_trim(afile)
        call getfile(ifail1,afile,bfile,1,nbb)
        if(ifail1.ne.0) go to 10
        afile=bfile

!       if(pm)then
!         aformat='f'
!       else
	  i=len_trim(afile)
	  extension=afile(max(i-3,1):i)
	  call casetrans(extension,'lo')
	  if(extension.eq.'.inu')then
	    aformat='u'
	  else if(extension.eq.'.inf')then
	    aformat='f'
	  else
20          if(aprompt(2:2).ne.' ') then
30            write(6,40,advance='no')
40            format(' Write a formatted or unformatted file? [f/u]: ')
	    else
	      write(6,50,advance='no')
50            format('   write a formatted or unformatted file? [f/u]: ')
	    end if
	    read(5,'(a)') aformat
	    if((aformat.eq.'e').or.(aformat.eq.'E')) then
	      write(6,*)
	      go to 10
	    end if
	    if(aformat.eq.' ') go to 20
	    call casetrans(aformat,'lo')
	    if((aformat.ne.'f').and.(aformat.ne.'u')) go to 20
	  end if
!       end if

	inquire (file=afile,opened=lopened)
	if(lopened) then
	  if(aprompt(2:2).ne.' ') then
	    write(6,51)
51          format(' file is already open  - try again')
	  else
	    write(6,52)
52          format('   file is already open  - try again')
	  end if
	  go to 10
	end if
!	inquire(file=afile,exist=lexist)
!	if(lexist) then
!55        if(aprompt(2:2).ne.' ') then
!	    write(6,56,advance='no')
!56          format(' File already exists: overwrite it?  [y/n]: ')
!	  else
!	    write(6,57,advance='no')
!57          format('   file already exists: overwrite it?  [y/n]: ')
!	  end if
!	  read(5,'(a)') aa
!	  if(aa.eq.' ') go to 55
!	  call casetrans(aa,'lo')
!	  if((index(eschar,aa).ne.0).or.(aa.eq.'n')) then
!	    write(6,*)
!	    go to 10
!	  end if
!	  if(aa.ne.'y') go to 55
!	end if

	iunit=nextunit()
	if(aformat.eq.'f') then
	  open(unit=iunit,file=afile,status='replace',err=100)
	  if(pm) then
	    write(iunit,60,err=700) columns, rows
60          format(1x,i7,1x,i7)
	  end if
	  do i=1,rows
	    write(iunit,65,err=700) array(:,i)
65          format(20i5)
	  end do
	else
	  open(unit=iunit,file=afile,form='binary',status='replace',err=100)
	  write(iunit,err=700)
	  write(iunit,err=700) array
	end if
	if(aprompt(2:2).ne.' ') then
	  if(aformat.eq.'f') then
	    write(amessage,70) trim(afile)
70          format('  - formatted integer array written to file ',a)
	  else
	    write(amessage,71) trim(afile)
71          format('  - unformatted integer array written to file ',a)
	  end if
	else
	  if(aformat.eq.'f')then
	    write(amessage,80) trim(afile)
80          format('    - formatted integer array written to file ',a)
	  else
	    write(amessage,81) trim(afile)
81          format('    - unformatted integer array written to file ',a)
	  end if
	end if
	call write_message
	aprompt=afile
	go to 999

100     write(amessage,110) trim(afile)
110     format('cannot open file ',a,' for output.')
	call write_message(error='yes',leadspace='yes')
	go to 995
700     write(amessage,710) trim(afile)
710     format(' Cannot write data to file ',a, &
	': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 995

995     ifail=1
999     close(unit=iunit,iostat=i)
	return

end subroutine write_integer_array
 