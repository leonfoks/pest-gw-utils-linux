!     Last change:  JD    9 May 2003   11:43 am
subroutine read_integer_array(ifail,aprompt,array,pm_header,rows,columns,defaultfile)

! -- Subroutine read_integer_array reads an integer array from a file.

! -- Arguments are as follows:-
!       ifail:         returned as non-zero if error condition is encountered.
!       aprompt:       prompt text for integer array filename
!       array:         integer array read from file
!       pm_header:     if "yes" number of columns and rows expected as first
!                      line of input file.
!       rows,columns:  number of rows and columns in finite-difference grid
!       defaultfile:   a default name for the integer array file


	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(inout)        :: aprompt
	integer, intent(out),dimension(:,:)     :: array
	character (len=*), intent(in), optional :: pm_header
	integer, intent(in), optional           :: rows,columns
        character (len=*), optional             :: defaultfile
	character (len=120)                     :: afile,bfile
	character (len=4)                       :: extension
	logical                                 :: pm
	character (len=1)                       :: aformat
	integer                                 :: i,iunit,nncol,nnrow,ifail1,nbb,idef

	ifail=0
	pm=.false.
	if(present(pm_header)) then
	  if(pm_header.ne.'yes') go to 5
	  if((.not.(present(rows))).or.(.not.(present(columns)))) &
	  call sub_error('READ_INTEGER_ARRAY')
	  pm=.true.
	end if

        idef=0
        if(present(defaultfile))then
          if(defaultfile.ne.' ')then
            idef=1
            aprompt=trim(aprompt)//' ['//trim(defaultfile)//']:'
          end if
        end if

5       imessage=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)')afile
	if(afile.eq.' ') then
          if(idef.eq.0) then
            go to 10
          else
            afile=defaultfile
            go to 15
          end if
        end if
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0)then
	  escset=1
	  return
	end if
        nbb=len_trim(afile)
        call getfile(ifail1,afile,bfile,1,nbb)
        if(ifail1.ne.0) go to 10
        afile=bfile
15      continue
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
30          write(6,40,advance='no')
40          format(' Is this a formatted or unformatted file? [f/u]: ')
	    read(5,'(a)') aformat
	    if((aformat.eq.'e').or.(aformat.eq.'E')) then
	      write(6,*)
	      go to 10
	    end if
	    if(aformat.eq.' ') go to 30
	    call casetrans(aformat,'lo')
	    if((aformat.ne.'f').and.(aformat.ne.'u')) go to 30
	  end if
!       end if

	iunit=nextunit()
	if(aformat.eq.'f') then
	  open(unit=iunit,file=afile,status='old',err=100)
	  if(pm) then
	    read(iunit,'(a)',end=150) cline
	    call linesplit(ifail,2)
	    if(ifail.ne.0) go to 200
	    nncol=char2int(ifail,1)
	    if(ifail.ne.0) go to 300
	    nnrow=char2int(ifail,2)
	    if(ifail.ne.0) go to 300
	    if((nncol.ne.columns).or.(nnrow.ne.rows)) go to 250
	    read(iunit,*,err=350,end=400) array
	  else
	    read(iunit,*,err=450,end=500) array
	  end if
	else
	  open(unit=iunit,file=afile,status='old',form='binary',err=120)
	  read(iunit,end=600)
	  read(iunit,err=550,end=600) array
	end if
	write(amessage,70) trim(afile)
70      format('  - integer array read from file ',a)
	call write_message
	aprompt=afile
	go to 999

100     if(imessage.eq.5) go to 995
	if(pm)then
	  write(amessage,110) trim(afile)
110       format('cannot open formatted integer array file ',&
	  a,' - try again.')
	else
	  write(amessage,111) trim(afile)
111       format('cannot open formatted integer array file ',a,' - try again.')
	end if
	call write_message(error='yes',increment=1,endspace='yes')
	go to 10

120     if(imessage.eq.5) go to 995
	write(amessage,130) trim(afile)
130     format('cannot open unformatted integer array file ',a,' - try again.')
	call write_message(error='yes',increment=1,endspace='yes')
	go to 10

150     write(amessage,160) trim(afile)
160     format('unexpected end to formatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

200     write(amessage,210) trim(afile)
210     format('column and row numbers are expected on first line of ',&
	'formatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

250     write(amessage,260) trim(afile)
260     format('the row and/or column number provided in header to ',&
	'integer array file ',a,' does not match that provided in ',&
	'grid specification file.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

300     write(amessage,310) trim(afile)
310     format('cannot read row and/or column number from first line of ',&
	'formatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

350     write(amessage,360) trim(afile)
360     format(' Data error encountered while reading formatted integer ',&
	'array file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 995

400     write(amessage,410) trim(afile)
410     format('unexpected end to formatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

450     write(amessage,460) trim(afile)
460     format(' Data error encountered while reading formatted integer ',&
	'array file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 995

500     write(amessage,510) trim(afile)
510     format('unexpected end to formatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

550     write(amessage,560) trim(afile)
560     format(' Error encountered while reading unformatted integer ',&
	'array file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 995

600     write(amessage,610) trim(afile)
610     format('unexpected end to unformatted integer array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

995     ifail=1
999     close(unit=iunit,iostat=i)
	return

end subroutine read_integer_array
 