!     Last change:  J    20 Dec 2004    1:04 pm
subroutine read_real_array(ifail,aprompt,array,pm_header,rows,columns)

! -- Subroutine read_real_array reads a real array from a file.

! -- Arguments are as follows:-
!       ifail:         returned as non-zero if error condition is encountered.
!       aprompt:       prompt text for real array filename
!       array:         real array read from file
!       pm_header:     if "yes" number of columns and rows expected as first
!                      line of input file.
!       rows,columns:  number of rows and columns in finite-difference grid

! -- Revision history:-
!       June-November, 1995:- version 1.

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(inout)        :: aprompt
	real, intent(out),dimension(:,:)        :: array
	character (len=*), intent(in), optional :: pm_header
	integer, intent(in), optional           :: rows,columns
	character (len=120)                     :: afile,bfile
	character (len=4)                       :: extension
	logical                                 :: pm
	character (len=1)                       :: aformat
	integer                                 :: i,iunit,nncol,nnrow,nbb,ifail1
        integer                                 :: kstp,kper,ilay
        real                                    :: pertim,totim
        character (len=16)                      :: text

	ifail=0
	pm=.false.
	if(present(pm_header)) then
	  if(pm_header.ne.'yes') go to 5
	  if((.not.(present(rows))).or.(.not.(present(columns)))) &
	  call sub_error('READ_REAL_ARRAY')
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
	  if(extension.eq.'.reu')then
	    aformat='u'
	  else if(extension.eq.'.ref')then
	    aformat='f'
	  else
20          if(aprompt(2:2).ne.' ') then
30            write(6,40,advance='no')
40            format(' Is this a formatted or unformatted file? [f/u]: ')
	    else
	      write(6,50,advance='no')
50            format('   is this a formatted or unformatted file? [f/u]: ')
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
!	  read(iunit,end=600) kstp,kper,pertim,totim,text,nncol,nnrow,ilay
	  read(iunit,end=600)
	  read(iunit,err=550,end=600) array
	end if
	if(aprompt(2:2).ne.' ') then
	  write(amessage,70) trim(afile)
70        format('  - real array read from file ',a)
	else
	  write(amessage,80) trim(afile)
80        format('    - real array read from file ',a)
	end if
	call write_message
	aprompt=afile
	go to 999

100     if(imessage.eq.5) go to 995
	if(pm)then
	  write(amessage,110) trim(afile)
110       format('cannot open formatted real array file ',&
	  a,' - try again.')
	else
	  write(amessage,111) trim(afile)
111       format('cannot open formatted real array file ',a,' - try again.')
	end if
	call write_message(error='yes',increment=1,endspace='yes')
	go to 10

120     if(imessage.eq.5) go to 995
	write(amessage,130) trim(afile)
130     format('cannot open unformatted real array file ',a,' - try again.')
	call write_message(error='yes',increment=1,endspace='yes')
	go to 10

150     write(amessage,160) trim(afile)
160     format('unexpected end to formatted real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

200     write(amessage,210) trim(afile)
210     format('column and row numbers are expected on first line of ',&
	'real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

250     write(amessage,260) trim(afile)
260     format('the row and/or column number provided in header to ',&
	'real array file ',a,' does not match that provided in ',&
	'grid specification file.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

300     write(amessage,310) trim(afile)
310     format('cannot read row and/or column number from first line of ',&
	'formatted real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

350     write(amessage,360) trim(afile)
360     format(' Data error encountered while reading formatted real ',&
	'array file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 995

400     write(amessage,410) trim(afile)
410     format('unexpected end to formatted real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

450     write(amessage,460) trim(afile)
460     format(' Data error encountered while reading formatted real ',&
	'array file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 995

500     write(amessage,510) trim(afile)
510     format('unexpected end to formatted real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

550     write(amessage,560) trim(afile)
560     format(' Error encountered while reading unformatted real ',&
	'array file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 995

600     write(amessage,610) trim(afile)
610     format('unexpected end to unformatted real array file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 995

995     ifail=1
999     close(unit=iunit,iostat=i)
	return

end subroutine read_real_array
 
