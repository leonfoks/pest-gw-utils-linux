!     Last change:  J    20 Dec 2004   12:46 pm
subroutine write_real_array(ifail,aprompt,array,pm_header,rows,columns,&
binary_header,atype,ntrans,kstp,kper,pertim,totim,text,ncol,nrow,ilay,istage, &
realfile,aaformat)

! -- Subroutine writ_real_array writes a real array to an output file.

! -- Arguments are as follows:-
!       ifail:         returned as non-zero if error condition arises
!       aprompt:       prompt requesting name of output file (contains output
!                      filename on exit)
!       array:         real array for output
!       pm_header:     if "yes" put number of columns and rows at start of file
!       rows,columns:  number of rows and columns in finite-difference grid
!       binary_header: if "yes" put MODFLOW/MT3D-type header at start of binary
!                      output file
!       atype:         determines whether MODFLOW or MT3D header
!       ntrans,kstp,kper,pertim,totim,text,ncol,nrow,ilay: used in header
!                                                          (if required)
!       istage:        used if making a different call to this subroutine to open
!                      real array file and to write it
!       realfile:      name of real array file if two separate calls to this subroutine
!       aaformat:      format status if two separate calls to this subroutine

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(inout)        :: aprompt
	real, intent(in),dimension(:,:)         :: array
	character (len=*), intent(in), optional :: pm_header
	integer, intent(in), optional           :: rows,columns
	character (len=*), optional             :: binary_header
	character (len=*), optional             :: atype
	integer, optional                       :: ntrans,kstp,kper
	real, optional                          :: pertim,totim
	character (len=16),optional             :: text
	integer, optional                       :: ncol,nrow,ilay
        integer, intent(in), optional           :: istage
        character (len=*), optional             :: realfile
        character (len=*), optional             :: aaformat
	character (len=120)                     :: afile,bfile
	character (len=4)                       :: extension
	logical                                 :: pm,lexist,lopened,mhead
	character (len=1)                       :: aformat,aa
	integer                                 :: i,iunit,nbb,ifail1,iistage,ierr
        integer                                 :: kkstp,kkper,iilay,nnrow,nncol
        real                                    :: ppertim,ttotim
        character (len=16)                      :: ttext

	ifail=0
	pm=.false.

	if(present(pm_header)) then
	  if(pm_header.ne.'yes') go to 5
	  if((.not.(present(rows))).or.(.not.(present(columns)))) &
	  call sub_error('WRITE_REAL_ARRAY')
	  pm=.true.
	end if

	mhead=.false.
	if(present(binary_header))then
	  if(binary_header.eq.'yes') then
	    if((.not.(present(atype))).or.(.not.(present(ntrans))).or. &
	       (.not.(present(kstp))).or.(.not.(present(kper))).or. &
	       (.not.(present(pertim))).or.(.not.(present(totim))).or. &
	       (.not.(present(text))).or.(.not.(present(ncol))).or. &
	       (.not.(present(nrow))).or.(.not.(present(ilay)))) &
	       call sub_error('WRITE_REAL_ARRAY')
	       mhead=.true.
	  else if(binary_header.ne.'no') then
	    call sub_error('WRITE_REAL_ARRAY')
	  end if
	end if

5       imessage=0

        if(present(istage)) then
          if((.not.(present(realfile))).or.(.not.(present(aaformat)))) &
          call sub_error('WRITE_REAL_ARRAY')
          iistage=istage
        else
          iistage=0
        end if
        if(iistage.eq.2) then
          afile=realfile
          aformat=aaformat
          go to 350
        end if

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
!       inquire(file=afile,exist=lexist)
!       if(lexist) then
!55        if(aprompt(2:2).ne.' ') then
!           write(6,56,advance='no')
!56          format(' File already exists: overwrite it?  [y/n]: ')
!         else
!           write(6,57,advance='no')
!57          format('   file already exists: overwrite it?  [y/n]: ')
!         end if
!         read(5,'(a)') aa
!         if(aa.eq.' ') go to 55
!         call casetrans(aa,'lo')
!         if((index(eschar,aa).ne.0).or.(aa.eq.'n')) then
!           write(6,*)
!           go to 10
!         end if
!         if(aa.ne.'y') go to 55
!       end if

        if(iistage.eq.1)then
          realfile=afile
          aaformat=aformat
          return
        end if

350     continue

	iunit=nextunit()
	if(aformat.eq.'f') then
	  open(unit=iunit,file=afile,status='replace',err=100)
	  if(pm) then
	    write(iunit,60,err=700) columns, rows
60          format(1x,i7,1x,i7)
	  end if
	  do i=1,rows
	    write(iunit,65,err=700) array(:,i)
65          format(7(1pe14.6))
	  end do
	else
          open(unit=iunit,file=afile,form='binary',iostat=ierr)
          if(ierr.ne.0)then
            open(unit=iunit,file=afile,form='unformatted',err=100)
	  end if
	  if(mhead)then
	    if(atype.eq.'modflow')then
	      write(iunit,err=700) kstp,kper,pertim,totim,text,ncol,nrow,ilay
	    else if(atype.eq.'mt3d')then
	      write(iunit,err=700) ntrans,kstp,kper,totim,text,ncol,nrow,ilay
	    else
	      call sub_error('WRITE_REAL_ARRAY')
	    end if
	  else
            kkstp=1
            kkper=1
            ppertim=1.0
            ttotim=1.0
            ttext=' '
            iilay=1
            nncol=columns
            nnrow=rows
!	    write(iunit,err=700)kkstp,kkper,ppertim,ttotim,ttext,nncol,nnrow,iilay
	    write(iunit,err=700)
	  end if
	  write(iunit,err=700) array
	end if
	if(aprompt(2:2).ne.' ') then
	  if(aformat.eq.'f') then
	    write(amessage,70) trim(afile)
70          format('  - formatted real array written to file ',a)
	  else
	    write(amessage,71) trim(afile)
71          format('  - unformatted real array written to file ',a)
	  end if
	else
	  if(aformat.eq.'f')then
	    write(amessage,80) trim(afile)
80          format('    - formatted real array written to file ',a)
	  else
	    write(amessage,81) trim(afile)
81          format('    - unformatted real array written to file ',a)
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

end subroutine write_real_array
 
