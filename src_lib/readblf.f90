!     Last change:  JD   11 Apr 2003    4:53 am
subroutine read_bore_list_file(ifail,aprompt,coord_check)

! -- Subroutine read_bore_list_file reads a bore listing file, checking
! -- it for errors.

! -- Arguments are as follows:-
!       ifail:       returned as zero unless an error condition is encountered
!       aprompt:     contains text used to prompt for filename
!       coord_check  if supplied as 'yes' a check is made that listed bores
!                    are cited in a previously-read bore coordinates file

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(in)           :: aprompt
	character (len=*), intent(in), optional :: coord_check

	integer                         :: iunit,ierr,iline,imsg,nbore,i,idup,nbb,ifail1
	character (len=20)              :: aline
	character (len=200)             :: amsg,atempf
	character (len=10)              :: atemp
	logical                         :: check

	if(associated(bore_list_id)) deallocate (bore_list_id)
	nullify (bore_list_id)
	check=.true.
	if(present(coord_check))then
	  if(coord_check.eq.'no')check=.false.
	end if    

	idup=0
	imessage=0
	ifail=0
10      write(6,'(a)',advance='no') aprompt
	read(5,'(a)') bore_list_file
	bore_list_file=adjustl(bore_list_file)
	if(index(eschar,bore_list_file(1:2)).ne.0) then
	  escset=1
	  return
	end if
	if(bore_list_file.eq.' ') go to 10
        nbb=len_trim(bore_list_file)
        call getfile(ifail1,bore_list_file,atempf,1,nbb)
        if(ifail1.ne.0) go to 10
        bore_list_file=atempf
	iunit=nextunit()
	open(unit=iunit,file=bore_list_file,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9990
	  write(amessage,20) trim(bore_list_file)
20        format('cannot open bore listing file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  go to 10
	end if
	write(initial_message,30) trim(bore_list_file)
30      format(' Errors in bore listing file ',a,' ----->')

! -- The bore listing file is first perused to count the number of entries.

	imessage=0
	iline=0
	nbore=0
	do
	  iline=iline+1
	  read(iunit,'(a)',err=9000,end=200) cline
	  call linesplit(ifail,1)
	  if(ifail.lt.0) cycle
!         if(ifail.gt.0) then
!           if(imessage.eq.0) call write_initial_message(leadspace='yes')
!           call num2char(iline,aline)
!           write(amessage,70) trim(aline)
!70          format('   Line ',a,': insufficient entries (minimum of one ', &
!           'entry required).')
!           call write_message(increment=1)
!           cycle
!         end if
	  nbore=nbore+1
	end do

200     if(nbore.eq.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,210)
210       format('   No bores read from file.')
	  call write_message
	  go to 9990
	end if
	if(imessage.ne.0) go to 9990

! -- The bore listing file is next re-read and data placed into memory.

	rewind(unit=iunit,err=9100)
	allocate(bore_list_id(nbore), stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,250)
250       format(' Cannot allocate sufficient memory to run program.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9990
	end if

	iline=0
	num_bore_list=0
	do
	  iline=iline+1
	  read(iunit,'(a)',end=365,err=9000) cline
	  call linesplit(ifail,1)
	  if(ifail.lt.0) cycle
	  num_bore_list=num_bore_list+1
	  call num2char(iline,aline)
	  write(amsg,280) trim(aline)
280       format('  Line ',a,':')
	  imsg=len_trim(amsg)+1
	  if(right_word(1)-left_word(1)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,310) amsg(1:imsg)
310         format(1x,a,'bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    bore_list_id(num_bore_list)='*'
	  else
	    bore_list_id(num_bore_list)=cline(left_word(1):right_word(1))
	    call casetrans(bore_list_id(num_bore_list),'hi')
	  end if
!!	  if(num_bore_list.ne.1)then
!!	    atemp=bore_list_id(num_bore_list)
!!	    if(atemp.eq.'*') go to 369
!!	    duplication: do i=1,num_bore_list-1
!!	      if(bore_list_id(i).eq.atemp) then
!               if(imessage.eq.0) call write_initial_message(leadspace='yes')
!               write(amessage,360) amsg(1:imsg), trim(atemp)
!360             format(1x,a,'bore ',a,' cited previously.')
!               call write_message(increment=1)
!!		idup=idup+1
!!		go to 367
!!	      end if
!!	    end do duplication
!!	  end if
369       continue
	  atemp=bore_list_id(num_bore_list)
	  if(atemp.eq.'*') go to 367
	  if(.not.check) go to 367
	  do i=1,num_bore_coord
	    if(bore_coord_id(i).eq.atemp) go to 367
	  end do
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,366) amsg(1:imsg),trim(atemp),trim(bore_coord_file)
366       format(1x,a,'no coordinates for bore ',a, &
	  ' found in bore coordinates file ',a,'.')
	  call write_message(increment=1)
367       continue
	end do

365     if(imessage.eq.0)then
	  call num2char(num_bore_list,aline)
	  write(amessage,370) trim(aline),trim(bore_list_file)
370       format('  - ',a,' bores read from bore listing file ',a)
	  call write_message
	  if(idup.ne.0)then
	    write(amessage,380)
380         format('    warning: duplicate listing of at least one bore.')
	    call write_message
	  end if
	  ifail=0
	  go to 9999
	end if
!       write(6,*)
	go to 9990

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(bore_list_file)
9010    format('Error reading line ',a,' of bore listing file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9100    write(amessage,9110) trim(bore_list_file)
9110    format(' Cannot rewind bore listing file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9990    ifail=1
9999    close(unit=iunit,iostat=ierr)
	imessage=0
	return

end subroutine read_bore_list_file
 