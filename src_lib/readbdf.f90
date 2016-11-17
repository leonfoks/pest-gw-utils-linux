subroutine read_bore_diff_file(ifail,aprompt,coord_check)

! -- Subroutine read_bore_diff_file reads a bore difference listing file, checking
! -- it for errors.

! -- Arguments are as follows:-
!       ifail:       returned as zero unless an error condition is encountered
!       aprompt:     contains text used to prompt for filename
!       coord_check  if supplied as 'yes' a check is made that listed bores
!                    are cited in a previously-read bore coordinates file

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(in)           :: aprompt
	character (len=*), intent(in), optional :: coord_check

	integer                         :: iunit,ierr,iline,imsg,nbore,i,nbb,ifail1
	character (len=20)              :: aline
	character (len=200)             :: amsg,atempf
	character (len=10)              :: atemp
	logical                         :: check

	if(associated(bore_list_id)) deallocate (bore_list_id)
	nullify (bore_list_id)
        if(associated(bore_diff_id)) deallocate (bore_diff_id)
        nullify (bore_diff_id)
	check=.true.
	if(present(coord_check))then
	  if(coord_check.eq.'no')check=.false.
	end if    

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
20        format('cannot open bore difference listing file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  go to 10
	end if
	write(initial_message,30) trim(bore_list_file)
30      format(' Errors in bore difference listing file ',a,' ----->')

! -- The bore difference listing file is first perused to count the number of entries.

	imessage=0
	iline=0
	nbore=0
	do
	  iline=iline+1
	  read(iunit,'(a)',err=9000,end=200) cline
          if(cline.eq.' ') cycle
	  call linesplit(ifail,3)
          if(ifail.gt.0) then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(iline,aline)
            write(amessage,70) trim(aline)
 70          format('   Line ',a,': insufficient entries (three ', &
            'entries required).')
            call write_message(increment=1)
            cycle
          end if
	  nbore=nbore+1
	end do
        if(imessage.ne.0) go to 9990

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
	allocate(bore_list_id(2*nbore),bore_diff_id(nbore),stat=ierr)
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
          if(cline.eq.' ') cycle
	  call linesplit(ifail,3)
	  num_bore_list=num_bore_list+1
	  call num2char(iline,aline)
	  write(amsg,280) trim(aline)
280       format('  Line ',a,':')
	  imsg=len_trim(amsg)+1
	  if(right_word(1)-left_word(1)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,310) amsg(1:imsg)
310         format(1x,a,'first bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    bore_list_id(num_bore_list)='*'
	  else
	    bore_list_id(num_bore_list)=cline(left_word(1):right_word(1))
	    call casetrans(bore_list_id(num_bore_list),'hi')
	  end if
	  if(right_word(2)-left_word(2)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,311) amsg(1:imsg)
311         format(1x,a,'second bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    bore_list_id(num_bore_list+nbore)='*'
	  else
	    bore_list_id(num_bore_list+nbore)=cline(left_word(2):right_word(2))
	    call casetrans(bore_list_id(num_bore_list+nbore),'hi')
	  end if
	  if(right_word(3)-left_word(3)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,312) amsg(1:imsg)
312         format(1x,a,'third bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    bore_diff_id(num_bore_list)='*'
	  else
	    bore_diff_id(num_bore_list)=cline(left_word(3):right_word(3))
	    call casetrans(bore_diff_id(num_bore_list),'hi')
	  end if
369       continue
	  if(.not.check) go to 367
	  atemp=bore_list_id(num_bore_list)
	  if(atemp.eq.'*') go to 3671
	  do i=1,num_bore_coord
	    if(bore_coord_id(i).eq.atemp) go to 3671
	  end do
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,366) amsg(1:imsg),trim(atemp),trim(bore_coord_file)
366       format(1x,a,'no coordinates for bore ',a, &
	  ' found in bore coordinates file ',a,'.')
	  call write_message(increment=1)
3671      continue
	  atemp=bore_list_id(num_bore_list+nbore)
	  if(atemp.eq.'*') go to 367
	  do i=1,num_bore_coord
	    if(bore_coord_id(i).eq.atemp) go to 367
	  end do
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,366) amsg(1:imsg),trim(atemp),trim(bore_coord_file)
	  call write_message(increment=1)
367       continue
	end do

365     if(imessage.eq.0)then
	  call num2char(num_bore_list,aline)
	  write(amessage,370) trim(aline),trim(bore_list_file)
370       format('  - ',a,' pairs read from bore diff listing file ',a)
	  call write_message
	  ifail=0
          num_bore_list=num_bore_list*2
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
        num_bore_list=num_bore_list*2
9999    close(unit=iunit,iostat=ierr)
	imessage=0
	return

end subroutine read_bore_diff_file
