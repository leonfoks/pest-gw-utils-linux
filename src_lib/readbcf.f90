!     Last change:  JD   23 Dec 2000    8:28 pm
subroutine read_bore_coord_file(ifail,aprompt)

! -- Subroutine read_bore_coord_file reads a bore coordinates file, checking
! -- it for errors.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       prompt:    contains text used to prompt for filename

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: aprompt
	integer                         :: iunit,ierr,iline,imsg,nbore,i,ilay,nbb,ifail1
	character (len=20)              :: aline
	character (len=200)             :: amsg,atempf
	character (len=10)              :: atemp

	if(associated(bore_coord_id)) deallocate (bore_coord_id)
	if(associated(bore_coord_east)) deallocate (bore_coord_east)
	if(associated(bore_coord_north)) deallocate (bore_coord_north)
	if(associated(bore_coord_layer)) deallocate (bore_coord_layer)
	nullify(bore_coord_id,bore_coord_east,bore_coord_north,bore_coord_layer)

	imessage=0
	ifail=0
10      if(bore_coord_file.eq.' ')then
	  write(6,'(a)',advance='no') aprompt
	else
	  write(6,'(a)',advance='no') aprompt(1:len(aprompt)-2)//' [' &
	  //trim(bore_coord_file)//']: '
	end if
	read(5,'(a)') amsg
	amsg=adjustl(amsg)
	if(index(eschar,amsg(1:2)).ne.0) then
	  escset=1
	  return
	end if
	if(amsg.eq.' ') then
	  if(bore_coord_file.eq.' ') go to 10
	  amsg=bore_coord_file
	else
          nbb=len_trim(amsg)
          call getfile(ifail1,amsg,atempf,1,nbb)
          if(ifail1.ne.0) go to 10
          amsg=atempf
        end if
	iunit=nextunit()
	open(unit=iunit,file=amsg,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9990
	  write(amessage,20) trim(amsg)
20        format('cannot open bore coordinates file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  if(amsg.eq.bore_coord_file) bore_coord_file=' '
	  go to 10
	end if
	bore_coord_file=amsg
	write(initial_message,30) trim(bore_coord_file)
30      format(' Errors in bore coordinates file ',a,' ----->')

! -- The bore coordinates file is first perused to count the number of entries.

	imessage=0
	iline=0
	nbore=0
	do
	  if(imessage.gt.40) go to 9990
	  iline=iline+1
	  read(iunit,'(a)',err=9000,end=200) cline
	  call linesplit(ifail,3)
	  if(ifail.lt.0) cycle
	  if(ifail.gt.0) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,70) trim(aline)
70          format('   Line ',a,': insufficient entries (minimum of three ', &
	    'entries required).')
	    call write_message(increment=1)
	    cycle
	  end if
	  nbore=nbore+1
	end do

200     if(nbore.eq.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,210)
210       format('   No bore coordinates read from file.')
	  call write_message
	  go to 9990
	end if
	if(imessage.ne.0) go to 9990

! -- The bore coordinates file is next re-read and data placed into memory.

	if(nbore.gt.1000) write(6,220)
220     format(' Working.....')
	rewind(unit=iunit,err=9100)
	allocate(bore_coord_id(nbore),bore_coord_east(nbore), &
	bore_coord_north(nbore),bore_coord_layer(nbore),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,250)
250       format(' Cannot allocate sufficient memory to run program.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9990
	end if

	iline=0
	num_bore_coord=0
	do
	  if(imessage.gt.40) go to 9990
	  iline=iline+1
	  read(iunit,'(a)',end=365,err=9000) cline
	  call linesplit(ifail,4)
	  if(ifail.lt.0) cycle
	  if(ifail.gt.0) then
	    ilay=0
	    call linesplit(ifail,3)
	  else
	    ilay=1
	  end if
	  num_bore_coord=num_bore_coord+1
	  call num2char(iline,aline)
	  write(amsg,280) trim(aline)
280       format('  Line ',a,':')
	  imsg=len_trim(amsg)+1
	  if(right_word(1)-left_word(1)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,310) amsg(1:imsg)
310         format(1x,a,'bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    bore_coord_id(num_bore_coord)='*'
	  else
	    bore_coord_id(num_bore_coord)=cline(left_word(1):right_word(1))
	    call casetrans(bore_coord_id(num_bore_coord),'hi')
	  end if
	  bore_coord_east(num_bore_coord)=char2double(ifail,2)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,330) amsg(1:imsg)
330         format(1x,a,'cannot read bore east coordinate.')
	    call write_message(increment=1)
	  end if
	  bore_coord_north(num_bore_coord)=char2double(ifail,3)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,350) amsg(1:imsg)
350         format(1x,a,'cannot read bore north coordinate.')
	    call write_message(increment=1)
	  end if
	  if(ilay.eq.0) then
	    bore_coord_layer(num_bore_coord)=-999
	  else
	    bore_coord_layer(num_bore_coord)=char2int(ifail,4)
	    if(ifail.ne.0) then
	      if(imessage.eq.0) call write_initial_message (leadspace='yes')
	      write(amessage,380) amsg(1:imsg)
380           format(1x,a,'cannot read bore layer number.')
	      call write_message(increment=1)
	    else
	      if(bore_coord_layer(num_bore_coord).le.0)then
		if(imessage.eq.0) call write_initial_message (leadspace='yes')
		write(amessage,390) amsg(1:imsg)
390             format(1x,a,'illegal bore layer number.')
		call write_message(increment=1)
	      end if
	    end if
	  end if
	  if(num_bore_coord.ne.1)then
	    atemp=bore_coord_id(num_bore_coord)
	    if(atemp.eq.'*') go to 369
	    duplication: do i=1,num_bore_coord-1
	      if(bore_coord_id(i).eq.atemp) then
		if(imessage.eq.0) call write_initial_message(leadspace='yes')
		write(amessage,360) amsg(1:imsg), trim(atemp)
360             format(1x,a,'bore ',a,' cited previously.')
		call write_message(increment=1)
		exit duplication
	      end if
	    end do duplication
	  end if
369       continue
	end do

365     if(imessage.eq.0)then
	  call num2char(num_bore_coord,aline)
	  write(amessage,370) trim(aline),trim(bore_coord_file)
370       format('  - ',a,' bores and coordinates read from bore ',&
	  'coordinates file ',a)
	  call write_message
	  ifail=0
	  go to 9999
	end if
	go to 9990

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(bore_coord_file)
9010    format('Error reading line ',a,' of bore coordinates file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9100    write(amessage,9110) trim(bore_coord_file)
9110    format(' Cannot rewind bore coordinates file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9990    ifail=1
9999    close(unit=iunit,iostat=ierr)
	imessage=0
	return

end subroutine read_bore_coord_file
 