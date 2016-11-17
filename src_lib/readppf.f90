!     Last change:  J    21 Apr 2002    5:18 pm
subroutine read_pilot_points_file(ifail,aprompt,accept_blank)

! -- Subroutine read_pilot_points_file reads a pilot points file, checking
! -- it for errors.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       prompt:    contains text used to prompt for filename
!       accept_blank:  "yes" if a blank response is accepted as indicating no file

	use defn
	use inter

	integer, intent(out)                      :: ifail
	character (len=*), intent(in)             :: aprompt
        character (len=*), intent(in), optional   :: accept_blank
	integer                         :: iunit,ierr,iline,imsg,npoint,i,iend
        character (len=3)               :: ablank
	character (len=20)              :: aline
	character (len=200)             :: amsg,aatemp
	character (len=12)              :: atemp

	if(associated(pilot_point_id)) deallocate (pilot_point_id)
	if(associated(pilot_point_east)) deallocate (pilot_point_east)
	if(associated(pilot_point_north)) deallocate (pilot_point_north)
	if(associated(pilot_point_zone)) deallocate (pilot_point_zone)
	if(associated(pilot_point_val)) deallocate (pilot_point_val)
	nullify(pilot_point_id,pilot_point_east,pilot_point_north,pilot_point_zone, &
        pilot_point_val)

	imessage=0
	ifail=0
        if(present(accept_blank))then
          ablank=accept_blank
        else
          ablank=' '
        end if

10      continue
        if(ablank(1:3).ne.'yes')then
          if(pilot_points_file.eq.' ')then
	    write(6,'(a)',advance='no') aprompt
	  else
	    write(6,'(a)',advance='no') aprompt(1:len(aprompt)-2)//' [' &
	    //trim(pilot_points_file)//']: '
	  end if
        else
          write(6,'(a)',advance='no') aprompt
        end if
	read(5,'(a)') amsg
	amsg=adjustl(amsg)
        if(ablank(1:3).eq.'yes')then
          if(amsg.eq.' ')then
            num_pilot_points=0
            return
          end if
        end if
	if(index(eschar,amsg(1:2)).ne.0) then
	  escset=1
	  return
	end if
	if(amsg.eq.' ') then
	  if(pilot_points_file.eq.' ') go to 10
	  amsg=pilot_points_file
        else
          iend=len_trim(amsg)
          call getfile(ifail,amsg,aatemp,1,iend)
          if(ifail.ne.0) go to 10
          amsg=aatemp
	end if
	iunit=nextunit()
	open(unit=iunit,file=amsg,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9990
	  write(amessage,20) trim(amsg)
20        format('cannot open pilot points file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  if(amsg.eq.pilot_points_file) pilot_points_file=' '
	  go to 10
	end if
	pilot_points_file=amsg
	write(initial_message,30) trim(pilot_points_file)
30      format(' Errors in pilot points file ',a,' ----->')

! -- The pilot points file is first perused to count the number of entries.

	imessage=0
	iline=0
	npoint=0
	do
	  if(imessage.gt.40) go to 9990
	  iline=iline+1
	  read(iunit,'(a)',err=9000,end=200) cline
          if(cline.eq.' ') cycle
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle
	  if(ifail.gt.0) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,70) trim(aline)
70          format('   Line ',a,': insufficient entries (five ', &
	    'entries required).')
	    call write_message(increment=1)
	    cycle
	  end if
	  npoint=npoint+1
	end do

200     if(npoint.eq.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,210)
210       format('   No pilot point data read from file.')
	  call write_message
	  go to 9990
	end if
	if(imessage.ne.0) go to 9990

! -- The pilot points file is next re-read and data placed into memory.

	if(npoint.gt.1000) write(6,220)
220     format(' Working.....')
	rewind(unit=iunit,err=9100)
	allocate(pilot_point_id(npoint),pilot_point_east(npoint), &
	pilot_point_north(npoint),pilot_point_zone(npoint), &
        pilot_point_val(npoint),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,250)
250       format(' Cannot allocate sufficient memory to run program.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9990
	end if

	iline=0
	num_pilot_points=0
	do
	  if(imessage.gt.40) go to 9990
	  iline=iline+1
	  read(iunit,'(a)',end=365,err=9000) cline
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle
	  num_pilot_points=num_pilot_points+1
	  call num2char(iline,aline)
	  write(amsg,280) trim(aline)
280       format('  Line ',a,':')
	  imsg=len_trim(amsg)+1
	  if(right_word(1)-left_word(1)+1.gt.12) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,310) amsg(1:imsg)
310         format(1x,a,'point identifier greater than 12 characters in length.')
	    call write_message(increment=1)
	    pilot_point_id(num_pilot_points)='*'
	  else
	    pilot_point_id(num_pilot_points)=cline(left_word(1):right_word(1))
	    call casetrans(pilot_point_id(num_pilot_points),'hi')
	  end if
	  pilot_point_east(num_pilot_points)=char2double(ifail,2)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,330) amsg(1:imsg)
330         format(1x,a,'cannot read point east coordinate.')
	    call write_message(increment=1)
	  end if
	  pilot_point_north(num_pilot_points)=char2double(ifail,3)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,350) amsg(1:imsg)
350         format(1x,a,'cannot read point north coordinate.')
	    call write_message(increment=1)
	  end if
	  pilot_point_zone(num_pilot_points)=char2int(ifail,4)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) call write_initial_message (leadspace='yes')
	    write(amessage,380) amsg(1:imsg)
380         format(1x,a,'cannot read point zone number.')
	    call write_message(increment=1)
	  else
!	    if(pilot_point_zone(num_pilot_points).eq.0)then
!	      if(imessage.eq.0) call write_initial_message (leadspace='yes')
!	      write(amessage,390) amsg(1:imsg)
!390           format(1x,a,'pilot point zone number cannot be zero.')
!	      call write_message(increment=1)
!	    end if
	  end if
	  pilot_point_val(num_pilot_points)=char2double(ifail,5)
	  if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,395) amsg(1:imsg)
395         format(1x,a,'cannot read value assigned to point.')
	    call write_message(increment=1)
	  end if
	  if(num_pilot_points.ne.1)then
	    atemp=pilot_point_id(num_pilot_points)
	    if(atemp.eq.'*') go to 369
	    duplication: do i=1,num_pilot_points-1
	      if(pilot_point_id(i).eq.atemp) then
		if(imessage.eq.0) call write_initial_message(leadspace='yes')
		write(amessage,360) amsg(1:imsg), trim(atemp)
360             format(1x,a,'point name "',a,'" cited previously.')
		call write_message(increment=1)
		exit duplication
	      end if
	    end do duplication
	  end if
369       continue
	end do

365     if(imessage.eq.0)then
	  call num2char(num_pilot_points,aline)
	  write(amessage,370) trim(aline),trim(pilot_points_file)
370       format('  - data for ',a,' pilot points read from pilot points file ',a)
	  call write_message
	  ifail=0
	  go to 9999
	end if
	go to 9990

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(pilot_points_file)
9010    format('Error reading line ',a,' of pilot points file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9100    write(amessage,9110) trim(pilot_points_file)
9110    format(' Cannot rewind pilot points file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9990    ifail=1
9999    close(unit=iunit,iostat=ierr)
	imessage=0
	return

end subroutine read_pilot_points_file
 