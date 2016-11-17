subroutine read_spec_data(ifail,gridspec)

! -- Subroutine read_spec_data reads grid specifications (except row and
!    column numbers) from a grid specification file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       gridspec:  defined type holding grid specifications

! -- Revision history:-
!       May-June, 1995: version 1.

	use defn
	use inter

	type (modelgrid)        :: gridspec
	integer, intent(out)    :: ifail
	character (len=8)       :: anum
	integer                 :: ierr

	ifail=0
10      gridspec%specline=gridspec%specline+1
	read(gridspec%specunit,'(a)',end=150) cline
	call linesplit(ifail,3)
	if(ifail.lt.0) go to 10
	if(ifail.gt.0) then
	  call num2char(gridspec%specline,anum)
	  write(amessage,20) trim(anum),trim(gridspec%specfile)
20        format('insufficient data items on line ',a,' of file ',a,'.')
	  call write_message(error='yes',leadspace='yes')
	  go to 500
	end if
	gridspec%east_corner=char2double(ifail,1)
	if(ifail.ne.0) go to 150
	gridspec%north_corner=char2double(ifail,2)
	if(ifail.ne.0) go to 150
	gridspec%rotation = char2double(ifail,3)
	if(ifail.ne.0) go to 170
	gridspec%cosang=cos(gridspec%rotation*3.14159265/180.0)
	gridspec%sinang=sin(gridspec%rotation*3.14159265/180.0)

	allocate(gridspec%delr(1:gridspec%ncol),&
	gridspec%delc(1:gridspec%nrow),stat=ierr)
	if(ierr.ne.0)then
	  write(6,50)
50        format(/,' *** Cannot allocate sufficient memory to run program ***',/)
	  ifail=1
	  return
	end if

	read(gridspec%specunit,*,err=400,end=400) gridspec%delr
	read(gridspec%specunit,*,err=400,end=400) gridspec%delc
	if(any(gridspec%delr.le.0.0)) go to 450
	if(any(gridspec%delc.le.0.0)) go to 450
	return

150     call num2char(gridspec%specline,anum)
	write(amessage,160) trim(anum),trim(gridspec%specfile)
160     format('cannot read grid corner coordinates from line ',a, &
	' of file ',a,'.')
	call write_message(error='yes',leadspace='yes')
	go to 500
170     call num2char(gridspec%specline,anum)
	write(amessage,180) trim(anum),trim(gridspec%specfile)
180     format('cannot read grid rotation from line ',a,' of file ',a,'.')
	call write_message(error='yes',leadspace='yes')
	go to 500
400     write(amessage,410) trim(gridspec%specfile)
410     format('cannot read DELR and/or DELC vector from grid specification ',&
	'file ',a,'.')
	call write_message(error='yes',leadspace='yes')
	go to 500
450     write(amessage,460) trim(gridspec%specfile)
460     format('grid cell row or column width less than or equal to zero in ',&
	'file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 500

500     ifail=1
	return

end subroutine read_spec_data
 