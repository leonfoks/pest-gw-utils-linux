subroutine read_spec_dim(ifail,gridspec)

! -- Subroutine read_spec_dim reads the grid row and column numbers from the
!    first line of a grid specification file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       gridspec:  defined type variable holding grid specifications

!    Revision history:-
!       June-November, 1995: verson 1.

	use defn
	use inter

	type (modelgrid)                :: gridspec
	integer, intent(out)            :: ifail
	character (len=8)               :: anum

	gridspec%specline=0
	ifail=0
10      gridspec%specline=gridspec%specline+1
	read(gridspec%specunit,'(a)',end=100) cline
	call linesplit(ifail,2)
	if(ifail.lt.0) go to 10
	if(ifail.gt.0) then
	  call num2char(gridspec%specline,anum)
	  write(amessage,20) trim(anum),trim(gridspec%specfile)
20        format('insufficient items on line ',a,' of file ',a,'.')
	  call write_message(error='yes',leadspace='yes')
	  go to 200
	end if
	gridspec%nrow=char2int(ifail,1)
	if(ifail.ne.0) go to 130
	gridspec%ncol=char2int(ifail,2)
	if(ifail.ne.0) go to 130
	return

100     write(amessage,110) trim(gridspec%specfile)
110     format('grid dimensions missing from file ',a,'.')
	call write_message(error='yes',leadspace='yes')
	go to 200
130     call num2char(gridspec%specline,anum)
	write(amessage,140) trim(anum),trim(gridspec%specfile)
140     format('cannot read grid dimensions from line ',a,' of file ',a,'.')
	call write_message(error='yes',leadspace='yes')
	go to 200

200     ifail=1
	return

end subroutine read_spec_dim
 