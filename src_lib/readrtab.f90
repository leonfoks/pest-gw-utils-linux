subroutine read_real_tab_file(ifail,tabunit,tabfile,realarray,gridspec)

! -- Subroutine read_real_tab_file reads a real array table file, checking
!    it for errors and transferring the results to a real array.

! -- Arguments are as follows:-
!       ifail:     returned as non-zero if error condition encountered.
!       tabunit:   unit number of real array table file
!       tabfile:   real array table fileneme
!       realarray: real array
!       gridspec:  defined type containing grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter
	implicit none

	integer, intent(out)                    :: ifail
	integer, intent(in)                     :: tabunit
	character(len=*), intent(in)            :: tabfile
	real, dimension(:,:),intent(inout)      :: realarray
	type(modelgrid), intent(in)             :: gridspec

	integer                                 :: icol,irow,iline
	real                                    :: rval
	character (len=10)                      :: aline

	write(initial_message,20) trim(tabfile)
20      format(' Errors in real array table file ',a,':-')

	write(6,25)
25      format(/,' Working.....')
	iline=0
	imessage=0
30      iline=iline+1
	read(tabunit,'(a)',err=9000,end=1000) cline
	call linesplit(ifail,3)
	if(ifail.lt.0) go to 30
	if(ifail.gt.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  call num2char(iline,aline)
	  write(amessage,50) trim(aline)
50        format('   Line ',a,': insufficient items in line  - 3 expected.')
	  call write_message(increment=1)
	  go to 30
	end if

	irow=char2int(ifail,1)
	if(ifail.ne.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  call num2char(iline,aline)
	  write(amessage,70) trim(aline)
70        format('   Line ',a,': cannot read grid row number from first column.')
	  call write_message(increment=1)
	else
	  if(irow.le.0) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,90) trim(aline)
90          format('   Line ',a,': grid row number less than one.')
	    call write_message(increment=1)
	  else if(irow.gt.gridspec%nrow) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,110) trim(aline)
110         format('   Line ',a,': grid row number exceeds grid row upper bound.')
	    call write_message(increment=1)
	  end if
	end if

	icol=char2int(ifail,2)
	if(ifail.ne.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  call num2char(iline,aline)
	  write(amessage,130) trim(aline)
130       format('   Line ',a,': cannot read grid column number from second ',&
	  'column.')
	  call write_message(increment=1)
	else
	  if(icol.le.0) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,150) trim(aline)
150         format('   Line ',a,': grid column number less than one.')
	    call write_message(increment=1)
	  else if(icol.gt.gridspec%ncol) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    call num2char(iline,aline)
	    write(amessage,160) trim(aline)
160         format('   Line ',a,': grid column number exceeds grid column ',&
	    'upper bound.')
	    call write_message(increment=1)
	  end if
	end if

	rval=char2real(ifail,3)
	if(ifail.ne.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  call num2char(iline,aline)
	  write(amessage,180) trim(aline)
180       format('   Line ',a,': cannot read real array value from ',&
	  'third column.')
	  call write_message(increment=1)
	else
	  if(imessage.eq.0) realarray(icol,irow)=rval
	end if
	go to 30

1000    if(imessage.eq.0) then
	  call num2char(iline-1,aline)
	  write(amessage,1030) trim(aline),trim(tabfile)
1030      format(' - ',a,' lines read from real array table file ',a)
	  call write_message(leadspace='yes',endspace='yes')
	  ifail=0
	  return
	else
	  ifail=1
	  return
	end if

9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(tabfile)
9010    format(' Cannot read line ',a,' of real array table file ',a)
	call write_message(leadspace='yes',endspace='yes')
	ifail=1
	return

end subroutine read_real_tab_file

 