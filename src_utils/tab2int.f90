!     Last change:  JD   16 Dec 2000    5:40 pm
program tab2int

! -- Program tab2int reads an integer array written in tabular form and 
!    rewrites it to a normal MODFLOW/MT3D-compatible integer array file.

	use defn
	use inter

	implicit none

	integer	:: ifail,ierr,idefault,tabunit,ncol,nrow,itemp,idate,iheader
	integer, allocatable, dimension(:,:)	:: intarray
	character (len=80)			:: aprompt,tabfile
	type (modelgrid) 			:: gridspec


	write(amessage,5)
5	format(' Program TAB2INT reads an integer array written in ',&
	'tabular form and rewrites it to a normal MODFLOW/MT3D-compatible ',&
	'integer array file.')
	call write_message(leadspace='yes',endspace='yes')

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message
	  go to 9900
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

	call readfig(gridspec%specfile)
10	call spec_open(ifail,gridspec)
	if((ifail.ne.0).or.(escset.eq.1)) then
	  write(6,*)
	  go to 9900
	end if
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) then
	  write(6,*)
	  go to 9900
	end if
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run TAB2INT.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

70	write(6,*)
	aprompt=' Enter name of integer array table file: '
	call open_input_file(ifail,aprompt,tabfile,tabunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(intarray,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,90)
90	    format(' Memory management error: cannot continue execution.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  write(6,*)
	  go to 10
	end if

110	write(6,120,advance='no')
120	format(' Enter value to assign to uncited cells: ')
	itemp=key_read(idefault)
	if(itemp.lt.0) go to 110
	if(escset.ne.0) then
	  escset=0
	  close(unit=tabunit,iostat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,140)
140	    format(' File handling error: cannot continue execution.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  go to 70
	end if
	if(itemp.gt.0)then
	  write(6,145)
145	  format(' Data input error  - try again.')
	  go to 110
	end if
	if((idefault.lt.-999).or.(idefault.gt.9999))then
	  write(6,150)
150	  format(' Number out of range  - try again.')
	  go to 110
	end if
	intarray=idefault

	call read_int_tab_file(ifail,tabunit,tabfile,intarray,gridspec)
	if(ifail.ne.0) then
	  write(6,*)
	  go to 9900
	end if

200     aprompt=' Enter name for output integer array file: '
        call write_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
        if(escset.eq.1) then
          escset=0
	  rewind(unit=tabunit,iostat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,220)
220	    format(' File handling error: cannot continue execution.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  write(6,*)
	  go to 110
	end if
	write(6,*)

9900	call close_files
	deallocate(intarray,stat=ierr)

end program tab2int

 