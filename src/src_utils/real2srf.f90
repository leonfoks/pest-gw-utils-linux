!     Last change:  JD   16 Dec 2000    3:20 pm
program real2srf

! -- Program real2srf translates a real array to surfer GRID format.

	use defn
	use inter

	implicit none


	integer					:: ifail,ierr,ncol,nrow,idate,iheader
	real					:: thresh,temp1,temp2,temp3,temp4
	real, allocatable, dimension(:,:)	:: realarray
	character (len=80)			:: aprompt
	type (modelgrid) 			:: gridspec


	write(amessage,5)
5	format(' Program REAL2SRF translates a real array to surfer GRID format.')
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
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	temp1=minval(gridspec%delr)
	temp2=maxval(gridspec%delr)
	temp3=minval(gridspec%delc)
	temp4=maxval(gridspec%delc)
	if((temp1+2.0*spacing(temp1).lt.temp2).or. &
	   (temp3+2.0*spacing(temp3).lt.temp4)) then
	   write(amessage,20) trim(gridspec%specfile)
20	   format(' The grid defined in grid specification file ',a,&
	   ' is not of uniform dimensions.')
	   call write_message(leadspace='yes')
	   go to 9900
	end if

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(realarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run REAL2SRF.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

55	write(6,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(realarray,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if

60	write(6,65,advance='no')
65	format(' Enter blanking threshold value for this array: ')
	if(key_read(thresh).ne.0) go to 60
	if(escset.eq.1)then
	  escset=0
	  go to 55
	end if
	if(pos_test(thresh,'blanking threshold').ne.0) go to 60
	
80	write(6,*)
	aprompt=' Enter name for SURFER grid file: '
	call write_surf_array(ifail,aprompt,realarray,thresh,gridspec)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 60
	end if

9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(realarray,stat=ierr)

end program real2srf

 