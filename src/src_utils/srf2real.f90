!     Last change:  JD   16 Dec 2000    5:06 pm
program srf2real

! -- Program srf2real rewrites a SURFER grid file as a MODFLOW-compatible
!    real array.

	use defn
	use inter

	implicit none


	integer                                 :: ifail,ierr,iheader,idate
	real                                    :: temp1,temp2,temp3,&
						   temp4,mblank,delta,sblank,&
						   slow,shigh
	real, allocatable, dimension(:,:)       :: realarray,transarray
	character (len=80)                      :: aprompt
	type (modelgrid)                        :: gridspec


	write(amessage,5)
5       format(' Program SRF2REAL rewrites a SURFER grid file as a ',&
	'MODFLOW-compatible real array.')
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



	sblank=1.701410e38

	call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
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
20         format(' The grid defined in grid specification file ',a,&
	   ' is not of uniform dimensions; hence a SURFER grid file',&
	   ' is not model-compatible.')
	   call write_message(leadspace='yes')
	   go to 9900
	end if

	allocate(realarray(gridspec%ncol,gridspec%nrow),&
	transarray(gridspec%ncol,gridspec%nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to run SRF2REAL.')
	  call write_message(leadspace='yes')
	  go to 9900
	end if

55      write(6,*)
	aprompt=' Enter name of SURFER grid file: '
	call read_surf_array(ifail,aprompt,realarray,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(realarray,transarray,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if

	write(6,*)
60      write(6,65,advance='no')
65      format(' Translate blanked grid values to what number? ')
	if(key_read(mblank).ne.0) go to 60
	if(escset.eq.1)then
	  escset=0
	  go to 55
	end if

	delta=3*spacing(sblank)
	slow=sblank-delta
	shigh=sblank+delta
	transarray=realarray
	where((realarray.ge.slow).and.(realarray.le.shigh)) transarray=mblank

	write(6,*)
80      aprompt=' Enter name for real array file: '
	call write_real_array(ifail,aprompt,transarray,pm_header=headerspec, &
	rows=gridspec%nrow,columns=gridspec%ncol)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 60
	end if

9900    call close_files
	call free_grid_mem(gridspec)
	deallocate(realarray,transarray,stat=ierr)
	write(6,*)

end program srf2real

 