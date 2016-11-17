!     Last change:  JD   16 Dec 2000    2:58 pm
program real2mif

! -- Program real2mif writes MAPINFO "MIF" and "MID" files of the finite-
!    difference grid based on the data contained in a real array and on
!    the contents of a window integer array.


	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer                                 :: ifail,ierr,ncol,nrow, &
						   outunit1,outunit2,itemp,&
						   izone,idate,iheader
	integer, allocatable, dimension(:,:)    :: intarray
	real, allocatable, dimension(:,:)       :: east,north,realarray
	character (len=80)                      :: aprompt,outfile1,outfile2


	write(amessage,5)
5       format(' Program REAL2MIF writes MAPINFO "MIF" and "MID" files ', &
	'of the finite difference grid based on the data contained in a ',&
	'real array and on the contents of a window integer array.')
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
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(intarray(ncol,nrow),realarray(ncol,nrow), &
	east(ncol,nrow),north(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to run REAL2MIF.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

52	write(6,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow, columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(intarray,realarray,east,north,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,55)
55	    format(' Memory management error: cannot continue execution.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  call free_grid_mem(gridspec)
	  go to 10
	end if

59	write(6,*)
60	aprompt=' Enter name of window integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow, columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  go to 52
	end if

80      write(6,*)
	aprompt=' Enter name for output "MIF" file: '
	call open_output_file(ifail,aprompt,outfile1,outunit1)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 60
	end if
	if(ifail.ne.0) go to 9900

90      aprompt=' Enter name for output "MID" file: '
	call open_output_file(ifail,aprompt,outfile2,outunit2)
	if(escset.eq.1) then
	  escset=0
	  close(unit=outunit1,iostat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,92)
92	    format(' File handling error: execution terminated.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  go to 80
	end if
	if(ifail.ne.0) go to 9900

	write(6,*)
95	write(6,96,advance='no')
96	format(' Enter AMG zone number of model area  [47 - 58]: ')
	itemp=key_read(izone)
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  close(unit=outunit2,iostat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,92)
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  go to 90
	end if
	if(itemp.lt.0) go to 95
	if(itemp.gt.0) then
	  write(6,97)
97	  format(' Data input error  - try again.')
	  go to 95
	end if
	if((izone.lt.47).or.(izone.gt.58))then
	  write(6,98)
98	  format(' Out of range  - try again.')
	  go to 95
	end if

	write(6,100)
100     format(/,' Working.....',/)
	
	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)
	call write_mif_file(ifail,outunit1,outfile1,intarray,east,north, &
	gridspec,izone,'r')
	if(ifail.ne.0) go to 9900

	write(amessage,150) trim(outfile1)
150     format('  - "MIF" file ',a,' written ok.')
	call write_message()

	call write_real_mid_file(ifail,outunit2,outfile2,intarray,realarray)
	if(ifail.ne.0) go to 9900

	write(amessage,160) trim(outfile2)
160     format('  - "MID" file ',a,' written ok.')
	call write_message(endspace='yes')

9900    call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,realarray,east,north,stat=ierr)

end program real2mif
 