!     Last change:  JD   16 Dec 2000    5:43 pm
program grid2pt

! -- Program GRID2PT tabulates the coordinates of the active cell centres 
!    of the finite difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer					:: ifail,ierr,ncol,nrow, &
						   outunit,layer,icol,irow,&
						   icellno,idate,iheader
	integer, allocatable, dimension(:,:)	:: intarray
	real, allocatable, dimension(:,:)	:: east,north
	character (len=80)			:: aprompt,outfile
	character (len=1)			:: aa
	logical					:: row_column, dummy_layer
	character (len=20)			:: atemp



	write(amessage,5)
5	format(' Program GRID2PT tabulates the coordinates of the active ',&
        'cell centres of the finite difference grid.')
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

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(intarray(ncol,nrow),east(ncol,nrow),north(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run GRID2PT')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

60	continue
	write(6,*)
	aprompt=' Enter name of window integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=gridspec%nrow, columns=gridspec%ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(intarray,east,north,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if
	
80	write(6,*)
	aprompt=' Enter name for output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  go to 60
	end if
	if(ifail.ne.0) go to 9900
	
	write(6,*)
100	write(6,110)
110	format(' In GRID2PT output file:-')
115	write(6,120,advance='no')
120	format(' Use row & column numbers or cell numbers?  [r/c] ')
	read(5,'(a)') aa
	call casetrans(aa,'lo')
	select case (aa)
	  case ('e')
	    close(unit=outunit)
	    go to 80
	  case ('r')
	    row_column= .true.
	  case('c')
	    row_column=.false.
	  case default
	    go to 115
	end select

150	write(6,160,advance='no')
160	format(' Include dummy layer number column? [y/n] ')
	read(5,'(a)') aa
	call casetrans(aa,'lo')
	select case (aa)
	  case('e')
	    write(6,*)
	    go to 100
	  case ('n')
	    dummy_layer=.false.
	  case('y')
	    dummy_layer=.true.
	  case default
	    go to 150
	end select

	if(dummy_layer) then
180	  write(6,190,advance='no')
190	  format(' Enter dummy layer number: ')
	  read(5,'(a)') atemp
	  if(atemp.eq.' ') go to 180
	  if(index(eschar,atemp(1:2)).ne.0) go to 150
	  call char2num(ifail,atemp,layer)
	  if(ifail.ne.0) go to 180
	end if

	write(6,220)
220	format(/,' Working.....')
	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)

	do irow=1,nrow
	  do icol=1,ncol
	    if(intarray(icol,irow).eq.0) cycle
	    if(row_column) then
	      write(outunit,250,advance='no',err=9000) irow,icol,&
	      east(icol,irow)+gridspec%east_corner,&
	      north(icol,irow)+gridspec%north_corner
250	      format(1x,i5,1x,i5,2x,f15.3,2x,f15.3)
	    else
	      call rc2cell(icellno,irow,icol,gridspec)
	      write(outunit,280,advance='no',err=9000)icellno,&
	      east(icol,irow)+gridspec%east_corner,&
	      north(icol,irow)+gridspec%north_corner
280	      format(1x,i11,2x,f15.3,2x,f15.3)
	    end if
	    if(dummy_layer)then
	      write(outunit,240) layer
240	      format(2x,i3)
	    else
	      write(outunit,*)
	    end if
	  end do
	end do

	call num2char(count(intarray.ne.0),atemp)
	write(amessage,260) trim(atemp),trim(outfile)
260	format('   - coordinates for ',a,' active cell centres written to ',&
        'file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9900

9000	write(amessage,9010) trim(outfile)
9010	format('cannot write cell centre information to file ',a,'.')
	call write_message(error='yes',leadspace='yes',endspace='yes')
	go to 9900

9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,east,north,stat=ierr)

end program grid2pt


 