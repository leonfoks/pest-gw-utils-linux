!     Last change:  JD   16 Dec 2000    4:36 pm
program grid2arc

! -- Program GRID2ARC produces a "line" and "points" file of the 
!    finite-difference grid in ARCINFO "generate" format.

	use defn
	use inter

	implicit none

interface

	subroutine linefile(iunit,iarray,east,north,gridspec)
	  use defn
	  integer, intent(in)                     :: iunit
	  integer, dimension(:,:), intent(in)     :: iarray
	  real, dimension(:,:), intent(in)        :: east,north
	  type(modelgrid), intent(in)             :: gridspec
	end subroutine linefile

	subroutine pointfile(iunit,iarray,east,north,gridspec)
	  use defn
	  integer, intent(in)                     :: iunit
	  integer, dimension(:,:), intent(in)     :: iarray
	  real, dimension(:,:), intent(in)        :: east,north
	  type(modelgrid), intent(in)             :: gridspec
	end subroutine pointfile

end interface


	type (modelgrid) gridspec

	integer                                 :: ifail,ierr,ncol,nrow, &
						   outunit1,outunit2,idate,iheader
	integer, allocatable, dimension(:,:)    :: intarray
	real, allocatable, dimension(:,:)       :: east,north
	character (len=80)                      :: aprompt,outfile1,outfile2
	character (len=1)                       :: aa
	character (len=20)                      :: atemp


	write(amessage,5)
5       format(' Program GRID2ARC produces a "line" and "points" file of the ',&
	'finite-difference grid in ARCINFO "generate" format.')
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
	allocate(intarray(ncol,nrow),east(ncol,nrow),north(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to run GRID2ARC')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

60      continue
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
	

80      write(6,*)
	aprompt=' Enter name for output "line" file: '
	call open_output_file(ifail,aprompt,outfile1,outunit1)
	if(escset.eq.1) then
	  escset=0
	  go to 60
	end if
	if(ifail.ne.0) go to 9900

90      aprompt=' Enter name for output "points" file: '
	call open_output_file(ifail,aprompt,outfile2,outunit2)
	if(escset.eq.1) then
	  escset=0
	  go to 80
	end if
	if(ifail.ne.0) go to 9900

	write(6,100)
100     format(/,' Working.....',/)
	
	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)

	call linefile(outunit1,intarray,east,north,gridspec)

	write(amessage,150) trim(outfile1)
150     format('  - ARCINFO generate "line" file ',a,' written ok.')
	call write_message()

	call pointfile(outunit2,intarray,east,north,gridspec)

	write(amessage,160) trim(outfile2)
160     format('  - ARCINFO generate "points" file ',a,' written ok.')
	call write_message(endspace='yes')

9900    call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,east,north,stat=ierr)

end program grid2arc




subroutine linefile(iunit,iarray,east,north,gridspec)

! -- Subroutine linefile writes an ARCINFO generate "line" file.

! -- Arguments are as follows:-
!       iunit:       file unit number for output
!       iarray:      "window" integer array
!       east, north: x and y coordinates of grid cell centres expressed in a
!                    coordinate system where the x direction is oriented
!                    east and the origin is at the top left corner of the grid
!       gridspec:    defined variable holding grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	implicit none

	integer, intent(in)                     :: iunit
	integer, dimension(:,:), intent(in)     :: iarray
	real, dimension(:,:), intent(in)        :: east,north
	type(modelgrid), intent(in)             :: gridspec
	integer                                 :: icol,irow,i,icellno
	real, dimension(4)                      :: ec,nc

	do irow=1,gridspec%nrow
	do icol=1,gridspec%ncol
	  if(iarray(icol,irow).eq.0) cycle
	  do i=1,4
	    call corner(i,ec(i),nc(i),icol,irow,gridspec,east,north)
	  end do
	  call rc2cell(icellno,irow,icol,gridspec)
	  write(iunit,50,err=100) icellno
50        format(1x,i7)
	  do i=1,2
	    write(iunit,70,err=100) gridspec%east_corner+ec(i), &
				    gridspec%north_corner+nc(i)
70          format(1x,f15.3,',  ',f15.3)
	  end do
	  do i=4,3,-1
	    write(iunit,70,err=100) gridspec%east_corner+ec(i), &
				    gridspec%north_corner+nc(i)
	  end do
	  write(iunit,70,err=100) gridspec%east_corner+ec(1), &
				  gridspec%north_corner+nc(1)
	  write(iunit,90)
90        format(1x,'end')
	end do
	end do
	write(iunit,90)
	return

100     write(6,110)
110     format(' Unable to write "line" file: file inaccessible or disk full.',/)
	stop

end subroutine linefile


subroutine pointfile(iunit,iarray,east,north,gridspec)

! -- Subroutine pointfile writes an ARCINFO generate "points" file.

! -- Arguments are as follows:-
!       iunit:       file unit number for output
!       iarray:      "window" integer array
!       east, north: x and y coordinates of grid cell centres expressed in a
!                    coordinate system where the x direction is oriented
!                    east and the origin is at the top left corner of the grid
!       gridspec:    defined variable holding grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	implicit none

	integer, intent(in)                     :: iunit
	integer, dimension(:,:), intent(in)     :: iarray
	real, dimension(:,:), intent(in)        :: east,north
	type(modelgrid), intent(in)             :: gridspec
	integer                                 :: icol,irow,icellno

	do irow=1,gridspec%nrow
	do icol=1,gridspec%ncol
	if(iarray(icol,irow).eq.0) cycle
	call rc2cell(icellno,irow,icol,gridspec)
	write(iunit,50,err=100) icellno,gridspec%east_corner+east(icol,irow), &
					gridspec%north_corner+north(icol,irow)
50      format(1x,i6,',  ',f15.3,',  ',f15.3)
	end do
	end do
	write(iunit,60)
60      format(1x,'end')
	return

100     write(6,110)
110     format(' Unable to write "points" file: file inaccessible or disk full.',/)
	stop
	
end subroutine pointfile
 