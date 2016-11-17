subroutine write_mif_file(ifail,outunit,outfile,intarray,east,north,gridspec, &
			izone,type)

! -- Subroutine write_mif_file writes a MAPINFO MIF file containing active
!    finite-difference grid cells as regions.

! -- Arguments are as folows:-
!       ifail:       returned as non-zero if error condition encountered
!       outunit:     unit number of output file
!       outfile:     name of output file
!       intarray:    "window" integer array
!       east, north  x and y coordinates of cell centres of finite-difference
!                    grid expressed in a coordinate system in which the x
!                    direction points east and the origin is at the top left
!                    of the finite difference grid
!       gridspec:    defined variable holding grid specifications
!       izone:       zone number (AMG only)
!       type:        "r" if corresponding MID file contains real array
!                    "i" if corresponding MID file contains integer array

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	implicit none

	integer, intent(out)                    :: ifail
	integer, intent(in)                     :: outunit
	character (len=*)                       :: outfile
	integer, dimension(:,:), intent(in)     :: intarray
	real, dimension(:,:), intent(in)        :: east,north
	type (modelgrid), intent(in)            :: gridspec
	integer, intent(in)                     :: izone
	character (len=1)                       :: type

	integer                                 :: irow,icol,i
	real, dimension(4)                      :: ec,nc
	double precision                        :: etmp,ntmp

	ifail=0
	write(outunit,20,err=900) 
20      format(' Version 2')
	write(outunit,40,err=900)
40      format(' Delimiter ","')
	write(outunit,60) 153-(56-izone)*6
60      format(' Coordsys Earth Projection 8, 12, "m", ',i3, &
	', 0, 0.9996, 500000, 10000000')
	write(outunit,80,err=900)
80      format(' Columns 3')
	write(outunit,90,err=900)
90      format('   Row Integer')
	write(outunit,100,err=900)
100     format('   Column Integer')
	if(type.eq.'r') then
	  write(outunit,110,err=900)
110       format('   Real_val Float')
	else
	  write(outunit,120,err=900)
120       format('   Int_val Integer')
	end if
	write(outunit,130,err=900)
130     format(' Data')
	write(outunit,*)

	do irow=1,gridspec%nrow
	do icol=1,gridspec%ncol
	  if(intarray(icol,irow).eq.0) cycle
	  do i=1,4
	    call corner(i,ec(i),nc(i),icol,irow,gridspec,east,north)
	  end do
	  write(outunit,150,err=900)
150       format(' Region 1')
	  write(outunit,160,err=900)
160       format(' 5')
	  do i=1,2
	    etmp=ec(i)+gridspec%east_corner
	    ntmp=nc(i)+gridspec%north_corner
	    write(outunit,170,err=900) etmp,ntmp
170         format(1x,f15.3,1x,f15.3)
	  end do
	  do i=4,3,-1
	    etmp=ec(i)+gridspec%east_corner
	    ntmp=nc(i)+gridspec%north_corner
	    write(outunit,170,err=900) etmp,ntmp
	  end do
	  etmp=ec(1)+gridspec%east_corner
	  ntmp=nc(1)+gridspec%north_corner
	  write(outunit,170,err=900) etmp,ntmp
	  etmp=east(icol,irow)+gridspec%east_corner
	  ntmp=north(icol,irow)+gridspec%north_corner
	  write(outunit,190,err=900) etmp,ntmp
190       format(' Center',1x,f15.3,1x,f15.3)
	end do
	end do
	return  

900     write(amessage,910) trim(outfile)
910     format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes', endspace='yes')
	ifail=1
	return

end subroutine write_mif_file
 