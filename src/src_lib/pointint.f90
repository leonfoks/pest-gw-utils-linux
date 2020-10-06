!     Last change:  JD   17 Dec 2000    0:04 am
subroutine point_interp(ncol,nrow,thresh,&
	fac1,fac2,fac3,fac4,icellno,jcellno,bhead,rarray,imethod)

! -- Subroutine INTERP multiplies the heads at grid cell centres by the
!    factors determined in subroutine FACTOR to carry out the interpolation
!    to bore locations. If a bore is near the edge of the grid, or borders
!    one or more dry or inactive cells, dummy heads are assigned to the
!    cell centres or pseudo cell centres surrounding the bore.

! -- Arguments are as follows:-
!       ncol, nrow:     number of columns and rows in finite-difference grid
!       thresh:         abs. value threshold above which cells are inactive
!       fac1,fac2,fac3,fac4:  interpolation factors from grid to point
!       icellno:        index of grid cell containing point
!       jcellno:        index of cell centre to top right of point
!       bhead:          interpolated value at point
!       rarray:         array from which interpolation takes place
!       imethod:        denotes whether to interpolate past last row/column
!                       of active cell centres


	integer, intent(in)                     :: ncol,nrow
	real, intent(in)                        :: thresh
	real, intent(in)                        :: fac1,fac2,fac3,fac4
	integer, intent(in)                     :: icellno,jcellno
	real, intent(out)                       :: bhead
	real, dimension(0:ncol+1,0:nrow+1),intent(in)  :: rarray
	character (len=*), intent(in), optional :: imethod

	integer :: i,j,irow,icol,jrow,jcol,k
	integer :: iact(2,2)
	real :: head(2,2)

	if(icellno.eq.-999) then
	  bhead=7.1e37
	  return
	end if
	irow=(icellno-1)/ncol+1
	icol=icellno-(irow-1)*ncol
	jrow=(jcellno-1)/(ncol+1)
	jcol=jcellno-jrow*(ncol+1)-1
	if(abs(rarray(icol,irow)).gt.thresh) then
	  bhead=5.1e37
	  return
	end if

	head(1,1)=rarray(jcol,jrow)
	head(1,2)=rarray(jcol+1,jrow)
	head(2,1)=rarray(jcol,jrow+1)
	head(2,2)=rarray(jcol+1,jrow+1)
	do k=1,2
	do j=1,2
	iact(j,k)=1
	end do
	end do
	if(abs(head(1,1)).gt.thresh) iact(1,1)=0
	if(abs(head(1,2)).gt.thresh) iact(1,2)=0
	if(abs(head(2,1)).gt.thresh) iact(2,1)=0
	if(abs(head(2,2)).gt.thresh) iact(2,2)=0
	if(imethod.eq.'inside')then
	  do k=1,2
	  do j=1,2
	    if(iact(j,k).ne.1) then
	      bhead=4.1e37
	      return
	    end if
	  end do
	  end do
	end if

	if(iact(1,1).eq.0)then
	  if(iact(1,2).eq.0)then
	    if(iact(2,1).eq.0)then
	      head(1,1)=head(2,2)
	      head(1,2)=head(2,2)
	      head(2,1)=head(2,2)
	    else if(iact(2,2).eq.0)then
	      head(1,1)=head(2,1)
	      head(1,2)=head(2,1)
	      head(2,2)=head(2,1)
	    else
	      head(1,1)=head(2,1)
	      head(1,2)=head(2,2)
	    end if
	  else if(iact(2,1).eq.0)then
	    if(iact(2,2).eq.0)then
	      head(1,1)=head(1,2)
	      head(2,1)=head(1,2)
	      head(2,2)=head(1,2)
	    else
	      head(1,1)=head(1,2)
	      head(2,1)=head(2,2)
	    end if
	  else if(iact(2,2).eq.0)then
	    if(icol.eq.jcol+1)then
	      head(1,1)=head(1,2)
	      head(2,1)=head(1,2)
	      head(2,2)=head(1,2)
	    else
	      head(1,1)=head(2,1)
	      head(1,2)=head(2,1)
	      head(2,2)=head(2,1)
	   end if
	  else
	    if(fac2+fac3.le.1.0e-6)then
	      head(1,1)=(head(1,2)+head(2,1))/2.0
	    else
	      head(1,1)=(head(1,2)*fac2+head(2,1)*fac3)/(fac2+fac3)
	    end if
	  end if
	else if(iact(1,2).eq.0)then
	  if(iact(2,2).eq.0)then
	    if(iact(2,1).eq.0)then
	      head(1,2)=head(1,1)
	      head(2,1)=head(1,1)
	      head(2,2)=head(1,1)
	    else
	      head(1,2)=head(1,1)
	      head(2,2)=head(2,1)
	    end if
	  else if(iact(2,1).eq.0)then
	    if(icol.eq.jcol)then
	      head(1,2)=head(1,1)
	      head(2,1)=head(1,1)
	      head(2,2)=head(1,1)
	    else
	      head(1,1)=head(2,2)
	      head(1,2)=head(2,2)
	      head(2,1)=head(2,2)
	    end if
	  else
	    if(fac1+fac4.le.1.0e-6)then
	      head(1,2)=(head(1,1)+head(2,2))/2.0
	    else
	      head(1,2)=(head(1,1)*fac1+head(2,2)*fac4)/(fac1+fac4)
	    end if
	  end if
	else if(iact(2,1).eq.0)then
	  if(iact(2,2).eq.0)then
	    head(2,1)=head(1,1)
	    head(2,2)=head(1,2)
	  else
	    if(fac1+fac4.le.1.0e-6)then
	      head(2,1)=(head(1,1)+head(2,2))/2.0
	    else
	      head(2,1)=(head(1,1)*fac1+head(2,2)*fac4)/(fac1+fac4)
	    end if
	  end if
	else if(iact(2,2).eq.0)then
	  if(fac2+fac3.le.1.0e-6)then
	    head(2,2)=(head(1,2)+head(2,1))/2.0
	  else
	    head(2,2)=(head(1,2)*fac2+head(2,1)*fac3)/(fac2+fac3)
	  end if
	end if

	bhead=head(1,1)*fac1+head(1,2)*fac2+head(2,1)*fac3+ &
	head(2,2)*fac4

	return
end subroutine point_interp
 