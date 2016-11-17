subroutine factor(gridspec,east,north,fac1,fac2,fac3,fac4,&
	icellno,jcellno)

! -- Subroutine factor calculates interpolation factors for a bore.
!    There are four factors, one for each of the four cell centres
!    surrounding the bore. When interpolation takes place, the heads at the
!    respective cell centres are multiplied by these factors and added
!    together to determine the head at the bore.

! -- Subroutine arguments are as follows:-
!       gridspec:            defined variable holding grid specifications
!       east, north:         east and north coordinates of bore
!       fac1,fac2,fac3,fac4  interpolation factors for bore
!       icellno              identifier for cell in which bore is situated
!       jcellno              identifier for the cell centre to the top left
!                            of the bore

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	type(modelgrid), intent(in)             :: gridspec
	double precision, intent(in)            :: east,north
	real, intent(out)                       :: fac1,fac2,fac3,fac4
	integer, intent(out)                    :: icellno,jcellno

	integer :: i,irow,icol,jrow,jcol,nrow,ncol
	real :: x,y,etemp,ntemp,rtemp1,rtemp2, &
	x1,y1,x2,y2,delx,dely
	real, dimension(:),pointer :: delr,delc

	nrow=gridspec%nrow
	ncol=gridspec%ncol
	delr=>gridspec%delr
	delc=>gridspec%delc

! -- First the bore coordinates are expressed as local grid coordinates.

	etemp=east-gridspec%east_corner
	ntemp=north-gridspec%north_corner
	x=etemp*gridspec%cosang+ntemp*gridspec%sinang
	y=ntemp*gridspec%cosang-etemp*gridspec%sinang
	if((x.lt.0.0).or.(y.gt.0.0))then
	  icellno=-999
	  return
	end if

! -- The location of the bore within the finite-difference grid is next
!    determined.

	rtemp1=0.0
	rtemp2=delr(1)/2.0
	do i=1,ncol+1
	if(x.le.rtemp2)then
	  x1=x-rtemp1
	  x2=rtemp2-x
	  delx=rtemp2-rtemp1
	  jcol=i-1
	  if(i.eq.1)then
	    icol=1
	  else if(i.eq.ncol+1)then
	    icol=ncol
	  else
	    if(x1.le.delr(i-1)/2.0)then
	      icol=i-1
	    else
	      icol=i
	    end if
	  end if
	  go to 100
	end if
	rtemp1=rtemp2
	if(i.lt.ncol)then
	  rtemp2=rtemp2+(delr(i)+delr(i+1))/2.0
	else if(i.eq.ncol) then
	  rtemp2=rtemp2+delr(i)/2.0
	end if
	end do
	icellno=-999
	return

100     rtemp1=0.0
	rtemp2=-delc(1)/2.0
	do i=1,nrow+1
	if(y.ge.rtemp2)then
	  y1=rtemp1-y
	  y2=y-rtemp2
	  dely=rtemp1-rtemp2
	  jrow=i-1
	  if(i.eq.1)then
	    irow=1
	  else if(i.eq.nrow+1)then
	    irow=nrow
	  else if(y1.le.delc(i-1)/2.0)then
	    irow=i-1
	  else
	    irow=i
	  end if
	  go to 200
	end if
	rtemp1=rtemp2
	if(i.lt.nrow)then
	  rtemp2=rtemp2-(delc(i)+delc(i+1))/2.0
	else if(i.eq.nrow) then
	  rtemp2=rtemp2-delc(i)/2.0
	end if
	end do
	icellno=-999
	return

! -- The interpolation factors are calculated.

200     jcellno=(ncol+1)*jrow+jcol+1
	icellno=ncol*(irow-1)+icol
	rtemp1=1.0/delx/dely
	fac1=x2*y2*rtemp1
	fac2=x1*y2*rtemp1
	fac3=x2*y1*rtemp1
	fac4=x1*y1*rtemp1

	return

end subroutine factor
 