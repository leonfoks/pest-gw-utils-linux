subroutine cell_coordinates(gridspec,east,north,icellrow,icellcol,roff,coff)

! -- Subroutine cell_coordinates calculates the location of a bore within the grid.
!    It also calculates ROFF and COFF as required by the MODFLOW observation
!    package.

! -- Subroutine arguments are as follows:-
!       gridspec:            defined variable holding grid specifications
!       east, north:         east and north coordinates of bore
!       icellrow, icellcol   the row and column within which the bore lies
!       roff, coff           the local cell coordinates of the bore as used by MODFLOW OBS package.

        use defn
        use inter

        type(modelgrid), intent(in)             :: gridspec
        double precision, intent(in)            :: east,north
        real, intent(out)                       :: roff,coff
        integer, intent(out)                    :: icellrow,icellcol

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
          icellrow=-999
          return
        end if

! -- The location of the bore within the finite-difference grid is next
!    determined.

        rtemp1=0.0
        do i=1,ncol
          rtemp2=rtemp1+delr(i)
          if(x.le.rtemp2)then
            icellcol=i
            x1=(rtemp1+rtemp2)*0.5
            coff=(x-x1)/delr(i)
            go to 100
          end if
          rtemp1=rtemp2
        end do
        icellrow=-999
        return

100     rtemp1=0.0
	do i=1,nrow
          rtemp2=rtemp1-delc(i)
          if(y.ge.rtemp2)then
            icellrow=i
            y1=(rtemp1+rtemp2)*0.5
            roff=-(y-y1)/delc(i)
            go to 200
          end if
          rtemp1=rtemp2
        end do
        icellrow=-999
        return

200     continue

        return

end subroutine cell_coordinates
