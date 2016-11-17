program genreal2srf

! -- Program genreal2srf interpolates from a MODFLOW array to a surfer grid file.

	use defn
	use inter

	implicit none

        integer  :: ifail,idate,iheader,ierr,ibeg,iend,icount,iunit
        integer  :: ncol,nrow,nx,ny,ix,iy,icellno,jcellno,irow,icol
        real     :: blankval,thresh,gt_thresh,rtemp
        real     :: fac1,fac2,fac3,fac4,zmin,zmax
        double precision :: xmin,ymin,deltax,deltay
        character*10  :: atemp
        character*30  :: anum
        character*200 :: aprompt,afile,outfile
        type (modelgrid) :: gridspec

        integer, allocatable :: intarray(:,:)
        real, allocatable    :: realarray(:,:),rarray(:,:),rval(:,:)
        double precision, allocatable :: xx(:),yy(:)


        write(amessage,5)
5	format(' Program GENREAL2SRF interpolates from a MODFLOW array to the nodes of ',  &
        'a SURFER grid, irrespective of the design of the MODFLOW grid. It then writes the SURFER grid file.')
	call write_message(leadspace='yes',endspace='yes')

        blankval=1.70141e38
        thresh=1.0e35

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

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(realarray(ncol,nrow),rarray(0:ncol+1,0:nrow+1),intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run GENREAL2SRF.')
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
	  deallocate(realarray,intarray,rarray,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if
        gt_thresh=thresh+2*spacing(thresh)
        rarray(0,:)=gt_thresh
        rarray(ncol+1,:)=gt_thresh
        rarray(:,0)=gt_thresh
        rarray(:,nrow+1)=gt_thresh
        do irow=1,nrow
          do icol=1,ncol
            rarray(icol,irow)=realarray(icol,irow)
          end do
        end do

60	continue
	write(6,*)
	aprompt=' Enter name of window integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=gridspec%nrow,columns=gridspec%ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  go to 55
	end if

        write(6,*)
        write(6,70)
70      format(' Enter specifications for SURFER grid: ')
75      write(6,80,advance='no')
80      format('    X direction grid minimum: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 75
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          go to 60
        end if
        call char2num(ifail,anum,xmin)
        if(ifail.ne.0) go to 75
85      write(6,90,advance='no')
90      format('    X direction spacing: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 85
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          write(6,*)
          go to 75
        end if
        call char2num(ifail,anum,deltax)
        if(ifail.ne.0) go to 85
95      write(6,100,advance='no')
100     format('    No. of X direction nodes: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 95
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          write(6,*)
          go to 85
        end if
        call char2num(ifail,anum,nx)
        if(ifail.ne.0) go to 95

        write(6,*)
105     write(6,110,advance='no')
110     format('    Y direction grid minimum: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 105
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          write(6,*)
          go to 95
        end if
        call char2num(ifail,anum,ymin)
        if(ifail.ne.0) go to 105
115     write(6,120,advance='no')
120     format('    Y direction spacing: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 115
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          write(6,*)
          go to 105
        end if
        call char2num(ifail,anum,deltay)
        if(ifail.ne.0) go to 115
125     write(6,130,advance='no')
130     format('    No. of Y direction nodes: ')
        read(5,'(a)') anum
        if(anum.eq.'') go to 125
        anum=adjustl(anum)
        call casetrans(anum,'lo')
        if(anum(1:2).eq.'e ') then
          write(6,*)
          go to 115
        end if
        call char2num(ifail,anum,ny)
        if(ifail.ne.0) go to 125

        write(6,*)
140     write(6,150,advance='no')
150     format(' Enter name for SURFER grid file: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 140
        afile=adjustl(afile)
        atemp=afile(1:2)
        call casetrans(atemp,'lo')
        if(atemp(1:2).eq.'e ')then
          write(6,*)
          go to 125
        end if
        ibeg=1
        iend=len_trim(afile)
        call getfile(ifail,afile,outfile,ibeg,iend)
        if(ifail.ne.0) go to 140

! -- The x and y coordinates of each point in the SURFER grid are now evaluated.

        allocate(xx(nx),yy(ny),rval(nx,ny),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,50)
	  call write_message
          go to 9900
        end if
        do ix=1,nx
          xx(ix)=(ix-1)*deltax+xmin
        end do
        do iy=1,ny
          yy(iy)=(iy-1)*deltay+ymin
        end do

! -- The real array is "blanked" using the integer array.

        do irow=1,nrow
          do icol=1,ncol
            if(intarray(icol,irow).eq.0) rarray(icol,irow)=gt_thresh
          end do
        end do

! -- Interpolation is now carried out and an array built.

        icount=0
        do ix=1,nx
          do iy=1,ny
            call factor(gridspec,xx(ix),yy(iy),fac1,fac2,fac3,fac4,icellno,jcellno)
            if(icellno.eq.-999)then
              rval(ix,iy)=blankval
            else
              icount=icount+1
              call point_interp(ncol,nrow,thresh,fac1,fac2,fac3,fac4,icellno,jcellno,rval(ix,iy),rarray)
            end if
          end do
        end do
	if(icount.eq.0) then
	  write(amessage,160)
160	  format(' Spatial interpolation could not take place to any cells in the SURFER grid from ', &
          'the MODFLOW real array.')
	  call write_message(leadspace='yes')
	  go to 9900
	end if

! -- The SURFER grid file is written.

	where (abs(rval).ge.thresh) rval=1.70141e38
	iunit=nextunit()
	open(unit=iunit,file=outfile)
	write(iunit,'(a)')'DSAA'
	write(iunit,'(1x,i5,1x,i5)') nx,ny
	write(iunit,240) xx(1),xx(nx)
240     format(1x,1pe18.10,2x,1pe18.10)
	write(iunit,240) yy(1),yy(ny)
	zmax=-1.1e35
	zmin=1.1e35
	do iy=1,ny
	  do ix=1,nx
	    rtemp=rval(ix,iy)
	    if(abs(rtemp).lt.thresh) then
	      if(rtemp.gt.zmax)zmax=rtemp
	      if(rtemp.lt.zmin)zmin=rtemp
	    end if
	  end do
	end do
	write(iunit,240) zmin,zmax
	do iy=1,ny
	  write(iunit,265) (rval(ix,iy),ix=1,nx)
265       format(7(1pe14.6))
	  write(iunit,*)
	end do
	write(6,270) trim(outfile)
270     format('  - SURFER grid file ',a,' written ok.')

9900	call close_files
        write(6,*)
	call free_grid_mem(gridspec)
	deallocate(intarray,realarray,rarray,rval,xx,yy,stat=ierr)

end program genreal2srf


