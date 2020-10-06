subroutine read_surf_array(ifail,aprompt,array,gridspec)

! -- Subroutine READ_SURF_ARRAY reads a real array from a SURFER grid file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless error condition encountered
!       aprompt:   promp text for SURFER grid filename
!       array:     real array
!       gridspec:  defined type containing grid specifications

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(in)           :: aprompt
	real, intent(out),dimension(:,:)        :: array
	type(modelgrid), intent(in)             :: gridspec
	integer                                 :: iunit,nx,ny,ierr,irow
	real                                    :: slength,delta
	double precision                        :: xlo,xhi,ylo,yhi
	character (len=80)                      :: afile

	ifail=0

5       imessage=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)') afile
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0)then
	  escset=1
	  return
	end if
	if(afile.eq.' ') go to 10

	iunit=nextunit()
	open(unit=iunit,file=afile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.gt.5) go to 9995
	  write(amessage,30) trim(afile)
30        format(' Cannot open SURFER grid file ',a,'  - try again.')
	  call write_message(increment=1)
	  go to 10
	end if
	read(iunit,'(a)',err=9000) cline
	if(cline(1:4).eq.'DSBB')then
	  write(amessage,50) trim(afile)
50        format(' File ',a,' is a binary SURFER grid file: program SRF2REAL ',&
	  'can only read an ASCII grid file. Click on the "Change" button ',&
	  'of the "Output Grid File" section of the SURFER ',&
	  '"Scattered Data Interpolation" dialogue box to make SURFER write ',&
	  'an ASCII grid file.')
	  call write_message(leadspace='yes')
	  go to 9995
	else if (cline(1:4).ne.'DSAA') then
	  write(amessage,70) trim(afile)
70        format(' File ',a,' is not a SURFER grid file  - try again.')
	  call write_message
	  close(unit=iunit)
	  go to 10
	end if
	read(iunit,*,iostat=ierr) nx,ny
	if(ierr.ne.0) then
	  write(amessage,90) trim(afile)
90        format(' Error reading number of grid rows and/or columns from SURFER ',&
	  'grid file ',a)
	  call write_message(leadspace='yes')
	  go to 9995
	end if
	if((nx.ne.gridspec%ncol).or.(ny.ne.gridspec%nrow)) then
	  write(amessage,100) trim(afile),trim(gridspec%specfile)
100       format(' Number of rows and/or columns in SURFER grid file ',a,&
	  ' does not agree with model grid dimensions as read from grid ',&
	  'specification file ',a)
	  call write_message(leadspace='yes')
	  go to 9995
	end if
	read(iunit,*,iostat=ierr) xlo,xhi
	if(ierr.ne.0)then
	  write(amessage,120) trim(afile)
120       format(' Error reading minimum and/or maximum x coordinate from ',&
	  'SURFER grid file ',a)
	  call write_message(leadspace='yes')
	  go to 9995
	end if
	slength=(xhi-xlo)/(gridspec%ncol-1)
	delta=3.0*spacing(slength)
	if((slength.lt.gridspec%delr(1)-delta).or. &
	   (slength.gt.gridspec%delr(1)+delta))then
	   write(amessage,140) trim(afile),trim(gridspec%specfile)
140        format(' Grid x-direction spacing in SURFER grid ',&
	   'file ',a,' differs from row-direction spacing represented in ',&
	   'grid specification file ',a)
	   call write_message(leadspace='yes')
	   go to 9995
	end if
	read(iunit,*,iostat=ierr) ylo,yhi
	if(ierr.ne.0)then
	  write(amessage,160) trim(afile)
160       format(' Error reading minimum and/or maximum y coordinate from ',&
	  'SURFER grid file ',a)
	  call write_message(leadspace='yes')
	  go to 9995
	end if
	slength=(yhi-ylo)/(gridspec%nrow-1)
	delta=3.0*spacing(slength)
	if((slength.lt.gridspec%delc(1)-delta).or. &
	   (slength.gt.gridspec%delc(1)+delta))then
	   write(amessage,170) trim(afile),trim(gridspec%specfile)
170        format(' Grid y-direction spacing in SURFER grid ',&
	   'file ',a,' differs from column-direction spacing represented in ',&
	   'grid specification file ',a)
	   call write_message(leadspace='yes')
	   go to 9995
	end if
	read(iunit,*,iostat=ierr) xlo,xhi
	if(ierr.ne.0)then
	  write(amessage,180) trim(afile)
180       format(' Error reading maximum and/or minimum z coordinate from ',&
	  'SURFER grid file ',a)
	  call write_message(leadspace='yes')
	  go to 9995
	end if

	do irow=gridspec%nrow,1,-1
	  read(iunit,*,iostat=ierr) array(:,irow)
	  if(ierr.ne.0)then
	    write(amessage,200) trim(afile)
200         format(' Error reading array from SURFER grid file ',a)
	    call write_message(leadspace='yes')
	    go to 9995
	  end if
	end do

	write(amessage,230) trim(afile)
230     format('  - array read from SURFER grid file ',a)
	call write_message
	go to 9999

9000    write(amessage,9010) trim(afile)
9010    format(' Error encountered while reading SURFER grid file ',a)
	call write_message(leadspace='yes')
	go to 9995

9995    ifail=1
9999    close(unit=iunit,iostat=ierr)
	return

end subroutine read_surf_array
 