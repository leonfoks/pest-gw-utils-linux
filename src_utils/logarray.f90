!     Last change:  JD   16 Dec 2000    1:31 am
program logarray

! -- Program LOGARRAY takes the log (to base 10) of a real array.

	use defn
	use inter

	implicit none

	integer					:: ifail,ierr,ncol,nrow,itemp,iheader,&
                                                   idate,irow,icol
        real					:: thresh1
        real, allocatable, dimension(:,:)	:: array1
        character (len=80)			:: aprompt
        type (modelgrid) 			:: gridspec
        character (len=1)			:: effect


	write(amessage,5)
5	format(' Program LOGARRAY takes the log (to base 10) of a real array.')
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
	if((ifail.ne.0).or.(escset.eq.1)) then
	  write(6,*)
	  go to 9900
	end if
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) then
	  write(6,*)
	  go to 9900
	end if
	call close_spec_file(gridspec,ok='yes')

        ncol=gridspec%ncol
        nrow=gridspec%nrow
        allocate(array1(ncol,nrow),stat=ierr)
        if(ierr.ne.0) then
          write(amessage,50)
50        format(' Cannot allocate sufficient memory to run LOGARRAY.')
          go to 9890
        end if

55	write(6,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,array1,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(array1,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
!	  call free_grid_mem(gridspec)
	  go to 10
	end if

60	write(6,65,advance='no')
65	format(' Enter inactive threshold for array (press ',&
        '<ENTER> if none): ')
	itemp=key_read(thresh1)
	if(escset.eq.1)then
	  escset=0
	  go to 55
	else if(itemp.lt.0) then
	  thresh1=huge(thresh1)
	else if(itemp.gt.0)then
	  write(6,70)
70	  format(' Data input error  - try again.')
	  go to 60
	endif
	if(pos_test(thresh1,'inactive threshold').ne.0) go to 60

! -- Log transformation is now carried out.

        do irow=1,nrow
          do icol=1,ncol
            if(abs(array1(icol,irow)).ge.thresh1) cycle
            if(array1(icol,irow).le.0.0)then
              write(amessage,72)
72            format(' Log transformation of array cannot take place because at least ',  &
              'one of its below-threshold elements is zero or negative.')
              go to 9890
            end if
            array1(icol,irow)=log10(array1(icol,irow))
          end do
        end do

        write(6,*)
        aprompt=' Enter name for output real array file: '
        call write_real_array(ifail,aprompt,array1,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 55
        end if
	write(6,*)
	go to 9900

9890	call write_message(leadspace='yes',endspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(array1,stat=ierr)

end program logarray

 