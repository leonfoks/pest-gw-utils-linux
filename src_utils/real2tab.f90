!     Last change:  JD   13 Jun 2002    1:55 pm
program real2tab

! -- Program REAL2TAB writes a real array as a three-column real array table
!    file.


	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer                                 :: ifail,ierr,ncol,nrow,icol,irow, &
						   outunit,itemp,idate,iheader
	integer, allocatable, dimension(:,:)    :: intarray
	real, allocatable, dimension(:,:)       :: realarray
	character (len=200)                     :: aprompt,outfile
        character (len=2)                       :: clog


	write(amessage,5)
5       format(' Program REAL2TAB writes a real array as a three-column ',&
	'real array table file.')
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
	allocate(intarray(ncol,nrow),realarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to run REAL2TAB.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call close_spec_file(gridspec,ok='yes')

52	write(6,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow, columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(intarray,realarray,stat=ierr)
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
	aprompt=' Enter name for output real array table file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 60
	end if
	if(ifail.ne.0) go to 9900

90      write(6,100,advance='no')
100     format(' Write native or log10 values to output file?  [n/l]: ')
        read(5,'(a)') clog
        if(clog.eq.' ') go to 90
        clog=adjustl(clog)
        call casetrans(clog,'lo')
        if(clog(1:1).eq.'e') then
          close(unit=outunit,iostat=ierr)
          if(ierr.ne.0)then
	    write(amessage,92)
92	    format(' File handling error: execution terminated.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  go to 80
	end if
        if((clog(1:1).ne.'l').and.(clog(1:1).ne.'n')) go to 90
        if(clog(1:1).eq.'l')then
          do irow=1,nrow
            do icol=1,ncol
              if(intarray(icol,irow).ne.0)then
                if(realarray(icol,irow).le.0.0)then
                  write(amessage,95)
95                format(' At least one element of the real array within the ', &
                  'integer array window is zero or negative - cannot log transform.')
                  call write_message(leadspace='yes',endspace='yes')
                  go to 9900
                end if
              end if
            end do
          end do
          do irow=1,nrow
            do icol=1,ncol
              if(intarray(icol,irow).ne.0) &
              realarray(icol,irow)=log10(realarray(icol,irow))
            end do
          end do
        end if

	call write_real_table_file(ifail,outunit,outfile,intarray,realarray)
	if(ifail.ne.0) go to 9900

	write(amessage,160) trim(outfile)
160     format('  - file ',a,' written ok.')
	call write_message(endspace='yes')

9900    call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,realarray,stat=ierr)

end program real2tab
