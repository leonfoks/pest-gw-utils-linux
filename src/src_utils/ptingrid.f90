!     Last change:  JD   16 Dec 2000    5:49 pm
program ptingrid

! -- Program PTINGRID locates a point with respect to the finite difference
!    grid and (optionally) provides the real or integer array value pertaining
!    to the cell in which the point lies.


	use defn
	use inter

	implicit none

	integer					:: ifail,i,ierr,outunit,j,irow,&
                                                   icol,icellno,jcellno,idate,iheader
	real					:: fac1,fac2,fac3,fac4
	type (modelgrid) gridspec
	integer, allocatable, dimension(:,:)	:: intarray
	real, allocatable, dimension(:,:)	:: realarray
	character (len=1)			:: aa
	character (len=80)			:: aprompt,outfile
	character (len=10)			:: atemp
	character (len=11)			:: aboreout
	character (len=12)			:: aeast,anorth
	character (len=6)			:: acol,arow


	write(amessage,5)
5       format(' Program PTINGRID locates a point with respect to the finite ',&
	'difference grid and (optionally) provides the real or integer array ',&
	'value pertaining to the cell in which the point lies.')
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

	call readfig(gridspec%specfile,bore_coord_file)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

30	call read_bore_coord_file(ifail, &
       ' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 10
	end if

50	call read_bore_list_file(ifail, &
       ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 30
	end if

60	write(6,*)
70	write(6,80,advance='no')
80	format(' Read integer array, real array or none [i/r/n]: ')
	read(5,'(a)') aa
	if(aa.eq.' ') go to 70
	if(index(eschar,aa).ne.0) then
	  write(6,*)
	  go to 50
	end if
	call casetrans(aa,'lo')
	if((aa.ne.'i').and.(aa.ne.'r').and.(aa.ne.'n')) go to 70
100	if(aa.eq.'i')then
	  allocate(intarray(gridspec%ncol,gridspec%nrow),stat=ierr)
	  if(ierr.ne.0) go to 9000
	  aprompt=' Enter name of integer array file: '
	  call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
          rows=gridspec%nrow,columns=gridspec%ncol)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    deallocate(intarray)
	    go to 60
	  end if
	else if(aa.eq.'r')then
	  allocate(realarray(gridspec%ncol,gridspec%nrow),stat=ierr)
	  if(ierr.ne.0) go to 9000
	  aprompt=' Enter name of real array file: '
	  call read_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
          rows=gridspec%nrow,columns=gridspec%ncol)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    deallocate(realarray)
	    go to 60
	  end if
	end if

110	write(6,*)
	aprompt=' Enter name for output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  if(aa.eq.'i') then
	    deallocate(intarray)
	    write(6,*)
	    escset=0
	    go to 100
	  else if(aa.eq.'r')then
	    deallocate(realarray)
	    write(6,*)
	    escset=0
	    go to 100
	  else
	    go to 60
	  end if
	end if
	if(ifail.ne.0) go to 9900

	write(outunit,150,advance='no',err=9200)
150	format(1x,'point_id',t14,'easting',t29,'northing',t44,'row',t54,&
        'column',t63,' ')
	if(aa.ne.'n')then
	  write(outunit,160,err=9200)
160	  format('array_value')
	else
	  write(outunit,*,err=9200)
	end if

	do i=1,num_bore_list
	  atemp=bore_list_id(i)
	  do j=1,num_bore_coord
	    if(bore_coord_id(j).eq.atemp) go to 190
	  end do
	  write(6,180)
180	  format(/,' Internal error: execution terminated.',/)
	  go to 9900
190	  call factor(gridspec,bore_coord_east(j),bore_coord_north(j), &
          fac1,fac2,fac3,fac4,icellno,jcellno)
	  call cell2rc(icellno,irow,icol,gridspec)
	  aboreout=atemp
	  call char_add(aboreout,',')
	  call num2char(bore_coord_east(j),aeast,11)
	  call char_add(aeast,',')
	  call num2char(bore_coord_north(j),anorth,11)
	  call char_add(anorth,',')
	  call num2char(irow,arow,5)
	  if(icellno.eq.-999) arow='--'
	  call char_add(arow,',')
	  call num2char(icol,acol,5)
	  if(icellno.eq.-999) acol='--'
	  if(aa.ne.'n') call char_add(acol,',')
	  write(outunit,210,err=9200,advance='no') trim(aboreout), &
          trim(aeast),trim(anorth),trim(arow),trim(acol)
210	  format(1x,a,t14,a,t29,a,t44,a,t54,a,t63,' ')
	  if(aa.eq.'i') then
	    if(icellno.eq.-999) then
	      atemp='--'
	    else
	      call num2char(intarray(icol,irow),atemp)
	    end if
	    write(outunit,220,err=9200) trim(atemp) 
	  else if(aa.eq.'r') then
	    if(icellno.eq.-999) then
	      atemp='--'
	    else
	      call num2char(realarray(icol,irow),atemp)
	    end if
	    write(outunit,220,err=9200) trim(atemp)
220	    format(a)
	  else
	    write(outunit,*,err=9200)
	  end if
	end do

	write(amessage,250) trim(outfile)
250	format('  - file ',a,' written ok')
	call write_message

	go to 9900

9000	write(amessage,9010)
9010	format(' Cannot allocate sufficient memory to run program PTINGRID.')
	call write_message
	go to 9900
9200	write(amessage,9210) trim(outfile)
9210	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message
9900    call close_files
	call free_grid_mem(gridspec)
	call free_bore_mem
        deallocate(realarray,intarray,stat=ierr)

end program ptingrid

 