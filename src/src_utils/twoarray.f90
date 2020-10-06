!     Last change:  JD   16 Dec 2000    1:31 am
program twoarray

! -- Program TWOARRAY combines two MODFLOW-compatible real arrays.


	use defn
	use inter

	implicit none

	integer					:: ifail,ierr,ncol,nrow,itemp,iheader,&
                                                   idate
	real					:: thresh1,thresh2
	real, allocatable, dimension(:,:)	:: array1,array2,array3
	character (len=80)			:: aprompt
	type (modelgrid) 			:: gridspec
	character (len=1)			:: effect


	write(amessage,5)
5	format(' Program TWOARRAY combines two MODFLOW-', &
	'compatible real arrays.')
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
	allocate(array1(ncol,nrow),array2(ncol,nrow), &
	array3(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run TWOARRAY.')
	  go to 9890
	end if

55	write(6,*)
	aprompt=' Enter name of first real array file: '
	call read_real_array(ifail,aprompt,array1,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(array1,array2,array3,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
!	  call free_grid_mem(gridspec)
	  go to 10
	end if

60	write(6,65,advance='no')
65	format(' Enter inactive threshold for first array (press ',&
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

80	write(6,*)
90	aprompt=' Enter name of second real array file: '
	call read_real_array(ifail,aprompt,array2,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 60
	end if

110	write(6,120,advance='no')
120	format(' Enter inactive threshold for second array (press ',&
        '<ENTER> if none): ')
	itemp=key_read(thresh2)
	if(escset.eq.1)then
	  escset=0
	  go to 80
	else if(itemp.lt.0) then
	  thresh2=huge(thresh2)
	else if(itemp.gt.0)then
	  write(6,70)
	  go to 110
	endif
	if(pos_test(thresh2,'inactive threshold').ne.0) go to 110

150	write(6,*)
160	write(6,170)
170	format(' Specify effect of second array on first array:-')
180	write(6,200,advance='no')
200	format(' (replace/add/subtract/multiply/divide)  [r/a/s/m/d]: ')
	read(5,'(a)') effect
	if(effect.eq.' ') go to 180
	if(index(eschar,effect).ne.0) then
	  write(6,*)
	  go to 110
	end if
	call casetrans(effect,'lo')
	if((effect.ne.'r').and.(effect.ne.'a').and.(effect.ne.'s').and. &
	   (effect.ne.'m').and.(effect.ne.'d')) go to 180

	array3=array1
	if(effect.eq.'r') then
	  array3=merge(array1,array2, &
	  ((abs(array1).ge.thresh1).or.(abs(array2).ge.thresh2)))
	else if(effect.eq.'a') then
	  where((abs(array1).lt.thresh1).and.(abs(array2).lt.thresh2)) &
	  array3=array1+array2
	else if(effect.eq.'s') then
	  where((abs(array1).lt.thresh1).and.(abs(array2).lt.thresh2)) &
	    array3=array1-array2
	else if(effect.eq.'m')then
	  where((abs(array1).lt.thresh1).and.(abs(array2).lt.thresh2)) &
	    array3=array1*array2
	else if(effect.eq.'d') then
	  itemp=count((abs(array1).lt.thresh1).and.(abs(array2).lt.thresh2).and. &
	  (array2.eq.0.0))
	  if(itemp.ne.0) go to 9100
	  where((abs(array1).lt.thresh1).and.(abs(array2).lt.thresh2)) &
	    array3=array1/array2
	end if

        write(6,*)
        aprompt=' Enter name for output real array file: '
        call write_real_array(ifail,aprompt,array3,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 150
        end if
	write(6,*)
	go to 9900

9100	write(amessage,9110)
9110	format(' Cannot divide second array into first array: at least one ',&
	'zero-valued element in active area.')
	go to 9890
9890	call write_message(leadspace='yes',endspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(array1,array2,array3,stat=ierr)

end program twoarray

 