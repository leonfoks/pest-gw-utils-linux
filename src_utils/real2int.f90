!     Last change:  JD   16 Dec 2000    2:39 pm
program real2int

! -- Program REAL2INT builds or modifies a MODFLOW-compatible integer array
!    based on an existing MODFLOW-compatible real array.


	use defn
	use inter

	implicit none

	integer	:: ifail,ierr,ncol,nrow,numreal,irow,icol,i,ikey,itemp,j, &
                   idate, iheader
	real	:: curreal,newreal,sp,rhuge
	integer :: iloc(1)
	integer, allocatable, dimension(:,:)	:: intarray
	real, allocatable, dimension(:,:)	:: realarray
	character (len=80)			:: aprompt,realfile
	type (modelgrid) 			:: gridspec
	character (len=1)			:: ause
	character (len=15)			:: anum1,anum2,anum3

	integer, parameter		:: MAXREAL=100
	integer				:: intval(MAXREAL)
	real				:: realval(MAXREAL),rwork(MAXREAL)
	character (len=15)		:: usernum(MAXREAL)


	write(amessage,5)
5	format(' Program REAL2INT builds or modifies a MODFLOW-compatible ',&
	'integer array based on an existing MODFLOW-compatible real array.')
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
	allocate(realarray(ncol,nrow),intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run REAL2INT.')
	  go to 9890
	end if

55	write(6,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(realarray,intarray,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
!	  call free_grid_mem(gridspec)
	  go to 10
	end if
	realfile=aprompt

80	write(6,*)
90	write(6,110,advance='no')
110	format(' Use ranges or individual values to build integer array? ', &
	' [r/i]: ')
	read(5,'(a)') ause
	if(ause.eq.' ') go to 90
	if(index(eschar,ause).ne.0) go to 55
	call casetrans(ause,'lo')
	if((ause.ne.'i').and.(ause.ne.'r')) go to 90
	if(ause.eq.'r') go to 800

! -- If necessary, the different numbers comprising the array are identified.

	numreal=1
	curreal=realarray(1,1)
	realval(1)=curreal
	row_travel: do irow=1,nrow
	  col_travel: do icol=1,ncol
	    newreal=realarray(icol,irow)
	    sp=spacing(curreal)
	    if((newreal.gt.curreal-sp).and.(newreal.lt.curreal+sp))&
	    cycle col_travel
	    curreal=newreal
	    prev_real: do i=1,numreal
	      sp=spacing(realval(i))
	      if((newreal.gt.realval(i)-sp).and. &
	         (newreal.lt.realval(i)+sp)) cycle col_travel
	    end do prev_real
	    numreal=numreal+1
	    if(numreal.gt.MAXREAL) then
	      call num2char(MAXREAL,anum1)
	      write(amessage,220) trim(anum1),trim(realfile)
220	      format(' To build an integer array on the basis of ',&
	      'individual real array values, the latter array can hold a ', &
	      'maximum of ',a,' different numbers. The array ', &
	      ' contained in file ',a,' holds more than this. To increase ',&
	      'this limit, edit REAL2INT source code, increase parameter ', &
	      'MAXREAL, and recompile.')
	      go to 9890
	    end if
	    realval(numreal)=newreal
	  end do col_travel
	end do row_travel

! -- The identified real numbers are next sorted.

        rwork=0.0
	rhuge=huge(realval(1))
	if(numreal.lt.MAXREAL) then
	  do i=numreal+1,MAXREAL
	    realval(i)=rhuge
	  end do
	end if
	do i=1,numreal
	  iloc=minloc(realval)
	  rwork(i)=realval(iloc(1))
	  realval(iloc(1))=rhuge  
	end do
	realval=rwork

! -- The real numbers are now presented to the user.

300 	write(6,*)
305	write(6,310)
310	format(' The following numbers have been detected in the real ',&
	'array:-')
	write(6,330)
330	format(' Enter corresponding integers.')
	i=0
355	i=i+1
	if(i.gt.numreal) go to 450
	write(anum1,'(1pg11.4)') realval(i)
	anum1=adjustl(anum1)
!	call num2char(realval(i),anum1)
360	write(6,370,advance='no') trim(anum1)
370	format('   enter integer corresponding to real number ',a,': ')
	ikey=key_read(itemp)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  if(i.gt.2) then
	    i=i-2
	    go to 355
	  else if(i.eq.2) then
	    go to 305
	  else
	    go to 90
	  end if
	end if
	if(ikey.eq.-1) then
	  go to 360
	else if(ikey.gt.0) then
	  write(6,380)
380	  format('   Data input error  - try again.')
	  go to 360
	end if
	intval(i)=itemp
	go to 355

450	do i=1,numreal
	  where((realarray.gt.realval(i)-spacing(realval(i))).and. &
	        (realarray.lt.realval(i)+spacing(realval(i)))) &
	        intarray=intval(i)
	end do
	go to 1200

800	write(6,*)
805	write(6,810)
810	format(' Enter the range boundaries:-')
	write(6,820)
820	format(' Use "+i" for "plus infinity" to terminate.')
	i=0
830	i=i+1
	if(i.gt.MAXREAL) then
	  call num2char(MAXREAL,anum1)
	  write(amessage,850) trim(anum1)
850	  format(' A maximum of ',a,' ranges can be defined. To increase ', &
	  'this upper limit edit REAL2INT source code, increase the value ', &
	  'of parameter MAXREAL, and recompile.')
	  go to 9890
	end if
	if(i.eq.1) then
	  anum1='minus infinity'
	else
	  anum1=usernum(i-1)
	end if
	anum1=adjustl(anum1)
	call num2char(i,anum2)
860	write(6,870,advance='no') trim(anum2),trim(anum1)
870	format('   range number ',a,':   ',a,' to ')
	read(5,'(a)') anum3
	if(anum3.eq.' ') go to 860
	anum3=adjustl(anum3)
	if(index(eschar,anum3(1:2)).ne.0) then
	  write(6,*)
	  if(i.eq.1)then
	    go to 90
	  else if(i.eq.2) then
	    go to 805
	  else
	    i=i-2
	    go to 830
	  endif
	end if
	if((anum3.eq.'+i').or.(anum3.eq.'+I').or.(anum3.eq.'i').or. &
	(anum3.eq.'I')) go to 930
	call char2num(ifail,anum3,realval(i))
	if(ifail.gt.0) then
	  write(6,890)
890	  format('   Data input error  - try again.')
	  go to 860
	end if
	if(i.ge.2) then
	  if(realval(i).le.realval(i-1)) then
	    write(6,900)
900	    format('   Range upper bound must exceed range lower bound  - ', &
	    'try again.')
	    go to 860
	  end if
	end if
	usernum(i)=anum3
	go to 830
930	numreal=i
	realval(i)=huge(realval(i))

	write(6,*)
940	write(6,950)
950	format(' Now enter the integers corresponding to these ranges:-')
	j=0
955	j=j+1
	if(j.gt.numreal) go to 1100
	if(j.eq.1) then
	  anum1='minus infinity'
	else
	  anum1=usernum(j-1)
	end if
	if(j.eq.numreal) then
	  anum2='plus infinity'
	else
	  anum2=usernum(j)
	end if
!	anum1=adjustl(anum1)
!	anum2=adjustl(anum2)
960	write(6,970,advance='no') trim(anum1),trim(anum2)
970	format('   enter integer for range ',a,' to ',a,': ')
	ikey=key_read(intval(j))
	if(ikey.lt.0) go to 960
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  if(j.eq.1) then
	    i=numreal-1
	    go to 830
	  else if(j.eq.2) then
	    go to 940
	  else
	    j=j-2
	    go to 955
	  end if
	end if
	if(ikey.gt.0) then
	  write(6,890)
	  go to 960
	end if
	if(intval(j).lt.-999) then
	  write(6,980)
980	  format('   Too small  - try again.')
	  go to 960
	else if(intval(j).gt.9999) then
	  write(6,990)
990	  format('   Too large  - try again.')
	  go to 960
	end if
	go to 955

1100	if(numreal.gt.1) then
	  do i=2,numreal
	    where((realarray.gt.realval(i-1)).and. &
		(realarray.le.realval(i))) intarray=intval(i)
	  end do
	end if
	where (realarray.le.realval(1)) intarray=intval(1)

1200    write(6,*)
        aprompt=' Enter name for output integer array file: '
        call write_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(escset.eq.1) then
          escset=0
          write(6,*)
	  if(ause.eq.'i') then
	    i=numreal-1
	    if(i.eq.0) then
	      go to 305
	    else
	      go to 355
	    end if
	  else
	    j=numreal-1
	    if(j.eq.0) then
	      go to 940
	    else
	      go to 955
	    end if
	  end if
        end if
        write(6,*)
        go to 9900

9890	call write_message(leadspace='yes',endspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,realarray,stat=ierr)

end program real2int

 