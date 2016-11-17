!     Last change:  JD   16 Dec 2000    2:00 am
program int2real

! -- Program INT2REAL builds or modifies a MODFLOW-compatible real array
!    based on an existing MODFLOW-compatible integer array.


	use defn
	use inter

	implicit none

	integer	:: ifail,ierr,ncol,nrow,itemp,numint,currint,irow,icol,newint, &
		   i,ihuge,corrunit,iline,j,inum,idate,iheader
	real	:: rtemp
	integer :: iloc(1)
	integer, allocatable, dimension(:,:)	:: intarray
	real, allocatable, dimension(:,:)	:: oldarray,newarray
	character (len=80)			:: aprompt,intfile,corrfile
	type (modelgrid) 			:: gridspec
	character (len=1)			:: action,aintreal
	character (len=15)			:: anum1,anum2

        integer, parameter :: MAXINT=50000
	integer		   :: intval(MAXINT),iwork(MAXINT)
	real		   :: realval(MAXINT)	


	write(amessage,5)
5	format(' Program INT2REAL builds or modifies a MODFLOW-compatible ',&
	'real array based on an existing MODFLOW-compatible integer array.')
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
	allocate(newarray(ncol,nrow),intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run INT2REAL.')
	  go to 9890
	end if

55	write(6,*)
	aprompt=' Enter name of integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(newarray,intarray,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
!	  call free_grid_mem(gridspec)
	  go to 10
	end if
	intfile=aprompt

80	write(6,*)
90	write(6,110,advance='no')
110	format(' Modify an existing real array or create a new one?  [m/c]: ')
	read(5,'(a)') action
	if(action.eq.' ') go to 90
	if(index(eschar,action).ne.0) go to 55
	call casetrans(action,'lo')
	if((action.ne.'m').and.(action.ne.'c')) go to 90
	if(action.eq.'m') then
	  allocate(oldarray(ncol,nrow),stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,50)
	    go to 9890
	  end if
	end if
150	if(action.eq.'m') then
	  aprompt=' Enter name of file holding existing real array: '
	  call read_real_array(ifail,aprompt,oldarray,pm_header=headerspec, &
          rows=nrow,columns=ncol)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    deallocate(oldarray,stat=ierr)
	    if(ierr.ne.0) then
	      write(amessage,57)
	      go to 9890
	    end if
	    go to 80
	  end if
	end if

! -- The different integers comprising the integer array are identified.

	numint=1
	currint=intarray(1,1)
	intval(1)=currint
	row_travel: do irow=1,nrow
	  col_travel: do icol=1,ncol
	    newint=intarray(icol,irow)
	    if(newint.eq.currint) cycle col_travel
	    currint=newint
	    prev_int: do i=1,numint
	      if(newint.eq.intval(i)) cycle col_travel
	    end do prev_int
	    numint=numint+1
	    if(numint.gt.MAXINT) then
	      call num2char(MAXINT,anum1)
	      write(amessage,220) trim(anum1),trim(intfile)
220	      format(' An integer array can hold a maximum of ',a, &
	      ' different integers for proper INT2REAL execution. The array ', &
	      ' contained in file ',a,' holds more than this. To increase ',&
	      'this limit, edit INT2REAL source code, increase parameter ', &
	      'MAXINT, and recompile.')
	      go to 9890
	    end if
	    intval(numint)=newint
	  end do col_travel
	end do row_travel

! -- The identified integers are next sorted.

        iwork=0
	ihuge=huge(intval(1))
	if(numint.lt.MAXINT) then
	  do i=numint+1,MAXINT
	    intval(i)=ihuge
	  end do
	end if
	do i=1,numint
	  iloc=minloc(intval)
	  iwork(i)=intval(iloc(1))
	  intval(iloc(1))=ihuge  
	end do
	intval=iwork

	write(6,*)
180	write(6,190,advance='no')
190	format(' Enter integer-real correspondence manually or using a file? ',&
	'[m/f]: ')
	read(5,'(a)') aintreal
	if(aintreal.eq.' ') go to 180
	if(index(eschar,aintreal).ne.0) then
	  if(action.eq.'m') then
	    write(6,*)
	    go to 150
	  else
	    go to 80
	  end if
	end if
	call casetrans(aintreal,'lo')
	if((aintreal.ne.'m').and.(aintreal.ne.'f')) go to 180
	if(aintreal.eq.'f') go to 600

! -- The integers are now presented to the user.

300 	write(6,*)
305	write(6,310)
310	format(' The following integers have been detected in the integer ',&
	'array:-')
	write(6,330)
330	format(' Enter corresponding real numbers.')
	if(action.eq.'m') then
	  write(6,350)
350	  format(' Press <ENTER> if integer effects no change to existing ', &
	  'real array.')
	end if
	i=0
	iwork=-1000
355	i=i+1
	if(i.gt.numint) go to 450
	call num2char(intval(i),anum1)
360	write(6,370,advance='no') trim(anum1)
370	format('   Enter real number corresponding to integer ',a,': ')
	itemp=key_read(rtemp)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  if(i.gt.2) then
	    i=i-2
	    go to 355
	  else if(i.eq.2) then
	    go to 305
	  else
	    go to 180
	  end if
	end if
	if(itemp.eq.-1) then
	  if(action.eq.'c') go to 360
	  go to 355
	else if(itemp.gt.0) then
	  write(6,380)
380	  format('   Data input error  - try again.')
	  go to 360
	end if
	realval(i)=rtemp
	iwork(i)=0
	go to 355
450	go to 1000

! -- An integer-real correspondence file is next read.

600	call open_input_file(ifail, &
	' Enter name of integer-real correspondence file: ',corrfile,corrunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 180
	end if

	iline=0
	imessage=0
	inum=0
	iwork=-1000
	read_a_line: do
	  iline=iline+1
	  read(corrunit,'(a)',err=9000,end=800) cline
	  call linesplit(ifail,2)
	  if(ifail.lt.0) cycle read_a_line
	  inum=inum+1
	  if(ifail.gt.0) then
	    if(imessage.eq.0) then
	      write(amessage,620) trim(corrfile)
620	      format(' Errors in integer-real correspondence file ',a,' ---->')
	      call write_message(leadspace='yes')
	    end if
	    call num2char(iline,anum1)
	    write(amessage,640) trim(anum1)
640	    format(' Line ',a,': insufficient items.')
	    call write_message(increment=1)
	    cycle read_a_line
	  end if
	  itemp=char2int(ifail,1)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) then
	      write(amessage,620) trim(corrfile)
	      call write_message(leadspace='yes')
	    end if
	    call num2char(iline,anum1)
	    write(amessage,650) trim(anum1)
650	    format(' Line ',a,': cannot read first item on line as an integer.')
	    call write_message(increment=1)
	  else
	    do i=1,numint
	      if(itemp.eq.intval(i)) then
	        if(iwork(i).eq.0) then
	          if(imessage.eq.0) then
	            write(amessage,620) trim(corrfile)
	            call write_message(leadspace='yes')
	          end if
	          call num2char(iline,anum1)
	          call num2char(itemp,anum2)
	          write(amessage,660) trim(anum1),trim(anum2)
660	          format(' Line ',a,': integer ',a,' cited previously.')
	          call write_message(increment=1)
	        end if
	        iwork(i)=0
	        go to 700
	      end if
	    end do
	    if(imessage.eq.0) then
	      write(amessage,620) trim(corrfile)
	      call write_message(leadspace='yes')
	    end if
	    call num2char(iline,anum1)
	    write(amessage,670) trim(anum1),trim(intfile)
670	    format(' Line ',a,': integer read as first item does not occur ',&
	    'in integer array read from file ',a)
	    call write_message(increment=1)
	  end if
700	  rtemp=char2real(ifail,2)
	  if(ifail.ne.0) then
	    if(imessage.eq.0) then
	      write(amessage,620) trim(corrfile)
	      call write_message(leadspace='yes')
	    end if
	    call num2char(iline,anum1)
	    write(amessage,720) trim(anum1)
720	    format(' Line ',a,': cannot read second item on line as a real ',&
	    'number.')
	    call write_message(increment=1)
	  else
	    realval(i)=rtemp
	  end if
	end do read_a_line
800	close (unit=corrunit,iostat=ierr)
	if(ierr.ne.0) then
	  write(amessage,810)
810	  format(' File management error: cannot continue execution.')
	  go to 9890
	end if
	if(imessage.ne.0) then
	  write(6,*)
	  go to 9900
	end if
	itemp=count(iwork.eq.0)
	if(itemp.eq.0) then
	  write(amessage,820) trim(corrfile)
820	  format(' No integer-real correspondence pairs read from file ',a)
	  go to 9890
	end if
	if(action.eq.'c') then
	  if(itemp.lt.numint) then
	    write(amessage,840) trim(corrfile),trim(intfile)
840	    format(' Integer-real corresondences were not supplied in file ',a,&
	    ' for all integers found in integer array file ',a)
	    call write_message(leadspace='yes')
	    write(amessage,850)
850	    format('   The following integer(s) were ommitted:-')
	    call write_message
	    cline=' '
	    j=1
	    do i=1,numint
	      if(iwork(i).eq.-1000) then
		call num2char(intval(i),anum1)
		anum1=adjustl(anum1)
	        cline(j+1:j+8)=anum1(1:8)
	        j=j+10
	        if(j.gt.72) then
		  write(6,'(a)') trim(cline)
	          j=1
		  cline=' '
	        end if
	      end if
	    end do
	    if(j.ne.1) write(6,'(a)') trim(cline)
	    write(6,*)
	    go to 9900
	  end if
	end if
	call num2char(inum,anum1)
	write(amessage,870) trim(anum1),trim(corrfile)
870	format(' - ',a,' integer-real pairs read from file ',a)
	call write_message

1000    if(action.eq.'m') newarray=oldarray
	do i=1,numint
	  if(iwork(i).eq.-1000) cycle
	  where(intarray.eq.intval(i)) newarray=realval(i)
	end do

! -- The output real array is next written.

        write(6,*)
        aprompt=' Enter name for output real array file: '
        call write_real_array(ifail,aprompt,newarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(escset.eq.1) then
          escset=0
          write(6,*)
	  if(aintreal.eq.'m') then
	    i=numint-1
	    if(i.eq.0) then
	      go to 305
	    else
	      go to 355
	    end if
	  else
	    go to 600
	  end if
        end if
        write(6,*)
        go to 9900

9000	call num2char(iline,anum1)
	write(amessage,9010) trim(anum1),trim(corrfile)
9010	format(' Error reading line ',a,' of integer-real ',&
	'correspondence file ',a)
	go to 9890
9890	call write_message(leadspace='yes',endspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,newarray,oldarray,stat=ierr)

end program int2real

