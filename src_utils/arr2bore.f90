program arr2bore

! -- Program ARR2BORE undertakes spatial interpolation from a MODFLOW/MT3D-compatible real array
!    to a set of points.

	use defn
	use inter

	implicit none

        integer            :: ifail,iheader,i,j,ierr,itemp,outunit,idate
        integer            :: ncol,nrow,icol,irow
        real               :: thresh,gt_thresh
        character*20       :: aval
        character*100      :: aprompt
        character*200      :: outfile

        integer, allocatable, dimension(:)         :: icellno,jcellno
        real, allocatable, dimension(:)            :: fac1,fac2,fac3,fac4,bore_val
        real, allocatable, dimension(:,:)          :: rarray,rarray1
        double precision, allocatable, dimension(:):: east,north
	type (modelgrid) gridspec


	write(amessage,5)
5       format(' Program ARR2BORE undertakes spatial interpolation from a MODFLOW/MT3D-compatible ', &
        'real array to a set of points.')
	call write_message(leadspace='yes',endspace='yes')

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          call write_message
          go to 9900
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          call write_message
          go to 9900
        endif
        if((iheader.ne.0).or.(headerspec.eq.' ')) then
          write(amessage,6)
6         format(' Cannot read array header specification from settings file ', &
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

	ncol=gridspec%ncol
	nrow=gridspec%nrow

        write(6,*)
30      call read_bore_coord_file(ifail, &
	' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 10
	end if

100     call read_bore_list_file(ifail, &
       ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 30
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
	        write(amessage,110) trim(bore_list_id(i)),trim(bore_list_file)
110		format(' Execution of program ARR2OBS cannot proceed as ',&
		'there are multiple occurrences of the same bore ',&
		'(i.e. bore ',a,') in bore listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	allocate(east(num_bore_list), north(num_bore_list), &
	fac1(num_bore_list), fac2(num_bore_list), fac3(num_bore_list), &
	fac4(num_bore_list), icellno(num_bore_list), jcellno(num_bore_list),  &
        bore_val(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,130)
130	  format(' Insufficient memory available to continue execution.')
	  go to 9890
	end if

	imessage=0
	list: do i=1,num_bore_list
	  do j=1,num_bore_coord
	    if(bore_coord_id(j).eq.bore_list_id(i)) then
	      east(i)=bore_coord_east(j)
	      north(i)=bore_coord_north(j)
	      cycle list
	    end if
	  end do
	  write(amessage,140) trim(bore_list_id(i)),trim(bore_list_file),&
	  trim(bore_coord_file)
140	  format(' No coordinates for bore ',a,' from bore listing ',&
	  'file ',a,' are provided in bore coordinates file ',a)
	  if(imessage.eq.0) write(6,*)
	  call write_message(increment=1)
	end do list
	if(imessage.ne.0) go to 9900

	do i=1,num_bore_list
	  call factor(gridspec,east(i),north(i),fac1(i),fac2(i),fac3(i), &
	  fac4(i),icellno(i),jcellno(i))
	end do
	if(all(icellno.eq.-999)) then
	  write(amessage,145) trim(bore_list_file),trim(gridspec%specfile)
145	  format(' None of the bores cited in bore listing file ',a, &
	  ' are within the bounds of the finite difference grid as defined ',&
	  'in grid specification file ',a)
	  go to 9890
	end if

        allocate(rarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200


155     write(6,*)
156     aprompt=' Enter name of real array file: '
        call read_real_array(ifail,aprompt,rarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) then
          escset=0
          write(6,*)
          deallocate(east,north,fac1,fac2,fac3,fac4,icellno,jcellno,bore_val,rarray,stat=ierr)
          if(ierr.ne.0) then
            write(amessage,157)
157         format(' Memory management error: cannot continue execution.')
            go to 9890
          end if
	  go to 100
        end if

160     write(6,165,advance='no')
165     format(' Enter inactive threshold for this array (<Enter> if 1E35): ')
        itemp=key_read(thresh)
        if(escset.eq.1)then
          escset=0
          go to 155
        else if(itemp.lt.0) then
          thresh=1.0E35
        else if(itemp.gt.0)then
	  write(6,170)
170	  format(' Data input error  - try again.')
	  go to 160
        endif
        if(pos_test(thresh,'inactive threshold').ne.0) go to 160


400     write(6,*)
410     call open_output_file(ifail, &
        ' Enter name for bore information output file: ',outfile,outunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 160
        end if

	allocate(rarray1(0:ncol+1,0:nrow+1),stat=ierr)
	if(ierr.ne.0) go to 9200
        rarray1=0.0
        gt_thresh=thresh+2*spacing(thresh)
        rarray1(0,:)=gt_thresh
        rarray1(ncol+1,:)=gt_thresh
        rarray1(:,0)=gt_thresh
        rarray1(:,nrow+1)=gt_thresh
        do irow=1,nrow
          do icol=1,ncol
            rarray1(icol,irow)=rarray(icol,irow)
          end do
        end do
        do i=1,num_bore_list
          call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
          fac3(i),fac4(i),icellno(i),jcellno(i),bore_val(i),   &
          rarray1)
        end do

        do i=1,num_bore_list
          if(bore_val(i).gt.7.0e37)then
            aval='not_in_grid'
          else if (bore_val(i).gt.4.0e37)then
            aval='inactive'
          else
            write(aval,'(1pg14.7)') bore_val(i)
          end if
          aval=adjustl(aval)
          write(outunit,1100) trim(bore_list_id(i)),east(i),north(i),trim(aval)
1100      format(1x,a,t14,2x,f15.2,2x,f15.2,2x,a)
	end do
1150	close(unit=outunit)
	write(6,1200) trim(outfile)
1200	format(' - file ',a,' written ok.')

1300	continue
	go to 9900



9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890

9890	call write_message(leadspace='yes')
9900	call close_files
	call free_bore_mem
	call free_grid_mem(gridspec)
	deallocate(icellno,jcellno,fac1,fac2,fac3,fac4,bore_val,rarray,rarray1,east,north,stat=ierr)
	write(6,*)

end program arr2bore

