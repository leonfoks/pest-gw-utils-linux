!     Last change:  J     2 May 2002    8:52 pm
program mod2smp

! -- Program MOD2SMP writes a bore sample file of model-generated heads over
!    time, interpolated to the sites of user-supplied bores.


	use defn
	use inter

	implicit none

	integer      :: ifail,ierr,ncol,nrow,i,j,modunit,dds,mms,yys,hhhs,&
	    		mmms,ssss,outunit,ntime,iarray,mcol,mrow,kstp,kper,&
			ntrans,ilay,kperold,ntransold,kstpold,itime,ddf,mmf, &
			yyf,hhhf,mmmf,sssf,icol,irow,idate,iheader
	real	     :: thresh,gt_thresh,day_convert,pertim,totim
	integer, allocatable, dimension(:)         :: layer,icellno,jcellno
	real, allocatable, dimension(:,:)          :: rarray,boreinfo
	real, allocatable, dimension(:)	           :: fac1,fac2,fac3,fac4, &
						      total_time
	double precision, allocatable, dimension(:):: east,north
	type (modelgrid) gridspec
	character (len=80)  :: modfile,outfile
	character (len=1)   :: af,at,ax
	character (len=15)  :: adate,atime,anum1,anum2
	character (len=16)  :: text


	write(amessage,5)
5       format(' Program MOD2SMP writes a bore sample file of model-generated ',&
	'heads over time, interpolated to the sites of user-supplied bores.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'
	
	kper=0
	ntrans=0
	kstp=0
	kperold=0
	ntransold=0
	kstpold=0

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
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
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

30      call read_bore_coord_file(ifail, &
	' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 10
	end if
	do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) go to 80
	end do
	go to 100
80	write(amessage,90) trim(bore_coord_file)
90	format(' The following bores were not provided with layer numbers in ',&
	'bore coordinates file ',a,' ----->')
	call write_message(leadspace='yes')
	j=2
	imessage=0
	write_bores: do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) then
	    write(amessage(j:),'(a)') trim(bore_coord_id(i))
	    j=j+11
	    if(j.ge.71) then
	      call write_message(increment=1)
	      if(imessage.gt.12) goto 9900
	      j=2
	    end if
	  end if
	end do write_bores
	if(j.ne.2) call write_message
	go to 9900

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
	        write(amessage,110) trim(bore_list_file)
110		format(' Execution of program MOD2SMP cannot proceed as ',&
		'there are multiple occurrences of the same bore in bore ',&
		'listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	allocate(east(num_bore_list), north(num_bore_list), &
	layer(num_bore_list),&
	fac1(num_bore_list), fac2(num_bore_list), fac3(num_bore_list), &
	fac4(num_bore_list), icellno(num_bore_list), jcellno(num_bore_list), &
	stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,130) trim(bore_list_file)
130	  format(' Too many bores cited in bore listing file ',a, &
	  ': insufficient memory available to continue execution.')
	  go to 9890
	end if

	imessage=0
	list: do i=1,num_bore_list
	  do j=1,num_bore_coord
	    if(bore_coord_id(j).eq.bore_list_id(i)) then
	      east(i)=bore_coord_east(j)
	      north(i)=bore_coord_north(j)
	      layer(i)=bore_coord_layer(j)
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

150	write(6,*)
	call open_input_file(ifail, &
	' Enter name of unformatted model-generated file: ',modfile,modunit, &
	file_format='unformatted')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  deallocate(east,north,layer,fac1,fac2,fac3,fac4,icellno,jcellno, &
	  stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 100
	end if

170	write(6,180,advance='no')
180	format(' Is this a MODFLOW or MT3D file?  [f/t]: ')
	read(5,'(a)') af
	if(af.eq.' ') go to 170
	call casetrans(af,'lo')
	if(index(eschar,af).ne.0) then
	  close(unit=modunit)
	  go to 150
	end if
	if((af.ne.'f').and.(af.ne.'t')) go to 170

190	write(6,192,advance='no')
192	format(' How many different output times are represented in this file? ')
	i=key_read(ntime)
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 170
	else if(i.eq.-1) then
	  go to 190
	else if(i.ne.0) then
	  write(6,220)
	  go to 190
	end if
	if(pos_test(ntime,'number of output times').ne.0) go to 190
	allocate(boreinfo(num_bore_list,ntime),total_time(ntime),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,195) trim(bore_list_file)
195	  format(' Cannot allocate sufficient memory to store information ', &
	  'for all bores represented in bore listing file ',a,' for all ',&
	  'of these output times.')
	  call write_message (leadspace='yes')
	  write(amessage,196)
196	  format(' Enter a smaller number "n" to read the first "n" output ',&
	  'times, or re-run model with fewer output times.')
	  call write_message (leadspace='yes',endspace='yes')
	  if(allocated(boreinfo))deallocate(boreinfo,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  if(allocated(total_time)) deallocate(total_time,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 190
	end if
	boreinfo=3.1e37
	
200	write(6,210,advance='no')
210	format(' Enter blanking threshold value for arrays in this file: ')
	i=key_read(thresh)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  deallocate(boreinfo,total_time,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 190
	else if(i.eq.-1) then
	  go to 200
	else if(i.ne.0) then
	  write(6,220)
220	  format(' Illegal input  - try again.')
	  go to 200
	end if
	if(thresh.gt.1.0e37) then
	  write(amessage,230)
230	  format(' Threshold must be less than 1.0E37  - try again.')
	  call write_message
	  go to 200
	end if

250	write(6,260,advance='no')
260	format(' Enter time units used by model (yr/day/hr/min/sec) [y/d/h/m/s]: ')
	read(5,'(a)') at
	if(at.eq.' ') go to 250
	if(index(eschar,at).ne.0) then
	  write(6,*)
	  go to 200
	end if
	call casetrans(at,'lo')
	if(at.eq.'s') then
	  day_convert=1.0/86400.0
	else if(at.eq.'m') then
	  day_convert=1.0/1440.0
	else if(at.eq.'h') then
	  day_convert=1.0/24.0
	else if(at.eq.'d') then
	  day_convert=1.0
	else if(at.eq.'y') then
	  day_convert=365.25
	else
	  go to 250
	end if

300	write(6,*)
310	if(datespec.eq.1) then
	  write(6,320,advance='no')
320	  format(' Enter simulation starting date [dd/mm/yyyy]: ')
	else
	  write(6,321,advance='no')
321	  format(' Enter simulation starting date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') adate
	if(adate.eq.' ') go to 310
	adate=adjustl(adate)
	if(index(eschar,adate(1:2)).ne.0) then
	  write(6,*)
	  go to 250
	end if
	call char2date(ifail,adate,dds,mms,yys)
	if(ifail.ne.0)then
	  write(6,340)
340	  format(' Illegal date  - try again.')
	  go to 310
	end if

360	write(6,370,advance='no')
370	format(' Enter simulation starting time [hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 360
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) go to 300
	call char2time(ifail,atime,hhhs,mmms,ssss)
	if(ifail.ne.0)then
	  write(6,380)
380	  format(' Illegal time  - try again.')
	  go to 360
	end if

400	write(6,*)
410	call open_output_file(ifail, &
	' Enter name for bore sample output file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 360
	end if

	deallocate(bore_coord_id,bore_coord_east,bore_coord_north, &
	bore_coord_layer,stat=ierr)
	if(ierr.ne.0) go to 9300
	nullify(bore_coord_id,bore_coord_east,bore_coord_north, &
	bore_coord_layer)	
	allocate(rarray(0:ncol+1,0:nrow+1),stat=ierr)
	if(ierr.ne.0) then
	  write(6,480)
480	  format(' Insufficient available memory to continue MOD2SMP ',&
	  'execution: try using fewer bores, model output times, or a ',&
	  'smaller model grid.')
	  go to 9890
	end if

        gt_thresh=thresh+2*spacing(thresh)
        rarray(0,:)=gt_thresh
        rarray(ncol+1,:)=gt_thresh
        rarray(:,0)=gt_thresh
        rarray(:,nrow+1)=gt_thresh

	iarray=0
	itime=0
	read_an_array: do
	  iarray=iarray+1
	  if(iarray.eq.10) write(6,460)
460	  format(/,' Working .....')
	  if(af.eq.'f')then
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=1000) kstp,kper,pertim,totim, &
	    text,mcol,mrow,ilay
	  else
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=1000) ntrans,kstp,kper,&
	    totim,text,mcol,mrow,ilay
	  end if
	  if((mcol.eq.0).or.(mrow.eq.0)) go to 9100
	  if((mrow.ne.nrow).or.(mcol.ne.ncol)) then
	    call num2char(iarray,anum1)
	    write(amessage,450) trim(anum1),trim(modfile),trim(gridspec%specfile)
450	    format(' Number of rows and columns read from header to array ',&
	    'number ',a,' from model output file ',a,' does not agree with ',&
	    'grid specifications as read from grid specification file ',a,   &
            '. Alternatively, the unformatted file protocol might be a problem - ', &
            'try using the alternative version of this program.')
	    go to 9890
	  end if
	  if(iarray.ne.1)then
	    if((kstp.ne.kstpold).or.(ntrans.ne.ntransold).or. &
	    (kper.ne.kperold)) then
	      itime=itime+1
	      if(itime.gt.ntime) then
		call num2char(ntime,anum1)
		write(amessage,500) trim(modfile),trim(anum1),trim(anum1)
500	        format(' Warning: more times are represented in file ',a, &
		' than the ',a,' times nominated by the user. Only arrays ',&
		'representing the first ',a,' times have been read.')
		call write_message(leadspace='yes')
		itime=itime-1
		go to 1000
	      end if
	      total_time(itime)=totim
	    end if
	  else
	    itime=itime+1
	    total_time(itime)=totim
	  end if
	  read(modunit,err=9200,end=9250) ((rarray(icol,irow),icol=1,ncol), &
	  irow=1,nrow)
	  kstpold=kstp
	  kperold=kper
	  ntransold=ntrans
	  do i=1,num_bore_list
	    if(layer(i).eq.ilay) then
	      call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),fac3(i), &
	      fac4(i),icellno(i),jcellno(i),boreinfo(i,itime),rarray)
	    end if
	  end do
	end do read_an_array

1000	if(iarray.eq.1) then
	  write(amessage,1010) trim(modfile)
1010	  format(' No arrays were found in model-generated unformatted ',&
	  'file ',a)
	  go to 9890
	end if
	write(6,*)
	close(unit=modunit)
	iarray=iarray-1
	ntime=itime
	call num2char(iarray,anum1)
	call num2char(itime,anum2)
	write(amessage,1050) trim(anum1),trim(anum2),trim(modfile)
1050	format(' - ',a,' arrays, covering ',a,' different model output ',&
	'times, read from file ',a)
	call write_message

	do i=1,num_bore_list
	  do itime=1,ntime
	    call elapsdate(total_time(itime),day_convert,dds,mms,yys,hhhs,&
	    mmms,ssss,ddf,mmf,yyf,hhhf,mmmf,sssf)
	    if(icellno(i).eq.-999) boreinfo(i,itime)=7.1e37
	    if(boreinfo(i,itime).gt.1.0e37) then
	      ax='x'
	    else
	      ax=' '
	    end if
	    if(datespec.eq.1) then
	      write(outunit,1100,err=9400) trim(bore_list_id(i)), ddf,mmf,yyf,&
	      hhhf,mmmf,sssf,boreinfo(i,itime),ax
	    else
	      write(outunit,1100,err=9400) trim(bore_list_id(i)), mmf,ddf,yyf,&
	      hhhf,mmmf,sssf,boreinfo(i,itime),ax
	    end if
1100	    format(1x,a,t14,i2.2,'/',i2.2,'/',i4.4,t28,i2.2,':',i2.2,':',i2.2,&
	    t38,1pg14.7,t55,a)
	  end do
	end do
	close(unit=outunit)
	write(amessage,1200) trim(outfile)
1200	format(' - bore sample file ',a,' written ok.')
	call write_message
	go to 9900

9100	call num2char(iarray,anum1)
	write(amessage,9110) trim(anum1),trim(modfile)
9110	format(' Error reading header to array number ',a,&
	' from model-generated file ',a)
	go to 9890
9200	call num2char(iarray,anum1)
	write(amessage,9210) trim(anum1),trim(modfile)
9210	format(' Error reading array number ',a,' from model ',&
	'output file ',a)
	go to 9890
9250	call num2char(iarray,anum1)
	write(amessage,9260) trim(anum1),trim(modfile)
9260	format(' Unexpected end of file encountered while reading array ',&
	'number ',a,' from model output file ',a)
	go to 9890
9300	write(amessage,9310)
9310	format(' Memory management error: cannot continue execution.')
	go to 9890
9400	write(amessage,9410) trim(outfile)
9410	format(' Error writing data to bore sample file ',a,&
	': file inaccessible or disk full.')
	go to 9890


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_bore_mem
	call free_grid_mem(gridspec)
	deallocate(rarray,east,north,fac1,fac2,fac3,fac4,icellno, &
	jcellno,boreinfo,total_time,stat=ierr)
	write(6,*)

end program mod2smp
 