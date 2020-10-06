!     Last change:  JD    9 May 2003   11:52 am
program mod2obs1

! -- Program MOD2OBS writes a bore sample file of model-generated heads over
!    time, interpolated to the sites and times of samples supplied in an 
!    existing bore sample file.

	use defn
	use inter

	implicit none

	logical      :: active
	integer      :: ifail,ierr,ncol,nrow,i,j,modunit,dds,mms,yys,hhhs,&
	    		mmms,ssss,outunit,iarray,mcol,mrow,kstp,kper,&
			ntrans,ilay,cols,ndays,nsecs,k,itry, &
			icol,irow,idate,sampunit,ilt0,iheader,   &
                        ndaystart,nsecstart,nlay,iline,numsamp1,jindex, &
                        nindex,day,mon,year,hour,sec,min,iflag,kk,kmin,ii,nind1
	real	     :: thresh,gt_thresh,day_convert,pertim,totim,  &
	                r_day_convert,bore_time,rtemp1,rtemp2,rlim,rrlim
	integer, allocatable, dimension(:)         :: layer,icellno,jcellno, &
                                                      numsamp,index_bore,iorder
	real, allocatable, dimension(:,:,:)        :: rarray1,rarray2
	real, allocatable, dimension(:)	           :: fac1,fac2,fac3,fac4
	integer, allocatable, dimension(:)         :: bore_days(:),bore_secs(:)
	real, allocatable, dimension(:)            :: bore_val(:)
	real, allocatable, dimension (:)           :: arrtime1(:),arrtime2(:)
	double precision                           :: dtemp
	double precision, allocatable, dimension(:):: east,north
	type (modelgrid) gridspec
        character (len=1)   :: abin,about
	character (len=10)  :: abore, aboreold, alay, alim
        character (len=12)  :: siteid
        character (len=23)  :: bcode
	character (len=80)  :: modfile,outfile,sampfile
	character (len=1)   :: af,at
	character (len=15)  :: adate,atime,anum1
	character (len=16)  :: text
	character (len=20)  :: aval


	write(amessage,5)
5       format(' Program MOD2OBS1 writes an ascii or binary  bore sample file of ', &
        'model-generated heads over time, interpolated to the sites and times ',   &
        'supplied in an existing ascii or binary bore sample file.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'
	
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

!	open(unit=75,file='debug.dat')		!debug

	call readfig(gridspec%specfile,bore_coord_file,sampfile)
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
110		format(' Execution of program MOD2OBS cannot proceed as ',&
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
	numsamp(num_bore_list), index_bore(num_bore_list), &
        iorder(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,130)
130	  format(' Insufficient memory available to continue execution.')
	  go to 9890
	end if
        index_bore=0

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

599     itry=0
600     continue
        if(itry.eq.5) go to 9900
        write(6,605,advance='no')
605     format(' Enter name of bore sample file: ')
        read(5,'(a)') sampfile
        if(sampfile.eq.' ') go to 600
        if(index(eschar,sampfile(1:2)).ne.0)then
          write(6,*)
	  deallocate(east,north,layer,fac1,fac2,fac3,fac4,icellno,jcellno,  &
          numsamp,index_bore,iorder,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 100
	end if
606     write(6,607,advance='no')
607     format(' Is this an ASCII or binary file?  [a/b]: ')
        read(5,'(a)') abin
        if(abin.eq.' ') go to 606
        if((abin.eq.'E').or.(abin.eq.'e'))then
          write(6,*)
          go to 599
        end if
        if((abin.eq.'A').or.(abin.eq.'a'))then
          abin='a'
        else if((abin.eq.'B').or.(abin.eq.'b'))then
          abin='b'
        else
          go to 606
        end if
        sampunit=nextunit()
        if(abin.eq.'a')then
          open(unit=sampunit,file=sampfile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,608) trim(sampfile)
608         format(' Cannot open file ',a,' as a formatted site sample file.')
            itry=itry+1
            go to 600
          end if
        else
          open(unit=sampunit,file=sampfile,status='old',form='binary',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,611) trim(sampfile)
611         format(' Cannot open file ',a,' as a binary site sample file.')
            itry=itry+1
            go to 600
          end if
          read(sampunit,iostat=ierr) bcode
          if(ierr.ne.0)then
            write(amessage,612) trim(sampfile)
612          format(' Cannot read header to binary site sample file ',a,  &
            ': are you sure that this is a file of the correct type? Try again.')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            close(unit=sampunit)
            go to 600
          end if
          if(bcode.ne.'binary_site_sample_file')then
            write(amessage,613) trim(sampfile)
613          format(' Incorrect header to file ',a,': are you sure that this is a ', &
            'binary site sample file? Try again.')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            close(unit=sampunit)
            go to 600
          end if
        end if

150	write(6,*)
	call open_input_file(ifail, &
	' Enter name of unformatted model-generated file: ',modfile,modunit, &
	file_format='unformatted')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  close(unit=sampunit)
	  go to 600
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

200	write(6,210,advance='no')
210	format(' Enter inactive threshold value for numbers in this file: ')
	i=key_read(thresh)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 170
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
	r_day_convert=1.0/day_convert

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
	ndaystart=numdays(1,1,1970,dds,mms,yys)
	nsecstart=numsecs(0,0,0,hhhs,mmms,ssss)
390	write(6,*)
395	write(6,396,advance='no')
396	format(' How many layers in the model? ')
	read(5,'(a)') alay
	if(alay.eq.' ') go to 395
	alay=adjustl(alay)
	if(index(eschar,alay(1:2)).ne.0) then
	  write(6,*)
	  go to 360
	end if
	call char2num(ifail,alay,nlay)
	if(ifail.ne.0) go to 395
	if(nlay.le.0) go to 395

440	write(amessage,445)
445	format(' If a sample time does not lie between model output ',  &
        'times, or if there is only one model output time, value at ', &
	'the sample time can equal that at nearest model output time:-')
	call write_message(leadspace='yes')
465	write(6,467,advance='no')
467	format('   Enter extrapolation limit in days (fractional if ', &
        'necessary): ')
	read(5,'(a)') alim
	if(alim.eq.' ') go to 465
	alim=adjustl(alim)
	if(index(eschar,alim(1:2)).ne.0) then
	  go to 390
	end if
	call char2num(ifail,alim,rlim)
	if(ifail.ne.0) go to 465
	if(rlim.lt.0.0) go to 465
	rrlim=rlim
	rlim=rlim*r_day_convert

400	write(6,*)
399     itry=0
401     continue
        if(itry.ge.5) go to 9900
402     write(6,403,advance='no')
403     format(' Enter name for bore sample output file: ')
        read(5,'(a)') outfile
        if(outfile.eq.' ') go to 402
        if(index(eschar,outfile(1:2)).ne.0)then
          write(6,*)
          go to 440
        end if
404     write(6,405,advance='no')
405     format(' Is this an ASCII or binary file?  [a/b]: ')
        read(5,'(a)') about
        if(about.eq.' ') go to 404
        if((about.eq.'E').or.(about.eq.'e'))then
          write(6,*)
          go to 399
        end if
        if((about.eq.'A').or.(about.eq.'a'))then
          about='a'
        else if((about.eq.'B').or.(about.eq.'b'))then
          about='b'
        else
          go to 404
        end if
        outunit=nextunit()
        if(about.eq.'a')then
          open(unit=outunit,file=outfile,iostat=ierr)
        else if(about.eq.'b')then
          open(unit=outunit,file=outfile,form='binary',iostat=ierr)
        end if
        if(ierr.ne.0)then
          write(amessage,407) trim(outfile)
407       format(' Cannot open file ',a,' for output - try again.')
          call write_message()
          itry=itry+1
          go to 401
        end if
        if(about.eq.'b')then
          write(outunit) 'binary_site_sample_file'
        end if

	deallocate(bore_coord_id,bore_coord_east,bore_coord_north, &
	bore_coord_layer,stat=ierr)
	if(ierr.ne.0) go to 9300
	nullify(bore_coord_id,bore_coord_east,bore_coord_north, &
	bore_coord_layer)	
	allocate(rarray1(0:ncol+1,0:nrow+1,nlay),                  &
        rarray2(0:ncol+1,0:nrow+1,nlay),arrtime1(nlay),arrtime2(nlay),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,480)
480	  format(' Insufficient available memory to continue MOD2OBS ',&
	  'execution: try using fewer bores, a smaller bore sample file, or a ',&
	  'smaller model grid.')
	  go to 9890
	end if

	aboreold=' '
	active=.false.
	iline=0
	imessage=0
	numsamp=0
	jindex=0
        read_sample: do
          iline=iline+1
          if(iline.eq.2000) write(6,590)
590       format(/,' Working.....')
          if(abin.eq.'a')then
            read(sampunit,'(a)',err=9600,end=700) cline
            cols=5
            call linesplit(ifail,5)
            if(ifail.lt.0) cycle read_sample
            if(ifail.gt.0) then
              call linesplit(ifail,4)
              if(ifail.ne.0) go to 9150
              cols=4
            end if
            if(right_word(1)-left_word(1).gt.9) go to 9500
            abore=cline(left_word(1):right_word(1))
          else
            read(sampunit,err=9600,end=700)  siteid,ndays,nsecs,dtemp
            abore=siteid
          end if
          call casetrans(abore,'hi')
	  if(abore.ne.aboreold)then
	    aboreold=abore
	    if(active) numsamp(i)=numsamp1
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) go to 620
	    end do
	    active=.false.
	    cycle read_sample
620	    numsamp1=0
	    active=.true.
	    index_bore(i)=jindex+1
	  else
	    if(.not.active) cycle read_sample
	  end if
          if(abin.eq.'a')then
            call read_rest_of_sample_line(ifail,cols,ndays,nsecs,dtemp, &
            iline,sampfile)
            if(ifail.ne.0) go to 9900
	    if(dtemp.lt.-1.0e38) cycle read_sample
          end if
          if(ndays.lt.ndaystart-rrlim-1) cycle read_sample
	  jindex=jindex+1
	  numsamp1=numsamp1+1
	end do read_sample

700	if(active) numsamp(i)=numsamp1
	nindex=jindex
	allocate(bore_days(nindex),bore_secs(nindex),bore_val(nindex), &
        stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,480)
	  go to 9890
	end if
	bore_val=-1.1e38
	nind1=nindex

! -- Now we re-read the bore sample file to read in data.

	rewind(unit=sampunit,iostat=ierr)
	if(ierr.ne.0)then
	  write(amessage,720) trim(sampfile)
720	  format(' Error re-winding bore sample file ',a)
	  go to 9890
	end if
        if(abin.eq.'b')then
          read(sampunit,iostat=ierr) bcode
        end if
	aboreold=' '
	jindex=0
	iline=0
	read_sample_1: do
	  iline=iline+1
          if(abin.eq.'a')then
	    read(sampunit,'(a)',err=9600,end=800) cline
	    cols=5
	    call linesplit(ifail,5)
	    if(ifail.lt.0) cycle read_sample_1
	    if(ifail.gt.0)then
	      call linesplit(ifail,4)
	      if(ifail.ne.0) go to 9150
	      cols=4
	    end if
	    abore=cline(left_word(1):right_word(1))
          else
            read(sampunit,err=9600,end=800)  siteid,ndays,nsecs,dtemp
            abore=siteid
          end if
	  call casetrans(abore,'hi')
	  if(abore.ne.aboreold)then
	    aboreold=abore
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) go to 820
	    end do
	    active=.false.
	    cycle read_sample_1
820	    active=.true.
	  else
	    if(.not.active) cycle read_sample_1
	  end if
          if(abin.eq.'a')then
	    call read_rest_of_sample_line(ifail,cols,ndays,nsecs,dtemp,iline, &
	    sampfile)
	    if(ifail.ne.0) go to 9900
	    if(dtemp.lt.-1.0e38) cycle read_sample_1
          end if
	  if(ndays.lt.ndaystart-rrlim-1) cycle read_sample_1
	  jindex=jindex+1
	  bore_days(jindex)=ndays
	  bore_secs(jindex)=nsecs
	end do read_sample_1
800	continue

!	do i=1,nindex					!debug
!	  write(75,*) bore_days(i),bore_secs(i)		!debug
!	end do						!debug
!	write(75,*)					!debug
!	do i=1,num_bore_list				!debug
!	  write(75,*) index_bore(i), numsamp(i)		!debug
!	end do						!debug

	rarray2=0.0
        gt_thresh=thresh+2*spacing(thresh)
        rarray2(0,:,:)=gt_thresh
        rarray2(ncol+1,:,:)=gt_thresh
        rarray2(:,0,:)=gt_thresh
        rarray2(:,nrow+1,:)=gt_thresh

	iarray=0
	arrtime1=0.0
	arrtime2=0.0
	read_an_array: do
	  iarray=iarray+1
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
	  if((mcol.le.0).or.(mrow.le.0).or.(ilay.le.0)) go to 9100
	  if((mrow.ne.nrow).or.(mcol.ne.ncol)) then
	    call num2char(iarray,anum1)
	    write(amessage,450) trim(anum1),trim(modfile),trim(gridspec%specfile)
450	    format(' Number of rows and columns read from header to array ',&
	    'number ',a,' from model output file ',a,' does not agree with ',&
	    'grid specifications as read from grid specification file ',a,   &
            '. Alternatively, unformatted file protocol might be a problem - ', &
            'try using the alternative version of this program.')
	    go to 9890
	  end if
	  if(ilay.gt.nlay)then
	    call num2char(iarray,anum1)
	    write(amessage,460) trim(anum1),trim(modfile)
460	    format(' Number of layers read from header to array number ', &
	    a,' from model output file ',a,' exceeds number of model ', &
	    'layers as supplied by user.')
	    go to 9890
	  end if

	  rarray1(:,:,ilay)=rarray2(:,:,ilay)
	  read(modunit,err=9200,end=9250)   &
          ((rarray2(icol,irow,ilay),icol=1,ncol),irow=1,nrow)
	  arrtime1(ilay)=arrtime2(ilay)
	  arrtime2(ilay)=totim
	  if(arrtime1(ilay).eq.0.0) then
	    interp: do i=1,num_bore_list
	      if(layer(i).ne.ilay) cycle interp
	      if(numsamp(i).eq.0) cycle interp
	      jindex=index_bore(i)
	      do j=jindex,jindex+numsamp(i)-1
	        bore_time=(bore_days(j)-ndaystart)*r_day_convert +       &
	        (bore_secs(j)-nsecstart)/86400.0*r_day_convert
	        if((bore_time.ge.arrtime2(ilay)-rlim).and.               &
	           (bore_time.le.arrtime2(ilay)+rlim))then
	           call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
	           fac3(i),fac4(i),icellno(i),jcellno(i),bore_val(j),   &
	           rarray2(0,0,ilay))
	        end if
	      end do
	    end do interp
	    cycle read_an_array
	  end if
	  interpolate: do i=1,num_bore_list
	    if(layer(i).ne.ilay) cycle interpolate
	    if(numsamp(i).eq.0) cycle interpolate
	    jindex=index_bore(i)
	    do j=jindex,jindex+numsamp(i)-1
!             if(bore_val(j).gt.-1.0e38) cycle 
	      bore_time=(bore_days(j)-ndaystart)*r_day_convert +       &
              (bore_secs(j)-nsecstart)/86400.0*r_day_convert
	      if(bore_time.gt.arrtime2(ilay)+rlim)then
	        continue
	      else if(bore_time.lt.arrtime1(ilay)-rlim)then
	        continue
	      else if(bore_time.lt.arrtime1(ilay))then
	        continue
!	        call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
!                fac3(i),fac4(i),icellno(i),jcellno(i),bore_val(j),   &
!	        rarray1(0,0,ilay))
	      else if(bore_time.ge.arrtime2(ilay))then
	        call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
                fac3(i),fac4(i),icellno(i),jcellno(i),bore_val(j),   &
	        rarray2(0,0,ilay))
	      else
	        call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
                fac3(i),fac4(i),icellno(i),jcellno(i),rtemp1,        &
	        rarray1(0,0,ilay))
	        call point_interp(ncol,nrow,thresh,fac1(i),fac2(i),  &
                fac3(i),fac4(i),icellno(i),jcellno(i),rtemp2,        &
	        rarray2(0,0,ilay))
	        bore_val(j)=rtemp1+(rtemp2-rtemp1)*(bore_time-arrtime1(ilay)) &
	        /(arrtime2(ilay)-arrtime1(ilay))
	      end if
	    end do
	  end do interpolate

	end do read_an_array

!note: this will not work if there is only one array per layer. There
!have to be many arrays per layers corresponding to different times.

1000	if(iarray.eq.1) then
	  write(amessage,1010) trim(modfile)
1010	  format(' No arrays were found in model-generated unformatted ',&
	  'file ',a)
	  go to 9890
	end if
	write(6,*)
	close(unit=modunit)
	iarray=iarray-1
	call num2char(iarray,anum1)
	write(amessage,1050) trim(anum1),trim(modfile)
1050	format(' - ',a,' arrays read from file ',a)
	call write_message

	iorder=index_bore
	do ii=1,num_bore_list

	  kmin=nind1+10
	  ilt0=0
          do kk=1,num_bore_list
	    if(numsamp(kk).eq.0) cycle
	    if(iorder(kk).lt.0) cycle
	    ilt0=ilt0+1
	    if(iorder(kk).lt.kmin)then
	      kmin=iorder(kk)
	      i=kk
	    end if
	  end do
	  if(ilt0.eq.0) go to 1150
!	  write(75,*) ilt0				!debug
	  iorder(i)=-1

	  if(numsamp(i).eq.0) cycle
	  jindex=index_bore(i)
	  do j=jindex,jindex+numsamp(i)-1
	  if(bore_val(j).lt.-1.0e38)then
	    cycle
!	    aval='cannot_interpolate'
	  else if(bore_val(j).gt.7.0e37)then
	    aval='not_in_grid'
	  else if (bore_val(j).gt.4.0e37)then
	    aval='dry_or_inactive'
	  else
	    write(aval,'(1pg14.7)') bore_val(j)
	  end if
	  aval=adjustl(aval)
	  call newdate(bore_days(j),1,1,1970,day,mon,year)
	  hour=bore_secs(j)/3600
	  sec=bore_secs(j)-hour*3600
	  min=sec/60
	  sec=sec-min*60
          if(about.eq.'a')then
	    if(datespec.eq.1)then
	      write(outunit,1100,err=9400) trim(bore_list_id(i)), &
	      day,mon,year,hour,min,sec,trim(aval)
	    else
	      write(outunit,1100,err=9400) trim(bore_list_id(i)), &
	      mon,day,year,hour,min,sec,trim(aval)
	    end if
1100	    format(1x,a,t14,i2.2,'/',i2.2,'/',i4.4,t28,i2.2,':',i2.2,':',i2.2,&
	    t40,a)
          else
            siteid=trim(bore_list_id(i))
            dtemp=bore_val(j)
            write(outunit)siteid,bore_days(j),bore_secs(j),dtemp
          end if
	  end do
	end do
1150	close(unit=outunit)
	write(amessage,1200) trim(outfile)
1200	format(' - bore sample file ',a,' written ok.')
	call write_message

	iflag=0
	if(any(numsamp.eq.0))then
	  iflag=1
	  go to 1249
	end if
	do i=1,num_bore_list
	  if(numsamp(i).ne.0) then
	    do j=index_bore(i),index_bore(i)+numsamp(i)-1
	      if(bore_val(j).gt.-1e38) goto 1230
	    end do
	    go to 1240
1230	  continue
	  end if
	end do
	go to 1249
1240	iflag=1
1249	continue
	if(iflag.eq.0) go to 1300

	write(amessage,1250) trim(bore_list_file),trim(outfile)
1250	format(' The following bores appearing in the bore listing ', &
	'file ',a,' will not appear in the bore sample file ',a,':-')
	call write_message(leadspace='yes')
	amessage=' '
        j=4
        imessage=0
        write_bores_1: do i=1,num_bore_list
          if(numsamp(i).eq.0) then
            write(amessage(j:),'(a)') trim(bore_list_id(i))
            j=j+11
            if(j.ge.69) then
              call write_message(increment=1)
              if(imessage.gt.12) goto 9900
              j=4
            end if
          end if
        end do write_bores_1
	write_bores_2: do i=1,num_bore_list
	  if(numsamp(i).ne.0)then
	    do k=index_bore(i),index_bore(i)+numsamp(i)-1
	      if(bore_val(k).gt.-1.0e38) go to 1290
	    end do
            write(amessage(j:),'(a)') trim(bore_list_id(i))
            j=j+11
            if(j.ge.69) then
              call write_message(increment=1)
              if(imessage.gt.12) goto 9900
              j=4
            end if
1290	    continue
	  endif
	end do write_bores_2
        if(j.ne.4) call write_message

1300	continue
	go to 9900

9100	call num2char(iarray,anum1)
	write(amessage,9110) trim(anum1),trim(modfile)
9110	format(' Error reading header to array number ',a,&
	' from model-generated file ',a)
	go to 9890
9150    call num2char(iline,anum1)
        write(amessage,9160) trim(anum1),trim(sampfile)
9160    format(' Error on line ',a,' of bore sample file ',a,': ',&
        'insufficient entries.')
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
9500    call num2char(iline,anum1)
        write(amessage,9510) trim(anum1),trim(sampfile)
9510    format(' Error on line ',a,' of bore sample file ',a,': bore ',&
        'identifier greater than 10 characters in length.')
        go to 9890
9600    call num2char(iline,anum1)
        write(amessage,9610) trim(anum1),trim(sampfile)
9610    format(' Error reading line ',a,' of bore sample file ',a)
        go to 9890


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_bore_mem
	call free_grid_mem(gridspec)
	deallocate(rarray1,rarray2,east,north,fac1,fac2,fac3,fac4,icellno, &
	jcellno,layer,numsamp,index_bore,arrtime1,arrtime2, &
	bore_days,bore_secs,bore_val,iorder,stat=ierr)
	write(6,*)

end program mod2obs1


 