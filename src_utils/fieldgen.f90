!     Last change:  JD    9 May 2003    1:56 pm
program fieldgen

! -- Program FIELDGEN generates 2D stochastic fields in model domain zones
!    using the Sequential Gaussian Simulation method.

	use defn
	use inter

	implicit none

        integer     :: ifail,idate,iheader,i,ncol,nrow,icol,irow,ierr,structunit, &
                       numint,currint,newint,ihuge,itemp,numstruct,numvario,      &
                       nbb,j,izone,ndat,istruct,k_ktype,n_nst,itrans,             &
                       irealization,isim,nparmcall,idim,maxtempdim,jerr,          &
                       irowmin,icolmin,irowmax,icolmax,n_nx,n_ny,iflag,uniform,   &
                       ind,ndig,outunit,ncond,n_nd,lt,npptzone,nz,ir,ic,ix,iy,nn, &
                       max_nx,max_ny,iunif,iitemp,icount,iseed,kstp,kper,ilay
        real        :: c_c0,rtemp1,rtemp2,eps,x_xmn,y_ymn,x_xsiz,y_ysiz,astdev,   &
                       delrmin,delcmin,xmax,ymax,xu,yu,xx,yy,rsum,ecrd,ncrd,      &
                       rmindist,rminval,dist,pertim,totim
        double precision :: minsep,minsep2,eastdiff,northdiff,distance,etemp,ntemp
        integer     :: iloc(1)
        integer     :: i_it(MAX_STRUCT_VARIO)
        real        :: c_cc(MAX_STRUCT_VARIO),a_ang(MAX_STRUCT_VARIO),  &
                       a_aa(MAX_STRUCT_VARIO),a_anis(MAX_STRUCT_VARIO)
        integer, allocatable, dimension(:)             :: ndmin,ndmax,ncnode,icc,irr
	integer, allocatable, dimension(:,:)	       :: intarray
        real, allocatable, dimension(:)                :: radmax
        real,allocatable, dimension(:)                 :: eastdat,northdat,valdat
        real,allocatable, dimension(:)                 :: temparray
        real, allocatable, dimension(:)                :: ec,nc
        real,allocatable,dimension(:,:)                :: rarray
        character (len=1),  allocatable, dimension(:)  :: akrig
        character (len=10), allocatable, dimension(:)  :: astructure
        character (len=1)                              :: as
        character (len=2)                              :: af
        character (len=10)                             :: ostring
        character (len=15)                             :: anum1,atemp,anum2
        character (len=16)                             :: text
	character (len=120)                            :: aprompt,intfile,structfile, &
                                                          atempf,arealfile,realbase

	type (modelgrid) gridspec
        type (geostructure), allocatable, dimension(:) :: structure(:)
        type (variogram), allocatable, dimension(:)    :: vario(:)


	integer, parameter :: MAXINT=100
	integer		   :: intval(MAXINT),iwork(MAXINT)
        real               :: avval(MAXINT)

!        open(unit=55,file='debug.dat')                     !debug

	write(amessage,5)
5	format(' Program FIELDGEN generates 2D stochastic fields in ',        &
               'model domain zones using the Sequential Gaussian ',           &
               'Simulation method.')

	call write_message(leadspace='yes',endspace='yes')

! -- Initialisation of variables used to write unformatted file.

        kstp=1
        kper=1
        pertim=1.0
        totim=1.0
        text=' '
        ilay=1

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig.')
	  call write_message
	  go to 9900
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig.')
	  call write_message
	  go to 9900
	end if

	call readfig(gridspec%specfile,pilotfile=pilot_points_file)
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
        delrmin=minval(gridspec%delr)
        delcmin=minval(gridspec%delc)
        rtemp1=delrmin
        eps=3.0*epsilon(rtemp1)
        do icol=1,ncol
          if(abs(rtemp1-gridspec%delr(icol)).gt.eps) go to 1310
        end do
        rtemp2=delcmin
        eps=3.0*epsilon(rtemp2)
        do irow=1,nrow
          if(abs(rtemp2-gridspec%delc(irow)).gt.eps) go to 1310
        end do
        iunif=1
        go to 1350
1310    iunif=0
1350    continue
        allocate(ec(ncol),nc(nrow),stat=ierr)
        if(ierr.ne.0)then
	  write(amessage,50)
	  go to 9890
	end if
        ec(1)=gridspec%delr(1)/2.0
        do icol=2,ncol
          ec(icol)=ec(icol-1)+(gridspec%delr(icol)+gridspec%delr(icol-1))/2.0
        end do
        nc(1)=-gridspec%delc(1)/2.0
        do irow=2,nrow
          nc(irow)=nc(irow-1)-(gridspec%delc(irow)+gridspec%delc(irow-1))/2.0
        end do

        ostring='0000000000'
        maxtempdim=ncol*nrow
	allocate(intarray(ncol,nrow),rarray(ncol,nrow),               &
        temparray(maxtempdim),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run FIELDGEN.')
	  go to 9890
	end if

30      call read_pilot_points_file(ifail, &
	' Enter name of conditioning pilot points file (<Enter> if none): ',  &
        accept_blank='yes')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(intarray,rarray,temparray,ec,nc,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
          call free_grid_mem(gridspec)
	  write(6,*)
	  go to 10
	end if

59      continue
        if(num_pilot_points.ne.0)then
!          write(6,*)
!61        write(6,63,advance='no')
!63        format(' Enter minimum allowable points separation: ')
!          if(key_read(minsep).ne.0) go to 61
!          if(escset.eq.1) then
!            write(6,*)
!            escset=0
!            go to 30
!          end if
          minsep=0.0d0
          imessage=0
          if(minsep.lt.0.0d0)minsep=0.0d0
!         if(minsep.gt.0.0d0)then
          minsep2=minsep*minsep
          if(num_pilot_points.gt.1)then
            do i=1,num_pilot_points-1
              do j=i+1,num_pilot_points
                eastdiff=pilot_point_east(i)-pilot_point_east(j)
                northdiff=pilot_point_north(i)-pilot_point_north(j)
                distance=eastdiff*eastdiff+northdiff*northdiff
                if((distance.le.minsep2).and. &
                   (pilot_point_zone(i).eq.pilot_point_zone(j)))then
                  imessage=imessage+1
                  if(imessage.gt.20) go to 9900
                  if(imessage.eq.1)then
                    write(amessage,65)
65                  format(' The following points are separated by zero distance ---->')
                    call write_message(leadspace='yes')
                  end if
                  write(amessage,66) trim(pilot_point_id(i)), &
                  trim(pilot_point_id(j))
66                format(1x,a,t15,a, t30, '(separation = 0.0)')
                  call write_message()
                end if
              end do
            end do
            end if
!          end if
          if(imessage.gt.0) go to 9900
        end if

        write(6,*)
70	continue
	aprompt=' Enter name of zonal integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
          write(6,*)
          go to 30
!          if(num_pilot_points.ne.0)then
!	    go to 59
!          else
!            write(6,*)
!            go to 30
!          end if
	end if
	intfile=aprompt

80	continue
	aprompt=' Enter name of structure file: '
	call open_input_file(ifail,aprompt,structfile,structunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
	  go to 70
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
	      write(amessage,150) trim(anum1),trim(intfile)
150	      format(' The zone integer array can hold a maximum of ',a, &
	      ' different integers for proper FIELDGEN execution. The array ', &
	      ' contained in file ',a,' holds more than this. To increase ',&
	      'this limit, edit FIELDGEN source code, increase parameter ', &
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

! -- The integers are now presented to the user. But first zonal arrays
!    are allocated.

        allocate(astructure(numint),radmax(numint), &
        akrig(numint),ndmin(numint),ndmax(numint),ncnode(numint),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

        radmax=0.0          ! radmax is an array
        akrig=' '
        astructure=' '

200 	write(6,*)
205	write(6,210)
210	format(' The following zones have been detected in the integer ',&
	'array:-')
	i=0
255	i=i+1
	if(i.gt.numint) go to 350
        write(6,*)
259	call num2char(intval(i),anum1)
260	write(6,270) trim(anum1)
270	format('    For zone characterised by integer value of ',a,':- ')
        atemp=' '
275     write(6,280,advance='no')
280     format('    Enter structure name (blank if no field generation for this zone): ')
        read(5,'(a)') atemp
        if(index(eschar,atemp(1:2)).ne.0)then
	  if(i.gt.2) then
	    i=i-2
	    go to 255
	  else if(i.eq.2) then
	    go to 200
	  else
            deallocate(astructure,radmax,akrig,ndmin,ndmax,ncnode,stat=ierr)
	    if(ierr.ne.0) then
	      write(amessage,57)
	      go to 9890
	    end if
            close(unit=structunit)
            write(6,*)
	    go to 80
	  end if
	end if
        atemp=adjustl(atemp)
        if(len_trim(atemp).gt.10)then
          write(amessage,285)
285       format('    Structure name must be 10 characters or less - try again.')
          call write_message()
          go to 275
        end if
        astructure(i)=atemp
        if(astructure(i).eq.' ') go to 255
        call casetrans(astructure(i),'lo')
289     write(6,290,advance='no')
290     format('    Use simple or ordinary kriging [s/o] in field generation: ')
	read(5,'(a)') akrig(i)
	if(akrig(i).eq.' ') go to 289
	if(index(eschar,akrig(i)).ne.0) then
	  write(6,*)
	  go to 259
	end if
	call casetrans(akrig(i),'lo')
	if((akrig(i).ne.'s').and.(akrig(i).ne.'o')) go to 289
!299     write(6,300,advance='no')
!300     format('    Enter search radius: ')
!        itemp=key_read(radmax(i))
!	if(escset.eq.1) then
!	  escset=0
!	  write(6,*)
!          go to 289
!        end if
!	if(itemp.ne.0) then
!          write(6,310)
310	  format('    Data input error  - try again.')
!	  go to 299
!	end if
!        if(radmax(i).le.0.0d0)then
!          write(6,311)
311       format('    Must be greater than zero - try again.')
!          go to 299
!        end if
        radmax(i)=1.0e15
        if(radmax(i).gt.1.0e15)radmax(i)=1.0e15

328     continue
        if(num_pilot_points.ne.0)then
          do j=1,num_pilot_points
            if(intval(i).eq.pilot_point_zone(j)) go to 336
          end do
          ndmax(i)=0
          ndmin(i)=0
          iflag=0
          go to 342
336       continue
          iflag=1
!329       write(6,330,advance='no')
!330       format('    Enter minimum number of conditioning points to use: ')
!          itemp=key_read(ndmin(i))
!	  if(escset.eq.1) then
!	    escset=0
!	    write(6,*)
!            go to 299
!          end if
! 	  if(itemp.ne.0) then
!            write(6,310)
!	    go to 329
!	  end if
!          if(ndmin(i).lt.0)ndmin(i)=0
          ndmin(i)=0           ! use instead of question
335       write(6,337,advance='no')
337       format('    Enter maximum number of conditioning points to use: ')
          itemp=key_read(ndmax(i))
	  if(escset.eq.1) then
	    escset=0
	    write(6,*)
            go to 289
          end if
	  if(itemp.ne.0) then
            write(6,310)
	    go to 335
	  end if
!          if(ndmax(i).lt.ndmin(i))then
          if(ndmax(i).le.0)then
            write(6,340)
340         format('    Must be greater than zero - try again.')
            go to 335
          end if
          if(ndmax(i).gt.100)then
            write(6,341)
341         format('    Must not be greater than 100 - try again.')
            go to 335
          end if
        else
          ndmin(i)=0
          ndmax(i)=0
        end if
342     continue
        write(6,343,advance='no')
343     format('    Enter maximum number of previously simulated nodes to use: ')
        itemp=key_read(ncnode(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          if(num_pilot_points.ne.0)then
            if(iflag.eq.1)then
              go to 328
            else
              go to 289
            end if
          else
            go to 289
          end if
        end if
        if(itemp.ne.0) then
          write(6,310)
          go to 342
        end if
        if(ncnode(i).le.0)then
          write(6,344)
344       format('    Must be greater than zero - try again.')
          go to 342
        end if
        if(ncnode(i).gt.100)then
          write(6,341)
          go to 342
        end if
        go to 255

350	continue

        do i=1,numint
          if(astructure(i).ne.' ')go to 351
        end do
        write(amessage,352)
352     format(' Field generation must take place for at least one zone in integer ', &
        'array - try again.')
        call write_message(leadspace='yes',endspace='yes')
        go to 200
351     continue

        write(6,*)
361     write(6,362,advance='no')
362     format(' How many realizations do you wish to generate? ')
        itemp=key_read(irealization)
        if(escset.ne.0)then
          escset=0
          i=numint-1
          go to 255
        end if
        if(itemp.ne.0)then
          write(6,310)
          go to 361
        end if
        if(irealization.lt.1)then
          write(6,311)
          go to 361
        end if
363     write(6,364,advance='no')
364     format(' Enter filename base for real array files: ')
        read(5,'(a)') atempf
        if(atempf.eq.' ') go to 363
        atempf=adjustl(atempf)
        if(index(eschar,atempf(1:2)).ne.0)then
          write(6,*)
          go to 361
        end if
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,realbase,1,nbb)
        if(ifail.ne.0) go to 363
365     write(6,366,advance='no')
366     format(' Write formatted or unformatted files? [f/u]: ')
        read(5,'(a)') af
        if(af.eq.' ') go to 365
        af=adjustl(af)
        call casetrans(af,'lo')
        if(af.eq.'e ') go to 363
        if((af.ne.'f ').and.(af.ne.'u '))go to 365

        as=' '
1369    if(iunif.eq.0)then
1370      write(6,1380,advance='no')
1380      format(/,' The model grid is non-uniform. To convert fields ', &
          'from stochastic subgrid to ',/,' model grid do you wish to average ', &
          'subgrid nodes or use closest subgrid ',/,' node?  [a/c]: ')
          read(5,'(a)') as
          if(as.eq.' ') go to 1370
          if((as.eq.'e').or.(as.eq.'E'))then
            write(6,*)
            go to 363
          end if
          if(as.eq.'A') as='a'
          if(as.eq.'C') as='c'
          if((as.ne.'a').and.(as.ne.'c')) go to 1370
        end if

! -- The geostatistical structure file is perused a first time to ascertain
!    array dimensions.

        call read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)
        if(ifail.ne.0) go to 9900
        if(numstruct.eq.0)then
          write(amessage,370) trim(structfile)
370       format(' No geostatistical structures found in file ',a)
          go to 9890
        end if
        if(numvario.eq.0)then
          write(amessage,380) trim(structfile)
380       format(' No variograms found in file ',a)
          go to 9890
        end if

! -- Memory is now allocated based on the contents of the structure file.

        allocate(structure(numstruct),vario(numvario),stat=ierr)
        if(ierr.ne.0)then
	  write(amessage,50)
	  go to 9890
	end if

! -- The remainder of the structure file is now read.

        call read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
        structure,vario)
        if(ifail.ne.0) go to 9900

! -- Next the east and north coordinate (relative to top left corner of the grid),
!    is calculated for every cell in the finite-difference grid.

        if(num_pilot_points.ne.0)then
          idim=num_pilot_points
        else
          idim=1
        end if
        allocate(eastdat(idim),northdat(idim),  &
        valdat(idim),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

! -- The coordinates of the pilot points have the grid left corner easting and
!    northing subtracted and are rotated.

        if(num_pilot_points.ne.0)then
          do i=1,num_pilot_points
            etemp=pilot_point_east(i)-gridspec%east_corner
            ntemp=pilot_point_north(i)-gridspec%north_corner
            pilot_point_east(i)=etemp*gridspec%cosang+ntemp*gridspec%sinang
            pilot_point_north(i)=ntemp*gridspec%cosang-etemp*gridspec%sinang
          end do
        end if

! Average field values are now requested.

389     write(6,*)
        icount=0
        write(6,394)
394     format(' Mean field values within each zone are now requested.')
        izone=0
395     izone=izone+1
        call num2char(intval(izone),anum1)
        if(izone.gt.numint) go to 397
        if(astructure(izone).eq.' ') go to 395
        do i=1,numstruct
          if(astructure(izone).eq.structure(i)%structname)then
            istruct=i
            go to 390
          end if
        end do
        write(amessage,402) trim(structfile),trim(astructure(izone)),trim(anum1)
        go to 9890
390     continue
        npptzone=0
        if(num_pilot_points.gt.0)then
          do j=1,num_pilot_points
            if(pilot_point_zone(j).eq.intval(izone)) go to 1210
          end do
          go to 1220
1210      npptzone=1
1220      continue
        end if
        itrans=structure(istruct)%transform
        write(6,*)
        if(npptzone.eq.0)then
391       if(itrans.eq.1)then
            write(6,392,advance='no') trim(anum1)
392         format('    Enter mean field value in zone with integer value ',a,': ')
          else
            write(6,393,advance='no') trim(anum1)
393         format('    Enter mean field value in zone with integer value ',a,': ')
          end if
          iitemp=key_read(avval(izone))
          if(escset.ne.0)then
            escset=0
          end if
          if(iitemp.ne.0) then
            write(6,310)
            go to 391
          end if
          if(itrans.eq.1)then
            if(avval(izone).le.0.0)then
              write(6,399)
399           format(/,'    Must be greater than zero due to log variogram - try again.',/)
              go to 391
            end if
            avval(izone)=log10(avval(izone))
          end if
        else
1299      if(itrans.eq.1)then
            write(6,1230) trim(anum1)
1230        format('    Enter mean field value in zone with integer value ',a)
          else
            write(6,1240) trim(anum1)
1240        format('    Enter mean field value in zone with integer value ',a)
          end if
          write(6,1250,advance='no')
1250      format('    (Hit <Enter> to obtain this value from conditioning ',  &
          'points in zone): ')
          read(5,'(a)') anum2
          if(anum2.ne.' ')then
            anum2=adjustl(anum2)
            if((anum2(1:2).eq.'e ').or.(anum2(1:2).eq.'E '))then
              escset=1
            else
              call char2num(ifail,anum2,avval(izone))
              if(ifail.ne.0) then
                write(6,310)
                go to 1299
              end if
            end if
            if(itrans.eq.1)then
              if(avval(izone).le.0.0)then
                write(6,399)
                go to 1299
              end if
              avval(izone)=log10(avval(izone))
            end if
          else
            avval(izone)=0.0
            nz=0
            do j=1,num_pilot_points
              if(pilot_point_zone(j).eq.intval(izone))then
                nz=nz+1
                if(itrans.eq.0)then
                  avval(izone)=avval(izone)+pilot_point_val(j)
                else
                  if(pilot_point_val(j).le.0.0) then
                    write(amessage,1260) trim(anum1)
1260                format(' According to the structure file, the variogram(s) ', &
                    'for zone with integer value ',a,' are log-transformed. However ', &
                    'a pilot point with a zero or negative value has been ',   &
                    'assigned to this zone.')
                    go to 9890
                  end if
                  avval(izone)=avval(izone)+log10(pilot_point_val(j))
                end if
              end if
            end do
            avval(izone)=avval(izone)/nz
          end if
        end if
        if(escset.eq.1)then
          escset=0
          if(icount.ge.1)then
            go to 389
          else
            jerr=0
            deallocate(eastdat,northdat,valdat,stat=ierr)
            if(ierr.ne.0) jerr=jerr+1
            deallocate(structure,vario,stat=ierr)
            if(ierr.ne.0) jerr=jerr+1
            if(jerr.ne.0)then
              write(amessage,57)
              go to 9890
            end if
            write(6,*)
            if(iunif.eq.0)then
              structunit=nextunit()
              open(unit=structunit,file=structfile,status='old')
              go to 1369
            else
              go to 363
            end if
          end if
        end if
        icount=icount+1
        go to 395

397     continue

        write(6,*)
1309    write(6,1315,advance='no')
1315    format(' Enter integer seed for random number generator [324853]: ')
        read(5,'(a)') anum2
        if(anum2.eq.' ') then
          iseed=324853
        else
          anum2=adjustl(anum2)
          if((anum2(1:2).eq.'E ').or.(anum2(1:2).eq.'e ')) then
            go to 389
          end if
          call char2num(ifail,anum2,iseed)
          if(ifail.ne.0)then
            write(6,1319)
1319        format(' Data input error - try again.')
            go to 1309
          end if
          if(iseed.le.0)then
            write(6,1320)
1320        format(' Must be greater than zero - try again.')
            go to 1309
          end if
        end if

        nparmcall=0
        do isim=1,irealization
          call num2char(isim,anum1)
          write(amessage,396) trim(anum1)
396       format(' Realization # ',a,' --->')
  	  call write_message(leadspace='yes')
          rarray=1.0e35          ! rarray is an array

          do izone=1,numint
            if(astructure(izone).eq.' ') go to 410
            itemp=intval(izone)
            ndat=0
            if(num_pilot_points.ne.0)then
              do i=1,num_pilot_points
                if(pilot_point_zone(i).eq.itemp)then
                  ndat=ndat+1
                  eastdat(ndat)=pilot_point_east(i)
                  northdat(ndat)=pilot_point_north(i)
                  valdat(ndat)=pilot_point_val(i)
                end if
              end do
            end if
            call num2char(intval(izone),anum1)
            write(6,401) trim(anum1)
401         format('  - carrying out field generation for integer array zone ',a,'....')
            do i=1,numstruct
              if(astructure(izone).eq.structure(i)%structname)then
                istruct=i
                go to 405
              end if
            end do
            call num2char(intval(izone),anum1)
            write(amessage,402) trim(structfile),trim(astructure(izone)),trim(anum1)
402         format(' Structure file ',a,' does not include specifications for ', &
            'structure "',a,'" needed for field generation of zone pertaining to ', &
            'integer array value of ',a)
            go to 9890
405         continue
            if(akrig(izone).eq.'o')then
              k_ktype=1
            else
              k_ktype=0
            end if
            n_nst=structure(istruct)%numvariogram
            c_c0=structure(istruct)%nugget
            itrans=structure(istruct)%transform
            do i=1,n_nst
              atemp=structure(istruct)%variogram_name(i)
              do j=1,numvario
                if(vario(j)%varname.eq.atemp)go to 420
              end do
              write(amessage,415) trim(atemp),trim(structure(istruct)%structname),  &
              trim(structfile)
415           format(' Specifications for variogram "',a,'" cited in structure "',a,  &
              '" in file ',a,' not found in this file.')
              go to 9890
420           continue
              i_it(i)=vario(j)%vartype
              if(i_it(i).eq.4)then
                write(amessage,422) trim(vario(j)%varname)
422             format(' Variogram "',a,'" uses the power model. Field generation ', &
                'cannot be carried out using a variogram of this type.')
                go to 9890
              end if
              if(i_it(i).eq.3)then
                write(amessage,423) trim(vario(j)%varname)
423             format(' Variogram "',a,'" uses the Gaussian model. Field generation ', &
                'cannot be carried out using a variogram of this type.')
                go to 9890
              end if

              c_cc(i)=structure(istruct)%variogram_contrib(i)
              a_ang(i)=vario(j)%angle
              a_aa(i)=vario(j)%a
              a_anis(i)=vario(j)%anis
            end do

            if(itrans.eq.1)then
              if(ndat.ne.0)then
                do j=1,ndat
                  if(valdat(j).le.0.0)then
                    call num2char(intval(izone),anum1)
                    write(amessage,1260) trim(anum1)
                    go to 9890
                  end if
                  valdat(j)=log10(valdat(j))
                end do
              end if
            end if

! -- Now a subgrid is designed for the use of the GSLIB stochasitic field generator.

! -- First we establish the nature of the subgrid, for example, whether it is uniform
!    or not.

            itemp=intval(izone)
            irowmin=99999999999
            icolmin=99999999999
            irowmax=0
            icolmax=0
            do irow=1,nrow
              do icol=1,ncol
                if(intarray(icol,irow).eq.itemp)then
                  if(icol.gt.icolmax)icolmax=icol
                  if(icol.lt.icolmin)icolmin=icol
                  if(irow.gt.irowmax)irowmax=irow
                  if(irow.lt.irowmin)irowmin=irow
                end if
              end do
            end do
            rtemp1=delrmin
            eps=3.0*epsilon(rtemp1)
            do icol=icolmin,icolmax
              if(abs(rtemp1-gridspec%delr(icol)).gt.eps) go to 1010
            end do
            rtemp2=delcmin
            eps=3.0*epsilon(rtemp2)
            do irow=irowmin,irowmax
              if(abs(rtemp2-gridspec%delc(irow)).gt.eps) go to 1010
            end do

! -- The subgrid is uniform. Hence the grid used by GSLIB will be a model
!    subgrid.

            uniform=1
            x_xmn=gridspec%delr(1)/2.0
            do j=2,icolmin
              x_xmn=x_xmn+(gridspec%delr(j)+gridspec%delr(j-1))/2.0
            end do
            y_ymn=-gridspec%delc(1)/2.0
            do j=2,irowmax
              y_ymn=y_ymn-(gridspec%delc(j)+gridspec%delc(j-1))/2.0
            end do
            n_nx=icolmax-icolmin+1
            n_ny=irowmax-irowmin+1
            x_xsiz=rtemp1
            y_ysiz=rtemp2

            go to 1050


1010        continue
            uniform=0
            x_xmn=0.0
            if(icolmin.gt.1)then
              do j=2,icolmin
                x_xmn=x_xmn+gridspec%delr(j-1)
              end do
            end if
            xmax=x_xmn+gridspec%delr(icolmin)
            if(icolmax.gt.icolmin)then
              do j=icolmin+1,icolmax
                xmax=xmax+gridspec%delr(j)
              end do
            end if
            y_ymn=0
            if(irowmin.gt.1)then
              do j=2,irowmin
                y_ymn=y_ymn-gridspec%delc(j-1)
              end do
            end if
            ymax=y_ymn
            y_ymn=y_ymn-gridspec%delc(irowmin)
            if(irowmax.gt.irowmin)then
              do j=irowmin+1,irowmax
                y_ymn=y_ymn-gridspec%delc(j)
              end do
            end if
            x_xsiz=delrmin
            y_ysiz=delcmin
            n_nx=(xmax-x_xmn)/delrmin+1
            n_ny=(ymax-y_ymn)/delcmin+1
!            write(55,*) 'zone ',itemp                       !debug
!            write(55,*) 'n_nx = ',n_nx                      !debug
!            write(55,*) 'n_ny = ',n_ny                      !debug
!            write(55,*) ' delrmin = ', delrmin              !debug
!            write(55,*) ' delcmin = ', delcmin              !debug
!            write(55,*) ' x_xsiz = ', x_xsiz                !debug
!            write(55,*) ' y_ysiz = ', y_ysiz                !debug
!            write(55,*) ' x_xmin, xmax = ',x_xmn,xmax      !debug
!            write(55,*) ' y_ymin, ymax = ',y_ymn,ymax      !debug

1050        continue
            nparmcall=nparmcall+1
            ncond=ndat
            n_nd=max(ncond,1)

            call readparm(nparmcall,k_ktype,radmax(izone),n_nx,x_xmn,x_xsiz, &
            n_ny,y_ymn,y_ysiz,ndmin(izone),ndmax(izone),ncnode(izone),       &
            n_nst,c_c0,i_it,c_cc,a_ang,a_aa,a_anis,                          &
            ncond,n_nd,eastdat,northdat,valdat,astdev,iseed)

            if(n_nx*n_ny.gt.maxtempdim)then
              if(allocated(temparray)) then
                deallocate(temparray,stat=ierr)
                if(ierr.ne.0)then
                  write(amessage,57)
	          go to 9890
	        end if
              end if
              maxtempdim=n_nx*n_ny
              allocate(temparray(maxtempdim),stat=ierr)
              if(ierr.ne.0)then
                write(amessage,57)
                go to 9890
              end if
            end if
            temparray=1.0e35          ! temparray is an array
            call sgsim(maxtempdim,temparray,avval(izone),astdev)

            itemp=intval(izone)
            if(uniform.eq.1)then
              ind=0
              do irow=irowmax,irowmin,-1
                do icol=icolmin,icolmax
                  ind=ind+1
                  if(intarray(icol,irow).eq.itemp)then
                    if(itrans.eq.0)then
                      rarray(icol,irow)=temparray(ind)
                    else
                      if(temparray(ind).lt.35)then
                        rarray(icol,irow)=10**(temparray(ind))
                      else
                        rarray(icol,irow)=1.0e35
                      end if
                    end if
                  end if
                end do
              end do
            else

! -- First we work out what MODFLOW cell each SGSIM cell lies in.

              if(allocated(irr))then
                if(n_nx.gt.max_nx)then
                  deallocate(irr,stat=ierr)
 	          if(ierr.ne.0) then
	            write(amessage,57)
	            go to 9890
	          end if
                  max_nx=n_nx
                  allocate(irr(max_nx),stat=ierr)
                  if(ierr.ne.0) then
                    write(amessage,57)
                    go to 9890
                  end if
                end if
                if(n_ny.gt.max_ny)then
                  deallocate(icc,stat=ierr)
 	          if(ierr.ne.0) then
	            write(amessage,57)
	            go to 9890
	          end if
                  max_ny=n_ny
                  allocate(icc(max_ny),stat=ierr)
                  if(ierr.ne.0) then
                    write(amessage,57)
                    go to 9890
                  end if
                end if
              else
                max_nx=n_nx
                max_ny=n_ny
                allocate(irr(max_nx),icc(max_ny),stat=ierr)
                if(ierr.ne.0) then
                  write(amessage,57)
                  go to 9890
                end if
              end if

! -- The following code fills the irr and icc arrays. These are the model cell
!    col and row numbers in which each stochastic cell lies.

              xu=gridspec%delr(1)
              if(icolmin.gt.1)then
                do ir=2,icolmin
                  xu=xu+gridspec%delr(ir)
                end do
              end if
              ir=icolmin
              xx=x_xmn
              irr(1)=ir
              do ix=2,n_nx
                xx=xx+x_xsiz
                if(xx.gt.xu)then
                  ir=ir+1
                  if(ir.le.ncol)then
                    xu=xu+gridspec%delr(ir)
                  else
                    xu=1e30
                    ir=-9999
                  end if
                end if
                irr(ix)=ir
              end do
              yu=0.0
              if(irowmax.gt.1)then
                do ic=2,irowmax
                  yu=yu-gridspec%delc(ic-1)
                end do
              end if
              ic=irowmax
              yy=y_ymn
              icc(1)=ic
              do iy=2,n_ny
                yy=yy+y_ysiz
                if(yy.gt.yu)then
                  ic=ic-1
                  if(ic.gt.0)then
                    yu=yu+gridspec%delc(ic)
                  else
                    yu=1e30
                    ic=-9999
                  end if
                end if
                icc(iy)=ic
              end do

! -- Now we evaluate values in model cells by averaging over SGSIM cells

              itemp=intval(izone)
              do irow=irowmin,irowmax
                ncrd=nc(irow)
                do icol=icolmin,icolmax
                  ecrd=ec(icol)
                  if(intarray(icol,irow).eq.itemp)then
                    nn=0
                    rsum=0.0
                    rmindist=1e30
                    do ind=1,n_nx*n_ny
                      iy=(ind-1)/n_nx+1
                      ix=ind-(iy-1)*n_nx
                      if((irr(ix).eq.icol).and.(icc(iy).eq.irow))then
                        nn=nn+1
                        rsum=rsum+temparray(ind)
!     write(55,*) ' ind, temparray, ix, iy, irr, icc = ',ind,temparray(ind),ix,iy,irr(ix),icc(iy)  !debug
                        if(as.eq.'c')then
                          xx=x_xmn+(ix-1)*x_xsiz
                          yy=y_ymn+(iy-1)*y_ysiz
                          dist=(xx-ecrd)*(xx-ecrd)+(yy-ncrd)*(yy-ncrd)
                          if(dist.lt.rmindist)then
                            rminval=temparray(ind)
                            rmindist=dist
                          end if
                        end if
                      end if
                    end do
                    if(nn.eq.0)then
                      write(6,*) ' Internal error - contact programmer.'
                      write(6,*) ' Row = ',irow,' Column = ',icol
                      stop
                    end if
                    if(as.eq.'a')then
                      rsum=rsum/nn
                    else
                      rsum=rminval
                    end if
                    if(itrans.eq.0)then
                      rarray(icol,irow)=rsum
                    else
                      if(rsum.lt.35)then
                        rarray(icol,irow)=10**rsum
                      else
                        rarray(icol,irow)=1.0e35
                      end if
                    end if
                  end if
                end do
              end do
            end if

410       continue
          end do

          ndig=log10(float(irealization))+1
          if(ndig.lt.2)ndig=2
          call num2char(isim,anum1)
          lt=len_trim(anum1)
          lt=ndig-lt
!          atempf=trim(realbase)//ostring(1:lt)//trim(anum1)//'.ref'
          if(af.eq.'f ')then
            atempf=trim(realbase)//trim(anum1)//'.ref'
          else
            atempf=trim(realbase)//trim(anum1)//'.reu'
          end if
          call addquote(atempf,arealfile)
          outunit=nextunit()
          if(af.eq.'f ')then
            open(unit=outunit,file=atempf)
            if(headerspec.eq.'yes') then
	      write(outunit,1150) ncol,nrow
1150          format(1x,i7,1x,i7)
 	    end if
	    do irow=1,nrow
	      write(outunit,1160) rarray(:,irow)
1160          format(7(1pe14.6))
  	    end do
	    write(amessage,1170) trim(arealfile)
1170        format('  - formatted real array written to file ',a)
          else
	    open(unit=outunit,file=atempf,form='binary',iostat=ierr)
	    if(ierr.ne.0)then
	    	  open(unit=outunit,file=atempf,form='unformatted')
	    end if
!            write(outunit) kstp,kper,pertim,totim,text,ncol,nrow,ilay
            write(outunit)
	    do irow=1,nrow
	      write(outunit) rarray(:,irow)
  	    end do
	    write(amessage,1171) trim(arealfile)
1171        format('  - unformatted real array written to file ',a)
          end if
          call write_message()
          close(unit=outunit)

        end do

        go to 9900

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
	call free_grid_mem(gridspec)
        call free_point_mem()
        deallocate(intarray,rarray,temparray,stat=ierr)
        deallocate(eastdat,northdat,valdat,stat=ierr)
        deallocate(astructure,radmax,akrig,ndmin,ndmax,ncnode,stat=ierr)
        deallocate(structure,vario,stat=ierr)
        deallocate(irr,icc,ec,nc,stat=ierr)

end program fieldgen


