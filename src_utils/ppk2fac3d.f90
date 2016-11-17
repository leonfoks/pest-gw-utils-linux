module geostat

      integer,allocatable :: nisb(:),ixsbtosr(:),iysbtosr(:),izsbtosr(:)
      real,allocatable    :: x(:),y(:),z(:),vr(:),ve(:),dh(:),tmp(:),                   &
               close(:),xa(:),ya(:),za(:),vra(:),vea(:),xdb(:),ydb(:),                  &
               zdb(:),cut(:),cdf(:)
      real*8,allocatable  :: r(:),rr(:),s(:),a(:)

      real, allocatable    :: sec3(:)                               !jd

end module



program ppk2fac3d

! -- Program PPK2FAC3D performs the first step in kriging from a set of
!    pilot points to a three-dimensional grid; it builds the factors by which
!    the value of each point contributes to the value at the grid cell.

	use defn
	use inter

	implicit none

        logical     :: lopened
        integer     :: ifail,idate,iheader,i,ncol,nrow,icol,irow,ierr,structunit, &
                       numint,currint,newint,ihuge,itemp,numstruct,numvario,outunit, &
                       nbb,j,izone,ndat,ngrid,istruct,k_ktype,n_nst,itrans
        integer     :: nlay,ilay,nncol,nnrow,izonegrid,iunit
        integer     :: icol_min,icol_max,irow_min,irow_max,ncol_range,nrow_range
        integer     :: MAXSBX,MAXSBY,MAXSBZ,MAXDIS
        integer     :: iparmcall
        integer     :: spectype
        real        :: elevgrid_min,elevgrid_max,eastgrid_min,eastgrid_max,northgrid_min,northgrid_max
        real        :: s_skmean,c_c0,pmx
        real        :: temin,temax,tnmin,tnmax,televmin,televmax
        real        :: rtemp
        real, allocatable :: thick(:),bot(:)
        double precision :: minsep,minsep2,eastdiff,northdiff,distance
        double precision :: elevdiff,eetemp,nntemp
        double precision :: dtemp
        integer     :: iloc(1)
        integer     :: i_it(MAX_STRUCT_VARIO)
        real        :: c_cc(MAX_STRUCT_VARIO),                                                      &
                       a_ang1(MAX_STRUCT_VARIO),a_ang2(MAX_STRUCT_VARIO),a_ang3(MAX_STRUCT_VARIO),  &
                       a_a_hmax(MAX_STRUCT_VARIO),a_a_hmin(MAX_STRUCT_VARIO),a_a_vert(MAX_STRUCT_VARIO)
        integer, allocatable, dimension(:)             :: icellno
        integer, allocatable, dimension(:)             :: inumdat
        integer, allocatable, dimension(:)             :: minpt,maxpt
        integer, allocatable, dimension(:,:,:)         :: intarray
        real, allocatable, dimension(:)                :: radmax_hmax,radmax_hmin,radmax_vert
        real, allocatable, dimension(:)                :: eastgrid,northgrid,elevgrid
        real,allocatable, dimension(:)                 :: eastdat,northdat,elevdat,valdat
        real,allocatable,dimension(:)                  :: east,north
        real, allocatable, dimension(:,:,:)            :: bottom,elev
        character (len=1)                              :: alogtrans,aoutfile,arealformat
        character (len=1),  allocatable, dimension(:)  :: akrig
        character (len=10)                             :: alay
        character (len=10), allocatable, dimension(:)  :: astructure
        character (len=15)                             :: anum1,atemp
	character (len=200)                            :: aprompt,structfile, &
                                                          atempf,outfile
        character (len=200)                            :: intbasename,elevbasename
        character (len=200)                            :: infile
        character (len=200), allocatable, dimension(:) :: elevfilename

	type (modelgrid) gridspec
        type (geostructure), allocatable, dimension(:) :: structure(:)
        type (variogram), allocatable, dimension(:)    :: vario(:)

	integer, parameter :: MAXINT=100
	integer		   :: intval(MAXINT),iwork(MAXINT)


        iparmcall=0

	write(amessage,5)
5	format(' Program PPK2FAC3D calculates point-to-cell factors by which kriging ',&
	'is undertaken from a set of pilot points to a 3D finite-difference grid.')
	call write_message(leadspace='yes',endspace='yes')

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  datespec=1
	  headerspec='no'
	  go to 1
!	  write(amessage,7)
!7	  format(' A settings file (settings.fig) was not found in the ', &
!	  'current directory.')
!	  call write_message
!	  go to 9900
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
1       continue

	call readfig(gridspec%specfile)
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

! -- This is a bit messy. We now re-read the grid spec file to see if it adopts the
!    new protocol where layer data is provided.

        spectype=-99999999
        iunit=nextunit()
        open(unit=iunit,file=gridspec%specfile,status='old',err=800)
        read(iunit,*,err=800,end=800) itemp,itemp,nlay
        if(nlay.le.0) go to 800
        read(iunit,*,err=800,end=800) dtemp,dtemp,dtemp
        read(iunit,*,err=800,end=800) (rtemp,icol=1,ncol)
        read(iunit,*,err=800,end=800) (rtemp,irow=1,nrow)
        read(iunit,*,err=800,end=800) spectype
        if(spectype.eq.0)then
          allocate(thick(nlay),bot(0:nlay),stat=ierr)
          if(ierr.ne.0) go to 800
          read(iunit,*,err=800,end=800) bot(0)
          read(iunit,*,err=800,end=800) (thick(ilay),ilay=1,nlay)
          do ilay=1,nlay
            bot(ilay)=bot(ilay-1)-thick(ilay)
         end do
        else if(spectype.eq.1)then
	  allocate(elevfilename(0:nlay),stat=ierr)
	  if(ierr.ne.0) go to 800
	  do ilay=0,nlay
	    read(iunit,*,err=800,end=800) elevfilename(ilay)
	  end do
	else
	  go to 800
	end if
	close(unit=iunit)
	go to 810
800     if(allocated(thick))deallocate(thick,stat=ierr)
        if(allocated(bot))deallocate(bot,stat=ierr)
        if(allocated(elevfilename)) deallocate(elevfilename,stat=ierr)
        inquire(unit=iunit,opened=lopened)
        if(lopened) close(unit=iunit,iostat=ierr)
        spectype=-99999999
810     continue

        if(spectype.ge.0) go to 850
20      write(6,22,advance='no')
22      format(' How many layers in model? ')
        if(key_read(nlay).ne.0) go to 20
	if(escset.eq.1) then
	  escset=0
          call free_grid_mem(gridspec)
	  write(6,*)
	  go to 10
        end if
        if(nlay.le.0) go to 20

        write(6,*)
570	continue
569     write(6,572,advance='no')
572     format(' Enter filename base of layer bottom elevation array files: ')
        read(5,'(a)') atempf
        atempf=adjustl(atempf)
        if((atempf(1:2).eq.'E ').or.(atempf(1:2).eq.'e '))then
          write(6,*)
          go to 20
        end if
        if(atempf.eq.' ') go to 570
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,elevbasename,1,nbb)
        if(ifail.ne.0) go to 569

850     continue

        allocate(bottom(ncol,nrow,0:nlay),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,50)
          go to 9890
        end if

	do ilay=0,nlay
	  if(spectype.eq.0)then
	    do irow=1,nrow
	      do icol=1,ncol
	        bottom(icol,irow,ilay)=bot(ilay)
	      end do
	    end do
	  else
	    if(spectype.eq.1)then
	      infile=elevfilename(ilay)
	    else
	      call num2char(ilay,alay)
	      infile=trim(elevbasename)//trim(alay)//'.ref'
	    end if
	    write(6,52) trim(infile)
	    iunit=nextunit()
	    open(unit=iunit,file=infile,status='old',iostat=ierr)
	    if(ierr.ne.0)then
	      write(amessage,460) trim(infile)
460	      format(' Cannot open real array file ',a,'.')
              go to 9890
            end if
	    if(headerspec.eq.'yes') then
	      read(iunit,'(a)',end=461) cline
	      call linesplit(ifail,2)
	      if(ifail.ne.0) go to 461
	      nncol=char2int(ifail,1)
	      if(ifail.ne.0) go to 461
	      nnrow=char2int(ifail,2)
	      if(ifail.ne.0) go to 461
	      if((nncol.ne.ncol).or.(nnrow.ne.nrow)) then
	        write(amessage,64) trim(infile)
                go to 9890
              end if
              go to 465
461           write(amessage,459) trim(infile)
459           format('Error reading NCOL/NROW header to real array contained in file ',a,'.')
              go to 9890
465           continue
            end if
            do irow=1,nrow
              read(iunit,*,err=466,end=467) (bottom(icol,irow,ilay),icol=1,ncol)
            end do
            go to 468
466         write(amessage,462) trim(infile)
462         format(' Error encountered in reading real array from file ',a,'.')
            go to 9890
467         write(amessage,464) trim(infile)
464         format(' Unexpected end encountered while reading real array from file ',a,'.')
            go to 9890
468         continue
            close(unit=iunit)
            write(6,681) trim(infile)
          end if
        end do

! -- Elevations of all cell centres are computed.

        allocate(elev(ncol,nrow,nlay),stat=ierr)
        if(ierr.ne.0) then
          write(amessage,50)
          go to 9890
        end if
        do ilay=1,nlay
          do icol=1,ncol
            do irow=1,nrow
              elev(icol,irow,ilay)=0.5*(bottom(icol,irow,ilay-1)+bottom(icol,irow,ilay))
            end do
          end do
        end do
        deallocate(bottom,stat=ierr)
        if(ierr.ne.0) then
          write(amessage,57)
          go to 9890
        end if

        write(6,*)
30      continue
        call read_3d_pilot_points_file(ifail, &
	' Enter name of 3D pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  deallocate(elev,stat=ierr)
	  if(ierr.ne.0)then
	    write(amessage,57)
57	    format(' Memory management error. Cannot continue execution.')
	    go to 9890
	  end if
	  if(spectype.eq.-99999999)then
	    go to 570
	  else
	    if(allocated(thick)) deallocate(thick,stat=ierr)
	    if(allocated(bot)) deallocate(bot,stat=ierr)
	    if(allocated(elevfilename)) deallocate(elevfilename,stat=ierr)
	    call free_grid_mem(gridspec)
	    go to 10
	  end if
	end if

!        write(6,*)
!61      write(6,63,advance='no')
!63      format(' Enter minimum allowable points separation: ')
!        if(key_read(minsep).ne.0) go to 61
!        if(escset.eq.1) then
!          write(6,*)
!          escset=0
!          go to 30
!        end if
!        imessage=0
!        if(minsep.lt.0.0d0)minsep=0.0d0
!!        if(minsep.gt.0.0d0)then
!          minsep2=minsep*minsep
!          if(num_pilot_points.gt.1)then
!          do i=1,num_pilot_points-1
!            do j=i+1,num_pilot_points
!              eastdiff=pilot_point_east(i)-pilot_point_east(j)
!              northdiff=pilot_point_north(i)-pilot_point_north(j)
!              elevdiff=pilot_point_elev(i)-pilot_point_elev(j)
!              distance=eastdiff*eastdiff+northdiff*northdiff+elevdiff*elevdiff
!!              if((distance.le.minsep2).and. &
!              if((distance.lt.minsep2).and. &
!                 (pilot_point_zone(i).eq.pilot_point_zone(j)))then
!                imessage=imessage+1
!                if(imessage.gt.20) go to 9900
!                if(imessage.eq.1)then
!                  write(amessage,65)
!65                format(' The following points are separated by less than the ', &
!                  'minimum separation ---->')
!                  call write_message(leadspace='yes')
!                end if
!                write(amessage,66) trim(pilot_point_id(i)), &
!                trim(pilot_point_id(j)),sqrt(distance)
!66              format(1x,a,t15,a, t30, '(separation = ',f12.3,')')
!                call write_message()
!              end if
!            end do
!          end do
!          end if
!!        end if
!        if(imessage.gt.0) go to 9900

        izonegrid=-99999999
        intbasename=' '
        write(6,*)
70	continue
69      write(6,72)
72      format(' Enter filename base of layer zonal integer array files.')
        write(6,73,advance='no')
73      format(' Press <Enter> if entire grid belongs to single zone: ')
        read(5,'(a)') atempf
        atempf=adjustl(atempf)
        if((atempf(1:2).eq.'E ').or.(atempf(1:2).eq.'e '))then
          write(6,*)
          go to 30
        end if
75      continue
        if(atempf.eq.' ')then
76        write(6,74,advance='no')
74        format(' Enter zone number for entire grid: ')
          if(key_read(izonegrid).ne.0) go to 76
	  if(escset.eq.1) then
	    escset=0
	    write(6,*)
            go to 70
          end if
          if(izonegrid.eq.-99999999) go to 76
          go to 78
        end if
        izonegrid=-99999999
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,intbasename,1,nbb)
        if(ifail.ne.0) then
          write(6,*)
          go to 69
        end if
	allocate(intarray(ncol,nrow,nlay),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run PPK2FAC3D.')
	  go to 9890
	end if
	do ilay=1,nlay
	  call num2char(ilay,alay)
	  infile=trim(intbasename)//trim(alay)//'.inf'
	  write(6,52) trim(infile)
52	  format(/,' - reading file ',a,'...')
	  iunit=nextunit()
	  open(unit=iunit,file=infile,status='old',iostat=ierr)
	  if(ierr.ne.0)then
	    write(amessage,60) trim(infile)
60	    format(' Cannot open integer array file ',a,'.')
            go to 9890
          end if
	  if(headerspec.eq.'yes') then
	    read(iunit,'(a)',end=618) cline
	    call linesplit(ifail,2)
	    if(ifail.ne.0) go to 618
	    nncol=char2int(ifail,1)
	    if(ifail.ne.0) go to 618
	    nnrow=char2int(ifail,2)
	    if(ifail.ne.0) go to 618
	    if((nncol.ne.ncol).or.(nnrow.ne.nrow)) then
	      write(amessage,64) trim(infile)
64	      format(' Number of columns and rows specified in NCOL/NROW header to file ',a,  &
              ' does not agree with grid dimensions in grid specification file.')
              go to 9890
            end if
            go to 651
618         write(amessage,611) trim(infile)
611         format(' Error reading NCOL/NROW header to integer array in file ',a,'.')
            go to 9890
651         continue
          end if
          do irow=1,nrow
            read(iunit,*,err=666,end=676) (intarray(icol,irow,ilay),icol=1,ncol)
          end do
          go to 68
666       write(amessage,62) trim(infile)
62        format(' Error encountered in reading integer array from file ',a,'.')
          go to 9890
676       write(amessage,671) trim(infile)
671       format(' Unexpected end encountered while reading integer array from file ',a,'.')
          go to 9890
68        continue
          close(unit=iunit)
          write(6,681) trim(infile)
681       format(' - file ',a,' read ok.')
        end do

78      continue
        write(6,*)
80	continue
	aprompt=' Enter name of structure file: '
	call open_input_file(ifail,aprompt,structfile,structunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
	  if(izonegrid.eq.-99999999)then
	    deallocate(intarray,stat=ierr)
	    if(ierr.ne.0)then
	      write(6,57)
	      go to 9890
	    end if
	    go to 70
	  else
	    atempf=' '
	    go to 75
	  end if

	end if

! -- The different integers comprising the integer array are identified.

        if(izonegrid.eq.-99999999)then
	  numint=1
	  currint=intarray(1,1,1)
	  intval(1)=currint
	  lay_travel: do ilay=1,nlay
            row_travel: do irow=1,nrow
              col_travel: do icol=1,ncol
	        newint=intarray(icol,irow,ilay)
	        if(newint.eq.currint) cycle col_travel
	        currint=newint
	        prev_int: do i=1,numint
	          if(newint.eq.intval(i)) cycle col_travel
	        end do prev_int
	        numint=numint+1
	        if(numint.gt.MAXINT) then
	          call num2char(MAXINT,anum1)
	          write(amessage,150) trim(anum1),trim(intbasename)
150	          format(' The zone integer arrays can hold a maximum of ',a, &
	          ' different integers for proper PPK2FAC3D execution. The arrays ', &
	          ' contained in files ',a,'*.inf holds more than this. To increase ',&
	          'this limit, edit PPK2FAC3D source code, increase parameter ', &
	          'MAXINT, and recompile.')
	          go to 9890
	        end if
	        intval(numint)=newint
	      end do col_travel
	    end do row_travel
	  end do lay_travel
	else
	  numint=1
	  intval(1)=izonegrid
	end if

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

        allocate(astructure(numint),    &
        akrig(numint),minpt(numint),maxpt(numint),stat=ierr)  !de-allocate below for "e"
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
	allocate(radmax_hmax(numint),radmax_hmin(numint),radmax_vert(numint),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
        radmax_hmax=0.0
        radmax_hmin=0.0
        radmax_vert=0.0
        akrig=' '
        astructure=' '

200     continue
        if(izonegrid.eq.-99999999)then
 	  write(6,*)
205	  write(6,210)
210	  format(' The following zones have been detected in the integer ',&
	  'arrays:-')
	end if
	i=0
255	i=i+1
	if(i.gt.numint) go to 350
        write(6,*)
259	call num2char(intval(i),anum1)
260     continue
        if(izonegrid.eq.-99999999)then
	  write(6,270) trim(anum1)
270	  format('    For zone characterised by integer value of ',a,':- ')
          atemp=' '
275       write(6,280,advance='no')
280       format('    Enter structure name (blank if no interpolation for this zone): ')
          read(5,'(a)') atemp
          if(index(eschar,atemp(1:2)).ne.0)then
	    if(i.gt.2) then
	      i=i-2
	      go to 255
	    else if(i.eq.2) then
	      go to 200
	    else
              deallocate(astructure,radmax_hmax,radmax_hmin,radmax_vert,akrig,minpt,maxpt,stat=ierr)
	      if(ierr.ne.0) then
	        write(amessage,57)
	        go to 9890
	      end if
              close(unit=structunit)
              write(6,*)
	      go to 80
	    end if
	  end if
	else
	  call num2char(intval(1),anum1)
282       write(6,281,advance='no') trim(anum1)
281       format(' Enter structure name for zone with integer value of ',a,': ')
          read(5,'(a)') atemp
          if(atemp.eq.' ') go to 282
          atemp=adjustl(atemp)
          if((atemp(1:2).eq.'E ').or.(atemp(1:2).eq.'e '))then
              deallocate(astructure,radmax_hmax,radmax_hmin,radmax_vert,akrig,minpt,maxpt,stat=ierr)
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
285       format(' Structure name must be 10 characters or less - try again.')
          if(izonegrid.eq.-99999999) amessage='   '//trim(amessage)
          call write_message()
          go to 260
        end if
        astructure(i)=atemp
        if(astructure(i).eq.' ') go to 255
        call casetrans(astructure(i),'lo')
289     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,290,advance='no')
290     format(' Perform simple or ordinary kriging [s/o]: ')
	read(5,'(a)') akrig(i)
	if(akrig(i).eq.' ') go to 289
	if(index(eschar,akrig(i)).ne.0) then
	  write(6,*)
	  go to 259
	end if
	call casetrans(akrig(i),'lo')
	if((akrig(i).ne.'s').and.(akrig(i).ne.'o')) go to 289
299     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,300,advance='no')
300     format(' Enter search radius in maximum horizontal elongation dirn: ')
        itemp=key_read(radmax_hmax(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 289
        end if
	if(itemp.ne.0) then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,310)
310	  format(' Data input error  - try again.')
	  go to 299
	end if
        if(radmax_hmax(i).le.0.0d0)then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,311)
311       format(' Must be greater than zero - try again.')
          go to 299
        end if
        if(radmax_hmax(i).gt.1.0e15)radmax_hmax(i)=1.0e15
312     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,313,advance='no')
313     format(' Enter search radius in minimum horizontal elongation dirn: ')
        itemp=key_read(radmax_hmin(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 299
        end if
	if(itemp.ne.0) then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,310)
	  go to 312
	end if
        if(radmax_hmin(i).le.0.0d0)then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,311)
          go to 312
        end if
        if(radmax_hmin(i).gt.1.0e15)radmax_hmin(i)=1.0e15
314     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,315,advance='no')
315     format(' Enter search radius in vertical dirn: ')
        itemp=key_read(radmax_vert(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 312
        end if
	if(itemp.ne.0) then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,310)
	  go to 314
	end if
        if(radmax_vert(i).le.0.0d0)then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,311)
          go to 314
        end if
        if(radmax_vert(i).gt.1.0e15)radmax_vert(i)=1.0e15

329     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,330,advance='no')
330     format(' Enter minimum number of points to use for interpolation: ')
        itemp=key_read(minpt(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 314
        end if
	if(itemp.ne.0) then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,310)
	  go to 329
	end if
        if(minpt(i).lt.1)then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,332)
332       format(' Must be greater than 1 - try again.')
          go to 329
        end if
335     continue
        if(izonegrid.eq.-99999999)then
          write(6,'(a3)',advance='no') '   '
        end if
        write(6,337,advance='no')
337     format(' Enter maximum number of pilot points to use for interpolation: ')
        itemp=key_read(maxpt(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 329
        end if
	if(itemp.ne.0) then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,310)
	  go to 335
	end if
        if(maxpt(i).lt.minpt(i))then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,340)
340       format(' Must be greater than min. no. of points - try again.')
          go to 335
        end if
        if(maxpt(i).gt.500)then
          if(izonegrid.eq.-99999999)then
            write(6,'(a3)',advance='no') '   '
          end if
          write(6,341)
341       format(' Must not be greater than 500 - try again.')
          go to 335
        end if
        go to 255

350	continue

        do i=1,numint
          if(astructure(i).ne.' ')go to 351
        end do
        write(amessage,352)
352     format(' Interpolation must take place for at least one zone in integer ', &
        'arrays - try again.')
        call write_message(leadspace='yes',endspace='yes')
        go to 200
351     continue

353     write(6,*)
354     write(6,355,advance='no')
355     format(' Enter name for interpolation factor file: ')
        read(5,'(a)') atempf
        atempf=adjustl(atempf)
        if(index(eschar,atempf(1:2)).ne.0) then
          i=numint-1
          go to 255
        end if
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,outfile,1,nbb)
        if(ifail.ne.0) go to 354
356     write(6,357,advance='no')
357     format(' Is this a formatted or unformatted file? [f/u]: ')
        read(5,'(a)') aoutfile
        if((aoutfile.eq.'e').or.(aoutfile.eq.'E'))then
          write(6,*)
          go to 354
        end if
        if(aoutfile.eq.' ') go to 356
        call casetrans(aoutfile,'lo')
        if((aoutfile.ne.'u').and.(aoutfile.ne.'f')) go to 356
        outunit=nextunit()
        if(aoutfile.eq.'f')then
          open(unit=outunit,file=outfile)
        else
          open(unit=outunit,file=outfile,form='binary')
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

        allocate(eastdat(num_pilot_points),northdat(num_pilot_points),elevdat(num_pilot_points),         &
        valdat(num_pilot_points),inumdat(num_pilot_points),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
	allocate(east(ncol),north(nrow),                                              &
        eastgrid(ncol*nrow*nlay),northgrid(ncol*nrow*nlay),elevgrid(ncol*nrow*nlay),            &
        icellno(ncol*nrow*nlay),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

        dtemp=gridspec%delr(1)/2.0
        east(1)=dtemp
        do i=2,ncol
          dtemp=dtemp+(gridspec%delr(i)+gridspec%delr(i-1))/2.0
          east(i)=dtemp
        end do
        dtemp=-gridspec%delc(1)/2.0
        north(1)=dtemp
        do i=2,nrow
          dtemp=dtemp-(gridspec%delc(i)+gridspec%delc(i-1))/2.0
          north(i)=dtemp
        end do

! -- The coordinates of the pilot points have the grid left corner easting and
!    northing subtracted.

! -- They are then expressed in a coordinate system where the origin is at the top
!    left corner of the grid and the row direction is presumed to be oriented easterly.

        do i=1,num_pilot_points
          eetemp=pilot_point_east(i)-gridspec%east_corner
          nntemp=pilot_point_north(i)-gridspec%north_corner
          pilot_point_east(i)=eetemp*gridspec%cosang+nntemp*gridspec%sinang
          pilot_point_north(i)=-eetemp*gridspec%sinang+nntemp*gridspec%cosang
        end do

! -- Kriging weights are now evaluated zone by zone for cell centres in
!    the finite-difference grid using GSLIB subroutines.

        call addquote(pilot_points_file,atempf)
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        if(intbasename.eq.' ')then
          atempf='nul'
        else
          call addquote(intbasename,atempf)
        end if
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        if(aoutfile.eq.'f')then
          write(outunit,*) ncol,nrow,nlay
        else
          write(outunit) ncol,nrow,nlay
        end if
        if(aoutfile.eq.'f')then
          write(outunit,*) num_pilot_points
          do i=1,num_pilot_points
            write(outunit,'(a)') trim(pilot_point_id(i))
          end do
        else
          write(outunit) num_pilot_points
          do i=1,num_pilot_points
            write(outunit) pilot_point_id(i)
          end do
        end if

        do izone=1,numint
          ngrid=0
          if(astructure(izone).eq.' ') go to 410
          itemp=intval(izone)
          elevgrid_min=1.0e35
          elevgrid_max=-1.0e35
          if(izonegrid.eq.-99999999)then
            icol_min=99999999
            icol_max=-99999999
            irow_min=99999999
            irow_max=-99999999
            do ilay=1,nlay
              do irow=1,nrow
                do icol=1,ncol
                  if(intarray(icol,irow,ilay).eq.itemp)then
                    if(icol.lt.icol_min)icol_min=icol
                    if(icol.gt.icol_max)icol_max=icol
                    if(irow.lt.irow_min)irow_min=irow
                    if(irow.gt.irow_max)irow_max=irow
                    ngrid=ngrid+1
                    eastgrid(ngrid)=east(icol)
                    northgrid(ngrid)=north(irow)
                    elevgrid(ngrid)=elev(icol,irow,ilay)
                    if(elevgrid(ngrid).gt.elevgrid_max)elevgrid_max=elevgrid(ngrid)
                    if(elevgrid(ngrid).lt.elevgrid_min)elevgrid_min=elevgrid(ngrid)
                    call rc2cell3d(icellno(ngrid),irow,icol,ilay,gridspec)
                  end if
                end do
              end do
            end do
          else
            icol_min=1
            icol_max=ncol
            irow_min=1
            irow_max=nrow
            do ilay=1,nlay
              do irow=1,nrow
                do icol=1,ncol
                  ngrid=ngrid+1
                  eastgrid(ngrid)=east(icol)
                  northgrid(ngrid)=north(irow)
                  elevgrid(ngrid)=elev(icol,irow,ilay)
                  if(elevgrid(ngrid).gt.elevgrid_max)elevgrid_max=elevgrid(ngrid)
                  if(elevgrid(ngrid).lt.elevgrid_min)elevgrid_min=elevgrid(ngrid)
                  call rc2cell3d(icellno(ngrid),irow,icol,ilay,gridspec)
                end do
              end do
            end do
          end if
          eastgrid_min=east(icol_min)
          eastgrid_max=east(icol_max)
          northgrid_min=north(irow_max)
          northgrid_max=north(irow_min)

          if(ngrid.eq.0) go to 410
          ndat=0
          do i=1,num_pilot_points
            if(pilot_point_zone(i).eq.itemp)then
              ndat=ndat+1
              eastdat(ndat)=pilot_point_east(i)
              northdat(ndat)=pilot_point_north(i)
              elevdat(ndat)=pilot_point_elev(i)
              valdat(ndat)=pilot_point_val(i)
              inumdat(ndat)=i
            end if
          end do
          if(ndat.eq.0) then
            call num2char(intval(izone),anum1)
            if(izonegrid.eq.-99999999)then
              write(6,399) trim(anum1)
399           format(/,' Warning: no pilot points assigned to integer array zone ',a)
              go to 410
            else
              write(amessage,398) trim(anum1)
398           format(' No points in pilot points file are assigned to integer array zone ',a,'.')
              go to 9890
            end if
          end if
          call num2char(intval(izone),anum1)
          write(6,401) trim(anum1)
401       format(/,' Carrying out interpolation for integer array zone ',a,'....')
          do i=1,numstruct
            if(astructure(izone).eq.structure(i)%structname)then
              istruct=i
              go to 405
            end if
          end do
          call num2char(intval(izone),anum1)
          write(amessage,402) trim(structfile),trim(astructure(izone)),trim(anum1)
402       format(' Structure file ',a,' does not include specifications for ', &
          'structure "',a,'" needed for interpolation of zone pertaining to ', &
          'an integer array value of ',a)
          go to 9890
405       continue
          if(akrig(izone).eq.'o')then
            k_ktype=1
            s_skmean=0.0
          else
            k_ktype=0
            if(structure(istruct)%mean.lt.-1.0e35)then
              write(amessage,407) trim(astructure(izone)),trim(structfile)
407           format(' In structure "',a,'" cited in structure file ',a,  &
              ', no mean value is provided; however a mean value is required if', &
              ' simple kriging (rather than ordinary kriging) is to be performed.')
              go to 9890
            else
              s_skmean=structure(istruct)%mean
            end if
          end if
          if(structure(istruct)%numcount_3d.ne.structure(istruct)%numvariogram)then
            write(amessage,408) trim(astructure(izone)),trim(structfile)
408         format(' Structure "',a,'" cited in structure file ',a,  &
            ' is required for 3D kriging; however at least one variogram which it cites ',  &
            'is not provided with three-dimensional specifications.')
            go to 9890
          end if
          n_nst=structure(istruct)%numvariogram
          c_c0=structure(istruct)%nugget
          pmx=structure(istruct)%maxpowercov
          itrans=structure(istruct)%transform
          do i=1,n_nst
            atemp=structure(istruct)%variogram_name(i)
            do j=1,numvario
              if(vario(j)%varname.eq.atemp)go to 420
            end do
            write(amessage,415) trim(atemp),trim(structure(istruct)%structname),  &
            trim(structfile)
415         format(' Specifications for variogram "',a,'" cited in structure "',a,  &
            '" in file ',a,' not found in this file.')
            go to 9890
420         continue
            i_it(i)=vario(j)%vartype
            if((k_ktype.eq.0).and.(i_it(i).eq.4))then
              write(amessage,422) trim(vario(j)%varname)
422           format(' Variogram "',a,'" uses the power model. Simple kriging ', &
              'cannot be carried out using a variogram of this type; use ordinary ', &
              'kriging.')
              go to 9890
            end if
            c_cc(i)=structure(istruct)%variogram_contrib(i)
            a_ang1(i)=vario(j)%ang1
            a_ang2(i)=vario(j)%ang2
            a_ang3(i)=vario(j)%ang3
            a_a_hmax(i)=vario(j)%a_hmax
            a_a_hmin(i)=vario(j)%a_hmin
            a_a_vert(i)=vario(j)%a_vert
          end do

! -- At this stage we adjust the maximum and minimum values to include pilot point values which may
!    be out of range of the grid.

          temin=1.0e35
          temax=-1.0e35
          tnmin=1.0e35
          tnmax=-1.0e35
          televmin=1.0e35
          televmax=-1.0e35
          do i=1,ndat
            if(eastdat(i).lt.temin)temin=eastdat(i)
            if(eastdat(i).gt.temax)temax=eastdat(i)
            if(northdat(i).lt.tnmin)tnmin=northdat(i)
            if(northdat(i).gt.tnmax)tnmax=northdat(i)
            if(elevdat(i).lt.televmin)televmin=elevdat(i)
            if(elevdat(i).gt.televmax)televmax=elevdat(i)
          end do

          eastgrid_min=min(temin,eastgrid_min)
          eastgrid_max=max(temax,eastgrid_max)
          northgrid_min=min(tnmin,northgrid_min)
          northgrid_max=max(tnmax,northgrid_max)
          elevgrid_min=min(televmin,elevgrid_min)
          elevgrid_max=max(televmax,elevgrid_max)

          ncol_range=icol_max-icol_min
          nrow_range=irow_max-irow_min
          iparmcall=iparmcall+1
          call readparm1(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,minpt(izone),maxpt(izone),   &
          radmax_hmax(izone),radmax_hmin(izone),radmax_vert(izone),k_ktype,s_skmean,       &
          n_nst,c_c0,i_it,c_cc,a_ang1,a_ang2,a_ang3,a_a_hmax,a_a_hmin,a_a_vert,            &
          ndat,eastdat,northdat,elevdat,valdat,ncol_range,nrow_range,iparmcall)
          call kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,                                &
          ngrid,ndat,inumdat,icellno,eastgrid,northgrid,elevgrid,                &
          outfile,aoutfile,outunit,pmx,itrans,                                   &
          ncol_range,nrow_range,eastgrid_min,eastgrid_max,northgrid_min,         &
          northgrid_max,elevgrid_min,elevgrid_max)

410     continue
        end do

        write(6,*)
        write(6,450) trim(outfile)
450     format('  - kriging factors written to file ',a)

        go to 9900


9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
	call free_grid_mem(gridspec)
        call free_point_mem()
        deallocate(intarray,east,north,eastgrid,northgrid,elevgrid,icellno,stat=ierr)
        deallocate(eastdat,northdat,elevdat,valdat,inumdat,stat=ierr)
        deallocate(astructure,radmax_hmax,radmax_hmin,radmax_vert,akrig,minpt,maxpt,stat=ierr)
        deallocate(structure,vario,stat=ierr)
        deallocate(bottom,elev,stat=ierr)
        deallocate(thick,bot,elevfilename,stat=ierr)

        call de_allocate_geostat()

end program ppk2fac3d



subroutine de_allocate_geostat

      use geostat

      if(allocated(nisb))deallocate(nisb)
      if(allocated(ixsbtosr))deallocate(ixsbtosr)
      if(allocated(iysbtosr))deallocate(iysbtosr)
      if(allocated(izsbtosr))deallocate(izsbtosr)
      if(allocated(x))deallocate(x)
      if(allocated(y))deallocate(y)
      if(allocated(z))deallocate(z)
      if(allocated(vr))deallocate(vr)
      if(allocated(ve))deallocate(ve)
      if(allocated(dh))deallocate(dh)
      if(allocated(tmp))deallocate(tmp)
      if(allocated(close))deallocate(close)
      if(allocated(xa))deallocate(xa)
      if(allocated(ya))deallocate(ya)
      if(allocated(za))deallocate(za)
      if(allocated(vra))deallocate(vra)
      if(allocated(vea))deallocate(vea)
      if(allocated(xdb))deallocate(xdb)
      if(allocated(ydb))deallocate(ydb)
      if(allocated(zdb))deallocate(zdb)
      if(allocated(cut))deallocate(cut)
      if(allocated(cdf))deallocate(cdf)
      if(allocated(r))deallocate(r)
      if(allocated(rr))deallocate(rr)
      if(allocated(s))deallocate(s)
      if(allocated(a))deallocate(a)
      if(allocated(sec3))deallocate(sec3)

end subroutine de_allocate_geostat




      subroutine readparm1(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,n_ndmin,n_ndmax,     &
      r_radius1,r_radius2,r_radius3,k_ktype,s_skmean,n_nst,c_c0,           &
      i_it,c_cc,a_ang1,a_ang2,a_ang3,a_hmax,a_hmin,a_vert,                 &
      n_nd,eastdat,northdat,elevdat,valdat,ncol_range,nrow_range,icall)

!      use       msflib
      use       geostat

      implicit integer(i-n), real(a-h, o-z)                    !jd

      include  'kt3d.inc'
      parameter(MV=100)
      real      var(MV)
      character datafl*512,jackfl*512,extfl*512,outfl*512,dbgfl*512,str*512,title*80
      logical   testfl

! -- Dimension statements added by myself for subroutine arguments.

      integer ncol_range,nrow_range
      integer i_it(n_nst)
      integer icall
      real c_c0
      real c_cc(n_nst),a_ang1(n_nst),a_ang2(n_nst),a_ang3(n_nst)
      real a_hmax(n_nst),a_hmin(n_nst),a_vert(n_nst)
      real valdat(n_nd),eastdat(n_nd),northdat(n_nd),elevdat(n_nd)

!
! FORTRAN Units:
!
      lin   = 81                       !jd
      ldbg  = 83                       !jd
      lout  = 84                       !jd
      lext  = 87                       !jd
      ljack = 88                       !jd

!
! Read Input Parameters:
!
!jd      read(lin,'(a512)',err=98) datafl
!jd      call chknam(datafl,512)
!jd      write(*,*) ' data file = ',datafl(1:40)

!jd      read(lin,*,err=98) idhl,ixl,iyl,izl,ivrl,iextv
!jd      write(*,*) ' columns = ',idhl,ixl,iyl,izl,ivrl,iextv

!jd      read(lin,*,err=98) tmin,tmax
!jd      write(*,*) ' trimming limits = ',tmin,tmax

!jd      read(lin,*,err=98) koption
!jd      write(*,*) ' kriging option = ',koption

      datafl=' '                 !jd
      idhl=-9999                 !jd
      ixl=-9999                  !jd
      iyl=-9999                  !jd
      izl=-9999                  !jd
      ivrl=-9999                 !jd
      iextv=-9999                !jd
      tmax=0.0                   !jd
      tmin=0.0                   !jd
      koption=0                  !jd
!
! This is an undocumented feature to have kt3d construct an IK-type
! distribution:
!
      iktype = 0
!jd      if(koption.lt.0) then
!jd            iktype  = 1
!jd            koption = -koption
!jd      end if
!jd      if(iktype.eq.1) then

!jd            read(lin,*,err=98) ncut
!jd            write(*,*) ' number of cutoffs = ',ncut
!
! Find the needed parameter:
!
!jd            MAXCUT = ncut
!
! Allocate the needed memory:
!21
!jd            allocate(cut(MAXCUT),stat = test)
!jd                  if(test.ne.0)then
!jd                        write(*,*)'ERROR: Allocation failed due to',
!jd     +                        ' insufficient memory.'
!jd                        stop
!jd                  end if
!22
!jd            allocate(cdf(MAXCUT),stat = test)
!jd                  if(test.ne.0)then
!jd                        write(*,*)'ERROR: Allocation failed due to',
!jd     +                        ' insufficient memory.'
!jd                        stop
!jd                  end if
!
!jd            read(lin,*,err=98) (cut(i),i=1,ncut)
!jd            write(*,*) ' cutoffs = ',(cut(i),i=1,ncut)

!jd      end if

!jd      read(lin,'(a512)',err=98) jackfl
!jd      call chknam(jackfl,512)
!jd      write(*,*) ' jackknife data file = ',jackfl(1:40)

         jackfl=' '                     !jd

!jd      read(lin,*,err=98) ixlj,iylj,izlj,ivrlj,iextvj
!jd      write(*,*) ' columns = ',ixlj,iylj,izlj,ivrlj,iextvj

       ixlj=0           !jd
       iylj=0           !jd
       izlj=0           !jd
       ivrlj=0          !jd
       iextvj=0         !jd

!jd      read(lin,*,err=98) idbg
!jd      write(*,*) ' debugging level = ',idbg

      idbg=0                                            !jd

!jd      read(lin,'(a512)',err=98) dbgfl
!jd      call chknam(dbgfl,512)
!jd      write(*,*) ' debugging file = ',dbgfl(1:40)

      dbgfl='debug_kt3d.dat'                            !jd

!jd      read(lin,'(a512)',err=98) outfl
!jd      call chknam(outfl,512)
!jd      write(*,*) ' output file = ',outfl(1:40)


      outfl='output_kt3d.dat'                           !jd

!jd      read(lin,*,err=98) nx,xmn,xsiz
!jd      write(*,*) ' nx, xmn, xsiz = ',nx,xmn,xsiz

      nx=1                                             !jd
      xmn=0.0                                          !jd
      xsiz=1.0                                         !jd

!jd      read(lin,*,err=98) ny,ymn,ysiz
!jd      write(*,*) ' ny, ymn, ysiz = ',ny,ymn,ysiz

      ny=1                                             !jd
      ymn=0.0                                          !jd
      ysiz=1.0                                         !jd

!jd      read(lin,*,err=98) nz,zmn,zsiz
!jd      write(*,*) ' nz, zmn, zsiz = ',nz,zmn,zsiz

      nz=1                                             !jd
      zmn=0.0                                          !jd
      zsiz=1.0                                         !jd

!jd      read(lin,*,err=98) nxdis,nydis,nzdis
!jd      write(*,*) ' block discretization:',nxdis,nydis,nzdis

      nxdis=1                                          !jd
      nydis=1                                          !jd
      nzdis=1                                          !jd

!jd      read(lin,*,err=98) ndmin,ndmax
!jd      write(*,*) ' ndmin,ndmax = ',ndmin,ndmax

      ndmin=n_ndmin                                    !jd
      ndmax=n_ndmax                                    !jd

!jd      read(lin,*,err=98) noct
!jd      write(*,*) ' max per octant = ',noct

      noct=0                                           !jd

!jd      read(lin,*,err=98) radius,radius1,radius2
!jd      write(*,*) ' search radii = ',radius,radius1,radius2
!jd      if(radius.lt.EPSLON) stop 'radius must be greater than zero'

      radius=r_radius1                                !jd
      radius1=r_radius2                                !jd
      radius2=r_radius3                                !jd

      radsqd = radius  * radius
      sanis1 = radius1 / radius
      sanis2 = radius2 / radius

!jd      read(lin,*,err=98) sang1,sang2,sang3
!jd      write(*,*) ' search anisotropy angles = ',sang1,sang2,sang3

!jd      read(lin,*,err=98) ktype,skmean
!jd      write(*,*) ' ktype, skmean =',ktype,skmean

      ktype=k_ktype                                    !jd
      skmean=s_skmean                                  !jd

!jd      read(lin,*,err=98) (idrif(i),i=1,9)
!jd      write(*,*) ' drift terms = ',(idrif(i),i=1,9)

      do i=1,9                                         !jd
        idrif(i)=0                                     !jd
      end do                                           !jd

!jd      read(lin,*,err=98) itrend
!jd      write(*,*) ' itrend = ',itrend

      itrend=0                                         !jd

!jd      read(lin,'(a512)',err=98) extfl
!jd      call chknam(extfl,40)
!jd      write(*,*) ' external drift file = ',extfl(1:40)

      extfl=' '                                        !jd

!jd      read(lin,*,err=98) iextve
!jd      write(*,*) ' variable in external drift file = ',iextve

      iextve=0                                         !jd

!jd      read(lin,*,err=98) nst(1),c0(1)
!jd      write(*,*) ' nst, c0 = ',nst(1),c0(1)

      nst(1)=n_nst                                     !jd
      c0(1)=c_c0                                       !jd

!jd      if(nst(1).le.0) then
!jd            write(*,9997) nst(1)
!jd 9997       format(' nst must be at least 1, it has been set to ',i4,/,
!jd     +             ' The c or a values can be set to zero')
!jd            stop
!jd      endif

      do i=1,nst(1)
!jd            read(lin,*,err=98) it(i),cc(i),ang1(i),ang2(i),ang3(i)

            it(i)=i_it(i)                               !jd
            cc(i)=c_cc(i)                               !jd
            ang1(i)=a_ang1(i)                           !jd
            ang2(i)=a_ang2(i)                           !jd
            ang3(i)=a_ang3(i)                           !jd

!jd            read(lin,*,err=98) aa(i),aa1,aa2

            aa(i)=a_hmax(i)                             !jd
            aa1  =a_hmin(i)                             !jd
            aa2  =a_vert(i)                             !jd

            anis1(i) = aa1 / max(aa(i),EPSLON)
            anis2(i) = aa2 / max(aa(i),EPSLON)
!jd            write(*,*) ' it,cc,ang[1,2,3]; ',it(i),cc(i),
!jd     +                   ang1(i),ang2(i),ang3(i)
!jd            write(*,*) ' a1 a2 a3: ',aa(i),aa1,aa2
            if(it(i).eq.4) then
                  if(aa(i).lt.0.0) stop ' INVALID power variogram'
                  if(aa(i).gt.2.0) stop ' INVALID power variogram'
            end if
      end do

      if(nst(1).eq.1)then                            !jd
        itemp=1                                      !jd
      else                                           !jd
        rtemp=-1.035                                 !jd
        do i=1,nst(1)                                !jd
          if(cc(i).gt.rtemp)then                     !jd
            rtemp=cc(i)                              !jd
            itemp=i                                  !jd
          end if                                     !jd
        end do                                       !jd
      end if                                         !jd
      sang1=ang1(itemp)                              !jd
      sang2=ang2(itemp)                              !jd
      sang3=ang3(itemp)                              !jd

!jd      close(lin)
!
! Find the needed parameters:
!
      MAXDIS = nxdis*nydis*nzdis
      MAXSAM = ndmax + 1
      MAXEQ = MAXSAM + MAXDT + 2
      MAXSBX = 1
!jd      if(nx.gt.1)then
!jd            MAXSBX = int(nx/2.00)
!jd            if(MAXSBX.gt.50)MAXSBX=50
!jd      end if

      if(ncol_range.gt.1)then                         !jd
            MAXSBX = int(ncol_range/2.00)             !jd
            if(MAXSBX.gt.50)MAXSBX=50                 !jd
      end if                                          !jd
!
      MAXSBY = 1
!jd      if(ny.gt.1)then
!jd            MAXSBY = int(ny/2.00)
!jd            if(MAXSBY.gt.50)MAXSBY=50
!jd      end if

      if(nrow_range.gt.1)then                         !jd
            MAXSBY = int(nrow_range/2.00)             !jd
            if(MAXSBY.gt.50)MAXSBY=50                 !jd
      end if                                          !jd
!
      MAXSBZ = 1
!jd      if(nz.gt.1)then
!jd            MAXSBZ = int(nz/2.00)
!jd            if(MAXSBZ.gt.50)MAXSBZ=50
!jd      end if
!
      MAXSB = MAXSBX*MAXSBY*MAXSBZ
      MXSXY = 4 * MAXSBX * MAXSBY
      MXSX  = 2 * MAXSBX
!
! Allocate the needed memory:
!1
      if(icall.eq.1)then
      allocate(nisb(MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!2
      if(icall.eq.1)then
      allocate(ixsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!3
      if(icall.eq.1)then
      allocate(iysbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!4
      if(icall.eq.1)then
      allocate(izsbtosr(8 * MAXSB),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!13
      if(icall.eq.1)then
      allocate(xa(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!14
      if(icall.eq.1)then
      allocate(ya(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!15
      if(icall.eq.1)then
      allocate(za(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!16
      if(icall.eq.1)then
      allocate(vra(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!17
      if(icall.eq.1)then
      allocate(vea(MAXSAM),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!18
      if(icall.eq.1)then
      allocate(xdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!19
      if(icall.eq.1)then
      allocate(ydb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!20
      if(icall.eq.1)then
      allocate(zdb(MAXDIS),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!23
      if(icall.eq.1)then
      allocate(r(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!24
      if(icall.eq.1)then
      allocate(rr(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!25

      if(icall.eq.1)then
      allocate(s(MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!26
      if(icall.eq.1)then
      allocate(a(MAXEQ * MAXEQ),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!
! Perform some quick error checking:
!
      if(ndmax.gt.MAXSAM) stop 'ndmax is too big - modify .inc file'
      if(ktype.eq.3.and.iextv.le.0) stop 'must have external variable'
      if(ixl.le.0.and.nx.gt.1) write(*,*) ' WARNING: ixl=0 and nx>1 ! '
      if(iyl.le.0.and.ny.gt.1) write(*,*) ' WARNING: iyl=0 and ny>1 ! '
      if(izl.le.0.and.nz.gt.1) write(*,*) ' WARNING: izl=0 and nz>1 ! '
!
! Check to make sure the data file exists, then either read in the
! data or write an error message and stop:
!

!jd      inquire(file=datafl,exist=testfl)
!jd      if(.not.testfl) then
!jd            write(*,*) 'ERROR data file ',datafl,' does not exist!'
!jd            stop
!jd      endif
!
! The data file exists so open the file and read in the header
! information. Initialize the storage that will be used to summarize
! the data found in the file:
!
      title(1:22) = 'KT3D ESTIMATES WITH: '
!jd      open(lin,file=datafl,status='OLD')
!jd      read(lin,*)
!jd      read(lin,*,err=99)       nvari
!jd      do i=1,nvari
!jd            read(lin,*)
!jd      end do
!jd      MAXDAT = 0
!jd 22   read(lin,*,end=33,err=99) (var(j),j=1,nvari)
!jd      if(var(ivrl).lt.tmin.or.var(ivrl).ge.tmax) go to 22
!jd      MAXDAT = MAXDAT + 1
!jd      go to 22
!jd 33   continue
!
! Allocate the needed memory:
!5

      MAXDAT=n_nd+1                                !jd

      if(icall.eq.1)then
      allocate(x(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!6
      if(icall.eq.1)then
      allocate(y(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!7
      if(icall.eq.1)then
      allocate(z(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if

      if(icall.eq.1)then
      allocate(sec3(MAXDAT),stat=ierr)                                                    !jd
            if(test.ne.0)then                                                             !jd
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'        !jd
                  stop                                                                    !jd
            end if                                                                        !jd
      end if
!8
      if(icall.eq.1)then
      allocate(vr(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!9

      if(icall.eq.1)then
      allocate(ve(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!10
      if(icall.eq.1)then
      allocate(dh(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!11
      if(icall.eq.1)then
      allocate(tmp(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!12
      if(icall.eq.1)then
      allocate(close(MAXDAT),stat = test)
            if(test.ne.0)then
                  write(*,*)'ERROR: Allocation failed due to insufficient memory.'
                  stop
            end if
      end if
!
!jd      rewind(lin)
!jd      read(lin,'(a58)') title(23:80)
!jd      read(lin,*,err=99)       nvari
      nd = 0
      av = 0.0
      ss = 0.0
!jd      do i=1,nvari
!jd            read(lin,'(a40)',err=99) str
!jd      end do
!
! Some tests on column numbers:
!
!jd      if(ixl.gt.nvari.or.iyl.gt.nvari.or.izl.gt.nvari.or.ivrl.gt.nvari)
!jd     +      then
!jd            write(*,*) 'There are only ',nvari,' columns in input data'
!jd            write(*,*) '  your specification is out of range'
!jd            stop
!jd      end if
!
! Read all the data until the end of the file:
!
!jd 2    read(lin,*,end=3,err=99) (var(j),j=1,nvari)
2     continue
!jd      vrt = var(ivrl)
!jd      if(vrt.lt.tmin.or.vrt.ge.tmax) go to 2
      nd = nd + 1
      if(nd.gt.n_nd) go to 3                         !jd
!jd      if(nd.gt.MAXDAT) then
!jd            write(*,*) ' ERROR: Exceeded available memory for data'
!jd            stop
!jd      end if
!
! Establish the location of this datum:
!
!jd      if(idhl.le.0) then
!jd            dh(nd) = -99
!jd      else
!jd            dh(nd) = var(idhl)
!jd      endif
!jd      if(ixl.le.0) then
!jd            x(nd) = xmn
!jd      else
!jd            x(nd) = var(ixl)
!jd      endif
!jd      if(iyl.le.0) then
!jd            y(nd) = ymn
!jd      else
!jd            y(nd) = var(iyl)
!jd      endif
!jd      if(izl.le.0) then
!jd            z(nd) = zmn
!jd      else
!jd            z(nd) = var(izl)
!jd      endif

         vrt=valdat(nd)                                 !jd
         x(nd)=eastdat(nd)                              !jd
         y(nd)=northdat(nd)                             !jd
         z(nd)=elevdat(nd)                              !jd
         dh(nd)=-99                                     !jd
         sec3(nd)=real(nd)                              !jd
!
! Establish the external drift variable (if needed):
!
      ve(nd) = 1.0
!jd      if(ktype.eq.3.or.ktype.eq.2) then
!jd            ve(nd) = var(iextv)
!jd            if(ve(nd).lt.tmin.or.ve(nd).ge.tmax) then
!jd                  write(*,*) ' External drift variable must be present at all data locations!'
!jd                  write(*,*) ' Encountered at data number ',nd
!jd                  stop
!jd            end if
!jd      end if
      vr(nd) = vrt
      av     = av + vrt
      ss     = ss + vrt*vrt
      go to 2
!jd 3    close(lin)
3     continue                     !jd
      nd=n_nd
!
! Compute the averages and variances as an error check for the user:
!
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
!jd      write(*,*) 'Data for KT3D: Variable number ',ivrl
!jd      write(*,*) '  Number   = ',nd
!jd      write(*,*) '  Average  = ',av
!jd      write(*,*) '  Variance = ',ss

       write(6,*)
       write(6,900) nd
900    format('   Number of pilot points for this zone     = ',i6)
       write(6,901) av
901    format('   Mean data value for these pilot points   = ',1pg12.5)
       write(6,902) sqrt(max(ss,0.0))
902    format('   Data standard deviation for these points = ',1pg12.5)
       write(6,903)
903    format('   Working....')

      if(nd.lt.1) then
            write(*,*) ' ERROR: there are no data'
            stop
      end if
!
! Open the debugging and output files:
!
!jd      open(ldbg,file=dbgfl,status='UNKNOWN')
!jd      open(lout,file=outfl,status='UNKNOWN')
!jd      write(lout,'(a80)') title

!jd      if(iktype.eq.0.and.koption.eq.0) then
!jd           write(lout,201) 2,nx,ny,nz
!jd           write(lout,102)
!jd 102       format('Estimate',/,'EstimationVariance')
!jd      end if
!jd      if(iktype.eq.0.and.koption.ge.1) then
!jd           write(lout,201) 7
!jd           write(lout,103)
!jd 103       format('X',/,'Y',/,'Z',/,'True',/,'Estimate',/,'EstimationVariance',/,'Error: est-true')
!jd      end if
!jd 201  format(4(1x,i4))

!jd      if(iktype.eq.1) then
!jd            if(koption.eq.0) then
!jd                  write(lout,201) ncut,nx,ny,nz
!jd            else
!jd                  write(lout,201) ncut+1
!jd            end if
!jd            do i=1,ncut
!jd                  write(lout,104) i,cut(i)
!jd 104              format('Threshold: ',i2,' = ',f12.5)
!jd            end do
!jd            if(koption.eq.1) write(lout,105)
!jd 105        format('true value')
!jd      end if
!
! Open the external drift file if needed and position it at the
! first grid node in the file:
!
!jd      if((ktype.eq.2.or.ktype.eq.3).and.koption.eq.0) then
!jd            inquire(file=extfl,exist=testfl)
!jd            if(.not.testfl) then
!jd                  write(*,*) 'ERROR file ',extfl,' does not exist!'
!jd                  stop
!jd            endif
!jd            open(lext,file=extfl,status='UNKNOWN')
!jd            read(lext,'(a40)',err=97) str
!jd            read(lext,*,err=97)       nvari
!jd            do i=1,nvari
!jd                  read(lext,'(a40)',err=97) str
!jd            end do
!jd            if(idbg.ge.3) write(ldbg,100) iextve
!jd 100        format('A secondary variable is being used.  The gridded '     &
!jd                   'file',/,'must have the same grid specifications '      &
!jd                   'as the grid you are kriging.',/,'The external '        &
!jd                   'drift variable was taken from column ',i2)
!jd      endif
!
! Set up for cross validation:
!
!jd      if(koption.eq.1) then
!jd            jackfl = datafl
!jd            idhlj  = idhl
!jd            ixlj   = ixl
!jd            iylj   = iyl
!jd            izlj   = izl
!jd            ivrlj  = ivrl
!jd            iextvj = iextv
!jd      end if
!
! Open the file with the jackknife data?
!
!jd      if(koption.gt.0) then
!jd            inquire(file=jackfl,exist=testfl)
!jd            if(.not.testfl) then
!jd                  write(*,*) 'ERROR file ',jackfl,' does not exist!'
!jd                  stop
!jd            endif
!jd            open(ljack,file=jackfl,status='OLD')
!jd            read(ljack,*,err=96)
!jd            read(ljack,*,err=96) nvarij
!jd            do i=1,nvarij
!jd                  read(ljack,*,err=96)
!jd            end do
!jd      end if
!
! Finished here:
!
      return
!
! Error in an Input File Somewhere:
!
 96   stop 'ERROR in jackknife file!'
 97   stop 'ERROR in external drift file!'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'
      end



      subroutine kt3d(MAXDIS,MAXSBX,MAXSBY,MAXSBZ,                                       &
                      n_npts,n_ndat,inumdat,icellno,epoint,npoint,zpoint,                &
                      outfile,aoutfile,outunit,pmx,itrans,                               &
                      ncol_range,nrow_range,e_min,e_max,n_min,n_max,elev_min,elev_max)
!-----------------------------------------------------------------------
!
!                Krige a 3-D Grid of Rectangular Blocks
!                **************************************
!
! This subroutine estimates point or block values of one variable by
! simple, ordinary, or kriging with a trend model.  It is also possible
! to estimate the trend directly.
!
!
!
! PROGRAM NOTES:
!
!   1. The data and parameters are passed in common blocks defined
!      in kt3d.inc.  Local storage is allocated in the subroutine
!      for kriging matrices, i.e.,
!         - xa,ya,za,vra   arrays for data within search neighborhood
!         - a,r,rr,s       kriging arrays
!         - xdb,ydb,zdb    relative position of discretization points
!         - cbb            block covariance
!   2. The kriged value and the kriging variance is written to Fortran
!      unit number "lout".
!
!
!
!
! Original:  A.G. Journel and C. Lemmer                             1981
! Revisions: A.G. Journel and C. Kostov                             1984
!-----------------------------------------------------------------------

      use geostat

      implicit integer(i-n), real(a-h, o-z)                          !jd

      include   'kt3d.inc'
      real*8     cbb
      real       var(20)
      logical    first,fircon,accept

      integer    n_npts,n_ndat,outunit,itrans                        !jd
      integer    inumdat(n_ndat),icellno(n_npts)                     !jd
      integer    ncol_range,nrow_range                               !jd
      real       epoint(n_npts),npoint(n_npts),zpoint(n_npts)        !jd
      real       e_min,e_max,n_min,n_max,elev_min,elev_max           !jd
      real       pmx                                                 !jd
      character*(*) outfile,aoutfile                                 !jd

      data       fircon/.true./
!
! Set up the rotation/anisotropy matrices that are needed for the
! variogram and search.  Also compute the maximum covariance for
! the rescaling factor:
!
!jd      write(*,*) 'Setting up rotation matrices for variogram and search'
      radsqd = radius * radius
!jd      PMX    = 999.0
      covmax = c0(1)
      do is=1,nst(1)
            call setrot_new(ang1(is),ang2(is),ang3(is),anis1(is),anis2(is),is,MAXROT,rotmat)
            if(it(is).eq.4) then
                  covmax = covmax + PMX
            else
                  covmax = covmax + cc(is)
            endif
      end do
      isrot = MAXNST + 1
      call setrot_new(sang1,sang2,sang3,sanis1,sanis2,isrot,MAXROT,rotmat)
!
! Finish computing the rescaling factor and stop if unacceptable:
!
      if(radsqd.lt.1.0) then
            resc = 2.0 * radius / max(covmax,0.0001)
      else
            resc =(4.0 * radsqd)/ max(covmax,0.0001)
      endif
      if(resc.le.0.0) then
            write(*,*) 'ERROR KT3D: The rescaling value is wrong ',resc
            write(*,*) '            Maximum covariance: ',covmax
            write(*,*) '            search radius:      ',radius
            stop
      endif
      resc = 1.0 / resc

!
! Set up for super block searching:
!
!jd      write(*,*) 'Setting up super block search strategy'
!jd      nsec = 2
         nsec = 3                                                               !jd
!jd      call setsupr(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,
!jd     +             vr,tmp,nsec,ve,dh,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb,
!jd     +             nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,
!jd     +             zmnsup,zsizsup)

      call setsupr_jd(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,  &          !jd
                   vr,tmp,nsec,ve,dh,sec3,MAXSBX,MAXSBY,MAXSBZ,nisb, &          !jd
                   nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,nzsup,  &          !jd
                   zmnsup,zsizsup,                                   &          !jd
                   ncol_range,nrow_range,e_min,e_max,n_min,n_max,    &          !jd
                   elev_min,elev_max)                                           !jd
      call picksup_new(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,        &
                   isrot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,      &
                   iysbtosr,izsbtosr)
!
! Compute the number of drift terms, if an external drift is being
! considered then it is one more drift term, if SK is being considered
! then we will set all the drift terms off and mdt to 0):
!
      mdt = 1
      do i=1,9
            if(ktype.eq.0.or.ktype.eq.2) idrif(i) = 0
            if(idrif(i).lt.0.or.idrif(i).gt.1) then
                  write(*,*) 'ERROR KT3D: invalid drift term',idrif(i)
                  stop
            endif
            mdt = mdt + idrif(i)
      end do
      if(ktype.eq.3) mdt = mdt + 1
      if(ktype.eq.0) mdt = 0
      if(ktype.eq.2) mdt = 0
!
! Set up the discretization points per block.  Figure out how many
! are needed, the spacing, and fill the xdb,ydb, and zdb arrays with
! the offsets relative to the block center (this only gets done once):
!
! In all cases the offsets are relative to the lower left corner.
! This is done for rescaling the drift terms in the kriging matrix.
!
      if(nxdis.lt.1) nxdis = 1
      if(nydis.lt.1) nydis = 1
      if(nzdis.lt.1) nzdis = 1
      ndb = nxdis * nydis * nzdis
      if(ndb.gt.MAXDIS) then
            write(*,*) 'ERROR KT3D: Too many discretization points',ndb
            write(*,*) '            Increase MAXDIS or lower n[xyz]dis'
            stop
      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      zdis = zsiz  / max(real(nzdis),1.0)
      i    = 0
      xloc = -0.5*(xsiz+xdis)
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  zloc = -0.5*(zsiz+zdis)
                  do iz=1,nzdis
                        zloc = zloc + zdis
                        i = i+1
                        xdb(i) = xloc + 0.5*xsiz
                        ydb(i) = yloc + 0.5*ysiz
                        zdb(i) = zloc + 0.5*zsiz
                  end do
            end do
      end do
!
! Initialize accumulators:
!
      nk    = 0
      xk    = 0.0
      vk    = 0.0
      xkmae = 0.0
      xkmse = 0.0
!
! Calculate Block Covariance. Check for point kriging.
!
      call cova3_jd(xdb(1),ydb(1),zdb(1),xdb(1),ydb(1),zdb(1),1,nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov,PMX)
!
! Set the ``unbias'' variable so that the matrix solution is more stable
!

      unbias = cov
      cbb    = dble(cov)
      if(ndb.gt.1) then
            cbb = 0.0
            do i=1,ndb
               do j=1,ndb
                  call cova3_jd(xdb(i),ydb(i),zdb(i),xdb(j),ydb(j),zdb(j),                  &
                     1,nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov,PMX)
                  if(i.eq.j) cov = cov - c0(1)
                  cbb = cbb + dble(cov)
               end do
            end do
            cbb = cbb/dble(real(ndb*ndb))
      end if
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      end if
!
! Mean values of the drift functions:
!
      do i=1,9
            bv(i) = 0.0
      end do
      do i=1,ndb
            bv(1) = bv(1) + xdb(i)
            bv(2) = bv(2) + ydb(i)
            bv(3) = bv(3) + zdb(i)
            bv(4) = bv(4) + xdb(i)*xdb(i)
            bv(5) = bv(5) + ydb(i)*ydb(i)
            bv(6) = bv(6) + zdb(i)*zdb(i)
            bv(7) = bv(7) + xdb(i)*ydb(i)
            bv(8) = bv(8) + xdb(i)*zdb(i)
            bv(9) = bv(9) + ydb(i)*zdb(i)
      end do
      do i=1,9
            bv(i) = (bv(i) / real(ndb)) * resc
      end do
!
! Report on progress from time to time:
!
      if(koption.eq.0) then
            nxy   = nx*ny
            nxyz  = nx*ny*nz
            nloop = nxyz
            irepo = max(1,min((nxyz/10),10000))
      else
            nloop = 10000000
            irepo = max(1,min((nd/10),10000))
      end if
      ddh = 0.0
!jd      write(*,*)
!jd      write(*,*) 'Working on the kriging '

      nloop=n_npts                 !jd
!
! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!
!jd      do index=1,nloop
      do i_ipts=1,n_npts
!jd      if((int(index/irepo)*irepo).eq.index) write(*,103) index
!jd103   format('   currently on estimate ',i9)
!
! Where are we making an estimate?
!
      if(koption.eq.0) then
!jd            iz   = int((index-1)/nxy) + 1
!jd            iy   = int((index-(iz-1)*nxy-1)/nx) + 1
!jd            ix   = index - (iz-1)*nxy - (iy-1)*nx
!jd            xloc = xmn + real(ix-1)*xsiz
!jd            yloc = ymn + real(iy-1)*ysiz
!jd            zloc = zmn + real(iz-1)*zsiz

          xloc=epoint(i_ipts)                     !jd
          yloc=npoint(i_ipts)                     !jd
          zloc=zpoint(i_ipts)                     !jd
          icell=icellno(i_ipts)                   !jd
      else
            read(ljack,*,err=96,end=2) (var(i),i=1,nvarij)
            ddh  = 0.0
            xloc = xmn
            yloc = ymn
            zloc = zmn
            true = UNEST
            secj = UNEST
            if(idhlj.gt.0)  ddh    = var(idhlj)
            if(ixlj.gt.0)   xloc   = var(ixlj)
            if(iylj.gt.0)   yloc   = var(iylj)
            if(izlj.gt.0)   zloc   = var(izlj)
            if(ivrlj.gt.0)  true   = var(ivrlj)
            if(iextvj.gt.0) extest = var(iextvj)
            if(true.lt.tmin.or.true.ge.tmax) true = UNEST
      end if

!
! Read in the external drift variable for this grid node if needed:
!
      if(ktype.eq.2.or.ktype.eq.3) then
            if(koption.eq.0) then
                  read(lext,*) (var(i),i=1,iextve)
                  extest = var(iextve)
            end if
            if(extest.lt.tmin.or.extest.ge.tmax) then
                  est  = UNEST
                  estv = UNEST
                  go to 1
            end if
            resce  = covmax / max(extest,0.0001)
      endif
!
! Find the nearest samples:
!
      call srchsupr_new(xloc,yloc,zloc,radsqd,isrot,MAXROT,rotmat,nsbtosr,            &
                    ixsbtosr,iysbtosr,izsbtosr,noct,nd,x,y,z,tmp,                 &
                    nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,ysizsup,               &
                    nzsup,zmnsup,zsizsup,nclose,close,infoct)

!       write(6,*) ' nclose = ',nclose     !debug
!
! Load the nearest data in xa,ya,za,vra,vea:
!
      na = 0
      do i=1,nclose
            ind    = int(close(i)+0.5)
            accept = .true.
            if(koption.ne.0.and.(abs(x(ind)-xloc)+abs(y(ind)-yloc)+ abs(z(ind)-zloc))  &
                                 .lt.EPSLON) accept = .false.
            if(koption.ne.0.and.(abs(dh(ind)-ddh)).lt.EPSLON) accept = .false.
            if(accept) then
                  if(na.lt.ndmax) then
                        na = na + 1
!jd                        xa(na)  = x(ind) - xloc + 0.5*xsiz
!jd                        ya(na)  = y(ind) - yloc + 0.5*ysiz
!jd                        za(na)  = z(ind) - zloc + 0.5*zsiz
                        xa(na)  = x(ind) - xloc
                        ya(na)  = y(ind) - yloc
                        za(na)  = z(ind) - zloc
                        vra(na) = vr(ind)
                        vea(na) = ve(ind)
                  end if
            end if
      end do
!
! Test number of samples found:
!
      if(na.lt.ndmin) then
            est  = UNEST
            estv = UNEST
            go to 1
      end if
!
! Test if there are enough samples to estimate all drift terms:
!
!jd      if(na.ge.1.and.na.le.mdt) then
!jd            if(fircon) then
!jd                  write(ldbg,999)
!jd                  write(6,999)
!jd                  fircon = .false.
!jd            end if
!jd            est  = UNEST
!jd            estv = UNEST
!jd            go to 1
!jd      end if
 999  format(' Encountered a location where there were too few data ',/,     &
             ' to estimate all of the drift terms but there would be',/,     &
             ' enough data for OK or SK.   KT3D currently leaves ',/,        &
             ' these locations unestimated.',/,                              &
             ' This message is only written once - the first time.',/)
!
! There are enough samples - proceed with estimation.
!
      if(na.le.1) then
!
! Handle the situation of only one sample:
!
            call cova3_jd(xa(1),ya(1),za(1),xa(1),ya(1),za(1),1,nst,MAXNST,     &
                       c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb1,PMX)
!
! Establish Right Hand Side Covariance:
!
            if(ndb.le.1) then
                  zero=0.0                                              !jd
!jd                  call cova3_jd(xa(1),ya(1),za(1),xdb(1),ydb(1),zdb(1),1,       &
!jd                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb,PMX)
                  call cova3_jd(xa(1),ya(1),za(1),zero,zero,zero,1,       &          !jd
                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb,PMX)          !jd

            else
                  cb  = 0.0
                  do i=1,ndb
                        call cova3_jd(xa(1),ya(1),za(1),xdb(i),ydb(i),          &
                                   zdb(i),1,nst,MAXNST,c0,it,cc,aa,1,        &
                                   MAXROT,rotmat,cmax,cov,PMX)
                        cb = cb + cov
                        dx = xa(1) - xdb(i)
                        dy = ya(1) - ydb(i)
                        dz = za(1) - zdb(i)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
!vd
!
! Early bug - always did OK in presence of one data.
!
!vd
            if(ktype.eq.2) skmean = extest
            if(ktype.eq.0.or.ktype.eq.2) then
                  wt   = cb / cb1
                  est  = wt * vra(1) + (1.0-wt) * skmean
                  estv = real(cbb) - wt*cb
            else
                  est  = vra(1)
                  estv = real(cbb) - 2.0*cb + cb1
            end if

            if(ktype.eq.0)then                                    !jd
              rrtemp=(1.0-wt)*skmean                              !jd
            else                                                  !jd
              rrtemp=0.0                                          !jd
            end if                                                !jd
            if(aoutfile.eq.'f')then                               !jd
              write(outunit,*) icell,itrans,1,rrtemp, &           !jd
              inumdat(nint(sec3(nint(close(1))))),wt                           !jd
            else                                                  !jd
              write(outunit)   icell,itrans,1,rrtemp, &           !jd
              inumdat(nint(sec3(nint(close(1))))),wt                           !jd
            end if                                                !jd

            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
            go to 1
      end if
!
! Go ahead and set up the OK portion of the kriging matrix:
!
      neq = mdt+na
!
! Initialize the main kriging matrix:
!
      first = .false.
      do i=1,neq*neq
            a(i) = 0.0
      end do
!
! Fill in the kriging matrix:
!
      do i=1,na
      do j=i,na
            call cova3_jd(xa(i),ya(i),za(i),xa(j),ya(j),za(j),1,nst,MAXNST,   &
                       c0,it,cc,aa,1,MAXROT,rotmat,cmax,cov,PMX)
            a(neq*(i-1)+j) = dble(cov)
            a(neq*(j-1)+i) = dble(cov)
      end do
      end do
!
! Fill in the OK unbiasedness portion of the matrix (if not doing SK):
!
      if(neq.gt.na) then
            do i=1,na
                  a(neq*(i-1)+na+1) = dble(unbias)
                  a(neq*na+i)       = dble(unbias)
            end do
      endif
!
! Set up the right hand side:
!
      do i=1,na
            if(ndb.le.1) then
                  zero=0.0
!jd                  call cova3_jd(xa(i),ya(i),za(i),xdb(1),ydb(1),zdb(1),1,     &
!jd                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb,PMX)
                  call cova3_jd(xa(i),ya(i),za(i),zero,zero,zero,1,     &
                       nst,MAXNST,c0,it,cc,aa,1,MAXROT,rotmat,cmax,cb,PMX)
            else
                  cb  = 0.0
                  do j=1,ndb
                        call cova3_jd(xa(i),ya(i),za(i),xdb(j),ydb(j),        &
                                   zdb(j),1,nst,MAXNST,c0,it,cc,aa,1,        &
                                   MAXROT,rotmat,cmax,cov,PMX)
                        cb = cb + cov
                        dx = xa(i) - xdb(j)
                        dy = ya(i) - ydb(j)
                        dz = za(i) - zdb(j)
                        if((dx*dx+dy*dy+dz*dz).lt.EPSLON) cb=cb-c0(1)
                  end do
                  cb = cb / real(ndb)
            end if
            r(i) = dble(cb)
      end do
      if(neq.gt.na) r(na+1) = dble(unbias)
!
! Add the additional unbiasedness constraints:
!
      im = na + 1
!
! First drift term (linear in "x"):
!
      if(idrif(1).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*resc)
            end do
            r(im) = dble(bv(1))
      endif
!
! Second drift term (linear in "y"):
!
      if(idrif(2).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*resc)
            end do
            r(im) = dble(bv(2))
      endif
!
! Third drift term (linear in "z"):
!
      if(idrif(3).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*resc)
            end do
            r(im) = dble(bv(3))
      endif
!
! Fourth drift term (quadratic in "x"):
!
      if(idrif(4).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*xa(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*xa(k)*resc)
            end do
            r(im) = dble(bv(4))
      endif
!
! Fifth drift term (quadratic in "y"):
!
      if(idrif(5).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(5))
      endif
!
! Sixth drift term (quadratic in "z"):
!
      if(idrif(6).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(za(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(za(k)*za(k)*resc)
            end do
            r(im) = dble(bv(6))
      endif
!
! Seventh drift term (quadratic in "xy"):
!
      if(idrif(7).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*ya(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*ya(k)*resc)
            end do
            r(im) = dble(bv(7))
      endif
!
! Eighth drift term (quadratic in "xz"):
!
      if(idrif(8).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(xa(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(xa(k)*za(k)*resc)
            end do
            r(im) = dble(bv(8))
      endif
!
! Ninth drift term (quadratic in "yz"):
!
      if(idrif(9).eq.1) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(ya(k)*za(k)*resc)
                  a(neq*(k-1)+im) = dble(ya(k)*za(k)*resc)
            end do
            r(im) = dble(bv(9))
      endif
!
! External drift term (specified by external variable):
!
      if(ktype.eq.3) then
            im=im+1
            do k=1,na
                  a(neq*(im-1)+k) = dble(vea(k)*resce)
                  a(neq*(k-1)+im) = dble(vea(k)*resce)
            end do
            r(im) = dble(extest*resce)
      endif
!
! Copy the right hand side to compute the kriging variance later:
!
      do k=1,neq
            rr(k) = r(k)
      end do
      kadim = neq * neq
      ksdim = neq
      nrhs  = 1
      nv    = 1
!
! If estimating the trend then reset all the right hand side terms=0.0:
!
      if(itrend.ge.1) then
            do i=1,na
                  r(i)  = 0.0
                  rr(i) = 0.0
            end do
      endif
!
! Write out the kriging Matrix if Seriously Debugging:
!
      if(idbg.eq.3) then
            write(ldbg,*) 'Estimating node index : ',ix,iy,iz
            is = 1 - neq
            do i=1,neq
                  is = 1 + (i-1)*neq
                  ie = is + neq - 1
                  write(ldbg,100) i,r(i),(a(j),j=is,ie)
 100              format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
            end do
      endif
!
! Solve the kriging system:
!
      call ktsol_new(neq,nrhs,nv,a,r,s,ising,maxeq)
!
! Compute the solution:
!
      if(ising.ne.0) then
            write(6,450) icell                                            !jd
450         format(' WARNING: singular kriging matrix for cell',i5)       !jd
            if(idbg.ge.3) write(ldbg,*) ' Singular Matrix ',ix,iy,iz
            est  = UNEST
            estv = UNEST
      else
            est  = 0.0
            estv = real(cbb)
            sumw=0.0                                                 !jd
            if(ktype.eq.2) skmean = extest
            do j=1,neq
                  estv = estv - real(s(j))*rr(j)
                  if(j.le.na) then
                        if(ktype.eq.0) then
                              est = est + real(s(j))*(vra(j)-skmean)
                              sumw=sumw+real(s(j))                   !jd
                        else if(ktype.eq.2) then
                              est = est + real(s(j))*(vra(j)-vea(j))
                        else
                              est = est + real(s(j))*vra(j)
                        endif
                  endif
            end do
            if(ktype.eq.0.or.ktype.eq.2) est = est + skmean
            nk   = nk + 1
            xk   = xk + est
            vk   = vk + est*est
            if(ktype.eq.0)then                                       !jd
              rrtemp=(1.0-sumw)*skmean                               !jd
            else                                                     !jd
              rrtemp=0.0                                             !jd
            end if                                                   !jd
            if(aoutfile.eq.'f')then                                  !jd
              write(outunit,*) icell,itrans,na,rrtemp, &             !jd
              ((inumdat(nint(sec3(nint(close(i))))),real(s(i))),i=1,na)           !jd
            else                                                     !jd
              write(outunit)   icell,itrans,na,rrtemp, &             !jd
              ((inumdat(nint(sec3(nint(close(i))))),real(s(i))),i=1,na)           !jd
            end if                                                   !jd

!
! Write the kriging weights and data if debugging level is above 2:
!
            if(idbg.ge.2) then
                  write(ldbg,*) '       '
                  write(ldbg,*) 'BLOCK: ',ix,iy,iz,' at ',xloc,yloc,zloc
                  write(ldbg,*) '       '
                  if(ktype.ne.0)write(ldbg,*) '  Lagrange : ',s(na+1)*unbias
                  write(ldbg,*) '  BLOCK EST: x,y,z,vr,wt '
                  do i=1,na
                        xa(i) = xa(i) + xloc - 0.5*xsiz
                        ya(i) = ya(i) + yloc - 0.5*ysiz
                        za(i) = za(i) + zloc - 0.5*zsiz
                        write(ldbg,'(5f12.3)') xa(i),ya(i),za(i),vra(i),s(i)
                  end do
                  write(ldbg,*) '  estimate, variance  ',est,estv
            endif
      endif
!
! END OF MAIN KRIGING LOOP:
!
 1          continue
            if(iktype.eq.0) then
                  if(koption.eq.0) then
!jd                        write(lout,'(g14.8,1x,g14.8)') est,estv
                  else
                        err = UNEST
                        if(true.ne.UNEST.and.est.ne.UNEST) then
                              err=est-true
                              xkmae = xkmae + abs(err)
                              xkmse = xkmse + err*err
                        end if
!jd                        write(lout,'(7(g14.8,1x))') xloc,yloc,zloc,true,est,estv,err
                  end if
            else
!
! Work out the IK-type distribution implicit to this data configuration
! and kriging weights:
!
                  do icut=1,ncut
                        cdf(icut) = -1.0
                  end do
                  wtmin = 1.0
                  do i=1,na
                        if(s(i).lt.wtmin) wtmin = s(i)
                  end do
                  sumwt = 0.0
                  do i=1,na
                        s(i)  = s(i) - wtmin
                        sumwt = sumwt + s(i)
                  end do
                  do i=1,na
                        s(i) = s(i) / max(0.00001,sumwt)
                  end do
                  if(na.gt.1.and.sumwt.gt.0.00001) then
                        do icut=1,ncut
                              cdf(icut) = 0.0
                              do i=1,na
                                    if(vra(i).le.cut(icut))cdf(icut)=cdf(icut)+s(i)
                              end do
                        end do
                  end if
                  if(koption.eq.0) then
!jd                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut)
                  else
!jd                        write(lout,'(30(f8.4))') (cdf(i),i=1,ncut),true
                  end if
            end if
      end do
 2    continue
      if(koption.gt.0) close(ljack)
!
! Write statistics of kriged values:
!

!jd      if(nk.gt.0.and.idbg.gt.0) then
!jd            xk    = xk/real(nk)
!jd            vk    = vk/real(nk) - xk*xk
!jd            xkmae = xkmae/real(nk)
!jd            xkmse = xkmse/real(nk)
!jd            write(ldbg,105) nk,xk,vk
!jd            write(*,   105) nk,xk,vk
!jd 105        format(/,'Estimated   ',i8,' blocks ',/,
!jd     +               '  average   ',g14.8,/,'  variance  ',g14.8,/)
!jd            if(koption.ne.0) then
!jd                  write(*,106) xkmae,xkmse
!jd 106              format(/,'  mean error',g14.8,/,'  mean sqd e',g14.8)
!jd            end if
!jd      endif
!
! All finished the kriging:
!
      return
 96   stop 'ERROR in jackknife file!'
      end




      subroutine setsupr_jd(nx,xmn,xsiz,ny,ymn,ysiz,nz,zmn,zsiz,nd,x,y,z,     &
                         vr,tmp,nsec,sec1,sec2,sec3,MAXSBX,MAXSBY,            &
                         MAXSBZ,nisb,nxsup,xmnsup,xsizsup,nysup,ymnsup,       &
                         ysizsup,nzsup,zmnsup,zsizsup,                        &
                         ncol_range,nrow_range,e_min,e_max,n_min,n_max,       &
                         elev_min,elev_max)
!-----------------------------------------------------------------------
!
!           Establish Super Block Search Limits and Sort Data
!           *************************************************
!
! This subroutine sets up a 3-D "super block" model and orders the data
! by super block number.  The limits of the super block is set to the
! minimum and maximum limits of the grid; data outside are assigned to
! the nearest edge block.
!
! The idea is to establish a 3-D block network that contains all the
! relevant data.  The data are then sorted by their index location in
! the search network, i.e., the index location is given after knowing
! the block index in each coordinate direction (ix,iy,iz):
!          ii = (iz-1)*nxsup*nysup + (iy-1)*nxsup + ix
! An array, the same size as the number of super blocks, is constructed
! that contains the cumulative number of data in the model.  With this
! array it is easy to quickly check what data are located near any given
! location.
!
!
!
! INPUT VARIABLES:
!
!   nx,xmn,xsiz      Definition of the X grid being considered
!   ny,ymn,ysiz      Definition of the Y grid being considered
!   nz,zmn,zsiz      Definition of the Z grid being considered
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   vr(nd)           Variable at each location.
!   tmp(nd)          Temporary storage to keep track of the super block
!                      index associated to each data (uses the same
!                      storage already allocated for the simulation)
!   nsec             Number of secondary variables to carry with vr
!   sec1(nd)         First secondary variable (if nsec >= 1)
!   sec2(nd)         Second secondary variable (if nsec >= 2)
!   sec3(nd)         Third secondary variable (if nsec = 3)
!   MAXSB[X,Y,Z]     Maximum size of super block network
!
!
!
! OUTPUT VARIABLES:
!
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the Y super block grid
!   nzsup,zmnsup,zsizsup  Definition of the Z super block grid
!
!
!
! EXTERNAL REFERENCES:
!
!   sortem           Sorting routine to sort the data
!
!
!
!-----------------------------------------------------------------------

      implicit integer(i-n), real(a-h, o-z)                    !jd

      real    x(*),y(*),z(*),vr(*),tmp(*),sec1(*),sec2(*),sec3(*)
      integer nisb(*)
      logical inflag

      integer ncol_range,nrow_range                            !jd
      real e_min,e_max,n_min,n_max,elev_min,elev_max           !jd
!
! Establish the number and size of the super blocks:
!
!jd      nxsup   = min(nx,MAXSBX)
!jd      nysup   = min(ny,MAXSBY)
      nxsup=min(ncol_range,MAXSBX)                             !jd
      nysup=min(nrow_range,MAXSBY)                             !jd
      nzsup   = min(nz,MAXSBZ)
!jd      xsizsup = real(nx)*xsiz/real(nxsup)
!jd      ysizsup = real(ny)*ysiz/real(nysup)
!jd      zsizsup = real(nz)*zsiz/real(nzsup)

      xsizsup = (e_max-e_min)*1.1/real(nxsup)                      !jd
      ysizsup = (n_max-n_min)*1.1/real(nysup)                      !jd
      zsizsup = (elev_max-elev_min)*1.1                            !jd

!jd      xmnsup  = (xmn-0.5*xsiz)+0.5*xsizsup
!jd      ymnsup  = (ymn-0.5*ysiz)+0.5*ysizsup
!jd      zmnsup  = (zmn-0.5*zsiz)+0.5*zsizsup

      xmnsup  = e_min-0.05*(e_max-e_min)                           !jd
      ymnsup  = n_min-0.05*(n_max-n_min)                           !jd
      zmnsup  = elev_min-0.05*(elev_max-elev_min)                  !jd

!
! Initialize the extra super block array to zeros:
!
      do i=1,nxsup*nysup*nzsup
            nisb(i) = 0
      end do
!
! Loop over all the data assigning the data to a super block and
! accumulating how many data are in each super block:
!
      do i=1,nd
            call getindx_new(nxsup,xmnsup,xsizsup,x(i),ix,inflag)
            call getindx_new(nysup,ymnsup,ysizsup,y(i),iy,inflag)
            call getindx_new(nzsup,zmnsup,zsizsup,z(i),iz,inflag)
            ii = ix + (iy-1)*nxsup + (iz-1)*nxsup*nysup
            tmp(i)   = ii
            nisb(ii) = nisb(ii) + 1
      end do
!
! Sort the data by ascending super block number:
!
      nsort = 4 + nsec
      call sortem_new(1,nd,tmp,nsort,x,y,z,vr,sec1,sec2,sec3)
!
! Set up array nisb with the starting address of the block data:
!
      do i=1,(nxsup*nysup*nzsup-1)
            nisb(i+1) = nisb(i) + nisb(i+1)
      end do
!
! Finished:
!
      return
      end



      subroutine cova3_jd(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,   &
                       irot,MAXROT,rotmat,cmax,cova,PMX)
!-----------------------------------------------------------------------
!
!                    Covariance Between Two Points
!                    *****************************
!
! This subroutine calculated the covariance associated with a variogram
! model specified by a nugget effect and nested varigoram structures.
! The anisotropy definition can be different for each nested structure.
!
!
!
! INPUT VARIABLES:
!
!   x1,y1,z1         coordinates of first point
!   x2,y2,z2         coordinates of second point
!   nst(ivarg)       number of nested structures (maximum of 4)
!   ivarg            variogram number (set to 1 unless doing cokriging
!                       or indicator kriging)
!   MAXNST           size of variogram parameter arrays
!   c0(ivarg)        isotropic nugget constant
!   it(i)            type of each nested structure:
!                      1. spherical model of range a;
!                      2. exponential model of parameter a;
!                           i.e. practical range is 3a
!                      3. gaussian model of parameter a;
!                           i.e. practical range is a*sqrt(3)
!                      4. power model of power a (a must be gt. 0  and
!                           lt. 2).  if linear model, a=1,c=slope.
!                      5. hole effect model
!   cc(i)            multiplicative factor of each nested structure.
!                      (sill-c0) for spherical, exponential,and gaussian
!                      slope for linear model.
!   aa(i)            parameter "a" of each nested structure.
!   irot             index of the rotation matrix for the first nested
!                    structure (the second nested structure will use
!                    irot+1, the third irot+2, and so on)
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!
!
! OUTPUT VARIABLES:
!
!   cmax             maximum covariance
!   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
!
!
!
! EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
!                      rotmat    computes rotation matrix for distance
!-----------------------------------------------------------------------

      implicit integer(i-n), real(a-h, o-z)                    !jd

!jd      parameter(PI=3.14159265,PMX=999.,EPSLON=1.e-5)
      parameter(PI=3.14159265,EPSLON=1.e-5)

      integer   nst(*),it(*)
      real      c0(*),cc(*),aa(*)
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist_new

      real      PMX                                         !jd
!
! Calculate the maximum covariance value (used for zero distances and
! for power model covariance):
!
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0(ivarg)
      do is=1,nst(ivarg)
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
!
! Check for "zero" distance, return with cmax if so:
!
      hsqd = sqdist_new(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
!
! Loop over all the structures:
!
      cova = 0.0
      do is=1,nst(ivarg)
            ist = istart + is - 1
!
! Compute the appropriate distance:
!
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist_new(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
!
! Spherical Variogram Model?
!
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
!
! Exponential Variogram Model?
!
            else if(it(ist).eq.2) then
!jd                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
                  cova = cova + cc(ist)*exp(-1.0*h/aa(ist))                   !jd
!
! Gaussian Variogram Model?
!
            else if(it(ist).eq.3) then
!jd                  cova = cova + cc(ist)*exp(-3.*(h/aa(ist))*(h/aa(ist)))
                  cova = cova + cc(ist)*exp(-1.*(h/aa(ist))*(h/aa(ist)))      !jd
!
! Power Variogram Model?
!
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
!
! Hole Effect Model?
!
            else if(it(ist).eq.5) then
!                 d = 10.0 * aa(ist)
!                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do
!
! Finished:
!
      return
      end


      subroutine ktsol_new(n,ns,nv,a,b,x,ktilt,maxeq)
!-----------------------------------------------------------------------
!
! Solution of a system of linear equations by gaussian elimination with
! partial pivoting.  Several right hand side matrices and several
! variables are allowed.
!
!
!         NOTE: All input matrices must be in double precision
!
!
! INPUT/OUTPUT VARIABLES:
!
!   n                Number of equations
!   ns               Number of right hand side matrices
!   nv               Number of variables.
!   a(n*n*nv)        left hand side matrices versus columnwise.
!   b(n*ns*nv)       input right hand side matrices.
!   x(n*ns*nv)       solution matrices.
!   ktilt            indicator of singularity
!                      =  0  everything is ok.
!                      = -1 n.le.1
!                      =  k  a null pivot appeared at the kth iteration.
!   tol              used in test for null pivot. depends on machine
!                      precision and can also be set for the tolerance
!                      of an ill-defined kriging system.
!
!
!-----------------------------------------------------------------------


!jd      implicit real*8 (a-h,o-z)

      implicit integer(i-n), real*8(a-h, o-z)                    !jd

      real*8 x(maxeq),a(maxeq*maxeq),b(maxeq)
!
! Make sure there are equations to solve:
!
      if(n.le.1) then
            ktilt = -1
            return
      endif
!
! Initialization:
!
      tol   = 0.1e-10
      ktilt = 0
      ntn   = n*n
      nm1   = n-1
!
! Triangulation is done variable by variable:
!
      do iv=1,nv
!
! Indices of location in vectors a and b:
!
            nva = ntn*(iv-1)
            nvb = n*ns*(iv-1)
!
! Gaussian elimination with partial pivoting:
!
            do k=1,nm1
                  kp1 = k+1
!
! Indice of the diagonal element in the kth row:
!
                  kdiag = nva+(k-1)*n+k
!
! Find the pivot - interchange diagonal element/pivot:
!
                  npiv = kdiag
                  ipiv = k
                  i1   = kdiag
                  do i=kp1,n
                        i1 = i1+1
                        if(abs(a(i1)).gt.abs(a(npiv))) then
                              npiv = i1
                              ipiv = i
                        endif
                  end do
                  t        = a(npiv)
                  a(npiv)  = a(kdiag)
                  a(kdiag) = t
!
! Test for singularity:
!
                  if(abs(a(kdiag)).lt.tol) then
                        ktilt=k
                        return
                  endif
!
! Compute multipliers:
!
                  i1 = kdiag
                  do i=kp1,n
                        i1    = i1+1
                        a(i1) = -a(i1)/a(kdiag)
                  end do
!
! Interchange and eliminate column per column:
!
                  j1 = kdiag
                  j2 = npiv
                  do j=kp1,n
                        j1    = j1+n
                        j2    = j2+n
                        t     = a(j2)
                        a(j2) = a(j1)
                        a(j1) = t
                        i1    = j1
                        i2    = kdiag
                        do i=kp1,n
                              i1    = i1+1
                              i2    = i2+1
                              a(i1) = a(i1)+a(i2)*a(j1)
                        end do
                  end do
!
! Interchange and modify the ns right hand matrices:
!
                  i1 = nvb+ipiv
                  i2 = nvb+k
                  do i=1,ns
                        t     = b(i1)
                        b(i1) = b(i2)
                        b(i2) = t
                        j1    = i2
                        j2    = kdiag
                        do j=kp1,n
                              j1    = j1+1
                              j2    = j2+1
                              b(j1) = b(j1)+b(i2)*a(j2)
                        end do
                        i1 = i1+n
                        i2 = i2+n
                  end do
            end do
!
! Test for singularity for the last pivot:
!
            kdiag = ntn*iv
            if(abs(a(kdiag)).lt.tol) then
                  ktilt = n
                  return
            endif
      end do
!
! End of triangulation. Now, solve back variable per variable:
!
      do iv=1,nv
!
! Indices of location in vectors a and b:
!
            nva  = ntn*iv
            nvb1 = n*ns*(iv-1)+1
            nvb2 = n*ns*iv
!
! Back substitution with the ns right hand matrices:
!
            do il=1,ns
                  do k=1,nm1
                        nmk = n-k
!
! Indice of the diagonal element of the (n-k+1)th row and of
! the (n-k+1)th element of the left hand side.
!
                        kdiag = nva-(n+1)*(k-1)
                        kb    = nvb2-(il-1)*n-k+1
                        b(kb) = b(kb)/a(kdiag)
                        t     = -b(kb)
                        i1    = kb
                        i2    = kdiag
                        do i=1,nmk
                              i1    = i1-1
                              i2    = i2-1
                              b(i1) = b(i1)+a(i2)*t
                        end do
                  end do
                  kdiag = kdiag-n-1
                  kb    = kb-1
                  b(kb) = b(kb)/a(kdiag)
            end do
!
! End of back substitution:
!
      end do
!
! Restitution of the solution:
!
      itot = n*ns*nv
      do i=1,itot
            x(i) = b(i)
      end do
!
! Finished:
!
      return
      end



      subroutine srchsupr_new(xloc,yloc,zloc,radsqd,irot,MAXROT,rotmat,    &
                          nsbtosr,ixsbtosr,iysbtosr,izsbtosr,noct,nd,      &
                          x,y,z,tmp,nisb,nxsup,xmnsup,xsizsup,             &
                          nysup,ymnsup,ysizsup,nzsup,zmnsup,zsizsup,       &
                          nclose,close,infoct)
!-----------------------------------------------------------------------
!
!              Search Within Super Block Search Limits
!              ***************************************
!
!
! This subroutine searches through all the data that have been tagged in
! the super block subroutine.  The close data are passed back in the
! index array "close".  An octant search is allowed.
!
!
!
! INPUT VARIABLES:
!
!   xloc,yloc,zloc   location of point being estimated/simulated
!   radsqd           squared search radius
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!   noct             If >0 then data will be partitioned into octants
!   nd               Number of data
!   x(nd)            X coordinates of the data
!   y(nd)            Y coordinates of the data
!   z(nd)            Z coordinates of the data
!   tmp(nd)          Temporary storage to keep track of the squared
!                      distance associated with each data
!   nisb()                Array with cumulative number of data in each
!                           super block.
!   nxsup,xmnsup,xsizsup  Definition of the X super block grid
!   nysup,ymnsup,ysizsup  Definition of the X super block grid
!   nzsup,zmnsup,zsizsup  Definition of the X super block grid
!
!
!
! OUTPUT VARIABLES:
!
!   nclose           Number of close data
!   close()          Index of close data
!   infoct           Number of informed octants (only computes if
!                      performing an octant search)
!
!
!
! EXTERNAL REFERENCES:
!
!   sqdist           Computes anisotropic squared distance
!   sortem           Sorts multiple arrays in ascending order
!
!
!
!-----------------------------------------------------------------------

      implicit integer(i-n), real(a-h, o-z)                    !jd

      real    x(*),y(*),z(*),tmp(*),close(*)
      real*8  rotmat(MAXROT,3,3),hsqd,sqdist_new
      integer nisb(*),inoct(8)
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
      logical inflag
!
! Determine the super block location of point being estimated:
!
      call getindx_new(nxsup,xmnsup,xsizsup,xloc,ix,inflag)
      call getindx_new(nysup,ymnsup,ysizsup,yloc,iy,inflag)
      call getindx_new(nzsup,zmnsup,zsizsup,zloc,iz,inflag)
!
! Loop over all the possible Super Blocks:
!
      nclose = 0
      do 1 isup=1,nsbtosr
!
! Is this super block within the grid system:
!
            ixsup = ix + ixsbtosr(isup)
            iysup = iy + iysbtosr(isup)
            izsup = iz + izsbtosr(isup)
            if(ixsup.le.0.or.ixsup.gt.nxsup.or.                      &
               iysup.le.0.or.iysup.gt.nysup.or.			     &
               izsup.le.0.or.izsup.gt.nzsup) go to 1
!
! Figure out how many samples in this super block:
!
            ii = ixsup + (iysup-1)*nxsup + (izsup-1)*nxsup*nysup
            if(ii.eq.1) then
                  nums = nisb(ii)
                  i    = 0
            else
                  nums = nisb(ii) - nisb(ii-1)
                  i    = nisb(ii-1)
            endif
!
! Loop over all the data in this super block:
!
            do 2 ii=1,nums
                  i = i + 1
!
! Check squared distance:
!
                  hsqd = sqdist_new(xloc,yloc,zloc,x(i),y(i),z(i),irot,    &
                                MAXROT,rotmat)
                  if(real(hsqd).gt.radsqd) go to 2
!
! Accept this sample:
!
                  nclose = nclose + 1
                  close(nclose) = real(i)
                  tmp(nclose)  = real(hsqd)
 2          continue
 1    continue
!
! Sort the nearby samples by distance to point being estimated:
!
      call sortem_new(1,nclose,tmp,1,close,c,d,e,f,g,h)
!
! If we aren't doing an octant search then just return:
!
      if(noct.le.0) return
!
! PARTITION THE DATA INTO OCTANTS:
!
      do i=1,8
            inoct(i) = 0
      end do
!
! Now pick up the closest samples in each octant:
!
      nt = 8*noct
      na = 0
      do j=1,nclose
            i  = int(close(j))
            h  = tmp(j)
            dx = x(i) - xloc
            dy = y(i) - yloc
            dz = z(i) - zloc
            if(dz.lt.0.) go to 5
            iq=4
            if(dx.le.0.0 .and. dy.gt.0.0) iq=1
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=2
            if(dx.lt.0.0 .and. dy.le.0.0) iq=3
            go to 6
 5          iq=8
            if(dx.le.0.0 .and. dy.gt.0.0) iq=5
            if(dx.gt.0.0 .and. dy.ge.0.0) iq=6
            if(dx.lt.0.0 .and. dy.le.0.0) iq=7
 6          continue
            inoct(iq) = inoct(iq) + 1
!
! Keep this sample if the maximum has not been exceeded:
!
            if(inoct(iq).le.noct) then
                  na = na + 1
                  close(na) = i
                  tmp(na)   = h
                  if(na.eq.nt) go to 7
            endif
      end do
!
! End of data selection. Compute number of informed octants and return:
!
 7    nclose = na
      infoct = 0
      do i=1,8
            if(inoct(i).gt.0) infoct = infoct + 1
      end do
!
! Finished:
!
      return
      end



      subroutine getindx_new(n,min,siz,loc,index,inflag)
!-----------------------------------------------------------------------
!
!     Gets the coordinate index location of a point within a grid
!     ***********************************************************
!
!
! n       number of "nodes" or "cells" in this coordinate direction
! min     origin at the center of the first cell
! siz     size of the cells
! loc     location of the point being considered
! index   output index within [1,n]
! inflag  true if the location is actually in the grid (false otherwise
!         e.g., if the location is outside then index will be set to
!         nearest boundary
!
!
!
!-----------------------------------------------------------------------
      integer   n,index
      real      min,siz,loc
      logical   inflag
!
! Compute the index of "loc":
!
      index = int( (loc-min)/siz + 1.5 )
!
! Check to see if in or out:
!
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
!
! Return to calling program:
!
      return
      end


      subroutine setrot_new(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
!-----------------------------------------------------------------------
!
!              Sets up an Anisotropic Rotation Matrix
!              **************************************
!
! Sets up the matrix to transform cartesian coordinates to coordinates
! accounting for angles and anisotropy (see manual for a detailed
! definition):
!
!
! INPUT PARAMETERS:
!
!   ang1             Azimuth angle for principal direction
!   ang2             Dip angle for principal direction
!   ang3             Third rotation angle
!   anis1            First anisotropy ratio
!   anis2            Second anisotropy ratio
!   ind              matrix indicator to initialize
!   MAXROT           maximum number of rotation matrices dimensioned
!   rotmat           rotation matrices
!
!
! NO EXTERNAL REFERENCES
!
!
!-----------------------------------------------------------------------
      implicit integer(i-n), real(a-h, o-z)                    !jd

      parameter(DEG2RAD=3.141592654/180.0,EPSLON=1.e-20)
      real*8    rotmat(MAXROT,3,3),afac1,afac2,sina,sinb,sint,cosa,cosb,cost
!
! Converts the input angles to three angles which make more
!  mathematical sense:
!
!         alpha   angle between the major axis of anisotropy and the
!                 E-W axis. Note: Counter clockwise is positive.
!         beta    angle between major axis and the horizontal plane.
!                 (The dip of the ellipsoid measured positive down)
!         theta   Angle of rotation of minor axis about the major axis
!                 of the ellipsoid.
!
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
!
! Get the required sines and cosines:
!
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
!
! Construct the rotation matrix in the required memory:
!
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
!
! Return to calling program:
!
      return
      end




      subroutine picksup_new(nxsup,xsizsup,nysup,ysizsup,nzsup,zsizsup,    &
                         irot,MAXROT,rotmat,radsqd,nsbtosr,ixsbtosr,	   &
                         iysbtosr,izsbtosr)
!-----------------------------------------------------------------------
!
!             Establish Which Super Blocks to Search
!             **************************************
!
! This subroutine establishes which super blocks must be searched given
! that a point being estimated/simulated falls within a super block
! centered at 0,0,0.
!
!
!
! INPUT VARIABLES:
!
!   nxsup,xsizsup    Definition of the X super block grid
!   nysup,ysizsup    Definition of the Y super block grid
!   nzsup,zsizsup    Definition of the Z super block grid
!   irot             index of the rotation matrix for searching
!   MAXROT           size of rotation matrix arrays
!   rotmat           rotation matrices
!   radsqd           squared search radius
!
!
!
! OUTPUT VARIABLES:
!
!   nsbtosr          Number of super blocks to search
!   ixsbtosr         X offsets for super blocks to search
!   iysbtosr         Y offsets for super blocks to search
!   izsbtosr         Z offsets for super blocks to search
!
!
!
! EXTERNAL REFERENCES:
!
!   sqdist           Computes anisotropic squared distance
!
!
!
!-----------------------------------------------------------------------
      implicit integer(i-n), real(a-h, o-z)                    !jd

      real*8  rotmat(MAXROT,3,3),hsqd,sqdist_new,shortest
      integer ixsbtosr(*),iysbtosr(*),izsbtosr(*)
!
! MAIN Loop over all possible super blocks:
!
      nsbtosr = 0
      do i=-(nxsup-1),(nxsup-1)
      do j=-(nysup-1),(nysup-1)
      do k=-(nzsup-1),(nzsup-1)
            xo = real(i)*xsizsup
            yo = real(j)*ysizsup
            zo = real(k)*zsizsup
!
! Find the closest distance between the corners of the super blocks:
!
            shortest = 1.0e21
            do i1=-1,1
            do j1=-1,1
            do k1=-1,1
                  do i2=-1,1
                  do j2=-1,1
                  do k2=-1,1
                        if(i1.ne.0.and.j1.ne.0.and.k1.ne.0.and.          &
                           i2.ne.0.and.j2.ne.0.and.k2.ne.0) then
                              xdis = real(i1-i2)*0.5*xsizsup + xo
                              ydis = real(j1-j2)*0.5*ysizsup + yo
                              zdis = real(k1-k2)*0.5*zsizsup + zo
                              hsqd = sqdist_new(0.0,0.0,0.0,xdis,ydis,zdis,   &
                                            irot,MAXROT,rotmat)
                              if(hsqd.lt.shortest) shortest = hsqd
                        end if
                  end do
                  end do
                  end do
            end do
            end do
            end do
!
! Keep this super block if it is close enoutgh:
!
            if(real(shortest).le.radsqd) then
                  nsbtosr = nsbtosr + 1
                  ixsbtosr(nsbtosr) = i
                  iysbtosr(nsbtosr) = j
                  izsbtosr(nsbtosr) = k
            end if
      end do
      end do
      end do
!
! Finished:
!
      return
      end



      real*8 function sqdist_new(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
!-----------------------------------------------------------------------
!
!    Squared Anisotropic Distance Calculation Given Matrix Indicator
!    ***************************************************************
!
! This routine calculates the anisotropic distance between two points
!  given the coordinates of each point and a definition of the
!  anisotropy.
!
!
! INPUT VARIABLES:
!
!   x1,y1,z1         Coordinates of first point
!   x2,y2,z2         Coordinates of second point
!   ind              The rotation matrix to use
!   MAXROT           The maximum number of rotation matrices dimensioned
!   rotmat           The rotation matrices
!
!
!
! OUTPUT VARIABLES:
!
!   sqdist           The squared distance accounting for the anisotropy
!                      and the rotation of coordinates (if any).
!
!
! NO EXTERNAL REFERENCES
!
!
!-----------------------------------------------------------------------
      implicit integer(i-n), real(a-h, o-z)                    !jd
      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
!
! Compute component distance vectors and the squared distance:
!
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist_new = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx      &
                   + rotmat(ind,i,2) * dy      &
                   + rotmat(ind,i,3) * dz
            sqdist_new = sqdist_new + cont * cont
      end do
      return
      end




      subroutine sortem_new(ib,ie,a,iperm,b,c,d,e,f,g,h)
!-----------------------------------------------------------------------
!
!                      Quickersort Subroutine
!                      **********************
!
! This is a subroutine for sorting a real array in ascending order. This
! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
! in collected algorithms of the ACM.
!
! The method used is that of continually splitting the array into parts
! such that all elements of one part are less than all elements of the
! other, with a third part in the middle consisting of one element.  An
! element with value t is chosen arbitrarily (here we choose the middle
! element). i and j give the lower and upper limits of the segment being
! split.  After the split a value q will have been found such that
! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
! performs operations on the two segments (i,q-1) and (q+1,j) as follows
! The smaller segment is split and the position of the larger segment is
! stored in the lt and ut arrays.  If the segment to be split contains
! two or fewer elements, it is sorted and another segment is obtained
! from the lt and ut arrays.  When no more segments remain, the array
! is completely sorted.
!
!
! INPUT PARAMETERS:
!
!   ib,ie        start and end index of the array to be sorteda
!   a            array, a portion of which has to be sorted.
!   iperm        0 no other array is permuted.
!                1 array b is permuted according to array a
!                2 arrays b,c are permuted.
!                3 arrays b,c,d are permuted.
!                4 arrays b,c,d,e are permuted.
!                5 arrays b,c,d,e,f are permuted.
!                6 arrays b,c,d,e,f,g are permuted.
!                7 arrays b,c,d,e,f,g,h are permuted.
!               >7 no other array is permuted.
!
!   b,c,d,e,f,g,h  arrays to be permuted according to array a.
!
! OUTPUT PARAMETERS:
!
!    a      = the array, a portion of which has been sorted.
!
!    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
!
! NO EXTERNAL ROUTINES REQUIRED:
!
!-----------------------------------------------------------------------
      implicit integer(i-n), real(a-h, o-z)                    !jd

      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
      integer   lt(64),ut(64),i,j,k,m,p,q
!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
!
! If this segment has more than two elements  we split it
!
 10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
!
! Start at the beginning of the segment, search for k such that a(k)>t
!
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
!
! Update q and search for another pair to interchange:
!
      q = q-1
      go to 20
 50   q = k-1
 60   continue
!
! The upwards search has now met the downwards search:
!
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
!
! Update m and split the new smaller segment
!
 80   m = m+1
      go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
!
! If lt and ut contain more segments to be sorted repeat process:
!
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end
