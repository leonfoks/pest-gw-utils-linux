!     Last change:  JD    9 May 2003   11:40 am

program ppk2fac

! -- Program PPK2FAC1 performs the first step in kriging from a set of
!    pilot points to a two-dimensional array; it builds the factors by which
!    the value of each point contributes to the value at the grid cell.

	use defn
	use inter

	implicit none

        integer     :: ifail,idate,iheader,i,ncol,nrow,icol,irow,ierr,structunit, &
                       numint,currint,newint,ihuge,itemp,numstruct,numvario,outunit, &
                       nbb,j,izone,ndat,ngrid,istruct,k_ktype,n_nst,regunit,itrans
        real        :: s_skmean,c_c0,pmx
        real        :: blankrad
        double precision :: minsep,minsep2,eastdiff,northdiff,distance
        integer     :: iloc(1)
        integer     :: i_it(MAX_STRUCT_VARIO)
        real        :: c_cc(MAX_STRUCT_VARIO),a_ang(MAX_STRUCT_VARIO),  &
                       a_aa(MAX_STRUCT_VARIO),a_anis(MAX_STRUCT_VARIO)
        integer, allocatable, dimension(:)             :: icellno
        integer, allocatable, dimension(:)             :: inumdat
        integer, allocatable, dimension(:)             :: minpt,maxpt
	integer, allocatable, dimension(:,:)	       :: intarray
        real, allocatable, dimension(:)                :: radmax,radmin,blankmax
        real, allocatable, dimension(:)                :: eastgrid,northgrid
        real,allocatable, dimension(:)                 :: eastdat,northdat,valdat
        real, allocatable, dimension(:)                :: totcov
        real, allocatable, dimension(:,:)              :: regarray
        real,allocatable,dimension(:,:)                :: east,north,variance
        character (len=1)                              :: alogtrans,aoutfile,arealformat
        character (len=1),  allocatable, dimension(:)  :: akrig
        character (len=10), allocatable, dimension(:)  :: astructure
        character (len=15)                             :: anum1,atemp
	character (len=200)                            :: aprompt,intfile,structfile, &
                                                          atempf,outfile,arealfile,regfile

	type (modelgrid) gridspec
        type (geostructure), allocatable, dimension(:) :: structure(:)
        type (variogram), allocatable, dimension(:)    :: vario(:)


	integer, parameter :: MAXINT=100
	integer		   :: intval(MAXINT),iwork(MAXINT)


	write(amessage,5)
5	format(' Program PPK2FAC1 calculates point-to-cell factors by which kriging ',&
	'is undertaken from a set of pilot points to the finite-difference grid.')
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
	allocate(intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run PPK2FAC.')
	  go to 9890
	end if

30      call read_pilot_points_file(ifail, &
	' Enter name of pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(intarray,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
          call free_grid_mem(gridspec)
	  write(6,*)
	  go to 10
	end if

        write(6,*)
61      write(6,63,advance='no')
63      format(' Enter minimum allowable points separation: ')
        if(key_read(minsep).ne.0) go to 61
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 30
        end if
        imessage=0
        if(minsep.lt.0.0d0)minsep=0.0d0
!        if(minsep.gt.0.0d0)then
          minsep2=minsep*minsep
          if(num_pilot_points.gt.1)then
          do i=1,num_pilot_points-1
            do j=i+1,num_pilot_points
              eastdiff=pilot_point_east(i)-pilot_point_east(j)
              northdiff=pilot_point_north(i)-pilot_point_north(j)
              distance=eastdiff*eastdiff+northdiff*northdiff
!              if((distance.le.minsep2).and. &
              if((distance.lt.minsep2).and. &
                 (pilot_point_zone(i).eq.pilot_point_zone(j)))then
                imessage=imessage+1
                if(imessage.gt.20) go to 9900
                if(imessage.eq.1)then
                  write(amessage,65)
65                format(' The following points are separated by less than the ', &
                  'minimum separation ---->')
                  call write_message(leadspace='yes')
                end if
                write(amessage,66) trim(pilot_point_id(i)), &
                trim(pilot_point_id(j)),sqrt(distance)
66              format(1x,a,t15,a, t30, '(separation = ',f12.3,')')
                call write_message()
              end if
            end do
          end do
          end if
!        end if
        if(imessage.gt.0) go to 9900

        write(6,*)
70	continue
	aprompt=' Enter name of zonal integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 61
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
	      ' different integers for proper PPK2FAC execution. The array ', &
	      ' contained in file ',a,' holds more than this. To increase ',&
	      'this limit, edit PPK2FAC source code, increase parameter ', &
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

        allocate(astructure(numint),radmax(numint),radmin(numint), &
        akrig(numint),minpt(numint),maxpt(numint),blankmax(numint),stat=ierr)  !de-allocate below for "e"
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

        radmax=0.0
        radmin=0.0
        blankmax=0.0
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
280     format('    Enter structure name (blank if no interpolation for this zone): ')
        read(5,'(a)') atemp
        if(index(eschar,atemp(1:2)).ne.0)then
	  if(i.gt.2) then
	    i=i-2
	    go to 255
	  else if(i.eq.2) then
	    go to 200
	  else
            deallocate(astructure,radmax,radmin,blankmax,akrig,minpt,maxpt,stat=ierr)
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
290     format('    Perform simple or ordinary kriging [s/o]: ')
	read(5,'(a)') akrig(i)
	if(akrig(i).eq.' ') go to 289
	if(index(eschar,akrig(i)).ne.0) then
	  write(6,*)
	  go to 259
	end if
	call casetrans(akrig(i),'lo')
	if((akrig(i).ne.'s').and.(akrig(i).ne.'o')) go to 289
299     write(6,300,advance='no')
300     format('    Enter search radius: ')
        itemp=key_read(radmax(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 289
        end if
	if(itemp.ne.0) then
          write(6,310)
310	  format('    Data input error  - try again.')
	  go to 299
	end if
        if(radmax(i).le.0.0d0)then
          write(6,311)
311       format('    Must be greater than zero - try again.')
          go to 299
        end if
        if(radmax(i).gt.1.0e15)radmax(i)=1.0e15

2990    write(6,3000,advance='no')
3000    format('    Enter blanking radius: ')
        itemp=key_read(blankmax(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 299
        end if
	if(itemp.ne.0) then
          write(6,310)
	  go to 2990
	end if
        if(blankmax(i).le.0.0d0)then
          write(6,311)
          go to 2990
        end if
        if(blankmax(i).gt.1.0e15)blankmax(i)=1.0e15

329     write(6,330,advance='no')
330     format('    Enter minimum number of pilot points to use for interpolation: ')
        itemp=key_read(minpt(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 2990
        end if
	if(itemp.ne.0) then
          write(6,310)
	  go to 329
	end if
        if(minpt(i).lt.1)then
          write(6,332)
332       format('    Must be greater than 1 - try again.')
          go to 329
        end if
335     write(6,337,advance='no')
337     format('    Enter maximum number of pilot points to use for interpolation: ')
        itemp=key_read(maxpt(i))
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
          go to 329
        end if
	if(itemp.ne.0) then
          write(6,310)
	  go to 335
	end if
        if(maxpt(i).lt.minpt(i))then
          write(6,340)
340       format('    Must be greater than min. no. of points - try again.')
          go to 335
        end if
        if(maxpt(i).gt.500)then
          write(6,341)
341       format('    Must not be greater than 500 - try again.')
          go to 335
        end if
        go to 255

350	continue

        do i=1,numint
          if(astructure(i).ne.' ')go to 351
        end do
        write(amessage,352)
352     format(' Interpolation must take place for at least one zone in integer ', &
        'array - try again.')
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
          open(unit=outunit,file=outfile,status='new',iostat=ierr)
          if(ierr.ne.0) then
            open(unit=outunit,file=outfile,status='replace',iostat=ierr)
          end if
        else
          open(unit=outunit,file=outfile,status='new',form='binary',iostat=ierr)
          if(ierr.ne.0) then
            open(unit=outunit,file=outfile,status='replace',form='binary',iostat=ierr)
          end if
        end if
        if(ierr.ne.0) then
          write(amessage,358) trim(outfile)
358       format(' Unable to open file ',a,' for output - try again.')
          go to 353
        end if

	allocate(variance(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

360     aprompt=' Enter name for output standard deviation array file: '
        call write_real_array(ifail,aprompt,variance,pm_header=headerspec, &
        rows=nrow,columns=ncol,istage=1,realfile=arealfile,aaformat=arealformat)
        if(escset.eq.1) then
          escset=0
          close(unit=outunit)
          deallocate(variance,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
	    go to 9890
	  end if
          go to 353
        end if


        aprompt = ' Enter name for regularisation information file: '
	call open_output_file(ifail,aprompt,regfile,regunit)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1)then
          escset=0
          write(6,*)
          go to 360
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

        allocate(eastdat(num_pilot_points),northdat(num_pilot_points),  &
        valdat(num_pilot_points),inumdat(num_pilot_points),  &
        regarray(num_pilot_points,num_pilot_points),totcov(num_pilot_points),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
	allocate(east(ncol,nrow),north(ncol,nrow),  &
        eastgrid(ncol*nrow),northgrid(ncol*nrow),icellno(ncol*nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
        variance=1.1e35
        regarray=-1.1e35           !regarray is an array
        totcov=-1.1e35             !totcov is an array

        call rel_centre_coords(east,north,gridspec)
        call grid2earth(east,north,gridspec)

! -- The coordinates of the pilot points have the grid left corner easting and
!    northing subtracted.

        do i=1,num_pilot_points
          pilot_point_east(i)=pilot_point_east(i)-gridspec%east_corner
          pilot_point_north(i)=pilot_point_north(i)-gridspec%north_corner
        end do

! -- Kriging weights are now evaluated zone by zone for cell centres in
!    the finite-difference grid using GSLIB subroutines.

        call addquote(pilot_points_file,atempf)
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        call addquote(intfile,atempf)
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        if(aoutfile.eq.'f')then
          write(outunit,*) ncol,nrow
        else
          write(outunit) ncol,nrow
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
          do irow=1,nrow
            do icol=1,ncol
              if(intarray(icol,irow).eq.itemp)then
                ngrid=ngrid+1
                eastgrid(ngrid)=east(icol,irow)
                northgrid(ngrid)=north(icol,irow)
                call rc2cell(icellno(ngrid),irow,icol,gridspec)
              end if
            end do
          end do
          if(ngrid.eq.0) go to 410
          ndat=0
          do i=1,num_pilot_points
            if(pilot_point_zone(i).eq.itemp)then
              ndat=ndat+1
              eastdat(ndat)=pilot_point_east(i)
              northdat(ndat)=pilot_point_north(i)
              valdat(ndat)=pilot_point_val(i)
              inumdat(ndat)=i
            end if
          end do
          if(ndat.eq.0) then
            call num2char(intval(izone),anum1)
            write(6,399) trim(anum1)
399         format(/,' Warning: no pilot points assigned to integer array zone ',a)
            go to 410
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
            a_ang(i)=vario(j)%angle
            a_aa(i)=vario(j)%a
            a_anis(i)=vario(j)%anis
          end do

          call readparm(minpt(izone),maxpt(izone),radmax(izone),k_ktype,s_skmean, &
          n_nst,c_c0,i_it,c_cc,a_ang,a_aa,a_anis,ndat,eastdat,northdat,valdat)
          blankrad=blankmax(izone)
          call kb2d(ngrid,ndat,ncol,nrow,inumdat,icellno,eastgrid,northgrid, &
          outfile,aoutfile,outunit,variance,pmx,itrans,num_pilot_points,regarray,totcov,  &
          blankrad)

410     continue
        end do

        write(6,*)
        write(6,450) trim(outfile)
450     format('  - kriging factors written to file ',a)
        call write_real_array(ifail,aprompt,variance,pm_header=headerspec, &
        rows=nrow,columns=ncol,istage=2,realfile=arealfile,aaformat=arealformat)

        write(regunit,'(i5)') num_pilot_points
        do i=1,num_pilot_points
          write(regunit,451) trim(pilot_point_id(i)), &
          pilot_point_east(i)+gridspec%east_corner, &
          pilot_point_north(i)+gridspec%north_corner
451       format(a,2x,f15.3,2x,f15.3)
        end do
        do i=1,num_pilot_points
          write(regunit,460) (regarray(j,i),j=1,num_pilot_points),totcov(i)
460       format(7(1x,1pg14.7))
        end do
        close(unit=regunit)
        write(6,470) trim(regfile)
470     format('  - regularisation information written to file ',a)


        go to 9900


9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
	call free_grid_mem(gridspec)
        call free_point_mem()
        deallocate(intarray,east,north,variance,eastgrid,northgrid,icellno,stat=ierr)
        deallocate(eastdat,northdat,valdat,inumdat,stat=ierr)
        deallocate(astructure,radmax,radmin,blankmax,akrig,minpt,maxpt,stat=ierr)
        deallocate(structure,vario,stat=ierr)
        deallocate(regarray,stat=ierr)
        deallocate(totcov,stat=ierr)

end program ppk2fac




      subroutine readparm(n_ndmin,n_ndmax,r_radius,k_ktype,s_skmean, &
      n_nst,c_c0,i_it,c_cc,a_ang,a_aa,a_anis,n_nd,eastdat,northdat,valdat)


      implicit integer(i-n), real(a-h, o-z)                    !jd

!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
!-----------------------------------------------------------------------
      include  'kb2d.inc'
      parameter(MV=20)
      real      var(MV)
      character datafl*40,outfl*40,dbgfl*40,str*40
      logical   testfl

! -- Dimension statements added by myself for subroutine arguments.

      integer i_it(n_nst)
      real c_cc(n_nst),a_ang(n_nst),a_aa(n_nst),a_anis(n_nst), &
      valdat(n_nd)
      real eastdat(n_nd),northdat(n_nd)
!
! Unit numbers:
!
      lout = 2
      ldbg = 3

!
! Read Input Parameters:
!

! -- Parameters that should have been read from input file.

      datafl=' '
      ixl=-9999
      iyl=-9999
      ivrl=-9999
      tmin=-1.1e35
      tmax=1.1e35
      idbg=0                    !debugging level
      dbgfl='debug1.dat'
      outfl='kb2d.out'
      nxdis=1
      nydis=1


! -- Here are the dummy grid specifications.

      nx=1
      xmn=0.0
      xsiz=1.0
      ny=1
      ymn=0.0
      ysiz=1.0


! -- The following variables were supplied through the subroutine argument.

      ndmin=n_ndmin
      ndmax=n_ndmax
      radius=r_radius
      ktype=k_ktype
      skmean=s_skmean
      nst=n_nst
      c0=c_c0

            do i=1,nst
                  it(i)=i_it(i)
                  cc(i)=c_cc(i)
                  ang(i)=a_ang(i)
                  aa(i)=a_aa(i)
!                  a2=a_a2(i)
!                  anis(i) = a2 / aa(i)
                  anis(i)=a_anis(i)
                  if(it(i).eq.4) then
                        if(aa(i).lt.0.0) stop ' INVALID power variogram'
                        if(aa(i).gt.2.0) stop ' INVALID power variogram'
                  end if
            end do

!      write(*,*)
      if(nst.gt.MAXNST)   stop ' nst is too big - contact programmer'
      if(ndmax.gt.MAXSAM) stop ' ndmax is too big - contact programmer'

      av = 0.0
      ss = 0.0
      nd = n_nd

      if(nd.gt.MAXDAT) then
            write(*,*) ' nd too big - contact programmer'
            stop
      end if
      do i=1,nd
        x(i)  = eastdat(i)
        y(i)  = northdat(i)
        vr(i) = valdat(i)
        av     = av + vr(i)
        ss     = ss + vr(i)*vr(i)
      end do

!
! Open the output files:
!
!      open(lout,file=outfl,status='UNKNOWN')
!      write(lout,300)
! 300  format('KB2D Output',/,'2',/,'Estimate',/,'Estimation Variance')
      open(ldbg,file=dbgfl,status='UNKNOWN')
!
! Compute the averages and variances as an error check for the user:
!
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
!
! Write Some of the Statistics to the screen:
!

       write(6,*)
       write(6,900) nd
900    format('   Number of pilot points for this zone     = ',i5)
       write(6,901) av
901    format('   Mean data value for these pilot points   = ',1pg12.5)
       write(6,902) sqrt(max(ss,0.0))
902    format('   Data standard deviation for these points = ',1pg12.5)
       write(6,903)
903    format('   Working....')
!      write(*,900) nd,av,sqrt(max(ss,0.0))
! 900  format(/' There are ',i8,' data with:',/,     &
!              '   mean value          = ',f12.5,/,  &
!              '   standard deviation  = ',f12.5,/)
      return
      end



      subroutine kb2d(n_npts,n_ndat,n_ncol,n_nrow,inumdat,icellno,epoint,npoint,outfile, &
      aoutfile,outunit,var_arr,pmx,itrans,npp,regarray,totcov,blankrad)

      implicit integer(i-n), real(a-h, o-z)                    !jd

!-----------------------------------------------------------------------
!
!           Ordinary/Simple Kriging of a 2-D Rectangular Grid
!           *************************************************
!
! This subroutine estimates point or block values of one variable by
! ordinary kriging.  All of the samples are rescanned for each block
! estimate; this makes the program simple but inefficient.  The data
! should NOT contain any missing values.  Unestimated points are
! returned as -1.0e21
!
!
!
! Original:  A.G. Journel                                           1978
! Revisions: B.E. Buxton                                       Apr. 1983
!-----------------------------------------------------------------------
      include  'kb2d.inc'
      real      xdb(MAXDIS),ydb(MAXDIS),xa(MAXSAM),ya(MAXSAM), &
                vra(MAXSAM),dist(MAXSAM)
      real*8    r(MAXSAM+1),rr(MAXSAM+1),s(MAXSAM+1),a(MAXKRG)
      real*8    dtemp                                               !jd
      integer   nums(MAXSAM)
      logical   first

! -- The following dimension statement was added by myself for the
!    subroutine arguments.

     integer n_npts,n_ndat,outunit,n_ncol,n_nrow,itrans             !jd
     integer icellno(n_npts),inumdat(n_ndat)                        !jd
     real pmx                                                       !jd
     real epoint(n_npts),npoint(n_npts)                             !jd
     real var_arr(n_ncol,n_nrow)                                    !jd
     real regarray(npp,npp)                                         !jd
     real totcov(npp)                                               !jd
     real blankrad                                                  !jd
     character*(*) outfile,aoutfile                                 !jd

     real dist_to_nearest                                           !jd

!jd      data      first/.true./,PMX/9999.0/

      first=.true.                                    !jd
!
! Echo the input parameters if debugging flag is >2:
!
      if(idbg.gt.2) then
            write(ldbg,*) 'KB2D Parameters'
            write(ldbg,*)
            write(ldbg,*) 'Variogram Parameters for ',nst,' structures:'
            write(ldbg,*) '  Nugget effect:         ',c0
            write(ldbg,*) '  Types of variograms:   ',(it(i),i=1,nst)
            write(ldbg,*) '  Contribution cc        ',(cc(i),i=1,nst)
            write(ldbg,*) '  Ranges:                ',(aa(i),i=1,nst)
            write(ldbg,*) '  Angle for Continuity:  ',(ang(i),i=1,nst)
            write(ldbg,*) '  Anisotropy Factors:    ',(anis(i),i=1,nst)
            write(ldbg,*) ' '
            write(ldbg,*) 'Grid for Kriging:'
            write(ldbg,*) '  Number of X and Y Blocks:',nx,ny
            write(ldbg,*) '  Origin of X and Y Blocks:',xmn,ymn
            write(ldbg,*) '  Size   of X and Y Blocks:',xsiz,ysiz
            write(ldbg,*) ' '
            write(ldbg,*) 'Discretization of blocks:  ',nxdis,nydis
            write(ldbg,*) 'Search Radius:             ',radius
            write(ldbg,*) 'Minimum number of samples: ',ndmin
            write(ldbg,*) 'Maximum number of samples: ',ndmax
            write(ldbg,*) ' '
      endif
!
! Echo the input data if debugging flag >1:
!
      if(idbg.ge.4) then
            do id=1,nd
                  write(ldbg,99) id,x(id),y(id),vr(id)
 99               format('Data: ',i5,' at ',2f12.3,' value: ',f12.5)
            end do
      endif
!
! Set up the discretization points per block.  Figure out how many
! are needed, the spacing, and fill the xdb and ydb arrays with the
! offsets relative to the block center (this only gets done once):
!
      ndb  = nxdis * nydis
      if(ndb.gt.MAXDIS) then
            write(*,*) ' ERROR KB2D: Too many discretization points '
            write(*,*) '             Increase MAXDIS or lower n[xy]dis'
            stop
      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      xloc = -0.5*(xsiz+xdis)
      i    = 0
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  i = i+1
                  xdb(i) = xloc
                  ydb(i) = yloc
            end do
      end do
!
! Initialize accumulators:
!
      cbb  = 0.0
      rad2 = radius*radius
!
! Calculate Block Covariance. Check for point kriging.
!
      cov   = cova2(xdb(1),ydb(1),xdb(1),ydb(1),nst,c0,PMX,cc, &
                    aa,it,ang,anis,first,passmaxcov)
!
! Keep this value to use for the unbiasedness constraint:
!
      unbias = cov
      first  = .false.
      if (ndb.le.1) then
            cbb = cov
      else
            do i=1,ndb
                  do j=1,ndb
                        cov = cova2(xdb(i),ydb(i),xdb(j),ydb(j),nst,c0, &
                                    PMX,cc,aa,it,ang,anis,first,passmaxcov)
                        if(i.eq.j) cov = cov - c0
                        cbb = cbb + cov
                  end do
            end do
            cbb = cbb/real(ndb*ndb)
      endif
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      endif
!
! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!
      nk = 0
      ak = 0.0
      vk = 0.0
!jd      do 4 iy=1,ny                                               !jd
!jd      yloc = ymn + (iy-1)*ysiz                                   !jd
!jd      do 4 ix=1,nx                                               !jd
!jd            xloc = xmn + (ix-1)*xsiz                             !jd

         do 4 i_ipts=1,n_npts                                       !jd
            xloc=epoint(i_ipts)                                     !jd
            yloc=npoint(i_ipts)                                     !jd
            icell=icellno(i_ipts)                                   !jd
!
! Find the nearest samples within each octant: First initialize
! the counter arrays:
!
            na = 0
            do isam=1,ndmax
                  dist(isam) = 1.0e+20
                  nums(isam) = 0
            end do
!
! Scan all the samples (this is inefficient and the user with lots of
! data should move to ktb3d):
!
            dist_to_nearest=1.0e30
            do 6 id=1,nd
                  dx = x(id) - xloc
                  dy = y(id) - yloc
                  h2 = dx*dx + dy*dy
                  if(h2.lt.dist_to_nearest)dist_to_nearest=h2
                  if(h2.gt.rad2) go to 6
!
! Do not consider this sample if there are enough close ones:
!
                  if(na.eq.ndmax.and.h2.gt.dist(na)) go to 6
!
! Consider this sample (it will be added in the correct location):
!
                  if(na.lt.ndmax) na = na + 1
                  nums(na)           = id
                  dist(na)           = h2
                  if(na.eq.1) go to 6
!
! Sort samples found thus far in increasing order of distance:
!
                  n1 = na-1
                  do ii=1,n1
                        k=ii
                        if(h2.lt.dist(ii)) then
                              jk = 0
                              do jj=k,n1
                                    j  = n1-jk
                                    jk = jk+1
                                    j1 = j+1
                                    dist(j1) = dist(j)
                                    nums(j1) = nums(j)
                              end do
                              dist(k) = h2
                              nums(k) = id
                              go to 6
                        endif
                  end do
 6          continue
!
! Is there enough samples?
!
            if(dist_to_nearest.gt.blankrad*blankrad) go to 1                             !jd
            if(na.lt.ndmin) then
                  if(idbg.ge.2) &
                  write(ldbg,*) 'Block ',ix,iy, 'not estimated'
                  est  = UNEST
                  estv = UNEST
                  go to 1
            endif
!
! Put coordinates and values of neighborhood samples into xa,ya,vra:
!
            do ia=1,na
                  jj      = nums(ia)
                  xa(ia)  = x(jj)
                  ya(ia)  = y(jj)
                  vra(ia) = vr(jj)
            end do
!
! Handle the situation of only one sample:
!
            if(na.eq.1) then
                  cb1 = cova2(xa(1),ya(1),xa(1),ya(1),nst,c0,   &
                              PMX,cc,aa,it,ang,anis,first,passmaxcov)
                  xx  = xa(1) - xloc
                  yy  = ya(1) - yloc
!
! Establish Right Hand Side Covariance:
!
                  if(ndb.le.1) then
                        cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0,  &
                                   PMX,cc,aa,it,ang,anis,first,passmaxcov)
                  else
                        cb  = 0.0
                        do i=1,ndb
                              cb = cb + cova2(xx,yy,xdb(i),ydb(i),nst, &
                                        c0,PMX,cc,aa,it,ang,anis,first,passmaxcov)
                              dx = xx - xdb(i)
                              dy = yy - ydb(i)
                              if((dx*dx+dy*dy).lt.EPSLON) &
                              cb = cb - c0
                        end do
                        cb = cb / real(ndb)
                  end if
                  if(ktype.eq.0) then
                        s(1) = cb/cbb
                        est  = s(1)*vra(1) + (1.0-s(1))*skmean
                        estv = cbb - s(1) * cb
                  else
                        s(1) = 1.0                                      !jd
                        est  = vra(1)
                        estv = cbb - 2.0*cb + cb1
                  end if

                  if(ktype.eq.0)then                                    !jd
                    rrtemp=(1.0-s(1))*skmean                            !jd
                  else                                                  !jd
                    rrtemp=0.0                                          !jd
                  end if                                                !jd
                  if(aoutfile.eq.'f')then                               !jd
                    write(outunit,*) icell,itrans,1,rrtemp, &           !jd
                    inumdat(nums(1)),s(1)                                 !jd
                  else                                                    !jd
                    write(outunit)   icell,itrans,1,rrtemp, &             !jd
                    inumdat(nums(1)),s(1)                                 !jd
                  end if                                                  !jd

            else
!
! Solve the Kriging System with more than one sample:
!
                  neq = na + ktype
                  nn  = (neq + 1)*neq/2
!
! Set up kriging matrices:
!
                  in=0
                  do j=1,na
!
! Establish Left Hand Side Covariance Matrix:
!
                        do i=1,j
                              in = in + 1
                              a(in) = dble( cova2(xa(i),ya(i),xa(j),  &
                                            ya(j),nst,c0,PMX,cc,aa,it, &
                                            ang,anis,first,passmaxcov) )
                              iitmp1=inumdat(nums(i))                           !jd
                              iitmp2=inumdat(nums(j))                           !jd
                              if(regarray(iitmp1,iitmp2).lt.-1.0e35)then        !jd
                                regarray(iitmp1,iitmp2)=passmaxcov -a(in)       !jd
                                regarray(iitmp2,iitmp1)=regarray(iitmp1,iitmp2) !jd
                                totcov(iitmp1)=passmaxcov                       !jd
                                totcov(iitmp2)=passmaxcov                       !jd
                              end if
                        end do
                        xx = xa(j) - xloc
                        yy = ya(j) - yloc
!
! Establish Right Hand Side Covariance:
!
                        if(ndb.le.1) then
                              cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0,  &
                                         PMX,cc,aa,it,ang,anis,first,passmaxcov)
                        else
                              cb  = 0.0
                              do j1=1,ndb
                                    cb = cb + cova2(xx,yy,xdb(j1),  &
                                         ydb(j1),nst,c0,PMX,cc,aa,  &
                                         it,ang,anis,first,passmaxcov)
                                    dx = xx - xdb(j1)
                                    dy = yy - ydb(j1)
                                    if((dx*dx+dy*dy).lt.EPSLON)  &
                                          cb = cb - c0
                              end do
                              cb = cb / real(ndb)
                        end if
                        r(j)  = dble(cb)
                        rr(j) = r(j)
                  end do
!
! Set the unbiasedness constraint:
!
                  if(ktype.eq.1) then
                        do i=1,na
                              in    = in + 1
                              a(in) = dble(unbias)
                        end do
                        in      = in + 1
                        a(in)   = 0.0
                        r(neq)  = dble(unbias)
                        rr(neq) = r(neq)
                  end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
                  if(idbg.ge.3) then
                        write(ldbg,101) ix,iy
                        is = 1
                        do i=1,neq
                              ie = is + i - 1
                              write(ldbg,102) i,r(i),(a(j),j=is,ie)
                              is = is + i
                        end do
 101                    format(/,'Kriging Matrices for Node: ',2i4,  &
                                 ' RHS first')
 102                    format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
                  endif
!
! Solve the Kriging System:
!
                  call ksol(1,neq,1,a,r,s,ising)
!
! Write a warning if the matrix is singular:
!
                  if(ising.ne.0) then
        write(6,450) icell                                            !jd
450     format(' WARNING: singular kriging matrix for cell',i5)      !jd
                        est  = UNEST
                        estv = UNEST
                        go to 1
                  endif
!
! Write the kriging weights and data if requested:
!
                  if(idbg.ge.2) then
                        write(ldbg,*) '       '
                        write(ldbg,*) 'BLOCK: ',ix,iy
                        write(ldbg,*) '       '
                        if(ktype.eq.1) write(ldbg,*)    &
                        '  Lagrange multiplier: ',s(neq)*unbias
                        write(ldbg,*) '  BLOCK EST: x,y,vr,wt '
                        do i=1,na
                        write(ldbg,'(4f12.3)') xa(i),ya(i),vra(i),s(i)
                        end do
                  endif
!
! Compute the estimate and the kriging variance:
!
                  est  = 0.0
                  estv = cbb
                  sumw = 0.0
                  if(ktype.eq.1) estv = estv - real(s(na+1))
                  do i=1,na
                        sumw = sumw + real(s(i))
                        est  = est  + real(s(i))*vra(i)
                        estv = estv - real(s(i)*rr(i))
                  end do
                  if(ktype.eq.0) est = est + (1.0-sumw)*skmean

                  if(ktype.eq.0)then                                    !jd
                    rrtemp=(1.0-sumw)*skmean                            !jd
                  else                                                  !jd
                    rrtemp=0.0                                          !jd
                  end if                                                !jd
                  if(aoutfile.eq.'f')then                               !jd
                    write(outunit,*) icell,itrans,na,rrtemp, &          !jd
                    ((inumdat(nums(i)),real(s(i))),i=1,na)                 !jd
                  else                                                     !jd
                    write(outunit)   icell,itrans,na,rrtemp, &             !jd
                    ((inumdat(nums(i)),real(s(i))),i=1,na)                 !jd
                  end if                                                   !jd

            endif
            if(idbg.ge.2) then
                  write(ldbg,*) '  est  ',est
                  write(ldbg,*) '  estv ',estv
                  write(ldbg,*) ' '
            endif
!
! Write the result to the output file:
!
! 1          write(lout,'(f8.3,1x,f8.3)') est,estv

1           continue

            if(est.gt.UNEST)then
              iirow=(icell-1)/n_ncol+1
              iicol=icell-((iirow-1)*n_ncol)
              rrtemp=estv
              if(rrtemp.lt.0.0)rrtemp=0.0
              var_arr(iicol,iirow)=sqrt(rrtemp)
            end if

            if(est.gt.UNEST) then
                  nk = nk + 1
                  ak = ak + est
                  vk = vk + est*est
            end if
!
! END OF MAIN LOOP OVER ALL THE BLOCKS:
!
 4    continue
      if(nk.ge.1) then
            ak = ak / real(nk)
            vk = vk/real(nk) - ak*ak
             write(6,105) nk
105          format('   No. of grid points for which factors were calculated = ',i5)
!             write(6,106) ak
!106          format('   Average interpolated value at these points           = ',1pg12.5)
!             if(vk.lt.0.0)vk=0.0
!             write(6,107) sqrt(vk)
!107          format('   Standard deviation of interpolated grid point values = ',1pg12.5)

!            write(ldbg,105) nk,ak,vk
!            write(*,   105) nk,ak,vk
! 105        format(/,'Estimated   ',i8,' blocks ',/,  &
!                     '  average   ',f9.4,/,'  variance  ',f9.4,/)
      else
        write(6,105) nk
      end if

! -- We add a little code to compute the diagonal elements of the regularisatoin
!    array and totcov for any pertinent blocks for which this has not been done.


      do i=1,n_ndat                                       !jd
        iitmp1=inumdat(i)                                !jd
        if(totcov(iitmp1).lt.-1.0e35)then                 !jd
          dtemp = dble( cova2(x(i),y(i),x(i),  &          !jd
                  y(i),nst,c0,PMX,cc,aa,it, &             !jd
                  ang,anis,first,passmaxcov))             !jd
          regarray(iitmp1,iitmp1)=passmaxcov-dtemp        !jd
          totcov(iitmp1)=passmaxcov                       !jd
        end if                                            !jd
      end do                                              !jd

      return
      end



      real function cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                          ang,anis,first,passmaxcov)              !jd

      implicit integer(i-n), real(a-h, o-z)                    !jd

!-----------------------------------------------------------------------
!
!              Covariance Between Two Points (2-D Version)
!              *******************************************
!
! This function returns the covariance associated with a variogram model
! that is specified by a nugget effect and possibly four different
! nested varigoram structures.  The anisotropy definition can be
! different for each of the nested structures (spherical, exponential,
! gaussian, or power).
!
!
!
! INPUT VARIABLES:
!
!   x1,y1            Coordinates of first point
!   x2,y2            Coordinates of second point
!   nst              Number of nested structures (max. 4).
!   c0               Nugget constant (isotropic).
!   PMX              Maximum variogram value needed for kriging when
!                      using power model.  A unique value of PMX is
!                      used for all nested structures which use the
!                      power model.  therefore, PMX should be chosen
!                      large enough to account for the largest single
!                      structure which uses the power model.
!   cc(nst)          Multiplicative factor of each nested structure.
!   aa(nst)          Parameter "a" of each nested structure.
!   it(nst)          Type of each nested structure:
!                      1. spherical model of range a;
!                      2. exponential model of parameter a;
!                           i.e. practical range is 3a
!                      3. gaussian model of parameter a;
!                           i.e. practical range is a*sqrt(3)
!                      4. power model of power a (a must be gt. 0  and
!                           lt. 2).  if linear model, a=1,c=slope.
!   ang(nst)         Azimuth angle for the principal direction of
!                      continuity (measured clockwise in degrees from Y)
!   anis(nst)        Anisotropy (radius in minor direction at 90 degrees
!                      from "ang" divided by the principal radius in
!                      direction "ang")
!   first            A logical variable which is set to true if the
!                      direction specifications have changed - causes
!                      the rotation matrices to be recomputed.
!
!
!
! OUTPUT VARIABLES: returns "cova2" the covariance obtained from the
!                   variogram model.
!
!
!
!-----------------------------------------------------------------------
      parameter(DTOR=3.14159265/180.0,EPSLON=0.0000001)
      real      aa(*),cc(*),ang(*),anis(*),rotmat(4,4),maxcov
      integer   it(*)
      logical   first
      save      rotmat,maxcov
!
! The first time around, re-initialize the cosine matrix for the
! variogram structures:
!
      if(first) then
            maxcov = c0
            do is=1,nst
                  azmuth       = (90.0-ang(is))*DTOR
                  rotmat(1,is) =  cos(azmuth)
                  rotmat(2,is) =  sin(azmuth)
                  rotmat(3,is) = -sin(azmuth)
                  rotmat(4,is) =  cos(azmuth)
                  if(it(is).eq.4) then
                        maxcov = maxcov + PMX
                  else
                        maxcov = maxcov + cc(is)
                  endif
            end do
      endif
      passmaxcov=maxcov                                            !jd
!
! Check for very small distance:
!
      dx = x2-x1
      dy = y2-y1
      if((dx*dx+dy*dy).lt.EPSLON) then
            cova2 = maxcov
            return
      endif
!
! Non-zero distance, loop over all the structures:
!
      cova2 = 0.0
      do is=1,nst
!
! Compute the appropriate structural distance:
!
            dx1 = (dx*rotmat(1,is) + dy*rotmat(2,is))
            dy1 = (dx*rotmat(3,is) + dy*rotmat(4,is))/anis(is)
            h   = sqrt(max((dx1*dx1+dy1*dy1),0.0))
            if(it(is).eq.1) then
!
! Spherical model:
!
                  hr = h/aa(is)
                  if(hr.lt.1.0) cova2 = cova2    &
                                      + cc(is)*(1.-hr*(1.5-.5*hr*hr))
            else if(it(is).eq.2) then
!
! Exponential model:
!
                  cova2 = cova2 +cc(is)*exp(-h/aa(is))
            else if(it(is).eq. 3) then
!
! Gaussian model:
!
                  hh=-(h*h)/(aa(is)*aa(is))
                  cova2 = cova2 +cc(is)*exp(hh)
            else
!
! Power model:
!
                  cov1  = PMX - cc(is)*(h**aa(is))
                  cova2 = cova2 + cov1
            endif
      end do
      return
      end





      subroutine ksol(nright,neq,nsb,a,r,s,ising)

!-----------------------------------------------------------------------
!
!                Solution of a System of Linear Equations
!                ****************************************
!
!
!
! INPUT VARIABLES:
!
!   nright,nsb       number of columns in right hand side matrix.
!                      for KB2D: nright=1, nsb=1
!   neq              number of equations
!   a()              upper triangular left hand side matrix (stored
!                      columnwise)
!   r()              right hand side matrix (stored columnwise)
!                      for kb2d, one column per variable
!
!
!
! OUTPUT VARIABLES:
!
!   s()              solution array, same dimension as  r  above.
!   ising            singularity indicator
!                      0,  no singularity problem
!                     -1,  neq .le. 1
!                      k,  a null pivot appeared at the kth iteration
!
!
!
! PROGRAM NOTES:
!
!   1. Requires the upper triangular left hand side matrix.
!   2. Pivots are on the diagonal.
!   3. Does not search for max. element for pivot.
!   4. Several right hand side matrices possible.
!   5. USE for ok and sk only, NOT for UK.
!
!
!-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
!
! If there is only one equation then set ising and return:
!
      if(neq.le.1) then
            ising = -1
            return
      endif
!
! Initialize:
!
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
!
! Start triangulation:
!
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
!
! Error checking - singular matrix:
!
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
!
! Finished triangulation, start solving back:
!
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
!
! Finished solving back, return:
!
      return
      end


!Notes:-
! standard deviations computed for gaussian case seem to be in error when a is large.


 
