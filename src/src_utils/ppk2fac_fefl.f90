!     Last change:  JD    9 May 2003   11:40 am

program ppk2fac_fefl

! -- Program PPK2FAC_FEFL performs the first step in kriging from a set of
!    pilot points to a two-dimensional FEFLOW mesh. Use of these factors is
!    made by the FAC2FEFL utility which performs the actual spatial interpolation.

        use defn
        use inter

        implicit none

        integer     :: ifail,i,ierr,structunit,                                      &
                       numint,currint,newint,ihuge,itemp,numstruct,numvario,outunit, &
                       nbb,j,izone,ndat,ngrid,istruct,k_ktype,n_nst,regunit,itrans,  &
                       ispace,jj,femunit
        real        :: s_skmean,c_c0,pmx
        double precision :: minsep,minsep2,eastdiff,northdiff,distance
        double precision :: east1,north1
        integer     :: iloc(1)
        integer     :: i_it(MAX_STRUCT_VARIO)
        real        :: c_cc(MAX_STRUCT_VARIO),a_ang(MAX_STRUCT_VARIO),  &
                       a_aa(MAX_STRUCT_VARIO),a_anis(MAX_STRUCT_VARIO)
        integer, allocatable, dimension(:)             :: icellno
        integer, allocatable, dimension(:)             :: inumdat
        integer, allocatable, dimension(:)             :: minpt,maxpt
        real, allocatable, dimension(:)                :: radmax,radmin
        real, allocatable, dimension(:)                :: eastgrid,northgrid
        real,allocatable, dimension(:)                 :: eastdat,northdat,valdat
        real, allocatable, dimension(:)                :: totcov
        real, allocatable, dimension(:,:)              :: regarray
        character (len=1)                              :: aa
        character (len=1)                              :: alogtrans,aoutfile
        character (len=5)                              :: premspace
        character (len=1),  allocatable, dimension(:)  :: akrig
        character (len=10), allocatable, dimension(:)  :: astructure
        character (len=15)                             :: anum1,atemp
        character (len=20)                             :: bin_unform
        character (len=200)                            :: aprompt,intfile,structfile, &
                                                          atempf,outfile,regfile,femfile

        type (geostructure), allocatable, dimension(:) :: structure(:)
        type (variogram), allocatable, dimension(:)    :: vario(:)


        integer, parameter :: MAXINT=100
        integer            :: intval(MAXINT),iwork(MAXINT)

        include 'bin_unform.inc'

        write(amessage,5)
5       format(' Program PPK2FAC_FEFL calculates point-to-element factors by which kriging ',&
       'is undertaken from a set of pilot points to a FEFLOW mesh.')
        call write_message(leadspace='yes',endspace='yes')

! -- The FEFLOW element property file is read.

6       write(6,7,advance='no')
7       format(' Enter name of FEFLOW FEM file for current project: ')
        read(5,'(a)',err=6) atempf
        if(atempf.eq.' ') go to 6
        atempf=adjustl(atempf)
        if(index(eschar,atempf(1:2)).ne.0) go to 9900
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,femfile,1,nbb)
        if(ifail.ne.0) go to 6
        femunit=nextunit()
        open(unit=femunit,file=femfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(6,8) trim(femfile)
8         format(/,' Cannot open file ',a,' - try again.',/)
          go to 6
        end if
        do
          read(femunit,'(a)',iostat=ierr) cline
          if(ierr.ne.0)then
            write(amessage,9) trim(femfile)
9           format(' Cannot find the "DIMENS" keyword in file ',a,'.')
            go to 9890
          end if
          cline=adjustl(cline)
          if(cline(1:7).eq.'DIMENS ') exit
        end do
        read(femunit,*,iostat=ierr) numnode_f,numelem_f
        if(ierr.ne.0)then
          write(amessage,11) trim(femfile)
11        format(' Cannot read mesh specifications from file ',a,'.')
          go to 9890
        end if
        close(unit=femunit)

        write(6,*)
10      continue
        call read_feflow_elem_prop_file(ifail,' Enter name of FEFLOW element property file: ')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          close(unit=femunit)
          write(6,*)
          go to 6
        end if

! -- The pilot points file is read.

        write(6,*)
30      call read_pilot_points_file(ifail, &
        ' Enter name of pilot points file: ')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          deallocate(eastelem_f,northelem_f,zoneelem_f,stat=ierr)
          if(ierr.ne.0) then
            write(amessage,57)
57          format(' Memory management error: cannot continue execution.')
            go to 9890
          end if
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
              if((distance.le.minsep2).and. &
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
80      continue
        aprompt=' Enter name of structure file: '
        call open_input_file(ifail,aprompt,structfile,structunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          write(6,*)
          escset=0
          go to 61
        end if

! -- The different integers comprising the index data are identified.

        numint=1
        currint=zoneelem_f(1)
        intval(1)=currint
        idat_travel: do i=1,numelem_f
          newint=zoneelem_f(i)
          if(newint.eq.currint) cycle idat_travel
          currint=newint
          prev_int: do j=1,numint
              if(newint.eq.intval(j)) cycle idat_travel
          end do prev_int
          numint=numint+1
          if(numint.gt.MAXINT) then
            call num2char(MAXINT,anum1)
            write(amessage,150) trim(epfile_f),trim(anum1),trim(epfile_f)
150         format(' PPK2FAC_FEFL expects zones for pilot point assignment ', &
            'to be defined in the sixth column of file ',a,' as the element property ',  &
            'value. A maximum of ',a,' such zones can be defined. There are more zones ', &
            'than this in file ',a,'. Edit PPK2FAC_FEFL source code, increase parameter ', &
            'MAXINT and re-compile program.')
            go to 9890
          end if
          intval(numint)=newint
        end do idat_travel

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
        akrig(numint),minpt(numint),maxpt(numint),stat=ierr)  !de-allocate below for "e"
        if(ierr.ne.0) go to 9200

        radmax=0.0
        radmin=0.0
        akrig=' '
        astructure=' '

199     continue
        premspace=' '
        ispace=4
200     write(6,*)
205     write(6,210)
210     format(' The following zones have been detected in the element property ',&
        'file:-')
        i=0
255     i=i+1
        if(i.gt.numint) go to 350
        write(6,*)
274     continue
259     call num2char(intval(i),anum1)
260     write(6,270) trim(anum1)
270     format('    For zone characterised by integer value of ',a,':- ')
        atemp=' '
275     write(6,280,advance='no')
280     format('    Enter structure name (blank if no interpolation for this zone): ')
        read(5,'(a)') atemp
        if(index(eschar,atemp(1:2)).ne.0)then
          if(i.gt.2) then
            i=i-2
            go to 255
          else if(i.eq.2) then
            go to 199
          else
            deallocate(astructure,radmax,radmin,akrig,minpt,maxpt,stat=ierr)
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
          write(amessage,285) premspace(1:ispace)
285       format(a,'Structure name must be 10 characters or less - try again.')
          call write_message()
          go to 274
        end if
        astructure(i)=atemp
        if(astructure(i).eq.' ') then
          go to 255
        end if
        call casetrans(astructure(i),'lo')
289     write(6,290,advance='no') premspace(1:ispace)
290     format(a,'Perform simple or ordinary kriging [s/o]: ')
        read(5,'(a)') akrig(i)
        if(akrig(i).eq.' ') go to 289
        if(index(eschar,akrig(i)).ne.0) then
          write(6,*)
          go to 274
        end if
        call casetrans(akrig(i),'lo')
        if((akrig(i).ne.'s').and.(akrig(i).ne.'o')) go to 289
299     write(6,300,advance='no') premspace(1:ispace)
300     format(a,'Enter search radius: ')
        itemp=key_read(radmax(i))
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 289
        end if
        if(itemp.ne.0) then
          write(6,310) premspace(1:ispace)
310       format(a,'Data input error  - try again.')
          go to 299
        end if
        if(radmax(i).le.0.0d0)then
          write(6,311) premspace(1:ispace)
311       format(a,'Must be greater than zero - try again.')
          go to 299
        end if
        if(radmax(i).gt.1.0e15)radmax(i)=1.0e15

329     write(6,330,advance='no') premspace(1:ispace)
330     format(a,'Enter minimum number of pilot points to use for interpolation: ')
        itemp=key_read(minpt(i))
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 299
        end if
        if(itemp.ne.0) then
          write(6,310) premspace(1:ispace)
          go to 329
        end if
        if(minpt(i).lt.1)then
          write(6,332) premspace(1:ispace)
332       format(a,'Must be greater than 1 - try again.')
          go to 329
        end if
335     write(6,337,advance='no') premspace(1:ispace)
337     format(a,'Enter maximum number of pilot points to use for interpolation: ')
        itemp=key_read(maxpt(i))
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 329
        end if
        if(itemp.ne.0) then
          write(6,310) premspace(1:ispace)
          go to 335
        end if
        if(maxpt(i).lt.minpt(i))then
          write(6,340) premspace(1:ispace)
340       format(a,'Must be greater than min. no. of points - try again.')
          go to 335
        end if
        if(maxpt(i).gt.500)then
          write(6,341) premspace(1:ispace)
341       format(a,'Must not be greater than 500 - try again.')
          go to 335
        end if
        go to 255

350     continue

        do i=1,numint
          if(astructure(i).ne.' ')go to 351
        end do
        write(amessage,352)
352     format(' Interpolation must take place for at least one zone ', &
        '- try again.')
        call write_message(leadspace='yes',endspace='yes')
        go to 199
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
          open(unit=outunit,file=outfile,status='new',form=bin_unform,iostat=ierr)
          if(ierr.ne.0) then
            open(unit=outunit,file=outfile,status='replace',form=bin_unform,iostat=ierr)
          end if
        end if
        if(ierr.ne.0) then
          write(amessage,358) trim(outfile)
358       format(' Unable to open file ',a,' for output - try again.')
          go to 353
        end if

        write(6,*)
        aprompt = ' Enter name for regularisation information file: '
        call open_output_file(ifail,aprompt,regfile,regunit)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1)then
          escset=0
          close(unit=outunit)
          go to 353
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
        if(ierr.ne.0) go to 9200

! -- The remainder of the structure file is now read.

        call read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
        structure,vario)
        if(ifail.ne.0) go to 9900

! -- Next the east and north coordinate (relative to the first node),
!    is calculated for every centroid in the mesh.

        east1=eastelem_f(1)
        north1=northelem_f(1)
        do i=1,numelem_f
          eastelem_f(i)=eastelem_f(i)-east1
          northelem_f(i)=northelem_f(i)-north1
        end do

        allocate(eastdat(num_pilot_points),northdat(num_pilot_points),  &
        valdat(num_pilot_points),inumdat(num_pilot_points),  &
        regarray(num_pilot_points,num_pilot_points),totcov(num_pilot_points),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(eastgrid(numelem_f),northgrid(numelem_f),icellno(numelem_f),stat=ierr)
        if(ierr.ne.0) go to 9200
        regarray=-1.1e35           !regarray is an array
        totcov=-1.1e35             !totcov is an array

! -- The coordinates of the pilot points have the first node easting and
!    northing subtracted.

        do i=1,num_pilot_points
          pilot_point_east(i)=pilot_point_east(i)-east1
          pilot_point_north(i)=pilot_point_north(i)-north1
        end do

! -- Kriging weights are now evaluated zone by zone for cell centres in
!    the finite-difference grid using GSLIB subroutines.

        call addquote(pilot_points_file,atempf)
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        call addquote(epfile_f,atempf)
        if(aoutfile.eq.'f')then
          write(outunit,'(a)') trim(atempf)
        else
          write(outunit) atempf
        end if
        if(aoutfile.eq.'f')then
          write(outunit,*) numelem_f
        else
          write(outunit) numelem_f
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

        jj=0
        do izone=1,numint
          ngrid=0
          if(astructure(izone).eq.' ') go to 410
          itemp=intval(izone)
          do i=1,numelem_f
            if(zoneelem_f(i).eq.itemp)then
              ngrid=ngrid+1
              eastgrid(ngrid)=eastelem_f(i)
              northgrid(ngrid)=northelem_f(i)
              icellno(ngrid)=i
            end if
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
399         format(/,' Warning: no pilot points assigned to zone ',a,'.')
            go to 410
          end if
          call num2char(intval(izone),anum1)
          write(6,401) trim(anum1)
401       format(/,' Carrying out interpolation for zone ',a,'....')
          do i=1,numstruct
            if(astructure(izone).eq.structure(i)%structname)then
              istruct=i
              go to 405
            end if
          end do
          call num2char(intval(izone),anum1)
          write(amessage,402) trim(structfile),trim(astructure(izone)),trim(anum1)
402       format(' Structure file ',a,' does not include specifications for ', &
          'structure "',a,'" needed for interpolation to zone ',a,'.')
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
          call kb2d(ngrid,ndat,inumdat,icellno,eastgrid,northgrid, &
          outfile,aoutfile,outunit,pmx,itrans,num_pilot_points,regarray,totcov)

410     continue
        end do

        write(6,*)
        write(6,450) trim(outfile)
450     format('  - kriging factors written to file ',a)

        write(regunit,'(i5)') num_pilot_points
        do i=1,num_pilot_points
          write(regunit,451) trim(pilot_point_id(i)), &
          pilot_point_east(i)+east1, &
          pilot_point_north(i)+north1
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

9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890

9890    call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        call free_point_mem()
        deallocate(eastgrid,northgrid,icellno,stat=ierr)
        deallocate(eastdat,northdat,valdat,inumdat,stat=ierr)
        deallocate(astructure,radmax,radmin,akrig,minpt,maxpt,stat=ierr)
        deallocate(structure,vario,stat=ierr)
        deallocate(regarray,stat=ierr)
        deallocate(totcov,stat=ierr)
        deallocate(eastelem_f,northelem_f,zoneelem_f,stat=ierr)

end program ppk2fac_fefl




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



      subroutine kb2d(n_npts,n_ndat,inumdat,icellno,epoint,npoint,outfile, &
      aoutfile,outunit,pmx,itrans,npp,regarray,totcov)

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

     integer n_npts,n_ndat,outunit,itrans                           !jd
     integer icellno(n_npts),inumdat(n_ndat)                        !jd
     real pmx                                                       !jd
     real epoint(n_npts),npoint(n_npts)                             !jd
     real regarray(npp,npp)                                         !jd
     real totcov(npp)                                               !jd
     character*(*) outfile,aoutfile                                 !jd

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
            do 6 id=1,nd
                  dx = x(id) - xloc
                  dy = y(id) - yloc
                  h2 = dx*dx + dy*dy
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
105          format('   No. of elements for which factors were calculated        = ',i6)
             if(n_npts-nk.eq.0)then
             write(6,1051) n_npts-nk
1051         format('   No. of elements beyond search radius of any pilot point  = ',i6)
             else
             write(6,1052) n_npts-nk
1052         format('   No. of elements beyond search radius of any pilot point  = ',i6,' : WARNING')
             end if
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



subroutine whichone_i(ifail,npar,ipar,ival,jval)

! -- Subroutine whichone_i locates an integer value in an array of integers.

        integer npar,ipar,i
        integer ifail
        integer jval
        integer ival(npar)

        ifail=0
        if((ipar.lt.1).or.(ipar.gt.npar)) ipar=1
        if(jval.eq.ival(ipar)) return
        if(ipar.ne.npar)then
          do 20 i=ipar+1,npar
          if(jval.eq.ival(i))then
            ipar=i
            return
          end if
20        continue
        end if
        if(ipar.ne.1)then
          do 40 i=ipar-1,1,-1
          if(jval.eq.ival(i)) then
            ipar=i
            return
          end if
40        continue
        end if
        ifail=1
        return

 end subroutine whichone_i



!Notes:-
! standard deviations computed for gaussian case seem to be in error when a is large.

