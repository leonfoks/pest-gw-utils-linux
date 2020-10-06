program ppcov

! -- Program PPCOV prepares a covariance matrix file for pilot point parameters based on a
!    geostatistical structure file.

	use defn
	use inter
	implicit none

        integer, parameter     :: MAXZONE=200
        integer, parameter     :: MAXVAR=10

        logical                :: first
        integer                :: ifail,ierr,i,j,structunit,numstruct,numvario,oldzone,numzone,  &
                                  izone,ihigh,ilow,j1,outunit,iizone,istruct,nst,ipp,jpp,itcount
        integer                :: ippzone(MAXZONE),iwork(MAXZONE),ippstruct(MAXZONE),   &
                                  it(MAXVAR)
        real                   :: c0,pmx,x1,y1,x2,y2,rtemp,passmaxcov,cova2
        real                   :: cc(MAXVAR),ang(MAXVAR),aa(MAXVAR),anis(MAXVAR)
        real, allocatable      :: covar(:,:)
        real(kind=kind(1.d0))  :: minsep,minsep2,eastdiff,northdiff,distance,e0,n0
        character*10           :: azone,aprefix,atemp
        character*20           :: astruct
        character*200          :: structfile,outfile,aprompt

	type (modelgrid) gridspec
        type (geostructure), allocatable, dimension(:)    :: structure(:)
        type (variogram),    allocatable, dimension(:)    :: vario(:)



        write(amessage,5)
5       format(' Program PP2COV prepares a covariance matrix file for pilot point parameters based on ', &
        'a geostatistical structure file.')
        call write_message(leadspace='yes',endspace='yes')

! -- The pilot points file is read.

        pilot_points_file=' '
	call readfig(gridspec%specfile,pilotfile=pilot_points_file)
30      call read_pilot_points_file(ifail, &
        ' Enter name of pilot points file: ')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

! -- A check is made that no two pilot points are too close together.

        write(6,*)
61      write(6,63,advance='no')
63      format(' Enter minimum allowable separation for points in same zone: ')
        if(key_read(minsep).ne.0) go to 61
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 30
        end if
        imessage=0
        if(minsep.lt.0.0d0) go to 61
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
65                format(' The following same-zone points are separated by less than this ', &
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
        if(imessage.gt.0) go to 9900

! -- The geostatistical structure file is read.

        write(6,*)
70      continue
	aprompt=' Enter name of structure file: '
	call open_input_file(ifail,aprompt,structfile,structunit)
	if(ifail.ne.0) go to 9890
	if(escset.ne.0) then
          close(unit=structunit)
	  write(6,*)
	  escset=0
	  go to 61
	end if

! -- The geostatistical structure file is perused a first time to ascertain
!    array dimensions.

        call read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)
        if(ifail.ne.0) go to 9900
        if(numstruct.eq.0)then
          write(amessage,370) trim(structfile)
370       format(' No geostatistical structures found in file ',a,'.')
          go to 9890
        end if
        if(numvario.eq.0)then
          write(amessage,380) trim(structfile)
380       format(' No variograms found in file ',a,'.')
          go to 9890
        end if

! -- Memory is now allocated based on the contents of the structure file.

        allocate(structure(numstruct),vario(numvario),stat=ierr)
        if(ierr.ne.0)then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to continue execution.')
	  go to 9890
	end if

! -- The remainder of the structure file is now read.

        call read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
        structure,vario)
        if(ifail.ne.0) go to 9900

! -- Zones in the pilot points file are numbered.

        do i=1,num_pilot_points
          izone=pilot_point_zone(i)
          if(i.eq.1)then
            numzone=1
            ippzone(numzone)=izone
            oldzone=izone
            go to 53
          else
            if(izone.eq.oldzone) cycle
            oldzone=izone
            do j=1,i-1
              if(ippzone(j).eq.izone) go to 53
            end do
            numzone=numzone+1
            if(numzone.gt.MAXZONE) then
              write(amessage,52)
52            format(' Too many pilot point zones cited in pilot points file. Increase ', &
              'MAXZONE and re-compile program.')
              go to 9890
            end if
            ippzone(numzone)=izone
          end if
53        continue
        end do

! -- The zones are now sorted.

        iwork=ippzone
        ihigh=-99999999
        do i=1,numzone
          if(ippzone(i).gt.ihigh)ihigh=ippzone(i)
        end do
        ihigh=ihigh+1
        do i=1,numzone
          ilow=ihigh
          do j=1,numzone
            if(iwork(j).lt.ilow)then
              j1=j
              ilow=iwork(j)
            end if
          end do
          ippzone(i)=ilow
          iwork(j1)=ihigh
        end do

! -- A geostatistical structure is assigned to each zone.

100     continue
        write(6,*)
        i=0
108       continue
          i=i+1
          if(i.gt.numzone) go to 140
109       call num2char(ippzone(i),azone)
          write(6,110,advance='no') trim(azone)
110       format(' Enter structure to use for pilot point zone ',a,': ')
          read(5,'(a)') astruct
          call casetrans(astruct,'lo')
          astruct=adjustl(astruct)
          if(astruct(1:2).eq.'e ')then
            if(i.eq.1)then
              deallocate(structure,vario,stat=ierr)
              if(ierr.ne.0)then
                write(amessage,120)
120             format(' Cannot deallocate memory to backtrack in program.')
                go to 9890
              end if
              write(6,*)
              go to 70
            else
              i=i-1
              write(6,*)
              go to 109
            end if
          else
            do j=1,numstruct
              if(astruct.eq.structure(j)%structname)then
                ippstruct(i)=j
                go to 108
              end if
            end do
            write(6,130) trim(structfile)
130         format(/,' This structure is not found in file ',a,' - try again.',/)
            go to 109
          end if
        go to 108

140     continue

! -- The name of the output matrix file is now acquired.

        write(6,*)
143     continue
        aprompt = ' Enter name for output matrix file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1)then
          escset=0
          write(6,*)
          i=numzone
          go to 109
        end if
        write(6,145,advance='no')
145     format(' Enter pilot point prefix for parameter name (<Enter> if none): ')
        read(5,'(a)') aprefix
        aprefix=adjustl(aprefix)
        call casetrans(aprefix,'lo')
        if(aprefix(1:2).eq.'e ')then
          close(unit=outunit)
          write(6,*)
          go to 143
        end if

! -- The covariance matrix is now allocated space.

        allocate(covar(num_pilot_points,num_pilot_points),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,150)
150       format(' Too many pilot points - insufficient memory to continue execution.')
          go to 9890
        end if
        covar=0.0   ! An array

! -- Elements of the covariance matrix are now evaluated.

        write(6,160)
160     format(/,' Filling covariance matrix....')
        itcount=0
        do izone=1,numzone
          iizone=ippzone(izone)
          istruct=ippstruct(izone)
          if(structure(istruct)%transform.eq.1) itcount=itcount+1
          nst=structure(istruct)%numvariogram
          if(nst.gt.MAXVAR)then
            write(amessage,413) trim(structfile)
413         format(' At least one structure in file ',a,' uses too many nested ',  &
            'variograms. Increase MAXVAR and re-compile program.')
            go to 9890
          end if
          c0=structure(istruct)%nugget
          pmx=structure(istruct)%maxpowercov
          do i=1,nst
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
            it(i)=vario(j)%vartype
            cc(i)=structure(istruct)%variogram_contrib(i)
            ang(i)=vario(j)%angle
            aa(i)=vario(j)%a
            anis(i)=vario(j)%anis
          end do

          do ipp=1,num_pilot_points
            if(pilot_point_zone(ipp).ne.iizone)cycle
            e0=pilot_point_east(ipp)
            n0=pilot_point_north(ipp)
            exit
          end do
          first=.true.
          do ipp=1,num_pilot_points
            if(pilot_point_zone(ipp).ne.iizone)cycle
            x1=pilot_point_east(ipp)-e0
            y1=pilot_point_north(ipp)-n0
            do jpp=1,num_pilot_points
              if(pilot_point_zone(jpp).ne.iizone) cycle
              x2=pilot_point_east(jpp)-e0
              y2=pilot_point_north(jpp)-n0
              rtemp = cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                            ang,anis,first,passmaxcov)
              first=.false.
              covar(ipp,jpp)=rtemp
              if(ipp.ne.jpp)covar(jpp,ipp)=rtemp
            end do
          end do
        end do

! -- The covariance matrix is now written to a matrix file.

        write(outunit,421) num_pilot_points,num_pilot_points,1
421     format(3i6)
        do i=1,num_pilot_points
          write(outunit,430) (covar(i,j),j=1,num_pilot_points)
430       format(8(1x,1pg14.7))
        end do
        write(outunit,440)
440     format('* row and column names')
        do i=1,num_pilot_points
          call casetrans(pilot_point_id(i),'lo')
          if(aprefix.eq.' ')then
            write(outunit,450) trim(pilot_point_id(i))
450         format(1x,a)
          else
            write(outunit,450) trim(aprefix)//trim(pilot_point_id(i))
          end if
        end do

        close(unit=outunit)
        write(6,490) trim(outfile)
490     format(' - file ',a,' written ok.')

! -- If any parameters are log-transformed a warning is issued.

        if(itcount.gt.0)then
          write(amessage,520)
520       format(' Warning: in any future processing of this covariance matrix, sensitivities ', &
          'for parameters with a log-variogram must be taken with respect to the log ', &
          'of the parameters.')
          call write_message(leadspace='yes',endspace='yes')
        end if


        go to 9900


9000    write(amessage,9010)
9010    format(' Too many pilot point zones - increase MAXZONE and re-compile program.')
        go to 9890


9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files()
        call free_point_mem()
        deallocate(structure,vario,covar,stat=ierr)

end program ppcov



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



