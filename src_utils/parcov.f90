program parcov

! -- Program PARCOV writes a covariance matrix file for a set of parameters whose names and coordinates are provided.

	use defn
	use inter
	implicit none

        logical                :: first
        integer, parameter     :: MAXVAR=10
        integer, parameter     :: MAXPAR=3000
        integer                :: ifail,ierr,i,j,structunit,numstruct,numvario,  &
                                  outunit,istruct,nst,itcount
        integer                :: ipar,jpar,npar,nbb,iunit,iline
        integer                :: it(MAXVAR)
        real                   :: c0,pmx,x1,y1,x2,y2,rtemp,passmaxcov,cova2
        real                   :: cc(MAXVAR),ang(MAXVAR),aa(MAXVAR),anis(MAXVAR)
        real, allocatable      :: covar(:,:)
        real(kind=kind(1.d0))  :: e0,n0,ee,nn
        real(kind=kind(1.d0))  :: ecoord(MAXPAR),ncoord(MAXPAR)
        character*12           :: aapar,aline
        character*12           :: apar(MAXPAR)
        character*20           :: astruct,atemp
        character*25           :: aastruct
        character*200          :: structfile,outfile,aprompt,afile,infile

        type (geostructure), allocatable, dimension(:)    :: structure(:)
        type (variogram),    allocatable, dimension(:)    :: vario(:)


        write(amessage,1)
1       format(' Program PARCOV writes a covariance matrix file for parameters whose ',  &
        'names and coordinates are provided.')
        call write_message(leadspace='yes',endspace='yes')

! -- The coordinates file is read

        infile=' '
4       write(6,5,advance='no')
5       format(' Enter name for parameter coordinates file: ')
        read(5,'(a)') afile
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E ')) go to 9900
        nbb=len_trim(afile)
        call getfile(ifail,afile,infile,1,nbb)
        if(ifail.ne.0) go to 4
        iunit=nextunit()
        open(unit=iunit,file=infile,status='old',iostat=ierr)
        call addquote(infile,afile)
        if(ierr.ne.0)then
          write(6,6) trim(afile)
6         format(' Cannot open file ',a,' - try again.')
          go to 4
        end if
        npar=0
        iline=0
        do
          npar=npar+1
          if(npar.gt.MAXPAR)then
            write(amessage,10) trim(afile)
10          format(' Too many parameters are provided in file ',a,'. Increase ',  &
            'MAXPAR and re-compile program.')
            go to 9890
          end if
9         continue
          iline=iline+1
          read(10,'(a)',end=50) cline
          if(cline.eq.' ') go to 9
	  call linesplit(ifail,3)
          if(ifail.ne.0) then
            call num2char(iline,aline)
            write(amessage,12) trim(aline),trim(afile)
12          format(' Insufficient entries found on line ',a,' of file ',a,'.')
            go to 9890
          end if
          ecoord(npar)=char2double(ifail,2)
          if(ifail.ne.0) go to 9100
          ncoord(npar)=char2double(ifail,3)
          if(ifail.ne.0) go to 9100
          if(right_word(1)-left_word(1).gt.11)then
            call num2char(iline,aline)
            write(amessage,15) trim(aline),trim(afile)
15          format(' Parameter name greater than 12 characters at line ',a,' of file ',a'.')
            go to 9890
          end if
          apar(npar)=cline(left_word(1):right_word(1))
          call casetrans(apar(npar),'lo')
        end do
50      continue
        npar=npar-1
        if(npar.eq.0)then
          write(amessage,55) trim(afile)
55        format(' No parameters provided in file ',a,'.')
          go to 9890
        end if
        close(unit=iunit)

! -- A check is made that parameter names are not duplicated.

        if(npar.gt.1)then
          do ipar=2,npar
            aapar=apar(ipar)
            do jpar=1,ipar-1
              if(apar(jpar).eq.aapar)then
                write(amessage,60) trim(aapar),trim(afile)
60              format(' Parameter name "',a,'" duplicated in file ',a,'.')
                go to 9890
              end if
            end do
          end do
          do ipar=2,npar
            ee=ecoord(ipar)
            nn=ncoord(ipar)
            do jpar=1,ipar-1
              if(equals(ee,ecoord(jpar)))then
                if(equals(nn,ncoord(jpar)))then
                  write(amessage,62) trim(apar(ipar)),trim(apar(jpar)),trim(afile)
62                format(' Parameters "',a,'" and "',a,'" are provided with identical ',  &
                  'coordinates in file ',a,'.')
                  go to 9890
                end if
              end if
            end do
          end do
        end if
        write(6,56) trim(afile)
56      format(' - file ',a,' read ok.')

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
	  go to 4
	end if

! -- The geostatistical structure file is perused a first time to ascertain
!    array dimensions.

        call read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)
        if(ifail.ne.0) go to 9900
        if(numstruct.eq.0)then
          call addquote(structfile,afile)
          write(amessage,370) trim(afile)
370       format(' No geostatistical structures found in file ',a,'.')
          go to 9890
        end if
        if(numvario.eq.0)then
          call addquote(structfile,afile)
          write(amessage,380) trim(afile)
380       format(' No variograms found in file ',a,'.')
          go to 9890
        end if

! -- Memory is now allocated based on the contents of the structure file.

        allocate(structure(numstruct),vario(numvario),stat=ierr)
        if(ierr.ne.0)then
	  write(amessage,51)
51        format(' Cannot allocate sufficient memory to continue execution.')
	  go to 9890
	end if

! -- The remainder of the structure file is now read.

        call read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
        structure,vario)
        if(ifail.ne.0) go to 9900
        call addquote(structfile,afile)
        write(6,56) trim(afile)

! -- A geostatistical structure is assigned to the parameter set.

100     continue
        write(6,*)
105     write(6,110,advance='no')
110     format(' Enter structure to use for parameters: ')
        read(5,'(a)') aastruct
        nbb=len_trim(aastruct)
        call getfile(ifail,aastruct,astruct,1,nbb)
        if(ifail.ne.0) go to 105
        call casetrans(astruct,'lo')
        astruct=adjustl(astruct)
        if(astruct(1:2).eq.'e ')then
          deallocate(structure,vario,stat=ierr)
          if(ierr.ne.0)then
            write(amessage,120)
120         format(' Cannot deallocate memory to backtrack in program.')
            go to 9890
          end if
          write(6,*)
          go to 70
        end if
        do j=1,numstruct
          if(astruct.eq.structure(j)%structname)then
            istruct=j
            go to 140
          end if
        end do
        call addquote(structfile,afile)
        write(6,130) trim(afile)
130     format(/,' This structure is not found in file ',a,' - try again.',/)
        go to 105

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
          go to 105
        end if

! -- The covariance matrix is now allocated space.

        allocate(covar(npar,npar),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,150)
150       format(' Too many parameters - insufficient memory to continue execution.')
          go to 9890
        end if
        covar=0.0   ! An array

! -- Elements of the covariance matrix are now evaluated.

        itcount=0
        write(6,160)
160     format(/,' Filling covariance matrix....')
        nst=structure(istruct)%numvariogram
        if(structure(istruct)%transform.eq.1) itcount=itcount+1
        if(nst.gt.MAXVAR)then
          write(amessage,413) trim(structfile)
413       format(' At least one structure in file ',a,' uses too many nested ',  &
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
415       format(' Specifications for variogram "',a,'" cited in structure "',a,  &
            '" in file ',a,' not found in this file.')
          go to 9890
420       continue
          it(i)=vario(j)%vartype
          cc(i)=structure(istruct)%variogram_contrib(i)
          ang(i)=vario(j)%angle
          aa(i)=vario(j)%a
          anis(i)=vario(j)%anis
        end do

        e0=ecoord(1)
        n0=ncoord(1)
        first=.TRUE.
        do ipar=1,npar
          x1=ecoord(ipar)-e0
          y1=ncoord(ipar)-n0
          do jpar=1,ipar
            x2=ecoord(jpar)-e0
            y2=ncoord(jpar)-n0
            rtemp = cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                            ang,anis,first,passmaxcov)
            covar(jpar,ipar)=rtemp
            if(ipar.ne.jpar)covar(ipar,jpar)=rtemp
          end do
        end do

! -- The covariance matrix is now written to a matrix file.

        write(outunit,421) npar,npar,1
421     format(3i6)
        do i=1,npar
          write(outunit,430) (covar(i,j),j=1,npar)
430       format(8(1x,1pg14.7))
        end do
        write(outunit,440)
440     format('* row and column names')
        do ipar=1,npar
          write(outunit,450) trim(apar(ipar))
450       format(1x,a)
        end do

        close(unit=outunit)
        call addquote(outfile,afile)
        write(6,490) trim(afile)
490     format(' - file ',a,' written ok.')

! -- If any parameters are log-transformed a warning is issued.

        if(itcount.gt.0)then
          call addquote(infile,afile)
          write(amessage,520) trim(afile)
520       format(' Warning: according to information in the structure file, this ',  &
          'covariance matrix pertains to the logs of parameters cited in file ',a,   &
          '. Future processing of this matrix must take this into account.')
          call write_message(leadspace='yes',endspace='yes')
        end if


        go to 9900

9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Cannot read parameter coordinates from line ',a,' of file ',a,'.')
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files()
        deallocate(structure,vario,covar,stat=ierr)

end program parcov



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



