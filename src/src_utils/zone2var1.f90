!     Last change:  JD   25 Jan 2007    7:18 am
program zone2var1

! -- Program ZONE2VAR1 computes a variogram of parameter values, and derivatives of this variogram,
!    based on ZONMDEF-produced parameterization.

	use defn
	use inter

	implicit none

        integer, parameter :: MAXLAY=80
        integer, parameter :: MAXSTRUCTZONE=10

        logical            :: first

        integer            :: n,i,j,m
        integer            :: jpar,izone,ipar,istructzone,oldzone,nstructzone, &
                              numstruct,numvario,itemp,istruct,maxsamprow,maxsampcol, &
                              is,iobs,nobs,npar
        integer            :: ifail,ierr,idate,iheader,iline
        integer            :: parnum
        integer            :: maxent,nument
        integer            :: ncol,nrow,mlay,ilay,icol,irow,nlay
        integer            :: d2punit,ircunit,structunit,varunit,insunit,derunit,pestunit, &
                              zonunit
        integer            :: nst,iv,ist,dflag,nb,zflag
        integer            :: nangle,iangle,maxsampang

        real               :: c0,PMX
        real               :: x1,x2,y1,y2
        real               :: cova2

        double precision   :: rtemp,contrib,offset_angle,tan_offset_angle,temp,rtemp2
        double precision   :: pi,one,ln10
        double precision   :: e0,n0,ediff,ndiff,covar_zero
        double precision   :: weight
        double precision   :: cosangle,sinangle,h,r,o,theta,phi,alpha

        character*1        :: at,ayn
        character*10       :: aline,alay,anum
        character*12       :: aapar,prefix,oldprefix
        character*20       :: atemp,azone,adist,agroup,aangle

        character*200      :: aprompt,afile,bfile,cfile
        character*200      :: d2pfile,ircfile,sintfile,structfile,varfile,insfile,derfile, &
                              pestfile,zonfile

        integer                       :: iszone(MAXSTRUCTZONE)
        integer                       :: it(4)

        real                          :: cc(4),aa(4),anis(4),ang(4)

        integer, allocatable          :: ival(:)
        integer, allocatable          :: layer(:)
        integer, allocatable          :: structnum(:)
        integer, allocatable          :: iarray(:,:),isarray(:,:)
        integer, allocatable          :: numsampcol(:),numsamprow(:)
        integer, allocatable          :: numcellpar(:)
        integer, allocatable          :: structzonepar(:)
        integer, allocatable          :: num_col_e(:,:),num_row_e(:,:)
        integer, allocatable          :: num_par_row(:,:),num_par_col(:,:)
        integer, allocatable          :: numsampang(:,:),num_ang_e(:,:,:),num_par_ang(:,:,:)

        double precision, allocatable :: pval(:),pvalkeep(:)
        double precision, allocatable :: sampintcol(:),sampintrow(:)
        double precision, allocatable :: east(:),north(:)
        double precision, allocatable :: eastpar(:),northpar(:)
        double precision, allocatable :: h_col(:,:),h_row(:,:)
        double precision, allocatable :: h_row_limit(:),h_col_limit(:)
        double precision, allocatable :: gamma_row_t(:,:),gamma_row_e(:,:),  &
                                         gamma_col_t(:,:),gamma_col_e(:,:)
        double precision, allocatable :: deriv_row(:,:),deriv_col(:,:)
        double precision, allocatable :: x_row(:)
        double precision, allocatable :: angle(:)
        double precision, allocatable :: sampintang(:,:),h_ang_limit(:,:),h_ang(:,:,:)
        double precision, allocatable :: gamma_ang_t(:,:,:),gamma_ang_e(:,:,:)
        double precision, allocatable :: deriv_ang(:,:,:)

        character*12, allocatable :: apar(:)
        character*20, allocatable :: bobs(:)

        type (modelgrid) gridspec
        type (geostructure), allocatable   :: structure(:)
        type (variogram), allocatable      :: vario(:)


        write(amessage,5)
5       format(' Program ZONE2VAR1 computes a parameter variogram based on ZONMDEF-generated ',  &
        'parameterization.')
        call write_message(leadspace='yes',endspace='yes')

! -- Initialization

        pi = 3.141592654
        one=1.0d0
        ln10=log(10.0d0)

! -- The settings file is read.

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          go to 9890
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          go to 9890
        endif
        if((iheader.ne.0).or.(headerspec.eq.' ')) then
          write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
          go to 9890
        end if

! -- Grid specifications are obtained.

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

        allocate(iarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200
        iarray=0   ! an array

! -- The distributed-to-pest parameter value file is opened.

50      write(6,*)
        call open_input_file(ifail, &
        ' Enter name of ZONMDEF distributed-to-PEST parameter file: ',d2pfile,d2punit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          deallocate(iarray)
          call free_grid_mem(gridspec)
          go to 10
        end if

! -- The distributed-to-PEST parameter file is read and a zonation array is built. The number
!    of model layers featured in it is also ascertained.

        write(6,60) trim(d2pfile)
60      format('  - reading file ',a,'...')
        iline=1
        read(d2punit,*,err=9000,end=9050)
        iline=2
        read(d2punit,*,err=9000,end=9050) mlay
        if(mlay.eq.0)then
          call num2char(iline,aline)
          write(amessage,70) trim(aline),trim(d2pfile)
70        format(' MLAY should not be zero at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(mlay.gt.0)then
          allocate(layer(mlay),stat=ierr)
          if(ierr.ne.0) go to 9200
          do ilay=1,mlay
            iline=iline+1
            read(d2punit,*,err=9000,end=9050) layer(ilay)
          end do
        else
          allocate(layer(MAXLAY),stat=ierr)
          if(ierr.ne.0) go to 9200
          layer=0      ! an array
          nlay=0
        end if
        iline=iline+1
        read(d2punit,*,err=9000,end=9050) npar
        allocate(apar(npar),ival(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        do ipar=1,npar
          iline=iline+1
          read(d2punit,*,err=9000,end=9050) apar(ipar)
          call casetrans(apar(ipar),'lo')
          n=index(apar(ipar),'_')
          if(n.le.1)then
            write(amessage,80) trim(d2pfile)
80          format(' Incorrect parameter naming convention in file ',a,'. Parameter names must be ',  &
            'a common prefix followed by a number.')
            go to 9890
          end if
          aapar=apar(ipar)
          go to 91
90        continue
          n=index(aapar,'_')
          if(n.eq.0) go to 100
91        aapar(n:n)=char(7)
          m=n
          go to 90
100       prefix=aapar(1:m-1)
          call replace_character(prefix,char(7),'_')
          if(ipar.gt.1)then
            if(prefix.ne.oldprefix) then
              write(amessage,80) trim(d2pfile)
              go to 9890
            end if
          end if
          oldprefix=prefix
          atemp=aapar(m+1:)
          call char2num(ifail,atemp,ival(ipar))
          if(ifail.ne.0) then
            write(amessage,80) trim(d2pfile)
            go to 9890
          end if
        end do

        iline=iline+1
        read(10,*,err=9000,end=9050) maxent
        if(maxent.ne.1)then
          call num2char(iline,aline)
          write(amessage,120) trim(aline),trim(d2pfile)
120       format(' MAXENT is expected to be 1 at line ',a,' of file ',a,   &
          ' if this file is to be used with this program.')
          go to 9890
        end if
        do
          iline=iline+1
          read(d2punit,*,err=9000,end=200) icol,irow,ilay,nument,parnum,contrib
          if(nument.ne.1)then
            call num2char(iline,aline)
            write(amessage,130) trim(aline),trim(d2pfile)
130         format(' NUMENT expected to be 1 at line ',a,' of file ',a,'.')
            go to 9890
          end if
          if(.not.equals(contrib,one)) then
            call num2char(iline,aline)
            write(amessage,140) trim(aline),trim(d2pfile)
140         format(' CONTRIB must be 1.0 at line ',a,' of file ',a,  &
            ' for grid parameterization to be suitable for processing by this program.')
            go to 9890
          end if
          if(mlay.gt.0)then
            do i=1,mlay
              if(ilay.eq.layer(i)) go to 150
            end do
            call num2char(iline,aline)
            write(amessage,160) trim(aline),trim(d2pfile)
160         format(' Illegal layer number at line ',a,' of file ',a,'.')
            go to 9890
150         continue
          else
            if(nlay.eq.0)then
              nlay=1
              layer(nlay)=ilay
            else
              do i=1,nlay
                if(ilay.eq.layer(i)) go to 180
              end do
              nlay=nlay+1
              if(nlay.gt.MAXLAY)then
                write(amessage,170)
170             format(' Increase MAXLAY and re-compile program.')
                go to 9890
              end if
              layer(nlay)=ilay
180           continue
            end if
          end if
          iarray(icol,irow)=parnum
        end do
200     continue
        if(mlay.lt.0)mlay=nlay
        close(unit=d2punit)
        write(6,210) trim(d2pfile)
210     format('  - file ',a,' read ok.')

!        open(unit=99,file='debug.dat')            !debug
!        write(99,*) 'mlay = ',mlay                !debug
!        do ilay=1,mlay                            !debug
!          write(99,*) layer(ilay)                 !debug
!        end do                                    !debug
!        call flush(99)                            !debug
!        open(unit=98,file='debug.inf')                      !debug
!        write(98,*) ncol,nrow                               !debug
!        do irow=1,nrow                                      !debug
!          write(98,211) (iarray(icol,irow),icol=1,ncol)     !debug
!211       format(20i5)                                      !debug
!        end do                                              !debug
!        close(unit=98)                                      !debug


! -- Now we need current parameter values so, to make sure that we have got them all, we read
!    integer-real-correspondence files for all layers.

        allocate(pval(npar),pvalkeep(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        pval=-1.0d300   ! an array

        write(6,*)
        ilay=0
204     continue
        ilay=ilay+1
        if(ilay.gt.mlay) go to 300
205       continue
          call num2char(layer(ilay),alay)
          aprompt=' Enter integer-real correspondence file for layer '//trim(alay)//': '
          call open_input_file(ifail,aprompt,ircfile,ircunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0)then
            escset=0
            if(ilay.eq.1)then
              deallocate(layer,apar,ival,pval,pvalkeep)
              go to 50
            else
              write(6,*)
              ilay=ilay-1
              go to 205
            end if
          end if
          iline=0
          jpar=1
          do
            iline=iline+1
            read(ircunit,'(a)',err=9100,end=250) cline
            if(cline.eq.' ') cycle
            call linesplit(ifail,2)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,220) trim(aline),trim(ircfile)
220           format(' Insufficient entries on line ',a,' of file ',a,'.')
              go to 9890
            end if
            izone=char2int(ifail,1)
            if(ifail.ne.0) go to 9100
            rtemp=char2double(ifail,2)
            if(ifail.ne.0) go to 9100
            call whichone_i(ifail,npar,jpar,ival,izone)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,230) trim(aline),trim(ircfile),trim(d2pfile)
230           format(' Integer cited in column 1 at line ',a,' of file ',a,    &
              ' is not linked to any parameter cited in file ',a,'.')
              go to 9890
            end if
            if(pval(jpar).gt.-1.0d299)then
              if(.not.equals(rtemp,pval(jpar)))then
                write(amessage,240) trim(apar(jpar)),trim(ircfile)
240             format(' Value provided to zone corresponding to parameter "',a,'" in file ',a,' differs ', &
                'from value provided in a previously-read integer-real correspondence file.')
                go to 9890
              end if
            else
              pval(jpar)=rtemp
            end if
          end do
250       continue
          close(unit=ircunit)
          write(6,270) trim(ircfile)
270       format('  - file ',a,' read ok.')
          go to 204
300     continue

        do ipar=1,npar
          if(pval(ipar).lt.-1.0d299)then
            write(amessage,310) trim(apar(ipar))
310         format(' No value has been assigned to a zone corresponding to parameter "',a,  &
            '" in any integer-real correspondence file.')
            go to 9890
          end if
        end do
        pvalkeep=pval    ! arrays

! -- The transformation status of parameters is now acquired.

        write(6,*)
320     write(6,330,advance='no')
330     format(' Is this parameter type log-transformed or untransformed? [l/u]: ')
        read(5,*) at
        call casetrans(at,'lo')
        if(at.eq.' ') go to 320
        if(at.eq.'e')then
          ilay=mlay
          write(6,*)
          go to 205
        end if
        if((at.ne.'l').and.(at.ne.'u')) go to 320

! -- Some aspects of the experimental variogram search are obtained.

        write(6,*)
1101    write(6,1100)
1100    format(' Search directions for experimental variogram construction will take place ')
        write(6,1120)
1120    format(' in the grid row and column directions.')
1129    write(6,1130,advance='no')
1130    format('   How many other directions for parameter search? ')
        if(key_read(nangle).ne.0) go to 1129
        if(escset.ne.0) then
          write(6,*)
          escset=0
          go to 320
        end if
        if(nangle.ne.0)then
          allocate(angle(nangle),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if
        iangle=0
1146    continue
        if(nangle.ne.0)then
1147      iangle=iangle+1
          if(iangle.gt.nangle) go to 1170
1148      call num2char(iangle,aangle)
1149      write(6,1150,advance='no') trim(aangle)
1150      format('   Enter angle number ',a,' from row direction (0 to 180 degrees): ')
          if(key_read(angle(iangle)).ne.0) go to 1149
          if(escset.ne.0) then
            write(6,*)
            escset=0
            if(iangle.eq.1)then
              deallocate(angle)
              go to 1129
            else
              iangle=iangle-1
              go to 1148
            end if
          end if
          if((angle(iangle).lt.0.0).or.(angle(iangle).gt.180.0)) go to 1148
!          if(iangle.gt.1)then
!            do i=1,iangle-1
!              if(equals(angle(iangle),angle(i)))then
!                write(6,1160)
!1160            format(/,'   *** Cannot equal a previous angle - try again.',/)
!                go to 1148
!              end if
!            end do
!          end if
          go to 1147
1170      continue
        end if

! -- The structural zonation integer array is read.

        allocate(isarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200

        write(6,*)
335     continue
        aprompt=' Enter name of structural zonation integer array file: '
        call read_integer_array(ifail,aprompt,isarray,  &
        pm_header=headerspec,rows=nrow,columns=ncol)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          deallocate(isarray)
          if(nangle.eq.0)then
            go to 1101
          else
            iangle=nangle-1
            go to 1146
          end if
        end if
        sintfile=aprompt

! -- The structural integer array defines the area over which experimental variograms
!    are to be computed. Using this device, some parameters can actually be excluded from
!    taking part in these variogram computations.
!    It also defines structural subdivision of this area.

        do irow=1,nrow
          do icol=1,ncol
            if((iarray(icol,irow).eq.0).or.(isarray(icol,irow).eq.0))then
              iarray(icol,irow)=0
              isarray(icol,irow)=0
            end if
          end do
        end do

 ! -- The number of zones in the structural integer array is defined.

        istructzone=0
        oldzone=0
        do irow=1,nrow
          do icol=1,ncol
            if(isarray(icol,irow).ne.0)then
              izone=isarray(icol,irow)
              if(izone.eq.oldzone) cycle
              if(istructzone.eq.0)then
                istructzone=1
                iszone(istructzone)=izone
                oldzone=izone
              else
                oldzone=izone
                do i=1,istructzone
                  if(iszone(i).eq.izone) go to 350
                end do
                istructzone=istructzone+1
                if(istructzone.gt.MAXSTRUCTZONE)then
                  write(amessage,345)
345               format(' Increase MAXTRUCTZONE and re-compile program.')
                  go to 9890
                end if
                iszone(istructzone)=izone
              end if
            end if
350         continue
          end do
        end do
        nstructzone=istructzone

! -- The zone numbers are now ordered.

        if(nstructzone.gt.1)then
360       continue
          do istruct=2,nstructzone
            if(iszone(istruct).lt.iszone(istruct-1))then
              itemp=iszone(istruct)
              iszone(istruct)=iszone(istruct-1)
              iszone(istruct-1)=itemp
              go to 360
            end if
          end do
        end if

!        write(99,*)                                           !debug
!        write(99,*) ' Ordered structural zone numbers'        !debug
!        write(99,*) ' numstructzone = ',nstructzone           !debug
!        do i=1,nstructzone                                    !debug
!          write(99,*) iszone(i)                               !debug
!        end do                                                !debug
!        call flush(99)                                        !debug

        allocate(structnum(nstructzone),sampintcol(nstructzone),numsampcol(nstructzone),   &
        sampintrow(nstructzone),numsamprow(nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(sampintang(nangle,nstructzone),numsampang(nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if

! -- The structure file is read.

        write(6,*)
370     continue
        aprompt=' Enter name of structure file: '
        call open_input_file(ifail,aprompt,structfile,structunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          write(amessage,375)
375       format(' Sorry, but you cannot backtrack from here as some data has been internally ',  &
          'altered and its original values cannot be retrieved.')
          go to 9890
        end if

! -- The geostatistical structure file is perused a first time to ascertain
!    array dimensions.

        call read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)
        if(ifail.ne.0) go to 9900
        if(numstruct.eq.0)then
          write(amessage,390) trim(structfile)
390       format(' No geostatistical structures found in file ',a,'.')
          go to 9890
        end if
        if(numvario.eq.0)then
          write(amessage,380) trim(structfile)
380       format(' No variograms found in file ',a,'.')
          go to 9890
        end if

! -- Memory is now allocated based on the contents of the structure file.

        allocate(structure(numstruct),vario(numvario),stat=ierr)
        if(ierr.ne.0) go to 9200

! -- The remainder of the structure file is now read.

        call read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
        structure,vario)
        if(ifail.ne.0) go to 9900
        write(6,282) trim(structfile)
282     format('  - structure file ',a,' read ok.')

! -- The names of the structures pertaining to different structure zones are acquired. Information
!    on experimental variogram construction is also generated.

400     continue
        i=0
455     i=i+1
          if(i.gt.nstructzone) go to 450
          write(6,*)
459       call num2char(iszone(i),anum)
460       write(6,470,advance='no') trim(anum)
470       format(' Enter structure name pertaining to structure zone ',a,': ')
          read(5,*) atemp
          if(atemp.eq.' ') go to 460
          atemp=adjustl(atemp)
          if(index(eschar,atemp(1:2)).ne.0)then
	    if(i.gt.1) then
              write(6,*)
	      i=i-1
	      go to 534
	    else if(i.eq.1) then
              close(unit=structunit)
              deallocate(structure,vario)
              write(6,*)
	      go to 370
	    end if
          end if
          if(len_trim(atemp).gt.10)then
            write(6,485)
485         format(/,' Structure name must be 10 characters or less - try again.',/)
            go to 460
          end if
          call casetrans(atemp,'lo')
          do j=1,numstruct
            if(atemp.eq.structure(j)%structname)then
              structnum(i)=j
              go to 510
            end if
          end do
          write(amessage,500) trim(structfile)
500       format(' There is no structure of this name in file ',a,'.')
          call write_message(leadspace='yes',endspace='yes')
          go to 460
510       continue
          if(at.eq.'l')then
            if(structure(j)%transform.ne.1)then
              write(amessage,502) trim(structfile),trim(atemp)
502           format(' According to file ',a,', structure "',a,'" does not have a ',  &
              'transformation status of "log".')
              call write_message(leadspace='yes',endspace='yes')
              go to 460
            end if
          else
            if(structure(j)%transform.ne.0)then
              write(amessage,503) trim(structfile),trim(atemp)
503           format(' According to file ',a,', structure "',a,'" does not have a ',  &
              'transformation status of "none".')
              call write_message(leadspace='yes',endspace='yes')
              go to 460
            end if
          end if
519       write(6,520,advance='no')
520       format(' Enter sampling length interval in row direction: ')
          if(key_read(sampintrow(i)).ne.0) go to 519
          if(escset.ne.0) then
            write(6,*)
            escset=0
            go to 459
          end if
          if(sampintrow(i).le.0.0d0) go to 519
524       write(6,525,advance='no')
525       format(' Enter number of sample intervals in row direction: ')
          if(key_read(numsamprow(i)).ne.0) go to 524
          if(escset.ne.0)then
            write(6,*)
            escset=0
            go to 519
          end if
          if(numsamprow(i).le.0) go to 524
530       write(6,531,advance='no')
531       format(' Enter sampling length interval in column direction: ')
          if(key_read(sampintcol(i)).ne.0) go to 530
          if(escset.ne.0) then
            write(6,*)
            escset=0
            go to 524
          end if
          if(sampintcol(i).le.0.0d0) go to 530
534       write(6,535,advance='no')
535       format(' Enter number of sample intervals in column direction: ')
          if(key_read(numsampcol(i)).ne.0) go to 534
          if(escset.ne.0)then
            write(6,*)
            escset=0
            go to 530
          end if
          if(numsampcol(i).le.0) go to 534

          iangle=0
1219      continue
          if(nangle.eq.0) go to 1280
            iangle=iangle+1
            if(iangle.gt.nangle) go to 1280
            write(aangle,'(f10.1)') angle(iangle)
            aangle=adjustl(aangle)
1220        write(6,1230,advance='no') trim(aangle)
1230        format(' Enter sampling length interval in ',a,' degrees direction: ')
            if(key_read(sampintang(iangle,i)).ne.0) go to 1220
            if(escset.ne.0) then
              write(6,*)
              escset=0
              if(iangle.eq.1)then
                go to 534
              else
                iangle=iangle-1
                write(aangle,'(f10.1)') angle(iangle)
                aangle=adjustl(aangle)
                go to 1240
              end if
            end if
            if(sampintang(iangle,i).le.0.0d0) go to 1220
1240        write(6,1250,advance='no') trim(aangle)
1250        format(' Enter number of samples in ',a,' degrees direction: ')
            if(key_read(numsampang(iangle,i)).ne.0) go to 1240
            if(escset.ne.0)then
              write(6,*)
              escset=0
              go to 1220
            end if
            if(numsampang(iangle,i).le.0) go to 1240
            go to 1219
1280      continue
          go to 455
450     continue

! -- The width of the sampling arc is obtained.

       write(6,*)
452    write(6,451,advance='no')
451    format(' Enter angular width of search arc in degrees: ')
       if(key_read(offset_angle).ne.0) go to 452
        if(escset.ne.0)then
          escset=0
          i=nstructzone
          write(6,*)
          if(nangle.eq.0)then
            go to 534
          else
            iangle=nangle
            write(aangle,'(f10.1)') angle(iangle)
            aangle=adjustl(aangle)
            go to 1240
          end if
        end if
        if(offset_angle.le.0.0d0) go to 452
        offset_angle=abs(offset_angle*0.5d0)
        tan_offset_angle=tan(offset_angle*pi/180.0d0)

!        write(99,*)             !debug
!        write(99,*) ' tan_offset_angle = ',tan_offset_angle   !debug
!        call flush(99)              !debug

! -- The names of model output files are now obtained.

        write(6,*)
550     aprompt = ' Enter name for experimental variogram output file: '
	call open_output_file(ifail,aprompt,varfile,varunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 452
        end if

        write(6,*)
554     write(6,555,advance='no')
555     format(' Write an instruction file to read this file?  [y/n]: ')
        read(5,*) ayn
        if(ayn.eq.' ') go to 554
        call casetrans(ayn,'lo')
556     continue
        if(ayn.eq.'e')then
          write(6,*)
          close(unit=varunit)
          go to 550
        else if(ayn.eq.'y')then
560       aprompt = ' Enter name for instruction file: '
	  call open_output_file(ifail,aprompt,insfile,insunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 554
          end if
        else if(ayn.eq.'n')then
          insfile=' '
        else
          go to 554
        end if

        write(6,*)
565     write(6,566,advance='no')
566     format(' Write a parameter derivatives file?  [y/n]: ')
        read(5,*) ayn
        if(ayn.eq.' ') go to 565
        call casetrans(ayn,'lo')
567     continue
        if(ayn.eq.'e')then
          write(6,*)
          if(insfile.eq.' ')then
            go to 554
          else
            ayn='y'
            close(unit=insunit)
            go to 556
          end if
        else if(ayn.eq.'y')then
570       aprompt = ' Enter name for parameter derivatives file: '
          call open_output_file(ifail,aprompt,derfile,derunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 565
          end if
        else if (ayn.eq.'n')then
          derfile=' '
        else
          go to 565
        end if
        dflag=0
        if(derfile.ne.' ')dflag=1

        write(6,*)
575     write(6,576,advance='no')
576     format(' Write a PEST building-block file?  [y/n]: ')
        read(5,*) ayn
        if(ayn.eq.' ') go to 575
        call casetrans(ayn,'lo')
577     continue
        if(ayn.eq.'e')then
          write(6,*)
          if(derfile.eq.' ')then
            go to 565
          else
            ayn='y'
            close(unit=derunit)
            go to 567
          end if
        else if(ayn.eq.'y')then
580       aprompt = ' Enter name for PEST building block file: '
	  call open_output_file(ifail,aprompt,pestfile,pestunit)
          if(ifail.ne.0) go to 9900
          if(escset.eq.1)then
            escset=0
            write(6,*)
            go to 575
          end if
        else if(ayn.eq.'n')then
          pestfile=' '
        else
          go to 575
        end if

        write(6,*)
585     write(6,586,advance='no')
586     format(' Write a ZONE2VAR2 input file?  [y/n]: ')
        read(5,*) ayn
        if(ayn.eq.' ') go to 585
        call casetrans(ayn,'lo')
        if(ayn.eq.'e')then
          write(6,*)
          if(pestfile.eq.' ')then
            go to 575
          else
            ayn='y'
            close(unit=pestunit)
            go to 577
          end if
        else if(ayn.eq.'y')then
          if(derfile.eq.' ')then
            write(6,587)
587         format(/,' This is not allowed unless a derivatives file is requested - try again.',/)
            go to 585
          end if
590       aprompt = ' Enter name for ZONE2VAR2 input file: '
	  call open_output_file(ifail,aprompt,zonfile,zonunit,file_format='unformatted')
          if(ifail.ne.0) go to 9900
          if(escset.eq.1)then
            escset=0
            write(6,*)
            go to 585
          end if
        else if(ayn.eq.'n')then
          zonfile=' '
        else
          go to 585
        end if
        zflag=0
        if(zonfile.ne.' ') zflag=1

! -- All input data has been acquired.

! -- Next the east and north coordinate (relative to top left corner of the grid),
!    is calculated for every cell in the finite-difference grid.

        write(6,601)
601     format(/,'  - computing theoretical semivariograms...')

        allocate(east(ncol),north(nrow),stat=ierr)
        if(ierr.ne.0) go to 9200

        temp=gridspec%delr(1)/2.0
        east(1)=temp
        do i=2,ncol
          temp=temp+(gridspec%delr(i)+gridspec%delr(i-1))/2.0
          east(i)=temp
        end do
        temp=-gridspec%delc(1)/2.0
        north(1)=temp
        do i=2,nrow
          temp=temp-(gridspec%delc(i)+gridspec%delc(i-1))/2.0
          north(i)=temp
        end do

!        write(99,*)                                 !debug
!        write(99,*) ' Cell centre easts --->'       !debug
!        do icol=1,ncol                              !debug
!          write(99,*) east(icol)                    !debug
!        end do                                      !debug
!        write(99,*)                                 !debug
!        write(99,*) ' Cell centre norths --->'      !debug
!        do irow=1,nrow                              !debug
!          write(99,*) north(irow)                   !debug
!        end do                                      !debug

! -- We now obtain the easting and northing of every zone that is involved in computation
!    of the variogram. This is done by averaging the centres of coordinates of cells
!    comprising each zone.

        allocate(eastpar(npar),northpar(npar),numcellpar(npar),structzonepar(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        eastpar=0.0d0         ! an array
        northpar=0.0d0        ! an array
        numcellpar=0          ! an array

        do irow=1,nrow
          do icol=1,ncol
            if(isarray(icol,irow).ne.0)then
              izone=iarray(icol,irow)
              call whichone_i(ifail,npar,jpar,ival,izone)
              if(ifail.ne.0)then
                write(amessage,602)
602             format(' Programming error type 2 - contact programmer.')
                go to 9890
              end if
              eastpar(jpar)=eastpar(jpar)+east(icol)
              northpar(jpar)=northpar(jpar)+north(irow)
              numcellpar(jpar)=numcellpar(jpar)+1
            end if
          end do
        end do
        do ipar=1,npar
          if(numcellpar(ipar).eq.0)then
            eastpar(ipar)=-1.0d300
            northpar(ipar)=-1.0d300
          else
            eastpar(ipar)=eastpar(ipar)/numcellpar(ipar)
            northpar(ipar)=northpar(ipar)/numcellpar(ipar)
          end if
        end do

!        write(99,*)                                         !debug
!        write(99,*) ' Parameter coordinates --->'           !debug
!        do ipar=1,npar                                      !debug
!          write(99,603) eastpar(ipar),northpar(ipar),numcellpar(ipar)   !debug
!603       format(1x,f12.3,2x,f12.3,2x,i4)                               !debug
!        end do                                                          !debug

! -- We now work out the structure zone in which each parameter lies.

        structzonepar=-99999999   ! an array
        jpar=1
        do irow=1,nrow
          do icol=1,ncol
            izone=iarray(icol,irow)
            if(izone.eq.0) cycle
            call whichone_i(ifail,npar,jpar,ival,izone)
            if(ifail.ne.0)then
              write(amessage,607)
607           format(' Programming error type 3 - contact programmer.')
              go to 9890
            end if
            itemp=isarray(icol,irow)
            if(structzonepar(jpar).eq.-99999999)then
              structzonepar(jpar)=itemp
            else if(structzonepar(jpar).ne.itemp)then
              write(amessage,605) trim(apar(jpar)),trim(d2pfile),trim(sintfile)
605           format(' The group of cells defining parameter "',a,  &
              '" as obtained from distributed-to-PEST parameter file ',a,     &
              ' spans more than one structure zone as defined ', &
              'in structure integer array file ',a,'.')
              go to 9890
            end if
          end do
        end do

!        write(99,*)                                         !debug
!        write(99,*) ' Parameter structure zones --->'      !debug
!        do ipar=1,npar                                      !debug
!          write(99,6061) apar(ipar),structzonepar(ipar)      !debug
!6061       format(1x,a,t15,i10)                              !debug
!        end do                                              !debug

! -- Each variogram angle cited in the structure file is corrected for the angle of rotation of the grid
!    (because grid row and column directions will form our coordinate system for this work).

        do i=1,numvario
          vario(i)%angle=vario(i)%angle+gridspec%rotation
        end do

! -- The sampling bin limits are now worked out.

        maxsamprow=0
        maxsampcol=0
        do istruct=1,nstructzone
          if(numsampcol(istruct).gt.maxsampcol) maxsampcol=numsampcol(istruct)
          if(numsamprow(istruct).gt.maxsamprow) maxsamprow=numsamprow(istruct)
        end do

        maxsampang=0
        do istruct=1,nstructzone
          do iangle=1,nangle
            if(numsampang(iangle,istruct).gt.maxsampang) maxsampang=numsampang(iangle,istruct)
          end do
        end do

        allocate(h_col(0:maxsampcol,nstructzone),h_row(0:maxsamprow,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(h_row_limit(nstructzone),h_col_limit(nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(h_ang(0:maxsampang,nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if
        allocate(h_ang_limit(nangle,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200

        do istruct=1,nstructzone
          h_col(0,istruct)=sampintcol(istruct)*0.5
          do i=1,numsampcol(istruct)
            h_col(i,istruct)=h_col(0,istruct)+sampintcol(istruct)*i
          end do
          h_col_limit(istruct)=h_col(numsampcol(istruct),istruct)
          h_row(0,istruct)=sampintrow(istruct)*0.5
          do i=1,numsamprow(istruct)
            h_row(i,istruct)=h_row(0,istruct)+sampintrow(istruct)*i
          end do
          h_row_limit(istruct)=h_row(numsamprow(istruct),istruct)
          if(nangle.ne.0)then
            do iangle=1,nangle
              h_ang(0,iangle,istruct)=sampintang(iangle,istruct)*0.5
              do i=1,numsampang(iangle,istruct)
                h_ang(i,iangle,istruct)=h_ang(0,iangle,istruct)+sampintang(iangle,istruct)*i
              end do
              h_ang_limit(iangle,istruct)=h_ang(numsampang(iangle,istruct),iangle,istruct)
            end do
          end if
        end do

! -- Theoretical semivariograms are now computed at the sample points.

        allocate(gamma_row_t(maxsamprow,nstructzone),gamma_row_e(maxsamprow,nstructzone),   &
                 gamma_col_t(maxsampcol,nstructzone),gamma_col_e(maxsampcol,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(num_col_e(maxsampcol,nstructzone),num_row_e(maxsamprow,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(gamma_ang_t(maxsampang,nangle,nstructzone),gamma_ang_e(maxsampang,nangle,nstructzone),  &
          stat=ierr)
          if(ierr.ne.0) go to 9200
          allocate(num_ang_e(maxsampang,nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if

        first=.true.
        do istruct=1,nstructzone
          i=structnum(istruct)
          c0=structure(i)%nugget
          PMX=structure(i)%maxpowercov
          nst=structure(i)%numvariogram
          if(nst.gt.4)then
            write(amessage,606) trim(structure(i)%structname),trim(structfile)
606         format(' Structure "',a,'" listed in file ',a,' cites more than 4 nested variograms. This is ',  &
            'not allowed.')
            go to 9890
          end if
          do ist=1,nst
            cc(ist)=structure(i)%variogram_contrib(ist)
            do iv=1,numvario
              if(vario(iv)%varname.eq.structure(i)%variogram_name(ist))then
                aa(ist)=vario(iv)%a
                it(ist)=vario(iv)%vartype
                anis(ist)=vario(iv)%anis
                ang(ist)=vario(iv)%angle
                go to 655
              end if
            end do
            write(amessage,652)
652         format(' Programming error type 1 - contact programmer.')
            go to 9890
655         continue
          end do
          y1=0.0
          y2=0.0
          x1=0.0
          x2=0.0
          covar_zero=cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,ang,anis,first)
          first=.false.
          do i=1,numsamprow(istruct)
            x2=(h_row(i,istruct)+h_row(i-1,istruct))*0.5
            gamma_row_t(i,istruct)=covar_zero-cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                          ang,anis,first)
          end do
          x2=0.0
          do i=1,numsampcol(istruct)
            y2=(h_col(i,istruct)+h_col(i-1,istruct))*0.5
            gamma_col_t(i,istruct)=covar_zero-cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                          ang,anis,first)
          end do
          if(nangle.ne.0)then
            do iangle=1,nangle
              cosangle=cos(angle(iangle)*pi/180.0d0)
              sinangle=sin(angle(iangle)*pi/180.0d0)
              do i=1,numsampang(iangle,istruct)
                h=(h_ang(i,iangle,istruct)+h_ang(i-1,iangle,istruct))*0.5
                x2=h*cosangle
                y2=h*sinangle
                gamma_ang_t(i,iangle,istruct)=covar_zero-cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                            ang,anis,first)
              end do
            end do
          end if
        end do
        write(6,657)
657     format('  - theoertical variograms computed ok.')

! -- Information is written to the start of the ZONE2VAR input file.

       if(zonfile.ne.' ')then
          write(6,1510) trim(zonfile)
1510      format(/,'  - writing header to ZON2VAR2 input file ',a,'...')
          write(zonunit) npar,mlay
          write(zonunit) (layer(ilay),ilay=1,mlay)
          write(zonunit) (ival(ipar),ipar=1,npar)
          write(zonunit) d2pfile
          do ipar=1,npar
            write(zonunit) apar(ipar)
          end do
          write(zonunit) at

          write(zonunit)  nstructzone,nangle
          write(zonunit)  maxsampcol,maxsamprow,maxsampang
          if(nangle.ne.0)then
            write(zonunit) (angle(iangle),iangle=1,nangle)
          end if

          do istruct=1,nstructzone
            write(zonunit) numsampcol(istruct),numsamprow(istruct)
            if(nangle.ne.0)then
              write(zonunit) (numsampang(iangle,istruct),iangle=1,nangle)
            end if
          end do

          do istruct=1,nstructzone
            write(zonunit) (h_row(i,istruct),i=0,numsamprow(istruct))
            write(zonunit) (h_col(i,istruct),i=0,numsampcol(istruct))
            if(nangle.ne.0)then
              do iangle=1,nangle
                write(zonunit) (h_ang(i,iangle,istruct),i=0,numsampang(iangle,istruct))
              end do
            end if
          end do

          write(zonunit) (iszone(istruct),istruct=1,nstructzone)
          write(zonunit) (structzonepar(ipar),ipar=1,npar)

          write(6,1531)
1531      format('  - header written ok.')
        end if

! -- Experimental variograms are now computed.

        write(6,658)
658     format(/,'  - computing observed semivariograms...')

        if(dflag.ne.0)then
          allocate(num_par_row(maxsamprow,npar),num_par_col(maxsampcol,npar),  &
                   deriv_row(maxsamprow,npar),deriv_col(maxsampcol,npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          num_par_row=0           ! an array
          num_par_col=0           ! an array
          deriv_row=0.0d0         ! an array
          deriv_col=0.0d0         ! an array
          if(nangle.ne.0)then
            allocate(num_par_ang(maxsampang,nangle,npar),                &
                     deriv_ang(maxsampang,nangle,npar),stat=ierr)
            if(ierr.ne.0) go to 9200
            num_par_ang=0           ! an array
            deriv_ang=0.0d0         ! an array
          end if
        end if

        gamma_row_e=0.0d0         ! an array
        gamma_col_e=0.0d0         ! an array
        num_row_e=0               ! an array
        num_col_e=0               ! an array
        if(nangle.ne.0)then
          gamma_ang_e=0.0d0       ! an array
          num_ang_e=0             ! an array
        end if

! -- If the geostatistical structure transform is log, then the log is taken of all parameter
!    values.

        if(at.eq.'l')then
          do ipar=1,npar
            if(pval(ipar).le.0.0d0)then
              write(amessage,662) trim(apar(ipar))
662           format(' Parameter "',a,'" has been provided with a zero or negative value. ', &
              'Log-transformation cannot apply to such parameters.')
              go to 9890
            end if
            pval(ipar)=log10(pval(ipar))
          end do
        end if

        do ipar=1,npar-1
          is=structzonepar(ipar)
          if((is.eq.0).or.(is.eq.-99999999)) cycle
          do istruct=1,nstructzone
            if(is.eq.iszone(istruct)) go to 720
          end do
720       continue
          e0=eastpar(ipar)
          n0=northpar(ipar)
          do jpar=ipar+1,npar
            if(structzonepar(jpar).ne.is) cycle
            ediff=abs(eastpar(jpar)-e0)
            ndiff=abs(northpar(jpar)-n0)
            if(ediff.lt.h_row(0,istruct)) go to 800
            if(ediff.gt.h_row_limit(istruct)) go to 800
            if(abs(ndiff).gt.ediff*tan_offset_angle) go to 800
            do i=1,numsamprow(istruct)
              if(ediff.le.h_row(i,istruct)) then
                rtemp=pval(ipar)-pval(jpar)
                rtemp2=rtemp*rtemp
                gamma_row_e(i,istruct)=gamma_row_e(i,istruct)+rtemp2
                num_row_e(i,istruct)=num_row_e(i,istruct)+1
                if(dflag.ne.0)then
                  num_par_row(i,ipar)=num_par_row(i,ipar)+1
                  num_par_row(i,jpar)=num_par_row(i,jpar)+1
                  deriv_row(i,ipar)=deriv_row(i,ipar)+rtemp
                  deriv_row(i,jpar)=deriv_row(i,jpar)-rtemp
                end if
                if(zflag.ne.0) write(zonunit) 1,i,istruct,ipar,jpar
                go to 810
              end if
            end do
            write(amessage,788)
788         format(' Programming error - contact programmer.')
            go to 9890
800         continue
            if(ndiff.lt.h_col(0,istruct)) go to 810
            if(ndiff.gt.h_col_limit(istruct)) go to 810
            if(abs(ediff).gt.ndiff*tan_offset_angle) go to 810
            do i=1,numsampcol(istruct)
              if(ndiff.le.h_col(i,istruct)) then
                rtemp=pval(ipar)-pval(jpar)
                rtemp2=rtemp*rtemp
                gamma_col_e(i,istruct)=gamma_col_e(i,istruct)+rtemp2
                num_col_e(i,istruct)=num_col_e(i,istruct)+1
                if(dflag.ne.0)then
                  num_par_col(i,ipar)=num_par_col(i,ipar)+1
                  num_par_col(i,jpar)=num_par_col(i,jpar)+1
                  deriv_col(i,ipar)=deriv_col(i,ipar)+rtemp
                  deriv_col(i,jpar)=deriv_col(i,jpar)-rtemp
                end if
                if(zflag.ne.0) write(zonunit) 2,i,istruct,ipar,jpar
                go to 810
              end if
            end do
            write(amessage,788)
            go to 9890
810         continue
            if(nangle.ne.0)then
              ediff=eastpar(jpar)-e0
              ndiff=northpar(jpar)-n0
              if(ndiff.lt.0.0d0)then   ! This keeps the points in the top 180 degrees of x-y coordinate plane.
                ndiff=-ndiff
                ediff=-ediff
              end if
              r=sqrt(ediff*ediff+ndiff*ndiff)
              if(ediff.eq.0.0d0)then
                theta=0.5*pi
              else
                theta=atan(ndiff/ediff)
                if(theta.lt.0.0d0)theta=theta+pi
              end if
              do iangle=1,nangle
                phi=angle(iangle)*pi/180.0d0
                alpha=theta-phi
                h=abs(r*cos(alpha))
                o=abs(r*sin(alpha))
                if(h.lt.h_ang(0,iangle,istruct)) go to 1300
                if(h.gt.h_ang_limit(iangle,istruct)) go to 1300
                if(o.gt.h*tan_offset_angle) go to 1300
                do i=1,numsampang(iangle,istruct)
                  if(h.le.h_ang(i,iangle,istruct)) then
                    rtemp=pval(ipar)-pval(jpar)
                    rtemp2=rtemp*rtemp
                    gamma_ang_e(i,iangle,istruct)=gamma_ang_e(i,iangle,istruct)+rtemp2
                    num_ang_e(i,iangle,istruct)=num_ang_e(i,iangle,istruct)+1
                    if(dflag.ne.0)then
                      num_par_ang(i,iangle,ipar)=num_par_ang(i,iangle,ipar)+1
                      num_par_ang(i,iangle,jpar)=num_par_ang(i,iangle,jpar)+1
                      deriv_ang(i,iangle,ipar)=deriv_ang(i,iangle,ipar)+rtemp
                      deriv_ang(i,iangle,jpar)=deriv_ang(i,iangle,jpar)-rtemp
                    end if
                    if(zflag.ne.0) write(zonunit) -iangle,i,istruct,ipar,jpar
                    go to 1300
                  end if
                end do
                write(amessage,788)
                go to 9890
1300            continue
              end do
            end if
          end do
        end do

        write(6,820)
820     format('  - observed semivariograms computed ok.')

! -- Observed semi-variograms are now written to a file.

        write(6,840) trim(varfile)
840     format(/,'  - writing experimental variogram file ',a,'...')
        nobs=0
        do istruct=1,nstructzone
          call num2char(iszone(istruct),anum)
          write(varunit,850) trim(anum)
850       format(//,' EXPERIMENTAL SEMIVARIOGRAMS FOR STRUCTURE ZONE CHARACTERIZED BY INTEGER VALUE OF ',A,'.')
          write(varunit,855)
855       format(/,' Row direction.')
          write(varunit,860)
860       format(/,' Separation          Gamma         Number_of_points')
          do i=1,numsamprow(istruct)
            nobs=nobs+1
            if(num_row_e(i,istruct).eq.0)then
              gamma_row_e(i,istruct)=0.0d0
            else
              gamma_row_e(i,istruct)=0.5*gamma_row_e(i,istruct)/num_row_e(i,istruct)
            end if
            write(varunit,870) 0.5*(h_row(i,istruct)+h_row(i-1,istruct)),gamma_row_e(i,istruct),  &
            num_row_e(i,istruct)
870         format(1x,1pg14.7,t20,1pg14.7,t40,i7)
          end do
          write(varunit,880)
880       format(/,' Column direction.')
          write(varunit,860)
          do i=1,numsampcol(istruct)
            nobs=nobs+1
            if(num_col_e(i,istruct).eq.0)then
              gamma_col_e(i,istruct)=0.0d0
            else
              gamma_col_e(i,istruct)=0.5*gamma_col_e(i,istruct)/num_col_e(i,istruct)
            end if
            write(varunit,870) 0.5*(h_col(i,istruct)+h_col(i-1,istruct)),gamma_col_e(i,istruct), &
            num_col_e(i,istruct)
          end do
          if(nangle.gt.0)then
            do iangle=1,nangle
              write(aangle,'(f10.1)') angle(iangle)
              aangle=adjustl(aangle)
              write(varunit,1350) trim(aangle)
1350          format(/,' Angular direction ',a,' degrees.')
              write(varunit,860)
              do i=1,numsampang(iangle,istruct)
                nobs=nobs+1
                if(num_ang_e(i,iangle,istruct).eq.0)then
                  gamma_ang_e(i,iangle,istruct)=0.0d0
                else
                  gamma_ang_e(i,iangle,istruct)=0.5*gamma_ang_e(i,iangle,istruct)/num_ang_e(i,iangle,istruct)
                end if
                write(varunit,870) 0.5*(h_ang(i,iangle,istruct)+h_ang(i-1,iangle,istruct)),  &
                gamma_ang_e(i,iangle,istruct),num_ang_e(i,iangle,istruct)
              end do
            end do
          end if
        end do
        close(unit=varunit)
        write(6,890) trim(varfile)
890     format('  - file ',a,' written ok.')

! -- If needed, observation names are computed.

        if((derfile.ne.' ').or.(pestfile.ne.' ').or.(insfile.ne.' '))then
          allocate(bobs(nobs),stat=ierr)
          if(ierr.ne.0) go to 9200
          iobs=0
          do istruct=1,nstructzone
            call num2char(iszone(istruct),azone)
            azone='z'//trim(azone)
            do i=1,numsamprow(istruct)
              iobs=iobs+1
              write(adist,'(f12.0)') (h_row(i,istruct)+h_row(i-1,istruct))*0.5
              bobs(iobs)=trim(azone)//'_'//trim(adjustl(adist))
              nb=len_trim(bobs(iobs))
              bobs(iobs)=bobs(iobs)(1:nb-1)//'_r'
            end do
            do i=1,numsampcol(istruct)
              iobs=iobs+1
              write(adist,'(f12.0)') (h_col(i,istruct)+h_col(i-1,istruct))*0.5
              bobs(iobs)=trim(azone)//'_'//trim(adjustl(adist))
              nb=len_trim(bobs(iobs))
              bobs(iobs)=bobs(iobs)(1:nb-1)//'_c'
            end do
            if(nangle.gt.0)then
              do iangle=1,nangle
                do i=1,numsampang(iangle,istruct)
                  iobs=iobs+1
                  write(adist,'(f12.0)') (h_ang(i,iangle,istruct)+h_ang(i-1,iangle,istruct))*0.5
                  bobs(iobs)=trim(azone)//'_'//trim(adjustl(adist))
                  write(aangle,'(f10.0)') angle(iangle)
                  aangle=adjustl(aangle)
                  nb=len_trim(aangle)
                  aangle=aangle(1:nb-1)
                  nb=len_trim(bobs(iobs))
                  bobs(iobs)=bobs(iobs)(1:nb-1)//'_'//trim(aangle)
                end do
              end do
            end if
          end do
        end if

        if(zflag.ne.0)then
          write(zonunit) 0,0,0,0,0
          do istruct=1,nstructzone
            write(zonunit) (num_row_e(i,istruct),i=1,maxsamprow)
            write(zonunit) (num_col_e(i,istruct),i=1,maxsampcol)
            if(nangle.ne.0)then
              do iangle=1,nangle
                write(zonunit) (num_ang_e(i,iangle,istruct),i=1,maxsampang)
              end do
            end if
          end do
          if(dflag.ne.0)then
            write(zonunit) ((num_par_row(i,ipar),i=1,maxsamprow),ipar=1,npar)
            write(zonunit) ((num_par_col(i,ipar),i=1,maxsampcol),ipar=1,npar)
            write(zonunit) (((num_par_ang(i,iangle,ipar),i=1,maxsampang),iangle=1,nangle),ipar=1,npar)
          end if
          write(zonunit) nobs
          do iobs=1,nobs
            write(zonunit) bobs(iobs)
          end do
        end if

! -- The instruction file to read the semi-variogram file is now written.

        if(insfile.ne.' ')then
          iobs=0
          write(6,895) trim(insfile)
895       format(/,'  - writing instruction file ',a,'...')
          write(insunit,896)
896       format('pif $')
          do istruct=1,nstructzone
            write(insunit,900)
900         format('l7')
            do i=1,numsamprow(istruct)
              iobs=iobs+1
              write(insunit,910) trim(bobs(iobs))
910           format('l1 ['a,']20:36')
            end do
            write(insunit,920)
920         format('l4')
            do i=1,numsampcol(istruct)
              iobs=iobs+1
              write(insunit,910) trim(bobs(iobs))
            end do
            if(nangle.ne.0)then
              do iangle=1,nangle
                write(insunit,920)
                do i=1,numsampang(iangle,istruct)
                  iobs=iobs+1
                  write(insunit,910) trim(bobs(iobs))
                end do
              end do
            end if
          end do
          close(unit=insunit)
          write(6,921) trim(insfile)
921       format('  - file ',a,' written ok.')
        end if

! -- The PEST building block file is now written.

        if(pestfile.ne.' ')then
          write(6,932) trim(pestfile)
932       format(/,'  - writing PEST building block file ',a,'...')
          write(pestunit,940)
940       format('* observation groups')
          do istruct=1,nstructzone
            call num2char(iszone(istruct),azone)
            azone=adjustl(azone)
            agroup='regulz'//trim(azone)//'_r'
            write(pestunit,950) trim(agroup)
950         format(a)
            agroup='regulz'//trim(azone)//'_c'
            write(pestunit,950) trim(agroup)
            if(nangle.ne.0)then
              do iangle=1,nangle
                write(aangle,'(f10.0)') angle(iangle)
                aangle=adjustl(aangle)
                nb=len_trim(aangle)
                aangle=aangle(1:nb-1)
                agroup='regulz'//trim(azone)//'_'//trim(aangle)
                write(pestunit,950) trim(agroup)
              end do
            end if
          end do
          write(pestunit,960)
960       format('* observation data')
          iobs=0
          do istruct=1,nstructzone
            call num2char(iszone(istruct),azone)
            azone=adjustl(azone)
            agroup='regulz'//trim(azone)//'_r'
            do i=1,numsamprow(istruct)
              iobs=iobs+1
              if(num_row_e(i,istruct).eq.0)then
                weight=0.0d0
              else
                weight=1.0d0/sqrt(float(num_row_e(i,istruct)))
              end if
              write(pestunit,970) trim(bobs(iobs)),gamma_row_t(i,istruct),weight,trim(agroup)
970           format(a,t15,1pg14.7,3x,1pg14.7,3x,a)
            end do
            agroup='regulz'//trim(azone)//'_c'
            do i=1,numsampcol(istruct)
              iobs=iobs+1
              if(num_col_e(i,istruct).eq.0)then
                weight=0.0d0
              else
                weight=1.0d0/sqrt(float(num_col_e(i,istruct)))
              end if
              write(pestunit,970) trim(bobs(iobs)),gamma_col_t(i,istruct),weight,trim(agroup)
            end do
            if(nangle.ne.0)then
              do iangle=1,nangle
                write(aangle,'(f10.0)') angle(iangle)
                aangle=adjustl(aangle)
                nb=len_trim(aangle)
                aangle=aangle(1:nb-1)
                agroup='regulz'//trim(azone)//'_'//trim(aangle)
                do i=1,numsampang(iangle,istruct)
                  iobs=iobs+1
                  if(num_ang_e(i,iangle,istruct).eq.0)then
                    weight=0.0d0
                  else
                    weight=1.0d0/sqrt(float(num_ang_e(i,iangle,istruct)))
                  end if
                  write(pestunit,970) trim(bobs(iobs)),gamma_ang_t(i,iangle,istruct),weight,trim(agroup)
                end do
              end do
            end if
          end do
          write(pestunit,990)
990       format('* model input/output')
          call addquote(insfile,afile)
          call addquote(varfile,bfile)
          write(pestunit,1000) trim(afile),trim(bfile)
1000      format(a,2x,a)
          close(unit=pestunit)
          write(6,921) trim(pestfile)
        end if

! -- Now derivatives are computed. Derivatives will be written in matrix file format.

!        write(99,*)                                          !debug
!        write(99,*) ' Parameter   structzonepar.'           !debug
!        do ipar=1,npar                                       !debug
!          write(99,*) trim(apar(ipar)), structzonepar(ipar)  !debug
!        end do                                               !debug

!        write(99,*)                                                                    !debug
!        write(99,*) ' Last observation'                                                !debug
!        do ipar=1,npar                                                                 !debug
!          write(99,*) trim(apar(ipar)),deriv_col(numsampcol(nstructzone),ipar)  !debug
!        end do                                                                         !debug



        if(derfile.ne.' ')then
          write(6,1010) trim(derfile)
1010      format(/,'  - writing derivatives matrix file ',a,'...')
          allocate(x_row(npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          write(derunit,1020) nobs,npar,2
          do istruct=1,nstructzone
            is=iszone(istruct)
            do i=1,numsamprow(istruct)
              do ipar=1,npar
                if(structzonepar(ipar).ne.is)then
                  x_row(ipar)=0.0d0
                else
                  itemp=num_par_row(i,ipar)
                  if(itemp.eq.0)then
                    x_row(ipar)=0.0d0
                  else
                    x_row(ipar)=deriv_row(i,ipar)/num_row_e(i,istruct)
                    if(at.eq.'l')then
                      x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                    end if
                  end if
                end if
              end do
              write(derunit,1020) (x_row(ipar),ipar=1,npar)
1020          format(8(1x,1pg14.7))
            end do
            do i=1,numsampcol(istruct)
              do ipar=1,npar
                if(structzonepar(ipar).ne.is)then
                  x_row(ipar)=0.0d0
                else
                  itemp=num_par_col(i,ipar)
                  if(itemp.eq.0)then
                    x_row(ipar)=0.0d0
                  else
                    x_row(ipar)=deriv_col(i,ipar)/num_col_e(i,istruct)
                    if(at.eq.'l')then
                      x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                    end if
                  end if
                end if
              end do
              write(derunit,1020) (x_row(ipar),ipar=1,npar)
            end do
            if(nangle.ne.0)then
              do iangle=1,nangle
                do i=1,numsampang(iangle,istruct)
                  do ipar=1,npar
                    if(structzonepar(ipar).ne.is)then
                      x_row(ipar)=0.0d0
                    else
                      itemp=num_par_ang(i,iangle,ipar)
                      if(itemp.eq.0)then
                        x_row(ipar)=0.0d0
                      else
                        x_row(ipar)=deriv_ang(i,iangle,ipar)/num_ang_e(i,iangle,istruct)
                        if(at.eq.'l')then
                          x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                        end if
                      end if
                    end if
                  end do
                  write(derunit,1020) (x_row(ipar),ipar=1,npar)
                end do
              end do
            end if
          end do
          write(derunit,1030)
1030      format('* row names')
          do iobs=1,nobs
            write(derunit,1035) trim(bobs(iobs))
1035        format(1x,a)
          end do
          write(derunit,1040)
1040      format('* column names')
          do ipar=1,npar
            write(derunit,1035) trim(apar(ipar))
          end do
          close(unit=derunit)
          write(6,921) trim(derfile)
       end if

       if(zonfile.ne.' ')then
         close(unit=zonunit)
         write(6,*)
         write(6,1041) trim(zonfile)
1041     format('  - ZONE2VAR2 input file ',a,' written ok.')
       end if

       go to 9900

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(d2pfile)
9010    format(' Error reading line ',a,' of distributed-to-PEST parameter file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(d2pfile)
9060    format(' Unexpected end encountered to distributed-to-PEST parameter file ',a,'.')
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(ircfile)
9110    format(' Error encountered in reading line ',a,' of integer-real ',  &
        'correspondence file ',a,'.')
        go to 9890



9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890


9890    call write_message(leadspace='yes')

9900    call close_files

        deallocate(ival,layer,structnum,iarray,isarray,numsampcol,numsamprow,numcellpar,  &
                   structzonepar,num_col_e,num_row_e,num_par_row,num_par_col,stat=ierr)

        deallocate(pval,sampintcol,sampintrow,east,north,eastpar,northpar,h_col,       &
                   h_row,h_row_limit,h_col_limit,gamma_row_t,gamma_row_e,gamma_col_t,  &
                   gamma_col_e,deriv_row,deriv_col,x_row,pvalkeep,stat=ierr)

        deallocate(apar,bobs,stat=ierr)

        deallocate(structure,vario,stat=ierr)

        deallocate(numsampang,num_ang_e,stat=ierr)
        deallocate(angle,sampintang,h_ang_limit,h_ang,gamma_ang_t,gamma_ang_e,stat=ierr)
        deallocate(num_par_ang,stat=ierr)
        deallocate(deriv_ang,stat=ierr)

	call free_grid_mem(gridspec)

end program zone2var1



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



 subroutine replace_character(astring,char1,char2)

 ! -- Subroutine replace_characer replaces char1 in string astring with char2.

      implicit none

      character*(*) astring
      character*1 char1,char2

      integer nb,i

      nb=len_trim(astring)
      if(nb.eq.0) return
      do i=1,nb
        if(astring(i:i).eq.char1) astring(i:i)=char2
      end do

      return

end subroutine replace_character



      real function cova2(x1,y1,x2,y2,nst,c0,PMX,cc,aa,it,  &
                          ang,anis,first)                             !jd

      implicit integer(i-n), real(a-h, o-z)                           !jd

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





