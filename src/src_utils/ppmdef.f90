program ppmdef

! -- Program PPMDEF builds a distributed-to-PEST parameter file for the use of ASENPROC where PEST parameters
!    are pilot points and MODFLOW employs distributed parameters through its adjoint solver.

	use defn
	use inter

	implicit none

        logical        :: lexist
        integer        :: nlay,nrow,ncol,nper,ixsec,ifrefm,nplpf,nclu,distribnclu,nzn,nml,  &
                          newnml,ii
        integer        :: ilay,irow,icol,icount,ifunction,mcol,mrow,npp,maxlen,ipp,nb,nbp, &
                          ndim,maxent,icellno,mlay,itrans,na,iconstant,itransold
        integer        :: ifail,ierr,iflag,ioc
        integer        :: iline,iunit,n,k,kk,i,j,itemp1,itemp2
        integer        :: idis,ibas,ilpf,izone,imult
        integer        :: locat_dis,locat_bas,locat_lpf,locat_zone,locat_mult
        integer        :: nameunit,facunit,defnunit,pestunit,disunit,basunit,lpfunit,multunit,  &
                          zoneunit,outunit,ppunit,templateunit,delunit
        integer        :: lpfunitkeep,multunitkeep
        integer        :: lpfkeepflag,multkeepflag,messageflag
        real           :: lbound,ubound,initval,rtemp
        real           :: distribparval
        double precision :: easting,northing

        character*1    :: facformat,aa
        character*10   :: aline,atype,arraytype,distribpartype,tempname,newmltarr,anum,aname
        character*20   :: distribpar,parname,zonename,tempmltarr,multname,atemp
        character*24   :: bname
        character*120  :: aprompt
        character*120  :: ppfile,intfile,atempf
        character*200  :: namefile,facfile,defnfile,pestfile,multfile,outfile,newmultfile,   &
                          disfile,basfile,lpffile,zonefile,facinfile,templatefile,modinfile, &
                          external_file
        character*200  :: lpffilekeep,multfilekeep
        character*200  :: bfile,afile

        integer, parameter :: MAXCLU=100
        integer            :: layer(MAXCLU),iz(10,MAXCLU)
        character*10       :: mltarr(MAXCLU),zonarr(MAXCLU)

        integer, allocatable      :: laycbd(:),iarray(:,:,:),itemparray(:,:),ipt(:)
        integer, allocatable      :: ivec(:),discard(:),ifound(:),ilaycount(:)
        real, allocatable         :: rvec(:),wt(:),rarray(:,:)
        character*12, allocatable :: ppname(:),ppname1(:)


        write(amessage,5)
5       format(' Program PPMDEF writes a distributed-to-PEST parameter file for pilot point ',  &
        'parameters for the subsequent use of ASENPROC.')
        call write_message(leadspace='yes',endspace='yes')


        lpfkeepflag=0
        multkeepflag=0
        messageflag=0

! -- The name of the name file is acquired.

100     continue
110     call open_input_file(ifail,          &
        ' Enter name of MODFLOW name file: ',namefile,nameunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

! -- The name of the krigging factor file is acquired.

        write(6,*)
150     aprompt=' Enter name of interpolation factor file: '
        call open_input_file(ifail,aprompt,facfile,facunit,form_prompt='yes', &
        fformat=facformat)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          close(unit=nameunit)
          write(6,*)
          go to 100
        end if
160     write(6,170,advance='no')
170     format(' To what MODFLOW distributed parameter does this pertain? ')
        read(5,'(a)',err=160) distribpar
        call casetrans(distribpar,'lo')
        if(distribpar(1:2).eq.'e ')then
          close(unit=facunit)
          write(6,*)
          go to 150
        end if
        call remquote(distribpar)

        write(6,*)
190     write(6,200,advance='no')
200     format(' Enter lower bound for pilot point parameters: ')
        if(key_read(lbound).ne.0) go to 190
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 160
        end if
220     write(6,230,advance='no')
230     format(' Enter upper bound for pilot point parameters: ')
        if(key_read(ubound).ne.0) go to 220
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 190
        end if
        if(ubound.le.lbound)then
          write(6,232)
232       format(' This must exceed lower bound - try again.')
          go to 220
        end if
235     write(6,236,advance='no')
236     format(' Enter initial value for pilot point parameters: ')
        if(key_read(initval).ne.0) go to 235
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 220
        end if
        if((initval.lt.lbound).or.(initval.gt.ubound))then
          write(6,237)
237       format(' This must lie between lower and upper bounds - try again.')
          go to 235
        end if

! -- The names of output files are now acquired.

        write(6,*)
250     call open_output_file(ifail, &
	' Enter name for distribted-to-PEST parameter file: ',defnfile,defnunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 235
	end if

270     call open_output_file(ifail, &
	' Enter name for PEST building block file: ',pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
          close(unit=defnunit)
	  write(6,*)
	  go to 250
	end if

! -- The name file is read.

        write(6,*)
        write(6,280) trim(namefile)
280     format(' - reading MODFLOW name file ',a,'...')
        iline=0
        idis=0
        ibas=0
        ilpf=0
        izone=0
        imult=0
        do
          iline=iline+1
          read(nameunit,'(a)',end=400) cline
          cline=adjustl(cline)
          if(cline(1:1).eq.'#') cycle
          if(cline.eq.' ') cycle
          call linesplit(ifail,3)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,320) trim(aline),trim(namefile)
320         format(' Insufficient entires on line ',a,' of MODFLOW name file ',a,'.')
            go to 9890
          end if
          atype=cline(left_word(1):right_word(1))
          call casetrans(atype,'lo')
          iunit=char2int(ifail,2)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,340) trim(aline),trim(namefile)
340         format(' Cannot read second entry on line ',a,  &
            ' of file ',a,'.')
            go to 9890
          end if
          bfile=cline(left_word(3):)
          aa=bfile(1:1)
          if((aa.eq.'"').or.(aa.eq.''''))then
            continue
          else
            aa=' '
          end if
          n=index(bfile(2:),aa)
          if(n.eq.0)then
            call num2char(iline,aline)
            write(amessage,350) trim(aline),trim(namefile)
350         format(' Improper filename at line ',a,' of name file ',a,'.')
            go to 9890
          end if
          if(aa.eq.' ')then
            bfile=bfile(1:n)
          else
            bfile=bfile(2:n)
          end if
          call casetrans(bfile,'lo')
          if(atype(1:3).eq.'dis')then
            if(idis.ne.0)then
              write(amessage,360) trim(namefile)
360           format(' More than one discretisation file is named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            disfile=bfile
            locat_dis=iunit
            idis=1
          else if(atype(1:3).eq.'bas')then
            if(ibas.ne.0)then
              write(amessage,370) trim(namefile)
370           format(' More than one BASIC package input file is named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            basfile=bfile
            locat_bas=iunit
            ibas=1
          else if(atype(1:3).eq.'lpf')then
            if(ilpf.ne.0)then
              write(amessage,380) trim(namefile)
380           format(' More than one LPF package input file is named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            lpffile=bfile
            locat_lpf=iunit
            ilpf=1
          else if(atype(1:4).eq.'zone')then
            if(izone.ne.0)then
              write(amessage,381) trim(namefile)
381           format(' More than one zone file is named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            zonefile=bfile
            locat_zone=iunit
            izone=1
          else if(atype(1:4).eq.'mult')then
            if(imult.ne.0)then
              write(amessage,382) trim(namefile)
382           format(' More than one multiplier file is named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            multfile=bfile
            locat_mult=iunit
            imult=1
          end if
        end do
400     continue
        if(ilpf.eq.0)then
          write(amessage,405) trim(namefile)
405       format(' No LPF file is cited in MODFLOW name file ',a,'; PPMDEF can presently only work with ',   &
          'LPF parameters.')
          go to 9890
        end if
        if(idis.eq.0)then
          write(amessage,406) trim(namefile)
406       format(' No discretisation file is cited in MODFLOW name file ',a,'.')
          go to 9890
        end if
        if(ibas.eq.0)then
          write(amessage,407) trim(namefile)
407       format(' No BASIC package input file is cited in MODFLOW name file ',a,'.')
          go to 9890
        end if
        if(imult.eq.0)then
          write(amessage,408) trim(namefile)
408       format(' No multiplier file is cited in MODFLOW name file ',a,'. A multiplier file ', &
          'must be cited in the name file, even if it is empty.')
          go to 9890
        end if

        close(unit=nameunit)
        write(6,410) trim(namefile)
410     format(' - file ',a,' read ok.')

! -- The discretisation file is now read in order to establish the dimensions of the grid.

        write(6,*)
        write(6,420) trim(disfile)
420     format(' - reading MODFLOW discretisation file ',a,'...')
        disunit=nextunit()
        open(unit=disunit,file=disfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,425) trim(disfile),trim(namefile)
425       format(' Cannot open discretization file ',a,' cited in MODFLOW ',  &
          'name file ',a,'.')
          go to 9890
        end if
        iline=0
        do
          iline=iline+1
          read(disunit,'(a)',err=9000,end=9050) cline
          cline=adjustl(cline)
          if(cline(1:1).ne.'#') exit
        end do
        call linesplit(ifail,4)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,426) trim(aline),trim(disfile)
426       format(' Insufficient entries on line ',a,' of MODFLOW discretisation file ',a,'.')
          go to 9890
        end if
        nlay=char2int(ifail,1)
        if(ifail.ne.0) go to 9070
        nrow=char2int(ifail,2)
        if(ifail.ne.0) go to 9070
        ncol=char2int(ifail,3)
        if(ifail.ne.0) go to 9070
        nper=char2int(ifail,4)
        if(ifail.ne.0) go to 9070
        if((nrow.le.0).or.(ncol.le.0).or.(nlay.le.0).or.(nper.le.0)) go to 9070
        allocate(laycbd(nlay),stat=ierr)
        if(ierr.ne.0) go to 9200
        read(disunit,*,iostat=ierr) (laycbd(i),i=1,nlay)
        if(ierr.ne.0)then
          write(amessage,440) trim(disfile)
440       format(' Error reading LAYCBD array from MODFLOW discretisation file ',a,'.')
          go to 9890
        end if
        close(unit=disunit)
        write(6,410) trim(disfile)

! -- The basic package input file is now read to obtain ibound arrays.

        allocate(iarray(ncol,nrow,nlay),rarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200

        write(6,*)
        write(6,500) trim(basfile)
500     format(' - reading MODFLOW BASIC package input file ',a,'...')
        basunit=nextunit()
        open(unit=basunit,file=basfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,510) trim(basfile),trim(namefile)
510       format(' Cannot open BASIC package input file ',a,' cited in MODFLOW ',  &
          'name file ',a,'.')
          go to 9890
        end if
        iline=0
        do
          iline=iline+1
          read(basunit,'(a)',err=9100,end=9150) cline
          cline=adjustl(cline)
          if(cline(1:1).ne.'#') exit
        end do
        ixsec=0
        ifrefm=0
        call casetrans(cline,'lo')
        if(index(cline,'xsection').ne.0) ixsec=1
        if(index(cline,'free').ne.0) ifrefm=1
        if(ixsec.ne.0)then
          write(amessage,520)
520       format(' MODFLOW model is running in cross-section mode; this is not permitted with ',  &
          'present version of PPMDEF.')
          go to 9890
        end if
        bname='BOUNDARY ARRAY          '
        do k=1,nlay
          kk=k
          call u2dint(ifail,iarray(1,1,kk),bname,nrow,ncol,kk,basunit,0,locat_bas)
          if(ifail.ne.0) go to 9890
        end do
        close(unit=basunit)
        write(6,410) trim(basfile)

! -- The LPF input file is now read to obtain data on the distributed parameter of interest.

        allocate(ivec(nlay),rvec(nlay),stat=ierr)
        if(ierr.ne.0) go to 9200

        write(6,*)
        write(6,550) trim(lpffile)
550     format(' - reading MODFLOW LPF package input file ',a,'...')
        lpfunit=nextunit()
        open(unit=lpfunit,file=lpffile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,560) trim(lpffile),trim(namefile)
560       format(' Cannot open LPF package input file ',a,' cited in MODFLOW ',  &
          'name file ',a,'.')
          go to 9890
        end if
        iline=0
        do
          iline=iline+1
          read(lpfunit,'(a)',err=9300,end=9350) cline
          cline=adjustl(cline)
          if(cline(1:1).ne.'#') exit
        end do
        call linesplit(ifail,3)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,570) trim(aline),trim(lpffile)
570       format(' Insufficient entries on line ',a,' of MODFLOW LPF package input file ',a,'.')
          go to 9890
        end if
        nplpf=char2int(ifail,3)
        if(ifail.ne.0) go to 9370
        if(nplpf.eq.0)then
          write(amessage,580) trim(lpffile)
580       format(' No MODFLOW parameters are cited in LPF package input file ',a,'.')
          go to 9890
        end if
        read(lpfunit,*,iostat=ierr) (ivec(ilay),ilay=1,nlay)
        if(ierr.ne.0)then
          arraytype='LAYTYP'
          go to 9390
        end if
        read(lpfunit,*,iostat=ierr) (ivec(ilay),ilay=1,nlay)
        if(ierr.ne.0)then
          arraytype='LAYAVG'
          go to 9390
        end if
        read(lpfunit,*,iostat=ierr) (rvec(ilay),ilay=1,nlay)
        if(ierr.ne.0)then
          arraytype='CHANI'
          go to 9390
        end if
        read(lpfunit,*,iostat=ierr) (ivec(ilay),ilay=1,nlay)
        if(ierr.ne.0)then
          arraytype='LAYVKA'
          go to 9390
        end if
        read(lpfunit,*,iostat=ierr) (ivec(ilay),ilay=1,nlay)
        if(ierr.ne.0)then
          arraytype='LAYWET'
          go to 9390
        end if
        do ilay=1,nlay
          if(ivec(ilay).ne.0) go to 590
        end do
        go to 600
590     read(lpfunit,'(a)',iostat=ierr) cline
        if(ierr.ne.0) then
          write(amessage,592) trim(lpffile)
592       format(' Error reading "[WETFCT IWETIT IHDWET]" line from MODFLOW LPF input file ',a,'.')
          go to 9890
        end if
        call linesplit(ifail,3)
        if(ifail.ne.0)then
          write(amessage,592) trim(lpffile)
          go to 9890
        end if
600     continue

! -- Now we look for data on the chosen parameter.

        do i=1,nplpf
          read(lpfunit,'(a)',err=9400,end=9420) cline
          if(cline.eq.' ') cycle
          call linesplit(ifail,4)
          parname=cline(left_word(1):right_word(1))
          call casetrans(parname,'lo')
          if(ifail.ne.0)then
            write(amessage,610) trim(parname),trim(lpffile)
610         format(' Lead line in definition of MODFLOW parameter "',a,'" has too few ',  &
            'entries in LPF package input file ',a,'.')
            go to 9890
          end if
          nclu=char2int(ifail,4)
          if(ifail.ne.0) then
            write(amessage,620) trim(parname),trim(lpffile)
620         format(' Cannot read NCLU variable for parameter "',a,'" from MODFLOW LPF ',  &
            'package input file ',a,'.')
            go to 9890
          end if
          if(parname.eq.distribpar) go to 700
          do j=1,nclu
630         continue
            read(lpfunit,'(a)',iostat=ierr) cline
            if(ierr.ne.0)then
              write(amessage,640) trim(parname),trim(lpffile)
640           format(' Error reading data for parameter "',a,'" from MODFLOW LPF package ', &
              'input file ',a,'.')
              go to 9890
            end if
            if(cline.eq.' ') go to 630
          end do
        end do
        write(amessage,650) trim(distribpar),trim(lpffile)
650     format(' Distributed parameter "',a,'" is not cited in MODFLOW LPF package input file ',a,'.')
        go to 9890

! -- Data pertaining to the distributed parameter of interest is now read.

700     distribpartype=cline(left_word(2):right_word(2))
        call casetrans(distribpartype,'lo')
        if((distribpartype.ne.'hk').and.                 &
           (distribpartype.ne.'hani').and.               &
           (distribpartype.ne.'vk').and.                 &
           (distribpartype.ne.'vani').and.               &
           (distribpartype.ne.'ss').and.                 &
           (distribpartype.ne.'sy').and.                 &
           (distribpartype.ne.'vkcb')) then
           write(amessage,710) trim(distribpar),trim(distribpartype),trim(lpffile)
710        format(' Parameter "',a,'" is an unknown parameter type - "',a,'" - in ',   &
           'MODFLOW LPF input file ',a,'.')
           go to 9890
        end if
        distribparval=char2real(ifail,3)
        if(ifail.ne.0)then
          write(amessage,720) trim(distribpar),trim(lpffile)
720       format(' Cannot read value assigned to parameter "',a,'" in LPF package input file ',a,'.')
          go to 9890
        end if
        if(nclu.gt.MAXCLU)then
          write(amessage,730)
730       format(' Increase MAXCLU and re-compile program.')
          go to 9890
        end if
        do i=1,nclu
          read(lpfunit,'(a)',err=9430,end=9450) cline
          call linesplit(ifail,3)
          if(ifail.ne.0) go to 9430
          layer(i)=char2int(ifail,1)
          if(ifail.ne.0) go to 9430
          mltarr(i)=cline(left_word(2):right_word(2))
          call casetrans(mltarr(i),'lo')
          call remquote(mltarr(i))
          if(mltarr(i).eq.'none')then
            write(amessage,731) trim(distribpar),trim(lpffile)
731         format(' PPMDEF does not permit "none" as a multiplier array name for a distributed parameter. ',  &
            'See parameter "',a,'" in LPF package input file ',a,'.')
            go to 9890
          end if
          zonarr(i)=cline(left_word(3):right_word(3))
          call casetrans(zonarr(i),'lo')
          call remquote(zonarr(i))
          if(zonarr(i).ne.'all')then
            cline=cline(right_word(3)+1:)
            cline=adjustl(cline)
            if(cline.eq.' ') go to 9430
            do j=10,1,-1
              call linesplit(ifail,j)
              if(ifail.eq.0) exit
              iz(j,i)=0
            end do
            do k=1,j
              iz(k,i)=char2int(ifail,k)
              if(ifail.ne.0) go to 9430
            end do
          end if
        end do
        close(unit=lpfunit)
        write(6,410) trim(lpffile)

! -- Before reading zone arrays, the integer ibound arrays are adjusted.

        do ilay=1,nlay
          do irow=1,nrow
            do icol=1,ncol
              if(iarray(icol,irow,ilay).ne.0) iarray(icol,irow,ilay)=1
            end do
          end do
        end do

! -- Now we read the zone array file if necessary.

        allocate(itemparray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200

        distribnclu=nclu
        do i=1,distribnclu
          if(zonarr(i).ne.'all') go to 750
        end do
        go to 795
750     continue
        if(izone.eq.0)then
          write(amessage,735) trim(lpffile),trim(namefile)
735       format(' At least one zone name is cited in MODFLOW LPF input file ',a,' yet no zone file ',  &
          'is cited in the MODFLOW name file ',a,'.')
          go to 9890
        end if
        write(6,*)
        write(6,740) trim(zonefile)
740     format(' - reading MODFLOW zone file ',a,'...')
        zoneunit=nextunit()
        open(unit=zoneunit,file=zonefile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,745) trim(zonefile),trim(namefile)
745       format(' Cannot open zone file ',a,' cited in MODFLOW ',  &
          'name file ',a,'.')
          go to 9890
        end if
        iline=0
        do
          iline=iline+1
          read(zoneunit,'(a)',err=9500,end=9520) cline
          cline=adjustl(cline)
          if(cline.eq.' ')cycle
          if(cline(1:1).ne.'#') exit
        end do
        call linesplit(ifail,1)
        nzn=char2int(ifail,1)
        if(ifail.ne.0) then
          write(amessage,770) trim(zonefile)
770       format(' Error reading NZN variable from MODFLOW zone file ',a,'.')
          go to 9890
        end if
        do izone=1,nzn
780       continue
          read(zoneunit,'(a)',err=9500,end=9520) cline
          if(cline.eq.' ') go to 780
          call linesplit(ifail,1)
          zonename=cline(left_word(1):right_word(1))
          call casetrans(zonename,'lo')
          if(zonename.eq.'all')then
            write(amessage,790) trim(zonefile)
790         format( '"ALL" is an illegal zone name in file ',a,'.')
            go to 9890
          end if
          bname='"'//trim(zonename)//'" ZONE ARRAY'
          call u2dint(ifail,itemparray,bname,nrow,ncol,0,zoneunit,0,locat_zone)
          if(ifail.ne.0) go to 9890
          do i=1,distribnclu
            if(zonarr(i).eq.zonename)then
              ilay=layer(i)
              do k=1,10
                if(iz(k,i).eq.0) cycle
                do irow=1,nrow
                  do icol=1,ncol
                    if(itemparray(icol,irow).eq.iz(k,i))iarray(icol,irow,ilay)=2
                  end do
                end do
              end do
            end if
          end do
        end do
        do i=1,distribnclu
          if(zonarr(i).eq.'all')then
            ilay=layer(i)
            do irow=1,nrow
              do icol=1,ncol
                if(iarray(icol,irow,ilay).ne.0) iarray(icol,irow,ilay)=2
              end do
            end do
          end if
        end do
        do ilay=1,nlay
          do irow=1,nrow
            do icol=1,ncol
              if(iarray(icol,irow,ilay).eq.2)then
                iarray(icol,irow,ilay)=1
              else
                iarray(icol,irow,ilay)=0
              end if
            end do
          end do
        end do
        close(unit=zoneunit)
        write(6,410) trim(zonefile)
795     continue

! -- So now our integer arrays define completely where interpolation is to take place.
! -- The next task is to re-write the LPF input file with a single multiplier array cited
!    for all incidences of the distributed parameter of interest. At the same time we
!    establish whether we can remove any multiplier arrays from the multiplier file.

! -- First we copy the existing LPF file to a new one.


        ii=0
        do
          ii=ii+1
          call num2char(ii,anum)
          anum=adjustl(anum)
          lpffilekeep=trim(lpffile)//'.k'//trim(anum)
          inquire(file=lpffilekeep,exist=lexist)
          if(.not.lexist)exit
        end do
        write(6,*)
        write(6,800) trim(lpffile),trim(lpffilekeep)
800     format(' - copying file ',a,' to ',a,' for safekeeping...')
        lpfunit=nextunit()
        open(unit=lpfunit,file=lpffile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,810) trim(lpffile)
810       format(' Cannot re-open file ',a,'.')
          go to 9890
        end if
        lpfunitkeep=nextunit()
        open(unit=lpfunitkeep,file=lpffilekeep)
        do
          read(lpfunit,'(a)',end=830) cline
          write(lpfunitkeep,'(a)') trim(cline)
        end do
830     continue
        close(unit=lpfunit)
        close(unit=lpfunitkeep)
        write(6,840) trim(lpffile)
840     format(' - file ',a,' copied ok.')
        lpfkeepflag=1

! -- Now we re-write the old one citing only the new multiplier file.

        write(6,*)
        write(6,850) trim(lpffile)
850     format(' - re-writing file ',a,' citing new multiplier array...')
        lpfunitkeep=nextunit()
        open(unit=lpfunitkeep,file=lpffilekeep,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,852) trim(lpffilekeep)
852       format(' Cannot re-open file ',a,'.')
          go to 9890
        end if
        lpfunit=nextunit()
        open(unit=lpfunit,file=lpffile)
        do
          read(lpfunitkeep,'(a)',err=9600,end=9600) cline
          write(lpfunit,'(a)') trim(cline)
          cline=adjustl(cline)
          if(cline(1:1).ne.'#') exit
        end do
        read(lpfunitkeep,*,err=9600) (ivec(ilay),ilay=1,nlay)
        write(lpfunit,860) (ivec(ilay),ilay=1,nlay)
860     format(100(i3))
        read(lpfunitkeep,*,err=9600) (ivec(ilay),ilay=1,nlay)
        write(lpfunit,860) (ivec(ilay),ilay=1,nlay)
        read(lpfunitkeep,*,err=9600) (rvec(ilay),ilay=1,nlay)
        write(lpfunit,870) (rvec(ilay),ilay=1,nlay)
870     format(8(1pg14.7,1x))
        read(lpfunitkeep,*,err=9600) (ivec(ilay),ilay=1,nlay)
        write(lpfunit,860) (ivec(ilay),ilay=1,nlay)
        read(lpfunitkeep,*,err=9600) (ivec(ilay),ilay=1,nlay)
        write(lpfunit,860) (ivec(ilay),ilay=1,nlay)
        do ilay=1,nlay
          if(ivec(ilay).ne.0) go to 890
        end do
        go to 900
890     read(lpfunitkeep,'(a)',err=9600) cline
        write(lpfunit,'(a)') trim(cline)
900     continue

        allocate(discard(distribnclu),ifound(distribnclu),stat=ierr)
        if(ierr.ne.0) go to 9200
        discard=0          ! an array
        ifound=0           ! an array
        newmltarr=trim(distribpar(1:8))

        do i=1,nplpf
          read(lpfunitkeep,'(a)',err=9600) cline
          if(cline.eq.' ') then
            write(lpfunit,*)
            cycle
          end if
          call linesplit(ifail,4)
          parname=cline(left_word(1):right_word(1))
          call casetrans(parname,'lo')
          nclu=char2int(ifail,4)
          if(ifail.ne.0)then
            write(amessage,905) trim(parname),trim(lpffilekeep)
905         format(' Cannot read NCLU for parameter "',a,'" from file ',a,'.')
            go to 9890
          end if
          if(parname.eq.distribpar) then
            write(lpfunit,906) cline(left_word(1):right_word(2)),cline(left_word(4):right_word(4))
906         format(a,' 1.0  ',a)
            do j=1,nclu
              read(lpfunitkeep,'(a)',err=9470,end=9490) cline
              call linesplit(ifail,3)
              write(lpfunit,920) cline(left_word(1):right_word(1)),                &
              trim(newmltarr),trim(cline(left_word(3):))
920           format(a,2x,a,2x,a)
              discard(j)=1
            end do
          else
            write(lpfunit,'(a)') trim(cline)
            do j=1,nclu
930           continue
              read(lpfunitkeep,'(a)',iostat=ierr) cline
              if(ierr.ne.0)then
                write(amessage,932) trim(parname),trim(lpffilekeep)
932             format(' Error reading data for parameter "',a,'" from LPF package ',  &
                'input file ',a,'.')
                go to 9890
              end if
              if(cline.eq.' ') then
                write(lpfunit,*)
                go to 930
              end if
              write(lpfunit,'(a)') trim(cline)
              call linesplit(ifail,3)
              tempmltarr=cline(left_word(2):right_word(2))
              call casetrans(tempmltarr,'lo')
              do k=1,distribnclu
                if(tempmltarr.eq.mltarr(k)) discard(k)=0
              end do
            end do
          end if
        end do
        do
          read(lpfunitkeep,'(a)',err=9600,end=950) cline
          write(lpfunit,'(a)') trim(cline)
        end do
950     continue
        write(6,960) trim(lpffile)
960     format(' - file ',a,' re-written ok.')

! -- First a copy of the existing multiplier file is made.

        ii=0
        do
          ii=ii+1
          call num2char(ii,anum)
          anum=adjustl(anum)
          multfilekeep=trim(multfile)//'.k'//trim(anum)
          inquire(file=multfilekeep,exist=lexist)
          if(.not.lexist)exit
        end do
        write(6,*)
        write(6,1000) trim(multfile),trim(multfilekeep)
1000    format(' - copying file ',a,' to ',a,' for safekeeping...')
        multunit=nextunit()
        open(unit=multunit,file=multfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,1010) trim(multfile)
1010       format(' Cannot open multiplier file ',a,'.')
          go to 9890
        end if
        multunitkeep=nextunit()
        open(unit=multunitkeep,file=multfilekeep)
        do
          read(multunit,'(a)',err=9700,end=1030) cline
          write(multunitkeep,'(a)') trim(cline)
        end do
1030    continue
        close(unit=multunit)
        close(unit=multunitkeep)
        write(6,1040) trim(multfile)
1040    format(' - file ',a,' copied ok.')
        multkeepflag=1

! -- Now the copied multiplier file is re-opened and read a first time.

        nml=0
        write(6,*)
        write(6,1100) trim(multfile)
1100    format(' - modifying multiplier file ',a,'...')
        multunitkeep=nextunit()
        open(unit=multunitkeep,file=multfilekeep,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,1110) trim(multfilekeep)
1110      format(' Cannot re-open file ',a,'.')
          go to 9890
        end if
        do
          read(multunitkeep,'(a)',end=1135) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).ne.'#') exit
        end do
        call linesplit(ifail,1)
        nml=char2int(ifail,1)
        if(ifail.ne.0)then
          write(amessage,1120) trim(multfile)
1120      format(' Cannot read number of multiplier arrays variable from first line ',  &
          'of multiplier file ',a,'.')
          go to 9890
        end if
        do i=1,nml
1130      continue
          read(multunitkeep,'(a)',err=9700,end=9720) cline
          if(cline.eq.' ') go to 1130
          call linesplit(ifail,1)
          multname=cline(left_word(1):right_word(1))
          call casetrans(multname,'lo')
          iflag=0
          do j=1,distribnclu
            if(multname.eq.mltarr(j)) then
              ifound(j)=1
              iflag=1
            end if
          end do
          call linesplit(ifail,2)
          if(ifail.eq.0)then
            atemp=cline(left_word(2):right_word(2))
            call casetrans(atemp,'lo')
            if(atemp.eq.'function')then
              read(multunitkeep,'(a)',err=9750,end=9750) cline
              call casetrans(cline,'lo')
              if(iflag.eq.0)then
                do j=1,distribnclu
                  if(index(cline,trim(mltarr(j))).ne.0) then
                    discard(j)=0
                  end if
                end do
              end if
            else
              bname='"'//trim(multname)//'" mult. array'
              call u2drel(ifail,rarray,bname,nrow,ncol,0,multunitkeep,0,locat_mult,ioc,cline,external_file)
              if(ifail.ne.0) go to 9890
            end if
          else
            bname='"'//trim(multname)//'" mult. array'
            call u2drel(ifail,rarray,bname,nrow,ncol,0,multunitkeep,0,locat_mult,ioc,cline,external_file)
            if(ifail.ne.0) go to 9890
          end if
        end do

! -- If necessary, modifications are now made.

1135    continue
        icount=0
        do i=1,distribnclu
          if((discard(i).eq.1).and.(ifound(i).eq.1)) icount=icount+1
        end do
        newnml=nml-icount+1
        rewind(unit=multunitkeep)
        multunit=nextunit()
        open(unit=multunit,file=multfile)
        do
          read(multunitkeep,'(a)',err=9790) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).eq.'#') then
            write(multunit,'(a)') trim(cline)
          else
            exit
          end if
        end do
        write(multunit,1150) newnml
1150    format(i5)
        do i=1,nml
1155      continue
          read(multunitkeep,'(a)',err=9790) cline
          if(cline.eq.' ') go to 1155
          call linesplit(ifail,1)
          multname=cline(left_word(1):right_word(1))
          call casetrans(multname,'lo')
          iflag=0
          do j=1,distribnclu
            if(multname.eq.mltarr(j))then
              if(discard(j).eq.1) iflag=1
            end if
          end do
          ifunction=0
          call linesplit(ifail,2)
          if(ifail.eq.0)then
            atemp=cline(left_word(2):right_word(2))
            call casetrans(atemp,'lo')
            if(atemp.eq.'function') ifunction=1
          end if
          if(iflag.eq.0)then
            write(multunit,'(a)') trim(cline)
          end if
          if(ifunction.eq.1)then
            read(multunitkeep,'(a)',err=9790) cline
            if(iflag.eq.0)then
              write(multunit,'(a)') trim(cline)
            end if
          else
            bname='"'//trim(multname)//'" mult. array'
            call u2drel(ifail,rarray,bname,nrow,ncol,0,multunitkeep,0,locat_mult,ioc,cline,external_file)
            if(ifail.ne.0) go to 9890
            if(iflag.eq.0)then
              if(ioc.eq.1)then
                write(multunit,'(a)') trim(cline)
              else
                iconstant=1
                rtemp=rarray(1,1)
                do irow=1,nrow
                  do icol=1,ncol
                    if(rarray(icol,irow).ne.rtemp)then
                      iconstant=0
                      go to 1169
                    end if
                  end do
                end do
1169            continue
                if(iconstant.eq.1)then
                  write(multunit,1168) rtemp
1168              format('constant    ',1pg14.7)
                else
                  outfile=trim(multname)//'_mult.ref'
                  outunit=nextunit()
                  open(unit=outunit,file=outfile)
                  do irow=1,nrow
                    write(outunit,1170) (rarray(icol,irow),icol=1,ncol)
1170                format(8(1pg14.7,1x))
                  end do
                  close(unit=outunit)
                  write(6,1180) trim(outfile)
1180              format(' - file ',a,' written ok.')
                  write(multunit,1190) trim(outfile)
1190              format('  OPEN/CLOSE  ''',a,'''     1.0 ''(FREE)'' -1')
                end if
              end if
            else
              if(ioc.eq.1)then
                if(index(external_file,'_mult').ne.0)then
                  delunit=nextunit()
                  open(unit=delunit,file=external_file,status='old',iostat=ierr)
                  if(ierr.eq.0)then
                    close(unit=delunit,status='delete',iostat=ierr)
                  end if
                end if
              end if
            end if
          end if
        end do

! -- Now we add the new multiplier file.

        newmultfile=trim(newmltarr)//'.ref'
        write(multunit,1200) trim(newmltarr)
1200    format(a)
        outfile=newmultfile
        write(multunit,1190) trim(outfile)

        close(unit=multunit)
        close(unit=multunitkeep)
        write(6,1205) trim(multfile)
1205    format(' - new multiplier file ',a,' written ok.')

! -- The first part of the factor file is read.

1300    continue

        write(6,*)
        write(6,1310) trim(facfile)
1310    format(' - processing interpolation factor file ',a,'...')
        if(facformat.eq.'f')then
          read(facunit,*,err=9770,end=9770) ppfile
          read(facunit,*,err=9770,end=9770) intfile
          read(facunit,*,err=9770,end=9770) mcol,mrow
        else
          read(facunit,err=9770,end=9770) ppfile
          read(facunit,err=9770,end=9770) intfile
          read(facunit,err=9770,end=9770) mcol,mrow
        end if
        call remquote(ppfile)
        call remquote(intfile)       !    (make sure this includes adjustl)
        if((mcol.ne.ncol).or.(mrow.ne.nrow))then
          write(amessage,1320) trim(facfile),trim(disfile)
1320      format(' Dimensions of finite difference grid as read from file ',a,' do not ',  &
          'agree with those as read from MODFLOW discretisation file ',a,'.')
          go to 9890
        end if
        if(facformat.eq.'f')then
          read(facunit,*,err=9770,end=9770) npp
        else
          read(facunit,err=9770,end=9770) npp
        end if
        allocate(ppname(npp),ppname1(npp),stat=ierr)
        if(ierr.ne.0) go to 9200
        do i=1,npp
          if(facformat.eq.'f')then
            read(facunit,*,err=9770,end=9770) ppname(i)
          else
            read(facunit,err=9770,end=9770) tempname
            ppname(i)=tempname
          end if
          call casetrans(ppname(i),'lo')
        end do

! -- We now work out how to form parameter names.

        maxlen=0
        do ipp=1,npp
          ppname(ipp)=adjustl(ppname(ipp))
          ppname1(ipp)=ppname(ipp)
          nb=len_trim(ppname(ipp))
          if(nb.gt.maxlen)maxlen=nb
        end do
        nbp=len_trim(distribpar)
        nbp=12-nb-1
        do ipp=1,npp
          ppname(ipp)=trim(distribpar(1:nbp))//'_'//trim(ppname(ipp))
        end do

! -- We now read the rest of the interpolation factor file to work out a value for NDIM.
! -- We also start writing the distributed-to-PEST parameter file.

        ndim=0
        maxent=0
        itransold=-999
        if(facformat.eq.'f')then
          do
            read(facunit,*,err=9770,end=1400) icellno,itrans,itemp2
            ndim=ndim+itemp2
            if(itemp2.gt.maxent)maxent=itemp2
            if(itransold.eq.-999)then
              itransold=itrans
            else
              if(itransold.ne.itrans) go to 9850
            end if
          end do
        else
          do
            read(facunit,err=9770,end=1400) icellno,itrans,itemp2
            ndim=ndim+itemp2
            if(itemp2.gt.maxent)maxent=itemp2
            if(itransold.eq.-999)then
              itransold=itrans
            else
              if(itransold.ne.itrans) go to 9850
            end if
          end do
        end if
1400    continue

        write(defnunit,1346) (ndim+1)*distribnclu,itrans
1346    format(2i10)

        allocate(ilaycount(nlay),stat=ierr)
        if(ierr.ne.0) go to 9200

        ilaycount=0                        ! an array
        do ilay=1,nlay
          do irow=1,nrow
            do icol=1,ncol
              if(iarray(icol,irow,ilay).ne.0) then
                ilaycount(ilay)=1
                go to 1350
              end if
            end do
          end do
1350      continue
        end do
        mlay=0
        do ilay=1,nlay
          if(ilaycount(ilay).ne.0) mlay=mlay+1
        end do
        write(defnunit,1360) mlay
        do ilay=1,nlay
          if(ilaycount(ilay).gt.0) write(defnunit,1360) ilay
        end do
1360    format(i5)
        write(defnunit,1360) npp
        do ipp=1,npp
          write(defnunit,1370) trim(ppname(ipp))
1370      format(a)
        end do
        write(defnunit,1360) maxent

        allocate(ipt(maxent),wt(maxent),stat=ierr)
        if(ierr.ne.0) go to 9200

        do ilay=1,nlay
          if(ilaycount(ilay).gt.0)then
            rewind(unit=facunit)
            if(facformat.eq.'f')then
              read(facunit,*,err=9770,end=9770) atempf
              read(facunit,*,err=9770,end=9770) atempf
              read(facunit,*,err=9770,end=9770) mcol,mrow
              read(facunit,*,err=9770,end=9770) itemp1
              do i=1,npp
                read(facunit,*,err=9770,end=9770) tempname
              end do
            else
              read(facunit,err=9770,end=9770) atempf
              read(facunit,err=9770,end=9770) atempf
              read(facunit,err=9770,end=9770) mcol,mrow
              read(facunit,err=9770,end=9770) itemp1
              do i=1,npp
                read(facunit,err=9770,end=9770) tempname
              end do
            end if
            do
              if(facformat.eq.'f')then         !          THERE IS AN ISSUE!!!! LOG INTERPOLATION
                read(facunit,*,err=9700,end=1391) icellno,itrans,na,rtemp,((ipt(i),wt(i)),i=1,na)
              else
                read(facunit,err=9700,end=1391) icellno,itrans,na,rtemp,((ipt(i),wt(i)),i=1,na)
              end if
              irow=(icellno-1)/ncol+1
              icol=icellno-((irow-1)*ncol)
              if((icol.gt.ncol).or.(icol.lt.1).or.(irow.gt.nrow).or.(irow.lt.1))then
                call num2char(icellno,anum)
                write(amessage,1389) trim(anum),trim(facfile)
1389            format(' A cell number of "',a,'" is provided in file ',a,'. This number ',   &
                'results in grid (row,column) numbers which are out of range.')
                go to 9890
              end if
              if(iarray(icol,irow,ilay).ne.0)then
                write(defnunit,1390) icol,irow,ilay,na,((ipt(i),wt(i)),i=1,na)
1390            format(4i6,100(i6,1x,1pg14.7))
              end if
            end do
1391        continue
          end if
        end do
        close(unit=facunit)
        close(unit=defnunit)
        write(6,1405) trim(facfile)
1405    format(' - file ',a,' read ok.')
        write(6,1420) trim(defnfile)
1420    format(' - file ',a,' written ok.')

! -- The pilot points template file is now written.

        templatefile=trim(distribpar)//'.tpl'
        modinfile=trim(distribpar)//'.pts'
        write(6,*)
        write(6,1430) trim(templatefile)
1430    format(' - writing pilot points template file ',a,'...')
        ppunit=nextunit()
        open(unit=ppunit,file=ppfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,1440) trim(ppfile),trim(facfile)
1440      format(' Cannot open pilot points file ',a,' cited in interpolation factor file ',a,'.')
          go to 9890
        end if
        templateunit=nextunit()
        open(unit=templateunit,file=templatefile)
        write(templateunit,1450)
1450    format('ptf $')
        do ipp=1,npp
          read(ppunit,*,err=9800,end=9820) aname,easting,northing,izone
          call casetrans(aname,'lo')
          if(aname.ne.ppname1(ipp))then
            write(amessage,1460) trim(ppfile),trim(facfile)
1460        format(' Pilot point names and ordering in file ',a,' do not correspond to ',  &
            'those in interpolation factor file ',a,'.')
            go to 9890
          end if
          write(templateunit,1470) trim(aname),easting,northing,izone,ppname(ipp)
1470      format(1x,a,t15,f12.3,3x,f12.3,2x,i5,3x,'$',a,'  $')
        end do
        close(unit=templateunit)
        close(unit=ppunit)
        write(6,1405) trim(ppfile)
        write(6,1420) trim(templatefile)


! -- Now the PEST file fragment must be written.

        write(6,*)
        write(6,1500) trim(pestfile)
1500    format(' - writing PEST building block file ',a,'....')
        write(pestunit,1530)
1530    format('* parameter groups')
        write(pestunit,1550) trim(distribpar)
1550    format(a,' relative 0.01 0.0 switch 2.0 parabolic')
        write(pestunit,1560)
1560    format('* parameter data')
        do ipp=1,npp
          write(pestunit,1570) trim(ppname(ipp)),initval,lbound,ubound,trim(distribpar)
1570      format(1x,a,t15,'log factor',3(2x,1pg14.7),2x,a,2x,' 1.0   0.0  0')
        end do
        write(pestunit,1575)
1575    format('* model input/output')
        write(pestunit,1576) trim(templatefile),trim(modinfile)
1576    format(a,3x,a)
        close(unit=pestunit)
        write(6,1580) trim(pestfile)
1580    format(' - file ',a,' written ok.')

! -- The FAC2REAL input re-direction file is now written.

        facinfile='fac2real_'//trim(distribpar)//'.in'
        write(6,*)
        write(6,1600) trim(facinfile)
1600    format(' - writing FAC2REAL keyboard input file ',a,'...')
        facunit=nextunit()
        open(unit=facunit,file=facinfile)
        call addquote(facfile,afile)
        write(facunit,'(a)') trim(afile)
        write(facunit,'(a)') trim(facformat)
        call addquote(modinfile,afile)
        write(facunit,'(a)') trim(afile)
        write(facunit,1620)
1620    format('s')
        write(facunit,1630)
1630    format('1.0e-20')
        write(facunit,1620)
        write(facunit,1640)
1640    format('1.0e20')
        call addquote(newmultfile,afile)
        write(facunit,'(a)') trim(afile)
        write(facunit,1650)
1650    format('1.0e35')
        close(unit=facunit)
        write(6,1660) trim(facinfile)
1660    format(' - file ',a,' written ok.')


        write(amessage,1602) trim(newmultfile)
1602    format(' Before the model can be run, FAC2REAL must be run to populate ', &
        'the external multiplier real array file ',a,'. This can be done using the command: ')
        call write_message(leadspace='yes')
        write(6,1622) trim(facinfile)
1622    format(/,'     fac2real < ',a)
        write(amessage,1634) trim(modinfile),trim(templatefile)
1634    format(' after the model input file ',a,' has been populated on the basis of ', &
        'the template file ',a,' and appropriate initial parameter values.')
        call write_message(leadspace='yes')
        write(amessage,1632)
1632    format(' The same command must be added to the model batch file.')
        call write_message(leadspace='yes')
        go to 9900

9000    write(amessage,9010) trim(disfile)
9010    format(' Error reading MODFLOW discretisation file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(disfile)
9060    format(' Premature end encountered to MODFLOW discretisation file ',a,'.')
        go to 9890
9070    call num2char(iline,aline)
        write(amessage,9080) trim(aline),trim(disfile)
9080    format(' Error in data at line ',a,' of MODFLOW discretisation file ',a,'.')
        go to 9890
9100    write(amessage,9110) trim(basfile)
9110    format(' Error reading MODFLOW BASIC package input file ',a,'.')
        go to 9890
9150    write(amessage,9160) trim(basfile)
9160    format(' Premature end encountered to MODFLOW BASIC package input file ',a,'.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890
9300    write(amessage,9310) trim(lpffile)
9310    format(' Error reading MODFLOW LPF package input file ',a,'.')
        go to 9890
9350    write(amessage,9360) trim(lpffile)
9360    format(' Premature end encountered to MODFLOW LPF package input file ',a,'.')
        go to 9890
9370    call num2char(iline,aline)
        write(amessage,9380) trim(aline),trim(lpffile)
9380    format(' Error in data at line ',a,' of MODFLOW LPF package input file ',a,'.')
        go to 9890
9390    write(amessage,9395) trim(arraytype),trim(lpffile)
9395    format(' Error reading ',a,' array from MODFLOW LPF package input file ',a,'.')
        go to 9890
9400    write(amessage,9410) trim(lpffile)
9410    format(' Error reading data pertaining to MODFLOW parameters from LPF package ',  &
        'input file ',a,'.')
        go to 9890
9420    write(amessage,9425) trim(lpffile)
9425    format(' Premature end encountered to MODFLOW LPF package input file ',a,' while reading ', &
        'MODFLOW parameter data.')
        go to 9890
9430    write(amessage,9440) trim(distribpar),trim(lpffile)
9440    format(' Error reading data for parameter "',a,'" from MODFLOW LPF package ',  &
        'input file ',a,'.')
        go to 9890
9450    write(amessage,9460) trim(lpffile),trim(distribpar)
9460    format(' Premature end encountered to MODFLOW LPF package input file ',a,   &
        ' when reading data for parameter "',a,'".')
        go to 9890
9470    write(amessage,9440) trim(distribpar),trim(lpffilekeep)
        go to 9890
9490    write(amessage,9460) trim(lpffilekeep),trim(distribpar)
9500    write(amessage,9510) trim(zonefile)
9510    format(' Error reading MODFLOW zone file ',a,'.')
        go to 9890
9520    write(amessage,9530) trim(zonefile)
9530    format(' Premature end encountered to MODFLOW zone file ',a,'.')
        go to 9890
9600    write(amessage,9610) trim(lpffilekeep)
9610    format(' Error re-reading file ',a,'.')
        go to 9890
9700    write(amessage,9710) trim(multfile)
9710    format(' Error reading multiplier array file ',a,'.')
        go to 9890
9720    write(amessage,9730) trim(multfile)
9730    format(' Premature end encountered to multiplier array file ',a,'.')
        go to 9890
9750    write(amessage,9760) trim(multname),trim(multfile)
9760    format(' Error reading function string for multiplier array "',a,'" from multiplier file ',a,'.')
        go to 9890
9770    write(amessage,9780) trim(facfile)
9780    format(' Error or premature end encountered in reading interpolation factor file ',a,'.')
        go to 9890
9790    write(amessage,9795) trim(multfilekeep)
9795    format(' Error re-reading file ',a,'.')
        go to 9890
9800    write(amessage,9810) trim(ppfile),trim(facfile)
9810    format(' Error encountered in reading pilot points file ',a,' cited in interpolation ',  &
        'factor file ',a,'.')
        go to 9890
9820    write(amessage,9830) trim(ppfile),trim(facfile)
9830    format(' Premature end encountered to pilot points file ',a,' cited in interpolation ',  &
        'factor file ',a,'.')
        go to 9890
9850    write(amessage,9860) trim(facfile)
9860    format(' According to interpolation factor file ',a,' some geostatistical structures ',  &
        'on which interpolation is based pertain to parameter logs while others pertain to native ', &
        'parameters. All interpolation must pertain to either log or native parameters, but not both.')
        go to 9890

9890    call write_message(leadspace='yes')
        messageflag=1
9900    call close_files

        if(messageflag.ne.0)then
          if(lpfkeepflag.eq.1)then
            write(6,9902) trim(lpffilekeep),trim(lpffile)
9902        format(' - copying file ',a,' back to file ',a,'...')
            call copyfile(ifail,lpffilekeep,lpffile)
          end if
          if(multkeepflag.eq.1)then
            write(6,9902) trim(multfilekeep),trim(multfile)
            call copyfile(ifail,multfilekeep,multfile)
          end if
        end if

        if(allocated(laycbd))deallocate(laycbd,stat=ierr)
        if(allocated(iarray))deallocate(iarray,stat=ierr)
        if(allocated(rarray))deallocate(rarray,stat=ierr)
        if(allocated(itemparray))deallocate(itemparray,stat=ierr)
        if(allocated(ipt))deallocate(ipt,stat=ierr)
        if(allocated(ivec))deallocate(ivec,stat=ierr)
        if(allocated(rvec))deallocate(rvec,stat=ierr)
        if(allocated(discard))deallocate(discard,stat=ierr)
        if(allocated(ifound))deallocate(ifound,stat=ierr)
        if(allocated(ilaycount))deallocate(ilaycount,stat=ierr)
        if(allocated(wt))deallocate(wt,stat=ierr)
        if(allocated(ppname))deallocate(ppname,stat=ierr)

end program ppmdef



subroutine remquote(astring)

        implicit none

        integer       :: n
        character*(*) :: astring

        do
          n=index(astring,'"')
          if(n.eq.0) exit
          astring(n:n)=' '
        end do
        continue
        do
          n=index(astring,'''')
          if(n.eq.0) exit
          astring(n:n)=' '
        end do

        astring=adjustl(astring)
        return

end subroutine remquote




      SUBROUTINE U2DREL(ifail,A,ANAME,II,JJ,K,IN,IOUT,locat_dis,ioc,cline,external_file)
!     ******************************************************************
!     ROUTINE TO INPUT 2-D REAL DATA MATRICES
!       A IS ARRAY TO INPUT
!       ANAME IS 24 CHARACTER DESCRIPTION OF A
!       II IS NO. OF ROWS
!       JJ IS NO. OF COLS
!       K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --)
!              IF K=0, NO LAYER IS PRINTED
!              IF K<0, CROSS SECTION IS PRINTED)
!       IN IS INPUT UNIT
!       IOUT IS OUTPUT UNIT
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------

      use defn, only : amessage
      implicit integer (i-n), real (a-h, o-z)

      integer ifail,locat_dis
      character*(*) cline
      CHARACTER*24 ANAME
      DIMENSION A(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*16 TEXT
      CHARACTER*200 FNAME
      character*(*) external_file
      DATA NUNOPN/99/
      INCLUDE 'openspec.inc'
!     ------------------------------------------------------------------
      ifail=0
      ioc=0
!
!1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)') CNTRL
!
!2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         write(amessage,8888)
8888     format(' "EXTERNAL" not allowed in array header when using PPMDEF.')
         ifail=1
         return
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
!!!!         WRITE(IOUT,15) LOCAT,FNAME
!!!!   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         ICLOSE=1
         ioc=1
         cline=cntrl
         external_file=fname
         return
      ELSE
!
!2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=500) LOCAT,CNSTNT,FMTIN,IPRN
    1    FORMAT(I10,F10.0,A20,I10)
         if(locat.eq.locat_dis)then
           locat=in
         else
           write(amessage,8889)
8889       format(' MODFLOW arrays must be read from same file as package input file or ',   &
           'from another file using OPEN/CLOSE functionality when using PPMDEF.')
           ifail=1
           return
         end if
      END IF
!
!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN)
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            IF(ICLOSE.NE.0) THEN
               IF(FMTIN.EQ.'(BINARY)') THEN
                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM=FORM,ACCESS=ACCESS,    &
                       ACTION=ACTION(1))
               ELSE
                  OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1))
               END IF
            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)') LOCAT=-LOCAT
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
         END IF
      END IF
!
!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.EQ.0) THEN
!
!4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
        DO 80 I=1,II
        DO 80 J=1,JJ
   80   A(J,I)=CNSTNT
!!!!        IF(K.GT.0) WRITE(IOUT,2) ANAME,CNSTNT,K
!!!!    2   FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR LAYER',I4)
!!!!        IF(K.LE.0) WRITE(IOUT,3) ANAME,CNSTNT
!!!!    3   FORMAT(1X,/1X,A,' =',1P,G14.6)
        RETURN
      ELSE IF(LOCAT.GT.0) THEN
!
!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
!!!!        IF(K.GT.0) THEN
!!!!           WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
!!!!   94      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        ELSE IF(K.EQ.0) THEN
!!!!           WRITE(IOUT,95) ANAME,LOCAT,FMTIN
!!!!   95      FORMAT(1X,///11X,A,/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        ELSE
!!!!           WRITE(IOUT,96) ANAME,LOCAT,FMTIN
!!!!   96      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        END IF
        DO 100 I=1,II
        IF(FMTIN.EQ.'(FREE)') THEN
           READ(LOCAT,*) (A(J,I),J=1,JJ)
        ELSE
           READ(LOCAT,FMTIN) (A(J,I),J=1,JJ)
        END IF
  100   CONTINUE
      ELSE
!
!4C-----LOCAT<0; READ UNFORMATTED ARRAY VALUES.
        LOCAT=-LOCAT
!!!!        IF(K.GT.0) THEN
!!!!           WRITE(IOUT,201) ANAME,K,LOCAT
!!!!  201      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        ELSE IF(K.EQ.0) THEN
!!!!           WRITE(IOUT,202) ANAME,LOCAT
!!!!  202      FORMAT(1X,///1X,A,/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        ELSE
!!!!           WRITE(IOUT,203) ANAME,LOCAT
!!!!  203      FORMAT(1X,///1X,A,' FOR CROSS SECTION',/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        END IF
        READ(LOCAT) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,NROW,ILAY
        READ(LOCAT) A
      END IF
!
!5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
      ZERO=0.
      IF(CNSTNT.EQ.ZERO) GO TO 320
      DO 310 I=1,II
      DO 310 J=1,JJ
      A(J,I)=A(J,I)*CNSTNT
  310 CONTINUE
!
!6------IF PRINT CODE (IPRN) >0 OR =0 THEN PRINT ARRAY VALUES.
!!!!  320 IF(IPRN.GE.0) CALL ULAPRW(A,ANAME,0,0,JJ,II,0,IPRN,IOUT)
320   continue
!
!7------RETURN
      RETURN
!
!8------CONTROL RECORD ERROR.
500   continue
      IF(K.GT.0) THEN
         WRITE(amessage,501) trim(ANAME),K
  501    FORMAT(' ERROR READING ARRAY CONTROL RECORD FOR ',A,' FOR LAYER',I4,':')
      ELSE
         WRITE(amessage,502) trim(ANAME)
  502    FORMAT(' ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
      END IF
      ifail=1
      return
!!!!      WRITE(IOUT,'(1X,A)') CNTRL
!!!!      CALL USTOP(' ')
      END


      SUBROUTINE U2DINT(ifail,IA,ANAME,II,JJ,K,IN,IOUT,locat_bas)
!     ******************************************************************
!     ROUTINE TO INPUT 2-D INTEGER DATA MATRICES
!       IA IS ARRAY TO INPUT
!       ANAME IS 24 CHARACTER DESCRIPTION OF IA
!       II IS NO. OF ROWS
!       JJ IS NO. OF COLS
!       K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --
!              IF K=0, NO LAYER IS PRINTED
!              IF K<0, CROSS SECTION IS PRINTED)
!       IN IS INPUT UNIT
!       IOUT IS OUTPUT UNIT
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------

      use defn, only : amessage
      implicit integer (i-n), real (a-h, o-z)

      integer ifail,locat_bas
      CHARACTER*24 ANAME
      DIMENSION IA(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*200 FNAME
      DATA NUNOPN/99/
      INCLUDE 'openspec.inc'
!     ------------------------------------------------------------------

      ifail=0
!
!1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)') CNTRL
!
!2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         write(amessage,8888)
8888     format(' "EXTERNAL" not allowed in array header when using PPMDEF.')
         ifail=1
         return
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
!!!!         WRITE(IOUT,15) LOCAT,FNAME
!!!!   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         ICLOSE=1
      ELSE
!
!2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=600) LOCAT,ICONST,FMTIN,IPRN
    1    FORMAT(I10,I10,A20,I10)
         if(locat.eq.locat_bas)then
           locat=in
         else
           write(amessage,8889)
8889       format(' MODFLOW arrays must be read from same file as package input file or ',   &
           'from another file using OPEN/CLOSE functionality when using PPMDEF.')
           ifail=1
           return
         end if
      END IF
!
!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,ICONST,R,IOUT,IN)
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            IF(ICLOSE.NE.0) THEN
               IF(FMTIN.EQ.'(BINARY)') THEN
                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM=FORM,ACCESS=ACCESS,    &
                       ACTION=ACTION(1))
               ELSE
                  OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1))
               END IF
            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)') LOCAT=-LOCAT
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
         END IF
      END IF
!
!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.EQ.0) THEN
!
!4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO ICONST. RETURN.
        DO 80 I=1,II
        DO 80 J=1,JJ
   80   IA(J,I)=ICONST
!!!!        IF(K.GT.0) WRITE(IOUT,82) ANAME,ICONST,K
!!!!   82   FORMAT(1X,/1X,A,' =',I15,' FOR LAYER',I4)
!!!!        IF(K.LE.0) WRITE(IOUT,83) ANAME,ICONST
!!!!   83   FORMAT(1X,/1X,A,' =',I15)
        RETURN
      ELSE IF(LOCAT.GT.0) THEN
!
!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
!!!!        IF(K.GT.0) THEN
!!!!           WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
!!!!   94      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        ELSE IF(K.EQ.0) THEN
!!!!           WRITE(IOUT,95) ANAME,LOCAT,FMTIN
!!!!   95      FORMAT(1X,///11X,A,/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        ELSE
!!!!           WRITE(IOUT,96) ANAME,LOCAT,FMTIN
!!!!   96      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!!!!     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
!!!!        END IF
        DO 100 I=1,II
        IF(FMTIN.EQ.'(FREE)') THEN
           READ(LOCAT,*) (IA(J,I),J=1,JJ)
        ELSE
           READ(LOCAT,FMTIN) (IA(J,I),J=1,JJ)
        END IF
  100   CONTINUE
      ELSE
!
!4C-----LOCAT<0; READ UNFORMATTED RECORD CONTAINING ARRAY VALUES.
        LOCAT=-LOCAT
!!!!        IF(K.GT.0) THEN
!!!!           WRITE(IOUT,201) ANAME,K,LOCAT
!!!!  201      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        ELSE IF(K.EQ.0) THEN
!!!!           WRITE(IOUT,202) ANAME,LOCAT
!!!!  202      FORMAT(1X,///11X,A,/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        ELSE
!!!!           WRITE(IOUT,203) ANAME,LOCAT
!!!!  203      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!!!!     1      1X,'READING BINARY ON UNIT ',I4)
!!!!        END IF
        READ(LOCAT)
        READ(LOCAT) IA
      END IF
!
!5------IF ICONST NOT ZERO THEN MULTIPLY ARRAY VALUES BY ICONST.
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
      IF(ICONST.EQ.0) GO TO 320
      DO 310 I=1,II
      DO 310 J=1,JJ
      IA(J,I)=IA(J,I)*ICONST
  310 CONTINUE
!
!6------IF PRINT CODE (IPRN) <0 THEN RETURN.
320   continue
!
!9------RETURN
      RETURN
!
!10-----CONTROL RECORD ERROR.
  600 IF(K.GT.0) THEN
         WRITE(amessage,601) trim(ANAME),K
  601    FORMAT(' ERROR READING ARRAY CONTROL RECORD FOR ',A,' FOR LAYER',I4,':')
      ELSE
         WRITE(amessage,602) trim(ANAME)
  602    FORMAT(' ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
      END IF
      ifail=1
      return
      END






      SUBROUTINE URWORD(LINE,ICOL,ISTART,ISTOP,NCODE,N,R,IOUT,IN)
!     ******************************************************************
!     ROUTINE TO EXTRACT A WORD FROM A LINE OF TEXT, AND OPTIONALLY
!     CONVERT THE WORD TO A NUMBER.
!        ISTART AND ISTOP WILL BE RETURNED WITH THE STARTING AND
!          ENDING CHARACTER POSITIONS OF THE WORD.
!        THE LAST CHARACTER IN THE LINE IS SET TO BLANK SO THAT IF ANY
!          PROBLEMS OCCUR WITH FINDING A WORD, ISTART AND ISTOP WILL
!          POINT TO THIS BLANK CHARACTER.  THUS, A WORD WILL ALWAYS BE
!          RETURNED UNLESS THERE IS A NUMERIC CONVERSION ERROR.  BE SURE
!          THAT THE LAST CHARACTER IN LINE IS NOT AN IMPORTANT CHARACTER
!          BECAUSE IT WILL ALWAYS BE SET TO BLANK.
!        A WORD STARTS WITH THE FIRST CHARACTER THAT IS NOT A SPACE OR
!          COMMA, AND ENDS WHEN A SUBSEQUENT CHARACTER THAT IS A SPACE
!          OR COMMA.  NOTE THAT THESE PARSING RULES DO NOT TREAT TWO
!          COMMAS SEPARATED BY ONE OR MORE SPACES AS A NULL WORD.
!        FOR A WORD THAT BEGINS WITH "'", THE WORD STARTS WITH THE
!          CHARACTER AFTER THE QUOTE AND ENDS WITH THE CHARACTER
!          PRECEDING A SUBSEQUENT QUOTE.  THUS, A QUOTED WORD CAN
!          INCLUDE SPACES AND COMMAS.  THE QUOTED WORD CANNOT CONTAIN
!          A QUOTE CHARACTER.
!        IF NCODE IS 1, THE WORD IS CONVERTED TO UPPER CASE.
!        IF NCODE IS 2, THE WORD IS CONVERTED TO AN INTEGER.
!        IF NCODE IS 3, THE WORD IS CONVERTED TO A REAL NUMBER.
!        NUMBER CONVERSION ERROR IS WRITTEN TO UNIT IOUT IF IOUT IS
!          POSITIVE; ERROR IS WRITTEN TO DEFAULT OUTPUT IF IOUT IS 0;
!          NO ERROR MESSAGE IS WRITTEN IF IOUT IS NEGATIVE.
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
      implicit integer (i-n), real (a-h, o-z)
      CHARACTER*(*) LINE
      CHARACTER*20 STRING
      CHARACTER*30 RW
      CHARACTER*1 TAB
!     ------------------------------------------------------------------
      TAB=CHAR(9)
!
!1------Set last char in LINE to blank and set ISTART and ISTOP to point
!1------to this blank as a default situation when no word is found.  If
!1------starting location in LINE is out of bounds, do not look for a
!1------word.
      LINLEN=LEN(LINE)
      LINE(LINLEN:LINLEN)=' '
      ISTART=LINLEN
      ISTOP=LINLEN
      LINLEN=LINLEN-1
      IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
!
!2------Find start of word, which is indicated by first character that
!2------is not a blank, a comma, or a tab.
      DO 10 I=ICOL,LINLEN
      IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.','     &
          .AND. LINE(I:I).NE.TAB) GO TO 20
10    CONTINUE
      ICOL=LINLEN+1
      GO TO 100
!
!3------Found start of word.  Look for end.
!3A-----When word is quoted, only a quote can terminate it.
20    IF(LINE(I:I).EQ.'''') THEN
         I=I+1
         IF(I.LE.LINLEN) THEN
            DO 25 J=I,LINLEN
            IF(LINE(J:J).EQ.'''') GO TO 40
25          CONTINUE
         END IF
!
!3B-----When word is not quoted, space, comma, or tab will terminate.
      ELSE
         DO 30 J=I,LINLEN
         IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.','      &
          .OR. LINE(J:J).EQ.TAB) GO TO 40
30       CONTINUE
      END IF
!
!3C-----End of line without finding end of word; set end of word to
!3C-----end of line.
      J=LINLEN+1
!
!4------Found end of word; set J to point to last character in WORD and
!-------set ICOL to point to location for scanning for another word.
40    ICOL=J+1
      J=J-1
      IF(J.LT.I) GO TO 100
      ISTART=I
      ISTOP=J
!
!5------Convert word to upper case and RETURN if NCODE is 1.
      IF(NCODE.EQ.1) THEN
         IDIFF=ICHAR('a')-ICHAR('A')
         DO 50 K=ISTART,ISTOP
            IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z')                     &
                   LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50       CONTINUE
         RETURN
      END IF
!
!6------Convert word to a number if requested.
100   IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
         RW=' '
         L=30-ISTOP+ISTART
         IF(L.LT.1) GO TO 200
         RW(L:30)=LINE(ISTART:ISTOP)
         IF(NCODE.EQ.2) READ(RW,'(I30)',ERR=200) N
         IF(NCODE.EQ.3) READ(RW,'(F30.0)',ERR=200) R
      END IF
      RETURN
!
!7------Number conversion error.
200   IF(NCODE.EQ.3) THEN
         STRING= 'A REAL NUMBER'
         L=13
      ELSE
         STRING= 'AN INTEGER'
         L=10
      END IF
!
!7A-----If output unit is negative, set last character of string to 'E'.
      IF(IOUT.LT.0) THEN
         N=0
         R=0.
         LINE(LINLEN+1:LINLEN+1)='E'
         RETURN
!
!7B-----If output unit is positive; write a message to output unit.
      ELSE IF(IOUT.GT.0) THEN
         IF(IN.GT.0) THEN
            WRITE(IOUT,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
         ELSE
            WRITE(IOUT,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
         END IF
201      FORMAT(1X,/1X,'FILE UNIT ',I4,' : ERROR CONVERTING "',A,      &
             '" TO ',A,' IN LINE:',/1X,A)
202      FORMAT(1X,/1X,'KEYBOARD INPUT : ERROR CONVERTING "',A,        &
             '" TO ',A,' IN LINE:',/1X,A)
!
!7C-----If output unit is 0; write a message to default output.
      ELSE
         IF(IN.GT.0) THEN
            WRITE(*,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
         ELSE
            WRITE(*,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
         END IF
      END IF
!
!7D-----STOP after writing message.
!      CALL USTOP(' ')
      END


      subroutine copyfile(ifail,file1,file2)

        use defn
        use inter
        implicit none

        integer       :: ifail,unit1,unit2,ierr
        character*(*) :: file1,file2

        ifail=0
        unit1=nextunit()
        open(unit=unit1,file=file1,status='old',iostat=ierr)
        if(ierr.ne.0) go to 9890
        unit2=nextunit()
        open(unit=unit2,file=file2)
        do
          read(unit1,'(a)',end=100) cline
          write(unit2,'(a)') trim(cline)
        end do
100     close(unit=unit1)
        close(unit=unit2)

        return

9890    ifail=1

        return
      end






! -- Do a test where significant parts of a layer are removed through absense of a value in an
!       integer array cited in a cluster. This tests that the iz() arrays have been properly read.
!       That is - check the stuff above line 795.
! -- Check that things work ok when the multiplier file is empty.

