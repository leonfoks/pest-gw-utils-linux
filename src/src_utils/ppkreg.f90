!     Last change:  JD    5 May 2003    2:53 pm


program ppkreg

! -- Program PPKREG adds a prior information and regularisation section to a PEST
!    control file whose parameters are based on pilot points. The regularisation
!    weights are assigned based on information contained within a PPK2FAC-generated
!    regularisation information file.

	use defn
	use inter

	implicit none



        integer                :: ifail,iline,npar,nobs,npargp,nprior,nobsgp, &
                                  ntplfle,ninsfle,ipar,numpoints,itype,numpartype, &
                                  j,jpar,ilastflag,jpoint,icount,ipoint,iflag, &
                                  npriornew,nobsgpnew,iprior,ierr,i,itemp,nbl,ii,  &
                                  lenchar
        integer                :: pestunit,regunit,pestoutunit
        integer, allocatable   :: itrans(:),ipartype(:),iparpoint(:),ipointnolink(:), &
                                  ius(:)
        real                   :: phimlim,wfinit,weight,wfmin,wfmax,rtemp
        real, allocatable      :: vario(:,:)
        real, allocatable      :: weightmul(:)
        character (len=1)      :: ayn,aus
        character (len=15)     :: aline,pestmode,anum,aobsgp,atemp,atemp1, &
                                  aprior,aweight,atrans,atempr
        character (len=120)    :: aprompt
        character (len=120)    :: pestfile,afile,regfile,pestoutfile
        character (len=300)    :: dline
        character (len=5), allocatable  :: prefix(:)
        character (len=12), allocatable :: point_id(:)
        character (len=12), allocatable :: apar(:),aparnopoint(:),apointnolink(:), &
                                           reggroup(:),piroot(:),oldgroup(:)

! -- The following variables must be adjusted together.

        integer, parameter         :: LENELINE=8000
        character (len=LENELINE)   :: eline
                                                                                      

	write(amessage,5)
5	format(' Program PPKREG adds a prior information and regularisation section ',&
        'to a PEST control file whose parameters are based on pilot points. ', &
        'Regularisation weights are assigned based on information contained ', &
        'within a PPK2FAC-generated regularisation information file.')
	call write_message(leadspace='yes',endspace='yes')


! -- The PEST control file is partially read to obtain the names and status
!    of the parameters featured in it.

25	aprompt=' Enter name of PEST control file: '
	call open_input_file(ifail,aprompt,pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

        afile=pestfile

        iline=1
        read(pestunit,700,err=9000,end=9050) cline
        cline=adjustl(cline)
        call casetrans(cline,'lo')
        if(cline(1:4).ne.'pcf ')then
          write(amessage,30) trim(pestfile)
30        format(' File ',a,' is not a PEST control file - try again.')
          call write_message(leadspace='yes',endspace='yes')
          close(unit=pestunit)
          go to 25
        end if
        iline=iline+1
        read(pestunit,700,err=9000,end=9050) cline
        iline=iline+1
        read(pestunit,700,err=9000,end=9050) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        pestmode=cline(left_word(2):right_word(2))
        call casetrans(pestmode,'lo')
!        if(pestmode(1:7).ne.'estimat')then
!          write(amessage,50) trim(pestfile)
!50        format(' PESTMODE in file ',a,' must be "estimation".')
!          go to 9890
!        end if
        iline=iline+1
        read(pestunit,*,err=9000,end=9050) npar,nobs,npargp,nprior,nobsgp
        if(npar.eq.1)then
          write(amessage,60) trim(pestfile)
60        format(' There is only one parameter cited in PEST control file ',a, &
          '. More parameters than this are required if regularisation is to be ', &
          'applied.')
          go to 9890
        end if
        allocate(apar(npar),itrans(npar),ipartype(npar),iparpoint(npar), &
        aparnopoint(npar),oldgroup(nobsgp),stat=ierr)
        if(ierr.ne.0) go to 9150
        iline=iline+1
        read(pestunit,*,err=9000,end=9050) ntplfle,ninsfle
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9070) cline
          call casetrans(cline,'lo')
          if(index(cline,'* parameter da').ne.0) go to 70
        end do
        continue
70      do ipar=1,npar
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9050) cline
          call linesplit(ifail,2)
          if(ifail.ne.0) go to 9100
          apar(ipar)=cline(left_word(1):right_word(1))
          call casetrans(apar(ipar),'lo')
          atrans=cline(left_word(2):right_word(2))
          call casetrans(atrans,'lo')
          if(atrans.eq.'log')then
            itrans(ipar)=1
          else if(atrans.eq.'none')then
            itrans(ipar)=0
          else if(atrans.eq.'fixed')then
            itrans(ipar)=-10000
          else if(atrans.eq.'tied')then
            itrans(ipar)=-1
          else
            call num2char(iline,aline)
            write(amessage,80) trim(aline),trim(pestfile)
80          format(' Incorrect PARTRANS entry at line ',a,' of PEST control file ',a,'.')
            go to 9890
          end if
        end do
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9250) cline
          call casetrans(cline,'lo')
          if(index(cline,'* observation gr').eq.0) cycle
          exit
        end do
        do i=1,nobsgp
          read(pestunit,'(a)',err=9000,end=9050) oldgroup(i)
          oldgroup(i)=adjustl(oldgroup(i))
          call casetrans(oldgroup(i),'lo')
        end do

! -- The regularisation information file is now read.

100	aprompt=' Enter name of regularisation information file: '
	call open_input_file(ifail,aprompt,regfile,regunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
          close(unit=pestunit)
          deallocate(apar,itrans,ipartype,iparpoint,aparnopoint,oldgroup,stat=ierr)
          if(ierr.ne.0) go to 9200
          write(6,*)
          go to 25
        end if
        afile=regfile
        iline=1
        read(regunit,*,err=9000,end=9050) numpoints
        allocate(point_id(numpoints),vario(numpoints,numpoints), &
        ipointnolink(numpoints),apointnolink(numpoints),stat=ierr)
        if(ierr.ne.0) go to 9150
        do i=1,numpoints
          iline=iline+1
          read(regunit,'(a)',err=9000,end=9050) point_id(i)
          call casetrans(point_id(i),'lo')
        end do
        do i=1,numpoints
          read(regunit,*,err=9450,end=9050) (vario(i,j),j=1,numpoints)
        end do
        close(unit=regunit)

! -- Information on the pilot point parameter families is now obtained.

        write(6,*)
140     write(6,150)
150     format(' How many pilot-point-based parameter families in PEST control file')
        write(6,149,advance='no')
149     format(' pertain to this regularisation information file: ')
        itemp=key_read(numpartype)
        if(escset.ne.0)then
          escset=0
          deallocate(point_id,vario,ipointnolink,apointnolink,stat=ierr)
          if(ierr.ne.0) go to 9200
          write(6,*)
          go to 100
        end if
        if(itemp.ne.0)then
          write(6,180)
180       format(' Data input error  - try again.')
          go to 140
        endif
        if(pos_test(numpartype,'no. of parameter families').ne.0) go to 140

        allocate(prefix(numpartype),weightmul(numpartype),ius(numpartype),  &
        reggroup(numpartype),piroot(numpartype),stat=ierr)
        if(ierr.ne.0) go to 9150

        itype=0
190      itype=itype+1
          if(itype.gt.numpartype) go to 250
195       call num2char(itype,anum)
          write(6,200) trim(anum)
200       format(/,' For family number ',a,':-')
          atemp=' '
209       write(6,210,advance='no')
210       format('   Enter parameter prefix (<Enter> if none): ')
          read(5,'(a)') atemp
          if(index(eschar,atemp(1:2)).ne.0) then
            if(itype.eq.1)then
              deallocate(prefix,weightmul,ius,reggroup,piroot,stat=ierr)
              if(ierr.ne.0) go to 9200
              write(6,*)
              go to 140
            else
              itype=itype-2
              go to 190
            end if
          end if
          prefix(itype)=atemp
          call casetrans(prefix(itype),'lo')
          if(itype.gt.1)then
            do j=1,itype-1
              if(prefix(j).eq.prefix(itype))then
                write(6,225)
225             format('   This prefix is already used - try again.')
                go to 209
              end if
            end do
          end if
220       write(6,230,advance='no')
230       format('   Apply uniform or geostatistical regularisation? [u/g]: ')
          read(5,'(a)') aus
          if(aus.eq.' ') go to 220
          if(index(eschar,aus).ne.0) then
            write(6,*)
            go to 209
          end if
          call casetrans(aus,'lo')
          if((aus.ne.'u').and.(aus.ne.'g')) go to 220
          if(aus.eq.'u')then
            ius(itype)=0
          else
            ius(itype)=1
          end if
240       write(6,245,advance='no')
245       format('   Enter weight multiplier: ')
          if(key_read(weightmul(itype)).ne.0) go to 240
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 220
          end if
          if(pos_test(weightmul(itype),'weight multiplier').ne.0) go to 240
248       write(6,246,advance='no')
246       format('   Enter new regularisation group name: ')
          read(5,'(a)') reggroup(itype)
          if(reggroup(itype).eq.' ') go to 248
          call casetrans(reggroup(itype),'lo')
          if(reggroup(itype)(1:2).eq.'e ')then
            write(6,*)
            go to 240
          end if
          if(reggroup(itype)(1:5).ne.'regul')then
            write(6,247)
247         format('   A regularisation group name must begin with "regul" - try again.')
            go to 248
          end if
          do ii=1,nobsgp
            if(reggroup(itype).eq.oldgroup(ii))then
              write(6,249)
249           format('   A group with this name is already defined in PEST control file',  &
              ' - try again.')
              go to 248
            end if
          end do
          if(itype.gt.1)then
            do ii=1,itype-1
              if(reggroup(itype).eq.reggroup(ii))then
                write(6,259)
259             format('   This name has already been supplied as a group name - try again.')
                go to 248
              end if
            end do
          end if
253       write(6,252,advance='no')
252       format('   Enter root name for new prior information: ')
          read(5,'(a)') atempr
          if(atempr.eq.' ') go to 253
          call casetrans(atempr,'lo')
          if(atempr(1:2).eq.'e ')then
            write(6,*)
            go to 248
          end if
          if(len_trim(atempr).gt.6)then
            write(6,258)
258         format('   Root name must be 6 characters or less - try again.')
            go to 253
          end if
          if(itype.gt.1)then
            do ii=1,itype-1
              if(piroot(ii).eq.atempr)then
                write(6,254)
254             format('   Name already used - try again.')
                go to 253
              end if
            end do
          end if
          piroot(itype)=atempr
        go to 190
250     continue

! -- Parameters in the PEST control file are now sorted into families.

        ipartype=0                     ! partype is an array
        iparpoint=0                    ! iparpoint is an array
        do itype=1,numpartype
          atemp=prefix(itype)
          if(atemp.ne.' ')then
            nbl=len_trim(atemp)
            do ipar=1,npar
              if(apar(ipar)(1:nbl).eq.atemp(1:nbl))then
                atemp1=apar(ipar)(nbl+1:)
                do j=1,numpoints
                  if(atemp1.eq.point_id(j))then
                    iparpoint(ipar)=j
                    ipartype(ipar)=itype
                    go to 270
                  end if
                end do
              end if
270           continue
            end do
          else
            do ipar=1,npar
              do j=1,numpoints
                if(apar(ipar).eq.point_id(j))then
                  iparpoint(ipar)=j
                  ipartype(ipar)=itype
                  go to 271
                end if
              end do
271         continue
            end do
          end if
        end do

        ilastflag=0
275     icount=count(iparpoint.eq.0)
        if(icount.ne.0)then
          write(amessage,280)
280       format(' The following parameters (and maybe others) in the PEST control ', &
          'file are not linked to a pilot point listed in the regularisation ', &
          'information file:-')
          call write_message(leadspace='yes')
          jpar=0
          lenchar=0
          do ipar=1,npar
            if(iparpoint(ipar).eq.0)then
              jpar=jpar+1
              lenchar=lenchar+1+len_trim(apar(ipar))
              if(lenchar.gt.LENELINE)then
                jpar=jpar-1
                go to 289
              end if
              aparnopoint(jpar)=apar(ipar)
            end if 
          end do
289       continue
          write(eline,300) (trim(aparnopoint(j)),j=1,jpar)
300       format(5000(1x,a))
          i=min(500,len_trim(eline))
          amessage=eline(1:i)
          call write_message()
325       write(6,320,advance='no')
320       format(' Is this ok? [y/n]: ')
          read(5,'(a)') ayn
          if(ayn.eq.' ') go to 325
          if(index(eschar,ayn).ne.0) then
            write(6,*)
            itype=numpartype-1
            go to 190
          end if
          call casetrans(ayn,'lo')
          if((ayn.ne.'y').and.(ayn.ne.'n')) go to 325
          if(ayn.eq.'n')then
            write(amessage,340)
340         format(' Edit input data prior to running this program again.')
            go to 9890
          end if
          ilastflag=1
        end if

! -- A check is made as to whether the regularisation information will be sufficient
!    for its purposes.

345     ipointnolink=0          !ipointnolink is an array
        do ipoint=1,numpoints
          do ipar=1,npar
            if(iparpoint(ipar).eq.ipoint) go to 346
          end do
          go to 350
346       continue
          icount=0
          do jpoint=1,numpoints
            if(ipoint.eq.jpoint)cycle
            if(vario(ipoint,jpoint).gt.-1.0e35)icount=icount+1
          end do
          if(icount.eq.0) ipointnolink(ipoint)=1
350     continue
        end do
        if(count(ipointnolink.ne.0).gt.0)then
          write(amessage,360) trim(regfile)
360       format(' The following pilot points in the regularisation information ', &
          'file are linked to parameters. However they are not linked to any other ', &
          'points in the variogram matrix contained in file ',a,'.')
          call write_message(leadspace='yes')
         jpoint=0
          do ipoint=1,numpoints
            if(ipointnolink(ipoint).ne.0)then
              jpoint=jpoint+1
              apointnolink(jpoint)=point_id(ipoint)
            end if
          end do
          write(eline,300) (trim(apointnolink(ipoint)),ipoint=1,jpoint)
          i=min(500,len_trim(eline))
          amessage=eline(1:i)
          call write_message()
370       write(6,320,advance='no')
          read(5,'(a)') ayn
          if(ayn.eq.' ') go to 370
          if(index(eschar,ayn).ne.0) then
            if(ilastflag.eq.1) then
              go to 275
            else
              write(6,*)
              itype=numpartype-1
              go to 190
            end if
          end if
          call casetrans(ayn,'lo')
          if((ayn.ne.'y').and.(ayn.ne.'n')) go to 370
          if(ayn.eq.'n')then
            write(amessage,380)
380         format(' Re-run PPK2FAC with a greater search radius to generate a ',&
            'new regularisation information file with more widespread point linkages.')
            go to 9890
          end if
          ilastflag=2
        end if

385     ipointnolink=0          !ipointnolink is an array
        do ipoint=1,numpoints
          do ipar=1,npar
            if(iparpoint(ipar).eq.ipoint) go to 390
          end do
          go to 400
390       continue
          icount=0
          do jpoint=1,numpoints
            if(ipoint.eq.jpoint)cycle
            if(vario(ipoint,jpoint).gt.-1.0e35)icount=icount+1
          end do
          if(icount.eq.1) ipointnolink(ipoint)=1
400     continue
        end do
        if(count(ipointnolink.ne.0).gt.0)then
          write(amessage,410) trim(regfile)
410       format(' The following pilot points in the regularisation information ', &
          'file are linked to parameters. However each is linked to only one other ', &
          'point in the variogram matrix contained in file ',a,'.')
          call write_message(leadspace='yes')
         jpoint=0
          do ipoint=1,numpoints
            if(ipointnolink(ipoint).ne.0)then
              jpoint=jpoint+1
              apointnolink(jpoint)=point_id(ipoint)
            end if
          end do
          write(eline,300) (trim(apointnolink(ipoint)),ipoint=1,jpoint)
          i=min(500,len_trim(eline))
          amessage=eline(1:i)
          call write_message()
420       write(6,320,advance='no')
          read(5,'(a)') ayn
          if(ayn.eq.' ') go to 420
          if(index(eschar,ayn).ne.0) then
            if(ilastflag.eq.1) then
              go to 275
            else if(ilastflag.eq.2) then
              go to 345
            else
              write(6,*)
              itype=numpartype-1
              go to 190
            end if
          end if
          call casetrans(ayn,'lo')
          if((ayn.ne.'y').and.(ayn.ne.'n')) go to 420
          if(ayn.eq.'n')then
            write(amessage,380)
            go to 9890
          end if
          ilastflag=3
        end if

! -- The parameters of the PEST control file are now checked.

425     iflag=0
        icount=count((itrans.eq.-10000).and.(iparpoint.ne.0))
        if(icount.ne.0) iflag=1
        icount=count((itrans.eq.-1).and.(iparpoint.ne.0))
        if(icount.ne.0) iflag=1

        if(iflag.ne.0)then
          write(amessage,430)
430       format(' There is at least one parameter in the PEST control file which ', &
          'is tied or fixed, yet is linked to a point in the regularisation ', &
          'information file. This will not be included in any regularisation ', &
          'prior information equations.')
          call write_message(leadspace='yes')
440       write(6,320,advance='no')
          read(5,'(a)') ayn
          if(ayn.eq.' ') go to 440
          if(index(eschar,ayn).ne.0) then
            if(ilastflag.eq.1) then
              go to 275
            else if(ilastflag.eq.2) then
              go to 345
            else if(ilastflag.eq.3)then
               go to 385
            else
              write(6,*)
              itype=numpartype-1
              go to 190
            end if
          end if
          call casetrans(ayn,'lo')
          if((ayn.ne.'y').and.(ayn.ne.'n')) go to 440
          if(ayn.eq.'n')then
            write(amessage,450)
450         format(' Then edit the PEST control file, changing the status of these ', &
            'parameters, prior to running PPKREG again.')
            go to 9890
          end if
          ilastflag=4
        end if

        do ipar=2,npar
          i=iparpoint(ipar)
          if(i.eq.0) go to 480
          do jpar=1,ipar-1
            j=iparpoint(jpar)
            if(j.eq.0) go to 475
            if(vario(i,j).gt.-1.0e35)then
              if(((itrans(ipar).eq.0).and.(itrans(jpar).eq.1)).or.  &
                 ((itrans(ipar).eq.1).and.(itrans(jpar).eq.0)))then
                   write(amessage,460)
460                format(' There are at least one pair of parameters cited ', &
                   'in the PEST control file which are linked by variogram ', &
                   'information in the regularisation information file for which ', &
                   'one parameter is log-transformed but the other is untransformed.')
                   go to 9890
              end if
            end if
475       continue
          end do
480     continue
        end do

! -- All data checks have now been made. Final user input is sought.

500     write(6,*)
        aprompt=' Enter name for new PEST control file: '
        call open_output_file(ifail,aprompt,pestoutfile,pestoutunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          if(ilastflag.eq.1) then
            go to 275
          else if(ilastflag.eq.2) then
            go to 345
          else if(ilastflag.eq.3) then
            go to 385
          else if(ilastflag.eq.4)then
            go to 425
          else
            write(6,*)
            itype=numpartype-1
            go to 190
          end if
        end if
520     write(6,530,advance='no')
530     format(' Enter target measurement objective function PHIMLIM: ')
        if(key_read(phimlim).ne.0) go to 520
        if(escset.eq.1) then
          write(6,*)
          escset=0
          close(unit=pestoutunit)
          go to 500
        end if
        if(pos_test(phimlim,'PHIMLIM').ne.0) go to 520

550     write(6,560,advance='no')
560     format(' Enter initial regularisation weight factor WFINIT: ')
        if(key_read(wfinit).ne.0) go to 550
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 520
        end if
        if(pos_test(wfinit,'WFINIT').ne.0) go to 550

        wfmin=wfinit*1.0e-10
        write(anum,565) wfmin
565     format(1pg10.4)
570     write(6,580,advance='no') trim(anum)
580     format(' Enter min. reg. weight factor WFMIN ( <Enter> if ',a,'): ')
        read(5,'(a)') atemp
        if(atemp.eq.' ') go to 590
        if(index(eschar,atemp(1:2)).ne.0) go to 550
        call char2num(ifail,atemp,rtemp)
        if(ifail.ne.0) then
          go to 570
        else
          if(pos_test(rtemp,'WFMIN').ne.0) go to 570
          wfmin=rtemp
        end if
590     continue
        if(wfmin.gt.wfinit)then
          write(6,592)
592       format(' WFMIN must be less than WFINIT - try again.')
          go to 570
        end if

        wfmax=wfinit*1.0e10
        write(anum,565) wfmax
600     write(6,610,advance='no') trim(anum)
610     format(' Enter max. reg. weight factor WFMAX ( <Enter> if ',a,'): ')
        read(5,'(a)') atemp
        if(atemp.eq.' ') go to 620
        if(index(eschar,atemp(1:2)).ne.0) go to 570
        call char2num(ifail,atemp,rtemp)
        if(ifail.ne.0) then
          go to 600
        else
          if(pos_test(rtemp,'WFMAX').ne.0) go to 600
          wfmax=rtemp
        end if
620     continue
        if(wfmax.lt.wfinit)then
          write(6,622)
622       format(' WFMAX must be greater than WFINIT - try again.')
          go to 600
        end if



! -- The PEST control file is now re-wound and the new one written.

         rewind(unit=pestunit,iostat=ierr)
         if(ierr.ne.0) then
           write(amessage,605) trim(pestfile)
605        format(' Error rewinding file ',a,'.')
           go to 9890
         end if

! -- The number of new items of prior information is worked out.

        icount=0
        do itype=1,numpartype
          do ipar=2,npar
            if(ipartype(ipar).ne.itype) go to 680
            if(itrans(ipar).lt.0) go to 680
            i=iparpoint(ipar)
            if(i.eq.0) go to 680
            do jpar=1,ipar-1
              if(ipartype(jpar).ne.itype) go to 675
              if(itrans(jpar).lt.0) go to 675
              j=iparpoint(jpar)
              if(j.eq.0) go to 675
              if(vario(i,j).lt.-1.0e35) go to 675
              icount=icount+1
675           continue
            end do
680       continue
          end do
685     continue
        end do

        npriornew=nprior+icount
        nobsgpnew=nobsgp+numpartype

! -- The first PEST control file is copied to the second PEST control file
! -- with changes made as appropriate.

        afile=pestfile
        do i=1,2
          read(pestunit,700,err=9000,end=9050) cline
700       format(a)
          write(pestoutunit,701,err=9300) trim(cline)
701       format(a)
        end do
        read(pestunit,700,err=9000,end=9050) cline
        write(pestoutunit,710,err=9300)
710     format(' restart regularisation')
        read(pestunit,700,err=9000,end=9050) cline
        write(pestoutunit,715) npar,nobs,npargp,npriornew,nobsgpnew
715     format(5i8)
        do
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
          call casetrans(cline,'lo')
          if(index(cline,'* observation gr').ne.0)go to 720
        end do
720     continue
        do itype=1,numpartype
          write(pestoutunit,701,err=9300) trim(reggroup(itype))
        end do
        if(nprior.ne.0)then
          do
            read(pestunit,700,err=9000,end=9350) cline
            write(pestoutunit,701,err=9300) trim(cline)
            call casetrans(cline,'lo')
            if(index(cline,'* prior info').ne.0) go to 730
          end do
730       continue
          do
            read(pestunit,700,err=9000,end=750) cline
            if(cline.eq.' ') go to 750
            dline=cline
            call casetrans(dline,'lo')
            dline=adjustl(dline)
            if(dline(1:6).eq.'* pred') go to 750
            if(dline(1:6).eq.'* regu') go to 750
            write(pestoutunit,701,err=9300) trim(cline)
          end do
750       continue
        else
          do
            read(pestunit,700,err=9000,end=9400) cline
            write(pestoutunit,701) trim(cline)
            call casetrans(cline,'lo')
            if(index(cline,'* model in').ne.0) go to 780
          end do
780       continue
          do i=1,ntplfle
            read(pestunit,700,err=9000,end=9050) cline
            write(pestoutunit,701,err=9300) trim(cline)
          end do
          do i=1,ninsfle
            read(pestunit,700,err=9000,end=9050) cline
            write(pestoutunit,701,err=9300) trim(cline)
          end do
          write(pestoutunit,790)
790       format('* prior information')
        end if

! -- The regularisation prior information is now written to the PEST control file.

        iprior=0
        do itype=1,numpartype
          do ipar=2,npar
            if(ipartype(ipar).ne.itype) go to 830
            if(itrans(ipar).lt.0) go to 830
            i=iparpoint(ipar)
            if(i.eq.0) go to 830
            do jpar=1,ipar-1
              if(ipartype(jpar).ne.itype) go to 820
              if(itrans(jpar).lt.0) go to 820
              j=iparpoint(jpar)
              if(j.eq.0) go to 820
              if(vario(i,j).lt.-1.0e35) go to 820
              iprior=iprior+1
              call num2char(iprior,aprior)
              aprior=trim(piroot(itype))//trim(aprior)
              if(ius(itype).eq.0)then
                weight=weightmul(itype)
              else
                weight=weightmul(itype)*sqrt(1.0/(2.0*vario(i,j)))
              end if
              if(itrans(ipar).eq.0)then
                write(pestoutunit,810) trim(aprior),trim(apar(ipar)),trim(apar(jpar)), &
                weight,trim(reggroup(itype))
810             format(a,2x,'1.0 * ',a,' - 1.0 * ',a,' = 0.0  ',1pg11.5,'  ',a)
              else
                write(pestoutunit,811) trim(aprior),trim(apar(ipar)),trim(apar(jpar)), &
                weight,trim(reggroup(itype))
811             format(a,2x,'1.0 * log(',a,') - 1.0 * log(',a,') = 0.0  ',1pg11.5,'  ',a)
              end if
820         continue
            end do
830       continue
          end do
        end do
        write(pestoutunit,870)
870     format('* regularisation')
        write(pestoutunit,880) phimlim,1.05*phimlim,0.1
880     format(3(1x,1pg14.7))
        write(pestoutunit,890) wfinit,wfmin,wfmax
890     format(3(1x,1pg12.5))
        write(pestoutunit,900)
900     format('   1.3 1.0e-2  1')
        close(unit=pestunit)
        close(unit=pestoutunit)
        write(6,840) trim(pestoutfile)
840     format(' - file ',a,' written ok.')
        go to 9900


9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(afile)
9060    format(' Unexpected end encountered to file ',a,'.')
        go to 9890
9070    write(amessage,9080) trim(pestfile)
9080    format(' Unexpected end to PEST control file ',a,' - looking for ', &
        '"* parameter data" section.')
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Insufficient entries on line ',a,' of file ',a,'.')
        go to 9890
9150    write(amessage,9160)
9160    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Memory management error - cannot de-allocate memory.')
        go to 9890
9250    write(amessage,9260) trim(pestfile)
9260    format(' Unexected end to PEST control file ',a,' - looking for ', &
        '"* observation groups" section.')
        go to 9890
9300    write(amessage,9310) trim(pestoutfile)
9310    format(' Error writing to new PEST control file ',a,'.')
        go to 9890
9350    write(amessage,9360) trim(pestfile)
9360    format(' Unexpected end to PEST control file ',a,' - looking for ', &
        '"* prior information" section.')
        go to 9890
9400    write(amessage,9410) trim(pestfile)
9410    format(' Unexpected end to PEST control file ',a,' - looking for ', &
        '"* model input/output" section.')
        go to 9890
9450    write(amessage,9460) trim(regfile)
9460    format(' Error reading variogram matrix from file ',a,'.')
        go to 9890



9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(itrans,apar,ipartype,iparpoint,aparnopoint,stat=ierr)
        deallocate(point_id,vario,ipointnolink,apointnolink,stat=ierr)
        deallocate(prefix,weightmul,ius,reggroup,piroot,stat=ierr)
        deallocate(iparpoint,oldgroup,stat=ierr)


end program ppkreg

! -- Ideally, if a parameter is log-transformed for estimation, then the variogram
! -- for regularisation should pertain to the log of the parameter - just through the
! -- way that the prior information is expressed. If a user has (for the sake of
! -- using MODFLOW 2000 to calculate derivatives, log-transformed K for estimation
! -- but not for regularisation, then he/she may consider running PPK2FAC again
! -- on the basis of variagrams of the logs of the parameters, just to build a
! -- regularisation information matrix based on log-parameters.


