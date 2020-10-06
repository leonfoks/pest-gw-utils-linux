!     Last change:  JD   18 Oct 2003    4:05 am


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
                                  lenchar,ipunit,dunit,ndata,ijustread,m,llo,l,cunit, &
                                  jj,nobsgpold,iend,ibeg
        integer                :: pestunit,regunit,pestoutunit
        integer, allocatable   :: itrans(:),ipartype(:),iparpoint(:),ipointnolink(:), &
                                  ius(:),isp(:),iuf(:),iud(:),n(:),used(:)
        real                   :: phimlim,wfinit,weight,wfmin,wfmax,rtemp,wmuli,wmulj
        real, allocatable      :: vario(:,:),prefvar(:,:),totcov(:),covmat(:,:)
        real, allocatable      :: weightmul(:)
        real, allocatable      :: a(:),b(:),c(:),maxweight(:),minweight(:),ww(:)
        double precision       :: ee,nn,sum,mindist,dtemp
        double precision, allocatable :: easting(:),northing(:),deast(:),dnorth(:)
        character (len=1)      :: ayn,aus,asp,auf,aug,aud
        character (len=15)     :: aline,pestmode,anum,aobsgp,atemp,atemp1, &
                                  aprior,aweight,atrans,atempr,aapoint
        character (len=200)    :: aprompt
        character (len=200)    :: pestfile,afile,regfile,pestoutfile,ppfile,datafile, &
                                  acovfile
        character (len=300)    :: dline
        character (len=5), allocatable  :: prefix(:)
        character (len=12), allocatable :: point_id(:)
        character (len=12), allocatable :: apar(:),aparnopoint(:),apointnolink(:), &
                                           reggroup(:),piroot(:),oldgroup(:)
        character (len=200), allocatable :: covfile(:),ecovfile(:)

! -- The following variables must be adjusted together.

        integer, parameter         :: LENELINE=8000
        character (len=LENELINE)   :: eline
                                                                                      

	write(amessage,5)
5	format(' Program PPKREG1 adds prior information and a regularisation section ',&
        'to a PEST control file whose parameters are based on pilot points. ', &
        'Regularisation weights are assigned based on information contained ', &
        'within a PPK2FAC1-generated regularisation information file.')
	call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        datafile=' '
        ijustread=0

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
        nobsgpold=nobsgp
        if(npar.eq.1)then
          write(amessage,60) trim(pestfile)
60        format(' There is only one parameter cited in PEST control file ',a, &
          '. More parameters than this are required if regularisation is to be ', &
          'applied.')
          go to 9890
        end if
        allocate(apar(npar),itrans(npar),ipartype(npar),iparpoint(npar), &
        aparnopoint(npar),oldgroup(nobsgp),ecovfile(nobsgp),stat=ierr)
        if(ierr.ne.0) go to 9150
        ecovfile=' '                      ! an array
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
          iline=iline+1
81        read(pestunit,'(a)',err=9000,end=9050) cline
          if(cline.eq.' ') go to 81
          call linesplit(ifail,1)
          oldgroup(i)=cline(left_word(1):right_word(1))
          call casetrans(oldgroup(i),'lo')
          call linesplit(ifail,2)
          if(ifail.eq.0)then
            cline=cline(left_word(2):)
            ibeg=1
            iend=len_trim(cline)
            call getfile(ifail,cline,ecovfile(i),ibeg,iend)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,82) trim(aline),trim(pestfile)
82            format(' Cannot read covariance matrix filename from line ',a,  &
              ' of file ',a)
              go to 9890
            end if
            call casetrans(ecovfile(i),'lo')      ! must not do in unix
          else
            ecovfile(i)=' '
          end if
        end do

! -- The regularisation information file is now read.

100	aprompt=' Enter name of regularisation information file: '
	call open_input_file(ifail,aprompt,regfile,regunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
          close(unit=pestunit)
          deallocate(apar,itrans,ipartype,iparpoint,aparnopoint,oldgroup,   &
          ecovfile,stat=ierr)
          if(ierr.ne.0) go to 9200
          write(6,*)
          go to 25
        end if
        afile=regfile
        iline=1
        read(regunit,*,err=9000,end=9050) numpoints
        allocate(point_id(numpoints),vario(numpoints,numpoints), &
        ipointnolink(numpoints),apointnolink(numpoints),  &
        easting(numpoints),northing(numpoints),ww(numpoints),totcov(numpoints),  &
        covmat(numpoints,numpoints),stat=ierr)
        if(ierr.ne.0) go to 9150
        do i=1,numpoints
          iline=iline+1
          read(regunit,'(a)',err=9000,end=9050) cline
          call linesplit(ifail,3)
          if(ifail.ne.0)then
            write(amessage,101) trim(regfile)
101         format(' Regularisation information file ',a,' does not appear to be ',  &
            'of the correct format. It appears to have been written by PPK2FAC rather ', &
            'than by PPK2FAC1.')
            go to 9890
          end if
          point_id(i)=cline(left_word(1):right_word(1))
          easting(i)=char2double(ifail,2)
          if(ifail.ne.0) go to 9000
          northing(i)=char2double(ifail,3)
          if(ifail.ne.0) go to 9000
          call casetrans(point_id(i),'lo')
        end do
        do i=1,numpoints
          read(regunit,*,err=9450,end=9050) (vario(i,j),j=1,numpoints),totcov(i)
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
          deallocate(point_id,vario,ipointnolink,apointnolink,easting,northing,ww,  &
          totcov,covmat,stat=ierr)
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
        reggroup(numpartype),piroot(numpartype),isp(numpartype),  &
        iuf(numpartype),prefvar(numpoints,numpartype),iud(numpartype),   &
        a(numpartype),b(numpartype),c(numpartype),n(numpartype),   &
        maxweight(numpartype),minweight(numpartype),covfile(numpartype),stat=ierr)
        if(ierr.ne.0) go to 9150

! -- Initialise prefvar.

        prefvar=-1.1e35

        itype=0
190     itype=itype+1
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
              deallocate(prefix,weightmul,ius,reggroup,piroot,isp,iuf,prefvar,iud,  &
              a,b,c,n,minweight,maxweight,covfile,stat=ierr)
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

1000      write(6,1010,advance='no')
1010      format('   Apply smoothness or preferred value regularisation? [s/p]: ')
          read(5,'(a)') asp
          if(asp.eq.' ') go to 1000
          if(index(eschar,asp).ne.0) then
            write(6,*)
            go to 209
          end if
          call casetrans(asp,'lo')
          if((asp.ne.'s').and.(asp.ne.'p')) go to 1000
          if(asp.eq.'s')then
            isp(itype)=0
          else
            isp(itype)=1
          end if

          if(isp(itype).eq.1)then
1030        write(6,1040,advance='no')
1040        format('   Use uniform preferred value, or read it from a file? [u/f]: ')
            read(5,'(a)') auf
            if(auf.eq.' ') go to 1030
            if(index(eschar,auf).ne.0) then
              write(6,*)
              go to 1000
            end if
            call casetrans(auf,'lo')
            if((auf.ne.'u').and.(auf.ne.'f')) go to 1030
            if(auf.eq.'u')then
              iuf(itype)=0
            else
              iuf(itype)=1
            end if

            if(iuf(itype).eq.0)then
1050          write(6,1060,advance='no')
1060          format('   Enter uniform preferred value: ')
              itemp=key_read(rtemp)
              if(escset.ne.0)then
                escset=0
                write(6,*)
                go to 1030
              end if
              if(itemp.ne.0)then
                write(6,1061)
1061            format(/,'   Data input error - try again.',/)
                go to 1050
              endif
              do i=1,numpoints
                prefvar(i,itype)=rtemp
              end do
            else
1070          write(6,1080,advance='no')
1080          format('   Enter pilot points file to read preferred values: ')
              read(5,'(a)') ppfile
              if(ppfile.eq.' ') go to 1070
              if((ppfile(1:2).eq.'e ').or.(ppfile(1:2).eq.'E '))then
                write(6,*)
                go to 1030
              end if
              ipunit=nextunit()
              open(unit=ipunit,file=ppfile,status='old',iostat=ierr)
              if(ierr.ne.0) go to 1070
              iline=0
              do
1090            iline=iline+1
                read(ipunit,'(a)',end=1150,err=9500) cline
                if(cline.eq.' ') go to 1090
                call linesplit(ifail,5)
                if(ifail.ne.0)then
                  call num2char(iline,aline)
                  write(amessage,1100) trim(aline),trim(ppfile)
1100              format(' Insufficient entries on line ',a,' of file ',a)
                  go to 9890
                end if
                aapoint=cline(left_word(1):right_word(1))
                call casetrans(aapoint,'lo')
                do i=1,numpoints
                  if(aapoint.eq.point_id(i)) go to 1110
                end do
                cycle
1110            continue
                prefvar(i,itype)=char2real(ifail,5)
                if(ifail.ne.0) go to 9500
              end do
1150          close(unit=ipunit)
              icount=0
              do i=1,numpoints
                if(prefvar(i,itype).lt.-1.0e35) icount=icount+1
              end do
              if(icount.gt.0)then
                write(amessage,1160) trim(ppfile)
1160            format(' No preferred values were supplied for the following ',  &
                'pilot points in file ',a,':-')
                call write_message(leadspace='yes')
                icount=0
                cline=' '
                do i=1,numpoints
                  if(prefvar(i,itype).lt.-1.0e35)then
                    icount=icount+1
                    cline=trim(cline)//'  '//trim(point_id(i))
                    if(len_trim(cline).gt.68)then
                      write(6,'(a)') trim(cline)
                      cline=' '
                      icount=0
                    end if
                  end if
                end do
                if(icount.ne.0) write(6,'(a)') trim(cline)
                go to 9900
              else
                write(6,1161) trim(ppfile)
1161            format('   - file ',a,' read ok.')
              end if
            end if
          end if

220       write(6,230,advance='no')
230       format('   Use weights of unity or geostatistical weighting? [u/g]: ')
          read(5,'(a)') aug
          if(aug.eq.' ') go to 220
          if(index(eschar,aug).ne.0) then
            write(6,*)
            go to 1000
          end if
          call casetrans(aug,'lo')
          if((aug.ne.'u').and.(aug.ne.'g')) go to 220
          if(aug.eq.'u')then
            ius(itype)=0
          else
            ius(itype)=1
          end if

1200      write(6,1210,advance='no')
1210      format('   Use uniform or data-density-dependent weight multiplier? [u/d]: ')
          read(5,'(a)') aud
          if(aud.eq.' ') go to 1200
          if(index(eschar,aud).ne.0) then
            write(6,*)
            go to 220
          end if
          call casetrans(aud,'lo')
          if((aud.ne.'u').and.(aud.ne.'d')) go to 1200
          if(aud.eq.'u')then
            iud(itype)=0
          else
            iud(itype)=1
          end if
          if(iud(itype).eq.0)then
240         write(6,245,advance='no')
245         format('   Enter uniform weight multiplier: ')
            if(key_read(weightmul(itype)).ne.0) go to 240
            if(escset.eq.1) then
              write(6,*)
              escset=0
              go to 1200
            end if
            if(pos_test(weightmul(itype),'weight multiplier').ne.0) go to 240
          else
1514        continue
            if(datafile.eq.' ')then
1511          write(6,1512,advance='no')
1512          format('   Enter name of data coordinates file: ')
              read(5,'(a)') datafile
              if(datafile.eq.' ') go to 1511
              if(index(eschar,datafile(1:2)).ne.0) then
                datafile=' '
                write(6,*)
                go to 1200
              end if
              dunit=nextunit()
              open(unit=dunit,file=datafile,status='old',iostat=ierr)
              if(ierr.ne.0)then
                datafile=' '
                go to 1511
              end if
              ndata=0
              do
                read(dunit,'(a)',end=1520) cline
                if(cline.eq.' ') cycle
                ndata=ndata+1
              end do
1520          continue
              if(ndata.eq.0)then
                write(amessage,1521) trim(datafile)
1521            format('   No data was found in file ',a,' - try again.')
                call write_message(leadspace='yes',endspace='yes')
                close(unit=dunit)
                datafile=' '
                go to 1511
              end if
              allocate(deast(ndata),dnorth(ndata),used(ndata),stat=ierr)
              if(ierr.ne.0) go to 9150
              rewind(unit=dunit)
              iline=0
              do i=1,ndata
1515            iline=iline+1
                read(dunit,'(a)') cline
                if(cline.eq.' ') go to 1515
                call linesplit(ifail,3)
                if(ifail.ne.0)then
                  call num2char(iline,aline)
                  write(amessage,1525) trim(aline),trim(datafile)
1525              format(' Insufficient entries on line ',a,' of file ',a)
                  go to 9890
                end if
                deast(i)=char2double(ifail,2)
                if(ifail.ne.0) go to 9600
                dnorth(i)=char2double(ifail,3)
                if(ifail.ne.0) go to 9600
              end do
              close(unit=dunit)
              call num2char(ndata,atemp)
              write(6,1530) trim(atemp),trim(datafile)
1530          format('   - ',a,' data points read from file ',a)
              ijustread=1
            else
              write(6,1531) trim(datafile)
1531          format('   Data read from file ',a,' is used for data density.')
              ijustread=0
            end if
            write(6,1560)
1560        format('   Weight multiplier function is "a+b*sum_over_closest_n(r**c)" :-')
1570        write(6,1580,advance='no')
1580        format('     Enter a: ')
            if(key_read(a(itype)).ne.0) go to 1570
            if(escset.ne.0)then
              escset=0
              write(6,*)
              if(ijustread.eq.1)then
                deallocate(deast,dnorth,used,stat=ierr)
                if(ierr.ne.0) go to 9200
                datafile=' '
                go to 1514
              else
                go to 1200
              end if
            end if
1590        write(6,1600,advance='no')
1600        format('     Enter b: ')
            if(key_read(b(itype)).ne.0) go to 1590
            if(escset.ne.0)then
              escset=0
              write(6,*)
              go to 1570
            end if
1610        write(6,1615,advance='no')
1615        format('     Enter n: ')
            if(key_read(n(itype)).ne.0) go to 1610
            if(escset.ne.0)then
              escset=0
              write(6,*)
              go to 1590
            end if
            if(n(itype).le.0)then
              write(6,1616)
1616          format(/,' Must be greater than zero - try again.',/)
              go to 1610
            end if
1620        write(6,1630,advance='no')
1630        format('     Enter c: ')
            if(key_read(c(itype)).ne.0) go to 1620
            if(escset.ne.0)then
              escset=0
              write(6,*)
              go to 1610
            end if
            if(c(itype).lt.0.0)then
              write(6,1676)
              go to 1620
            end if
1640        write(6,1650,advance='no')
1650        format('     Enter maximum allowable weight factor: ')
            if(key_read(maxweight(itype)).ne.0) go to 1640
            if(escset.ne.0)then
              escset=0
              write(6,*)
              go to 1620
            end if
            if(maxweight(itype).le.0.0)then
              write(6,1616)
              go to 1640
            end if
1660        write(6,1670,advance='no')
1670        format('     Enter minimum allowable weight factor: ')
            if(key_read(minweight(itype)).ne.0) go to 1660
	    if(escset.ne.0)then
              escset=0
              write(6,*)
              go to 1640
            end if
            if(minweight(itype).ge.maxweight(itype))then
              write(6,1675)
1675          format(/,' Minimum weight must be less than maximum weight - try again.',/)
              go to 1660
            end if
            if(minweight(itype).lt.0.0)then
              write(6,1676)
1676          format(/,' Must not be negative - try again.',/)
              go to 1660
            end if
          end if

248       write(6,246,advance='no')
246       format('   Enter new regularisation group name: ')
          read(5,'(a)') reggroup(itype)
          if(reggroup(itype).eq.' ') go to 248
          call casetrans(reggroup(itype),'lo')
          if(reggroup(itype)(1:2).eq.'e ')then
            write(6,*)
            if(datafile.ne.' ')then
              deallocate(deast,dnorth,used,stat=ierr)
              if(ierr.ne.0) go to 9200
              datafile=' '
            end if
            itype=itype-1
            go to 190
          end if
          if(reggroup(itype)(1:5).ne.'regul')then
            write(6,247)
247         format(/,'   A regularisation group name must begin with "regul" - try again.',/)
            go to 248
          end if
          do ii=1,nobsgpold
            if(reggroup(itype).eq.oldgroup(ii))then
              write(6,249)
249           format(/,'   This group is already used in PEST control file',  &
              ' - try again.',/)
              go to 248
            end if
          end do
          if(itype.gt.1)then
            do ii=1,itype-1
              if(reggroup(itype).eq.reggroup(ii))then
                write(6,259)
259             format(/,'   This name has already been supplied as a ',  &
                'group name - try again.',/)
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
          if((isp(itype).eq.1).and.(ius(itype).eq.1))then
1310        write(6,1320,advance='no')
1320        format('   Enter covariance matrix file for this prior information: ')
            read(5,'(a)') cline
            if(cline.eq.' ') go to 1310
            if((cline(1:2).eq.'E ').or.(cline(1:2).eq.'e '))then
              write(6,*)
              go to 253
            end if
            cline=adjustl(cline)
            ibeg=1
            iend=len_trim(cline)
            call getfile(ifail,cline,covfile(itype),ibeg,iend)
            call casetrans(covfile(itype),'lo')           ! not in unix
            do jj=1,nobsgpold
              if(covfile(itype).eq.ecovfile(jj))then
                write(amessage,1321) trim(pestfile)
1321            format(' This is the same as the name of a covariance matrix ', &
                'file already cited in PEST control file ',a,' - try again.')
                call write_message(leadspace='yes',endspace='yes')
                go to 1310
              end if
            end do
            if(itype.gt.1)then
              do jj=1,itype-1
                if(covfile(itype).eq.covfile(jj))then
                  write(amessage,1322)
1322              format(' This name is already used - try again.')
                  call write_message(leadspace='yes',endspace='yes')
                  go to 1310
                end if
              end do
            end if
          end if
          
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

        do itype=1,numpartype
          do ipar=2,npar
            if(ipartype(ipar).ne.itype) cycle
            i=iparpoint(ipar)
            if(i.eq.0) go to 480
            do jpar=1,ipar-1
              if(ipartype(jpar).ne.itype) cycle
              j=iparpoint(jpar)
              if(j.eq.0) go to 475
              if(vario(i,j).gt.-1.0e35)then
                if(((itrans(ipar).eq.0).and.(itrans(jpar).eq.1)).or.  &
                   ((itrans(ipar).eq.1).and.(itrans(jpar).eq.0)))then
                     write(amessage,460)
460                  format(' There are at least one pair of parameters cited ', &
                     'in the PEST control file which are linked by variogram ', &
                     'information in the regularisation information file for which ', &
                     'one parameter is log-transformed but the other is untransformed.')
                     go to 9890
                end if
              end if
475         continue
            end do
480       continue
          end do
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
          if(isp(itype).eq.0)then
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
675             continue
              end do
680         continue
            end do
          else
            do ipar=1,npar
              if(ipartype(ipar).ne.itype) cycle
              if(itrans(ipar).lt.0) cycle
              i=iparpoint(ipar)
              if(i.eq.0) cycle
              icount=icount+1
            end do
          end if
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
          if((isp(itype).eq.1).and.(ius(itype).eq.1))then
            call addquote(covfile(itype),acovfile)
            write(pestoutunit,702,err=9300) trim(reggroup(itype)),trim(acovfile)
702         format(1x,a,2x,a)
          else
            write(pestoutunit,702,err=9300) trim(reggroup(itype))
          end if
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

        do itype=1,numpartype
        
! -- If data density-dependent weighting is activated, we now work out the weight for each
!    pilot point.

          if(iud(itype).eq.0)then
            ww=1.0                 ! ww is an array
          else
            do i=1,numpoints
            
              do ipar=1,npar
                if(ipartype(ipar).ne.itype) cycle
                if(iparpoint(ipar).ne.i) cycle
                go to 1798
              end do
              cycle
            
1798          continue
              ee=easting(i)
              nn=northing(i)
              used=0            ! used is an array
              sum=0.0d0
              do m=1,min(n(itype),ndata)
                mindist=1.0d300
                llo=0
                do l=1,ndata
                  if(used(l).ne.0) cycle
                  dtemp=(ee-deast(l))*(ee-deast(l)) +    &
                        (nn-dnorth(l))*(nn-dnorth(l))
                  if(dtemp.lt.mindist)then
                    llo=l
                    mindist=dtemp
                  end if
                end do
                if(llo.ne.0)then
                  used(llo)=l
                  if(mindist.eq.0.0d0)then
                    continue
                  else 
                    if(c(itype).gt.0.0)then
                      sum=sum+sqrt(mindist)**c(itype)
                    else
                      sum=sum+1.0
                    end if
                  end if
                else
                  go to 1799
                end if
              end do
1799          continue
              ww(i)=a(itype)+b(itype)*sum
              if(ww(i).gt.maxweight(itype)) ww(i)=maxweight(itype)
              if(ww(i).lt.minweight(itype)) ww(i)=minweight(itype)
            end do
          end if
          
          iprior=0    !  make change in ppkreg too? chek that basenames are different
          if(isp(itype).eq.0)then
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
                  if(iud(itype).eq.0)then
                    weight=weightmul(itype)
                  else
                    weight=sqrt(ww(i)*ww(j))
                  end if
                else
                  if(iud(itype).eq.0)then
                    weight=weightmul(itype)*sqrt(1.0/(2.0*vario(i,j)))
                  else
                    weight=sqrt(1.0/(2.0*vario(i,j)))*sqrt(ww(i)*ww(j))
                  end if
                end if
                if(itrans(ipar).eq.0)then
                  write(pestoutunit,810) trim(aprior),trim(apar(ipar)),trim(apar(jpar)), &
                  weight,trim(reggroup(itype))
810               format(a,2x,'1.0 * ',a,' - 1.0 * ',a,' = 0.0  ',1pg11.5,'  ',a)
                else
                  write(pestoutunit,811) trim(aprior),trim(apar(ipar)),trim(apar(jpar)), &
                  weight,trim(reggroup(itype))
811               format(a,2x,'1.0 * log(',a,') - 1.0 * log(',a,') = 0.0  ',1pg11.5,'  ',a)
                end if
820           continue
              end do
830         continue
            end do
          else
            do ipar=1,npar
              if(ipartype(ipar).ne.itype) cycle
              if(itrans(ipar).lt.0) cycle
              i=iparpoint(ipar)
              if(i.eq.0) cycle
              if(itrans(ipar).eq.0)then
                rtemp=prefvar(i,itype)
              else
                rtemp=prefvar(i,itype)
                if(rtemp.le.0.0)then
                  write(amessage,1710)
1710              format(' A negative or zero preferred regularisation value was supplied for ',  &
                  'a log-transformed '  &
                  'pilot point parameter. PPKREG1 cannot continue execution.')
                  go to 9890
                end if
                rtemp=log10(rtemp)
              end if
              iprior=iprior+1
              call num2char(iprior,aprior)
              aprior=trim(piroot(itype))//trim(aprior)
              if(ius(itype).eq.0)then
                if(iud(itype).eq.0)then
                  weight=weightmul(itype)
                else
                  weight=ww(i)
                end if
              else
                weight=1.0
              end if
              if(itrans(ipar).eq.0)then
                write(pestoutunit,1720) trim(aprior),trim(apar(ipar)),rtemp,   &
                weight,trim(reggroup(itype))
1720            format(a,2x,'1.0 * ',a,' = ',1pg12.5,'  ',1pg12.5,'  ',a)
              else
                write(pestoutunit,1721) trim(aprior),trim(apar(ipar)),rtemp,  &
                weight,trim(reggroup(itype))
1721            format(a,2x,'1.0 * log(',a,')  = ',1pg12.5,'  ',1pg12.5,'  ',a)
              end if
            end do

! -- The covariance matrix file is written if appropriate.

            if((ius(itype).eq.1).and.(isp(itype).eq.1))then
              cunit=nextunit()
              open(unit=cunit,file=covfile(itype),action='write',iostat=ierr)
              if(ierr.ne.0)then
                write(amessage,1810) trim(covfile(itype))
1810            format(' Cannot open file ',a,' to write covariance matrix.')
                go to 9890
              end if
              ii=0
              do ipar=1,npar
                if(ipartype(ipar).ne.itype) cycle
                if(itrans(ipar).lt.0) cycle
                i=iparpoint(ipar)
                if(i.eq.0) cycle
                ii=ii+1
                if(iud(itype).eq.0)then
                  wmuli=weightmul(itype)
                else
                  wmuli=ww(i)
                end if
                jj=0
                do jpar=1,ipar
                  if(ipartype(jpar).ne.itype) cycle
                  if(itrans(jpar).lt.0) cycle
                  j=iparpoint(jpar)
                  if(j.eq.0) cycle
                  jj=jj+1
                  if(iud(itype).eq.0)then
                    wmulj=weightmul(itype)
                  else
                    wmulj=ww(j)
                  end if
                  if(vario(i,j).gt.-1.0e35)then
                    covmat(ii,jj)=(totcov(i)-vario(i,j))/(wmuli*wmulj)
                  else
                    covmat(ii,jj)=0.0
                  end if
                  covmat(jj,ii)=covmat(ii,jj)
                end do
              end do
              do i=1,ii
                write(cunit,1830) (covmat(i,j),j=1,ii)
1830            format(8(1x,1pg14.7))
              end do
              write(6,1831) trim(covfile(itype))
1831          format(' - file ',a,' written ok.')
              close(unit=cunit)
            end if
          end if
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
9500    call num2char(iline,aline)
        write(amessage,9510) trim(aline),trim(ppfile)
9510    format(' Error reading line ',a,' of file ',a)
        go to 9890
9600    call num2char(iline,aline)
        write(amessage,9610) trim(aline),trim(datafile)
9610    format(' Cannot read data point coordinates from line ',a,' of file ',a)
        go to 9890


9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(itrans,apar,ipartype,iparpoint,aparnopoint,stat=ierr)
        deallocate(point_id,vario,ipointnolink,apointnolink,totcov,covmat,stat=ierr)
        deallocate(easting,northing,stat=ierr)
        deallocate(prefix,weightmul,ius,reggroup,piroot,isp,iuf,iud,stat=ierr)
        deallocate(iparpoint,oldgroup,stat=ierr)
        deallocate(prefvar,stat=ierr)
        deallocate(deast,dnorth,used,stat=ierr)
        deallocate(a,b,c,n,minweight,maxweight,covfile,ecovfile,ww,stat=ierr)


end program ppkreg

! -- Ideally, if a parameter is log-transformed for estimation, then the variogram
! -- for regularisation should pertain to the log of the parameter - just through the
! -- way that the prior information is expressed. If a user has (for the sake of
! -- using MODFLOW 2000 to calculate derivatives, log-transformed K for estimation
! -- but not for regularisation, then he/she may consider running PPK2FAC again
! -- on the basis of variagrams of the logs of the parameters, just to build a
! -- regularisation information matrix based on log-parameters.








