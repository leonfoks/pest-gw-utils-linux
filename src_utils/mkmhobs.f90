!     Last change:  JD   15 Sep 2006    0:17 am
program mkmhobs

! -- Program MKMHOBS writes a MODFLOW heads observation input file, a MODFLOW heads DATA instruction
!    file, and a PEST building block file based on measurements supplied in a bore sample file.

	use defn
	use inter

	implicit none

        integer, parameter :: MAXUNIT=100
        integer            :: iunit(MAXUNIT)

        logical  :: active,lexist
        integer  :: ifail,idate,iheader,i,j,k,iline,idis,iu,n,ibas,ierr,itemp,kk,jj,  &
                    istp,ntotdays,jindex,nindex,numsamp1
        integer  :: icount,jcount,kcount
        integer  :: lloc,istart,istop,ixsec,ifrefm,cols,itt
        integer  :: dds,mms,yys,hhhs,mmms,ssss,ndaystart,nsecstart,ndays,nsecs,ndayfin,nsecfin
        integer  :: itotstep,ntotstep,irow,icol,nrow,ncol,nnrow,nncol,nh,nnh,nlay,nper,iper
        integer  :: iu1,iunit1,iunit2,hobunit,datunit,basunit,sampunit,nameunit,obsunit,obsdatunit,  &
                    insunit,pestunit,disunit,locat_dis,locat_bas
        real     :: weight,dweight,tweight,tdweight,rtemp,hobdry
        double precision :: day_convert,r_day_convert
        double precision :: totim,delt,one,totdays,dtemp,dtemp1,toffset,refval,ptotim
	type (modelgrid) gridspec

        character*1   :: aso,aad,at,aa
        character*10  :: aboreold,abore,anum,acount
        character*15  :: adate,atime,aline,aper
        character*20  :: obsname
        character*30  :: atype
        character*200 :: sampfile,namefile,obsfile,obsdatfile,insfile,pestfile,disfile,  &
                         namefilekeep,basfile
        character*200 :: cfile,bfile

        integer, allocatable          :: layer(:),numsamp(:),index_bore(:),icellrow(:),icellcol(:)
        integer, allocatable          :: nstp(:),issflg(:),laycbd(:)
        integer, allocatable          :: iarray(:,:)
        real, allocatable             :: rvector(:),rarray(:,:)
        real, allocatable             :: perlen(:),tsmult(:)
        real, allocatable             :: roff(:),coff(:)
        real, allocatable             :: svalue(:)
        double precision, allocatable :: east(:),north(:)
        double precision, allocatable :: ttotim(:)
        double precision, allocatable :: soffset(:)

        character*24 bname
        CHARACTER*24 ANAME(5)
        DATA ANAME(1) /'                    DELR'/
        DATA ANAME(2) /'                    DELC'/
        DATA ANAME(3) /'TOP ELEVATION OF LAYER 1'/
        DATA ANAME(4) /'  MODEL LAYER BOTTOM EL.'/
        DATA ANAME(5) /'BOT. EL. OF QUASI-3D BED'/


        write(amessage,5)
5       format(' Program MKMHOBS writes a MODFLOW input heads observation file, ', &
       'an instruction file for a MODFLOW heads DATA output file, and a PEST ',  &
       'building block file based on measurements supplied in a bore sample file.')
        call write_message(leadspace='yes',endspace='yes')

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          call write_message
          go to 9900
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          call write_message
          go to 9900
        endif
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
          'settings.fig')
          call write_message
          go to 9900
        end if

        call readfig(gridspec%specfile,bore_coord_file,sampfile)
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

        write(6,*)
30      call read_bore_coord_file(ifail, &
        ' Enter name of bore coordinates file: ')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 10
        end if
	do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) go to 80
	end do
	go to 99
80	write(amessage,90) trim(bore_coord_file)
90	format(' The following bores were not provided with layer numbers in ',&
	'bore coordinates file ',a,' ----->')
	call write_message(leadspace='yes')
	j=2
	imessage=0
	write_bores: do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) then
	    write(amessage(j:),'(a)') trim(bore_coord_id(i))
	    j=j+11
	    if(j.ge.71) then
	      call write_message(increment=1)
	      if(imessage.gt.12) goto 9900
	      j=2
	    end if
	  end if
	end do write_bores
	if(j.ne.2) call write_message
	go to 9900

99      write(6,*)
100     call read_bore_list_file(ifail, &
       ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 30
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
	        write(amessage,110) trim(bore_list_file)
110		format(' Execution of program MKMHOBS cannot proceed as ',&
		'there are multiple occurrences of the same bore in bore ',&
		'listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	allocate(east(num_bore_list), north(num_bore_list),           &
	layer(num_bore_list),roff(num_bore_list),coff(num_bore_list), &
	numsamp(num_bore_list), index_bore(num_bore_list),            &
        icellrow(num_bore_list),icellcol(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,130)
130	  format(' Insufficient memory available to continue execution.')
	  go to 9890
	end if
        index_bore=0            ! an array
        icellrow=0              ! an array
        icellcol=0              ! an array
        layer=0                 ! an array

	imessage=0
	list: do i=1,num_bore_list
	  do j=1,num_bore_coord
	    if(bore_coord_id(j).eq.bore_list_id(i)) then
	      east(i)=bore_coord_east(j)
	      north(i)=bore_coord_north(j)
	      layer(i)=bore_coord_layer(j)
	      cycle list
	    end if
	  end do
	  write(amessage,140) trim(bore_list_id(i)),trim(bore_list_file),&
	  trim(bore_coord_file)
140	  format(' No coordinates for bore ',a,' from bore listing ',&
	  'file ',a,' are provided in bore coordinates file ',a)
	  if(imessage.eq.0) write(6,*)
	  call write_message(increment=1)
	end do list
	if(imessage.ne.0) go to 9900

        do i=1,num_bore_list
          call cell_coordinates(gridspec,east(i),north(i),icellrow(i),icellcol(i),roff(i),coff(i))
        end do
        if(any(icellrow.eq.-999))then
          write(6,141)
141       format(/,' At least one bore named in bore listing file is outside of model grid.')
          write(6,142)
142       format(' A list of offending bores follows:-')
          do i=1,num_bore_list
            if(icellrow(i).eq.-999)then
              write(6,143) trim(bore_list_id(i))
143           format(5x,a)
            end if
          end do
          write(6,*)
          go to 9900
        end if

!        open(unit=75,file='debug.dat')                                             !debug
!        do i=1,num_bore_list                                                       !debug
!          write(75,155) bore_list_id(i), icellrow(i),icellcol(i),roff(i),coff(i)   !debug
!155       format(1x,a,t15,i5,2x,i5,2x,f10.3,2x,f10.3)                              !debug
!        end do                                                                     !debug
!        call flush(75)                                                             !debug
!        write(75,*)                                                                !debug

        write(6,*)
146     call open_named_input_file(ifail,          &
        ' Enter name of bore sample file: ',sampfile,sampunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          write(6,*)
          deallocate(east,north,layer,roff,coff,numsamp,index_bore,icellrow,icellcol,stat=ierr)
          if(ierr.ne.0) go to 9200
          go to 100
        end if

160     write(6,165,advance='no')
165     format(' Assign observations to end of stress period or occurrence time? [s/o]: ')
        read(5,'(a)') aso
        call casetrans(aso,'lo')
        if(aso.eq.'e')then
          close(unit=sampunit)
          write(6,*)
          go to 146
        end if
        if((aso.ne.'s').and.(aso.ne.'o')) go to 160

170     write(6,175,advance='no')
175     format(' Use absolutes or differences (i.e. drawdowns) in calibration process? [a/d]: ')
        read(5,'(a)') aad
        call casetrans(aad,'lo')
        if(aad.eq.'e')then
          write(6,*)
          go to 160
        end if
        if((aad.ne.'a').and.(aad.ne.'d')) go to 170
        if(aad.eq.'d')then
          write(amessage,176)
176       format(' MKMHOBS does not presently allow this option.')
          go to 9890
        end if

200     continue
        if(aad.eq.'a')then
          write(6,210,advance='no')
210       format(' Enter value for all weights: ')
        else
          write(6,211,advance='no')
211       format(' Enter value for head weights: ')
        end if
        i=key_read(weight)
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 170
        else if(i.eq.-1) then
          go to 200
        else if(i.ne.0) then
          write(6,220)
220       format(' Illegal input  - try again.')
          go to 200
        end if
        if(weight.lt.0.0) then
          write(6,230)
230       format(' Weight must not be less than zero  - try again.')
          go to 200
        end if
        dweight=weight

214     continue
        if(aad.eq.'d')then
          write(6,212,advance='no')
212       format(' Enter value for drawdown weights: ')
          i=key_read(dweight)
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 200
          else if(i.eq.-1) then
            go to 214
          else if(i.ne.0) then
            write(6,220)
            go to 214
          end if
          if(dweight.lt.0.0) then
            write(6,230)
            go to 214
          end if
        end if


        write(6,*)
150     call open_input_file(ifail,          &
        ' Enter name of MODFLOW name file: ',namefile,nameunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
	  escset=0
	  write(6,*)
          if(aad.eq.'a')then
            go to 200
          else
            go to 214
          end if
        end if

310     if(datespec.eq.1) then
          write(6,320,advance='no')
320       format(' Enter simulation starting date [dd/mm/yyyy]: ')
        else
          write(6,321,advance='no')
321       format(' Enter simulation starting date [mm/dd/yyyy]: ')
        end if
        read(5,'(a)') adate
        if(adate.eq.' ') go to 310
        adate=adjustl(adate)
        if(index(eschar,adate(1:2)).ne.0) then
          write(6,*)
          close(unit=nameunit)
          go to 150
        end if
        call char2date(ifail,adate,dds,mms,yys)
        if(ifail.ne.0)then
          write(6,340)
340       format(' Illegal date  - try again.')
          go to 310
        end if

360	write(6,370,advance='no')
370	format(' Enter simulation starting time [hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 360
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) then
          write(6,*)
          go to 310
        end if
	call char2time(ifail,atime,hhhs,mmms,ssss)
	if(ifail.ne.0)then
	  write(6,380)
380	  format(' Illegal time  - try again.')
	  go to 360
	end if
	ndaystart=numdays(1,1,1970,dds,mms,yys)
	nsecstart=numsecs(0,0,0,hhhs,mmms,ssss)


250     write(6,260,advance='no')
260     format(' Enter time units used by model (yr/day/hr/min/sec) [y/d/h/m/s]: ')
        read(5,'(a)') at
        if(at.eq.' ') go to 250
        if(index(eschar,at).ne.0) then
          write(6,*)
          go to 360
        end if
        call casetrans(at,'lo')
        if(at.eq.'s') then
          day_convert=1.0d0/86400.0d0
        else if(at.eq.'m') then
          day_convert=1.0d0/1440.0d0
        else if(at.eq.'h') then
          day_convert=1.0d0/24.0d0
        else if(at.eq.'d') then
          day_convert=1.0d0
        else if(at.eq.'y') then
          day_convert=365.25d0
        else
          go to 250
        end if
        r_day_convert=1.0d0/day_convert

! -- The names of MKMHOBS output files are now acquired.

400	write(6,*)
410	call open_output_file(ifail, &
	' Enter name for MODFLOW OBS file: ',obsfile,obsunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 250
	end if

420	call open_output_file(ifail, &
	' Enter name for MODFLOW observation data file: ',obsdatfile,obsdatunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
          close(unit=obsunit)
	  go to 400
	end if

430	call open_output_file(ifail, &
	' Enter name for corresponding PEST instruction file: ',insfile,insunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
          close(unit=obsdatunit)
	  go to 420
	end if

450	call open_output_file(ifail, &
	' Enter name for PEST building block file: ',pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
          close(unit=insunit)
	  go to 430
	end if

        write(6,*)
451     write(6,452,advance='no')
452     format(' Enter value for MODFLOW to use for dry observation cells: ')
        i=key_read(hobdry)
        if(escset.ne.0)then
          escset=0
          write(6,*)
          close(unit=pestunit)
          go to 450
        else if(i.eq.-1) then
          go to 451
        else if(i.ne.0) then
          write(6,220)
          go to 451
        end if

! -- The name file is read and the names of some other MODFLOW input files obtained.

        write(6,505) trim(namefile)
505     format(/,' - processing MODFLOW name file ',a,'...')
        iline=0
        idis=0
        ibas=0
        iu=0
        do
          iline=iline+1
          read(nameunit,'(a)',end=600) cline
          cline=adjustl(cline)
          if(cline(1:1).eq.'#') cycle
          if(cline.eq.' ') cycle
          call linesplit(ifail,3)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,510) trim(aline),trim(namefile)
510         format(' Insufficient entires on line ',a,' of MODFLOW name file ',a,'.')
            go to 9890
          end if
          atype=cline(left_word(1):right_word(1))
          call casetrans(atype,'lo')
          if(atype.eq.'hob')then
            write(amessage,520) trim(namefile)
520         format(' Name file ',a,' already cites a HOB file. Use of MKMHOBS assumes that ',  &
            'no HOB file is already cited in the name file.')
            go to 9890
          end if
          iu=iu+1
          if(iu.gt.MAXUNIT)then
            write(amessage,523)
523         format(' Increase MAXUNIT and re-compile program.')
            go to 9890
          end if
          iunit(iu)=char2int(ifail,2)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,521) trim(aline),trim(namefile)
521         format(' Cannot read unit number from second data columun at line ',a,  &
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
            write(amessage,530) trim(aline),trim(namefile)
530         format(' Improper filename at line ',a,' of name file ',a,'.')
            go to 9890
          end if
          if(aa.eq.' ')then
            bfile=bfile(1:n)
          else
            bfile=bfile(2:n)
          end if
          call casetrans(bfile,'lo')
          cfile=obsdatfile
          call casetrans(cfile,'lo')
          if(bfile.eq.cfile)then
            call num2char(iline,aline)
            write(amessage,540) trim(obsdatfile),trim(aline),trim(namefile)
540         format(' Requested new observation data file ',a,' already cited at line ',a,' of ',  &
            'MODFLOW name file ',a,'.')
            go to 9890
          end if
          if(atype.eq.'dis')then
            if(idis.ne.0)then
              write(amessage,545) trim(namefile)
545           format(' More than one discretisation file named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            disfile=bfile
            locat_dis=iunit(iu)
            idis=1
          else if(atype(1:3).eq.'bas')then
            if(ibas.ne.0)then
              write(amessage,546) trim(namefile)
546           format(' More than one basic package input file named in MODFLOW name file ',a,'.')
              go to 9890
            end if
            basfile=bfile
            locat_bas=iunit(iu)
            ibas=1
          end if
        end do
600     continue
        close(unit=nameunit)
        write(6,602) trim(namefile)
602     format(' - file ',a,' processed ok.')

! -- The discretisation file is now opened and read.

        write(6,605) trim(disfile)
605     format(/,' - reading MODFLOW discretisation file ',a,'...')
        disunit=nextunit()
        open(unit=disunit,file=disfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,610) trim(disfile),trim(namefile)
610       format(' Cannot open discretisation file ',a,' cited in MODFLOW ',  &
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
          write(amessage,611) trim(aline),trim(disfile)
611       format(' Insufficient entries on line ',a,' of MODFLOW discretisation file ',a,'.')
          go to 9890
        end if
        nlay=char2int(ifail,1)
        if(ifail.ne.0) go to 9070
        nnrow=char2int(ifail,2)
        if(ifail.ne.0) go to 9070
        nncol=char2int(ifail,3)
        if(ifail.ne.0) go to 9070
        nper=char2int(ifail,4)
        if(ifail.ne.0) go to 9070
        if((nnrow.le.0).or.(nncol.le.0).or.(nlay.le.0).or.(nper.le.0)) go to 9070
        if((nnrow.ne.nrow).or.(nncol.ne.ncol))then
          write(amessage,622)
622       format(' Grid dimensions in discretisation file do not match those in ',  &
          'grid specification file.')
          go to 9890
        end if
        allocate(laycbd(nlay),stat=ierr)
        if(ierr.ne.0) go to 9200
        read(disunit,*,iostat=ierr) (laycbd(i),i=1,nlay)
        if(ierr.ne.0)then
          write(amessage,630) trim(disfile)
630       format(' Error reading LAYCBD array from MODFLOW discretisation file ',a,'.')
          go to 9890
        end if
        itemp=max(ncol,nrow)
        allocate(rvector(itemp),rarray(ncol,nrow),iarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(perlen(nper),nstp(nper),tsmult(nper),issflg(nper))
        if(ierr.ne.0) go to 9200

! -- Read the delr and delc vectors.

        call u1drel(ifail,rvector,aname(1),ncol,disunit,0,locat_dis)
        if(ifail.ne.0) go to 9890
        call u1drel(ifail,rvector,aname(2),nrow,disunit,0,locat_dis)
        if(ifail.ne.0) go to 9890

! -- Read the top elevation of layer 1.

        call u2drel(ifail,rarray,aname(3),nrow,ncol,0,disunit,0,locat_dis)
        if(ifail.ne.0) go to 9890

! -- Read the bottom elevations.

        do k=1,nlay
          kk=k
          call u2drel(ifail,rarray,aname(4),nrow,ncol,kk,disunit,0,locat_dis)
          if(ifail.ne.0) go to 9890
          if(laycbd(k).ne.0)then
            kk=k
            call u2drel(ifail,rarray,aname(5),nrow,ncol,kk,disunit,0,locat_dis)
            if(ifail.ne.0) go to 9890
          end if
        end do

! -- Read stress period timing data.

        do iper=1,nper
          read(disunit,'(a)',err=9400,end=9400) cline
          lloc=1
          call urword(cline,lloc,istart,istop,3,itemp,perlen(iper),0,disunit)
          call urword(cline,lloc,istart,istop,2,nstp(iper),rtemp,0,disunit)
          call urword(cline,lloc,istart,istop,3,itemp,tsmult(iper),0,disunit)
          call urword(cline,lloc,istart,istop,1,itemp,rtemp,0,disunit)
          if(cline(istart:istop).eq.'TR')then
            issflg(iper)=0
          else if(cline(istart:istop).eq.'SS')then
            issflg(iper)=1
          else
            call num2char(iper,aper)
            write(amessage,640) trim(aper),trim(disfile)
640         format(' SSFLAG must be "TR" or "SS" for stress period ',a,' in MODFLOW discretisation file ',a,'.')
            go to 9890
          end if
          if(nstp(iper).le.0)then
            call num2char(iper,aper)
            write(amessage,650) trim(aper),trim(disfile)
650         format(' Stress period ',a,' has no time time steps according ',   &
            'to information supplied in file ',a,'.')
            go to 9890
          end if
          if(perlen(iper).le.0.0)then
            call num2char(iper,aper)
            write(amessage,660) trim(aper),trim(disfile)
660         format(' Period length must be greater than zero for stress period ',a,    &
            ' cited in file ',a,'.')
            go to 9890
          end if
          if(tsmult(iper).le.0.0)then
            call num2char(iper,aper)
            write(amessage,670) trim(aper),trim(disfile)
670         format(' Time period multiplier must be greater than zero for stress ', &
            'period ',a,' cited in file ',a,'.')
            go to 9890
          end if
        end do
        close(unit=disunit)

! -- TOTIM for the end of each time step is now evaluated.

        ntotstep=0
        do iper=1,nper
          ntotstep=ntotstep+nstp(iper)
        end do
        allocate(ttotim(ntotstep),stat=ierr)
        if(ierr.ne.0) go to 9200
        totim=0.0d0
        ptotim=0.0d0
        itotstep=0
        do iper=1,nper
          delt=perlen(iper)/dble(nstp(iper))
          one=1.0d0
          if(tsmult(iper).ne.one)then
             delt=perlen(iper)*(one-tsmult(iper))/(one-tsmult(iper)**nstp(iper))
          end if
          do istp=1,nstp(iper)
            if(istp.ne.1) delt=tsmult(iper)*delt
            totim=totim+delt
            itotstep=itotstep+1
            ttotim(itotstep)=totim
          end do
          ptotim=ptotim+perlen(iper)
          ttotim(itotstep)=ptotim
        end do
        write(6,680) trim(disfile)
680     format(' - discretisation file ',a,' processed ok.')

!        write(75,*) ntotstep                           !debug
!        do itotstep=1,ntotstep                         !debug
!          write(75,*) ttotim(itotstep)              !debug
!        end do                                         !debug
!        call flush(75)                            !debug


! -- We now evaluate the days and seconds of the end of the simulation relative to 1/1/1970.

        totdays=ttotim(ntotstep)*day_convert
        ntotdays=totdays
        ndayfin=ndaystart+ntotdays
        nsecfin=nsecstart+nint(dble(totdays-ntotdays)*86400.0d0)
        if(nsecfin.ge.86400)then
          ndayfin=ndayfin+1
          nsecfin=nsecfin-86400
        end if

!        write(75,*)                                  !debug
!        write(75,*) ' ndaystart = ',ndaystart        !debug
!        write(75,*) ' nsecstart = ',nsecstart        !debug
!        write(75,*) ' ndayfin   = ',ndayfin          !debug
!        write(75,*) ' nsecfin   = ',nsecfin          !debug
!        call flush(75)                               !debug

! -- The next task is to read IBOUND arrays to see if any bore is in an inactive cell.

        write(6,705) trim(basfile)
705     format(/,' - reading MODFLOW basic package input file ',a,'...')
        basunit=nextunit()
        open(unit=basunit,file=basfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,710) trim(basfile),trim(namefile)
710       format(' Cannot open basic package input file ',a,' cited in MODFLOW ',  &
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
          write(amessage,720)
720       format(' MODFLOW model is running in cross-section mode; this is not permitted with ',  &
          'present version of MKMHOBS.')
          go to 9890
        end if
        bname='          BOUNDARY ARRAY'
        do k=1,nlay
          kk=k
          call u2dint(ifail,iarray,bname,nrow,ncol,kk,basunit,0,locat_bas)
          if(ifail.ne.0) go to 9890
          do i=1,num_bore_list
            if(layer(i).eq.k)then
              irow=icellrow(i)
              icol=icellcol(i)
              if(iarray(icol,irow).eq.0)then
                icellrow(i)=-999999
              end if
            end if
          end do
        end do
        if(any(icellrow.eq.-999999))then
          write(6,730)
730       format(/,' At least one bore named in bore listing file is in an inactive cell.')
          write(6,740)
740       format(' A list of offending bores follows:-')
          do i=1,num_bore_list
            if(icellrow(i).eq.-999999)then
              write(6,743) trim(bore_list_id(i))
743           format(5x,a)
            end if
          end do
          write(6,*)
          go to 9900
        end if
        close(unit=basunit)
        write(6,750) trim(basfile)
750     format(' - basic package file ',a,' read ok.')

! -- The bore sample file must now be read a first time to work out the number of head observations.

        write(6,805) trim(sampfile)
805     format(/,' - reading bore sample file ',a,'...')
        nh=0
        aboreold=' '
        active=.false.
        iline=0
        imessage=0
        numsamp=0     ! an array
        jindex=0
        read_sample: do
          iline=iline+1
          read(sampunit,'(a)',err=9300,end=830) cline
          cols=5
          call linesplit(ifail,5)
          if(ifail.lt.0) cycle read_sample
          if(ifail.gt.0) then
            call linesplit(ifail,4)
            if(ifail.ne.0) go to 9350
            cols=4
          end if
          if(right_word(1)-left_word(1).gt.9) go to 9370
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if(abore.ne.aboreold)then
	    aboreold=abore
	    if(active) numsamp(i)=numsamp1
	    do i=1,num_bore_list
	      if(abore.eq.bore_list_id(i)) go to 620
	    end do
	    active=.false.
	    cycle read_sample
620	    numsamp1=0
	    active=.true.
	    index_bore(i)=jindex+1
	  else
	    if(.not.active) cycle read_sample
	  end if
          call read_rest_of_sample_line(ifail,cols,ndays,nsecs,dtemp, &
          iline,sampfile)
          if(ifail.ne.0) go to 9900
	  if(dtemp.lt.-1.0e38) cycle read_sample
          if((ndays.lt.ndaystart).or.                                  &
             ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) cycle read_sample
          if((ndays.gt.ndayfin).or.                                    &
             ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin))) cycle read_sample
	  jindex=jindex+1
	  numsamp1=numsamp1+1
          nh=nh+1
	end do read_sample

830     if(active) numsamp(i)=numsamp1
	nindex=jindex
        allocate(soffset(nindex),svalue(nindex),stat=ierr)
        if(ierr.ne.0) go to 9200

!        write(75,*)                                                     !debug
!        do i=1,num_bore_list                                            !debug
!          write(75,831) trim(bore_list_id(i)), numsamp(i), index_bore(i)  !debug
!831       format(1x,a,t20,i5,i6)                                          !debug
!        end do                                                          !debug
!        write(75,*) ' nindex = ', nindex                                !debug
!        write(75,*) ' nh     = ', nh                                    !debug
!        call flush(75)                                                  !debug

! -- The bore sample file is now re-read and the MODFLOW observation file simultaneously re-written.

        rewind(unit=sampunit)
        write(6,840) trim(obsfile)
840     format(' - writing MODFLOW head observation file ',a,'...')
        write(6,850) trim(insfile)
850     format(' - writing instruction file ',a,'...')
        write(6,851) trim(pestfile)
851     format(' - writing PEST building block file ',a,'...')
        write(pestunit,852)
852     format('* observation groups')
        write(pestunit,853)
853     format('head')
        if(aad.eq.'d')then
          write(pestunit,854)
854       format('headdiff')
        end if
        write(pestunit,855)
855     format('* observation data')

	aboreold=' '
	active=.false.
	iline=0
	imessage=0
	jindex=0
        read_sample1: do
          iline=iline+1
          read(sampunit,'(a)',err=9300,end=2000) cline
          cols=5
          call linesplit(ifail,5)
          if(ifail.lt.0) cycle read_sample1
          if(ifail.gt.0) then
            call linesplit(ifail,4)
            cols=4
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if(abore.ne.aboreold)then
	    aboreold=abore
            do i=1,num_bore_list
              if(abore.eq.bore_list_id(i)) go to 920
            end do
            active=.false.
            cycle read_sample1
920         numsamp1=0
            active=.true.
          else
            if(.not.active) cycle read_sample1
          end if
          call read_rest_of_sample_line(ifail,cols,ndays,nsecs,dtemp, &
          iline,sampfile)
          if(dtemp.lt.-1.0e38) cycle read_sample1
          if((ndays.lt.ndaystart).or.                                  &
             ((ndays.eq.ndaystart).and.(nsecs.lt.nsecstart))) cycle read_sample1
          if((ndays.gt.ndayfin).or.                                    &
             ((ndays.eq.ndayfin).and.(nsecs.gt.nsecfin))) cycle read_sample1
          toffset=(ndays-ndaystart)+dble(nsecs-nsecstart)/86400.0d0
          toffset=toffset*r_day_convert
	  numsamp1=numsamp1+1
          jindex=jindex+1
          if(aso.eq.'o')then
            soffset(jindex)=toffset
          else
            dtemp1=1.0d30
            jj=0
            do j=1,ntotstep
              if(abs(toffset-ttotim(j)).lt.dtemp1) then
                dtemp1=abs(toffset-ttotim(j))
                jj=j
              end if
            end do
            soffset(jindex)=ttotim(jj)
          end if
          svalue(jindex)=dtemp
        end do read_sample1

! -- We now evaluate the total number of observations.

2000    continue
        if(aso.eq.'o')then
          nnh=nh
        else
          nnh=0
          do i=1,num_bore_list
            if(numsamp(i).eq.1)then
              nnh=nnh+1
            else if(numsamp(i).gt.1)then
              nnh=nnh+1
              do j=index_bore(i)+1,index_bore(i)+numsamp(i)-1
                if(.not.equals(soffset(j),soffset(j-1))) nnh=nnh+1
              end do
            end if
          end do
        end if

! -- Now we find unit numbers for the new head observation file and for the observation
!    data output file.

        icount=0
        do i=10,100
          do j=1,iu
            if(iunit(j).eq.i) go to 1010
          end do
          icount=icount+1
          if(icount.eq.1)then
            iu1=i
          else
            exit
          end if
1010      continue
        end do

        hobunit=iu1
        datunit=i

! -- The observation file is now written.

        write(obsunit,870)
870     format('# MODFLOW head observation file prepared by MKMHOBS.')
        write(obsunit,880)
880     format('#')
        write(obsunit,890) nnh,0,2,datunit,hobdry
890     format(4i10,1x,1pg14.7)
        write(obsunit,900)
900     format('     1.0     1.0')
        write(insunit,901)
901     format('pif #')
        write(insunit,902)
902     format('l1')

        do i=1,num_bore_list
          if(numsamp(i).eq.0) cycle
          if(numsamp(i).eq.1)then
            obsname=trim(bore_list_id(i))//'_1'
            call casetrans(obsname,'lo')
            toffset=soffset(index_bore(i))
            write(obsunit,930) trim(obsname),layer(i),icellrow(i),icellcol(i),1,toffset,   &
            roff(i),coff(i),svalue(index_bore(i)),1.0/weight,1,1
930         format(1x,a,t23,i6,i6,i6,i6,1x,1pg13.6,1x,1pg14.7,1x,1pg14.7,1x,1pg14.7,1x,1pg13.6,i6,i6)
            write(insunit,935) trim(obsname)
935         format('l1  !',a,'!')
            write(pestunit,936) trim(obsname),svalue(index_bore(i)),weight
936         format(1x,a,t23,1pg14.7,2x,1pg13.6,2x,'head')
          else
            icount=1
            do j=index_bore(i)+1,index_bore(i)+numsamp(i)-1
              if(.not.equals(soffset(j),soffset(j-1)))icount=icount+1
            end do
            if(icount.eq.1)then
              obsname=trim(bore_list_id(i))//'_1'
              call casetrans(obsname,'lo')
              jcount=1
              dtemp=svalue(index_bore(i))
              tweight=weight
              do j=index_bore(i)+1,index_bore(i)+numsamp(i)-1
                if(equals(soffset(j),soffset(j-1)))then
                  dtemp=dtemp+svalue(j)
                  tweight=tweight+weight
                  jcount=jcount+1
                end if
              end do
              dtemp=dtemp/jcount
              toffset=soffset(index_bore(i))
              write(obsunit,930) trim(obsname),layer(i),icellrow(i),icellcol(i),1,toffset,   &
              roff(i),coff(i),dtemp,1.0/tweight,1,1
              write(insunit,935) trim(obsname)
              write(pestunit,936) trim(obsname),dtemp,tweight
            else
              obsname=bore_list_id(i)
              call casetrans(obsname,'lo')
              write(obsunit,930) trim(obsname),layer(i),icellrow(i),icellcol(i),-icount,   &
              soffset(index_bore(i)),roff(i),coff(i),svalue(index_bore(i)),1.0/weight,1,1
              if(aad.eq.'a')then
                itt=1
              else
                itt=2
              end if
              write(obsunit,940) itt
940           format(i6)
              jcount=0
              do j=index_bore(i),index_bore(i)+numsamp(i)-1
                if((j.eq.index_bore(i)).or.    &
                  ((j.gt.index_bore(i)).and.(.not.equals(soffset(j),soffset(j-1)))))then
                  kcount=0
                  tweight=0.0d0
                  tdweight=0.0d0
                  dtemp=0.0d0
                  do k=index_bore(i),index_bore(i)+numsamp(i)-1
                    if(equals(soffset(k),soffset(j)))then
                      kcount=kcount+1
                      dtemp=dtemp+svalue(k)
                      tweight=tweight+weight
                      tdweight=tdweight+dweight
                    end if
                  end do
                  dtemp=dtemp/kcount
                  jcount=jcount+1
                  call num2char(jcount,acount)
                  obsname=trim(bore_list_id(i))//'_'//trim(acount)
                  call casetrans(obsname,'lo')
                  write(obsunit,950) trim(obsname),1,soffset(j),dtemp,    &
                  1.0/tweight,1.0/tdweight,1,1
950               format(1x,a,t23,i6,1x,1pg13.6,1x,1pg14.7,1x,1pg13.6,1x,1pg13.6,i6,i6)
                  if(aad.eq.'a')then
                    write(pestunit,936) trim(obsname),dtemp,tweight
                    write(insunit,935) trim(obsname)
                  else
                    if(jcount.eq.1)then
                      write(pestunit,936) trim(obsname),dtemp,tweight
                      refval=dtemp
                      write(insunit,935) trim(obsname)
                    else
                      write(pestunit,951) trim(obsname),dtemp-refval,tdweight
951                   format(1x,a,t23,1pg14.7,2x,1pg13.6,2x,'headdiff')
                      write(insunit,952) trim(obsname)
952                   format('l1  !dum!  !dum! w w !',a,'!')
                    end if
                  end if
                end if
              end do
            end if
          end if
        end do

        close(unit=insunit)
        close(unit=obsunit)
        write(6,960) trim(obsfile)
960     format(' - file ',a,' written ok.')
        write(6,960) trim(insfile)

! -- Now we finish off the PEST building block file with the observation file and instruction name.

        write(pestunit,970)
970     format('* model input/output')
        call addquote(insfile,bfile)
        call addquote(obsdatfile,cfile)
        write(pestunit,980) trim(bfile),trim(cfile)
980     format(1x,a,3x,a)
        close(unit=pestunit)
        write(6,990) trim(pestfile)
990     format(' - file ',a,' written ok.')

! -- The existing name file is copied to another file for safe keeping.

        i=0
        do
          i=i+1
          call num2char(i,anum)
          anum=adjustl(anum)
          namefilekeep=trim(namefile)//'.k'//trim(anum)
          inquire(file=namefilekeep,exist=lexist)
          if(.not.lexist)exit
        end do
        call system('copy '//trim(namefile)//' '//trim(namefilekeep)//' > nul')

! -- The new name file is now written

        iunit1=nextunit()
        open(unit=iunit1,file=namefilekeep,status='old')
        iunit2=nextunit()
        open(unit=iunit2,file=namefile)
        do
          read(iunit1,'(a)',end=1100) cline
          write(iunit2,'(a)') trim(cline)
        end do
1100    continue
        close(unit=iunit1)
        write(iunit2,'(a)') '#'
        write(iunit2,1110)
1110    format('# The following two lines were added by MKMHOBS.')
        call addquote(obsfile,bfile)
        write(iunit2,1120) hobunit,trim(bfile)
1120    format('hob',t22,i3,5x,a)
        call addquote(obsdatfile,bfile)
        write(iunit2,1130) datunit,trim(bfile)
1130    format('DATA',t22,i3,5x,a)
        close(unit=iunit2)
        write(6,1140) trim(namefile)
1140    format(' - new name file ',a,' written ok.')
        write(6,1150) trim(namefilekeep)
1150    format(' - old name file copied to ',a,'.')

        go to 9900

9000    write(amessage,9010) trim(disfile)
9010    format(' Error reading MODFLOW discretisation file ',a,'.')
        go to 9890
9020    write(amessage,9030)
9030    format(' Error in de-allocating memory.')
        go to 9890
9050    write(amessage,9060) trim(disfile)
9060    format(' Premature end encountered to MODFLOW discretisation file ',a,'.')
        go to 9890
9070    call num2char(iline,aline)
        write(amessage,9080) trim(aline),trim(disfile)
9080    format(' Error in data at line ',a,' of MODFLOW discretisation file ',a,'.')
        go to 9890
9100    write(amessage,9110) trim(basfile)
9110    format(' Error reading MODFLOW basic package input file ',a,'.')
        go to 9890
9150    write(amessage,9160) trim(basfile)
9160    format(' Premature end encountered to MODFLOW basic package input file ',a,'.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890
9300    call num2char(iline,aline)
        write(amessage,9310) trim(aline),trim(sampfile)
9310    format(' Error encountered in reading line ',a,' of bore sample file ',a,'.')
        go to 9890
9350    call num2char(iline,aline)
        write(amessage,9360) trim(aline),trim(sampfile)
9360    format(' Insufficient entries on line ',a,' of bore sample file ',a,'.')
        go to 9890
9370    call num2char(iline,aline)
        write(amessage,9380) trim(aline),trim(sampfile)
9380    format(' Bore name has more than 10 characters at line ',a,' of bore sample file ',a,'.')
        go to 9890
9400    write(amessage,9410) trim(disfile)
9410    format(' Error reading stress period timing data from discretisation file ',a,'.')
        go to 9890


9890    call write_message(leadspace='yes')
9900    call close_files
        call free_bore_mem
        call free_grid_mem(gridspec)

        if(allocated(laycbd)) deallocate(laycbd,stat=ierr)
        if(allocated(layer)) deallocate(layer,stat=ierr)
        if(allocated(numsamp)) deallocate(numsamp,stat=ierr)
        if(allocated(index_bore)) deallocate(index_bore,stat=ierr)
        if(allocated(icellrow)) deallocate(icellrow,stat=ierr)
        if(allocated(icellcol)) deallocate(icellcol,stat=ierr)
        if(allocated(nstp)) deallocate(nstp,stat=ierr)
        if(allocated(issflg)) deallocate(issflg,stat=ierr)
        if(allocated(rvector)) deallocate(rvector,stat=ierr)
        if(allocated(rarray)) deallocate(rarray,stat=ierr)
        if(allocated(iarray)) deallocate(iarray,stat=ierr)
        if(allocated(perlen)) deallocate(perlen,stat=ierr)
        if(allocated(tsmult)) deallocate(tsmult,stat=ierr)
        if(allocated(roff)) deallocate(roff,stat=ierr)
        if(allocated(coff)) deallocate(coff,stat=ierr)
        if(allocated(svalue)) deallocate(svalue,stat=ierr)
        if(allocated(east)) deallocate(east,stat=ierr)
        if(allocated(north)) deallocate(north,stat=ierr)
        if(allocated(ttotim)) deallocate(ttotim,stat=ierr)
        if(allocated(soffset)) deallocate(soffset,stat=ierr)
	write(6,*)

end program mkmhobs




      SUBROUTINE U1DREL(ifail,A,ANAME,JJ,IN,IOUT,locat_dis)
!     ******************************************************************
!     ROUTINE TO INPUT 1-D REAL DATA MATRICES
!       A IS ARRAY TO INPUT
!       ANAME IS 24 CHARACTER DESCRIPTION OF A
!       JJ IS NO. OF ELEMENTS
!       IN IS INPUT UNIT
!       IOUT IS OUTPUT UNIT
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------

      use defn, only : amessage
      implicit integer (i-n), real (a-h, o-z)

      integer ifail,locat_dis
      CHARACTER*24 ANAME
      DIMENSION A(JJ)
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
8888     format(' "EXTERNAL" not allowed in array header when using MKMHOBS.')
         ifail=1
         return
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
!!!!         WRITE(IOUT,15) LOCAT,FNAME
!!!!   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1))
         ICLOSE=1
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
           'from another file using OPEN/CLOSE functionality when using MKMHOBS.')
           ifail=1
           return
         end if
      END IF
!
!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN)
         IF(LOCAT.GT.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
         END IF
      END IF
!
!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.GT.0) GO TO 90
!
!4A-----LOCAT <0 OR =0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
      DO 80 J=1,JJ
   80 A(J)=CNSTNT
!!!!      WRITE(IOUT,3) ANAME,CNSTNT
!!!!    3 FORMAT(1X,/1X,A,' =',1P,G14.6)
      RETURN
!
!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
   90 CONTINUE
!!!!      WRITE(IOUT,5) ANAME,LOCAT,FMTIN
!!!!    5 FORMAT(1X,///11X,A,/
!!!!     1       1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A20)
      IF(FMTIN.EQ.'(FREE)') THEN
      READ(LOCAT,*) (A(J),J=1,JJ)
      ELSE
         READ(LOCAT,FMTIN) (A(J),J=1,JJ)
      END IF
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
!
!5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
      ZERO=0.
      IF(CNSTNT.EQ.ZERO) GO TO 120
      DO 100 J=1,JJ
  100 A(J)=A(J)*CNSTNT
!
!6------IF PRINT CODE (IPRN) =0 OR >0 THEN PRINT ARRAY VALUES.
120   CONTINUE
!!!!      IF(IPRN.EQ.0) THEN
!!!!         WRITE(IOUT,1001) (A(J),J=1,JJ)
!!!!1001     FORMAT((1X,1PG12.5,9(1X,G12.5)))
!!!!      ELSE IF(IPRN.GT.0) THEN
!!!!         WRITE(IOUT,1002) (A(J),J=1,JJ)
!!!!1002     FORMAT((1X,1PG12.5,4(1X,G12.5)))
!!!!      END IF
!
!7------RETURN
      RETURN
!
!8------CONTROL RECORD ERROR.
500   WRITE(amessage,502) ANAME
502   FORMAT(' ERROR READING ARRAY CONTROL RECORD FOR ',A,'.')
      ifail=1
      return
!!!!      WRITE(IOUT,'(1X,A)') CNTRL
!!!!      CALL USTOP(' ')
      END




      SUBROUTINE U2DREL(ifail,A,ANAME,II,JJ,K,IN,IOUT,locat_dis)
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
      CHARACTER*24 ANAME
      DIMENSION A(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*16 TEXT
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
8888     format(' "EXTERNAL" not allowed in array header when using MKMHOBS.')
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
         READ(CNTRL,1,ERR=500) LOCAT,CNSTNT,FMTIN,IPRN
    1    FORMAT(I10,F10.0,A20,I10)
         if(locat.eq.locat_dis)then
           locat=in
         else
           write(amessage,8889)
8889       format(' MODFLOW arrays must be read from same file as package input file or ',   &
           'from another file using OPEN/CLOSE functionality when using MKMHOBS.')
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
8888     format(' "EXTERNAL" not allowed in array header when using MKMHOBS.')
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
           'from another file using OPEN/CLOSE functionality when using MKMHOBS.')
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




! -- Point out in documentation that EXTERNAL files not allowed when reading discretisation
!    or basic package input file. Also, we use MODFLOW's reading functionality which does
!    not give as good error messages.


