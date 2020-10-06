!     Last change:  J     8 Mar 2005   10:04 pm
program genreg

! -- Program GENREG adds many different kinds of regularisation to pilot point
!    disposed throughout a two or three-dimensional model domain.

	use defn
	use inter

	implicit none

        integer, parameter    :: MAXOBS=10000
        integer, parameter    :: MAXREGPAR=5000
        integer, parameter    :: MAXNEWGROUP=100
        integer, parameter    :: MAXLINKAGE=100000
        integer, parameter    :: MAXWARN=1000
        integer               :: ifail,pestunit,ierr,npar,nobs,npargp,nprior,  &
                                 nobsgp,ipar,i,genunit,outunit,iline,regpresent, &
                                 fpos,lpos,iblock,nblock,rblock,iprior,foundprior, &
                                 ntpfle,ninsfle,tempunit,blocktype,nb,nt,nrpar,iobs, &
                                 lt1,ibeg,iend,nf,irpar,numnewgroup,anynewprior,  &
                                 mobs,ifound,ilog,iparcount,jrpar,jmindist,ii, &
                                 nrpar_1,nrpar_2,ilog_1,ilog_2,nr1,nr2,ir1,ir2, &
                                 jprior,jp,iii,nf_1,nf_2,iwarn,paramindex,paramtrans
        integer               :: unassigned_int
        integer               :: linkage1(MAXLINKAGE),linkage2(MAXLINKAGE)

        real                  :: unassigned_real,unassigned_thresh
        real                  :: ww,diffrat

        double precision      :: distmin,dist,anis_angle,cosang,sinang,sepx,sepy, &
                                 asepx,asepy,mindist,eee,nnn,distmin1,distmin2,wtd,wt, &
                                 smalldist,east1,north1,east2,north2

        character (len=10)    :: aline,atrans
        character (len=12)    :: priorname,aapar1,aapar2
        character (len=20)    :: aprior,pestmode,ablock
        character (len=30)    :: atemp1,atemp
        character (len=200)   :: aprompt,pestfile,afile,genfile,outfile,tempfile, &
                                 bfile

        integer, allocatable               :: itrans(:)
        character (len=12), allocatable    :: apar(:),oldgroup(:)

! -- Regularisation variables from PEST control file.

        integer    :: memsave, conjgrad, iregadj, linreg, cgitnlim
        real       :: phimlim, phimaccept, fracphim, cgrtol, &
                      wfinit, wfmin, wfmax, wffac, wftol

! -- Entries in REGULARISATION block.

        integer    :: reg_memsave, reg_conjgrad, reg_iregadj, reg_linreg, &
                      reg_cgitnlim
        real       :: reg_phimlim, reg_phimaccept, reg_fracphim, reg_cgrtol, &
                      reg_wfinit, reg_wfmin, reg_wfmax, reg_wffac, reg_wftol

! -- Entries in REGSPEC block.

        integer    :: spec_type, reg_type, max_pilot_points, min_pilot_points,   &
                      val_type, weight_type, weight_obs_dist, equality_type,     &
                      diffrat_type,diffrat_val_type,warn_less_min
        real       :: value, weight, weight_multiplier, obs_dist_a, obs_dist_b, &
                      obs_dist_c, obs_dist_minwt, obs_dist_maxwt, search_radius, &
                      weight_sep_a, weight_sep_b, weight_sep_c, weight_sep_anis_bearing, &
                      weight_sep_anis_ratio, weight_sep_maxwt, weight_sep_minwt
        character (len=8)    :: family_prefix, family_prefix_1, family_prefix_2
        character (len=12)   :: reg_group,paramname
        character (len=200)  :: pilot_points_filename, obs_coord_filename, &
                                pilot_points_filename_1,pilot_points_filename_2, &
                                pilot_points_filename_dr

! -- Pilot point and observation coordinate data; note that observations are only used
!    from the point of view of regularisation weights calculation.

        integer            :: parindex(MAXREGPAR),parindex_1(MAXREGPAR),  &
                              parindex_2(MAXREGPAR)
        real               :: val_par(MAXREGPAR),val_par_1(MAXREGPAR),  &
                              val_par_2(MAXREGPAR)
        double precision   :: east_obs(MAXOBS),north_obs(MAXOBS)
        double precision   :: east_par(MAXREGPAR), north_par(MAXREGPAR),  &
                              ddist(MAXREGPAR),east_par_1(MAXREGPAR),     &
                              east_par_2(MAXREGPAR),north_par_1(MAXREGPAR), &
                              north_par_2(MAXREGPAR)
        character (len=12) :: arpar(MAXREGPAR),newgroup(MAXNEWGROUP),  &
                              arpar_1(MAXREGPAR),arpar_2(MAXREGPAR),  &
                              parwarn(MAXWARN)

! -- Processing begins


	write(amessage,5)
5       format(' Program GENREG adds regularisation prior information to a ', &
        'PEST control file in order to assist in the estimation of pilot ', &
        'point parameters for a two- or three-dimensional model domain.')
        call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        iprior=0
        foundprior=0
        unassigned_real=-1.1e37
        unassigned_thresh=unassigned_real*(1.0 - 2.0*epsilon(1.0))
        unassigned_int=-99999999
        regpresent=0
        rblock=0
        nblock=0
        numnewgroup=0
        tempfile='t###.###'

! -- Default values for regularisation control variables.

        phimlim=unassigned_real
        phimaccept=unassigned_real
        fracphim=0.0
        memsave=0
        conjgrad=0
        cgrtol=1.0e-5
        cgitnlim=500
        wfinit=1.0
        wfmin=1.0e-10
        wfmax=1.0e10
        wffac=1.3
        wftol=1.0e-2
        linreg=0
        iregadj=1

        reg_phimlim=unassigned_real
        reg_phimaccept=unassigned_real
        reg_fracphim=unassigned_real
        reg_memsave=unassigned_int
        reg_conjgrad=unassigned_int
        reg_cgrtol=unassigned_real
        reg_cgitnlim=unassigned_int
        reg_wfinit=unassigned_real
        reg_wfmin=unassigned_real
        reg_wfmax=unassigned_real
        reg_wffac=unassigned_real
        reg_wftol=unassigned_real
        reg_linreg=unassigned_int
        reg_iregadj=unassigned_int

! -- The PEST control file is partially read to obtain the names and status
!    of the parameters featured in it.

10	aprompt=' Enter name of existing PEST control file: '
	call open_input_file(ifail,aprompt,pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

! -- The GENREG control file is now opened.

25	aprompt=' Enter name of GENREG control file: '
	call open_input_file(ifail,aprompt,genfile,genunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          close(unit=pestunit)
          write(6,*)
          go to 10
        end if

! -- The new PEST control file is opened.

131     aprompt = ' Enter name for new PEST control file: '
        call open_output_file(ifail,aprompt,outfile,outunit)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1)then
          escset=0
          close(unit=genunit)
          write(6,*)
          go to 25
        end if

! -- The existing PEST control file is now read.

        call addquote(pestfile,afile)
        iline=1
        read(pestunit,26,err=9000,end=9050) cline
26      format(a)
        cline=adjustl(cline)
        call casetrans(cline,'lo')
        if(cline(1:3).ne.'pcf')then
          write(amessage,30) trim(afile)
30        format(' File ',a,' is not a PEST control file - ',   &
          '"pcf" header line not found.')
          go to 9890
        end if
        iline=iline+1
        read(pestunit,26,err=9000,end=9050) cline
        iline=iline+1
        read(pestunit,26,err=9000,end=9050) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        pestmode=cline(left_word(2):right_word(2))
        call casetrans(pestmode,'lo')
        if(pestmode.eq.'prediction')then
          write(amessage,50) trim(afile)
50        format(' In file ',a,' PEST is asked to run in predictive ',  &
          'analysis mode; GENREG requires that PEST be run in regularisation ', &
          'mode.')
          go to 9890
        else if(index(pestmode,'regul').eq.0)then
          write(amessage,51) trim(afile)
51        format(' PEST must be asked to run in regularisation mode ', &
          'in PEST control file ',a,'.')
          go to 9890
        end if
        iline=iline+1
        read(pestunit,*,err=9000,end=9050) npar,nobs,npargp,nprior,nobsgp
        if(npar.eq.1)then
          write(amessage,60) trim(afile)
60        format(' NPAR is supplied as 1 on line 4 of PEST control file ',a, &
          '. More parameters than this are required if regularisation is to be ', &
          'applied.')
          go to 9890
        end if
        allocate(apar(npar),itrans(npar),oldgroup(nobsgp),stat=ierr)
        if(ierr.ne.0) go to 9150
        iline=iline+1
        read(pestunit,*,err=9000,end=9050) ntpfle,ninsfle
        if((ntpfle.le.0).or.(ninsfle.le.0)) go to 9170
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9070) cline
          call casetrans(cline,'lo')
          if(index(cline,'* parameter da').ne.0) go to 70
        end do
70      continue
        do ipar=1,npar
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
            write(amessage,80) trim(aline),trim(afile)
80          format(' Incorrect PARTRANS entry at line ',a,' of PEST control file ',a,'.')
            go to 9890
          end if
        end do
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9090) cline
          call casetrans(cline,'lo')
          if(index(cline,'* observation gr').eq.0) cycle
          exit
        end do
        do i=1,nobsgp
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9050) cline
          call linesplit(ifail,1)
          oldgroup(i)=cline(left_word(1):right_word(1))
          call casetrans(oldgroup(i),'lo')
        end do

! -- The contents of the "regularisation" section are read (if this section is present).

        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9200) cline
          call casetrans(cline,'lo')
          if(index(cline,'* regul').ne.0) go to 130
        end do
130     continue
        regpresent=1
        iline=iline+1
        call num2char(iline,aline)
        read(pestunit,'(a)',err=9000,end=9050) cline
        call casetrans(cline,'lo')
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        phimlim=char2real(ifail,1)
        if(ifail.ne.0) go to 9250
        phimaccept=char2real(ifail,2)
        if(ifail.ne.0) go to 9250
        fpos=-1
        lpos=-1
        call linesplit(ifail,3)
        if(ifail.eq.0)then
          fracphim=char2real(ifail,3)
          if(ifail.ne.0)then
            atemp=cline(left_word(3):right_word(3))
            call casetrans(atemp,'lo')
            if(atemp.eq.'nomemsave')then
              memsave=0
              lpos=3
            else if(atemp.eq.'memsave')then
              memsave=1
              lpos=3
            else
              write(amessage,150) trim(aline),trim(afile)
150           format(' Third variable on line ',a,' of file ',a,' must be ',  &
              'FRACPHIM (a number) or MEMSAVE (a character variable either ', &
              '"memsave" or "nomemsave").')
              go to 9890
            end if
          else
            fpos=1
          end if
        end if
        if(fpos.ne.-1)then
          call linesplit(ifail,4)
          if(ifail.eq.0)then
            atemp=cline(left_word(4):right_word(4))
            call casetrans(atemp,'lo')
            if(atemp.eq.'nomemsave')then
              memsave=0
              lpos=4
            else if(atemp.eq.'memsave')then
              memsave=1
              lpos=4
            else
              write(amessage,160) trim(aline),trim(afile)
160           format(' MEMSAVE must be "memsave" or "nomemsave" at line ',a,  &
              ' of file ',a,'.')
              go to 9890
            end if
          end if
        end if
        if(lpos.ne.-1)then
          call linesplit(ifail,lpos+1)
          if(ifail.ne.0) go to 200
          atemp=cline(left_word(lpos+1):right_word(lpos+1))
          call casetrans(atemp,'lo')
          if(atemp.eq.'cg')then
            conjgrad=1
          else if(atemp.eq.'nocg')then
            conjgrad=0
          else
            go to 9250
          end if
          call linesplit(ifail,lpos+3)
          if(ifail.ne.0)then
            write(amessage,170) trim(aline),trim(afile)
170         format(' If CONJGRAD is assigned a value (either "cg" or "nocg") ', &
            'then values must be supplied ', &
            'for CGRTOL and CGITNLIM on line ',a,' of file ',a,'.')
            go to 9890
          end if
          cgrtol=char2real(ifail,lpos+2)
          if(ifail.ne.0) go to 9250
          cgitnlim=char2int(ifail,lpos+3)
          if(ifail.ne.0) go to 9250
        end if
200     continue
        iline=iline+1
        call num2char(iline,aline)
        read(pestunit,'(a)',err=9000,end=9250) cline
        call linesplit(ifail,3)
        if(ifail.ne.0) go to 9100
        wfinit=char2real(ifail,1)
        if(ifail.ne.0) go to 9250
        wfmin=char2real(ifail,2)
        if(ifail.ne.0) go to 9250
        wfmax=char2real(ifail,3)
        if(ifail.ne.0) go to 9250
        call linesplit(ifail,4)
        if(ifail.eq.0)then
          atemp=cline(left_word(4):right_word(4))
          call casetrans(atemp,'lo')
          if((atemp.ne.'linreg').and.(atemp.ne.'nonlinreg'))go to 9250
          if(atemp.eq.'linreg')then
            linreg=1
          else
            linreg=0
          end if
        end if
        iline=iline+1
        call num2char(iline,aline)
        read(pestunit,'(a)',err=9000,end=9250) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        wffac=char2real(ifail,1)
        if(ifail.ne.0) go to 9250
        wftol=char2real(ifail,2)
        if(ifail.ne.0) go to 9250
        call linesplit(ifail,3)
        if(ifail.eq.0)then
          iregadj=char2int(ifail,3)
          if(ifail.ne.0) go to 9250
        end if

300     continue
        rewind(unit=pestunit)

! -- The scratch file is opened.

        tempunit=nextunit()
        open(unit=tempunit,file=tempfile,action='write',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,302) trim(tempfile)
302       format(' Cannot open temporary scratch file ',a,'.')
          go to 9890
        end if

! -- Blocks of the GENREG control file are now read and processed.

! -- First the reading of this file is iniitalised.

        call addquote(genfile,afile)
        iline=0
        nblock=0
        rblock=0

! -- A new block is now searched for.

        do
          iline=iline+1
          call num2char(iline,aline)
          read(genunit,'(a)',err=9000,end=8000) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).eq.'#') cycle
          call linesplit(ifail,2)
          if(ifail.ne.0) go to 9300
          atemp=cline(left_word(1):right_word(1))
          call casetrans(atemp,'lo')
          if(atemp.ne.'start')then
            write(amessage,320) trim(aline),trim(afile)
320         format(' First entry on line ',a,' of file ',a,' expected ',  &
            'to be "START".')
            go to 9890
          end if
          atemp=cline(left_word(2):right_word(2))
          call casetrans(atemp,'lo')
          if(atemp.eq.'regspec')then
            blocktype=1
            nblock=nblock+1
          else if(atemp(1:7).eq.'regular')then
            blocktype=2
            rblock=rblock+1
            if(rblock.eq.2)then
              write(amessage,325) trim(afile)
325           format(' More than one REGULARISATION block found in file ',a,'.')
              go to 9890
            end if
          else
            write(amessage,330) trim(aline),trim(afile)
330         format(' Second entry on line ',a,' of file ',a,' must be ', &
            '"REGSPEC" or "REGULARISATION".')
            go to 9890
          end if

! -- All entries in the block are now initialised.

          spec_type=unassigned_int
          reg_type=unassigned_int
          diffrat_type=unassigned_int
          diffrat_val_type=unassigned_int
          equality_type=unassigned_int
          val_type=unassigned_int
          weight_type=unassigned_int
          weight_obs_dist=unassigned_int
          max_pilot_points=unassigned_int
          min_pilot_points=unassigned_int
          warn_less_min=unassigned_int

          value=unassigned_real
          weight=unassigned_real
          weight_multiplier=unassigned_real
          obs_dist_a=unassigned_real
          obs_dist_b=unassigned_real
          obs_dist_c=unassigned_real
          obs_dist_minwt=unassigned_real
          obs_dist_maxwt=unassigned_real
          search_radius=unassigned_real
          weight_sep_a=unassigned_real
          weight_sep_b=unassigned_real
          weight_sep_c=unassigned_real
          weight_sep_anis_bearing=unassigned_real
          weight_sep_anis_ratio=unassigned_real
          weight_sep_minwt=unassigned_real
          weight_sep_maxwt=unassigned_real

          family_prefix=' '
          family_prefix_1=' '
          family_prefix_2=' '
          reg_group=' '
          paramname=' '
          pilot_points_filename=' '
          pilot_points_filename_dr=' '
          pilot_points_filename_1=' '
          pilot_points_filename_2=' '
          obs_coord_filename=' '

          anynewprior=0

          if(blocktype.eq.1)then
            ablock='REGSPEC'
          else
            ablock='REGULARISATION'
          end if
          write(6,332) trim(ablock)
332       format(/,' Reading data for ',a,' block...')

! -- Entries are now read and stored until the end of the block is encountered.

          do
            iline=iline+1
            call num2char(iline,aline)
            read(genunit,'(a)',end=9050,err=9000) cline
            if(cline.eq.' ') cycle
            if(cline(1:1).eq.'#') cycle
            cline=adjustl(cline)
            atemp=cline(1:4)
            call casetrans(atemp,'lo')
            if(atemp(1:4).ne.'end ')write(6,333) trim(cline)
333         format(6x,a)
            call linesplit(ifail,2)
            if(ifail.ne.0) go to 9300
            atemp=cline(left_word(1):right_word(1))
            call casetrans(atemp,'lo')
            ibeg=left_word(2)
            iend=len_trim(cline)

            select case(atemp)

            case('phimlim')
              if(blocktype.ne.2) go to 9400
              if(reg_phimlim.gt.unassigned_thresh)go to 9850
              reg_phimlim=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_phimlim.le.0.0) go to 9450

            case('phimaccept')
              if(blocktype.ne.2)go to 9400
              if(reg_phimaccept.gt.unassigned_thresh)go to 9850
              reg_phimaccept=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_phimaccept.le.0.0) go to 9450

            case('fracphim')
              if(blocktype.ne.2)go to 9400
              if(reg_fracphim.gt.unassigned_thresh)go to 9850
              reg_fracphim=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('memsave')
              if(blocktype.ne.2)go to 9400
              if(reg_memsave.ne.unassigned_int)go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'memsave')then
                reg_memsave=1
              else if(atemp1.eq.'nomemsave')then
                reg_memsave=0
              else
                go to 9500
              end if

            case('conjgrad')
              if(blocktype.ne.2)go to 9400
              if(reg_conjgrad.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'cg')then
                reg_conjgrad=1
              else if(atemp1.eq.'nocg')then
                reg_conjgrad=0
              else
                go to 9500
              end if

            case('cgrtol')
              if(blocktype.ne.2)go to 9400
              if(reg_cgrtol.gt.unassigned_thresh) go to 9850
              reg_cgrtol=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_cgrtol.le.0.0) go to 9450

            case('cgitnlim')
              if(blocktype.ne.2)go to 9400
              if(reg_cgitnlim.ne.unassigned_int) go to 9850
              reg_cgitnlim=char2int(ifail,2)
              if(ifail.ne.0) go to 9370
              if(reg_cgitnlim.le.0) go to 9450

            case('wfinit')
              if(blocktype.ne.2)go to 9400
              if(reg_wfinit.gt.unassigned_thresh) go to 9850
              reg_wfinit=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_wfinit.le.0.0) go to 9450

            case('wfmin')
              if(blocktype.ne.2)go to 9400
              if(reg_wfmin.gt.unassigned_thresh) go to 9850
              reg_wfmin=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_wfmin.le.0.0) go to 9450

            case('wfmax')
              if(blocktype.ne.2)go to 9400
              if(reg_wfmax.gt.unassigned_thresh) go to 9850
              reg_wfmax=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_wfmax.le.0.0) go to 9450

            case('wffac')
              if(blocktype.ne.2)go to 9400
              if(reg_wffac.gt.unassigned_thresh) go to 9850
              reg_wffac=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_wffac.le.0.0) go to 9450

            case('wftol')
              if(blocktype.ne.2)go to 9400
              if(reg_wftol.gt.unassigned_thresh) go to 9850
              reg_wftol=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(reg_wftol.le.0.0) go to 9450

            case('linreg')
              if(blocktype.ne.2)go to 9400
              if(reg_linreg.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'linreg')then
                reg_linreg=1
              else if(atemp1.eq.'nolinreg')then
                reg_linreg=0
              else
                go to 9500
              end if

            case('iregadj')
              if(blocktype.ne.2)go to 9400
              if(reg_iregadj.ne.unassigned_int) go to 9850
              reg_iregadj=char2int(ifail,2)
              if(ifail.ne.0) go to 9370

            case('spec_type')
              if(blocktype.ne.1) go to 9550
              if(spec_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'within_family')then
                spec_type=1
              else if(atemp1.eq.'between_family')then
                spec_type=2
              else if(atemp1.eq.'indiv_family')then
                spec_type=3
              else
                go to 9500
              end if

            case('equality_type')
              if(blocktype.ne.1) go to 9550
              if(equality_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'next_pcf')then
                equality_type=1
              else if(atemp1.eq.'spatial')then
                equality_type=2
              else
                go to 9500
              end if

            case('reg_type')
              if(blocktype.ne.1) go to 9550
              if(reg_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'specified_value')then
                reg_type=1
              else if(atemp1.eq.'difference')then
                reg_type=2
              else if(atemp1.eq.'ratio')then
                reg_type=3
              else if(atemp1.eq.'equality')then
                reg_type=4
              else
                go to 9500
              end if

            case('diffrat_type')
              if(blocktype.ne.1) go to 9550
              if(diffrat_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'name')then
                diffrat_type=1
              else if(atemp1.eq.'spatial')then
                diffrat_type=2
              else
                go to 9500
              end if

            case('diffrat_val_type')
              if(blocktype.ne.1) go to 9550
              if(diffrat_val_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'uniform')then
                diffrat_val_type=1
              else if(atemp1.eq.'file')then
                diffrat_val_type=2
              else
                go to 9500
              end if

            case('weight_type')
              if(blocktype.ne.1) go to 9550
              if(weight_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'uniform')then
                weight_type=1
              else if((atemp1.eq.'sep_power').or.(atemp1.eq.'sep_pow'))then
                weight_type=2
              else if(atemp1.eq.'sep_exp')then
                weight_type=3
              else if(atemp1.eq.'sep_log')then
                weight_type=4
              else
                go to 9500
              end if

            case('weight_obs_dist')
              if(blocktype.ne.1) go to 9550
              if(weight_obs_dist.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if((atemp1.eq.'yes').or.(atemp1.eq.'y')) then
                weight_obs_dist=1
              else if((atemp1.eq.'no').or.(atemp1.eq.'n'))then
                weight_obs_dist=0
              else
                go to 9500
              end if

            case('warn_less_min')
              if(blocktype.ne.1) go to 9550
              if(warn_less_min.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if((atemp1.eq.'yes').or.(atemp1.eq.'y')) then
                warn_less_min=1
              else if((atemp1.eq.'no').or.(atemp1.eq.'n'))then
                warn_less_min=0
              else
                go to 9500
              end if

            case('val_type')
              if(blocktype.ne.1) go to 9550
              if(val_type.ne.unassigned_int) go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              if(atemp1.eq.'uniform')then
                val_type=1
              else if(atemp1.eq.'file')then
                val_type=2
              else
                go to 9500
              end if

            case('family_prefix')
              if(blocktype.ne.1) go to 9550
              if(family_prefix.ne.' ') go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              family_prefix=atemp1
              if(family_prefix.eq.' ')then
                write(amessage,334) trim(aline),trim(afile)
334             format(' Blank FAMILY_PREFIX not allowed at line ',a,' of file ',a,'.')
                go to 9890
              end if

            case('family_prefix_1')
              if(blocktype.ne.1) go to 9550
              if(family_prefix_1.ne.' ') go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              family_prefix_1=atemp1
              if(family_prefix_1.eq.' ')then
                write(amessage,334) trim(aline),trim(afile)
                go to 9890
              end if

            case('family_prefix_2')
              if(blocktype.ne.1) go to 9550
              if(family_prefix_2.ne.' ') go to 9850
              call getfile(ifail,cline,atemp1,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(atemp1,'lo')
              family_prefix_2=atemp1
              if(family_prefix_2.eq.' ')then
                write(amessage,334) trim(aline),trim(afile)
                go to 9890
              end if

            case('reg_group')
              if(blocktype.ne.1) go to 9550
              if(reg_group.ne.' ') go to 9850
              call getfile(ifail,cline,reg_group,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(reg_group,'lo')

            case('parameter')
              if(blocktype.ne.1) go to 9550
              if(paramname.ne.' ') go to 9850
              call getfile(ifail,cline,paramname,ibeg,iend)
              if(ifail.ne.0) go to 9660
              call casetrans(paramname,'lo')

            case('pilot_points_filename')
              if(blocktype.ne.1) go to 9550
              if(pilot_points_filename.ne.' ') go to 9850
              call getfile(ifail,cline,pilot_points_filename,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,350) trim(aline),trim(afile)
350             format(' Error reading pilot points filename from line ',a,  &
                ' of file ',a,'.')
                go to 9890
              end if
              call casetrans(pilot_points_filename,'lo')       ! not in unix

            case('pilot_points_filename_dr')
              if(blocktype.ne.1) go to 9550
              if(pilot_points_filename_dr.ne.' ') go to 9850
              call getfile(ifail,cline,pilot_points_filename_dr,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,350) trim(aline),trim(afile)
                go to 9890
              end if
              call casetrans(pilot_points_filename_dr,'lo')       ! not in unix

            case('pilot_points_filename_1')
              if(blocktype.ne.1) go to 9550
              if(pilot_points_filename_1.ne.' ') go to 9850
              call getfile(ifail,cline,pilot_points_filename_1,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,350) trim(aline),trim(afile)
                go to 9890
              end if
              call casetrans(pilot_points_filename_1,'lo')       ! not in unix

            case('pilot_points_filename_2')
              if(blocktype.ne.1) go to 9550
              if(pilot_points_filename_2.ne.' ') go to 9850
              call getfile(ifail,cline,pilot_points_filename_2,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,350) trim(aline),trim(afile)
                go to 9890
              end if
              call casetrans(pilot_points_filename_2,'lo')       ! not in unix

            case('obs_coord_filename')
              if(blocktype.ne.1) go to 9550
              if(obs_coord_filename.ne.' ') go to 9850
              call getfile(ifail,cline,obs_coord_filename,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,351) trim(aline),trim(afile)
351             format(' Error reading observation coordinates filename from ', &
                'line ',a,' of file ',a,'.')
                go to 9890
              end if
              call casetrans(obs_coord_filename,'lo')       ! not in unix

            case('min_pilot_points')
              if(blocktype.ne.1) go to 9550
              if(min_pilot_points.ne.unassigned_int) go to 9850
              min_pilot_points=char2int(ifail,2)
              if(ifail.ne.0) go to 9370
              if(min_pilot_points.le.0) go to 9450

            case('max_pilot_points')
              if(blocktype.ne.1) go to 9550
              if(max_pilot_points.ne.unassigned_int) go to 9850
              max_pilot_points=char2int(ifail,2)
              if(ifail.ne.0) go to 9370
              if(max_pilot_points.le.0) go to 9450

            case('value')
              if(blocktype.ne.1) go to 9550
              if(value.gt.unassigned_thresh) go to 9850
              value=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('weight')
              if(blocktype.ne.1) go to 9550
              if(weight.gt.unassigned_thresh) go to 9850
              weight=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(weight.lt.0.0) go to 9750

            case('weight_multiplier')
              if(blocktype.ne.1) go to 9550
              if(weight_multiplier.gt.unassigned_thresh) go to 9850
              weight_multiplier=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(weight_multiplier.lt.0.0) go to 9750

            case('search_radius')
              if(blocktype.ne.1) go to 9550
              if(search_radius.gt.unassigned_real) go to 9850
              search_radius=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(search_radius.le.0.0) go to 9450

            case('obs_dist_a')
              if(blocktype.ne.1) go to 9550
              if(obs_dist_a.gt.unassigned_thresh) go to 9850
              obs_dist_a=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('obs_dist_b')
              if(blocktype.ne.1) go to 9550
              if(obs_dist_b.gt.unassigned_thresh) go to 9850
              obs_dist_b=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('obs_dist_c')
              if(blocktype.ne.1) go to 9550
              if(obs_dist_c.gt.unassigned_thresh) go to 9850
              obs_dist_c=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('obs_dist_minwt')
              if(blocktype.ne.1) go to 9550
              if(obs_dist_minwt.gt.unassigned_thresh) go to 9850
              obs_dist_minwt=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(obs_dist_minwt.lt.0.0) go to 9750

            case('obs_dist_maxwt')
              if(blocktype.ne.1) go to 9550
              if(obs_dist_maxwt.gt.unassigned_thresh) go to 9850
              obs_dist_maxwt=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(obs_dist_maxwt.lt.0.0) go to 9750

            case('weight_sep_a')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_a.gt.unassigned_thresh) go to 9850
              weight_sep_a=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('weight_sep_b')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_b.gt.unassigned_thresh) go to 9850
              weight_sep_b=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('weight_sep_c')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_c.gt.unassigned_thresh) go to 9850
              weight_sep_c=char2real(ifail,2)
              if(ifail.ne.0) go to 9350

            case('weight_sep_anis_bearing')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_anis_bearing.gt.unassigned_thresh) go to 9850
              weight_sep_anis_bearing=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if((weight_sep_anis_bearing.lt.0.0).or.          &
                 (weight_sep_anis_bearing.gt.180.0))then
                 write(amessage,352)
352              format(' WEIGHT_SEP_ANIS_BEARING must be between 0.0 and 180.0 ', &
                 'degrees inclusive.')
                 go to 9890
              end if
              if(weight_sep_anis_bearing.eq.360.0) weight_sep_anis_bearing=0.0

            case('weight_sep_anis_ratio')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_anis_ratio.gt.unassigned_thresh) go to 9850
              weight_sep_anis_ratio=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(weight_sep_anis_ratio.lt.1.0)then
                write(amessage,353)
353             format(' WEIGHT_SEP_ANIS_RATIO must not be less than 1.0')
                go to 9890
              end if

            case('weight_sep_minwt')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_minwt.gt.unassigned_thresh) go to 9850
              weight_sep_minwt=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(weight_sep_minwt.lt.0.0) go to 9750

            case('weight_sep_maxwt')
              if(blocktype.ne.1) go to 9550
              if(weight_sep_maxwt.gt.unassigned_thresh) go to 9850
              weight_sep_maxwt=char2real(ifail,2)
              if(ifail.ne.0) go to 9350
              if(weight_sep_maxwt.lt.0.0) go to 9750

            case('end')
              go to 400

            case default
              call casetrans(atemp,'hi')
              write(amessage,335) trim(atemp),trim(aline),trim(afile)
335           format(' Unrecognised keyword - "',a,'" - at line ',a,' of file ',a,'.')
              go to 9890

            end select

          end do

! -- All items in the block have been read. Now the block is processed.

400       continue
          if(blocktype.eq.1)then
            ablock='REGSPEC'
          else
            ablock='REGULARISATION'
          end if
          write(6,336) trim(ablock)
336       format(' Processing ',a,' block...')

! -- The REGULARISATION block is processed.

          if(blocktype.eq.2)then
            if(reg_phimlim.gt.unassigned_thresh)then
              phimlim=reg_phimlim
            end if
            if(reg_phimaccept.gt.unassigned_thresh)then
              phimaccept=reg_phimaccept
            end if
            if(reg_phimlim.gt.unassigned_thresh)then
              if(reg_phimaccept.lt.unassigned_thresh)then
                phimaccept=phimlim*1.05
              end if
            end if

            if(phimlim.lt.unassigned_thresh)then
              write(amessage,405) trim(afile)
405           format(' A value has not been assigned to the target measurement ', &
              'objective function PHIMLIM either in the existing PEST control ', &
              'file or in a REGULARISATION block of file ',a,'.')
              go to 9890
            end if
            if(reg_fracphim.gt.unassigned_thresh)then
              fracphim=reg_fracphim
            end if
            if(reg_memsave.ne.unassigned_int)then
              memsave=reg_memsave
            end if
            if(reg_conjgrad.ne.unassigned_int)then
              conjgrad=reg_conjgrad
            end if
            if(reg_cgrtol.gt.unassigned_thresh)then
              cgrtol=reg_cgrtol
            end if
            if(reg_cgitnlim.ne.unassigned_int)then
              cgitnlim=reg_cgitnlim
            end if
            if(reg_wfinit.gt.unassigned_thresh)then
              wfinit=reg_wfinit
            end if
            if(reg_wfmin.gt.unassigned_thresh)then
              wfmin=reg_wfmin
            end if
            if(reg_wfmax.gt.unassigned_thresh)then
              wfmax=reg_wfmax
            end if
            if(reg_wffac.gt.unassigned_thresh)then
              wffac=reg_wffac
            end if
            if(reg_wftol.gt.unassigned_thresh)then
              wftol=reg_wftol
            end if
            if(reg_linreg.ne.unassigned_int)then
              linreg=reg_linreg
            end if
            if(reg_iregadj.ne.unassigned_int)then
              iregadj=reg_iregadj
            end if

! -- Now the REGSPEC block

          else if (blocktype.eq.1)then

            if(spec_type.eq.unassigned_int)then
              write(amessage,415)
415           format(' SPEC_TYPE keyword missing from REGSPEC block.')
              go to 9890
            end if
            if(weight_type.eq.unassigned_int)then
              write(amessage,414)
414           format(' WEIGHT_TYPE keyword missing from REGSPEC block.')
              go to 9890
            end if
            if(spec_type.eq.2)then
              if((weight_type.eq.2).or.(weight_type.eq.3).or.(weight_type.eq.4))then
                if(pilot_points_filename.eq.' ')then
                  write(amessage,412)
412               format(' If WEIGHT_TYPE is "sep_power", "sep_exp" or "sep_log" then ', &
                  'a PILOT_POINTS_FILENAME keyword must be supplied in REGSPEC block.')
                  go to 9890
                end if
              end if
            end if
            if(weight_obs_dist.eq.1)then
              if((spec_type.ne.2).and.(spec_type.ne.3))then   ! This error condition handled later.
                if((obs_coord_filename.eq.' ').or.(pilot_points_filename.eq.' '))then
                  write(amessage,418)
418               format(' If WEIGHT_OBS_DIST is "yes" and SPEC_TYPE is "within_family" then ', &
                  'both an OBS_COORD_FILENAME and a PILOT_POINTS_FILENAME must be ', &
                  'supplied in REGSPEC block.')
                  go to 9890
                end if
              end if
            end if
            if(weight_type.eq.1)then
              if(weight.lt.unassigned_thresh)then
                write(amessage,421)
421             format(' If WEIGHT_TYPE is "uniform" then a WEIGHT keyword must be ', &
                'supplied in REGSPEC block.')
                go to 9890
              end if
            end if
            if(weight_obs_dist.eq.1)then
              if((obs_dist_a.lt.unassigned_thresh).or.   &
                 (obs_dist_b.lt.unassigned_thresh).or.   &
                 (obs_dist_c.lt.unassigned_thresh).or.   &
                 (obs_dist_minwt.lt.unassigned_thresh).or. &
                 (obs_dist_maxwt.lt.unassigned_thresh))then
                 write(amessage,435)
435              format(' If WEIGHT_OBS_DIST is "yes"  then the following ', &
                 'keywords must be supplied in REGSPEC block - OBS_DIST_A ', &
                 'OBS_DIST_B, OBS_DIST_C, OBS_DIST_MINWT, OBS_DIST_MAXWT.')
                 go to 9890
              end if
              if(obs_dist_minwt.ge.obs_dist_maxwt)then
                write(amessage,436)
436             format(' OBS_DIST_MAXWT must be greater than OBD_DIST_MINWT ', &
                'in REGSPEC block.')
                go to 9890
              end if
            end if
            if(reg_group.eq.' ')then
              write(amessage,390)
390           format(' REG_GROUP keyword missing from REGSPEC block.')
              go to 9890
            end if

            if(weight_obs_dist.eq.unassigned_int) weight_obs_dist=0

! -- WITHIN_FAMILY regularisation

            if(spec_type.eq.1)then

              if(reg_type.eq.unassigned_int)then
                write(amessage,416)
416             format(' If SPEC_TYPE is "within_family" then a REG_TYPE keyword ', &
                'must be supplied in REGSPEC block.')
                go to 9890
              end if
              if(family_prefix.eq.' ')then
                write(amessage,420)
420             format(' If SPEC_TYPE is "within_family" then a FAMILY_PREFIX ', &
                'keyword must be supplied in REGSPEC block.')
                go to 9890
              end if
              nf=len_trim(family_prefix)
              nrpar=0
              do ipar=1,npar
                if(apar(ipar)(1:nf).eq.family_prefix(1:nf))then
                  if(itrans(ipar).ge.0)then
                    nrpar=nrpar+1
                    if(nrpar.gt.MAXREGPAR)then
                      write(amessage,408)
408                   format(' Out of memory - increase MAXREGPAR and re-compile ', &
                      'program.')
                      go to 9890
                    end if
                    arpar(nrpar)=apar(ipar)
                    parindex(nrpar)=ipar
                  end if
                end if
              end do
              if(nrpar.eq.0)then
                call addquote(pestfile,bfile)
                write(amessage,419) trim(bfile),trim(family_prefix)
419             format(' No adjustable (i.e. non-tied and nonfixed) parameters ', &
                'in PEST control file ',a,   &
                ' have a prefix of "',a,'", this being the FAMILY_PREFIX ', &
                'supplied in REGSPEC block.')
                go to 9890
              end if
              if((reg_type.ne.1).and.(reg_type.ne.4))then
                write(amessage,424)
424             format(' If SPEC_TYPE is "within_family" then REG_TYPE must ', &
                'be either "specified_value" or "equality".')
                go to 9890
              end if

! -- Specified value regularisation

              if(reg_type.eq.1)then

                if(val_type.eq.unassigned_int)then
                  write(amessage,417)
417               format(' If REG_TYPE is "specified_value" and SPEC_TYPE is ',  &
                  '"within_family", then a VAL_TYPE keyword must be supplied ',  &
                  'in REGSPEC block.')
                  go to 9890
                end if

                if(weight_type.ne.1)then
                  write(amessage,430)
430               format(' If REG_TYPE is "specified_value" then WEIGHT_TYPE must ', &
                  '"uniform" in REGSPEC block.')
                  go to 9890
                end if

                if(val_type.eq.2)then
                  if(pilot_points_filename.eq.' ')then
                    write(amessage,438)
438                 format(' If VAL_TYPE is "file" then a PILOT_POINTS_FILENAME ', &
                    'must be supplied in REGSPEC block.')
                    go to 9890
                  end if
                end if

                if((val_type.eq.2).or.(weight_obs_dist.ne.0))then
                  call read_pp_file(ifail,nrpar,arpar,east_par,north_par,  &
                  val_par,parindex,family_prefix,pilot_points_filename)
                  if(ifail.ne.0) go to 9890
                end if

                if(weight_obs_dist.ne.0)then
                  call read_oc_file(ifail,MAXOBS,mobs,east_obs,north_obs,  &
                  obs_coord_filename)
                  if(ifail.ne.0) go to 9890
                end if

                if(val_type.eq.1)then
                  if(value.lt.unassigned_thresh)then
                    write(amessage,439)
439                 format(' When SPEC_TYPE is "within_family", REG_TYPE is ', &
                    '"specified_value" and VAL_TYPE is "uniform", the VALUE keyword must ', &
                    'be supplied.')
                    go to 9890
                  end if
                  do irpar=1,nrpar
                    val_par(irpar)=value
                  end do
                end if

                do irpar=1,nrpar
                  iprior=iprior+1
                  anynewprior=anynewprior+1
                  call num2char(iprior,aprior)
                  priorname='gr'//trim(aprior)
                  wt=weight
                  if(weight_obs_dist.eq.1)then
                    distmin=1.0d300
                    do iobs=1,mobs
                      dist=(east_par(irpar)-east_obs(iobs))*    &
                           (east_par(irpar)-east_obs(iobs))+    &
                           (north_par(irpar)-north_obs(iobs))*  &
                           (north_par(irpar)-north_obs(iobs))
                      if(dist.lt.distmin)distmin=dist
                    end do
                    distmin=sqrt(distmin)
                    ww=obs_dist_a+obs_dist_b*(distmin**obs_dist_c)
                    if(ww.gt.obs_dist_maxwt)ww=obs_dist_maxwt
                    if(ww.lt.obs_dist_minwt)ww=obs_dist_minwt
                    wt=ww*wt
                  end if
                  ww=wt
                  if(weight_multiplier.gt.unassigned_thresh)then
                    ww=ww*weight_multiplier
                  end if
                  if(itrans(parindex(irpar)).eq.0)then
                    write(tempunit,440) trim(priorname),trim(arpar(irpar)), &
                    val_par(irpar),ww,trim(reg_group)
440                 format(1x,a,t15,'1.0 * ',a,' = ',1pg14.7,1x,1pg13.6,1x,a)
                  else
                    if(val_par(irpar).le.0.0)then
                      write(amessage,450) trim(arpar(irpar))
450                   format(' Parameter "',a,'" is log transformed in PEST ',   &
                      'control file; specified value for regularisation prior ', &
                      'information must be positive for this parameter.')
                      go to 9890
                    end if
                    write(tempunit,460) trim(priorname),trim(arpar(irpar)), &
                    log10(val_par(irpar)),ww,trim(reg_group)
460                 format(1x,a,t15,'1.0 * log(',a,') = ',1pg14.7,1x, &
                    1pg13.6,1x,a)
                  end if
                end do

! -- Equality regularisation.

              else if(reg_type.eq.4)then

                 if(equality_type.eq.unassigned_int)then
                   write(amessage,610)
610                format(' IF REG_TYPE is specified as "equality" then an ', &
                   'EQUALITY_TYPE keyword must be present in REGSPEC block.')
                   go to 9890
                 end if
                 call log_test(ilog,npar,nrpar,parindex,itrans)
                 if(ilog.lt.0)then
                   write(amessage,630)
630                format(' If REG_TYPE is set to "equality" then all non-tied ', &
                   'and non-fixed parameters ', &
                   'in selected family must be either log transformed or ', &
                   'untransformed.')
                   go to 9890
                 end if
                 if(nrpar.eq.1)then
                   write(amessage,614)
614                format(' If REG_TYPE is set to "equality", there must be more ', &
                   'than one adjustable parameter in specified family.')
                   go to 9890
                 end if
                 if(warn_less_min.eq.1)then
                   write(amessage,611)
611                format(' WARN_LESS_MIN must be set to "no" or omitted if ', &
                   'SPEC_TYPE is set to "within_family" in REGSPEC block.')
                   go to 9890
                 end if
                 if(equality_type.eq.1)then
                   if(weight_type.ne.1)then
                     write(amessage,615)
615                  format(' If EQUALITY_TYPE is specified as "next_pcf" then ', &
                     'WEIGHT_TYPE must be set to "uniform" in REGSPEC block.')
                     go to 9890
                   end if
                   if(weight.lt.unassigned_thresh)then
                     write(amessage,620)
620                  format(' If WEIGHT_TYPE is set to "uniform" then a WEIGHT ', &
                     'keyword must be supplied in REGSPEC block.')
                     go to 9890
                   end if
                   if(weight_obs_dist.ne.0)then
                     write(amessage,621)
621                  format(' If EQUALITY_TYPE is specified as "next_pcf" then ', &
                     'WEIGHT_OBS_DIST must be set to "no" (or ommitted) in REGSPEC block.')
                     go to 9890
                   end if
                   ww=weight
                   if(weight_multiplier.gt.unassigned_thresh)then
                     ww=ww*weight_multiplier
                   end if
                   do irpar=1,nrpar
                     aapar1=arpar(irpar)
                     if(irpar.lt.nrpar)then
                       aapar2=arpar(irpar+1)
                     else
                       aapar2=arpar(1)
                     end if
                     iprior=iprior+1
                     anynewprior=anynewprior+1
                     call num2char(iprior,aprior)
                     priorname='gr'//trim(aprior)
                     if(ilog.eq.0)then
                       write(tempunit,640) trim(priorname),trim(aapar1), &
                       trim(aapar2),ww,trim(reg_group)
640                    format(1x,a,t15,'1.0 * ',a,' - 1.0 * ',a,   &
                       ' = 0.0 ',1pg13.6,1x,a)
                     else
                       write(tempunit,650) trim(priorname),trim(aapar1), &
                       trim(aapar2),ww,trim(reg_group)
650                    format(1x,a,t15,'1.0 * log(',a,') - 1.0 * log(',a,   &
                       ') = 0.0 ',1pg13.6,1x,a)
                     end if
                   end do
                 else if(equality_type.eq.2)then
                   if(pilot_points_filename.eq.' ')then
                     write(amessage,660)
660                  format(' If EQUALITY_TYPE is "spatial" then a ', &
                     'PILOT_POINTS_FILENAME must be supplied in REGSPEC block.')
                     go to 9890
                   end if
                   if((min_pilot_points.eq.unassigned_int).or.   &
                      (max_pilot_points.eq.unassigned_int).or.   &
                      (search_radius.lt.unassigned_thresh))then
                      write(amessage,670)
670                   format(' If EQUALITY_TYPE is "spatial" then all ', &
                      'of the MIN_PILOT_POINTS, MAX_PILOT_POINTS and ', &
                      'SEARCH_RADIUS keywords must be supplied in REGSPEC block.')
                      go to 9890
                   end if
                   if(min_pilot_points.le.0)then
                     write(amessage,680)
680                  format(' MIN_PILOT_POINTS must be greater than zero in ', &
                     'REGSPEC block.')
                     go to 9890
                   end if
                   if(max_pilot_points.lt.min_pilot_points)then
                     write(amessage,681)
681                  format(' MAX_PILOT_POINTS must not be less than MIN_PILOT_POINTS ', &
                     'in REGSPEC block.')
                     go to 9890
                   end if
                   if(search_radius.le.0.0d0)then
                     write(amessage,690)
690                  format(' SEARCH_RADIUS must be greater than zero in REGSPEC block.')
                     go to 9890
                   end if

                   if((weight_type.eq.2).or.(weight_type.eq.3).or.(weight_type.eq.4))then
                     if((weight_sep_a.lt.unassigned_thresh).or.     &
                        (weight_sep_b.lt.unassigned_thresh).or.     &
                        (weight_sep_c.lt.unassigned_thresh).or.     &
                        (weight_sep_anis_bearing.lt.unassigned_thresh).or. &
                        (weight_sep_anis_ratio.lt.unassigned_thresh).or.   &
                        (weight_sep_minwt.lt.unassigned_thresh).or.        &
                        (weight_sep_maxwt.lt.unassigned_thresh)) then
                        write(amessage,700)
700                     format(' If WEIGHT_TYPE is "sep_power", "sep_exp" or "sep_log" ', &
                        'then the WEIGHT_SEP_A, WEIGHT_SEP_B, WEIGHT_SEP_C ', &
                        'WEIGHT_SEP_MINWT, WEIGHT_SEP_MAXWT, ',                &
                        'WEIGHT_SEP_ANIS_BEARING and WEIGHT_SEP_ANIS_RATIO ', &
                        'keywords must all be present in REGSPEC block.')
                        go to 9890
                     end if
                     if(weight_sep_maxwt.lt.weight_sep_minwt)then
                       write(amessage,702)
702                    format(' If WEIGHT_TYPE is "sep_power", "sep_exp" or "sep_log" ', &
                       'then WEIGHT_SEP_MAXWT must not be smaller than WEIGHT_SEP_MINWT.')
                       go to 9890
                     end if
                   end if
                   if((weight_type.eq.3).and.(weight_sep_c.lt.0.0))then
                     write(amessage,701)
701                  format(' If WEIGHT_TYPE is "sep_exp" then WEIGHT_SEP_C must ', &
                     'be positive.')
                     go to 9890
                   end if
                   call read_pp_file(ifail,nrpar,arpar,east_par,north_par,  &
                   val_par,parindex,family_prefix,pilot_points_filename)
                   if(ifail.ne.0) go to 9890
                   if(weight_obs_dist.eq.1)then
                     call read_oc_file(ifail,MAXOBS,mobs,east_obs,north_obs,  &
                     obs_coord_filename)
                     if(ifail.ne.0) go to 9890
                   end if
                   if(weight_sep_anis_ratio.lt.unassigned_thresh)then
                     weight_sep_anis_ratio=1.0d0
                     weight_sep_anis_bearing=0.0d0
                   end if
                   anis_angle=90.0d0-weight_sep_anis_bearing
                   cosang=cos(anis_angle*3.14159265/180.0)
                   sinang=sin(anis_angle*3.14159265/180.0)
                   do irpar=1,nrpar
                     aapar1=arpar(irpar)
                     iparcount=0
                     do jrpar=1,nrpar
                       if(irpar.eq.jrpar)cycle
                       sepx=east_par(jrpar)-east_par(irpar)
                       sepy=north_par(jrpar)-north_par(irpar)
                       if(weight_sep_anis_ratio.ne.1.0d0)then
                         asepx=sepx*cosang+sepy*sinang
                         asepy=-sepx*sinang+sepy*cosang
                         asepy=asepy*weight_sep_anis_ratio
                       else
                         asepx=sepx
                         asepy=sepy
                       end if
                       ddist(jrpar)=sqrt(asepx*asepx+asepy*asepy)
                     end do
                     do i=1,max_pilot_points
                       mindist=1.0d300
                       jmindist=0
                       do jrpar=1,nrpar
                         if(irpar.eq.jrpar) cycle
                         if(ddist(jrpar).gt.search_radius)cycle
                         if(ddist(jrpar).lt.mindist)then
                           mindist=ddist(jrpar)
                           jmindist=jrpar
                         end if
                       end do
                       if(jmindist.eq.0) go to 750
                       aapar2=arpar(jmindist)
                       ddist(jmindist)=1.0d301
                       if(weight_type.eq.1)then
                         ww=weight
                       else
                         if(weight_type.eq.2)then
                           ww=weight_sep_a+weight_sep_b*mindist**weight_sep_c
                         else if(weight_type.eq.3)then
                           ww=weight_sep_a+weight_sep_b*exp(-mindist*weight_sep_c)
                         else if(weight_type.eq.4)then
                           if(mindist.le.0.0d0)then
                             write(amessage,710)
710                          format(' WEIGHT_TYPE must not be "sep_log" if two parameters ', &
                             'are at the same location.')
                             go to 9890
                           end if
                           ww=weight_sep_a+weight_sep_b*((log10(mindist))**weight_sep_c)
                         end if
                         if(ww.gt.weight_sep_maxwt)ww=weight_sep_maxwt
                         if(ww.lt.weight_sep_minwt)ww=weight_sep_minwt
                       end if
                       if(irpar.lt.jmindist)then
                         if(weight_obs_dist.eq.1)then
                           do ii=1,2
                             if(ii.eq.1)then
                               eee=east_par(irpar)
                               nnn=north_par(irpar)
                             else
                               eee=east_par(jmindist)
                               nnn=north_par(jmindist)
                             end if
                             smalldist=1.0d300
                             do iobs=1,mobs
                               dist=(eee-east_obs(iobs))*    &
                                    (eee-east_obs(iobs))+    &
                                    (nnn-north_obs(iobs))*  &
                                    (nnn-north_obs(iobs))
                               if(dist.lt.smalldist)smalldist=dist
                              end do
                             if(ii.eq.1)then
                               distmin1=sqrt(smalldist)
                             else
                               distmin2=sqrt(smalldist)
                             end if
                           end do
                           distmin=max(distmin1,distmin2)
                           wtd=obs_dist_a+obs_dist_b*(distmin**obs_dist_c)
                           if(wtd.gt.obs_dist_maxwt)wtd=obs_dist_maxwt
                           if(wtd.lt.obs_dist_minwt)wtd=obs_dist_minwt
                           ww=ww*wtd
                         end if
                       end if
                       if(weight_multiplier.gt.unassigned_thresh)ww=ww*weight_multiplier
                       iparcount=iparcount+1
                       if(irpar.lt.jmindist)then
                         iprior=iprior+1
                         anynewprior=anynewprior+1
                         call num2char(iprior,aprior)
                         priorname='gr'//trim(aprior)
                         if(ilog.eq.0)then
                           write(tempunit,640) trim(priorname),trim(aapar1), &
                           trim(aapar2),ww,trim(reg_group)
                         else
                           write(tempunit,650) trim(priorname),trim(aapar1), &
                           trim(aapar2),ww,trim(reg_group)
                         end if
                       end if
                     end do
750                  continue
                     if(iparcount.lt.min_pilot_points)then
                       if(weight_sep_anis_ratio.eq.1.0d0)then
                         write(amessage,771) trim(arpar(irpar))
771                      format(' Less than MIN_PILOT_POINTS parameters lie ', &
                         'within one SEARCH_RADIUS of parameter "',a,'".')
                         go to 9890
                       else
                         write(amessage,772) trim(arpar(irpar))
772                      format(' Less than MIN_PILOT_POINTS parameters lie ', &
                         'within one anisotropy-adjusted search radius of ', &
                         'parameter "',a,'".')
                         go to 9890
                       end if
                     end if
                   end do
                 end if         ! end of equality_type=2
              end if           ! end of reg_type=4

            else if (spec_type.eq.2)then

! -- Between-family regularisation.

              if((family_prefix_1.eq.' ').or.(family_prefix_2.eq.' '))then
                write(amessage,805)
805             format(' If SPEC_TYPE is "between_family" then both a ', &
                '"FAMILY_PREFIX_1" and a "FAMILY_PREFIX_2 keyword must ', &
                'be supplied in REGSPEC block.')
                go to 9890
              end if
              if(family_prefix_1.eq.family_prefix_2)then
                write(amessage,806)
806             format(' FAMILY_PREFIX_1 must be different from FAMILY_PREFIX_2 ', &
                'in REGSPEC block.')
                go to 9890
              end if
              if((reg_type.ne.2).and.(reg_type.ne.3))then
                write(amessage,810)
810             format(' If SPEC_TYPE is "between_family" then REG_TYPE must ', &
                'be "difference" or "ratio" in REGSPEC block.')
                go to 9890
              end if
              if(diffrat_type.eq.unassigned_int)then
                write(amessage,820)
820             format(' If SPEC_TYPE is "between_family" then a DIFFRAT_TYPE ', &
                'keyword must be supplied in REGSPEC block.')
                go to 9890
              end if
              if(diffrat_val_type.eq.unassigned_int)then
                write(amessage,830)
830             format(' If SPEC_TYPE is "between_family" then a ',  &
                'DIFFRAT_VAL_TYPE keyword must be supplied in REGSPEC block.')
                go to 9890
              end if
              if(diffrat_type.eq.2)then
                if(diffrat_val_type.eq.2)then
                  write(amessage,840)
840               format(' If DIFFRAT_TYPE is "spatial" then DIFFRAT_VAL_TYPE ', &
                  'must be "uniform" in REGSPEC block.')
                  go to 9890
                end if
              end if
              if(diffrat_type.eq.2)then
                if((pilot_points_filename_1.eq.' ').or.     &
                   (pilot_points_filename_2.eq.' ').or.     &
                   (search_radius.lt.unassigned_thresh).or.  &
                   (min_pilot_points.eq.unassigned_int).or.  &
                   (max_pilot_points.eq.unassigned_int))then
                   write(amessage,850)
850                format(' If DIFFRAT_TYPE is "spatial" then all of the ', &
                   'following keywords must be supplied in the REGSPEC block - ', &
                   'PILOT_POINTS_FILENAME_1, PILOT_POINTS_FILENAME_2, '  &
                   'SEARCH_RADIUS, MIN_PILOT_POINTS and MAX_PILOT_POINTS.')
                   go to 9890
                end if
              end if
              if(diffrat_val_type.eq.1)then
                if(value.lt.unassigned_thresh)then
                  write(amessage,860)
860               format(' If DIFFRAT_VAL_TYPE is "uniform", then a VALUE keyword ', &
                  'must be provided in REGSPEC block.')
                  go to 9890
                end if
              else if(diffrat_val_type.eq.2)then
                if(pilot_points_filename_dr.eq.' ')then
                  write(amessage,870)
870               format(' If DIFFRAT_VAL_TYPE is "file" then a ',  &
                  'PILOT_POINTS_FILENAME_DR keyword must be provided in REGSPEC block.')
                  go to 9890
                end if
              end if

              if(weight_type.ne.1)then
                write(amessage,871)
871             format(' If SPEC_TYPE is "between_family" then WEIGHT_TYPE ', &
                'must be set to "uniform" in REGSPEC block.')
                go to 9890
              end if
              if(weight.lt.unassigned_thresh)then
                write(amessage,872)
872             format(' If WEIGHT_TYPE is set to "uniform" then a WEIGHT ', &
                'keyword must be present in REGSPEC block.')
                go to 9890
              end if

              call get_regularised_parameters(ifail,npar,MAXREGPAR,nrpar_1,  &
              family_prefix_1,itrans,parindex_1,apar,arpar_1,pestfile)
              if(ifail.ne.0) go to 9890
              call get_regularised_parameters(ifail,npar,MAXREGPAR,nrpar_2,  &
              family_prefix_2,itrans,parindex_2,apar,arpar_2,pestfile)
              if(ifail.ne.0) go to 9890

              call log_test(ilog_1,npar,nrpar_1,parindex_1,itrans)
              if(ilog_1.lt.0)then
                write(amessage,875)
875             format(' If SPEC_TYPE is "between family" then non-tied ', &
                'and non-fixed members of both families must all be log-transformed ', &
                'or must all be untransformed.')
                go to 9890
              end if
              call log_test(ilog_2,npar,nrpar_2,parindex_2,itrans)
              if(ilog_2.lt.0)then
                write(amessage,875)
                go to 9890
              end if
              if(ilog_1.ne.ilog_2)then
                write(amessage,875)
                go to 9890
              end if
              if(ilog_1.eq.1)then
                if(reg_type.eq.2)then
                  write(amessage,876)
876               format(' If REG_TYPE is set to "difference" then parameters ', &
                  'in both families must be untransformed - use "ratio" instead.')
                  go to 9890
                end if
              else if(ilog_1.eq.0)then
                if(reg_type.eq.3)then
                  write(amessage,877)
877               format(' If REG_TYPE is set to "ratio" then parameters ', &
                  'in both families must be log-transformed - use "difference" ', &
                  'instead.')
                  go to 9890
                end if
              end if

              if(weight_obs_dist.ne.0)then
                if(obs_coord_filename.eq.' ')then
                  write(amessage,878)
878               format(' If WEIGHT_OBS_DIST is "yes" then an OBS_COORD_FILENAME ', &
                  'keyword must be provided in REGSPEC block.')
                  go to 9890
                end if
                call read_oc_file(ifail,MAXOBS,mobs,east_obs,north_obs,  &
                obs_coord_filename)
                if(ifail.ne.0) go to 9890
              end if

              if(diffrat_type.eq.1)then

                if(nrpar_1.ne.nrpar_2)then
                  write(amessage,880)
880               format(' If DIFFRAT_TYPE is "name" then the number of adjustable ', &
                  '(i.e. non-fixed and non-tied) parameters and the roots of their ', &
                  'names must be the same for the two parameter families.')
                  go to 9890
                end if

                if(weight_obs_dist.eq.1)then
                  if(pilot_points_filename_dr.eq.' ')then
                    write(amessage,881)
881                 format(' DIFFRAT_TYPE is "name" and WEIGHT_OBS_DIST is "yes" ', &
                    'then a PILOT_POINTS_FILENAME_DR keyword must be present in ', &
                    'REGSPEC block.')
                    go to 9890
                  end if
                end if

                if((diffrat_val_type.eq.2).or.(weight_obs_dist.eq.1))then
                   call read_pp_file(ifail,nrpar_1,arpar_1,east_par,north_par,  &
                   val_par,parindex_1,family_prefix_1,pilot_points_filename_dr)
                   if(ifail.ne.0) go to 9890
                end if

                nf_1=len_trim(family_prefix_1)
                nf_2=len_trim(family_prefix_2)
                do irpar=1,nrpar_1
                  do jrpar=1,nrpar_2
                    if(arpar_2(jrpar)(nf_2:).eq.arpar_1(irpar)(nf_1:)) go to 900
                  end do
                  write(amessage,890) trim(arpar_1(irpar))
890               format(' No match for parameter "',a,'" from first parameter ', &
                  'family was found in second parameter family.')
                  go to 9890
900               continue
                  if(diffrat_val_type.eq.1)then
                    diffrat=value
                  else
                    diffrat=val_par(irpar)
                  end if
                  if(ilog_1.eq.1)then
                    if(diffrat.le.0.0)then
                      if(diffrat_val_type.eq.1)then
                        write(amessage,910)
910                     format(' VALUE must be greater than zero if REG_TYPE ', &
                        'is set to "ratio".')
                        go to 9890
                      else
                        write(amessage,920)
920                     format(' If REG_SEC is set to "ratio" then ratio values ', &
                        'supplied in the PILOT_POINTS_FILENAME_DR file must be ', &
                        'greater than zero.')
                        go to 9890
                      end if
                    end if
                  end if
                  ilog=ilog_1                            ! for safety
                  if(ilog_1.eq.1) diffrat=log10(diffrat)
                  wt=weight
                  if(weight_obs_dist.eq.1)then
                    distmin=1.0d300
                    do iobs=1,mobs
                      dist=(east_par(irpar)-east_obs(iobs))*    &
                           (east_par(irpar)-east_obs(iobs))+    &
                           (north_par(irpar)-north_obs(iobs))*  &
                           (north_par(irpar)-north_obs(iobs))
                      if(dist.lt.distmin)distmin=dist
                    end do
                    distmin=sqrt(distmin)
                    ww=obs_dist_a+obs_dist_b*(distmin**obs_dist_c)
                    if(ww.gt.obs_dist_maxwt)ww=obs_dist_maxwt
                    if(ww.lt.obs_dist_minwt)ww=obs_dist_minwt
                    wt=ww*wt
                  end if
                  if(weight_multiplier.gt.unassigned_thresh) wt=wt*weight_multiplier
                  iprior=iprior+1
                  anynewprior=anynewprior+1
                  call num2char(iprior,aprior)
                  priorname='gr'//trim(aprior)
                  aapar1=arpar_1(irpar)
                  aapar2=arpar_2(jrpar)
                  if(ilog_1.eq.0)then
                    write(tempunit,930) trim(priorname),trim(aapar1), &
                    trim(aapar2),diffrat,wt,trim(reg_group)
930                 format(1x,a,t15,'1.0 * ',a,' - 1.0 * ',a,   &
                    ' = ',1pg13.6,1x,1pg13.6,1x,a)
                  else
                    write(tempunit,940) trim(priorname),trim(aapar1), &
                    trim(aapar2),diffrat,wt,trim(reg_group)
940                 format(1x,a,t15,'1.0 * log(',a,') - 1.0 * log(',a,   &
                    ') = ',1pg13.6,1x,1pg13.6,1x,a)
                  end if
                end do

              else if(diffrat_type.eq.2)then

                 call read_pp_file(ifail,nrpar_1,arpar_1,east_par_1,  &
                 north_par_1,val_par_1,parindex_1,family_prefix_1,    &
                 pilot_points_filename_1)
                 if(ifail.ne.0) go to 9890
                 call read_pp_file(ifail,nrpar_2,arpar_2,east_par_2,  &
                 north_par_2,val_par_2,parindex_2,family_prefix_2,    &
                 pilot_points_filename_2)
                 if(ifail.ne.0) go to 9890

                 iwarn=0
                 jprior=0
                 do ii=1,2
                   if(ii.eq.1)then
                     nr1=nrpar_1
                     nr2=nrpar_2
                   else
                     nr1=nrpar_2
                     nr2=nrpar_1
                   end if
                   do ir1=1,nr1
                     if(ii.eq.1)then
                       east1=east_par_1(ir1)
                       north1=north_par_1(ir1)
                       aapar1=arpar_1(ir1)
                     else
                       east1=east_par_2(ir1)
                       north1=north_par_2(ir1)
                       aapar1=arpar_2(ir1)
                     end if
                     do ir2=1,nr2
                       if(ii.eq.1)then
                         east2=east_par_2(ir2)
                         north2=north_par_2(ir2)
                       else
                         east2=east_par_1(ir2)
                         north2=north_par_1(ir2)
                       end if
                       ddist(ir2)=(east1-east2)*(east1-east2)+   &
                                  (north1-north2)*(north1-north2)
                       ddist(ir2)=sqrt(ddist(ir2))
                     end do
                     iparcount=0
                     do i=1,max_pilot_points
                       mindist=1.0d300
                       jmindist=0
                       do ir2=1,nr2
                         if(ddist(ir2).gt.search_radius)cycle
                         if(ddist(ir2).lt.mindist)then
                           mindist=ddist(ir2)
                           jmindist=ir2
                         end if
                       end do
                       if(jmindist.eq.0) go to 1050
                       iparcount=iparcount+1
                       if(ii.eq.1)then
                         jprior=jprior+1
                         if(jprior.gt.MAXLINKAGE)then
                           write(amessage,950)
950                        format(' Too many regularisation linkages - ', &
                           'increase MAXLINKAGE and re-compile program.')
                           go to 9890
                         end if
                         linkage1(jprior)=ir1
                         linkage2(jprior)=jmindist
                       else
                         if(jprior.gt.0)then
                           do jp=1,jprior
                             if(linkage1(jp).eq.jmindist)then
                               if(linkage2(jp).eq.ir1)then
                                 go to 990
                               end if
                             end if
                           end do
                         end if
                       end if
                       if(ii.eq.2)then
                         aapar2=arpar_1(jmindist)
                       else
                         aapar2=arpar_2(jmindist)
                       end if
                       ddist(jmindist)=1.0d301
                       ww=weight
                       if(weight_obs_dist.eq.1)then
                         do iii=1,2
                           if(iii.eq.1)then
                             eee=east1
                             nnn=north1
                           else
                             if(ii.eq.1)then
                               eee=east_par_2(jmindist)
                               nnn=north_par_2(jmindist)
                             else
                               eee=east_par_1(jmindist)
                               nnn=north_par_1(jmindist)
                             end if
                           end if
                           smalldist=1.0d300
                           do iobs=1,mobs
                             dist=(eee-east_obs(iobs))*    &
                                  (eee-east_obs(iobs))+    &
                                  (nnn-north_obs(iobs))*  &
                                  (nnn-north_obs(iobs))
                             if(dist.lt.smalldist)smalldist=dist
                           end do
                           if(iii.eq.1)then
                             distmin1=sqrt(smalldist)
                           else
                             distmin2=sqrt(smalldist)
                           end if
                         end do
                         distmin=max(distmin1,distmin2)
                         wtd=obs_dist_a+obs_dist_b*(distmin**obs_dist_c)
                         if(wtd.gt.obs_dist_maxwt)wtd=obs_dist_maxwt
                         if(wtd.lt.obs_dist_minwt)wtd=obs_dist_minwt
                         ww=ww*wtd
                       end if
                       if(weight_multiplier.gt.unassigned_thresh)ww=ww*weight_multiplier
                       iprior=iprior+1
                       anynewprior=anynewprior+1
                       call num2char(iprior,aprior)
                       priorname='gr'//trim(aprior)
                       diffrat=value
                       if(ilog_1.eq.1)then
                         if(diffrat.le.0)then
                           write(amessage,910)
                           go to 9890
                         end if
                         diffrat=log10(diffrat)
                       end if
                       if(ii.eq.1)then
                         if(ilog_1.eq.0)then
                           write(tempunit,930) trim(priorname),trim(aapar1), &
                           trim(aapar2),diffrat,ww,trim(reg_group)
                         else
                           write(tempunit,940) trim(priorname),trim(aapar1), &
                           trim(aapar2),diffrat,ww,trim(reg_group)
                         end if
                       else
                         if(ilog_1.eq.0)then
                           write(tempunit,930) trim(priorname),trim(aapar2), &
                           trim(aapar1),diffrat,ww,trim(reg_group)
                         else
                           write(tempunit,940) trim(priorname),trim(aapar2), &
                           trim(aapar1),diffrat,ww,trim(reg_group)
                         end if
                       end if
990                    continue
                     end do
1050                 continue
                     if(iparcount.lt.min_pilot_points)then
                        if(warn_less_min.ne.1)then
                          write(amessage,960) trim(aapar1)
960                       format(' Less than MIN_PILOT_POINTS parameters from ', &
                          'other parameter family lie ', &
                          'within one SEARCH_RADIUS of parameter "',a,'".')
                          go to 9890
                        else
                          iwarn=iwarn+1
                          if(iwarn.gt.MAXWARN)then
                            write(amessage,965)
965                         format(' Too many unmatched parameters: increase ', &
                            'MAXWARN and re-compile program.')
                            go to 9890
                          end if
                          parwarn(iwarn)=aapar1
                        end if
                     end if

                   end do
                 end do
                 if(iwarn.ne.0)then
                   write(amessage,970)
970                format(' Warning: the following parameters cannot be matched ', &
                   'to any parameters from other family within the ',  &
                   'user-supplied SEARCH_RADIUS:-')
                   call write_message(leadspace='yes')
                   write(6,980) (trim(parwarn(i)),i=1,iwarn)
980                format(5(a13))
                   write(6,*)
                 end if

              end if

            else if(spec_type.eq.3)then

! - "Indiv-family" regularisation.

              if(paramname.eq.' ')then
                write(amessage,1051)
1051            format(' If SPEC_TYPE is set to "indiv_family" then a PARAMETER ', &
                'keyword must be supplied in REGSPEC block.')
                go to 9890
              else
                do ipar=1,npar
                  if(paramname.eq.apar(ipar)) go to 1060
                end do
                write(amessage,1055) trim(paramname)
1055            format(' No PARAMETER named "',a,'" is present in PEST control file.')
                go to 9890
1060            continue
                paramindex=ipar
                paramtrans=itrans(ipar)
                if(paramtrans.lt.0)then
                  write(amessage,1070)
1070              format(' PARAMETER specified in REGSPEC block must not be tied ', &
                  'or fixed in PEST control file.')
                  go to 9890
                end if
              end if
              if(family_prefix.eq.' ')then
                write(amessage,1080)
1080            format(' If SPEC_TYPE is set to "indiv_family" then a FAMILY_PREFIX ', &
                'keyword must be present in REGSPEC block.')
                go to 9890
              else
                call get_regularised_parameters(ifail,npar,MAXREGPAR,nrpar,  &
                family_prefix,itrans,parindex,apar,arpar,pestfile)
                if(ifail.ne.0) go to 9890
              end if
              if((reg_type.ne.2).and.(reg_type.ne.3))then
                write(amessage,1090)
1090            format(' If SPEC_TYPE is set to "indiv_family" then REG_TYPE must be ', &
                'set to "difference" or "ratio" in REGSPEC block.')
                go to 9890
              end if
              call log_test(ilog,npar,nrpar,parindex,itrans)
              if(ilog.lt.0)then
                write(amessage,1100)
1100            format(' If SPEC_TYPE is set to "indiv_family" then all non-tied ', &
                'and non-fixed parameters ', &
                'belonging to selected parameter family must be either log-', &
                'transformed or untransformed.')
                go to 9890
              else
                if(ilog.eq.0)then
                  if(reg_type.eq.3)then
                    write(amessage,1110)
1110                format(' If REG_TYPE is set to "ratio" then all adjustable ', &
                    'parameters in selected parameter family must be log-transformed.')
                    go to 9890
                  end if
                else if(ilog.eq.1)then
                  if(reg_type.eq.2)then
                    write(amessage,1120)
1120                format(' If REG_TYPE is set to "difference" then all adjustable ', &
                    'parameters in selected parameter family must be ', &
                    'untransformed.')
                    go to 9890
                  end if
                end if
              end if
              if(reg_type.eq.2)then
                if(paramtrans.eq.1)then
                  write(amessage,1122)
1122              format(' If REG_TYPE is set to "difference" then specified PARAMETER ', &
                  'in REGSPEC block must NOT be log-transformed in PEST control file.')
                  go to 9890
                end if
              else if(reg_type.eq.3)then
                if(paramtrans.eq.0)then
                  write(amessage,1123)
1123              format(' If REG_TYPE is set to "ratio" then specified PARAMETER ', &
                  'in REGSPEC block must be log-transformed in PEST control file.')
                  go to 9890
                end if
                if(diffrat_val_type.eq.1)then
                  if(value.lt.unassigned_thresh)then
                    write(amessage,1150)
1150                format(' If DIFFRAT_VAL_TYPE is set to "uniform" then ', &
                    'a VALUE keyword must be present in REGSPEC block.')
                    go to 9890
                  end if
                  if(value.le.0.0)then
                    write(amessage,1125)
1125                format(' If REG_TYPE is set to "ratio" and DIFFRAT_VAL_TYPE ', &
                    'is set to "uniform", then VALUE must be a positive number ', &
                    'in REGSPEC block.')
                    go to 9890
                  end if
                end if
              end if
              if(diffrat_val_type.eq.unassigned_int)then
                write(amessage,1140)
1140            format(' If SPEC_TYPE is set to "indiv_family" then a ', &
                'DIFFRAT_VAL_TYPE keyword must be present in REGSPEC block.')
                go to 9890
              else
                if(diffrat_val_type.eq.2)then
                  if(pilot_points_filename_dr.eq.' ')then
                    write(amessage,1160)
1160                format(' If DIFFRAT_VAL_TYPE is set to "file" then ', &
                    'a PILOT_POINTS_FILENAME_DR keyword must be present in ', &
                    'REGSPEC block.')
                    go to 9890
                  end if
                end if
              end if
              if(weight_type.ne.1)then
                write(amessage,1162)
1162            format(' If SPEC_TYPE is set to "indiv_family" then WEIGHT_TYPE ', &
                'must be set to "uniform" in REGSPEC block.')
                go to 9890
              end if
              if(weight_obs_dist.eq.1)then
                if(pilot_points_filename_dr.eq.' ')then
                  write(amessage,1170)
1170              format(' If SPEC_TYPE is "indiv_family" and WEIGHT_OBS_DIST is "yes" ', &
                  'then a PILOT_POINTS_FILENAME_DR keyword must be present in ', &
                  'REGSPEC block.')
                  go to 9890
                end if
              end if

              if((weight_obs_dist.eq.1).or.(diffrat_val_type.eq.2))then
                call read_pp_file(ifail,nrpar,arpar,east_par,north_par,  &
                val_par,parindex,family_prefix,pilot_points_filename_dr)
                if(ifail.ne.0) go to 9890
              end if
              if(weight_obs_dist.eq.1)then
                if(obs_coord_filename.eq.' ')then
                  write(amessage,1175)
1175              format(' If WEIGHT_OBS_DIST is set to "yes" then an ',  &
                  'OBS_COORD_FILENAME must be supplied in ', &
                  'REGSPEC block.')
                  go to 9890
                end if
                call read_oc_file(ifail,MAXOBS,mobs,east_obs,north_obs,  &
                obs_coord_filename)
                if(ifail.ne.0) go to 9890
              end if

              do irpar=1,nrpar
                iprior=iprior+1
                anynewprior=anynewprior+1
                call num2char(iprior,aprior)
                priorname='gr'//trim(aprior)
                wt=weight
                if(weight_obs_dist.eq.1)then
                  distmin=1.0d300
                  do iobs=1,mobs
                    dist=(east_par(irpar)-east_obs(iobs))*    &
                         (east_par(irpar)-east_obs(iobs))+    &
                         (north_par(irpar)-north_obs(iobs))*  &
                         (north_par(irpar)-north_obs(iobs))
                    if(dist.lt.distmin)distmin=dist
                  end do
                  distmin=sqrt(distmin)
                  ww=obs_dist_a+obs_dist_b*(distmin**obs_dist_c)
                  if(ww.gt.obs_dist_maxwt)ww=obs_dist_maxwt
                  if(ww.lt.obs_dist_minwt)ww=obs_dist_minwt
                  wt=ww*wt
                end if
                ww=wt
                if(weight_multiplier.gt.unassigned_thresh)then
                  ww=ww*weight_multiplier
                end if
                if(diffrat_val_type.eq.1) val_par(irpar)=value
                if(paramtrans.eq.0)then
                  write(tempunit,1180) trim(priorname),trim(arpar(irpar)), &
                  trim(paramname),val_par(irpar),ww,trim(reg_group)
1180              format(1x,a,t15,'1.0 * ',a,' - 1.0 * ',a,' = ',     &
                  1pg14.7,1x,1pg13.6,1x,a)
                else
                  if(val_par(irpar).le.0.0)then
                    write(amessage,1190) trim(arpar(irpar))
1190                format(' Value supplied for regularisation ratio for parameter "',  &
                    a,'" in ', &
                    'PILOT_POINTS_FILENAME_DR file must be positive.')
                    go to 9890
                  end if
                  write(tempunit,1200) trim(priorname),trim(arpar(irpar)), &
                  trim(paramname),log10(val_par(irpar)),ww,trim(reg_group)
1200              format(1x,a,t15,'1.0 * log(',a,') - 1.0 * log(',a,') = ', &
                  1pg14.7,1x,1pg13.6,1x,a)
                end if
              end do

            end if

            if(anynewprior.ne.0)then
              ifound=0
              do i=1,nobsgp
                if(reg_group.eq.oldgroup(i)) go to 560
              end do
              if(numnewgroup.ne.0)then
                do i=1,numnewgroup
                  if(reg_group.eq.newgroup(i)) go to 560
                end do
              end if
              numnewgroup=numnewgroup+1
              if(numnewgroup.gt.MAXNEWGROUP)then
                write(amessage,552)
552             format(' Increase MAXNEWGROUP and re-compile program.')
                go to 9890
              end if
              newgroup(numnewgroup)=reg_group
560           continue
            else
              write(amessage,565)
565           format(' No new items of prior information were added by current ', &
              'REGSPEC block: check its contents.')
              go to 9890
            end if

          end if
          if(blocktype.eq.1)then
            ablock='REGSPEC'
          else
            ablock='REGULARISATION'
          end if
          write(6,551) trim(ablock)
551       format(1x,a,' block processed ok.')

        end do


8000    continue
        if(nblock.eq.0)then
          write(amessage,8010) trim(afile)
8010      format(' No REGSPEC blocks have been found in file ',a,'.')
          go to 9890
        end if
        close(unit=tempunit)

! -- All that remains to be done now is to write the new pest control file.

        call addquote(outfile,bfile)
        write(6,8015) trim(bfile)
8015    format(/,' Writing PEST control file ',a,'...')
        do i=1,3
          read(pestunit,'(a)') cline
          write(outunit,'(a)',err=9800) trim(cline)
        end do
        read(pestunit,'(a)') cline
        nprior=nprior+iprior
        write(outunit,8020,err=9800) npar,nobs,npargp,nprior,nobsgp+numnewgroup
8020    format(5i8)
        iline=4
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000) cline
          write(outunit,'(a)',err=9800) trim(cline)
          if(index(cline,'* observation gr').ne.0) go to 8030
        end do
8030    continue
        do i=1,nobsgp
          read(pestunit,'(a)',err=9000) cline
          write(outunit,'(a)',err=9800) trim(cline)
        end do
        if(numnewgroup.ne.0)then
          do i=1,numnewgroup
            write(outunit,'(a)',err=9800) trim(newgroup(i))
          end do
        end if
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9050) cline
          if(cline(1:1).eq.'*')then
            if(index(cline,'* prior information').ne.0)foundprior=1
          end if
          if(cline(1:6).eq.'* regu')go to 8080
          write(outunit,'(a)',err=9800) trim(cline)
        end do

8080    continue
        tempunit=nextunit()
        open(unit=tempunit,file=tempfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,8090) trim(tempfile)
8090      format(' Error re-opening temporary file ',a,'.')
          go to 9890
        end if
        if(foundprior.eq.0)then
          write(outunit,8095,err=9800)
8095      format('* prior information')
        end if
        do i=1,iprior
          read(tempunit,'(a)',end=8200) cline
          write(outunit,'(a)',err=9800) trim(cline)
        end do
8200    continue
        close(unit=tempunit)
        write(outunit,8210,err=9800)
8210    format('* regularisation')
        write(cline,8220) phimlim,phimaccept,fracphim
8220    format(1x,1pg14.6,1x,1pg14.6,1x,1pg14.6)
        if(memsave.eq.1)then
          cline=trim(cline)//' memsave'
        end if
        if(conjgrad.eq.1)then
          if(memsave.eq.0)then
            cline=trim(cline)//' nomemsave'
          end if
          lt1=len_trim(cline)+1
          write(cline(lt1:),8230) cgrtol,cgitnlim
8230      format(' cg ',1pg14.6,1x,i5)
        end if
        write(outunit,'(a)',err=9800) trim(cline)
        write(cline,8240) wfinit,wfmin,wfmax
8240    format(3(1x,1pg13.6))
        if(linreg.eq.1)then
          cline=trim(cline)//' linreg'
        end if
        write(outunit,'(a)',err=9800) trim(cline)
        write(outunit,8250) wffac,wftol,iregadj
8250    format(2(1pg13.6),i5)

        close(unit=tempunit)
        close(unit=outunit)
        write(6,8260) trim(outfile)
8260    format(' - file ',a,' written ok.')


        go to 9900


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(afile)
9060    format(' Unexpected end encountered to file ',a,'.')
        go to 9890
9070    write(amessage,9080) trim(afile)
9080    format(' Cannot find "parameter data" section of PEST control file ',a,'.')
        go to 9890
9090    write(amessage,9095) trim(afile)
9095    format(' Cannot find "observation groups" section of PEST control file ',a,'.')
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Insufficient entries on line ',a,' of file ',a)
        go to 9890
9150    write(amessage,9160)
9160    format(' Cannot allocate sufficient memory to continue ',&
        'execution.')
        go to 9890
9170    write(amessage,9180) trim(afile)
9180    format(' File ',a,' is an illegal PEST control file. Check it with ', &
        'PESTCHEK.')
        go to 9890
9200    write(amessage,9210) trim(afile)
9210    format(' No "regularisation" section found in PEST control file ',a,'.')
        go to 9890
9250    write(amessage,9260) trim(aline),trim(afile)
9260    format(' Error reading "regularisation" data from existing PEST ', &
        'control file; error occurs at line ',a,' of file ',a,'.')
        go to 9890
9300    write(amessage,9310) trim(aline),trim(afile)
9310    format(' Two entries are exected on line ',a,' of file ',a,'.')
        go to 9890
9350    write(amessage,9360) trim(aline),trim(afile)
9360    format(' Cannot read real number as second entry on line ',a,  &
        ' of file ',a,'.')
        go to 9890
9370    write(amessage,9380) trim(aline),trim(afile)
9380    format(' Cannot read integer as second entry on line ',a,  &
        ' of file ',a,'.')
        go to 9890
9400    write(amessage,9410) trim(aline),trim(afile)
9410    format(' Keyword at line ',a,' of file ',a,' is only allowed in REGULARISATION block.')
        go to 9890
9450    write(amessage,9460) trim(aline),trim(afile)
9460    format(' Second entry on line ',a' of file ',a,' must be a ',  &
        'positive number.')
        go to 9890
9500    write(amessage,9510) trim(aline),trim(afile)
9510    format(' Second item is illegal on line ',a,' of file ',a,'.')
        go to 9890
9550    write(amessage,9560) trim(aline),trim(afile)
9560    format(' Keyword at line ',a,' of file ',a,' is only allowed in REGSPEC block.')
        go to 9890
9660    write(amessage,9670) trim(aline),trim(afile)
9670    format(' Error reading text string from line ',a,' of file ',a,'.')
        go to 9890
9700    write(amessage,9710) trim(aline),trim(afile)
9710    format(' Second entry on line ',a,' of file ',a,' must be "yes" or "no".')
        go to 9890
9750    write(amessage,9760) trim(aline),trim(afile)
9760    format(' Second entry on line ',a' of file ',a,' must be a ',  &
        'non-negative number.')
        go to 9890
9800    write(amessage,9810) trim(bfile)
9810    format(' Error writing to file ',a,'.')
        go to 9890
9850    write(amessage,9860) trim(aline),trim(afile)
9860    format(' Keyword repeated in block at line ',a,' of file ',a,'.')
        go to 9890

9890    call write_message(leadspace='yes')

9900    call close_files
        deallocate(apar,itrans,oldgroup,stat=ierr)
        write(6,*)

end program genreg



subroutine read_oc_file(ifail,MAXOBS,nobs,east_obs,north_obs,obs_coord_filename)

! -- Subroutine read_oc_file reads the eastings and northings from an
!    observation coordinates file.

        use defn
        use inter

        implicit none

        integer, intent(out)          :: ifail
        integer, intent(in)           :: MAXOBS
        integer, intent(out)          :: nobs
        double precision, intent(out) :: east_obs(MAXOBS),north_obs(MAXOBS)
        character (len=*)             :: obs_coord_filename

        integer                       :: ierr,iunit,iline
        character (len=10)            :: aline
        character (len=200)           :: afile

        ifail=0
        iunit=nextunit()
        call addquote(obs_coord_filename,afile)
        write(6,10) trim(afile)
10      format( '      - reading file ',a,'...')
        open(unit=iunit,file=obs_coord_filename,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,30) trim(afile)
30        format(' Cannot open observation coordinates file ',a,'.')
          go to 9890
        end if

        nobs=0
        iline=0
        do
          iline=iline+1
          read(iunit,'(a)',end=100,err=9000) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).eq.'#') cycle
          nobs=nobs+1
          if(nobs.gt.MAXOBS)then
            write(amessage,50) trim(afile)
50          format(' Insufficient memory; too many observations ', &
            'in file ',a,' - increase MAXOBS and re-compile program.')
            go to 9890
          end if
          call linesplit(ifail,3)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,60) trim(aline),trim(afile)
60          format(' Insufficient entries on line ',a,' of file ',a,'.')
            go to 9890
          end if
          east_obs(nobs)=char2double(ifail,2)
          if(ifail.ne.0) go to 9000
          north_obs(nobs)=char2double(ifail,3)
          if(ifail.ne.0) go to 9000
        end do
100     close(unit=iunit)
        call num2char(nobs,aline)
        write(6,110) trim(aline),trim(afile)
110     format('      - ',a,' observation well coordinates read from file ',a)
        return

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890

9890    ifail=1

        return

end subroutine read_oc_file



subroutine read_pp_file(ifail,nrpar,arpar,east_par,north_par,    &
                    val_par,parindex,family_prefix, pilot_points_filename)

! -- Subroutine read_pp_file reads a pilot points file. It extracts coordinates
!    and values from that file pertaining to identified parameters.

        use defn
        use inter
        implicit none

        integer, intent(out)             :: ifail
        integer, intent(inout)           :: nrpar
        character (len=*), intent(inout) :: arpar(nrpar)
        double precision, intent(out)    :: east_par(nrpar),north_par(nrpar)
        real, intent(out)                :: val_par(nrpar)
        integer, intent(inout)           :: parindex(nrpar)
        character (len=*)                :: family_prefix
        character (len=*), intent(in)    :: pilot_points_filename

        integer                       :: iunit,ierr,icount,j,lt1,ibeg,irpar,iline
        character (len=10)            :: aline
        character (len=20)            :: app
        character (len=200)           :: afile


! -- Initialisation

        ifail=0
        icount=0
        do irpar=1,nrpar
          val_par(irpar)=-1.1e37
        end do

! -- The file is opened and processed line by line.

        iunit=nextunit()
        call addquote(pilot_points_filename,afile)
        write(6,10) trim(afile)
10      format( '      - reading file ',a,'...')
        open(unit=iunit,file=pilot_points_filename,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,30) trim(afile)
30        format(' Cannot open pilot points file ',a,'.')
          go to 9890
        end if

        lt1=len_trim(family_prefix)+1
        iline=0
        do
          iline=iline+1
          read(iunit,'(a)',end=100,err=9000) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).eq.'#') cycle
          call linesplit(ifail,5)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,60) trim(aline),trim(afile)
60          format(' Insufficient entries on line ',a,' of file ',a,'; five ',  &
            'entries per line are expected in a pilot points file.')
            go to 9890
          end if
          app=cline(left_word(1):right_word(1))
          call casetrans(app,'lo')
          do irpar=1,nrpar
            if((arpar(irpar)(lt1:).eq.app).or.(arpar(irpar).eq.app))then
              east_par(irpar)=char2double(ifail,2)
              if(ifail.ne.0) go to 9000
              north_par(irpar)=char2double(ifail,3)
              if(ifail.ne.0) go to 9000
              val_par(irpar)=char2real(ifail,5)
              if(ifail.ne.0) go to 9000
              icount=icount+1
              go to 80
            end if
          end do
80        continue
        end do

! -- The pilot points file has been read to completion. Now the number of
!    pertinent pilot point parameters is reduced if necessary.
!    Note that if pilot point names are duplicated, the above reading scheme
!    only reads data pertaining to the first pilot point.

100     close(unit=iunit)
        if(icount.eq.0)then
          write(amessage,110) trim(afile),trim(family_prefix)
110       format(' There are no pilot points in file ',a,' which match ', &
          'any parameters whose prefix is "',a,'" from PEST control file.')
          go to 9890
        end if
        if(icount.ne.nrpar)then
          ibeg=1
120       continue
          do irpar=ibeg,nrpar
            if(val_par(irpar).lt.-1.0e37) go to 130
          end do
          go to 150
130       continue
          write(amessage,131) trim(arpar(irpar)),trim(afile)
131       format(' No match found for parameter "',a,'" in pilot points file ',a,'.')
          go to 9890
! -- The following code existed as an alternative to the above error message.
          if(irpar.eq.nrpar)then
            nrpar=nrpar-1
            go to 150
          else
            do j=irpar,nrpar-1
              arpar(j)=arpar(j+1)
              val_par(j)=val_par(j+1)
              east_par(j)=east_par(j+1)
              north_par(j)=north_par(j+1)
              parindex(j)=parindex(j+1)
            end do
            nrpar=nrpar-1
            ibeg=irpar
            go to 120
          end if
        end if

150     continue
        call num2char(nrpar,aline)
        write(6,170) trim(aline),trim(afile)
170     format('      - ',a,' matching pilot points read from file ',a)
        return


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890

9890    continue
        ifail=1
        return

end subroutine read_pp_file



subroutine log_test(ilog,npar,nrpar,parindex,itrans)

        implicit none

        integer, intent(out) :: ilog
        integer, intent(in)  :: npar
        integer, intent(in)  :: nrpar
        integer, intent(in)  :: parindex(nrpar)
        integer, intent(in)  :: itrans(npar)

        integer              :: icountl,icountn,irpar

        icountl=0
        icountn=0
        do irpar=1,nrpar
          if(itrans(parindex(irpar)).eq.1)then
            icountl=icountl+1
          else if(itrans(parindex(irpar)).eq.0)then
            icountn=icountn+1
          end if
        end do
        if((icountl.eq.0).and.(icountn.eq.0))then
          ilog=-10
          return
        end if
        if((icountl.gt.0).and.(icountn.gt.0))then
          ilog=-1
          return
        end if
        if(icountl.gt.0)then
          ilog=1
        else
          ilog=0
        end if

        return

end subroutine log_test



subroutine get_regularised_parameters(ifail,npar,MAXREGPAR,nrpar,family_prefix,  &
               itrans,parindex,apar,arpar,pestfile)


! -- Subroutine GET_REGULARISED_PARAMETERS selects parameter values on the basis
!    of a parameter prefix.

        use defn
        use inter
        implicit none

        integer, intent(out)           :: ifail
        integer, intent(in)            :: npar
        integer, intent(in)            :: MAXREGPAR
        integer, intent(out)           :: nrpar
        character (len=*), intent(in)  :: family_prefix
        integer, intent(in)            :: itrans(npar)
        integer, intent(out)           :: parindex(MAXREGPAR)
        character (len=*), intent(in)  :: apar(npar)
        character (len=*), intent(out) :: arpar(MAXREGPAR)
        character (len=*), intent(in)  :: pestfile


        integer             :: nf,ipar
        character (len=200) :: bfile

        ifail=0

        nf=len_trim(family_prefix)
        nrpar=0
        do ipar=1,npar
          if(apar(ipar)(1:nf).eq.family_prefix(1:nf))then
            if(itrans(ipar).ge.0)then
              nrpar=nrpar+1
              if(nrpar.gt.MAXREGPAR)then
                write(amessage,40)
40              format(' Too many parameters fit prefix - increase MAXREGPAR ', &
                'and re-compile program.')
                go to 9890
              end if
              arpar(nrpar)=apar(ipar)
              parindex(nrpar)=ipar
            end if
          end if
        end do
        if(nrpar.eq.0)then
          call addquote(pestfile,bfile)
          write(amessage,100) trim(bfile),trim(family_prefix)
100       format(' No adjustable (i.e. non-tied and nonfixed) parameters ', &
          'in PEST control file ',a,   &
          ' have a prefix of "',a,'".')
          go to 9890
        end if
        return

9890    continue
        ifail=1
        return

end subroutine get_regularised_parameters

