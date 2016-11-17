!     Last change:  JD   28 Aug 2006   11:34 pm

program fac2mf2k

! -- Program FAC2MF2K adds parameters to a MODFLOW 2000 input dataset, based on a set of
!    pilot points. Note that PPK2FAC must have been run first.

	use defn
	use inter

	implicit none
        logical       :: lexist,lopened
        integer       :: ifail,idate,iheader,iline,ll,ncol,nrow,nbb,ierr,itemp, &
                         j,i,nlay,mrow,mcol,ip,nclu,icl,layer,jreplace,mclu,mplpf, &
                         ipoint,inum,ilay,ncopylpf,nplpf,iwdflg,irep,ipp,na,icellno,izone, &
                         newnplpf,iirep,iz,nzn,iwrite,newnzn,irflag,ncopyzon,irow,icol, &
                         im,nml,newnml,lt,newmul,ncopymul,iarr,itemp1,ifunct,iprn,nplist, &
                         iparout,newnplist,ncopysen,ln,isens,maxiter,maxchange, &
                         ibeflg,iycflg,iostar,nopt,nfit,iap,iprcov,iprint,lprint,lastx, &
                         npng,ipr,mpr,iwtp,ii,ie,jj,numarray,itrans,itemptrans,ncopyint, &
                         nrp,nar,inew,itr,newnar,iflagarr,ncopynam,nnppt,iipt,ipestint, &
                         iinterp,nostop,hdrybot,limop,jjflag,kk,nper,iper,ncopyrch, &
                         rchoutunit,numrchpar,numlpfpar,nprch,nrchop,irchcb,newnprch, &
                         jrepce,inrechold,inrech,inirch,rchunit,itrexist
        integer       :: nameunit,facunit,tunit,disunit,lpfunit,zonunit,modzonunit, &
                         modmulunit,mulunit,senunit,pesunit,maxminunit
        integer, allocatable, dimension(:)     :: mdat,nnclu,llayer,ipt,iused,found, &
                                                  removed_par,removed_zon,removed_mul
        integer, allocatable, dimension(:,:)   :: intarray,ia
        real                                   :: rtemp,parval,b,bl,bu,bscal,tol,sosc, &
                                                  sosr,rmar,rmarm,csa,fconv,minthick
        real, allocatable, dimension(:)        :: wt
        real, allocatable, dimension(:)        :: rmdat
        real, allocatable, dimension(:,:)      :: ra
        real, allocatable, dimension(:,:,:)    :: amul
        character (len=1)                      :: facformat,aao
        character (len=3)                      :: abase
        character (len=10)                     :: aline,atemp,parnam,partyp,anum, &
                                                  mltarr,zonarr,zonnam,mltnam,afunction, &
                                                  arpt
        character (len=15)                     :: atempp,arrmin,arrmax,armi,arma
        character (len=120)                    :: namefile,afile,senfile,lpffile, &
                                                  pesfile,intfile,aprompt,facfile, &
                                                  disfile,bfile,zonfile,mulfile, &
                                                  maxminfile,interpfile,pestctlfile, &
                                                  rchfile
        character (len=1000)                   :: ccline
        character (len=15),allocatable, dimension(:)   :: oldarrayname(:), &
                                                          oldarrmin(:),oldarrmax(:)
        character (len=12), allocatable, dimension(:)  :: wpoints



	write(amessage,5)
5	format(' Program FAC2MF2K adds parameters to a MODFLOW 2000 input dataset, ', &
        'based on a set of pilot points. Note that PPK2FAC must have been run first.')
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
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message
	  go to 9900
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if


! -- Initialisation

        lpffile=' '
        rchfile=' '
        senfile=' '
        pesfile=' '
        intfile=' '
        interpfile=' '
        pestctlfile=' '
        disfile=' '
        zonfile=' '
        mulfile=' '


        ncopylpf=0
        ncopyrch=0
        ncopyzon=0
        ncopymul=0
        ncopysen=0
        ncopyint=0
        ncopynam=0

        zonunit=0
        mulunit=0

        iprn=0
        numparamreplace=0

        ipestint=0
        iinterp=1
        nostop=0
        hdrybot=0
        limop=0
        minthick=0.0
        tunit=1
        itrexist=0

10      aprompt=' Enter name of MODFLOW 2000 name file: '
        call open_input_file(ifail,aprompt,namefile,nameunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900


! -- The MODFLOW 2000 name file is read and pertinent filenames are extracted.

        afile=namefile
        iline=0
100     iline=iline+1
        read(nameunit,'(a)',end=120,err=9000) cline
        if(cline.eq.' ') go to 100
        if(cline(1:1).eq.'#') go to 100
        call linesplit(ifail,3)
        atemp=cline(left_word(1):right_word(1))
        call casetrans(atemp,'lo')
        if(atemp.eq.'lpf')then
          if(ifail.ne.0)go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,lpffile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
        else if(atemp.eq.'rch')then
          if(ifail.ne.0)go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,rchfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
          rchoutunit=char2int(ifail,2)
        else if(atemp.eq.'sen')then
          if(ifail.ne.0)go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,senfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
        else if(atemp.eq.'pes')then
          if(ifail.ne.0)go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,pesfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
        else if(atemp.eq.'dis')then
          if(ifail.ne.0) go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,disfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
        else if(atemp.eq.'zone')then
          if(ifail.ne.0) go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,zonfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
          modzonunit=char2int(ifail,2)
          if(ifail.ne.0)then
            write(amessage,102) trim(namefile)
102         format(' Cannot read unit number for zone file from MODFLOW name file ',a)
            go to 9890
          end if
        else if(atemp.eq.'mult')then
          if(ifail.ne.0) go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,mulfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
          modmulunit=char2int(ifail,2)
          if(ifail.ne.0)then
            write(amessage,103) trim(namefile)
103         format(' Cannot read unit number for multiplier file from MODFLOW name file ',a)
            go to 9890
          end if
        else if(atemp.eq.'asp')then
          if(ifail.ne.0) go to 9050
          ll=len_trim(cline)
          call getfile(ifail,cline,interpfile,left_word(3),ll)
          if(ifail.ne.0) go to 9100
        end if
        go to 100
120     continue
        close(unit=nameunit)
        if((lpffile.eq.' ').and.(rchfile.eq.' '))then
          write(amessage,140) trim(namefile)
140       format(' An input file for neither the LPF package nor the RCH package is ', &
          'cited in MODFLOW 2000 name file ',a)
          go to 9890
        end if
!        if(senfile.eq.' ')then
!          write(amessage,150) trim(namefile)
!150       format(' An input file for the MODFLOW 2000 sensitivity process is not ', &
!          'cited in MODFLOW 2000 name file ',a)
!          go to 9890
!        end if
        if(disfile.eq.' ')then
          write(amessage,155) trim(namefile)
155       format(' The name of a MODFLOW 2000 discretisation file is not ', &
          'cited in MODFLOW 2000 name file ',a)
          go to 9890
        end if
        if(zonfile.eq.' ')then
          write(amessage,156) trim(namefile)
156       format(' The name of a MODFLOW 2000 zone file is not ', &
          'cited in MODFLOW 2000 name file ',a)
          go to 9890
        end if
        if(mulfile.eq.' ')then
          write(amessage,157) trim(namefile)
157       format(' The name of a MODFLOW 2000 multiplier file is not ', &
          'cited in MODFLOW 2000 name file ',a)
          go to 9890
        end if

        write(6,*)
200     aprompt=' Enter name of PPK2FAC-generated interpolation factor file: '
        call open_input_file(ifail,aprompt,facfile,facunit,form_prompt='yes', &
        fformat=facformat)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          write(6,*)
          escset=0
          go to 10
        end if

        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9150,end=9200) afile
        else
          read(facunit,err=9150,end=9200) afile
        end if
        afile=adjustl(afile)
        nbb=len_trim(afile)
        call getfile (ifail,afile,pilot_points_file,1,nbb)
        if(ifail.ne.0)pilot_points_file=' '

        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9150,end=9200) afile
        else
          read(facunit,err=9150,end=9200) afile
        end if
        afile=adjustl(afile)
        nbb=len_trim(afile)
        call getfile (ifail,afile,intfile,1,nbb)
        if(ifail.ne.0)intfile=' '

        if(facformat.eq.'f')then
          read(facunit,*,err=9150,end=9200) ncol,nrow
        else
          read(facunit,err=9150,end=9200) ncol,nrow
        end if
        allocate(intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,250)
250	  format(' Cannot allocate sufficient memory to run FAC2REAL.')
	  go to 9890
	end if


        if(facformat.eq.'f')then
          read(facunit,*,err=9150,end=9200) nnppt
        else
          read(facunit,err=9150,end=9200) nnppt
        end if
        allocate(wpoints(nnppt),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,250)
	  go to 9890
	end if
        if(facformat.eq.'f')then
          do i=1,nnppt
            read(facunit,'(a)',err=9000,end=9200) wpoints(i)
          end do
        else
          do i=1,nnppt
            read(facunit,err=9000,end=9200) wpoints(i)
          end do
        end if
        do i=1,nnppt
          call casetrans(wpoints(i),'hi')
        end do


! -- The pilot points file is read.

280     call read_pilot_points_file(ifail, &
	' Enter name of pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(intarray,wpoints,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,290)
290	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
          close(unit=facunit)
	  write(6,*)
	  go to 200
	end if
        if(nnppt.ne.num_pilot_points)then
          write(amessage,59)
          go to 9890
        else
          do i=1,num_pilot_points
            if(pilot_point_id(i).ne.wpoints(i))then
              write(amessage,59)
59            format(' The pilot points in the pilot points file are not the same, ', &
              'or are not arranged in the same order, as the pilot points in the ', &
              'pilot points file read by PPF2FAC when it calculated the ', &
              'factors contained in the factor file.')
              go to 9890
            end if
          end do
        end if

! -- The integer array file is read.

300     aprompt=' Enter name of integer array zonation file '
        call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=nrow,columns=ncol,defaultfile=intfile)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 280
        end if
        intfile=aprompt

! -- The parameter replacement file is now read.

        write(6,*)
320     call read_parameter_replacement_file(ifail, &
	' Enter name of parameter replacement file: ',ncol,nrow)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 300
	end if


        write(6,*)
324     write(6,325,advance='no')
325     format(' Enter name base for new multiplier and zonation arrays ', &
        '[3 chars or less]: ')
        abase=' '
        read(5,'(a)') abase
        if(abase.eq.' ') go to 324
        if(index(eschar,abase(1:2)).ne.0) then
          write(6,*)
          go to 320
        end if
        call casetrans(abase,'hi')


        write(6,*)
326     maxminfile=' '
        afile=' '
321     continue
        if(interpfile.eq.' ')then
          write(6,318,advance='no')
318       format(' Enter name for ASP input file: ')
          read(5,'(a)') afile
          if(afile.eq.' ') go to 321
        else
          write(6,327,advance='no') trim(interpfile)
327       format(' Enter name for ASP input file [',a,']: ')
          read(5,'(a)') afile
        end if
        if(afile.eq.' ')then
          maxminfile=interpfile
        else
          if (index(eschar,afile(1:2)).ne.0) then
            write(6,*)
            go to 324
          end if
          call casetrans(afile,'lo')
          nbb=len_trim(afile)
          call getfile(ifail,afile,maxminfile,1,nbb)
          if(ifail.ne.0) go to 326
        end if
        if(interpfile.ne.' ')then
          if(maxminfile.ne.interpfile)then
            write(amessage,319) trim(interpfile)
319         format(' An ASP input file named ',a,' is already ', &
            'cited in the MODFLOW name file. You must use this same filename ', &
            'or remove the "ASP" entry from the MODFLOW-2000 name file.')
            call write_message(leadspace='yes')
            write(6,*)
            go to 326
          end if
        end if
        inquire(file=maxminfile,exist=lexist)
        aao='o'
        if(lexist)then
328       write(6,329,advance='no')
329       format(' This file already exists; append or overwrite [a/o]: ')
          read(5,'(a)') aao
          if(aao.eq.' ') go to 328
          if(index(eschar,aao).ne.0) then
            write(6,*)
            go to 321
          end if
          call casetrans(aao,'lo')
          if((aao.ne.'a').and.(aao.ne.'o')) go to 328
        end if

        if(aao.eq.'o')then
          maxminunit=nextunit()
          open(unit=maxminunit,file=maxminfile,iostat=ierr)
          if(ierr.ne.0) go to 326
        end if

! -- Before any parameters are manipulated, the number of layers in the current
!    model is ascertained. This is done by reading the discretisation file.

       disunit=nextunit()
       open(unit=disunit,file=disfile,status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,350) trim(disfile)
350      format(' Cannot open MODFLOW 2000 discretisation file ',a)
         go to 9890
       end if
360    read(disunit,'(a)',err=9350,end=9400) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 360
       call linesplit(ifail,4)
       if(ifail.ne.0)then
         write(amessage,370) trim(disfile)
370      format(' Insufficient entries on first data line of MODFLOW 2000 ', &
         'discretisation file ',a)
         go to 9890
       end if
       nlay=char2int(ifail,1)
       if(ifail.ne.0)then
         write(amessage,375) trim(disfile)
375      format(' Cannot read value for NLAY from MODFLOW discretisation file ',a)
         go to 9890
       end if
       mrow=char2int(ifail,2)
       if(ifail.ne.0)then
         write(amessage,376) trim(disfile)
376      format(' Cannot read value for NROW from MODFLOW discretisation file ',a)
         go to 9890
       end if
       mcol=char2int(ifail,3)
       if(ifail.ne.0)then
         write(amessage,377) trim(disfile)
377      format(' Cannot read value for NCOL from MODFLOW discretisation file ',a)
         go to 9890
       end if
       if((mrow.ne.nrow).or.(mcol.ne.ncol))then
         write(amessage,380)
380      format(' Grid dimensions as read from MODFLOW discretisation file do not ', &
         'match those used when generating pilot point interpolation factors using ', &
         'program PPK2FAC. ')
         go to 9890
       end if
       nper=char2int(ifail,4)
       if(ifail.ne.0)then
         write(amessage,378) trim(disfile)
378      format(' Cannot read value for NPER from MODFLOW discretisation file ',a)
         go to 9890
       end if
       close(unit=disunit)

       do i=1,numparamreplace
         do j=1,replace(i)%numout
           replace(i)%outparamlay1(j)=1
           replace(i)%outparamlay2(j)=nlay
         end do
       end do

! -- Some array allocation is carried out.

       itemp=max(nlay,nper)
       allocate(mdat(itemp),llayer(nlay),rmdat(nlay),stat=ierr)
       if(ierr.ne.0)then
         write(amessage,250)
         go to 9890
       end if

      allocate(ipt(num_pilot_points),wt(num_pilot_points),iused(num_pilot_points), &
      found(numparamreplace),stat=ierr)
      if(ierr.ne.0)then
        write(amessage,250)
        go to 9890
      end if
      iused=0         !iused is an array
      found=0         !found is an array

      iirep=0
      do irep=1,numparamreplace
        do j=1,replace(irep)%numout
          iirep=iirep+1
        end do
      end do
      allocate(removed_par(iirep),removed_zon(numzoneremove), &
      removed_mul(nummultremove),stat=ierr)
      if(ierr.ne.0)then
        write(amessage,250)
        go to 9890
      end if
      removed_par=0      !removed_par is an array
      removed_zon=0      !removed_zon is an array
      removed_mul=0      !removed_mul is an array

      allocate(ia(ncol,nrow),ra(ncol,nrow),stat=ierr)
      if(ierr.ne.0)then
        write(amessage,250)
        go to 9890
      end if

! -- The interpolation factor file is read in order to ascertain whether any pilot
!    points are not needed as parameters.

      itrans=-9999
      afile=facfile
      do
        if(facformat.eq.'f')then
          read(facunit,*,err=9850,end=384) icellno,itemptrans, &
          na,rtemp,((ipt(i),wt(i)),i=1,na)
        else
          read(facunit,err=9850,end=384)   icellno,itemptrans, &
          na,rtemp,((ipt(i),wt(i)),i=1,na)
        end if
        if(rtemp.ne.0.0)then
          write(amessage,383) trim(facfile)
383       format(' It is apparent from the interpolation factor file ',a, &
          ' that PPK2FAC undertook SIMPLE kriging to generate these factors. However ',&
          'only ORDINARY kriging is allowed if these factors are to be used by ', &
          'MODFLOW 2000.')
          go to 9890
        end if
        if(itrans.eq.-9999)then
          itrans=itemptrans
        else
          if(itrans.ne.itemptrans)then
            write(amessage,381) trim(facfile)
381         format(' According to the data presented in file ',a,' spatial ', &
            'interpolation is based on logarithms of values pertaining ', &
            'to some pilot points, and on native values for others. For any ', &
            'parameter type, it must be all of one or the other when using FAC2MF2K.')
            go to 9890
          end if
        end if
        do i=1,na
          iused(ipt(i))=1
        end do
      end do
384   continue
      rewind(unit=facunit,iostat=ierr)
      if(ierr.ne.0)then
        write(amessage,385) trim(facfile)
385     format(' Cannot rewind interpolation factor file ',a)
        go to 9890
      end if

      if(itrans.eq.1)then
        do irep=1,numparamreplace
          if(replace(irep)%transformtype.ne.1)then
            write(amessage,386) trim(facfile),trim(repfile)
386         format(' According to the information in file ',a,', spatial ', &
            'interpolation is based on the logarithm of parameter values rather ', &
            'than on native parameter values. Accordingly, TRANSFORMTYPE in the ', &
            'replacement file ',a,' must also be "log" for all parameters ', &
            'interpolated from these pilot points.')
            go to 9890
          end if
        end do
      end if

! -- The number of LPF and RCH replacements is now evaluated.

      numrchpar=0
      numlpfpar=0
      do irep=1,numparamreplace
        if(replace(irep)%replacetype.eq.'RCH')then
          numrchpar=numrchpar+1
        else
          numlpfpar=numlpfpar+1
        end if
      end do

! -- Information is now altered in the MODFLOW 2000 LPF file. First the existing LPF file
!    is copied to a temporary file.

       if(numlpfpar.eq.0) go to 4000

       inquire(file=lpffile,exist=lexist)
       if(.not.lexist)then
         write(amessage,400) trim(lpffile)
400      format(' Cannot find MODFLOW 2000 LPF input file ',a)
         go to 9890
       end if
       call addquote(lpffile,afile)
       call system('copy '//trim(afile)//' t###.lpf > nul')
       tunit=nextunit()
       open(unit=tunit,file='t###.lpf',status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,410)
410      format(' Cannot open temporary work file t###.lpf')
         go to 9890
       end if
405    read(tunit,'(a)',err=9250,end=9300) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 405
       call linesplit(ifail,3)
       if(ifail.ne.0)then
         write(amessage,406) trim(lpffile)
406      format(' Insufficient entries on first data line of file ',a)
         go to 9890
       end if
       nplpf=char2int(ifail,3)
       if(ifail.ne.0)then
         write(amessage,420) trim(lpffile)
420      format(' Error reading value of NPLPF from MODFLOW 2000 LPF ', &
         'input file ',a)
         go to 9890
       end if
       if(nplpf.gt.0)then
         allocate(nnclu(nplpf),stat=ierr)
         if(ierr.ne.0)then
           write(amessage,250)
           go to 9890
         end if
       else
         write(amessage,425) trim(lpffile)
425      format(' No parameters are defined in MODFLOW 2000 LPF input file ',a)
         go to 9890
       end if
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !laytyp
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !layavg
       read(tunit,*,err=9250,end=9300) (rmdat(i),i=1,nlay)    !chani
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !layvka
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !laywet


       iwdflg=0
       do i=1,nlay
         if(mdat(i).ne.0) iwdflg=iwdflg+1
       end do
       if(iwdflg.ne.0) then
         read(tunit,'(a)',err=9250,end=9300) cline
       end if
       mplpf=0
       do ip=1,nplpf
         read(tunit,'(a)',err=9250,end=9300) cline
         call linesplit(ifail,4)
         if(ifail.ne.0)then
           write(amessage,430) trim(lpffile)
430        format(' Error reading parameter data from file ',a,': insufficient ', &
           'entries on parameter data line.')
           go to 9890
         end if
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         partyp=cline(left_word(2):right_word(2))
         call casetrans(partyp,'hi')
         parval=char2real(ifail,3)
         if(ifail.ne.0)then
           write(amessage,440) trim(parnam),trim(lpffile)
440        format(' Cannot read value for parameter "',a,'" from MODFLOW 2000 ', &
           'LPF input file ',a)
           go to 9890
         end if
         nclu=char2int(ifail,4)
         if(ifail.ne.0)then
           write(amessage,450) trim(parnam),trim(lpffile)
450        format(' Cannot read number of clusters for parameter "',a,'" from ', &
           'MODFLOW 2000 input file ',a)
           go to 9890
         end if
         jreplace=0
         do irep=1,numparamreplace
           if(replace(irep)%replacetype.eq.'RCH') cycle
           do j=1,replace(irep)%numout
             if(replace(irep)%outparamname(j).eq.parnam)then
               if(replace(irep)%replacetype.ne.partyp)then
                 write(amessage,455) trim(parnam),trim(lpffile),trim(repfile), &
                 trim(partyp),trim(lpffile),trim(repfile)
455              format(' Parameter "',a,'" from file ',a,' has been marked ', &
                 'for removal in file ',a,'; however it is of a different type ("',a, &
                 '") in file ',a,' from that expected from the information in file ',a)
                 go to 9890
               end if
               jreplace=1
             end if
           end do
         end do

         do irep=1,numparamreplace
           if(partyp.eq.replace(irep)%replacetype) then
             found(irep)=1
           end if
         end do

         if(jreplace.eq.0)then
           do icl=1,nclu
             read(tunit,'(a)',err=9250,end=9300) cline
           end do
         else
           mclu=0
           do icl=1,nclu
             read(tunit,'(a)',err=9250,end=9300) cline
             call linesplit(ifail,4)
             if(ifail.ne.0)then
               write(amessage,470) trim(parnam),trim(lpffile)
470            format(' Insufficient entries on one of the cluster lines for parameter "', &
               a,'" in file ',a)
               go to 9890
             end if
             layer=char2int(ifail,1)
             if(ifail.ne.0)then
               write(amessage,480) trim(parnam),trim(lpffile)
480            format(' Cannot read layer number from one of the cluster lines for ', &
               'parameter "',a,'" in file ',a)
               go to 9890
             end if
             do irep=1,numparamreplace
               if(replace(irep)%replacetype.eq.'RCH') cycle
               do j=1,replace(irep)%numout
                 if(replace(irep)%outparamname(j).eq.parnam)then
                   if((layer.ge.replace(irep)%outparamlay1(j)).and. &
                      (layer.le.replace(irep)%outparamlay2(j)))then
                        mclu=mclu+1
                        jjflag=0
                        do kk=1,replace(irep)%numin
                          if((layer.ge.replace(irep)%newparamlay1(kk)).and. &
                             (layer.le.replace(irep)%newparamlay2(kk))) jjflag=1
                        end do
                        if(jjflag.eq.0)then
                          write(amessage,481) trim(parnam)
481                       format(' Parameter ',a,' has been identified for ', &
                          'removal in parameter replacement file; however it ', &
                          'is associated with more layers in the model than the ' &
                          'parameter replacing it in the same replacement set.')
                          go to 9890
                        end if
                        go to 490
                   end if
                 end if
               end do
             end do
490          continue
           end do
           nclu=nclu-mclu
           if(nclu.eq.0)then
             mplpf=mplpf+1
           end if
           nnclu(ip)=nclu
         end if
       end do

       do irep=1,numparamreplace
         if(found(irep).eq.0)then
           if(replace(irep)%replacetype.ne.'RCH')then
           call num2char(irep,anum)
           write(amessage,505) trim(lpffile),trim(replace(irep)%replacetype), &
           trim(anum),trim(repfile)
505        format(' No parameters were found in LPF input file ',a,' of type "', &
           a,'"; parameters of this type must be defined in this file if parameters ', &
           'defined in replacement set ',a,' from file ',a,' are to be inserted into ', &
           'the MODFLOW 2000 dataset for estimation.')
           go to 9890
           end if
         end if
       end do


       newnplpf=nplpf-mplpf+count(iused.ne.0)*numlpfpar


! -- The temporary copy of the LPF file is now re-wound and the new LPF file written.

       rewind(unit=tunit,iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,510)
510      format(' Cannot rewind temporary file t###.lpf')
         go to 9890
       end if

       ncopylpf=1
       lpfunit=nextunit()
       open(unit=lpfunit,file=lpffile)
520    read(tunit,'(a)',err=9250,end=9300) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
         write(lpfunit,'(a)') trim(cline)
         go to 520
       end if
       call linesplit(ifail,3)
       call num2char(newnplpf,anum)
       write(lpfunit,530) cline(1:right_word(2)),trim(anum)
530    format(a,2x,a)
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !laytyp
       write(lpfunit,540) (mdat(i),i=1,nlay)
540    format(80(i3))
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !layavg
       write(lpfunit,540) (mdat(i),i=1,nlay)
       read(tunit,*,err=9250,end=9300) (rmdat(i),i=1,nlay)    !chani
       write(lpfunit,541) (rmdat(i),i=1,nlay)
541    format(80(1x,1pg13.6))       
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !layvka
       write(lpfunit,540) (mdat(i),i=1,nlay)
       read(tunit,*,err=9250,end=9300) (mdat(i),i=1,nlay)    !laywet
       write(lpfunit,540) (mdat(i),i=1,nlay)
       if(iwdflg.ne.0) then
         read(tunit,'(a)',err=9250,end=9300) cline
         write(lpfunit,'(a)') trim(cline)
       end if

! -- Parameters are removed as directed.

       do ip=1,nplpf
         read(tunit,'(a)',err=9250,end=9300) cline
         call linesplit(ifail,4)
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         partyp=cline(left_word(2):right_word(2))
         call casetrans(partyp,'hi')
         parval=char2real(ifail,3)
         nclu=char2int(ifail,4)
         jreplace=0
         iirep=0
         do irep=1,numparamreplace
           do j=1,replace(irep)%numout
             iirep=iirep+1
             if(replace(irep)%replacetype.ne.'RCH') then
               if(replace(irep)%outparamname(j).eq.parnam)then
                 jreplace=1
                 removed_par(iirep)=1
               end if
             end if
           end do
         end do
         if(jreplace.eq.0)then
           write(lpfunit,'(a)') trim(cline)
           do icl=1,nclu
             read(tunit,'(a)',err=9250,end=9300) cline
             write(lpfunit,'(a)') trim(cline)
           end do
         else
           if(nnclu(ip).gt.0)then
             write(lpfunit,560) cline(1:right_word(3)), nnclu(ip)
560          format(a,2x,i5)
           end if
           do icl=1,nclu
             read(tunit,'(a)',err=9250,end=9300) cline
             call linesplit(ifail,4)
             layer=char2int(ifail,1)
             do irep=1,numparamreplace
               if(replace(irep)%replacetype.eq.'RCH') cycle
               do j=1,replace(irep)%numout
                 if(replace(irep)%outparamname(j).eq.parnam)then
                   if((layer.ge.replace(irep)%outparamlay1(j)).and. &
                      (layer.le.replace(irep)%outparamlay2(j)))then
                        go to 570
                   end if
                 end if
               end do
             end do
             write(lpfunit,'(a)') trim(cline)
570          continue
           end do
         end if
       end do

! -- New parameters are added.

       do irep=1,numparamreplace
         if(replace(irep)%replacetype.eq.'RCH') cycle
         llayer=0          !llayer is an array
         nclu=0
         do i=1,nlay
           do j=1,replace(irep)%numin
             if((i.ge.replace(irep)%newparamlay1(j)).and.  &
                (i.le.replace(irep)%newparamlay2(j)))then
               llayer(i)=1
             end if
           end do
         end do
         nclu=count(llayer.ne.0)
         if(nclu.eq.0)then
           call num2char(irep,anum)
           write(amessage,580) trim(anum),trim(repfile)
580        format(' The layers defined using the "NEWPARAMLAYS" keyword for ', &
           'replacement set ',a,' in file ',a,' are not in accord with the number ', &
           'of model layers.')
           go to 9890
         end if
         ipp=0
         do ipoint=1,num_pilot_points
           if(iused(ipoint).eq.0) cycle
           ipp=ipp+1
           parnam=trim(replace(irep)%newparamprefix)//trim(pilot_point_id(ipoint))
           parnam=adjustl(parnam)
           if(len_trim(parnam).gt.10) parnam=parnam(1:10)
           partyp=replace(irep)%replacetype
           parval=pilot_point_val(ipoint)
           write(lpfunit,600) trim(parnam),trim(partyp),parval,nclu
600        format(a,t12,a,t24,1pg13.6,2x,i3)
           inum=0
           izone=pilot_point_zone(ipoint)
           do j=1,ipoint
             if(iused(j).eq.0) cycle
             if(pilot_point_zone(j).eq.izone)inum=inum+1
           end do
           call num2char(inum,anum)
           mltarr=trim(abase)//trim(anum)
           zonarr=abase
           do ilay=1,nlay
             if(llayer(ilay).ne.0)then
               write(lpfunit,620) ilay,trim(mltarr),trim(zonarr),pilot_point_zone(ipoint)
620            format(i3,t12,a,t24,a,t36,i4)
             end if
           end do
         end do
       end do

! -- The remainder of the LPF input file is now written.

      do
        read(tunit,'(a)',err=9250,end=650) ccline
        write(lpfunit,'(a)') trim(ccline)
      end do

650   continue
      close(unit=tunit)
      close(unit=lpfunit)
      bfile=trim(lpffile)//'.old'
      call addquote(bfile,afile)
      call system('copy t###.lpf '//trim(afile)//' > nul')

4000  continue
      if(numrchpar.eq.0) go to 690

! -- Information is now altered in the MODFLOW 2000 RCH file.

!    First the existing RCH file is copied to a temporary file.


       if(rchfile.eq.' ') then
         write(amessage,4350)
4350     format(' Parameters of type RCH were cited in the parameter replacement ', &
         'file; however the RCH process is not cited in the MODFLOW 2000 name ', &
         'file.')
         go to 9890
       end if
       inquire(file=rchfile,exist=lexist)
       if(.not.lexist)then
         write(amessage,4400) trim(rchfile)
4400     format(' Cannot find MODFLOW 2000 RCH input file ',a)
         go to 9890
       end if
       call addquote(rchfile,afile)
       call system('copy '//trim(afile)//' t###.rch > nul')
       tunit=nextunit()
       open(unit=tunit,file='t###.rch',status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,4410)
4410      format(' Cannot open temporary work file t###.rch')
         go to 9890
       end if
4405   read(tunit,'(a)',err=9270,end=9320) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 4405
       call casetrans(cline,'lo')
       call linesplit(ifail,2)
       if(ifail.ne.0)then
         write(amessage,4406) trim(rchfile)
4406      format(' Insufficient entries on first data line of MODFLOW 2000 ', &
          'RCH package input file ',a)
         go to 9890
       end if
       if(cline(left_word(1):right_word(1)).ne.'parameter')then
         write(amessage,4407) trim(rchfile)
4407     format(' The word "PARAMETER" is expected on the first ', &
         'data line of MODFLOW 2000 RCH package input file ',a)
         go to 9890
       end if
       nprch=char2int(ifail,2)
       if(ifail.ne.0)then
         write(amessage,4420) trim(rchfile)
4420     format(' Error reading value of NPRCH from MODFLOW 2000 RCH ', &
         'input file ',a)
         go to 9890
       end if
       if(nprch.le.0)then
         write(amessage,4425) trim(rchfile)
4425     format(' No parameters are defined in MODFLOW 2000 RCH input file ',a)
         go to 9890
       end if
       read(tunit,*,err=9270,end=9320) nrchop,irchcb
       jrepce=0
       do ip=1,nprch
         read(tunit,'(a)',err=9270,end=9320) cline
         call linesplit(ifail,4)
         if(ifail.ne.0)then
           write(amessage,430) trim(rchfile)
           go to 9890
         end if
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         partyp=cline(left_word(2):right_word(2))
         call casetrans(partyp,'hi')
         parval=char2real(ifail,3)
         if(ifail.ne.0)then
           write(amessage,4440) trim(parnam),trim(rchfile)
4440       format(' Cannot read value for parameter "',a,'" from MODFLOW 2000 ', &
           'RCH input file ',a)
           go to 9890
         end if
         nclu=char2int(ifail,4)
         if(ifail.ne.0)then
           write(amessage,450) trim(parnam),trim(rchfile)
           go to 9890
         end if
         do irep=1,numparamreplace
           do j=1,replace(irep)%numout
             if(replace(irep)%outparamname(j).eq.parnam)then
               if(replace(irep)%replacetype.ne.partyp)then
                 write(amessage,455) trim(parnam),trim(rchfile),trim(repfile), &
                 trim(partyp),trim(rchfile),trim(repfile)
                 go to 9890
               end if
               jrepce=jrepce+1
             end if
           end do
         end do

         do icl=1,nclu
           read(tunit,'(a)',err=9270,end=9320) cline
         end do
       end do

       newnprch=nprch-jrepce+count(iused.ne.0)*numrchpar


! -- The temporary copy of the RCH file is now re-wound and the new RCH file written.

       rewind(unit=tunit,iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,4510)
4510      format(' Cannot rewind temporary file t###.rch')
         go to 9890
       end if

       ncopyrch=1
       rchunit=nextunit()
       open(unit=rchunit,file=rchfile)
4520   read(tunit,'(a)',err=9270,end=9320) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
         write(rchunit,'(a)') trim(cline)
         go to 4520
       end if
       call num2char(newnprch,anum)
       write(rchunit,4530) trim(anum)
4530   format('PARAMETER',2x,a)
       read(tunit,'(a)',err=9270,end=9320) cline
       write(rchunit,'(a)') trim(cline)

! -- Parameters are removed as directed.

       do ip=1,nprch
         read(tunit,'(a)',err=9270,end=9320) cline
         call linesplit(ifail,4)
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         partyp=cline(left_word(2):right_word(2))
         call casetrans(partyp,'hi')
         parval=char2real(ifail,3)
         nclu=char2int(ifail,4)
         jreplace=0
         iirep=0
         do irep=1,numparamreplace
           do j=1,replace(irep)%numout
             iirep=iirep+1
             if(replace(irep)%replacetype.eq.'RCH')then
               if(replace(irep)%outparamname(j).eq.parnam)then
                 jreplace=1
                 removed_par(iirep)=1
               end if
             end if
           end do
         end do
         if(jreplace.eq.0)then
           write(rchunit,'(a)') trim(cline)
           do icl=1,nclu
             read(tunit,'(a)',err=9270,end=9320) cline
             write(rchunit,'(a)') trim(cline)
           end do
         else
           do icl=1,nclu
             read(tunit,'(a)',err=9270,end=9320) cline
           end do
         end if
       end do

! -- New parameters are added.

       do irep=1,numparamreplace
         if(replace(irep)%replacetype.ne.'RCH') cycle
         do j=1,replace(irep)%numin
           if((replace(irep)%newparamlay1(j).le.0).or.   &
              (replace(irep)%newparamlay2(j).gt.nper))then
              write(amessage,4565) trim(repfile)
4565          format(' NEWPARAMPERS values for at least one replacement block ', &
              'of type RCH in file ',a,' are not compatible with number of model ', &
              'stress periods.')
              go to 9890
           end if
         end do
         found(irep)=1
         ipp=0
         do ipoint=1,num_pilot_points
           if(iused(ipoint).eq.0) cycle
           ipp=ipp+1
           parnam=trim(replace(irep)%newparamprefix)//trim(pilot_point_id(ipoint))
           parnam=adjustl(parnam)
           if(len_trim(parnam).gt.10) parnam=parnam(1:10)
           partyp=replace(irep)%replacetype
           parval=pilot_point_val(ipoint)
           write(rchunit,600) trim(parnam),trim(partyp),parval,1
           inum=0
           izone=pilot_point_zone(ipoint)
           do j=1,ipoint
             if(iused(j).eq.0) cycle
             if(pilot_point_zone(j).eq.izone)inum=inum+1
           end do
           call num2char(inum,anum)
           mltarr=trim(abase)//trim(anum)
           zonarr=abase
           write(rchunit,4580) trim(mltarr),trim(zonarr),pilot_point_zone(ipoint)
4580       format(1x,a,t24,a,t40,i4)
         end do
4670     continue
       end do

! -- The remainder of the RCH input file is now written.

      do iper=1,nper
        read(tunit,*,err=9270,end=9320) inrech,inirch
        jreplace=0
        do irep=1,numparamreplace
          if(replace(irep)%replacetype.ne.'RCH') cycle
          do j=1,replace(irep)%numin
            if((replace(irep)%newparamlay1(j).le.iper).and.   &
               (replace(irep)%newparamlay2(j).ge.iper))then
               if(jreplace.ne.0)then
                 write(amessage,4675) trim(repfile)
4675             format(' At least two new RCH parameter sets have ', &
                 'been assigned to the same stress period in file ',a,  &
                 ' or overlapping NEWPARAMPERS designations have been used; ', &
                 'only one parameter set replacement per stress period is allowed.')
                 go to 9890
               end if
               inrechold=inrech
               inrech=ipp
               jreplace=irep
            end if
          end do
        end do
4680    continue
        write(rchunit,4690) inrech,inirch
4690    format(2i10)
        if(jreplace.eq.0)then
          do i=1,inrech
            read(tunit,'(a)',err=9270,end=9320) cline
            call linesplit(ifail,2)
            if(ifail.ne.0)then
              call num2char(iper,anum)
              write(amessage,4700) trim(anum),trim(rchfile)
4700          format(' Insufficient entries in parameter name line for stress ', &
              'period ',a,' in MODFLOW 2000 recharge file ',a)
              go to 9890
            end if
            parnam=cline(left_word(1):right_word(1))
            call casetrans(parnam,'hi')
            do irep=1,numparamreplace
              if(replace(irep)%replacetype.ne.'RCH') cycle
              do j=1,replace(irep)%numout
                if(replace(irep)%outparamname(j).eq.parnam)then
                  call num2char(iper,anum)
                  write(amessage,4710) trim(parnam),trim(rchfile),trim(repfile), &
                  trim(anum)
4710              format(' Parameter ',a,' of type RCH was designated for removal ', &
                  'from MODFLOW 2000 recharge file ',a,  &
                  ' in parameter replacement file ',a,'; however it is still used ', &
                  'in stress period ',a,' even after parameter replacement has ', &
                  'taken place in this file.')
                  go to 9890
                end if
              end do
            end do
            write(rchunit,'(a)') trim(cline)
          end do
        else
          do i=1,inrechold
            read(tunit,'(a)',err=9270,end=9320) cline
          end do
          do ipoint=1,num_pilot_points
            if(iused(ipoint).eq.0) cycle
            parnam=trim(replace(jreplace)%newparamprefix)//trim(pilot_point_id(ipoint))
            parnam=adjustl(parnam)
            if(len_trim(parnam).gt.10) parnam=parnam(1:10)
            write(rchunit,4720) trim(parnam),-1
4720        format(1x,a,t20,i5)
          end do
        end if
        if((nrchop.eq.2).and.(inirch.gt.0))then
          call u2dint(ifail,1,nrow,ncol,ia,tunit,rchoutunit,rchunit,iprn)
        end if
      end do

! -- A check is made that all recharge parameters have been used.

      do irep=1,numparamreplace
        if(replace(irep)%replacetype.ne.'RCH') cycle
        if(found(irep).eq.0)then
          write(amessage,4740) trim(repfile)
4740      format(' At least one recharge replacement set provided in file ',a,  &
          ' could not be included in the input dataset because it did not ', &
          'pertain to any stress period undertaken by the model.')
          go to 9890
        end if
      end do

      close(unit=tunit)
      close(unit=rchunit)
      bfile=trim(rchfile)//'.old'
      call addquote(bfile,afile)
      call system('copy t###.rch '//trim(afile)//' > nul')

690   continue

! -- Attemtion must now be turned to the zone file.
! -- First it is copied to a temporary working file.

       inquire(file=zonfile,exist=lexist)
       if(.not.lexist)then
         write(amessage,700) trim(zonfile)
700      format(' Cannot find MODFLOW 2000 zone file ',a)
         go to 9890
       end if
       call addquote(zonfile,afile)
       call system('copy '//trim(afile)//' t###.zon > nul')
       tunit=nextunit()
       open(unit=tunit,file='t###.zon',status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,710)
710      format(' Cannot open temporary work file t###.zon')
         go to 9890
       end if

! -- The temporary zone file is read a first time in order to ensure that the
!    names requested for removal are indeed present within that file.

740    read(tunit,'(a)',err=9500,end=9550) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 740
       call linesplit(ifail,1)
       nzn=char2int(ifail,1)
       if(ifail.ne.0)then
         write(amessage,750) trim(zonfile)
750      format(' Cannot read NZN variable from MODFLOW 2000 zone file ',a)
         go to 9890
       end if
       do iz=1,nzn
         read(tunit,'(a)',err=9600,end=9550) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 9600
         zonnam=cline(left_word(1):right_word(1))
         call casetrans(zonnam,'hi')
         irflag=0
         do i=1,numzoneremove
           if(zonnam(1:10).eq.removezone(i)(1:10))then
             removed_zon(i)=1
             irflag=1
           end if
         end do
         atemp=abase
         if(zonnam.eq.atemp)then
           if(irflag.eq.0)then
             write(amessage,755) trim(abase),trim(zonfile)
755          format(' Your choice of name base ("',a,'") for new multiplier and ', &
             'zonation arrays will not result in a unique zone array name in file ',a)
             go to 9890
           end if
         end if
         iwrite=0
         call u2dint(ifail,iwrite,nrow,ncol,ia,tunit,modzonunit,zonunit,iprn)
         if(ifail.eq.1) then
           go to 9600
         else if(ifail.eq.2) then
           go to 9550
         end if
       end do
       newnzn=nzn-count(removed_zon.ne.0)+1

       rewind(unit=tunit,iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,760)
760      format(' Cannot rewind temporary file t###.zon')
         go to 9890
       end if

! -- The new zone file is now written - first the existing files.

       ncopyzon=1
       zonunit=nextunit()
       open(unit=zonunit,file=zonfile)
780    read(tunit,'(a)',err=9500,end=9550) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
         write(zonunit,'(a)') trim(cline)
         go to 780
       end if
       write(zonunit,'(i5)') newnzn
       do iz=1,nzn
         read(tunit,'(a)',err=9600,end=9550) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 9600
         zonnam=cline(left_word(1):right_word(1))
         call casetrans(zonnam,'hi')
         iwrite=1
         do i=1,numzoneremove
           if(zonnam(1:10).eq.removezone(i)(1:10))then
             iwrite=0
           end if
         end do
         if(iwrite.ne.0)then
           write(zonunit,'(a)') trim(cline)
         end if
         call u2dint(ifail,iwrite,nrow,ncol,ia,tunit,modzonunit,zonunit,iprn)
       end do

! -- Now the new zone is added.

       write(zonunit,'(a)') trim(abase)
       write(zonunit,800) iprn
800    format('INTERNAL  1  ''(FREE)'' ',i5)
       do irow=1,nrow
         write(zonunit,810) (intarray(icol,irow),icol=1,ncol)
810      format(20(i5))
       end do
       

       close(unit=tunit)
       close(unit=zonunit)
       bfile=trim(zonfile)//'.old'
       call addquote(bfile,afile)
       call system('copy t###.zon '//trim(afile)//' > nul')


! -- Attemtion must now be turned to the multiplier file.
! -- First it is copied to a temporary working file.

       inquire(file=mulfile,exist=lexist)
       if(.not.lexist)then
         write(amessage,830) trim(mulfile)
830      format(' Cannot find MODFLOW 2000 multiplier file ',a)
         go to 9890
       end if
       call addquote(mulfile,afile)
       call system('copy '//trim(afile)//' t###.mul > nul')
       tunit=nextunit()
       open(unit=tunit,file='t###.mul',status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,840)
840      format(' Cannot open temporary work file t###.mul')
         go to 9890
       end if

! -- The temporary multiplier file is read a first time in order to ensure that the
!    names requested for removal are indeed present within that file.

850    read(tunit,'(a)',err=9700,end=9750) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 850
       call linesplit(ifail,1)
       nml=char2int(ifail,1)
       if(ifail.ne.0)then
         write(amessage,860) trim(mulfile)
860      format(' Cannot read NML variable from MODFLOW 2000 multiplier file ',a)
         go to 9890
       end if
       do im=1,nml
         read(tunit,'(a)',err=9800,end=9750) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 9800
         mltnam=cline(left_word(1):right_word(1))
         call casetrans(mltnam,'hi')
         irflag=0
         do i=1,nummultremove
           if(mltnam(1:10).eq.removemult(i)(1:10))then
             removed_mul(i)=1
             irflag=1
           end if
         end do
         lt=len_trim(abase)
         do i=1,lt
           if(mltnam(i:i).ne.abase(i:i)) go to 880
         end do
         atemp=mltnam(lt+1:)
         if(atemp.eq.' ') go to 880
         call char2num(ifail,atemp,itemp)
         if(ifail.ne.0) go to 880
         if(itemp.le.num_pilot_points) then
           if(irflag.eq.0)then
             write(amessage,865) trim(abase),trim(mulfile)
865          format(' Your choice of name base ("',a,'") for new multiplier and ', &
             'zonation arrays will not result in unique multiplier array names in file ',a)
             go to 9890
           end if
         end if
880      continue
         call linesplit(ifail,2)
         if(ifail.eq.0)then
           afunction=cline(left_word(2):right_word(2))
           call casetrans(afunction,'hi')
           if(afunction.eq.'FUNCTION')then
             read(tunit,'(a)',err=9800,end=9750) cline
             go to 885
           end if
         end if
         iwrite=0
         call u2drel(ifail,iwrite,nrow,ncol,ra,tunit,modmulunit,mulunit,iprn)
         if(ifail.eq.1) then
           go to 9800
         else if(ifail.eq.2) then
           go to 9750
         end if
885      continue
       end do

       rewind(unit=tunit,iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,890)
890      format(' Cannot rewind temporary file t###.mul')
         go to 9890
       end if


! -- Next we must work out how many new multiplier arrays are needed.

       newmul=0
       do i=1,num_pilot_points
         itemp=pilot_point_zone(i)
         if(i.gt.1)then
           do j=1,i-1
             if(pilot_point_zone(j).eq.itemp) go to 900
           end do
         end if
!         if(iused(i).eq.0) go to 900
         itemp1=count((pilot_point_zone.eq.itemp).and.(iused.ne.0))
         if(itemp1.gt.newmul)newmul=itemp1
900    continue
       end do

       newnml=nml-count(removed_mul.ne.0)+newmul

! -- The new multiplier file is now written - first the existing arrays.

       ncopymul=1
       mulunit=nextunit()
       open(unit=mulunit,file=mulfile)
910    read(tunit,'(a)',err=9700,end=9750) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
         write(mulunit,'(a)') trim(cline)
         go to 910
       end if
       write(mulunit,'(i5)') newnml
       do im=1,nml
         ifunct=0
         read(tunit,'(a)',err=9800,end=9750) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 9800
         mltnam=cline(left_word(1):right_word(1))
         call casetrans(mltnam,'hi')
         call linesplit(ifail,2)
         if(ifail.eq.0)then
           afunction=cline(left_word(2):right_word(2))
           call casetrans(afunction,'hi')
           if(afunction.eq.'FUNCTION') ifunct=1
         end if
         iwrite=1
         do i=1,nummultremove
           if(mltnam(1:10).eq.removemult(i)(1:10))then
             iwrite=0
           end if
         end do
         if(ifunct.eq.1)then
           if(iwrite.eq.0)then
             read(tunit,'(a)',err=9800,end=9750) cline
           else
             write(mulunit,'(a)') trim(cline)
             read(tunit,'(a)',err=9800,end=9750) cline
             write(mulunit,'(a)') trim(cline)
           end if
         else
           if(iwrite.ne.0)then
             write(mulunit,'(a)') trim(cline)
           end if
           call u2drel(ifail,iwrite,nrow,ncol,ra,tunit,modmulunit,mulunit,iprn)
         end if
       end do

! -- Now the new multiplier arrays are added; first memory is allocated.

       allocate(amul(ncol,nrow,newmul),stat=ierr)
       if(ierr.ne.0) then
         write(amessage,250)
         go to 9890
       end if
       amul=0.0      !amul is an array

! -- The weights array is now read line by line and weights are allocated to the
!    respective array.

       if(facformat.eq.'f')then
         read(facunit,'(a)',err=9150,end=9200) afile
         read(facunit,'(a)',err=9150,end=9200) afile
         read(facunit,*,err=9150,end=9200) ncol,nrow
         read(facunit,*,err=9150,end=9200) nnppt
         do iipt=1,nnppt
           read(facunit,'(a)',err=9000,end=9100) wpoints(iipt)
          end do
       else
         read(facunit,err=9150,end=9200) afile
         read(facunit,err=9150,end=9200) afile
         read(facunit,err=9150,end=9200) ncol,nrow
         read(facunit,err=9150,end=9200) nnppt
         do iipt=1,nnppt
           read(facunit,err=9000,end=9100) wpoints(iipt)
         end do
       end if

       afile=facfile
       do
         if(facformat.eq.'f')then
           read(facunit,*,err=9850,end=950) icellno,itemptrans, &
           na,rtemp,((ipt(i),wt(i)),i=1,na)
         else
           read(facunit,err=9850,end=950)   icellno,itemptrans, &
           na,rtemp,((ipt(i),wt(i)),i=1,na)
         end if
         do i=1,na
           irow=(icellno-1)/ncol+1
           icol=icellno-((irow-1)*ncol)
           itemp=pilot_point_zone(ipt(i))
           iarr=0
           do j=1,ipt(i)
             if(iused(j).eq.0) cycle
             if(pilot_point_zone(j).eq.itemp) iarr=iarr+1
           end do
           amul(icol,irow,iarr)=wt(i)
         end do
       end do

950    continue
       close(unit=facunit)

! -- The new weights arrays are now written.

       do i=1,newmul
         call num2char(i,atemp)
         atemp=trim(abase)//trim(atemp)
         write(mulunit,'(a)') trim(atemp)
         write(mulunit,960) iprn
960      format('INTERNAL  1.0  ''(FREE)'' ',i5)
         do irow=1,nrow
           write(mulunit,970) (amul(icol,irow,i),icol=1,ncol)
970        format(7(1x,1pg14.7))
         end do
       end do

       close(unit=tunit)
       close(unit=mulunit)
       bfile=trim(mulfile)//'.old'
       call addquote(bfile,afile)
       call system('copy t###.mul '//trim(afile)//' > nul')


! -- If a MODFLOW 2000 sensitivity file is present, then that is altered.

       if(senfile.eq.' ') go to 3000

! -- The sensitivity file is read a first time and a count is made of the number
!    of parameters


       inquire(file=senfile,exist=lexist)
       if(.not.lexist)then
         write(amessage,1010) trim(senfile)
1010      format(' Cannot find MODFLOW 2000 sensitivity input file ',a)
         go to 9890
       end if
       call addquote(senfile,afile)
       call system('copy '//trim(afile)//' t###.sen > nul')
       tunit=nextunit()
       open(unit=tunit,file='t###.sen',status='old',iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,1020)
1020      format(' Cannot open temporary work file t###.sen')
         go to 9890
       end if
1030   read(tunit,'(a)',err=7000,end=7050) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 1030
       call linesplit(ifail,4)
       if(ifail.ne.0)then
         write(amessage,1050) trim(senfile)
         go to 9890
       end if
       nplist=char2int(ifail,1)
       if(ifail.ne.0)then
         write(amessage,1050) trim(senfile)
1050    format(' Error reading value of NPLIST from MODFLOW 2000 sensitivity ', &
         'input file ',a)
         go to 9890
       end if
       read(tunit,'(a)',err=7000,end=7050) cline
       iparout=0
       do ip=1,nplist
         read(tunit,'(a)',err=7100,end=7050) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 7100
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         do irep=1,numparamreplace
           do j=1,replace(irep)%numout
             if(replace(irep)%outparamname(j).eq.parnam)then
               iparout=iparout+1
               go to 1055
             end if
           end do
         end do
1055     continue
       end do
       newnplist=nplist-iparout+count(iused.ne.0)*numparamreplace

! -- The temporary copy of the sensitivity file is now re-wound and the new
!    sensitivity file written.

       rewind(unit=tunit,iostat=ierr)
       if(ierr.ne.0)then
         write(amessage,1080)
1080     format(' Cannot rewind temporary file t###.sen')
         go to 9890
       end if

       ncopysen=1
       senunit=nextunit()
       open(unit=senunit,file=senfile)
1090   read(tunit,'(a)',err=7000,end=7050) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
         write(senunit,'(a)') trim(cline)
         go to 1090
       end if
       call linesplit(ifail,4)
       call num2char(newnplist,anum)
       write(senunit,1100) trim(anum),cline(left_word(2):right_word(3)),trim(anum)
1100   format(a,2x,a,2x,a)
       read(tunit,'(a)',err=7000,end=7050) cline
       write(senunit,'(a)') trim(cline)
       do ip=1,nplist
         read(tunit,'(a)',err=7100,end=7050) cline
         call linesplit(ifail,1)
         if(ifail.ne.0) go to 7100
         parnam=cline(left_word(1):right_word(1))
         call casetrans(parnam,'hi')
         iirep=0
         do irep=1,numparamreplace
           do j=1,replace(irep)%numout
             iirep=iirep+1
             if(replace(irep)%outparamname(j).eq.parnam)then
               if(removed_par(iirep).ne.0)then
                 go to 1105
               end if
             end if
           end do
         end do
         write(senunit,'(a)') trim(cline)
1105   continue
       end do

! -- The extra, pilot-point-based parameters are now added.

       do irep=1,numparamreplace
         do i=1,num_pilot_points
           if(iused(i).eq.0) cycle
           atempp=trim(replace(irep)%newparamprefix)//trim(pilot_point_id(i))
           parnam=atempp(1:10)
           parnam=adjustl(parnam)
           call casetrans(parnam,'hi')
           isens=1
           ln=replace(irep)%transformtype
           b=pilot_point_val(i)
           bl=replace(irep)%minpval
           bu=replace(irep)%maxpval
           bscal=min(abs(bl),abs(bu))
           bscal=bscal*0.1
           write(senunit,1120) trim(parnam),isens,ln,b,bl,bu,bscal
1120       format(1x,a,t13,1x,i2,1x,i2,4(1x,1pg14.7))
         end do
       end do
       close(unit=senunit)
       close(unit=tunit)
       bfile=trim(senfile)//'.old'
       call addquote(bfile,afile)
       call system('copy t###.sen '//trim(afile)//' > nul')

! -- If the PES file is present it is now read in order to check whether
!    any prior information in that file cites parameters which have been removed.

       if(pesfile.eq.' ') go to 1210
       pesunit=nextunit()
       open(unit=pesunit,file=pesfile,status='old',iostat=ierr)
       if(ierr.ne.0) go to 1200
1160   read(pesunit,'(a)',err=1200,end=1200) cline
       if((cline.eq.' ').or.(cline(1:1).eq.'#')) go to 1160
!       read(pesunit,*,err=1200,end=1200) maxiter,maxchange,tol,sosc
       read(pesunit,*,err=1200,end=1200) ibeflg,iycflg,iostar,nopt,nfit,sosr,rmar, &
       rmarm,iap
       read(pesunit,*,err=1200,end=1200) iprcov,iprint,lprint
       read(pesunit,*,err=1200,end=1200) csa, fconv,lastx
       read(pesunit,*,err=1200,end=1200) npng,ipr,mpr
       if((ipr.eq.0 ).and.(mpr.eq.0)) then
         close(unit=pesunit)
         go to 1200
       end if
       if(npng.gt.0)then
         read(pesunit,*,err=1200,end=1200) (rtemp,i=1,npng)
       end if
       if(ipr.gt.0)then
         do i=1,ipr
           read(pesunit,'(a)',err=1200,end=1200) cline
           call linesplit(ifail,1)
           if(ifail.ne.0) go to 1200
           parnam=cline(left_word(1):right_word(1))
           call casetrans(parnam,'hi')
           iirep=0
           do irep=1,numparamreplace
             do j=1,replace(irep)%numout
               iirep=iirep+1
               if(parnam.eq.replace(irep)%outparamname(j))then
                 if(removed_par(iirep).eq.1)go to 7300
               end if
             end do
           end do
         end do
         read(pesunit,*,err=1200,end=1200) iwtp
         do i=1,ipr
           read(pesunit,*,err=1200,end=1200) (rtemp,j=1,ipr)
         end do
       end if
       if(mpr.gt.0)then
         do i=1,mpr
           read(pesunit,'(a)',err=1200,end=1200) cline
           call casetrans(cline,'hi')
           ie=index(cline,'=')
           if(ie.eq.0) cycle
           iirep=0
           do irep=1,numparamreplace
             do j=1,replace(irep)%numout
               iirep=iirep+1
               if(removed_par(iirep).eq.0) go to 1180
               ii=index(cline,trim(replace(irep)%outparamname(j)))
               if(ii.ne.0)then
                 if(ii.lt.ie) go to 1180
                 atempp=trim(replace(irep)%outparamname(j))
                 lt=len_trim(atempp)
                 jj=ii+lt
                 if((cline(ii-1:ii-1).eq.'+').or.(cline(ii-1:ii-1).eq.'-').or. &
                    (cline(ii-1:ii-1).eq.'=').or.(cline(ii-1:ii-1).eq.' '))then
                    if((cline(jj:jj).eq.' ').or.(cline(jj:jj).eq.'-').or. &
                       (cline(jj:jj).eq.'+').or.(cline(jj:jj+3).eq.'STAT')) go to 7300
                 end if
               end if
1180         continue
             end do
           end do
         end do
       end if
1200   continue
       close(unit=pesunit,iostat=ierr)
1210   continue

3000   continue


! -- The ASP input file is now written.


        numarray=0
        do irep=1,numparamreplace
          if(replace(irep)%minival.lt.-1.0e35)then
            numarray=numarray+1
          end if
          if(replace(irep)%maxival.lt.-1.0e35)then
            numarray=numarray+1
          end if
        end do

1299    continue
        if((aao.eq.'o').or.(aao.eq.'x'))then
          if(aao.eq.'x') go to 1298
          write(maxminunit,1280) ipestint,iinterp
1280      format(1x,2i5,t40,'IPESTINT  INTERP')
          write(maxminunit,1290) nostop,hdrybot,limop,minthick
1290      format(1x,3i5,2x,1pg13.6,t40,'NOSTOP  HDRYBOT   LIMOP  MINTHICK')
1298      continue
          write(maxminunit,1300) numparamreplace,numarray
1300      format(1x,i5,1x,i5,t40,'NUMPARAMREPLACE   NUMARRAY')
          iarr=0
          do irep=1,numparamreplace
            mdat=0                !mdat is an array
            do i=1,replace(irep)%numin
              if(replace(irep)%replacetype.ne.'RCH')then
                do j=1,nlay
                  if((j.ge.replace(irep)%newparamlay1(i)).and.  &
                     (j.le.replace(irep)%newparamlay2(i))) mdat(j)=1
                end do
              else
                do j=1,nper
                  if((j.ge.replace(irep)%newparamlay1(i)).and.  &
                     (j.le.replace(irep)%newparamlay2(i))) mdat(j)=1
                end do
              end if
            end do
            if(replace(irep)%minival.lt.-1.0e35)then
              iarr=iarr+1
              call num2char(iarr,anum)
              arrmin='ARRAY'//trim(anum)
            else
              call num2char(replace(irep)%minival,arrmin)
            end if
            if(replace(irep)%maxival.lt.-1.0e35)then
              iarr=iarr+1
              call num2char(iarr,anum)
              arrmax='ARRAY'//trim(anum)
            else
              call num2char(replace(irep)%maxival,arrmax)
            end if
            if(replace(irep).replacetype.ne.'RCH')then
              write(maxminunit,1310) replace(irep)%replacetype, &
              itrans, trim(arrmin), trim(arrmax), &
              (mdat(i),i=1,nlay)
1310          format(1x,a,t10,i3,t15,a,t30,a,t45,1000(i3))
            else
              write(maxminunit,1310) replace(irep)%replacetype, &
              itrans, trim(arrmin), trim(arrmax), &
              (mdat(i),i=1,nper)
            end if
          end do

          if(numarray.gt.0)then
            iarr=0
            do irep=1,numparamreplace
              if(replace(irep)%minival.lt.-1.0e35)then
                iarr=iarr+1
                call num2char(iarr,anum)
                arrmin='ARRAY'//trim(anum)
                write(maxminunit,1320) trim(arrmin)
1320            format(1x,a)
                do irow=1,nrow
                  write(maxminunit,1330) &
                  (replace(irep)%miniarray(icol,irow),icol=1,ncol)
1330              format(7(1x,1pg14.7))
                end do
              end if
              if(replace(irep)%maxival.lt.-1.0e35)then
                iarr=iarr+1
                call num2char(iarr,anum)
                arrmax='ARRAY'//trim(anum)
                write(maxminunit,1320) trim(arrmax)
                do irow=1,nrow
                  write(maxminunit,1330) &
                  (replace(irep)%maxiarray(icol,irow),icol=1,ncol)
                end do
              end if
            end do
          end if
          close(unit=maxminunit)
        else
          call addquote(maxminfile,afile)
          call system('copy '//trim(afile)//' t###.asp > nul')
          tunit=nextunit()
          open(unit=tunit,file='t###.asp',status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,1340)
1340        format(' Cannot open temporary work file t###.asp')
            go to 9890
          end if
1341      read(tunit,'(a)',err=9650,end=9670) cline
          if(cline.eq.' ') go to 1341
          if(cline(1:1).eq.'#') go to 1341
          call linesplit(ifail,2)
          if(ifail.ne.0) go to 9650
          ipestint=char2int(ifail,1)
          if(ifail.ne.0) go to 9650
          iinterp=char2int(ifail,2)
          if(ifail.ne.0) go to 9650
          read(tunit,*,err=9650,end=9670) nostop,hdrybot,limop,minthick
          if(ipestint.ne.0)then
            read(tunit,*,err=9650,end=9670) pestctlfile
          end if
          if(iinterp.ne.0) then
            read(tunit,*,err=9650,end=9670) nrp,nar
            if(nrp.eq.0)then
              nar=0
              iinterp=0
              go to 1369
            end if
            allocate(oldarrmin(nrp),oldarrmax(nrp),stat=ierr)
            if(ierr.ne.0)then
              write(amessage,250)
              go to 9890
            end if
            if(nar.ne.0)then
              allocate(oldarrayname(nar),stat=ierr)
              if(ierr.ne.0)then
                write(amessage,250)
                go to 9890
              end if
            end if
            do i=1,nrp
              read(tunit,*,err=9650,end=9670) arpt,itr,oldarrmin(i),oldarrmax(i)
              if(itr.ne.0) itrexist=1
              call casetrans(arpt,'hi')
              call casetrans(oldarrmin(i),'hi')
              call casetrans(oldarrmax(i),'hi')
            end do
            if(nar.ne.0)then
              do i=1,nar
                read(tunit,'(a)',err=9650,end=9670) oldarrayname(i)
                oldarrayname(i)=adjustl(oldarrayname(i))
                call casetrans(oldarrayname(i),'hi')
                do irow=1,nrow
                  read(tunit,*,err=9650,end=9670) (ra(icol,irow),icol=1,ncol)
                end do
              end do
            end if
          else
            nrp=0
            nar=0
          end if
1369      continue
          rewind(unit=tunit,iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,1370)
1370        format(' Error rewinding temporary file t###.asp.')
            go to 9890
          end if
          ncopyint=1
          maxminunit=nextunit()
          open(unit=maxminunit,file=maxminfile)
1371      read(tunit,'(a)',err=9650,end=9670) cline
          if((cline.eq.' ').or.(cline(1:1).eq.'#')) then
            write(maxminunit,'(a)') trim(cline)
            go to 1371
          end if
          read(tunit,*,err=9650,end=9670) nostop,hdrybot,limop,minthick
          if(ipestint.ne.0)then
            read(tunit,*,err=9650,end=9670) pestctlfile
          end if
          write(maxminunit,1280) ipestint,1
          write(maxminunit,1290) nostop,hdrybot,limop,minthick
          if(ipestint.ne.0)then
            call addquote(pestctlfile,bfile)
            write(maxminunit,'(1x,a)') trim(bfile)
          end if
          if((nrp.eq.0).and.(nar.eq.0))then
            close(unit=tunit,iostat=ierr)
            bfile=trim(maxminfile)//'.old'
            call addquote(bfile,afile)
            call system('copy t###.asp '//trim(afile)//' > nul')
            aao='x'
            go to 1299
          end if
          read(tunit,*,err=9650,end=9670) nrp,nar
          write(maxminunit,1300) numparamreplace+nrp,numarray+nar
          do i=1,nrp
            read(tunit,'(a)',err=9650,end=9670) cline
            call linesplit(ifail,1)
            call casetrans(cline(left_word(1):right_word(1)),'hi')
            write(maxminunit,'(a)') trim(cline)
          end do
          iarr=0
          do irep=1,numparamreplace
            mdat=0                !mdat is an array
            do i=1,replace(irep)%numin
              if(replace(irep)%replacetype.ne.'RCH')then
                do j=1,nlay
                  if((j.ge.replace(irep)%newparamlay1(i)).and.  &
                     (j.le.replace(irep)%newparamlay2(i))) mdat(j)=1
                end do
              else
                do j=1,nper
                  if((j.ge.replace(irep)%newparamlay1(i)).and.  &
                     (j.le.replace(irep)%newparamlay2(i))) mdat(j)=1
                end do
              end if
            end do
            if(replace(irep)%minival.lt.-1.0e35)then
1390          iarr=iarr+1
              call num2char(iarr,anum)
              arrmin='ARRAY'//trim(anum)
              do i=1,nar
                if(arrmin.eq.oldarrayname(i)) go to 1390
              end do
            else
              call num2char(replace(irep)%minival,arrmin)
            end if
            if(replace(irep)%maxival.lt.-1.0e35)then
1400          iarr=iarr+1
              call num2char(iarr,anum)
              arrmax='ARRAY'//trim(anum)
              do i=1,nar
                if(arrmax.eq.oldarrayname(i)) go to 1400
              end do
            else
              call num2char(replace(irep)%maxival,arrmax)
            end if
            if(replace(irep)%replacetype.ne.'RCH')then
              write(maxminunit,1310) replace(irep)%replacetype, &
              itrans, trim(arrmin), trim(arrmax), &
              (mdat(i),i=1,nlay)
            else
              write(maxminunit,1310) replace(irep)%replacetype, &
              itrans, trim(arrmin), trim(arrmax), &
              (mdat(i),i=1,nper)
            end if
          end do
          if(nar.ne.0)then
            do i=1,nar
              read(tunit,'(a)',err=9650,end=9670) oldarrayname(i)
              oldarrayname(i)=adjustl(oldarrayname(i))
              call casetrans(oldarrayname(i),'hi')
              write(maxminunit,1320) trim(oldarrayname(i))
              do irow=1,nrow
                read(tunit,*,err=9650,end=9670) (ra(icol,irow),icol=1,ncol)
              end do
              do irow=1,nrow
                write(maxminunit,1330) (ra(icol,irow),icol=1,ncol)
              end do
            end do
          end if
          if(numarray.gt.0)then
            iarr=0
            do irep=1,numparamreplace
              if(replace(irep)%minival.lt.-1.0e35)then
1410            iarr=iarr+1
                call num2char(iarr,anum)
                arrmin='ARRAY'//trim(anum)
                do i=1,nar
                  if(arrmin.eq.oldarrayname(i)) go to 1410
                end do
                write(maxminunit,1320) trim(arrmin)
                do irow=1,nrow
                  write(maxminunit,1330) &
                  (replace(irep)%miniarray(icol,irow),icol=1,ncol)
                end do
              end if
              if(replace(irep)%maxival.lt.-1.0e35)then
1420            iarr=iarr+1
                call num2char(iarr,anum)
                arrmax='ARRAY'//trim(anum)
                do i=1,nar
                  if(arrmax.eq.oldarrayname(i)) go to 1420
                end do
                write(maxminunit,1320) trim(arrmax)
                do irow=1,nrow
                  write(maxminunit,1330) &
                  (replace(irep)%maxiarray(icol,irow),icol=1,ncol)
                end do
              end if
            end do
          end if
          close(unit=tunit,iostat=ierr)
          close(unit=maxminunit)
          bfile=trim(maxminfile)//'.old'
          call addquote(bfile,afile)
          call system('copy t###.asp '//trim(afile)//' > nul')
        end if

! Finally, the MODFLOW-2000 name file is altered if appropriate.

        if(interpfile.eq.' ')then
          call addquote(namefile,afile)
          call system('copy '//trim(afile)//' t###.nam > nul')
          ncopynam=1
          tunit=nextunit()
          open(unit=tunit,file='t###.nam',status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,1430)
1430        format(' Cannot open temporary file ',a,'.')
            go to 9890
          end if
          nameunit=nextunit()
          open(unit=nameunit,file=namefile)
          iline=0
          afile='t###.nam'
          do
            iline=iline+1
            read(tunit,'(a)',err=9000,end=1500) cline
            write(nameunit,'(a)') trim(cline)
          end do
1500      write(nameunit,1510) trim(maxminfile)
1510      format('ASP  50  ',a)
          close(unit=tunit)
          close(unit=nameunit)
          bfile=trim(namefile)//'.old'
          call addquote(bfile,afile)
          call system('copy t###.nam '//trim(afile)//' > nul')
        end if

! -- A few messages are now written to the user to finish off.

        if(numlpfpar.ne.0)then
          write(amessage,3010) trim(lpffile)
3010      format(' - new LPF input file ',a,' written ok.')
          call write_message()
          afile=trim(lpffile)//'.old'
          write(amessage,3020) trim(afile)
3020      format(' - previous LPF input file copied to file ',a)
          call write_message()
        end if
        if(numrchpar.ne.0)then
          write(amessage,3022) trim(rchfile)
3022      format(' - new RCH input file ',a,' written ok.')
          call write_message()
          afile=trim(rchfile)//'.old'
          write(amessage,3023) trim(afile)
3023      format(' - previous RCH input file copied to file ',a)
          call write_message()
        end if
        iirep=0
        do irep=1,numparamreplace
          do j=1,replace(irep)%numout
            iirep=iirep+1
            if(removed_par(iirep).eq.0)then
              write(amessage,3050) trim(replace(irep)%outparamname(j))
3050          format('   Warning: parameter "',a,'" not found ', &
              'in MODFLOW 2000 LPF or RCH input file file (whichever is appropriate ', &
              'for the nominated parameter type).')
              call write_message()
            end if
          end do
        end do

        write(6,*)
        write(amessage,3060) trim(zonfile)
3060    format(' - new zone file ',a,' written ok.')
        call write_message()
        afile=trim(zonfile)//'.old'
        write(amessage,3070) trim(afile)
3070    format(' - previous zone file copied to file ',a)
        call write_message()
        do i=1,numzoneremove
          if(removed_zon(i).eq.0)then
            write(amessage,3080) trim(removezone(i)),trim(zonfile)
3080        format('   Warning: zone "',a,'" not found in file ',a)
            call write_message()
          end if
        end do

        write(6,*)
        write(amessage,3100) trim(mulfile)
3100    format(' - new multiplier file ',a,' written ok.')
        call write_message()
        afile=trim(mulfile)//'.old'
        write(amessage,3110) trim(afile)
3110    format(' - previous multiplier file copied to file ',a)
        call write_message()
        do i=1,nummultremove
          if(removed_mul(i).eq.0)then
            write(amessage,3120) trim(removemult(i)),trim(mulfile)
3120        format('   Warning: multiplier array "',a,'" not found in file ',a)
            call write_message()
          end if
        end do

        if(senfile.ne.' ')then
          write(6,*)
          write(amessage,3150) trim(senfile)
3150      format(' - new sensitivity file ',a,' written ok.')
          call write_message()
          afile=trim(senfile)//'.old'
          write(amessage,3160) trim(afile)
3160      format(' - previous sensitivity file copied to file ',a)
          call write_message()
        end if

        write(6,*)
        write(amessage,3170) trim(maxminfile)
3170    format(' - ASP input file ',a,' written ok.')
        call write_message()
        if((aao.eq.'a').or.(aao.eq.'x'))then
          afile=trim(maxminfile)//'.old'
          write(amessage,3180) trim(afile)
3180      format(' - previous ASP input file copied to file ',a)
          call write_message()
        end if

        if(interpfile.eq.' ')then
          write(6,*)
          write(amessage,3190) trim(namefile)
3190      format(' - MODFLOW name file ',a,' written ok.')
          call write_message()
          afile=trim(namefile)//'.old'
          write(amessage,3200) trim(afile)
3200      format(' - previous MODFLOW name file copied to file ',a)
          call write_message()
        end if

        write(6,*)
        write(amessage,3220)
3220    format(' Note: since the MODFLOW name file has been altered to include ', &
        'reference to the ASP process, ', &
        'you must use MF2KASP rather than MODFLOW-2000, unless you remove this line', &
        ' from the name file.')
        call write_message()
        if((itrans.eq.1).or.(itrexist.eq.1))then
          write(amessage,3230)
3230      format(' However as kriging interpolation in the present case is based on ', &
          'a log-variogram, use of the ASP process (and hence MF2KASP) is mandatory.')
          call write_message()
        end if

        go to 9900


7000    write(amessage,7010) trim(senfile)
7010    format(' Error encountered while reading data from early part of ', &
        'MODFLOW 2000 sensitivity input file ',a)
        go to 9890

7050    write(amessage,7060) trim(senfile)
7060    format(' Unexpected end encountered to sensitivity input file ',a)
        go to 9890

7100    write(amessage,7110) trim(senfile)
7110    format(' Error reading parameter data from sensitivity input file ',a)
        go to 9890

7300    write(amessage,7310) trim(pesfile)
7310    format(' There is prior information cited in the MODFLOW 2000 parameter ', &
        'estimation input file ',a,' which includes parameters which have just been ', &
        'removed from the parameter estimation process.')
        go to 9890


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a)
        go to 9890

9050    call num2char(iline,aline)
        write(amessage,9060) trim(aline),trim(afile)
9060    format(' Insufficient entries on line ',a,' of file ',a)
        go to 9890

9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Error reading filename from line ',a,' of file ',a)
        go to 9890

9150    write(amessage,9160) trim(facfile)
9160    format(' Error reading header information in file ',a)
        go to 9890

9200    write(amessage,9210) trim(facfile)
9210    format(' Unexpected end encountered to file ',a)
        go to 9890

9250    write(amessage,9260) trim(lpffile)
9260    format(' Error encountered while reading data from early part of ', &
        'MODFLOW 2000 LPF input file ',a)
        go to 9890

9270    write(amessage,9280) trim(rchfile)
9280    format(' Error encountered while reading data from ', &
        'MODFLOW 2000 RCH input file ',a)
        go to 9890

9300    write(amessage,9310) trim(lpffile)
9310    format(' Unexpected end encountered to LPF input file ',a)
        go to 9890

9320    write(amessage,9330) trim(rchfile)
9330    format(' Unexpected end encountered to RCH input file ',a)
        go to 9890

9350    write(amessage,9360) trim(disfile)
9360    format(' Error encountered while reading data from early part of ', &
        'MODFLOW 2000 discretisation input file ',a)
        go to 9890

9400    write(amessage,9410) trim(disfile)
9410    format(' Unexpected end encountered to MODFLOW 2000 discretisation input ', &
        'file ',a)
        go to 9890

9450    write(amessage,9460) trim(facfile)
9460    format(' Error reading factor data from interpolation factor file ',a)
        go to 9890

9500    write(amessage,9510) trim(zonfile)
9510    format(' Error encountered while reading first part of MODFLOW 2000 zone file ',a)
        go to 9890

9550    write(amessage,9560) trim(zonfile)
9560    format(' Unexpected end encountered to MODFLOW 2000 zone file ',a)
        go to 9890

9600    call num2char(iz,anum)
        write(amessage,9610) trim(anum),trim(zonfile)
9610    format(' Error reading zone array number ',a,' from MODFLOW 2000 zone ', &
        'file ',a)
        go to 9890

9650    write(amessage,9660) trim(maxminfile)
9660    format(' Error reading existing ASP input file ',a,'.')
        go to 9890

9670    write(amessage,9680) trim(maxminfile)
9680    format(' Unexpected end encountered to existing ASP input ', &
        'file ',a,',')
        go to 9890

9700    write(amessage,9710) trim(mulfile)
9710    format(' Error encountered while reading first part of MODFLOW 2000 ', &
        'multiplier file ',a)
        go to 9890

9750    write(amessage,9760) trim(mulfile)
9760    format(' Unexpected end encountered to MODFLOW 2000 multiplier file ',a)
        go to 9890

9800    call num2char(im,anum)
        write(amessage,9810) trim(anum),trim(mulfile)
9810    format(' Error reading multiplier array number ',a,' from MODFLOW 2000 ', &
        'multiplier file ',a)
        go to 9890

9850    write(amessage,9860) trim(facfile)
9860    format(' Error reading data from weights file ',a,'.')
        go to 9890



9890    call write_message(leadspace='yes',endspace='yes')
9891    continue
        inquire(unit=tunit,opened=lopened)
        if(lopened) close(unit=tunit)
        if(ncopylpf.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=lpfunit,iostat=ierr)
          call addquote(lpffile,afile)
          call system('copy t###.lpf '//trim(afile)//' > nul')
        end if
        if(ncopyrch.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=rchunit,iostat=ierr)
          call addquote(rchfile,afile)
          call system('copy t###.rch '//trim(afile)//' > nul')
        end if
        if(ncopyzon.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=zonunit,iostat=ierr)
          call addquote(zonfile,afile)
          call system('copy t###.zon '//trim(afile)//' > nul')
        end if
        if(ncopymul.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=mulunit,iostat=ierr)
          call addquote(mulfile,afile)
          call system('copy t###.mul '//trim(afile)//' > nul')
        end if
        if(ncopysen.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=senunit,iostat=ierr)
          call addquote(senfile,afile)
          call system('copy t###.sen '//trim(afile)//' > nul')
        end if
        if(ncopyint.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=maxminunit,iostat=ierr)
          call addquote(maxminfile,afile)
          call system('copy t###.asp '//trim(afile)//' > nul')
        end if
        if(ncopynam.eq.1)then
          close(unit=tunit,iostat=ierr)
          close(unit=nameunit,iostat=ierr)
          call addquote(namefile,afile)
          call system('copy t###.nam '//trim(afile)//' > nul')
        end if
        inquire(file='t###.lpf',exist=lexist)
        if(lexist) call system('del t###.lpf > nul')
        inquire(file='t###.zon',exist=lexist)
        if(lexist) call system('del t###.zon > nul')
        inquire(file='t###.mul',exist=lexist)
        if(lexist) call system('del t###.mul > nul')
        inquire(file='t###.sen',exist=lexist)
        if(lexist) call system('del t###.sen > nul')
        inquire(file='t###.asp',exist=lexist)
        if(lexist) call system('del t###.asp > nul')
        inquire(file='t###.nam',exist=lexist)
        if(lexist) call system('del t###.nam > nul')

9900    call close_files
        call free_point_mem()
        call free_replace_mem()
        deallocate(intarray,mdat,replace,removezone,removemult,nnclu,llayer, &
        ipt,iused,wt,found,removed_par,removed_zon,ia,ra,removed_mul,amul, &
        oldarrayname,oldarrmin,oldarrmax,wpoints,rmdat,stat=ierr)

end program fac2mf2k



      SUBROUTINE U2DINT(ifail,iwrite,ii,jj,IA,IN,inmod,IOUT,iprn)

      implicit integer (i-n), real(a-h,o-z)

      DIMENSION IA(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*200 FNAME
      DATA NUNOPN/99/


      ifail=0
!1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)',err=9000,end=9100) CNTRL

!2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN,ifail1)
      if(ifail1.ne.0) go to 9000
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=INmod
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
!         WRITE(IOUT,15) LOCAT,FNAME
!   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT',I4,':',/1X,A)
!         ICLOSE=1
      ELSE

!2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=9000,end=9100) LOCAT,ICONST,FMTIN,IPRN
    1    FORMAT(I10,I10,A20,I10)
      END IF

!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,ICONST,R,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN,ifail1)
            if(ifail1.ne.0) go to 9000
            FMTIN=CNTRL(ISTART:ISTOP)
!            IF(ICLOSE.NE.0) THEN
!               IF(FMTIN.EQ.'(BINARY)') THEN
!                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM='UNFORMATTED')
!               ELSE
!                  OPEN(UNIT=LOCAT,FILE=FNAME)
!               END IF
!            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'BINARY') LOCAT=-LOCAT
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN,ifail1)
            if(ifail1.ne.0) go to 9000
         END IF
      END IF

!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT) 200,50,90

!4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO ICONST. RETURN.

50    continue
!   50 DO 80 I=1,II
!      DO 80 J=1,JJ
!   80 IA(J,I)=ICONST
!      IF(K.GT.0) WRITE(IOUT,82) ANAME,ICONST,K
!   82 FORMAT(1X,/1X,A,' =',I15,' FOR LAYER',I4)
!      IF(K.LE.0) WRITE(IOUT,83) ANAME,ICONST
!   83 FORMAT(1X,/1X,A,' =',I15)

      if(iwrite.eq.1)then
        write(iout,'(a)') trim(cntrl)
      end if
      RETURN

!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
   90 CONTINUE
!      IF(K.GT.0) THEN
!         WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
!   94    FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      ELSE IF(K.EQ.0) THEN
!         WRITE(IOUT,95) ANAME,LOCAT,FMTIN
!   95    FORMAT(1X,///11X,A,/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      ELSE
!         WRITE(IOUT,96) ANAME,LOCAT,FMTIN
!   96    FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      END IF

      if(locat.eq.inmod)then
      DO 100 I=1,II
      IF(FMTIN.EQ.'(FREE)') THEN
         READ(in,*,err=9000,end=9100) (IA(J,I),J=1,JJ)
      ELSE
         READ(in,FMTIN,err=9000,end=9100) (IA(J,I),J=1,JJ)
      END IF
  100 CONTINUE
      end if

      if(iwrite.eq.1)then
        if(locat.ne.inmod)then
          write(iout,'(a)') trim(cntrl)
        else
          write(iout,8000) iconst,iprn
8000      format('INTERNAL ',i10,' ''(FREE)'' ',i10)
          do i=1,ii
            write(iout,8010) (ia(j,i),j=1,jj)
8010        format(20(i5))
          end do
        end if
      end if
      GO TO 300

!C4C-----LOCAT<0; READ UNFORMATTED RECORD CONTAINING ARRAY VALUES.
  200 LOCAT=-LOCAT
!      IF(K.GT.0) THEN
!         WRITE(IOUT,201) ANAME,K,LOCAT
!  201    FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!     1    1X,'READING BINARY ON UNIT',I4)
!      ELSE IF(K.EQ.0) THEN
!         WRITE(IOUT,202) ANAME,LOCAT
!  202    FORMAT(1X,///11X,A,/
!     1    1X,'READING BINARY ON UNIT',I4)
!      ELSE
!         WRITE(IOUT,203) ANAME,LOCAT
!  203    FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!     1    1X,'READING BINARY ON UNIT',I4)
!      END IF
!      READ(LOCAT)
!      READ(LOCAT) IA
       if(iwrite.eq.1)then
         write(iout,'(a)') trim(cntrl)
       end if

!C5------IF ICONST NOT ZERO THEN MULTIPLY ARRAY VALUES BY ICONST.

300    continue
!  300 IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
!      IF(ICONST.EQ.0) GO TO 320
!      DO 310 I=1,II
!      DO 310 J=1,JJ
!      IA(J,I)=IA(J,I)*ICONST
!  310 CONTINUE

       return

9000   ifail=1
       return
9100   ifail=2
       return

!10-----CONTROL RECORD ERROR.

!  600 IF(K.GT.0) THEN
!         WRITE(IOUT,601) ANAME,K
!  601    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,
!     1     ' FOR LAYER',I4,':')
!      ELSE
!         WRITE(IOUT,602) ANAME
!  602    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
!      END IF
!      WRITE(IOUT,'(1X,A)') CNTRL
!      STOP
      END


      SUBROUTINE URWORD(LINE,ICOL,ISTART,ISTOP,NCODE,N,R,IOUT,IN,ifail)

!-----VERSION 1003 05AUG1992 URWORD
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

      implicit integer (i-n), real(a-h,o-z)

      CHARACTER*(*) LINE
      CHARACTER*20 RW,STRING
!     ------------------------------------------------------------------
!
!1------Set last char in LINE to blank and set ISTART and ISTOP to point
!1------to this blank as a default situation when no word is found.  If
!1------starting location in LINE is out of bounds, do not look for a
!1------word.

      ifail=0

      LINLEN=LEN(LINE)
      LINE(LINLEN:LINLEN)=' '
      ISTART=LINLEN
      ISTOP=LINLEN
      LINLEN=LINLEN-1
      IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
!
!2------Find start of word, which is indicated by first character that
!2------is not a blank and not a comma.
      DO 10 I=ICOL,LINLEN
      IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.',') GO TO 20
10    CONTINUE
      ICOL=LINLEN+1
      GO TO 100

!3------Found start of word.  Look for end.
!3A-----When word is quoted, only a quote can terminate it.
20    IF(LINE(I:I).EQ.'''') THEN
         I=I+1
         IF(I.LE.LINLEN) THEN
            DO 25 J=I,LINLEN
            IF(LINE(J:J).EQ.'''') GO TO 40
25          CONTINUE
         END IF

!3B-----When word is not quoted, space or comma will terminate.
      ELSE
         DO 30 J=I,LINLEN
         IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.',') GO TO 40
30       CONTINUE
      END IF

!3C-----End of line without finding end of word; set end of word to
!3C-----end of line.
      J=LINLEN+1

!4------Found end of word; set J to point to last character in WORD and
!-------set ICOL to point to location for scanning for another word.
40    ICOL=J+1
      J=J-1
      IF(J.LT.I) GO TO 100
      ISTART=I
      ISTOP=J

!5------Convert word to upper case and RETURN if NCODE is 1.
      IF(NCODE.EQ.1) THEN
         IDIFF=ICHAR('a')-ICHAR('A')
         DO 50 K=ISTART,ISTOP
            IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z')  &
                   LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50       CONTINUE
         RETURN
      END IF

!6------Convert word to a number if requested.
100   IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
         RW=' '
         L=20-ISTOP+ISTART
         IF(L.LT.1) GO TO 200
         RW(L:20)=LINE(ISTART:ISTOP)
         IF(NCODE.EQ.2) READ(RW,'(I20)',ERR=200) N
         IF(NCODE.EQ.3) READ(RW,'(F20.0)',ERR=200) R
      END IF
      RETURN

!7------Number conversion error.

200   continue
      ifail=1
      return

!200   IF(NCODE.EQ.3) THEN
!         STRING= 'A REAL NUMBER'
!         L=13
!      ELSE
!         STRING= 'AN INTEGER'
!         L=10
!      END IF

!7A-----If output unit is negative, set last character of string to 'E'.
!      IF(IOUT.LT.0) THEN
!         N=0
!         R=0.
!         LINE(LINLEN+1:LINLEN+1)='E'
!         RETURN

!7B-----If output unit is positive; write a message to output unit.
!      ELSE IF(IOUT.GT.0) THEN
!         IF(IN.GT.0) THEN
!            WRITE(IOUT,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
!         ELSE
!            WRITE(IOUT,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
!         END IF
!201      FORMAT(1X,/1X,'FILE UNIT',I4,' : ERROR CONVERTING "',A,
!     1       '" TO ',A,' IN LINE:',/1X,A)
!202      FORMAT(1X,/1X,'KEYBOARD INPUT : ERROR CONVERTING "',A,
!     1       '" TO ',A,' IN LINE:',/1X,A)

!7C-----If output unit is 0; write a message to default output.
!      ELSE
!         IF(IN.GT.0) THEN
!            WRITE(*,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
!         ELSE
!            WRITE(*,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
!         END IF
!      END IF

!7D-----STOP after writing message.
!      STOP
      END



      SUBROUTINE U2DREL(ifail,iwrite,II,JJ,A,IN,inmod,IOUT,iprn)

      implicit integer (i-n), real(a-h,o-z)

      DIMENSION A(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*200 FNAME
      DATA NUNOPN/99/


      ifail=0
!1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)',err=9000,end=9100) CNTRL

!2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN,ifail1)
      if(ifail1.ne.0) go to 9000
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=INmod
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
!         WRITE(IOUT,15) LOCAT,FNAME
!   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT',I4,':',/1X,A)
!         ICLOSE=1
      ELSE

!2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,err=9000,end=9100) LOCAT,CNSTNT,FMTIN,IPRN
    1    FORMAT(I10,F10.0,A20,I10)
      END IF

!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN,ifail1)
         if(ifail1.ne.0) go to 9000
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN,ifail1)
            if(ifail1.ne.0) go to 9000
            FMTIN=CNTRL(ISTART:ISTOP)
!            IF(ICLOSE.NE.0) THEN
!               IF(FMTIN.EQ.'(BINARY)') THEN
!                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM='UNFORMATTED')
!               ELSE
!                  OPEN(UNIT=LOCAT,FILE=FNAME)
!               END IF
!            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)') LOCAT=-LOCAT
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN,ifail1)
            if(ifail1.ne.0) go to 9000
         END IF
      END IF

!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT) 200,50,90

!4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.

   50 continue
!  50  DO 80 I=1,II
!      DO 80 J=1,JJ
!   80 A(J,I)=CNSTNT
!      IF(K.GT.0) WRITE(IOUT,2) ANAME,CNSTNT,K
!    2 FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR LAYER',I4)
!      IF(K.LE.0) WRITE(IOUT,3) ANAME,CNSTNT
!    3 FORMAT(1X,/1X,A,' =',1P,G14.6)

      if(iwrite.eq.1)then
        write(iout,'(a)') trim(cntrl)
      end if
      RETURN

!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
   90 CONTINUE
!      IF(K.GT.0) THEN
!         WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
!   94    FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      ELSE IF(K.EQ.0) THEN
!         WRITE(IOUT,95) ANAME,LOCAT,FMTIN
!   95    FORMAT(1X,///11X,A,/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      ELSE
!         WRITE(IOUT,96) ANAME,LOCAT,FMTIN
!   96    FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
!     1    1X,'READING ON UNIT',I4,' WITH FORMAT: ',A)
!      END IF

      if(locat.eq.inmod)then
      DO 100 I=1,II
      IF(FMTIN.EQ.'(FREE)') THEN
         READ(in,*,err=9000,end=9100) (A(J,I),J=1,JJ)
      ELSE
         READ(in,FMTIN,err=9000,end=9100) (A(J,I),J=1,JJ)
      END IF
  100 CONTINUE
      end if

      if(iwrite.eq.1)then
        if(locat.ne.inmod)then
          write(iout,'(a)') trim(cntrl)
        else
          write(iout,8000) cnstnt,iprn
8000      format('INTERNAL ',1pg14.7,' ''(FREE)'' ',i10)
          do i=1,ii
            write(iout,8010) (a(j,i),j=1,jj)
8010        format(7(1pg14.7))
          end do
        end if
      end if
      GO TO 300

!4C-----LOCAT<0; READ UNFORMATTED ARRAY VALUES.
  200 LOCAT=-LOCAT
!      IF(K.GT.0) THEN
!         WRITE(IOUT,201) ANAME,K,LOCAT
!  201    FORMAT(1X,///11X,A,' FOR LAYER',I4,/
!     1    1X,'READING BINARY ON UNIT',I4)
!      ELSE IF(K.EQ.0) THEN
!         WRITE(IOUT,202) ANAME,LOCAT
!  202    FORMAT(1X,///1X,A,/
!     1    1X,'READING BINARY ON UNIT',I4)
!      ELSE
!         WRITE(IOUT,203) ANAME,LOCAT
!  203    FORMAT(1X,///1X,A,' FOR CROSS SECTION',/
!     1    1X,'READING BINARY ON UNIT',I4)
!      END IF
!      READ(LOCAT) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,NROW,ILAY
!      READ(LOCAT) A
       if(iwrite.eq.1)then
         write(iout,'(a)') trim(cntrl)
       end if

!5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.

300   continue
!  300 IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
!      ZERO=0.
!      IF(CNSTNT.EQ.ZERO) GO TO 320
!      DO 310 I=1,II
!      DO 310 J=1,JJ
!      A(J,I)=A(J,I)*CNSTNT
!  310 CONTINUE

      return

9000   ifail=1
       return
9100   ifail=2
       return

!C8------CONTROL RECORD ERROR.
!  500 IF(K.GT.0) THEN
!         WRITE(IOUT,501) ANAME,K
!  501    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,
!     1     ' FOR LAYER',I4,':')
!      ELSE
!         WRITE(IOUT,502) ANAME
!  502    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
!      END IF
!      WRITE(IOUT,'(1X,A)') CNTRL
!      STOP
      END





