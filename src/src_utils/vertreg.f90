program vertreg

! -- Program VERTREG adds vertical regularisation prior information to a PEST
!    control file whose parameters are based on pilot points. 

	use defn
	use inter

	implicit none



        integer, parameter     :: MAXPOINTS=2000
        integer, parameter     :: MAXZONES=20
        integer                :: ifail,iline,npar,nobs,npargp,nprior,nobsgp,ierr, &
                                  ipar,i,ifile,punit,ipoint,numzones,np,jj, &
                                  icount,ilog,iz,dunit,ndata,newprior,jlo,iparam,j,kk
        integer                :: pestunit,nfile,iprior,ipp,izz,k,m,l,llo,n,itemp,outunit
        integer                :: npoints(2),ipzone(MAXPOINTS,2),lpoint(MAXPOINTS,2)
        integer                :: izone(MAXZONES)
        integer, allocatable   :: itrans(:),ipartype(:),used(:)
        real                   :: diff,weight,a,b,c,maxweight,minweight,rhs,weighting,sum
        real                   :: zdiff(MAXZONES),ww(2)
        double precision       :: exdist,exdistsq,eee,nnn,mindist,dtemp,rdist
        double precision       :: easting(MAXPOINTS,2),northing(MAXPOINTS,2), &
                                  ee(2),nn(2)
        double precision, allocatable :: deast(:),dnorth(:)
        character (len=1)      :: az,aw
        character (len=15)     :: aline,pestmode,anum,aobsgp,atemp,atemp1, &
                                  aprior,aweight,atrans,atempr,azone,aregulgp, &
                                  priorbase,app,apar1,apar2,priorname
        character (len=200)    :: aprompt
        character (len=200)    :: pestfile,afile,datafile,outfile
        character (len=200)    :: pointsfile(2)
        character (len=10)     :: prefix(2)
        character (len=12)     :: apoint(MAXPOINTS,2)
        character (len=12), allocatable :: apar(:),oldgroup(:)
                                                                                      

	write(amessage,5)
5       format(' Program VERTREG adds regularisation prior information ', &
        'to a PEST control file whose parameters are based on pilot points. ', &
        'Points are linked in the vertical direction according to a user-supplied ', &
        'preferred parameter ratio or difference.')
        call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        nfile=2	
	

! -- The PEST control file is partially read to obtain the names and status
!    of the parameters featured in it.

25	aprompt=' Enter name of existing PEST control file: '
	call open_input_file(ifail,aprompt,pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

        afile=pestfile

        iline=1
        read(pestunit,26,err=9000,end=9050) cline
26      format(a)        
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
        read(pestunit,26,err=9000,end=9050) cline
        iline=iline+1
        read(pestunit,26,err=9000,end=9050) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        pestmode=cline(left_word(2):right_word(2))
        call casetrans(pestmode,'lo')
        if(pestmode.eq.'prediction')then
          write(amessage,50) trim(pestfile)
50        format(' In file ',a,' PEST is asked to run in predictive ',  &
          'analysis mode; VERTREG requires that PEST be run in regularisation ', &
          'mode.')
          go to 9890
        else if(index(pestmode,'regul').eq.0)then
          write(amessage,51) trim(pestfile)
51        format(' PEST must be asked to run in regularisation mode ', &
          'in PEST control file ',a)
          go to 9890
        end if
        iline=iline+1
        read(pestunit,*,err=9000,end=9050) npar,nobs,npargp,nprior,nobsgp
        if(npar.eq.1)then
          write(amessage,60) trim(pestfile)
60        format(' There is only one parameter cited in PEST control file ',a, &
          '. More parameters than this are required if regularisation is to be ', &
          'applied.')
          go to 9890
        end if
        allocate(apar(npar),itrans(npar),ipartype(npar),oldgroup(nobsgp),stat=ierr)
        if(ierr.ne.0) go to 9150
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
80          format(' Incorrect PARTRANS entry at line ',a,' of PEST control file ',a)
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
          read(pestunit,'(a)',err=9000,end=9050) oldgroup(i)
          oldgroup(i)=adjustl(oldgroup(i))
          call casetrans(oldgroup(i),'lo')
        end do
        rewind(unit=pestunit)

! -- The pilot points files are now read.

        ifile=0
90        ifile=ifile+1
          if(ifile.gt.2) go to 210
          write(6,*)
100       continue
          if(ifile.eq.1)then
            write(6,110,advance='no')
110         format(' Enter name of pilot points file for one model layer: ')
          else
            write(6,111,advance='no')
111         format(' Enter name of pilot points file for another model layer: ')
          end if
          read(5,'(a)') pointsfile(ifile)
          if(pointsfile(ifile).eq.' ') go to 100
          if(index(eschar,pointsfile(ifile)(1:2)).ne.0) then
            if(ifile.eq.1)then
              close(unit=pestunit)
              deallocate(apar,itrans,ipartype,oldgroup,stat=ierr)
              if(ierr.ne.0) go to 9200
              write(6,*)
              go to 25
            else
              ifile=0
              go to 90
            end if
          end if
          npoints(ifile)=0
          punit=nextunit()
          open(unit=punit,file=pointsfile(ifile),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,112) trim(pointsfile(ifile))
112         format(' Cannot open file ',a,' - try again.')
            call write_message(leadspace='yes',endspace='yes')
            go to 100
          end if
          iline=0
          ipoint=0
          do
            iline=iline+1
            read(punit,'(a)',end=200) cline
            if(cline.eq.' ') cycle
            call linesplit(ifail,4)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,120) trim(aline),trim(pointsfile(ifile))
120           format(' Insufficient entries on line ',a,' of file ',a)
              go to 9890
            end if
            ipoint=ipoint+1
            if(ipoint.gt.MAXPOINTS)then
              write(amessage,130)
130           format(' Increase MAXPOINTS and re-compile program.')
              go to 9890
            end if
            apoint(ipoint,ifile)=cline(left_word(1):right_word(1))
            call casetrans(apoint(ipoint,ifile),'lo')
            easting(ipoint,ifile)=char2double(ifail,2)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,140) trim(aline),trim(pointsfile(ifile))
140           format(' Cannot read easting from line ',a,' of file ',a)
              go to 9890
            end if
            northing(ipoint,ifile)=char2double(ifail,3)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,150) trim(aline),trim(pointsfile(ifile))
150           format(' Cannot read northing from line ',a,' of file ',a)
              go to 9890
            end if
            ipzone(ipoint,ifile)=char2int(ifail,4)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,155) trim(aline),trim(pointsfile(ifile))
155           format(' Cannot read zone number from line ',a,' of file ',a)
              go to 9890
            end if
          end do
200       npoints(ifile)=ipoint
          if(npoints(ifile).eq.0)then
            write(amessage,205) trim(pointsfile(ifile))
205         format(' No pilot points were found in file ',a)
            go to 9890
          end if            
          call num2char(ipoint,aline)
          write(6,160) trim(aline),trim(pointsfile(ifile))
160       format(' - data for ',a,' points read from file ',a)
          close(unit=punit)
165       write(6,170,advance='no')
170       format(' Enter prefix for corresponding parameter names (<Enter> if none): ')
          read(5,'(a)') prefix(ifile)
          call casetrans(prefix(ifile),'lo')
          prefix(ifile)=adjustl(prefix(ifile))
          if(prefix(ifile)(1:2).eq.'e ')then
            write(6,*)
            go to 100
          end if
          go to 90
210     continue

! -- It is now verified that at least some parameters for each set of pilot points
!    is present in the PEST control file. In doing this, parameters are placed into
!    the two pilot point groups.

        ipartype=0      ! ipartype is an array
        do ifile=1,nfile
          icount=0
          do ipoint=1,npoints(ifile)
            atemp=trim(prefix(ifile))//trim(apoint(ipoint,ifile))
            do ipar=1,npar
              if(apar(ipar).eq.atemp)then
                ipartype(ipar)=ifile
                icount=icount+1
                go to 220
              end if
            end do
220         continue
          end do
          if(icount.eq.0)then
            write(amessage,230) trim(pestfile),trim(pointsfile(ifile))
230         format(' There appear to be no parameters in the PEST control file ', &
            a,' which are derived from pilot points cited in the pilot points ', &
            'file ',a)
            go to 9890
          end if
        end do
        
! -- We now verify that all pertinent adjustable parameters are log transformed
!    or that they are all not log-transformed.

        do ifile=1,2
          icount=0
          do ipar=1,npar
            if(ipartype(ipar).eq.ifile)then
              if(itrans(ipar).ge.0)then
                icount=icount+1
              end if
            end if
          end do
          if(icount.eq.0)then
            write(amessage,240) trim(pestfile),trim(pointsfile(ifile))
240         format(' All pilot point parameters in PEST control file file ',a,  &
            ' derived from pilot points in file ',a,' are either tied or fixed.')
            go to 9890
          end if
        end do
        
        ilog=-999
        do ipar=1,npar
          if(ipartype(ipar).ne.0)then
            if(itrans(ipar).ge.0)then
              if(ilog.eq.-999)then
                if(itrans(ipar).eq.1)then
                  ilog=1
                else
                  ilog=0
                end if
              else
                if(itrans(ipar).ne.ilog)then
                  write(amessage,250) trim(pestfile),trim(pointsfile(1)),trim(pointsfile(2))
250               format(' Pilot point parameters in file ',a,' which are derived ', &
                  'from pilot points in files ',a,' and ',a,' and which are neither ',  &
                  'tied nor fixed must be either ALL log transformed or ALL untransformed.')
                  go to 9890
                end if
              end if
            end if
          end if
        end do   
             
! -- Zonal data are now obtained.

        write(6,*)
260     continue
        if(ilog.eq.1)then
          write(6,270,advance='no')
270       format(' Will regularised vertical property ratio depend on zones? [y/n]: ')
        else
          write(6,280,advance='no')
280       format(' Will regularised vertical property difference depend on zones? [y/n]: ')
        end if
        read(5,'(a)') az
        call casetrans(az,'lo')
        if(az.eq.' ') go to 260
        if(az.eq.'e') then
          write(6,*)
          ifile=2
          go to 100
        end if
        if((az.ne.'y').and.(az.ne.'n')) go to 260
290     continue
        iz=0
        if(az.eq.'y')then
299       write(6,300)
300       format(' From which pilot points file is zonation taken?')
          write(6,310) trim(pointsfile(1))
310       format('     enter 1 if from file ',a,':')
          write(6,320) trim(pointsfile(2))
320       format('     enter 2 if from file ',a,':')
329       write(6,330,advance='no')
330       format(' Enter your choice: ')
          read(5,'(a)') atemp
          if(atemp.eq.' ') go to 329
          if(index(eschar,atemp(1:2)).ne.0) then
            write(6,*)
            go to 260
          end if
          call char2num(ifail,atemp,iz)
          if(ifail.ne.0) go to 329
          if((iz.ne.1).and.(iz.ne.2)) go to 329
        end if
        
! -- If the user has indicated that differences or ratios are zone-specific, then
!    the zones in the pertinent pilot points file are sorted.

        if(iz.gt.0)then
          numzones=0
          np=npoints(iz)
          itemp=999999999
          do i=1,np
            if(ipzone(i,iz).lt.itemp)then
              itemp=ipzone(i,iz)
            end if
          end do
          numzones=1
          izone(numzones)=itemp
          do
            icount=0
            itemp=999999999
            do i=1,np
              if(ipzone(i,iz).lt.itemp)then
                if(ipzone(i,iz).gt.izone(numzones))then
                  itemp=ipzone(i,iz)
                  icount=icount+1
                end if
              end if
            end do
            if(icount.gt.0)then
              numzones=numzones+1
              if(numzones.gt.MAXZONES)then
                write(amessage,340)
340             format(' Increase MAXZONES and re-compile program.')
                go to 9890
              end if
              izone(numzones)=itemp
            else
              go to 347
            end if
          end do
        end if

! -- We now ascertain the parameter difference or ratio to use for regularisation.

347     continue
        if(iz.eq.0)then
348       if(ilog.eq.1)then
349         write(6,350,advance='no')
350         format(' Enter (first/second) parameter ratio for vertical regularisation: ')
          else
359         write(6,360,advance='no')
360         format(' Enter (first-second) parameter difference for vertical regularisation: ')
          end if
          if(key_read(diff).ne.0) go to 347
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 260
          end if
          if(ilog.eq.1)then
            if(diff.le.0.0) then
              write(6,361)
361           format(/,' Ratio must be positive - try again.',/)
              go to 348
            end if
            diff=log10(diff)
          end if
	else
          i=0
370       i=i+1
          if(i.gt.numzones)go to 400
          call num2char(izone(i),azone)
          azone=adjustl(azone)
388       if(ilog.eq.1)then
            write(6,390,advance='no') trim(azone)
390         format(' Enter (first/second) regularisation parameter ratio for zone ',a,': ')
          else
            write(6,395,advance='no') trim(azone)
395         format(' Enter (first-second) regularisation parameter difference for zone ',a,': ')
          end if
          if(key_read(zdiff(i)).ne.0) go to 388
          if(escset.ne.0)then
            escset=0
            write(6,*)
            if(i.ge.2)then
              i=i-2
              go to 370
            else
              go to 290
            end if
          end if
          if(ilog.eq.1)then
            if(zdiff(i).le.0.0) then
              write(6,361)
              go to 388
            end if
            zdiff(i)=log10(zdiff(i))
          end if
          go to 370
400       continue          
        end if  
            
! The exclusion distance is read.

        write(6,*)
420     write(6,430,advance='no')
430     format(' Enter lateral exclusion distance: ')
        if(key_read(exdist).ne.0) go to 420
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 347
        end if
        if(exdist.le.0.0) then
          write(6,440)
440       format(/,' Must be a postive number - try again.',/)
          go to 420
        end if
        exdistsq=exdist*exdist

! -- And now the regularisation subgroup name.

        write(6,*)
450     write(6,460,advance='no')
460     format(' Enter name for new regularisation subgroup: ')
        read(5,'(a)') aregulgp
        if(aregulgp.eq.' ') go to 450
        call casetrans(aregulgp,'lo')
        if(aregulgp(1:2).eq.'e ')then
          write(6,*)
          go to 420
        end if
        if(aregulgp(1:5).ne.'regul')then
          write(6,470)
470       format(/,' This name must begin with "regul" - try again.',/)
          go to 450
        end if
        if(len_trim(aregulgp).gt.12)then
          write(6,475)
475       format(/,' This must be less than 12 characters - try again.',/)
          go to 450
        end if
        do i=1,nobsgp
          if(aregulgp.eq.oldgroup(i))then
            write(amessage,490)
490         format(' Observation group name is already used in PEST control ', &
            'file - try again.')
            call write_message(leadspace='yes',endspace='yes')
            go to 450
          end if
        end do
        
! -- The basename for the new prior information is now supplied.

491     write(6,492,advance='no')
492     format(' Enter basename for new prior information: ')
        read(5,'(a)') priorbase
        if(priorbase.eq.' ') go to 491
        call casetrans(priorbase,'lo')
        if(priorbase(1:2).eq.'e ')then
          write(6,*)
          go to 450
        end if
        if(len_trim(priorbase).gt.6)then
          write(6,495)
495       format(/,' This name must be less than 6 characters in length - try again',/)
          go to 491
        end if

! -- Do we need data-density-dependent weighting?

        write(6,*)
500     write(6,510,advance='no')
510     format(' Do you require data-density-dependent weighting?  [y/n]: ')
        read(5,'(a)') aw
        if(aw.eq.' ') go to 500
        call casetrans(aw,'lo')
        if(aw.eq.'e')then
          write(6,*)
          go to 491
        end if
        if((aw.ne.'y').and.(aw.ne.'n')) go to 500
        if(aw.eq.'n')then
493       write(6,494,advance='no')
494       format(' Enter weight for all new prior information: ')
          if(key_read(weight).ne.0) go to 493
	  if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 500
          end if
          if(weight.lt.0.0)then
            write(6,496)
496         format(/,' Must be positive - try again.',/)
            go to 493
          end if
          go to 700
        end if
511     write(6,512,advance='no')
512     format(' Enter name of data coordinates file: ')
        read(5,'(a)') datafile
        if(datafile.eq.' ') go to 511
        if(index(eschar,datafile(1:2)).ne.0) then
          write(6,*)
          go to 500
        end if
        dunit=nextunit()
        open(unit=dunit,file=datafile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(6,513) trim(datafile)
513       format(/,' Cannot open file ',a,' - try again.',/)
          go to 511
        end if
        ndata=0
        do
          read(dunit,'(a)',end=520) cline
          if(cline.eq.' ') cycle
          ndata=ndata+1
        end do
520     continue
        if(ndata.eq.0)then
          write(amessage,521) trim(datafile)
521       format(' No data was found in file ',a,' - try again.')
          call write_message(leadspace='yes',endspace='yes')
          close(unit=dunit)
          go to 511
        end if
        allocate(deast(ndata),dnorth(ndata),used(ndata),stat=ierr)
        if(ierr.ne.0) go to 9150
        rewind(unit=dunit)
        iline=0
        do i=1,ndata
515       iline=iline+1
          read(dunit,'(a)') cline
          if(cline.eq.' ') go to 515
          call linesplit(ifail,3)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,525) trim(aline),trim(datafile)
525         format(' Insufficient entries on line ',a,' of file ',a)
            go to 9890
          end if
          deast(i)=char2double(ifail,2)
          if(ifail.ne.0) go to 9600
          dnorth(i)=char2double(ifail,3)
          if(ifail.ne.0) go to 9600
        end do
        close(unit=dunit)
        call num2char(ndata,atemp)
        write(6,530) trim(atemp),trim(datafile)
530     format(' - ',a,' data points read from file ',a)
        write(6,560)
560     format(' Weights function is "a+b*sum_over_closest_n(r**c)" :-')
570     write(6,580,advance='no')
580     format('     Enter a: ')
        if(key_read(a).ne.0) go to 570
        if(escset.ne.0)then
          escset=0
          write(6,*)
          deallocate(deast,dnorth,used,stat=ierr)
          if(ierr.ne.0) go to 9200
          go to 511
        end if
590     write(6,600,advance='no')
600     format('     Enter b: ')
        if(key_read(b).ne.0) go to 590
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 570
        end if
610     write(6,615,advance='no')
615     format('     Enter n: ')
        if(key_read(n).ne.0) go to 610
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 590
        end if
        if(n.le.0)then
          write(6,616)
616       format(/,' Must be greater than zero - try again.',/)
          go to 610
        end if 
620     write(6,630,advance='no')
630     format('     Enter c: ')
        if(key_read(c).ne.0) go to 620
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 610
        end if
        if(c.lt.0.0)then
          write(6,676)
          go to 620
        end if
640     write(6,650,advance='no')
650     format('     Enter maximum allowable weight: ')
        if(key_read(maxweight).ne.0) go to 640
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 620
        end if
        if(maxweight.le.0.0)then
          write(6,616)
          go to 640
        end if 
660     write(6,670,advance='no')
670     format('     Enter minimum allowable weight: ')
        if(key_read(minweight).ne.0) go to 660
	if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 640
        end if
        if(minweight.ge.maxweight)then
          write(6,675)
675       format(/,' Minimum weight must be less than maximum weight - try again.',/)
          go to 660
        end if
        if(minweight.lt.0.0)then
          write(6,676)
676       format(/,' Must not be negative - try again.',/)
          go to 660
        end if
       
700     continue
        write(6,*)
710     write(6,720,advance='no')
720     format(' Enter name for new PEST control file: ')
        read(5,'(a)') outfile
        if(outfile.eq.' ') go to 710
        if((outfile(1:2).eq.'e ').or.(outfile(1:2).eq.'E '))then
          write(6,*)
          if(aw.eq.'y')then
            go to 660
          else
            go to 493
          end if
        end if
        outunit=nextunit()
        open(unit=outunit,file=outfile,action='write',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,730) trim(outfile)
730       format(' Cannot open file ',a,' for output - try again.')
          call write_message(leadspace='yes',endspace='yes')
          go to 710
        end if
        
        write(6,740)
740     format(/,' Working....')

! -- First we evaluate how many new items of prior information will be generated.

        lpoint=0                              ! lpoint is an array.
        
        newprior=0
        do i=1,npoints(1)
          call is_adj_param(npar,apar,itrans,apoint(i,1),prefix(1),iparam)
          if(iparam.eq.0) cycle
          eee=easting(i,1)
          nnn=northing(i,1)
          mindist=1.0d300
          jlo=0
          do j=1,npoints(2)
            call is_adj_param(npar,apar,itrans,apoint(j,2),prefix(2),iparam)
            if(iparam.eq.0) go to 745
            rdist=(eee-easting(j,2))*(eee-easting(j,2))+     &
                  (nnn-northing(j,2))*(nnn-northing(j,2))
                  if(rdist.lt.mindist)then
                    jlo=j
                    mindist=rdist
                  end if
745               continue
          end do
          if(mindist.le.exdistsq)then
            newprior=newprior+1
            lpoint(i,1)=jlo
            lpoint(jlo,2)=-i
          else
            lpoint(i,1)=-999999999
          end if
        end do
        do i=1,npoints(2)
          if(lpoint(i,2).ne.0) cycle
          call is_adj_param(npar,apar,itrans,apoint(i,2),prefix(2),iparam)
          if(iparam.eq.0) cycle
          eee=easting(i,2)
          nnn=northing(i,2)
          mindist=1.0d300
          jlo=0
          do j=1,npoints(1)
            call is_adj_param(npar,apar,itrans,apoint(j,1),prefix(1),iparam)
            if(iparam.eq.0) cycle
            rdist=(eee-easting(j,1))*(eee-easting(j,1))+     &
                  (nnn-northing(j,1))*(nnn-northing(j,1))
                  if(rdist.lt.mindist)then
                    jlo=j
                    mindist=rdist
                  end if
          end do
          if(mindist.le.exdistsq)then
            newprior=newprior+1
            lpoint(i,2)=jlo
            if(lpoint(jlo,1).eq.0) lpoint(jlo,1)=-i
          else
            lpoint(i,2)=-999999999
          end if
        end do
        if(newprior.eq.0)then
          write(amessage,746)
746       format(' No new prior information can be added based on vertical ',  &
          'regularisation. Is the exclusion distance set too low? ')
          go to 9890
        end if

! -- The existing PEST control file is copied up until the point where
!    some control variables need to be incremented.

        afile=pestfile
        iline=0
        do i=1,3
          iline=iline+1
          read(pestunit,750,err=9000,end=9050) cline
          write(outunit,750,err=9550) trim(cline)
750       format(a)
        end do
        iline=iline+1
        read(pestunit,750,err=9000,end=9050) cline
        write(outunit,751,err=9550) npar,nobs,npargp,nprior+newprior,nobsgp+1
751     format(1x,5(i6,1x))
        
! -- The existing PEST control file is now copied to the new one up until the stage 
!    where the new observation group must be added.

        do
          iline=iline+1
          read(pestunit,750,err=9000,end=9090) cline
          write(outunit,750,err=9550) trim(cline)
          if(index(cline,'* observation gr').ne.0) go to 760
        end do
760     continue
        do i=1,nobsgp
          iline=iline+1
          read(pestunit,750,err=9000,end=9050) cline
          write(outunit,750,err=9550) trim(cline)
        end do
        write(outunit,750,err=9550) ' '//trim(aregulgp)
        
! -- We now look for a place to write the new prior information. This will be just above
!    the regularisation section, as we assume that the existing file for a PEST run is in 
!    regularisation mode.

        iprior=0
        do
          iline=iline+1
          read(pestunit,750,err=9000,end=9400) cline
          if(index(cline,'* regul').ne.0) go to 770
          if(index(cline,'* prior informati').ne.0) iprior=1
          write(outunit,750,err=9550) trim(cline)
        end do
770     continue

! -- The extra prior information is now added.

        if(iprior.eq.0)then
          write(outunit,780,err=9550)
780       format('* prior information')
        end if
        ipp=0
        do ifile=1,2
          do i=1,npoints(ifile)
            if(lpoint(i,ifile).le.0) cycle
            if(iz.eq.0)then
              rhs=diff
            else
              if(iz.eq.ifile)then
                izz=ipzone(i,ifile)
              else
                jj=lpoint(i,ifile)
                izz=ipzone(jj,3-ifile)
              end if
              do kk=1,numzones
                if(izz.eq.izone(kk))then
                  rhs=zdiff(kk)
                  go to 781
                end if
              end do
781           continue
            end if
            if(aw.eq.'n')then
              weighting=weight
            else
              ee(1)=easting(i,ifile)
              nn(1)=northing(i,ifile)
              ee(2)=easting(lpoint(i,ifile),3-ifile)
              nn(2)=northing(lpoint(i,ifile),3-ifile)
              do k=1,2
                if(k.eq.2)then
                  if((ee(2).eq.ee(1)).and.(nn(2).eq.nn(1)))then
                    ww(2)=ww(1)
                    go to 800
                  end if
                end if
                used=0            ! used is an array
                sum=0.0d0
                do m=1,min(n,ndata)
                  mindist=1.0d300
                  llo=0
                  do l=1,ndata
                    if(used(l).ne.0) cycle
                    dtemp=(ee(k)-deast(l))*(ee(k)-deast(l)) +    &
                         (nn(k)-dnorth(l))*(nn(k)-dnorth(l))
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
                      if(c.gt.0.0)then
                        sum=sum+sqrt(mindist)**c
                      else
                        sum=sum+1.0
                      end if
                    end if
                  else
                    go to 799
                  end if
                end do
799             continue
                ww(k)=a+b*sum
              end do
800           weighting=max(ww(1),ww(2))
              if(weighting.gt.maxweight) weighting=maxweight
              if(weighting.lt.minweight) weighting=minweight
            end if
            ipp=ipp+1
            call num2char(ipp,app)
            priorname=trim(priorbase)//trim(app)
            if(ifile.eq.1)then
              apar1=trim(prefix(ifile))//trim(apoint(i,ifile))
              apar2=trim(prefix(3-ifile))//trim(apoint(lpoint(i,ifile),3-ifile))
            else
              apar2=trim(prefix(ifile))//trim(apoint(i,ifile))
              apar1=trim(prefix(3-ifile))//trim(apoint(lpoint(i,ifile),3-ifile))
            end if
            if(ilog.eq.0)then
              write(outunit,810) trim(priorname),trim(apar1),trim(apar2),rhs,  &
              weighting,trim(aregulgp)
810           format(a,t14,'1.0 * ',a,' - 1.0 * ',a,' = ',1pg13.5,2x,1pg13.5,2x,a)
            else
              write(outunit,820) trim(priorname),trim(apar1),trim(apar2),rhs,  &
              weighting,trim(aregulgp)
820           format(a,t14,'1.0 * log(',a,') - 1.0 * log(',a,') = ',1pg13.5,2x,1pg13.5,2x,a)
            end if
          end do
        end do
        write(outunit,830,err=9550)
830     format('* regularisation')
        do
          read(pestunit,'(a)',end=900) cline
          write(outunit,'(a)',err=9550) trim(cline)
        end do
900     continue
        write(6,910) trim(outfile)
910     format(' - file ',a,' written ok.')

! -- We now list which parameters have not been linked.

        icount=0
        do ifile=1,2
          do i=1,npoints(ifile)
            if(lpoint(i,1).eq.-999999999) icount=icount+1
          end do
        end do
        if(icount.ne.0)then
          write(amessage,920) trim(pointsfile(1)),trim(pointsfile(2))
920       format(' The following non-fixed and non-tied parameters, ',  &
          'built from pilot points in files ',  &
          a,' and ',a,' are not featured in new prior information. Increase ',  &
          'the exclusion distance to include these points if you wish:-')
          call write_message(leadspace='yes',endspace='yes')
          cline=' '
          icount=0
          do ifile=1,2
            do i=1,npoints(ifile)
              if(lpoint(i,ifile).eq.-999999999) then
                cline=trim(cline)//'  '//trim(prefix(ifile))//trim(apoint(i,ifile))
                icount=icount+1
                if(len_trim(cline).gt.68)then
                  write(6,'(a)') trim(cline)
                  cline=' '
                  icount=0
                end if
              end if
            end do
          end do
          if(icount.ne.0) write(6,'(a)') trim(cline)
        end if
        go to 9900
       
9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a)
        go to 9890
9050    write(amessage,9060) trim(afile)
9060    format(' Unexpected end encountered to file ',a)
        go to 9890
9070    write(amessage,9080) trim(pestfile)
9080    format(' Cannot find "parameter data" section of PEST control file ',a)
        go to 9890
9090    write(amessage,9095) trim(pestfile)
9095    format(' Cannot find "observation groups" section of PEST control file ',a)
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Insufficient entries on line ',a,' of file ',a)
        go to 9890
9150    write(amessage,9160)
9160    format(' Cannot allocate sufficient memory to continue ',&
        'execution.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Memory management error.')
        go to 9890
9400    write(6,9410) trim(pestfile)
9410    format(' Cannot find "regularisation" section of PEST control file ',a)
        go to 9890
9550    write(amessage,9560) trim(outfile)
9560    format(' Cannot write to file ',a)
        go to 9890
9600    call num2char(iline,aline)
        write(amessage,9610) trim(aline),trim(datafile)
9610    format(' Cannot read coordinates from line ',a,' of file ',a)
        go to 9890


9890    call write_message(leadspace='yes')
9900    call close_files
        deallocate(apar,itrans,ipartype,oldgroup,deast,dnorth,used,stat=ierr)
        write(6,*)

end program vertreg
 

subroutine is_adj_param(npar,apar,itrans,apoint,prefix,iparam)

        implicit none
        
        integer            :: npar,iparam,i
        integer            :: itrans(npar)
        character (len=*)  :: apar(npar)
        character (len=*)  :: prefix,apoint
        
        iparam=0
        do i=1,npar
          if(itrans(i).lt.0) cycle
          if(trim(prefix)//apoint.eq.apar(i))then
            iparam=i
            return
          end if
        end do
        
        return

end subroutine is_adj_param 
