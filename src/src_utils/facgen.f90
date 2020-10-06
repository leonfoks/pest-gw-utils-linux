!     Last change:  J     3 Mar 2005    1:10 am

program facgen

! -- Program FACGEN is part of the PEST-TETRAD interface. It calculates kriging
!    factors on the basis of information contained in a pilot points parameterization
!    file and keyboard entries provided by the user.

	use defn
	use inter

	implicit none

! -- Variables are declared.

        integer, parameter :: MAXPP=3000

        logical   :: lexist
        integer   :: ifail,pppunit,rdfunit,rpfunit,isunit,iline,nx,ny,nz,  &
                     itemp,nregion,ierr,iregion,nline,il,i1,i2,j1,j2,k1,k2, &
                     i,j,k,nproperties,ix,iy,iz,iblock,ipfile,iproperty,iffile, &
                     ibeg,iend,ilog,ivar,maxpoints,minpoints,punit, &
                     jline,ipp,iloop,ncol,nrow,outunit,ngrid,minpt,maxpt, &
                     numhorcell,npp,noest,ilb,iub,itr,iiv,ipptrans

        integer   :: k_ktype,n_nst,itrans,ndat
        integer   :: i_it(1)
        real      :: s_skmean,c_c0,pmx,radmax
        real      :: c_cc(1),a_aa(1),a_ang(1),a_anis(1)

        double precision :: e_corner,n_corner,angle,a,anis,anis_angle,radius

        double precision :: lppbound,uppbound
        character*1   :: aa,aoutfile
        character*15  :: aregion,aproperty
        character*15  :: aline,atrans
        character*40  :: akey,atemp,anum
        character*120 :: aprompt
        character*200 :: pppfile,rdffile,rpffile,isfile,afile,bfile,pfile, &
                         ffile,cfile,apfile,affile

        type(modelgrid)  :: gridspec

        integer          :: inumdat(MAXPP)
        real             :: valdat(MAXPP),eastdat(MAXPP),northdat(MAXPP)
        double precision :: ep(MAXPP),np(MAXPP)

        integer,allocatable           :: icellno(:),igrid(:,:,:)
        real,allocatable              :: east(:,:),north(:,:),eastgrid(:),  &
                                         northgrid(:),variance(:,:),dx(:),dy(:)
        character*12, allocatable     :: aregname(:),apropname(:)

! -- An initial screen message is written.


	write(amessage,5)
5	format(' Program FACGEN generates interpolation factor files ',    &
        'corresponding to pilot points files named in a pilot points ',  &
        'parameterization file.')
	call write_message(leadspace='yes',endspace='yes')

! -- Initialization

        iloop=0

! -- Information is acquired from the user.

20      continue
        aprompt=' Enter name of pilot points parameterization file : '
        call open_input_file(ifail,aprompt,pppfile,pppunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

40	continue
	aprompt=' Enter name of regions definition file            : '
	call open_input_file(ifail,aprompt,rdffile,rdfunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
          close(unit=pppunit)
	  go to 20
	end if

50	continue
        aprompt=' Enter name of region properties file             : '
        call open_input_file(ifail,aprompt,rpffile,rpfunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          write(6,*)
          escset=0
          close(unit=rdfunit)
          go to 40
        end if

60      continue
        aprompt=' Enter name of IS file for current model          : '
        call open_input_file(ifail,aprompt,isfile,isunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
          close(unit=rpfunit)
	  go to 50
	end if

        write(6,*)
80      continue
        write(6,90,advance='no')
90      format(' Enter easting of top left grid corner : ')
	itemp=key_read(e_corner)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
          close(unit=isunit)
	  go to 60
	else if(itemp.eq.-1) then
	  go to 80
	else if(itemp.ne.0) then
	  write(6,100)
100	  format(' Illegal input  - try again.')
	  go to 80
	end if

110     continue
        write(6,120,advance='no')
120     format(' Enter northing top left grid corner   : ')
	itemp=key_read(n_corner)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 80
	else if(itemp.eq.-1) then
	  go to 110
	else if(itemp.ne.0) then
	  write(6,100)
	  go to 110
	end if

130     continue
        write(6,140,advance='no')
140     format(' Enter bearing of grid row direction   : ')
	itemp=key_read(angle)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 110
	else if(itemp.eq.-1) then
	  go to 130
	else if(itemp.ne.0) then
	  write(6,100)
	  go to 130
	end if
        if((angle.lt.0.0d0).or.(angle.gt.180.0d0))then
          write(6,150)
150       format(' Bearing must be between 0.0 and 180 degrees - try again.')
          go to 130
        end if
        angle=90.0d0-angle

! -- The region definitions file is first read.

        write(6,*)
        call addquote(rdffile,afile)
        write(6,155) trim(afile)
155     format(' - reading regions definition file ',a,'...')
        iline=1
        read(rdfunit,'(a)',err=9000,end=9050) cline
        call linesplit(ifail,3)
        if(ifail.ne.0) go to 9150
        nx=char2int(ifail,1)
        if(ifail.ne.0) go to 9000
        ny=char2int(ifail,2)
        if(ifail.ne.0) go to 9000
        nz=char2int(ifail,3)
        if(ifail.ne.0) go to 9000
        if((nx.le.0).or.(ny.le.0).or.(nz.le.0)) then
          call num2char(iline,aline)
          write(amessage,160) trim(aline),trim(afile)
160       format(' Grid dimension cannot be zero or negative at line ',a,  &
          ' of file ',a,'.')
          go to 9890
        end if

        iline=2
        read(rdfunit,'(a)',err=9000,end=9050) cline
        call linesplit(ifail,1)
        if(ifail.ne.0) go to 9150
        nregion=char2int(ifail,1)
        if(ifail.ne.0) go to 9000

        allocate(aregname(nregion),igrid(nx,ny,nz),variance(nx,ny),   &
        dx(nx),dy(ny),stat=ierr)
        if(ierr.ne.0) go to 9100
        igrid=0    ! an array
        numhorcell=nx*ny
        allocate(eastgrid(numhorcell),northgrid(numhorcell),icellno(numhorcell), &
        stat=ierr)
        if(ierr.ne.0) go to 9100
        allocate(east(nx,ny),north(nx,ny),stat=ierr)
        if(ierr.ne.0) go to 9100

        do iregion=1,nregion
          iline=iline+1
          read(rdfunit,'(a)',err=9000,end=9050) cline
          call linesplit(ifail,2)
          if(ifail.ne.0) go to 9150
          aregname(iregion)=cline(left_word(1):right_word(1))
          call casetrans(aregname(iregion),'lo')
          nline=char2int(ifail,2)
          if(ifail.ne.0) go to 9000

          do il=1,nline
            iline=iline+1
            read(rdfunit,'(a)',err=9000,end=9050) cline
            call linesplit(ifail,6)
            if(ifail.ne.0) go to 9150
            i1=char2int(ifail,1)
            if(ifail.ne.0) go to 9000
            i2=char2int(ifail,2)
            if(ifail.ne.0) go to 9000
            j1=char2int(ifail,3)
            if(ifail.ne.0) go to 9000
            j2=char2int(ifail,4)
            if(ifail.ne.0) go to 9000
            k1=char2int(ifail,5)
            if(ifail.ne.0) go to 9000
            k2=char2int(ifail,6)
            if(ifail.ne.0) go to 9000
            if((i1.le.0).or.(j1.le.0).or.(k1.le.0).or.       &
               (i2.gt.nx).or.(j2.gt.ny).or.(k2.gt.nz))then
               call num2char(iline,aline)
               write(amessage,170) trim(aline),trim(afile)
170            format(' Index out of range at line ',a,' of file ',a,'.')
               go to 9890
            end if
            if((i2.lt.i1).or.(j2.lt.j1).or.(k2.lt.k1))then
              call num2char(iline,aline)
              write(amessage,171) trim(aline),trim(afile)
171           format(' Second index less than first index for either ', &
              'IX, IY, IZ at line ',a,' of file ',a,'.')
              go to 9890
            end if
            do k=k1,k2
              do j=j1,j2
                do i=i1,i2
                  igrid(i,j,k)=iregion
                end do
              end do
            end do
          end do
        end do

! -- A check is made that all grid blocks are assigned to a region.

       if(any(igrid.eq.0))then
         write(amessage,200) trim(afile)
200      format(' Some grid blocks are not assigned to a region in ',   &
         'region definitions file ',a,'.')
         go to 9890
       end if

! -- The RDF file can now be closed.

       close(unit=rdfunit)
       write(6,220) trim(afile)
220    format(' - file ',a,' read ok.')


! -- The names of estimable properties are now read from the region
!    properties file.

       write(6,*)
       call addquote(rpffile,afile)
       write(6,225) trim(afile)
225    format(' - reading region properties file ',a,'...')
       iline=1
       call num2char(iline,aline)
       read(rpfunit,'(a)',err=9500,end=9550) cline
       if(cline.eq.' ') go to 9500
       call linesplit(ifail,2)
       if(ifail.ne.0) go to 9150
       itemp=char2int(ifail,1)
       if(ifail.ne.0) go to 9500
       if(itemp.le.0)then
         write(amessage,224) trim(aline),trim(afile)
224      format(' Number of regions must be positive at line ',a,   &
         ' of file ',a,'.')
         go to 9890
       end if
       if(itemp.ne.nregion)then
         call addquote(rdffile,bfile)
         write(amessage,222) trim(bfile),trim(afile)
222      format(' Number of regions cited in regions definition file ',a,   &
         ' differs from number of regions cited in region properties file ',a,'.')
         go to 9890
       end if
       nproperties=char2int(ifail,2)
       if(ifail.ne.0) go to 9500
       if(nproperties.le.0)then
         write(amessage,226) trim(aline),trim(afile)
226      format(' Number of properties must be positive at line ',a,   &
         ' of file ',a,'.')
         go to 9890
       end if
       allocate(apropname(nproperties),stat=ierr)
       if(ierr.ne.0) go to 9100


       iline=2
       call num2char(iline,aline)
       read(rpfunit,'(a)',err=9500,end=9550) cline
       if(cline.eq.' ') go to 9500
       if(nproperties.gt.NUM_WORD_DIM-1)then
         write(amessage, 223)
223      format(' Too many properties - increase NUM_WORD_DIM and re-compile ', &
         'program.')
         go to 9890
       end if
       call linesplit(ifail,nproperties+1)
       if(ifail.ne.0) go to 9150
       atemp=cline(left_word(1):right_word(1))
       call casetrans(atemp,'lo')
       if(index(atemp,'region').eq.0)then
         write(amessage,227) trim(aline),trim(afile)
227      format(' First entry on line ',a,' of file ',a,' must be "REGION".')
         go to 9890
       end if
       do i=1,nproperties
         apropname(i)=cline(left_word(i+1):right_word(i+1))
         call casetrans(apropname(i),'lo')
       end do

! -- Enough of the region properties file has been read - so close it now.

       close(unit=rpfunit)
       write(6,220) trim(afile)


! -- Grid DXs and DYs are now read from the IS file.

       write(6,*)
       call addquote(isfile,afile)
       write(6,229) trim(afile)
229    format(' - reading IS file ',a,'...')
       do i=1,2
         read(isunit,'(a)',err=9200,end=9250) cline
       end do
       read(isunit,*,err=9200,end=9250) ix
       read(isunit,*,err=9200,end=9250) iy
       read(isunit,*,err=9200,end=9250) iz
       if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz))then
         call addquote(rdffile,bfile)
         write(amessage,230) trim(afile),trim(bfile)
230      format(' Grid dimensions in IS file ',a,' do not match those supplied ', &
         'in regions definitions file ',a,'.')
         go to 9890
       end if
       do
         read(isunit,'(a)',err=9200,end=250) cline
         cline=adjustl(cline)
         call casetrans(cline,'lo')
         if(cline(1:3).eq.'dx ') go to 270
         cycle
250      write(amessage,260) trim(afile)
260      format(' Cannot find DX array in IS file ',a,'.')
         go to 9890
       end do
270    continue
       read(isunit,*,iostat=ierr) (dx(i),i=1,nx)
       if(ierr.ne.0)then
         write(amessage,280) trim(isfile)
280      format(' Error encountered while reading DX array from IS file ',a,'.')
         go to 9890
       end if
       do
         read(isunit,'(a)',err=9200,end=290) cline
         cline=adjustl(cline)
         call casetrans(cline,'lo')
         if(cline(1:3).eq.'dy ') go to 310
         cycle
290      write(amessage,300) trim(afile)
300      format(' Cannot find DY array in IS file ',a,'.')
         go to 9890
       end do
310    continue
       do iy=1,ny
         read(isunit,*,iostat=ierr) (variance(ix,iy),ix=1,nx)
       end do
       do iy=1,ny
         dy(iy)=variance(1,iy)
       end do
       if(ierr.ne.0)then
         write(amessage,320) trim(isfile)
320      format(' Error encountered when reading DY array from IS file ',a,'.')
         go to 9890
       end if

! -- The IS file can now be closed.

       close(unit=isunit)
       write(6,220) trim(afile)

! -- The pilot points parameterisation file is now read and, where necessary,
!    factors are generated for each pilot points file found therein.

! -- First the pilot point filename and interpolation factor filename are obtained.

       write(6,*)
       call addquote(pppfile,afile)
       write(6,350) trim(afile)
350    format(' - reading pilot points parameterization file ',a,'...')
       iblock=0
       iline=0

! -- We look for the start of a block.

       do
         iline=iline+1
         call num2char(iline,aline)
         read(pppunit,'(a)',err=9300,end=1000) cline
         if(cline.eq.' ') cycle
         cline=adjustl(cline)
         if(cline(1:1).eq.'#') cycle
         call linesplit(ifail,1)
         atemp=cline(left_word(1):right_word(1))
         call casetrans(atemp,'lo')
         if(atemp.ne.'start') then
           write(amessage,358) trim(aline),trim(afile)
358        format(' "START" keyword expected at line ',a,' of pilot ',  &
           'points parameterization file ',a,'.')
           go to 9890
         end if
         call linesplit(ifail,2)
         if(ifail.ne.0)then
           write(amessage,360) trim(aline),trim(afile)
360        format(' Block name must follow "START" keyword at line ',a,  &
           ' of file ',a,'.')
           go to 9890
         end if
         atemp=cline(left_word(2):right_word(2))
         call casetrans(atemp,'lo')
         if(atemp(1:11).ne.'spatial_par')then
           write(amessage,370) trim(aline),trim(afile)
370        format(' "SPATIAL_PARAMETER_DEFINITION" string must follow "START" ', &
           'keyword at line ',a,' of pilot points parameterization file ',a,'.')
           go to 9890
         end if

! -- We now read the block.

         iblock=iblock+1
         ipfile=0
         iproperty=0
         iffile=0
         iregion=0
         iub=0
         ilb=0
         iiv=0
         itr=0

         do
           iline=iline+1
           call num2char(iline,aline)
           read(pppunit,'(a)',err=9300,end=9350) cline
           if(cline.eq.' ') cycle
           cline=adjustl(cline)
           if(cline(1:1).eq.'#') cycle
           call linesplit(ifail,2)
           if(ifail.ne.0) go to 9150
           akey=cline(left_word(1):right_word(1))
           call casetrans(akey,'lo')
           ibeg=left_word(2)
           iend=len_trim(cline)

           if(akey.eq.'region')then
             if(iregion.ne.0)then
               write(amessage,379) trim(aline),trim(afile)
379            format(' Repeated REGION keyword at line ',a,' of file ',a,'.')
               go to 9890
             end if
             call getfile(ifail,cline,aregion,ibeg,iend)
             if(ifail.ne.0) go to 9300
             call casetrans(aregion,'lo')
             do i=1,nregion
               if(aregion.eq.aregname(i)) go to 390
             end do
             call addquote(rdffile,bfile)
             write(amessage,380) trim(aregion),trim(aline),trim(afile),trim(bfile)
380          format(' Region "',a,'" cited at line ',a,' of pilot points ',  &
             'parameterization file ',a,' is not cited in regions definition file ',a,'.')
             go to 9890
390          continue
             iregion=i

           else if(akey.eq.'property')then
             if(iproperty.ne.0)then
               write(amessage,391) trim(aline),trim(afile)
391            format(' Repeated PROPERTY keyword at line ',a,' of file ',a,'.')
               go to 9890
             end if
             call getfile(ifail,cline,aproperty,ibeg,iend)
             if(ifail.ne.0) go to 9300
             call casetrans(aproperty,'lo')
             do i=1,nproperties
               if(aproperty.eq.apropname(i)) go to 420
             end do
             call addquote(rpffile,bfile)
             write(amessage,400) trim(aproperty),trim(aline),trim(afile),trim(bfile)
400          format(' Property "',a,'" cited at line ',a,' of pilot points ',  &
             'parameterization file ',a,' is not cited in region property file ',a,'.')
             go to 9890
420          continue
             iproperty=i

           else if(akey.eq.'pilot_points_file')then
             if(ipfile.ne.0)then
               write(amessage,421) trim(aline),trim(afile)
421            format(' Repeated PILOT_POINTS_FILE keyword at line ',a,' of file ',a,'.')
               go to 9890
             end if
             call getfile(ifail,cline,pfile,ibeg,iend)
             if(ifail.ne.0) go to 9300
             ipfile=1

           else if(index(akey,'hor_interp').ne.0)then
             if(iffile.ne.0)then
               write(amessage,422) trim(aline),trim(afile)
422            format(' Repeated HORIZONTAL_INTERPOLATION_FACTOR_FILE keyword at line ',a, &
               ' of file ',a,'.')
               go to 9890
             end if
             call getfile(ifail,cline,ffile,ibeg,iend)
             if(ifail.ne.0) go to 9300
             iffile=1

           else if(akey.eq.'transform')then
             if(itr.ne.0)then
               write(amessage,241) trim(aline),trim(afile)
241            format(' Repeated TRANSFORM keyword at line ',a,        &
               ' of pilot points parameterization file ',a)
               go to 9890
             end if
             call getfile(ifail,cline,atrans,ibeg,iend)
             if(ifail.ne.0)then
               write(amessage,2422) trim(aline),trim(afile)
2422           format(' Cannot read TRANSFORM at line ',a,' of pilot points ', &
               'parameterization file ',a)
               go to 9890
             end if
             call casetrans(atrans,'lo')
             if(atrans.eq.'log')then
               ipptrans=1
             else if(atrans.eq.'none')then
               ipptrans=0
             else
               write(amessage,242) trim(aline),trim(afile)
242            format(' TRANSFORM type must be "log" or "none" at line ',a, &
               ' of pilot points parameterization file ',a)
               go to 9890
             end if
             itr=1

           else if(akey.eq.'initial_value')then
             if(iiv.ne.0)then
               write(amessage,243) trim(aline),trim(afile)
243            format(' Repeated INITIAL_VALUE keyword at line ',a,        &
               ' of pilot points parameterization file ',a)
               go to 9890
             end if
             iiv=1

           else if(akey.eq.'lower_bound')then
             if(ilb.ne.0)then
               write(amessage,245) trim(aline),trim(afile)
245            format(' Repeated LOWER_BOUND keyword at line ',a,        &
               ' of pilot points parameterization file ',a)
               go to 9890
             end if
             call getfile(ifail,cline,anum,ibeg,iend)
             if(ifail.ne.0)then
               write(amessage,246) trim(aline),trim(afile)
246            format(' Cannot read LOWER_BOUND from line ',a,' of ', &
               'pilot points parameterization file ',a)
               go to 9890
             else
               call char2num(ifail,anum,lppbound)
               if(ifail.ne.0)then
                 write(amessage,246) trim(aline),trim(afile)
                 go to 9890
               end if
             end if
             ilb=1

           else if(akey.eq.'upper_bound')then
             if(iub.ne.0)then
               write(amessage,247) trim(aline),trim(afile)
247            format(' Repeated UPPER_BOUND keyword at line ',a,        &
               ' of pilot points parameterization file ',a)
               go to 9890
             end if
             call getfile(ifail,cline,anum,ibeg,iend)
             if(ifail.ne.0)then
               write(amessage,248) trim(aline),trim(afile)
248            format(' Cannot read UPPER_BOUND from line ',a,' of ', &
               'pilot points parameterization file ',a)
               go to 9890
             else
               call char2num(ifail,anum,uppbound)
               if(ifail.ne.0)then
                 write(amessage,248) trim(aline),trim(afile)
                 go to 9890
               end if
             end if
             iub=1

           else if(akey.eq.'end')then

             atemp=cline(left_word(2):right_word(2))
             call casetrans(atemp,'lo')
             if(atemp(1:17).ne.'spatial_parameter')then
               write(amessage,440) trim(aline),trim(afile)
440            format(' The string "SPATIAL_PARAMETER_DEFINITION" should follow ',  &
               '"END" at line ',a,' of pilot points parameterization file ',a,'.')
               go to 9890
             end if
             if(ipfile.eq.0)then
               akey='PILOT_POINTS_FILE'
               go to 9400
             end if
             if(iffile.eq.0)then
               akey='HOR_INTERPOLATION_FACTOR_FILE'
               go to 9400
             end if
             if(iproperty.eq.0)then
               akey='PROPERTY'
               go to 9400
             end if
             if(iregion.eq.0)then
               akey='REGION'
               go to 9400
             end if
             if(ilb.eq.0)then
               akey='LOWER_BOUND'
               go to 9400
             end if
             if(iub.eq.0)then
               akey='UPPER_BOUND'
               go to 9400
             end if
             if(itr.eq.0)then
               akey='TRANSFORM'
               go to 9400
             end if

             if((ipptrans.eq.1).and.(lppbound.le.0.0d0))then
               write(amessage,361) trim(aline),trim(afile)
361            format(' LOWER_BOUND must be positive if TRANSFORM is "log" for ', &
               'block that ends at line ',a,' of pilot points parameterization ', &
               'file ',a)
               go to 9890
             end if
             if(lppbound.ge.uppbound)then
               write(amessage,382) trim(aline),trim(afile)
382            format(' UPPER_BOUND must be higher than LOWER_BOUND for ', &
              'block that ends at line ',a,' of pilot points parameterization ', &
              'file ',a)
              go to 9890
            end if

! -- A series of questions is now asked of the user.

             write(6,*)
             call addquote(pfile,bfile)
449          write(6,500,advance='no') trim(bfile)
500          format(' Generate interpolation factors for pilot points file ',a,  &
             '? (y/n): ')
             read(5,'(a)') aa
             if(aa.eq.' ') go to 449
             call casetrans(aa,'lo')
             if(aa.eq.'e')then
               write(6,510)
510            format(' Sorry - no backtracking allowed for this prompt - try again.')
               go to 449
             else if(aa.eq.'n')then
               write(6,512) trim(afile)
512            format(/,' - continuing to read pilot points parameterization file ',a,'...')
               go to 930
             else if(aa.eq.'y')then
               continue
             else
               write(6,513)
513            format(' Response must be "y" or "n" - try again.')
               go to 449
             end if
             inquire(file=pfile,exist=lexist)
             if(.not.lexist)then
               write(6,520) trim(bfile)
520            format(' Cannot find pilot points file ',a,' - try again.')
               go to 449
             end if

521          continue

             inquire(file=ffile,exist=lexist)
525          continue
             if(lexist)then
               call addquote(ffile,cfile)
523            write(6,522,advance='no') trim(cfile)
522            format(' Overwrite existing interpolation factor file ',a,'? (y/n): ')
               read(5,'(a)') aa
               if(aa.eq.' ') go to 523
               call casetrans(aa,'lo')
               if((aa.eq.'e').or.(aa.eq.'n'))then
                 write(6,*)
                 go to 449
               end if
               if((aa.ne.'y').and.(aa.ne.'n')) then
                 write(6,513)
                 go to 523
               end if
             end if

             ilog=ipptrans

550          write(6,560,advance='no')
560          format(' Base kriging on exponential/spherical/power variogram? (x/s/p) : ')
             read(5,'(a)') aa
             if(aa.eq.' ') go to 550
             call casetrans(aa,'lo')
             if(aa.eq.'e')then
               write(6,*)
               if(lexist)then
                 go to 525
               else
                go to 449
              end if
             end if
             if(aa.eq.'x')then
               ivar=1
             else if(aa.eq.'s')then
               ivar=2
             else if(aa.eq.'p')then
               ivar=3
             else
               go to 550
             end if

570          write(6,580,advance='no')
580          format(' Enter "a"                                                      : ')
             i=key_read(a)
             if(escset.ne.0)then
	       escset=0
               write(6,*)
               go to 550
             else if(i.eq.-1) then
               go to 570
             else if(i.ne.0) then
               write(6,590)
590            format(' Illegal input  - try again.')
               go to 570
             end if
             if(a.le.0.0d0)then
               write(6,600)
600            format(' Negative "a" not allowed - try again.')
               go to 570
             end if
             if((ivar.eq.3).and.(a.gt.2.0))then
               write(6,601)
601            format(' Exponent in power variogram must be 2.0 or less - try again.')
               go to 570
             end if

610          write(6,620,advance='no')
620          format(' Enter anisotropy                                               : ')
             i=key_read(anis)
             if(escset.ne.0)then
               escset=0
               write(6,*)
               go to 570
             else if(i.eq.-1) then
               go to 610
             else if(i.ne.0) then
               write(6,590)
               go to 610
             end if
             if(anis.lt.1.0d0)then
               write(6,630)
630            format(' Anisotropy must be 1.0 or greater - try again.')
               go to 610
             end if

640          continue
             if(anis.eq.1.0d0)then
               anis_angle=0.0d0
             else
650            write(6,660,advance='no')
660            format(' Enter bearing of major anisotropy axis                         : ')
               i=key_read(anis_angle)
               if(escset.ne.0)then
                 escset=0
                 write(6,*)
                 go to 610
               else if(i.eq.-1) then
                 go to 650
               else if(i.ne.0) then
                 write(6,590)
                 go to 650
               end if
               if((anis_angle.lt.0.0d0).or.(anis_angle.gt.180.0d0))then
                 write(6,150)
                 go to 650
               end if
             end if

670          write(6,680,advance='no')
680          format(' Enter pilot point search radius                                : ')
             i=key_read(radius)
             if(escset.ne.0)then
               escset=0
               write(6,*)
               if(anis.eq.1.0d0)then
                 go to 610
               else
                 go to 640
               end if
             else if(i.eq.-1) then
               go to 670
             else if(i.ne.0) then
               write(6,590)
               go to 670
             end if
             if(radius.le.0.0d0)then
               write(6,690)
690            format(' Search radius must be greater than zero - try again.')
               go to 670
             end if

700          write(6,710,advance='no')
710          format(' Enter maximum pilot points from which to interpolate           : ')
             i=key_read(maxpoints)
             if(escset.ne.0)then
               escset=0
               write(6,*)
               go to 670
             else if(i.eq.-1) then
               go to 700
             else if(i.ne.0) then
               write(6,590)
               go to 700
             end if
             if(maxpoints.lt.1)then
               write(6,720)
720            format(' Maximum pilot points must be greater than zero - try again.')
               go to 700
             end if

730          write(6,740,advance='no')
740          format(' Enter minimum pilot points from which to interpolate           : ')
             i=key_read(minpoints)
             if(escset.ne.0)then
               escset=0
               write(6,*)
               go to 700
             else if(i.eq.-1) then
               go to 730
             else if(i.ne.0) then
               write(6,590)
               go to 730
             end if
             if(minpoints.lt.1)then
               write(6,750)
750            format(' Minimum pilot points must be greater than zero - try again.')
               go to 730
             end if
             if(minpoints.gt.maxpoints)then
               write(6,760)
760            format(' Minimum pilot points must not exceed maximum ', &
               'pilot points - try again.')
               go to 730
             end if

! -- Screen data has been gathered. Now the pilot points file is read.

             write(6,*)
             call addquote(pfile,apfile)
             write(6,770) trim(apfile)
770          format(' - reading pilot points file ',a,'...')
             punit=nextunit()
             open(unit=punit,file=pfile,status='old',iostat=ierr)
             if(ierr.ne.0)then
               write(amessage,780) trim(apfile)
780            format(' Cannot open pilot points file ',a,'.')
               go to 9890
             end if
             jline=0
             ipp=0
             do
790            jline=jline+1
               read(punit,'(a)',err=9600,end=850) cline
               if(cline.eq.' ') go to 790
               call linesplit(ifail,3)
               if(ifail.ne.0) then
                 call num2char(jline,aline)
                 write(amessage,800)  trim(aline),trim(apfile)
800              format(' Three entries expected on line ',a,' of pilot points ', &
                 'file ',a,'.')
                 go to 9890
               end if
               ipp=ipp+1
               if(ipp.gt.MAXPP)then
                 write(amessage,810) trim(apfile)
810              format(' Too many pilot points in file ',a,'; increase MAXPP ',  &
                 'and re-compile program.')
                 go to 9890
               end if
               ep(ipp)=char2double(ifail,1)
               if(ifail.ne.0)go to 9600
               np(ipp)=char2double(ifail,2)
               if(ifail.ne.0)go to 9600
               valdat(ipp)=char2real(ifail,3)
               if(ifail.ne.0) go to 9600
               inumdat(ipp)=ipp
             end do
850          continue
             npp=ipp
             if(npp.eq.0)then
               write(amessage,860) trim(apfile)
860            format(' No pilot points were found in pilot points file ',a,'.')
               go to 9890
             end if
             close(unit=punit)
             call num2char(npp,aline)
             write(6,870) trim(aline),trim(apfile)
870          format(' - ',a,' pilot points read from file ',a,'.')

             write(6,880)
880          format(' - calculating kriging factors...')

! -- If this is our first loop we set up some arrays.

             iloop=iloop+1
             if(iloop.eq.1)then
               ncol=nx
               nrow=ny

! -- The pseudo-MODFLOW grid specifications are set.

               gridspec%nrow=nrow
               gridspec%ncol=ncol
               gridspec%east_corner=e_corner
               gridspec%north_corner=n_corner
               gridspec%rotation=angle
               gridspec%cosang=cos(gridspec%rotation*3.14159265/180.0)
               gridspec%sinang=sin(gridspec%rotation*3.14159265/180.0)
               allocate(gridspec%delr(1:gridspec%ncol),&
               gridspec%delc(1:gridspec%nrow),stat=ierr)
               if(ierr.ne.0) go to 9100
               do i=1,ncol
                 gridspec%delr(i)=dx(i)
               end do
               do i=1,nrow
                 gridspec%delc(i)=dy(i)
               end do
               gridspec%specfile=' '

! -- Eastings and norhings relative to top left of grid are evaluated for all
!    grid cells.

               call rel_centre_coords(east,north,gridspec)
               call grid2earth(east,north,gridspec)

             end if


! -- The coordinates of the pilot points have the grid left corner easting and
!    northing subtracted and are placed into single precision arrays.

             do i=1,npp
               eastdat(i)=ep(i)-gridspec%east_corner
               northdat(i)=np(i)-gridspec%north_corner
             end do

! -- For the specified region we now accumulate cell centre coordinates.

             ngrid=0
             do iy=1,ny
               do ix=1,nx
                 do iz=1,nz
                   if(igrid(ix,iy,iz).eq.iregion)then
                     ngrid=ngrid+1
                     eastgrid(ngrid)=east(ix,iy)
                     northgrid(ngrid)=north(ix,iy)
                     call rc2cell(icellno(ngrid),iy,ix,gridspec)
                     go to 890
                   end if
                 end do
890              continue
               end do
             end do
             if(ngrid.eq.0) then
               call addquote(pppfile,afile)
               write(amessage,900) trim(aregname(iregion)),trim(afile)
900            format(' No grid cells belong to region "',a,'" cited in ',  &
               'SPATIAL_PARAMETER_DEFINITION block of pilot points ',      &
               'parameterization file ',a,'.')
               go to 9890
             end if

! -- We now set up for ordinary kriging.

             k_ktype=1
             s_skmean=0.0
             n_nst=1
             c_c0=0.0
             pmx=1.0e4
             itrans=ilog       !check
             if(ivar.eq.2)then
               i_it(1)=1
             else if(ivar.eq.1)then
               i_it(1)=2
             else if(ivar.eq.3)then
               i_it(1)=4
             end if
             c_cc(1)=1.0
             a_ang(1)=anis_angle
             a_aa(1)=a
             a_anis(1)=1.0/anis
             minpt=minpoints
             maxpt=maxpoints
             radmax=radius
             ndat=npp
             aoutfile='f'

! -- The kriging factor file is now opened.

             outunit=nextunit()
             call addquote(ffile,affile)
             open(unit=outunit,file=ffile,iostat=ierr)
             if(ierr.ne.0)then
               write(amessage,915) trim(affile)
915            format(' Cannot open interpolation factor file ',a,' for writing.')
               go to 9890
             end if

! -- The first line of the file is written.

             write(outunit,916) npp,ngrid,itrans,lppbound,uppbound
916          format(3i6,1x,1pg14.7,1x,1pg14.7)

! -- Krigging is carried out.

             call readparm(minpt,maxpt,radmax,k_ktype,s_skmean, &
             n_nst,c_c0,i_it,c_cc,a_ang,a_aa,a_anis,ndat,eastdat,northdat,valdat)
             call kb2d(ngrid,ndat,ncol,nrow,inumdat,icellno,eastgrid,northgrid, &
             ffile,aoutfile,outunit,variance,pmx,itrans,noest)

             close(unit=outunit)
             if(noest.eq.0)then
               write(6,920) trim(affile)
920            format(' - file ',a,' written ok.')
             else
               call num2char(noest,aline)
               call addquote(pfile,afile)
               write(amessage,922) trim(aline),trim(aregname(iregion)),trim(afile)
922            format(' Kriging factors could not be calculated for ',a,' cells ',  &
               'in region "',a,'" based on pilot points in file ',a,'; increase ', &
               'search radius and/or decrease minimum ', &
               'interpolation points setting.')
               go to 9890
             end if
             go to 930

           else
             call casetrans(akey,'hi')
             write(amessage,381) trim(akey),trim(aline),trim(afile)
381          format(' Unrecognised keyword "',a,'" at line ',a,' of pilot ',  &
             'points parameterization file ',a,'.')
             go to 9890
           end if

         end do
930      continue

       end do


1000   continue
       close(unit=pppunit)

       if(iblock.eq.0)then
         write(amessage,1010) trim(afile)
1010     format(' No SPATIAL_PARAMETER_DEFINITION blocks were found in pilot ', &
         'points parameterization file ',a,'.')
         go to 9890
       else
!         call num2char(iblock,aline)
!         write(amessage,1020) trim(aline),trim(afile)
!1020     format(' - ',a,' blocks read from file ',a,'.')
!         call write_message()
       end if


       go to 9900


9000    continue
        call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of regions definition file ',a,'.')
        go to 9890
9050    continue
        write(amessage,9060) trim(afile)
9060    format(' Unexpected end encountered to regions defintion file ', &
        a,'.')
        go to 9890
9100    write(amessage,9110)
9110    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890
9150    call num2char(iline,aline)
        write(amessage,9160) trim(aline),trim(afile)
9160    format(' Insufficient entries on line ',a,' of file ',a,'.')
        go to 9890
9200    write(amessage,9210) trim(afile)
9210    format(' Error encountered while reading IS file ',a,'.')
        go to 9890
9250    write(amessage,9260) trim(afile)
9260    format(' Unexpected end encountered to IS file ',a,'.')
        go to 9890
9300    write(amessage,9310) trim(aline),trim(afile)
9310    format(' Error encountered while reading line ',a,   &
        ' of pilot points parameterization file ',a,'.')
        go to 9890
9350    write(amessage,9360) trim(afile)
9360    format(' Unexpected end encountered to pilot points parameterization ', &
        'file ',a,' before end of SPATIAL_PARAMETER_DEFINITION block.')
        go to 9890
9400    write(amessage,9410) trim(akey),trim(aline),trim(afile)
9410    format(' No ',a,' keyword found ',         &
        'in SPATIAL_PARAMETER_DEFINITION block that ends ',   &
        'at line ',a,' of file ',a)
        go to 9890
9500    write(amessage,9510) trim(aline),trim(afile)
9510    format(' Error reading line ',a,' of region properties file ',a,'.')
        go to 9890
9550    write(amessage,9560) trim(afile)
9560    format(' Unexpected end encountered to region properties file ',a,'.')
        go to 9890
9600    call num2char(jline,aline)
        write(amessage,9610) trim(aline),trim(apfile)
9610    format(' Error reading line ',a,' of pilot points file ',a,'.')
        go to 9890

9890	call write_message(leadspace='yes')

9900    continue

        call close_files
        call free_grid_mem(gridspec)
        if(allocated(icellno)) deallocate(icellno,stat=ierr)
        if(allocated(igrid)) deallocate(igrid,stat=ierr)
        if(allocated(east)) deallocate(east,stat=ierr)
        if(allocated(north)) deallocate(north,stat=ierr)
        if(allocated(eastgrid)) deallocate(eastgrid,stat=ierr)
        if(allocated(northgrid)) deallocate(northgrid,stat=ierr)
        if(allocated(variance)) deallocate(variance,stat=ierr)
        if(allocated(aregname)) deallocate(aregname,stat=ierr)
        if(allocated(apropname)) deallocate(apropname,stat=ierr)
        if(allocated(dx)) deallocate(dx,stat=ierr)
        if(allocated(dy)) deallocate(dy,stat=ierr)

end program facgen


 
      subroutine readparm(n_ndmin,n_ndmax,r_radius,k_ktype,s_skmean, &
      n_nst,c_c0,i_it,c_cc,a_ang,a_aa,a_anis,n_nd,eastdat,northdat,valdat)


      implicit integer(i-n), real(a-h, o-z)                    !jd

!-----------------------------------------------------------------------
!
!                  Initialization and Read Parameters
!                  **********************************
!
! The input parameters and data are read in from their files. Some quick
! error checking is performed and the statistics of all the variables
! being considered are written to standard output.
!
!
!
!-----------------------------------------------------------------------
      include  'kb2d.inc'
      parameter(MV=20)
      real      var(MV)
      character datafl*40,outfl*40,dbgfl*40,str*40
      logical   testfl

! -- Dimension statements added by myself for subroutine arguments.

      integer i_it(n_nst)
      real c_cc(n_nst),a_ang(n_nst),a_aa(n_nst),a_anis(n_nst), &
      valdat(n_nd)
      real eastdat(n_nd),northdat(n_nd)
!
! Unit numbers:
!
      lout = 2
      ldbg = 3

!
! Read Input Parameters:
!

! -- Parameters that should have been read from input file.

      datafl=' '
      ixl=-9999
      iyl=-9999
      ivrl=-9999
      tmin=-1.1e35
      tmax=1.1e35
      idbg=0                    !debugging level
      dbgfl='debug1.dat'
      outfl='kb2d.out'
      nxdis=1
      nydis=1


! -- Here are the dummy grid specifications.

      nx=1
      xmn=0.0
      xsiz=1.0
      ny=1
      ymn=0.0
      ysiz=1.0


! -- The following variables were supplied through the subroutine argument.

      ndmin=n_ndmin
      ndmax=n_ndmax
      radius=r_radius
      ktype=k_ktype
      skmean=s_skmean
      nst=n_nst
      c0=c_c0

            do i=1,nst
                  it(i)=i_it(i)
                  cc(i)=c_cc(i)
                  ang(i)=a_ang(i)
                  aa(i)=a_aa(i)
!                  a2=a_a2(i)
!                  anis(i) = a2 / aa(i)
                  anis(i)=a_anis(i)
                  if(it(i).eq.4) then
                        if(aa(i).lt.0.0) stop ' INVALID power variogram'
                        if(aa(i).gt.2.0) stop ' INVALID power variogram'
                  end if
            end do

!      write(*,*)
      if(nst.gt.MAXNST)   stop ' nst is too big - contact programmer'
      if(ndmax.gt.MAXSAM) stop ' ndmax is too big - contact programmer'

      av = 0.0
      ss = 0.0
      nd = n_nd

      if(nd.gt.MAXDAT) then
            write(*,*) ' nd too big - contact programmer'
            stop
      end if
      do i=1,nd
        x(i)  = eastdat(i)
        y(i)  = northdat(i)
        vr(i) = valdat(i)
        av     = av + vr(i)
        ss     = ss + vr(i)*vr(i)
      end do

!
! Open the output files:
!
!      open(lout,file=outfl,status='UNKNOWN')
!      write(lout,300)
! 300  format('KB2D Output',/,'2',/,'Estimate',/,'Estimation Variance')
      if(idbg.gt.0)open(ldbg,file=dbgfl,status='UNKNOWN')
!
! Compute the averages and variances as an error check for the user:
!
      av = av / max(real(nd),1.0)
      ss =(ss / max(real(nd),1.0)) - av * av
!
! Write Some of the Statistics to the screen:
!

!!!       write(6,*)
!!!       write(6,900) nd
!!!900    format('   Number of pilot points for this region   = ',i5)
!!!       write(6,901) av
!!!901    format('   Mean data value for these pilot points   = ',1pg12.5)
!!!       write(6,902) sqrt(max(ss,0.0))
!!!902    format('   Data standard deviation for these points = ',1pg12.5)
!!!       write(6,903)
!!!903    format('   Working....')
!      write(*,900) nd,av,sqrt(max(ss,0.0))
! 900  format(/' There are ',i8,' data with:',/,     &
!              '   mean value          = ',f12.5,/,  &
!              '   standard deviation  = ',f12.5,/)
      return
      end



      subroutine kb2d(n_npts,n_ndat,n_ncol,n_nrow,inumdat,icellno,epoint,npoint,outfile, &
      aoutfile,outunit,var_arr,pmx,itrans,i_noest)

      implicit integer(i-n), real(a-h, o-z)                    !jd

!-----------------------------------------------------------------------
!
!           Ordinary/Simple Kriging of a 2-D Rectangular Grid
!           *************************************************
!
! This subroutine estimates point or block values of one variable by
! ordinary kriging.  All of the samples are rescanned for each block
! estimate; this makes the program simple but inefficient.  The data
! should NOT contain any missing values.  Unestimated points are
! returned as -1.0e21
!
!
!
! Original:  A.G. Journel                                           1978
! Revisions: B.E. Buxton                                       Apr. 1983
!-----------------------------------------------------------------------
      include  'kb2d.inc'
      real      xdb(MAXDIS),ydb(MAXDIS),xa(MAXSAM),ya(MAXSAM), &
                vra(MAXSAM),dist(MAXSAM)
      real*8    r(MAXSAM+1),rr(MAXSAM+1),s(MAXSAM+1),a(MAXKRG)
      integer   nums(MAXSAM)
      logical   first

! -- The following dimension statement was added by myself for the
!    subroutine arguments.

     integer n_npts,n_ndat,outunit,n_ncol,n_nrow,itrans             !jd
     integer i_noest                                                !jd
     integer icellno(n_npts),inumdat(n_ndat)                        !jd
     real pmx                                                       !jd
     real epoint(n_npts),npoint(n_npts)                             !jd
     real var_arr(n_ncol,n_nrow)                                    !jd
!jd     real regarray(npp,npp)                                         !jd
     character*(*) outfile,aoutfile                                 !jd

!jd      data      first/.true./,PMX/9999.0/

      first=.true.                                    !jd
      i_noest=0                                       !jd
!
! Echo the input parameters if debugging flag is >2:
!
      if(idbg.gt.2) then
            write(ldbg,*) 'KB2D Parameters'
            write(ldbg,*)
            write(ldbg,*) 'Variogram Parameters for ',nst,' structures:'
            write(ldbg,*) '  Nugget effect:         ',c0
            write(ldbg,*) '  Types of variograms:   ',(it(i),i=1,nst)
            write(ldbg,*) '  Contribution cc        ',(cc(i),i=1,nst)
            write(ldbg,*) '  Ranges:                ',(aa(i),i=1,nst)
            write(ldbg,*) '  Angle for Continuity:  ',(ang(i),i=1,nst)
            write(ldbg,*) '  Anisotropy Factors:    ',(anis(i),i=1,nst)
            write(ldbg,*) ' '
            write(ldbg,*) 'Grid for Kriging:'
            write(ldbg,*) '  Number of X and Y Blocks:',nx,ny
            write(ldbg,*) '  Origin of X and Y Blocks:',xmn,ymn
            write(ldbg,*) '  Size   of X and Y Blocks:',xsiz,ysiz
            write(ldbg,*) ' '
            write(ldbg,*) 'Discretization of blocks:  ',nxdis,nydis
            write(ldbg,*) 'Search Radius:             ',radius
            write(ldbg,*) 'Minimum number of samples: ',ndmin
            write(ldbg,*) 'Maximum number of samples: ',ndmax
            write(ldbg,*) ' '
      endif
!
! Echo the input data if debugging flag >1:
!
      if(idbg.ge.4) then
            do id=1,nd
                  write(ldbg,99) id,x(id),y(id),vr(id)
 99               format('Data: ',i5,' at ',2f12.3,' value: ',f12.5)
            end do
      endif
!
! Set up the discretization points per block.  Figure out how many
! are needed, the spacing, and fill the xdb and ydb arrays with the
! offsets relative to the block center (this only gets done once):
!
      ndb  = nxdis * nydis
      if(ndb.gt.MAXDIS) then
            write(*,*) ' ERROR KB2D: Too many discretization points '
            write(*,*) '             Increase MAXDIS or lower n[xy]dis'
            stop
      endif
      xdis = xsiz  / max(real(nxdis),1.0)
      ydis = ysiz  / max(real(nydis),1.0)
      xloc = -0.5*(xsiz+xdis)
      i    = 0
      do ix =1,nxdis
            xloc = xloc + xdis
            yloc = -0.5*(ysiz+ydis)
            do iy=1,nydis
                  yloc = yloc + ydis
                  i = i+1
                  xdb(i) = xloc
                  ydb(i) = yloc
            end do
      end do
!
! Initialize accumulators:
!
      cbb  = 0.0
      rad2 = radius*radius
!
! Calculate Block Covariance. Check for point kriging.
!
      cov   = cova2(xdb(1),ydb(1),xdb(1),ydb(1),nst,c0,PMX,cc, &
                    aa,it,ang,anis,first,passmaxcov)
!
! Keep this value to use for the unbiasedness constraint:
!
      unbias = cov
      first  = .false.
      if (ndb.le.1) then
            cbb = cov
      else
            do i=1,ndb
                  do j=1,ndb
                        cov = cova2(xdb(i),ydb(i),xdb(j),ydb(j),nst,c0, &
                                    PMX,cc,aa,it,ang,anis,first,passmaxcov)
                        if(i.eq.j) cov = cov - c0
                        cbb = cbb + cov
                  end do
            end do
            cbb = cbb/real(ndb*ndb)
      endif
      if(idbg.gt.1) then
            write(ldbg,*) ' '
            write(ldbg,*) 'Block Covariance: ',cbb
            write(ldbg,*) ' '
      endif
!
! MAIN LOOP OVER ALL THE BLOCKS IN THE GRID:
!
      nk = 0
      ak = 0.0
      vk = 0.0
!jd      do 4 iy=1,ny                                               !jd
!jd      yloc = ymn + (iy-1)*ysiz                                   !jd
!jd      do 4 ix=1,nx                                               !jd
!jd            xloc = xmn + (ix-1)*xsiz                             !jd

         do 4 i_ipts=1,n_npts                                       !jd
            xloc=epoint(i_ipts)                                     !jd
            yloc=npoint(i_ipts)                                     !jd
            icell=icellno(i_ipts)                                   !jd
!
! Find the nearest samples within each octant: First initialize
! the counter arrays:
!
            na = 0
            do isam=1,ndmax
                  dist(isam) = 1.0e+20
                  nums(isam) = 0
            end do
!
! Scan all the samples (this is inefficient and the user with lots of
! data should move to ktb3d):
!
            do 6 id=1,nd
                  dx = x(id) - xloc
                  dy = y(id) - yloc
                  h2 = dx*dx + dy*dy
                  if(h2.gt.rad2) go to 6
!
! Do not consider this sample if there are enough close ones:
!
                  if(na.eq.ndmax.and.h2.gt.dist(na)) go to 6
!
! Consider this sample (it will be added in the correct location):
!
                  if(na.lt.ndmax) na = na + 1
                  nums(na)           = id
                  dist(na)           = h2
                  if(na.eq.1) go to 6
!
! Sort samples found thus far in increasing order of distance:
!
                  n1 = na-1
                  do ii=1,n1
                        k=ii
                        if(h2.lt.dist(ii)) then
                              jk = 0
                              do jj=k,n1
                                    j  = n1-jk
                                    jk = jk+1
                                    j1 = j+1
                                    dist(j1) = dist(j)
                                    nums(j1) = nums(j)
                              end do
                              dist(k) = h2
                              nums(k) = id
                              go to 6
                        endif
                  end do
 6          continue
!
! Is there enough samples?
!
            if(na.lt.ndmin) then
                  i_noest=i_noest+1
                  if(idbg.ge.2) &
                  write(ldbg,*) 'Block ',ix,iy, 'not estimated'
                  est  = UNEST
                  estv = UNEST
                  go to 1
            endif
!
! Put coordinates and values of neighborhood samples into xa,ya,vra:
!
            do ia=1,na
                  jj      = nums(ia)
                  xa(ia)  = x(jj)
                  ya(ia)  = y(jj)
                  vra(ia) = vr(jj)
            end do
!
! Handle the situation of only one sample:
!
            if(na.eq.1) then
                  cb1 = cova2(xa(1),ya(1),xa(1),ya(1),nst,c0,   &
                              PMX,cc,aa,it,ang,anis,first,passmaxcov)
                  xx  = xa(1) - xloc
                  yy  = ya(1) - yloc
!
! Establish Right Hand Side Covariance:
!
                  if(ndb.le.1) then
                        cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0,  &
                                   PMX,cc,aa,it,ang,anis,first,passmaxcov)
                  else
                        cb  = 0.0
                        do i=1,ndb
                              cb = cb + cova2(xx,yy,xdb(i),ydb(i),nst, &
                                        c0,PMX,cc,aa,it,ang,anis,first,passmaxcov)
                              dx = xx - xdb(i)
                              dy = yy - ydb(i)
                              if((dx*dx+dy*dy).lt.EPSLON) &
                              cb = cb - c0
                        end do
                        cb = cb / real(ndb)
                  end if
                  if(ktype.eq.0) then
                        s(1) = cb/cbb
                        est  = s(1)*vra(1) + (1.0-s(1))*skmean
                        estv = cbb - s(1) * cb
                  else
                        est  = vra(1)
                        estv = cbb - 2.0*cb + cb1
                  end if

                  if(ktype.eq.0)then                                    !jd
                    rrtemp=(1.0-s(1))*skmean                            !jd
                  else                                                  !jd
                    rrtemp=0.0                                          !jd
                  end if                                                !jd
                  if(aoutfile.eq.'f')then                               !jd
                    write(outunit,*) icell,1, &                         !jd
                    inumdat(nums(1)),1.0                                  !jd
                  else                                                    !jd
                    write(outunit)   icell,1, &                           !jd
                    inumdat(nums(1)),1.0                                  !jd
                  end if                                                  !jd

            else
!
! Solve the Kriging System with more than one sample:
!
                  neq = na + ktype
                  nn  = (neq + 1)*neq/2
!
! Set up kriging matrices:
!
                  in=0
                  do j=1,na
!
! Establish Left Hand Side Covariance Matrix:
!
                        do i=1,j
                              in = in + 1
                              a(in) = dble( cova2(xa(i),ya(i),xa(j),  &
                                            ya(j),nst,c0,PMX,cc,aa,it, &
                                            ang,anis,first,passmaxcov) )
                              iitmp1=inumdat(nums(i))                           !jd
                              iitmp2=inumdat(nums(j))                           !jd
!jd                              if(regarray(iitmp1,iitmp2).lt.-1.0e35)then        !jd
!jd                                regarray(iitmp1,iitmp2)=passmaxcov -a(in)       !jd
!jd                                regarray(iitmp2,iitmp1)=regarray(iitmp1,iitmp2) !jd
!jd                              end if
                        end do
                        xx = xa(j) - xloc
                        yy = ya(j) - yloc
!
! Establish Right Hand Side Covariance:
!
                        if(ndb.le.1) then
                              cb = cova2(xx,yy,xdb(1),ydb(1),nst,c0,  &
                                         PMX,cc,aa,it,ang,anis,first,passmaxcov)
                        else
                              cb  = 0.0
                              do j1=1,ndb
                                    cb = cb + cova2(xx,yy,xdb(j1),  &
                                         ydb(j1),nst,c0,PMX,cc,aa,  &
                                         it,ang,anis,first,passmaxcov)
                                    dx = xx - xdb(j1)
                                    dy = yy - ydb(j1)
                                    if((dx*dx+dy*dy).lt.EPSLON)  &
                                          cb = cb - c0
                              end do
                              cb = cb / real(ndb)
                        end if
                        r(j)  = dble(cb)
                        rr(j) = r(j)
                  end do
!
! Set the unbiasedness constraint:
!
                  if(ktype.eq.1) then
                        do i=1,na
                              in    = in + 1
                              a(in) = dble(unbias)
                        end do
                        in      = in + 1
                        a(in)   = 0.0
                        r(neq)  = dble(unbias)
                        rr(neq) = r(neq)
                  end if
!
! Write out the kriging Matrix if Seriously Debugging:
!
                  if(idbg.ge.3) then
                        write(ldbg,101) ix,iy
                        is = 1
                        do i=1,neq
                              ie = is + i - 1
                              write(ldbg,102) i,r(i),(a(j),j=is,ie)
                              is = is + i
                        end do
 101                    format(/,'Kriging Matrices for Node: ',2i4,  &
                                 ' RHS first')
 102                    format('    r(',i2,') =',f7.4,'  a= ',9(10f7.4))
                  endif
!
! Solve the Kriging System:
!
                  call ksol(1,neq,1,a,r,s,ising)
!
! Write a warning if the matrix is singular:
!
                  if(ising.ne.0) then
        write(6,450) icell                                            !jd
450     format(' WARNING: singular kriging matrix for cell',i5)      !jd
                        est  = UNEST
                        estv = UNEST
                        go to 1
                  endif
!
! Write the kriging weights and data if requested:
!
                  if(idbg.ge.2) then
                        write(ldbg,*) '       '
                        write(ldbg,*) 'BLOCK: ',ix,iy
                        write(ldbg,*) '       '
                        if(ktype.eq.1) write(ldbg,*)    &
                        '  Lagrange multiplier: ',s(neq)*unbias
                        write(ldbg,*) '  BLOCK EST: x,y,vr,wt '
                        do i=1,na
                        write(ldbg,'(4f12.3)') xa(i),ya(i),vra(i),s(i)
                        end do
                  endif
!
! Compute the estimate and the kriging variance:
!
                  est  = 0.0
                  estv = cbb
                  sumw = 0.0
                  if(ktype.eq.1) estv = estv - real(s(na+1))
                  do i=1,na
                        sumw = sumw + real(s(i))
                        est  = est  + real(s(i))*vra(i)
                        estv = estv - real(s(i)*rr(i))
                  end do
                  if(ktype.eq.0) est = est + (1.0-sumw)*skmean

                  if(ktype.eq.0)then                                    !jd
                    rrtemp=(1.0-sumw)*skmean                            !jd
                  else                                                  !jd
                    rrtemp=0.0                                          !jd
                  end if                                                !jd
                  if(aoutfile.eq.'f')then                               !jd
                    write(outunit,*) icell,na, &                        !jd
                    ((inumdat(nums(i)),real(s(i))),i=1,na)                 !jd
                  else                                                     !jd
                    write(outunit)   icell,na, &                           !jd
                    ((inumdat(nums(i)),real(s(i))),i=1,na)                 !jd
                  end if                                                   !jd

            endif
            if(idbg.ge.2) then
                  write(ldbg,*) '  est  ',est
                  write(ldbg,*) '  estv ',estv
                  write(ldbg,*) ' '
            endif
!
! Write the result to the output file:
!
! 1          write(lout,'(f8.3,1x,f8.3)') est,estv

1           continue

            if(est.gt.UNEST)then
              iirow=(icell-1)/n_ncol+1
              iicol=icell-((iirow-1)*n_ncol)
              rrtemp=estv
              if(rrtemp.lt.0.0)rrtemp=0.0
              var_arr(iicol,iirow)=sqrt(rrtemp)
            end if

            if(est.gt.UNEST) then
                  nk = nk + 1
                  ak = ak + est
                  vk = vk + est*est
            end if
!
! END OF MAIN LOOP OVER ALL THE BLOCKS:
!
 4    continue
      if(nk.ge.1) then
            ak = ak / real(nk)
            vk = vk/real(nk) - ak*ak
!!!             write(6,105) nk
!!!105          format('   No. of grid points for which factors were calculated = ',i5)
!             write(6,106) ak
!106          format('   Average interpolated value at these points           = ',1pg12.5)
!             if(vk.lt.0.0)vk=0.0
!             write(6,107) sqrt(vk)
!107          format('   Standard deviation of interpolated grid point values = ',1pg12.5)

!            write(ldbg,105) nk,ak,vk
!            write(*,   105) nk,ak,vk
! 105        format(/,'Estimated   ',i8,' blocks ',/,  &
!                     '  average   ',f9.4,/,'  variance  ',f9.4,/)
      else
!!!        write(6,105) nk
      end if
      return
      end
 
 
 
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





      subroutine ksol(nright,neq,nsb,a,r,s,ising)

!-----------------------------------------------------------------------
!
!                Solution of a System of Linear Equations
!                ****************************************
!
!
!
! INPUT VARIABLES:
!
!   nright,nsb       number of columns in right hand side matrix.
!                      for KB2D: nright=1, nsb=1
!   neq              number of equations
!   a()              upper triangular left hand side matrix (stored 
!                      columnwise)
!   r()              right hand side matrix (stored columnwise)
!                      for kb2d, one column per variable
!
!
!
! OUTPUT VARIABLES:
!
!   s()              solution array, same dimension as  r  above.
!   ising            singularity indicator
!                      0,  no singularity problem
!                     -1,  neq .le. 1
!                      k,  a null pivot appeared at the kth iteration
!
!
!
! PROGRAM NOTES:
!
!   1. Requires the upper triangular left hand side matrix.
!   2. Pivots are on the diagonal.
!   3. Does not search for max. element for pivot.
!   4. Several right hand side matrices possible.
!   5. USE for ok and sk only, NOT for UK.
!
!
!-----------------------------------------------------------------------
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      real*8   a(*),r(*),s(*)
!
! If there is only one equation then set ising and return:
!
      if(neq.le.1) then
            ising = -1
            return
      endif
!
! Initialize:
!
      tol   = 0.1e-06
      ising = 0
      nn    = neq*(neq+1)/2
      nm    = nsb*neq
      m1    = neq-1
      kk    = 0
!
! Start triangulation:
!
      do k=1,m1
            kk=kk+k
            ak=a(kk)
            if(abs(ak).lt.tol) then
                  ising=k
                  return
            endif
            km1=k-1
            do iv=1,nright
                  nm1=nm*(iv-1)
                  ii=kk+nn*(iv-1)
                  piv=1./a(ii)
                  lp=0
                  do i=k,m1
                        ll=ii
                        ii=ii+i
                        ap=a(ii)*piv
                        lp=lp+1
                        ij=ii-km1
                        do j=i,m1
                              ij=ij+j
                              ll=ll+j
                              a(ij)=a(ij)-ap*a(ll)
                        end do
                        do llb=k,nm,neq
                              in=llb+lp+nm1
                              ll1=llb+nm1
                              r(in)=r(in)-ap*r(ll1)
                        end do
                  end do
            end do
      end do
!
! Error checking - singular matrix:
!
      ijm=ij-nn*(nright-1)
      if(abs(a(ijm)).lt.tol) then
            ising=neq
            return
      endif
!
! Finished triangulation, start solving back:
!
      do iv=1,nright
            nm1=nm*(iv-1)
            ij=ijm+nn*(iv-1)
            piv=1./a(ij)
            do llb=neq,nm,neq
                  ll1=llb+nm1
                  s(ll1)=r(ll1)*piv
            end do
            i=neq
            kk=ij
            do ii=1,m1
                  kk=kk-i
                  piv=1./a(kk)
                  i=i-1
                  do llb=i,nm,neq
                        ll1=llb+nm1
                        in=ll1
                        ap=r(in)
                        ij=kk
                        do j=i,m1
                              ij=ij+j
                              in=in+1
                              ap=ap-a(ij)*s(in)
                        end do
                        s(ll1)=ap*piv
                  end do
            end do
      end do
!
! Finished solving back, return:
!
      return
      end


! -- Check variable grid spacing
 