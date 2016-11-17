!     Last change:  J     2 May 2002    8:49 pm
program conc2elev

! -- Program CONC2ELEV determines the location of the "freshwater/saltwater interface"
!    from a sequence of concentration arrays. It provides the position of this interface
!    in a number of different formats.

        use defn
        use inter

	implicit none

        integer, parameter   :: MAXCROSS=100
        integer              :: nbb,iunit,n,ifail,idate,iheader,i,ierr,icross
        integer              :: ucnunit,outunit
        integer              :: iarray,icount
        integer              :: ilay,irow,icol,ncol,nrow,ilay1,ilay2,mrow,mcol
        integer              :: ntrans,kstp,kper
        integer              :: itype(MAXCROSS)
        integer, allocatable :: filled(:)

        real                 :: thresh_interface,thresh_inactive,totim,dummy1,dummy2,rtemp
        real                 :: elev1,elev2,elev3,mid1,mid2,conc1,conc2,frac
        real                 :: dist(MAXCROSS),eastdist(MAXCROSS),northdist(MAXCROSS),z(MAXCROSS)
        real, allocatable    :: bottom(:,:,:),conc(:,:,:),elev(:,:)
        real, allocatable    :: east(:,:),north(:,:),x(:,:),y(:,:)

        character*1          :: adm
        character*10         :: alay1,alay2,atemp,alay,anum
        character*16         :: text
        character*100        :: aprompt
        character*200        :: afile,basename,infile,ucnfile,outfile,bfile

	type (modelgrid) gridspec


	write(amessage,5)
5       format(' Program CONC2ELEV computes the elevation of the saltwater/freshwater ', &
        'interface from a set of set of layer-specific concentration arrays read from ', &
        'a MT3D UCN file.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'
	
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
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

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

        write(6,*)
30      write(6,40,advance='no')
40      format(' Enter upper layer for processing: ')
	i=key_read(ilay1)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 10
	else if(i.eq.-1) then
	  go to 30
	else if(i.ne.0) then
	  write(6,50)
50	  format(' Illegal input  - try again.')
	  go to 30
        else if(ilay1.le.0) then
          write(6,50)
          go to 30
	end if
60      write(6,70,advance='no')
70      format(' Enter lower layer for processing: ')
	i=key_read(ilay2)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 30
	else if(i.eq.-1) then
	  go to 60
	else if(i.ne.0) then
	  write(6,50)
	  go to 60
	end if
        if(ilay2.lt.ilay1)then
          write(6,71)
71        format(' *** Must be greater than upper layer - try again ***')
          go to 60
        end if

        write(6,*)
80      call num2char(ilay1-1,alay1)
        call num2char(ilay2,alay2)
        write(6,90,advance='no') trim(alay1),trim(alay2)
90      format(' Enter filename base for layer ',a,' to ',a,' bottom elev arrays: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 80
        atemp=afile(1:2)
        call casetrans(atemp,'lo')
        if(atemp(1:2).eq.'e ') then
          write(6,*)
          go to 60
        end if
        nbb=len_trim(afile)
        call getfile(ifail,afile,basename,1,nbb)
        if(ifail.ne.0) go to 80
        allocate(bottom(ncol,nrow,ilay1-1:ilay2),stat=ierr)
        if(ierr.ne.0) go to 9200
        do ilay=ilay1-1,ilay2
          call num2char(ilay,alay)
          infile=trim(basename)//trim(alay)//'.ref'
          iunit=nextunit()
          open(unit=iunit,file=infile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            call addquote(infile,afile)
            write(amessage,110) trim(afile)
110         format(' Cannot open file ',a,'.')
            go to 9890
          end if
          if(headerspec.eq.'yes')then
            read(iunit,*,err=9000,end=9050) mcol,mrow
            if((mrow.ne.nrow).or.(mcol.ne.ncol))then
              call addquote(infile,afile)
              call addquote(gridspec%specfile,bfile)
              write(amessage,111) trim(afile),trim(bfile)
111           format(' Number of columns or number of rows in header to array file ',a, &
              ' is not the same as that in finite-difference grid as read from from ', &
              'grid specification file ',a,'.')
              go to 9890
            end if
          end if
          do irow=1,nrow
            read(iunit,*,err=9000,end=9050) (bottom(icol,irow,ilay),icol=1,ncol)
          end do
          close(unit=iunit)
          write(6,120) trim(infile)
120       format(' - file ',a,' read ok.')
        end do

        write(6,*)
150     continue
	call open_input_file(ifail, &
	' Enter name of unformatted MT3D concentration file: ',ucnfile,ucnunit, &
	file_format='unformatted')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  deallocate(bottom)
	  go to 80
	end if
180     write(6,190,advance='no')
190     format(' Enter simulation time to read arrays for: ')
	i=key_read(totim)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
          close(unit=ucnunit)
	  go to 150
	else if(i.eq.-1) then
	  go to 180
	else if(i.ne.0) then
	  write(6,195)
195	  format(' Illegal input  - try again.')
	  go to 180
        else if(totim.lt.0.0)then
          write(6,195)
          go to 180
	end if
210     write(6,220,advance='no')
220     format(' Enter threshold concentration defining interface: ')
	i=key_read(thresh_interface)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
          go to 180
	else if(i.eq.-1) then
	  go to 210
	else if(i.ne.0) then
	  write(6,195)
	  go to 210
	end if
        if(thresh_interface.le.0.0)then
          write(6,221)
221       format(' *** Must be positive - try again ***')
          go to 210
        end if
230     write(6,240,advance='no')
240     format(' Enter threshold concentration defining inactive cell: ')
	i=key_read(thresh_inactive)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
          go to 210
	else if(i.eq.-1) then
	  go to 230
	else if(i.ne.0) then
	  write(6,195)
	  go to 230
	end if
        if(thresh_inactive.le.thresh_interface)then
          write(6,*)
          write(6,241)
241       format(' *** Must be greater than interface concentration. Try these again ***')
          write(6,*)
          go to 210
        end if

350     write(6,360)
360     format(' Use dummy value or cell midpoint elev when interface above/below top/bottom')
        write(6,370,advance='no')
370     format('   active cell centre? [d/m]: ')
        read(5,'(a)') adm
        if(adm.eq.' ') go to 350
        call casetrans(adm,'lo')
        if(adm.eq.'e') then
          write(6,*)
          go to 230
        end if
        if((adm.ne.'d').and.(adm.ne.'m')) go to 350
373     continue
        if(adm.eq.'d')then
371       write(6,372,advance='no')
372       format(' Enter dummy value for above top cell centre: ')
          i=key_read(dummy1)
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 350
          else if(i.eq.-1) then
            go to 371
          else if(i.ne.0) then
            write(6,195)
            go to 371
          end if
374       write(6,375,advance='no')
375       format(' Enter dummy value for below bottom cell centre: ')
          i=key_read(dummy2)
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 371
          else if(i.eq.-1) then
            go to 374
          else if(i.ne.0) then
            write(6,195)
            go to 374
          end if
        end if

        allocate(conc(ncol,nrow,ilay1:ilay2),filled(ilay1:ilay2),elev(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200
        filled=0
        iarray=0.0
        call addquote(ucnfile,afile)
        write(6,379) trim(afile)
379     format(/,' - reading file ',a,'...')
        do
          iarray=iarray+1
          mrow=0
          mcol=0
	  read(ucnunit,err=9100,end=300) ntrans,kstp,kper,   &
          rtemp,text,mcol,mrow,ilay
          if((mcol.le.0).or.(mrow.le.0).or.(ilay.le.0)) go to 9100
          if((mrow.ne.nrow).or.(mcol.ne.ncol)) then
            call num2char(iarray,anum)
            write(amessage,250) trim(anum),trim(afile),trim(gridspec%specfile)
250         format(' Number of rows and columns read from header to array ',&
	    'number ',a,' from MT3D concentration file ',a,' does not agree with ',&
	    'grid specifications as read from grid specification file ',a,   &
            '. Alternatively, unformatted file protocol might be a problem - ', &
            'try using an alternative version of this program.')
	    go to 9890
	  end if
          read(ucnunit,err=9130,end=9150) ((elev(icol,irow),icol=1,ncol),irow=1,nrow)
          if(equals(totim,rtemp))then
            if((ilay.ge.ilay1).and.(ilay.le.ilay2))then
              do irow=1,nrow
                do icol=1,ncol
                  conc(icol,irow,ilay)=elev(icol,irow)
                end do
              end do
              filled(ilay)=1
            end if
          end if
        end do
300     continue
        do ilay=ilay1,ilay2
          if(filled(ilay).ne.0) go to 320
        end do
        write(amessage,310) trim(afile)
310     format(' No arrays were found in file ',a,' corresponding to the nominated ',  &
        'simulation time.')
        go to 9890
320     do ilay=ilay1,ilay2
          if(filled(ilay).eq.0)then
            call num2char(ilay,alay)
            write(amessage,330) trim(afile),trim(alay)
330         format(' No concentration array was found in file ',a,' for layer ',a,  &
            ' at the nominated simulation time.')
            go to 9890
          end if
        end do
        write(6,340) trim(afile)
340     format(' - file ',a,' read ok.')

! -- First we evaluate the elevation of the interface in each vertical column of the grid.

        elev=1.0e35
        do irow=1,nrow
          do icol=1,ncol
            icount=0
            do ilay=ilay1,ilay2
              if(abs(conc(icol,irow,ilay)).gt.thresh_inactive) then
                icount=0    ! restart the count
                cycle
              end if
              icount=icount+1
              if(conc(icol,irow,ilay).ge.thresh_interface)then
                if(icount.eq.1)then
                  if(adm.eq.'d')then
                    elev(icol,irow)=dummy1
                  else
                    elev(icol,irow)=0.5*(bottom(icol,irow,ilay-1)+bottom(icol,irow,ilay))
                  end if
                else
                  elev1=bottom(icol,irow,ilay-2)
                  elev2=bottom(icol,irow,ilay-1)
                  elev3=bottom(icol,irow,ilay)
                  mid1=0.5*(elev1+elev2)
                  mid2=0.5*(elev2+elev3)
                  conc1=conc(icol,irow,ilay-1)
                  conc2=conc(icol,irow,ilay)
                  rtemp=mid1+(thresh_interface-conc1)/(conc2-conc1)*(mid2-mid1)
                  elev(icol,irow)=rtemp
                end if
                go to 363
              else
                mid2=0.5*(bottom(icol,irow,ilay)+bottom(icol,irow,ilay-1))
              end if
            end do
            if(adm.eq.'d')then
              elev(icol,irow)=dummy2
            else
              elev(icol,irow)=mid2
            end if
363         continue
          end do
        end do

! -- The elevation array is recorded.

        write(6,*)
1350    continue
        aprompt=' Enter name for interface elevation real array output file: '
        call write_real_array(ifail,aprompt,elev,pm_header=headerspec, &
        rows=nrow,columns=ncol)
        if(escset.eq.1) then
          escset=0
          write(6,*)
          deallocate(conc,elev,filled)
          rewind(unit=ucnunit)
          if(adm.eq.'d')then
            go to 373
          else
            go to 350
          end if
        end if

! -- The row/column intersection file is opened.

1400    write(6,*)
1410    call open_output_file(ifail, &
	' Enter name for row/column intersection file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  go to 1350
	end if

! -- First we obtain coordinates of grid cell centres in grid coordinates.

        allocate(east(ncol,nrow),north(ncol,nrow),x(ncol,nrow),y(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200

	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)
        call rel_centre_coords(x,y,gridspec)

        do ilay=ilay1,ilay2
          call num2char(ilay,alay)
          write(outunit,1460) trim(alay)
1460      format(/,' Row intersections of interface for layer ',a,' ---->')
          do irow=1,nrow
            icross=0
            do icol=2,ncol
              conc1=conc(icol-1,irow,ilay)
              conc2=conc(icol,irow,ilay)
              if(abs(conc1).gt.thresh_inactive) cycle
              if(abs(conc2).gt.thresh_inactive) cycle
              if(((conc1.lt.thresh_interface).and.(conc2.ge.thresh_interface)).or.   &
                 ((conc1.ge.thresh_interface).and.(conc2.lt.thresh_interface)))then
                icross=icross+1
                if(icross.gt.MAXCROSS)then
                  write(amessage,1470)
1470              format(' Increase MAXCROSS and re-compile program.')
                  go to 9890
                end if
                frac=(thresh_interface-conc1)/(conc2-conc1)
                dist(icross)=frac*(x(icol,irow)-x(icol-1,irow))+x(icol-1,irow)
                eastdist(icross)=frac*(east(icol,irow)-east(icol-1,irow))+east(icol-1,irow)
                northdist(icross)=frac*(north(icol,irow)-north(icol-1,irow))+north(icol-1,irow)
                mid2=0.5*(bottom(icol,irow,ilay)+bottom(icol,irow,ilay-1))
                mid1=0.5*(bottom(icol-1,irow,ilay)+bottom(icol-1,irow,ilay-1))
                z(icross)=frac*(mid2-mid1)+mid1
                if(conc2.gt.conc1)then
                  itype(icross)=1
                else
                  itype(icross)=-1
                end if
              end if
            end do
            if(icross.eq.0)then
              write(outunit,1480) irow,icross
1480          format(i7,i7,1000(2x,i3,2x,e15.8,2x,f15.3,2x,f15.3,2x,e15.8))
            else
              write(outunit,1480) irow,icross,(itype(i),dist(i),  &
              eastdist(i)+gridspec%east_corner,northdist(i)+gridspec%north_corner,z(i),i=1,icross)
            end if
          end do
          call num2char(ilay,alay)
          write(outunit,1461) trim(alay)
1461      format(/' Column intersections of interface for layer ',a,' ---->')
          do icol=1,ncol
            icross=0
            do irow=2,nrow
              conc1=conc(icol,irow-1,ilay)
              conc2=conc(icol,irow,ilay)
              if(abs(conc1).gt.thresh_inactive) cycle
              if(abs(conc2).gt.thresh_inactive) cycle
              if(((conc1.lt.thresh_interface).and.(conc2.ge.thresh_interface)).or.   &
                 ((conc1.ge.thresh_interface).and.(conc2.lt.thresh_interface)))then
                icross=icross+1
                if(icross.gt.MAXCROSS)then
                  write(amessage,1470)
                end if
                frac=(thresh_interface-conc1)/(conc2-conc1)
                dist(icross)=frac*(y(icol,irow)-y(icol,irow-1))+y(icol,irow-1)
                dist(icross)=-dist(icross)
                eastdist(icross)=frac*(east(icol,irow)-east(icol,irow-1))+east(icol,irow-1)
                northdist(icross)=frac*(north(icol,irow)-north(icol,irow-1))+north(icol,irow-1)
                mid2=0.5*(bottom(icol,irow,ilay)+bottom(icol,irow,ilay-1))
                mid1=0.5*(bottom(icol,irow-1,ilay)+bottom(icol,irow-1,ilay-1))
                z(icross)=frac*(mid2-mid1)+mid1
                if(conc2.gt.conc1)then
                  itype(icross)=1
                else
                  itype(icross)=-1
                end if
              end if
            end do
            if(icross.eq.0)then
              write(outunit,1480) icol,icross
            else
              write(outunit,1480) icol,icross,(itype(i),dist(i),  &
              eastdist(i)+gridspec%east_corner,northdist(i)+gridspec%north_corner,z(i),i=1,icross)
            end if
          end do
        end do
        call addquote(outfile,afile)
        write(6,340) trim(afile)

	go to 9900

9000    call addquote(infile,afile)
        write(amessage,9010) trim(afile)
9010    format(' Error encountered in reading real array from file ',a,'.')
        go to 9890
9050    call addquote(infile,afile)
        write(amessage,9060) trim(afile)
9060    format(' Premature end to file encountered when reading array from file ',a,'.')
        go to 9890
9070    call addquote(infile,afile)
        write(amessage,9080) trim(afile)
9080    format(' File ',a,' appears to have a number-of-columns, number-of-rows header. ', &
        'This is not allowed for this sequence of arrays.')
        go to 9890
9100	call num2char(iarray,anum)
        call addquote(ucnfile,afile)
	write(amessage,9110) trim(anum),trim(afile)
9110	format(' Error reading header to array number ',a,' found in file ',a,'.')
	go to 9890
9130    call num2char(iarray,anum)
        call addquote(ucnfile,afile)
        write(amessage,9140) trim(anum),trim(afile)
9140    format(' Error encountered while reading array number ',a,' from file ',a,'.')
        go to 9890
9150    call num2char(iarray,anum)
        call addquote(ucnfile,afile)
        write(amessage,9160) trim(afile),trim(anum)
9160    format(' Premature end to file ',a,' encountered while reading array number ',a,'.')
        go to 9890
9200	write(amessage,9210)
9210	format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)

        if(allocated(filled)) deallocate(filled,stat=ierr)
        if(allocated(bottom)) deallocate(bottom,stat=ierr)
        if(allocated(conc))   deallocate(conc,stat=ierr)
        if(allocated(elev))   deallocate(elev,stat=ierr)
        if(allocated(east))   deallocate(east,stat=ierr)
        if(allocated(north))  deallocate(north,stat=ierr)
        if(allocated(x))      deallocate(x,stat=ierr)
        if(allocated(y))      deallocate(y,stat=ierr)

end program conc2elev

