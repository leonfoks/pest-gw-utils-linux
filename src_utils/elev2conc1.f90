program elev2conc1

! -- Program ELEV2CONC1 is identical to ELEV2CONC except that it also writes formatted "zero flow head" files.

       use defn
       use inter

       implicit none

       integer  :: ncol,nrow,laylo,layhi,ilay,irow,icol,ntrans,kstp,kper
       integer  :: ierr,ichoice,ifail,ibeg,iend,nlay
       integer  :: idate,iheader,nbb,iunit,mcol,mrow,i
       integer  :: outunit,uoutunit
       real     :: rtemp,cellelev,diff,totim,wid,pi,a
       real     :: rof,cs,slope,rtemp1
       character*1   :: asw
       character*10  :: atemp
       character*10  :: alay
       character*16  :: text
       character*20  :: anum
       character*100 :: aprompt
       character*200 :: specfile,actbase,infile,elevbase,saltfile,  &
                        afile,outbase,uoutfile,outfile,concbase,bfile,houtbase

       integer, allocatable :: active(:,:,:)
       real, allocatable    :: top(:,:,:),interface(:,:),width(:,:),fresh(:,:),salt(:,:),conc(:,:,:)

       type (modelgrid) gridspec


       pi=3.1415927

       write(amessage,1)
1      format(' Program ELEV2CONC1 writes MT3D and/or SEAWAT initial concentration ', &
       'arrays based on the spatial disposition of the saltwater/freshwater interface as ',  &
       'supplied by the user. It also writes "zero flow initial head" arrays.')
       call write_message(leadspace='yes',endspace='yes')

       include 'unformat1.inc'
	
       call read_settings(ifail,idate,iheader)
       if(ifail.eq.1) then
         write(amessage,7)
7        format(' A settings file (settings.fig) was not found in the ', &
         'current directory.')
         call write_message
         go to 9900
       else if(ifail.eq.2) then
         write(amessage,8)
8        format(' Error encountered while reading settings file settings.fig')
         call write_message
         go to 9900
       endif
       if((iheader.ne.0).or.(headerspec.eq.' ')) then
         write(amessage,6)
6        format(' Cannot read array header specification from settings file ', &
         'settings.fig')
         call write_message
         go to 9900
       end if

       call readfig(gridspec%specfile)
10     call spec_open(ifail,gridspec)
       if(ifail.ne.0) go to 9900
       if(escset.eq.1) go to 9900
       call read_spec_dim(ifail,gridspec)
       if(ifail.ne.0) go to 9900
       call close_spec_file(gridspec,ok='yes')

       ncol=gridspec%ncol
       nrow=gridspec%nrow

       write(6,*)
140    write(6,142,advance='no')
142    format(' Enter number of layers in model: ')
       i=key_read(nlay)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 10
       else if(i.eq.-1) then
         go to 140
       else if(i.ne.0) then
         write(6,144)
144      format(' Illegal input  - try again.')
         go to 140
       end if
       if(nlay.le.0) then
         write(6,144)
         go to 140
       end if

150    write(6,155,advance='no')
155    format(' Enter smaller layer number for present analysis: ')
       i=key_read(laylo)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 140
       else if(i.eq.-1) then
         go to 150
       else if(i.ne.0) then
         write(6,144)
         go to 150
       end if
       if((laylo.le.0).or.(laylo.gt.nlay)) then
         write(6,144)
         go to 150
       end if
160    write(6,165,advance='no')
165    format(' Enter larger layer number for present analysis: ')
       i=key_read(layhi)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 150
       else if(i.eq.-1) then
         go to 160
       else if(i.ne.0) then
         write(6,144)
         go to 160
       end if
       if((layhi.le.0).or.(layhi.gt.nlay)) then
         write(6,144)
         go to 160
       end if
       if(layhi.lt.laylo)then
         write(6,167)
167      format(' Must not be less than smaller layer - try again.')
         go to 160
       end if

! -- Arrays are allocated.

       allocate(active(ncol,nrow,laylo:layhi),                                   &
                top(ncol,nrow,laylo:layhi+1),                                    &
                interface(ncol,nrow),                                            &
                fresh(ncol,nrow),                                                &
                salt(ncol,nrow),                                                 &
                width(ncol,nrow),                                                &
                conc(ncol,nrow,nlay),stat=ierr)
       if(ierr.ne.0) go to 9200

! -- Layer activity arrays are read.

       write(6,*)
170    write(6,180,advance='no')
180    format(' Enter filename base for activity arrays: ')
       read(5,'(a)') afile
       if(afile.eq.' ') go to 170
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         deallocate(active,top,interface,fresh,salt,width,conc)
         write(6,*)
         go to 160
       end if
       nbb=len_trim(afile)
       call getfile(ifail,afile,actbase,1,nbb)
       if(ifail.ne.0) go to 170
       do ilay=laylo,layhi
         call num2char(ilay,alay)
         infile=trim(actbase)//trim(alay)//'.inf'
         call addquote(infile,afile)
         iunit=nextunit()
         open(unit=iunit,file=infile,status='old',iostat=ierr)
         if(ierr.ne.0)then
           write(amessage,190) trim(afile)
190        format(' Cannot open file ',a,'.')
           go to 9890
         end if
         if(headerspec.eq.'yes')then
           read(iunit,*,err=9000,end=9000) mcol,mrow
           if((mrow.ne.nrow).or.(mcol.ne.ncol))then
             call addquote(gridspec%specfile,bfile)
             write(amessage,191) trim(afile),trim(bfile)
191          format(' Number of columns or number of rows in header to array file ',a,   &
             ' is not the same as that in finite-difference grid as read from from ',    &
             'grid specification file ',a,'.')
             go to 9890
           end if
         end if
         do irow=1,nrow
           read(iunit,*,iostat=ierr) (active(icol,irow,ilay),icol=1,ncol)
           if(ierr.ne.0)then
             write(amessage,200) trim(afile)
200          format(' Error encountered while reading array from file ',a,'.')
             go to 9890
           end if
         end do
         close(unit=iunit)
         write(6,210) trim(afile)
210      format('  - file ',a,' read ok.')
       end do

! -- Array top elevations are now read.

       write(6,*)
230    write(6,240,advance='no')
240    format(' Enter filename base of layer bottom elevation arrays: ')
       read(5,'(a)') afile
       if(afile.eq.' ') go to 230
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         write(6,*)
         go to 170
       end if
       nbb=len_trim(afile)
       call getfile(ifail,afile,elevbase,1,nbb)
       if(ifail.ne.0) go to 230
       do ilay=laylo,layhi+1
         call num2char(ilay-1,alay)
         infile=trim(elevbase)//trim(alay)//'.ref'
         iunit=nextunit()
         call addquote(infile,afile)
         open(unit=iunit,file=infile,status='old',iostat=ierr)
         if(ierr.ne.0)then
           write(amessage,190) trim(afile)
           go to 9890
         end if
         if(headerspec.eq.'yes')then
           read(iunit,*,err=9000,end=9000) mcol,mrow
           if((mrow.ne.nrow).or.(mcol.ne.ncol))then
             call addquote(gridspec%specfile,bfile)
             write(amessage,191) trim(afile),trim(bfile)
             go to 9890
           end if
         end if
         do irow=1,nrow
           read(iunit,*,err=9050,end=9050)(top(icol,irow,ilay),icol=1,ncol)
         end do
         close(unit=iunit)
         write(6,210) trim(afile)
       end do

! -- The interface elevation array file is now read.

       write(6,*)
300    continue
       aprompt=' Enter name of interface elevation array file: '
       call read_real_array(ifail,aprompt,interface,pm_header=headerspec, &
       rows=nrow,columns=ncol)
       if(ifail.ne.0) go to 9900
       if(escset.eq.1) then
         escset=0
         write(6,*)
         go to 230
       end if

! -- The nature of the interface is now acquired.

       write(6,*)
350    write(6,360)
360    format(' Enter nature of concentration variation across interface:-')
       write(6,380)
380    format('      if linear      - enter 1')
       write(6,390)
390    format('      if sigmoidal   - enter 2')
       write(6,391)
391    format('      if exponential - enter 3')
400    write(6,410,advance='no')
410    format(' Enter your choice: ')
       i=key_read(ichoice)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 300
       else if(i.eq.-1) then
         go to 400
       else if(i.ne.0) then
         write(6,144)
         go to 400
       else if((ichoice.ne.1).and.(ichoice.ne.2).and.(ichoice.ne.3)) then
         write(6,144)
         go to 400
       end if

! -- Concentration data is now acquired.

       write(6,*)
420    write(6,430,advance='no')
430    format(' Enter interface width value or formatted real array file: ')
       read(5,'(a)') afile
       if(afile.eq.' ') go to 420
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         write(6,*)
         go to 350
       end if
       anum=afile(1:20)
       call char2num(ifail,anum,rtemp)
       if(ifail.eq.0)then
         if(rtemp.lt.0.0)then
           write(6,431)
431        format(' Must not be negative - try again.')
           go to 420
         end if
         width=rtemp    ! an array
       else
         ibeg=1
         iend=len_trim(afile)
         call getfile(ifail,afile,infile,ibeg,iend)
         if(ifail.ne.0) go to 420
         call addquote(infile,afile)
         iunit=nextunit()
         open(unit=iunit,file=infile,status='old',iostat=ierr)
         if(ierr.ne.0)then
           write(6,435) trim(afile)
435        format(' Cannot open file ',a,' - try again.')
           go to 420
         end if
         if(headerspec.eq.'yes')then
           read(iunit,*,err=9000,end=9000) mcol,mrow
           if((mrow.ne.nrow).or.(mcol.ne.ncol))then
             call addquote(gridspec%specfile,bfile)
             write(amessage,191) trim(afile),trim(bfile)
             go to 9890
           end if
         end if
         do irow=1,nrow
           read(iunit,*,iostat=ierr) (width(icol,irow),icol=1,ncol)
           if(ierr.ne.0)then
             write(amessage,440) trim(afile)
440          format(' Error reading array from file ',a,'.')
             go to 9890
           end if
         end do
         close(unit=iunit)
         write(6,450) trim(afile)
450      format('  - file ',a,' read ok.')
         do irow=1,nrow
           do icol=1,ncol
             if(width(icol,irow).lt.0.0) width(icol,irow)=0.0
           end do
         end do
       end if

470    write(6,480,advance='no')
480    format(' Enter fresh water concentration value or formatted real array file: ')
       read(5,'(a)') afile
       if(afile.eq.' ') go to 470
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         write(6,*)
         go to 420
       end if
       anum=afile(1:20)
       call char2num(ifail,anum,rtemp)
       if(ifail.eq.0)then
         fresh=rtemp    ! an array
       else
         ibeg=1
         iend=len_trim(afile)
         call getfile(ifail,afile,infile,ibeg,iend)
         call addquote(infile,afile)
         if(ifail.ne.0) go to 470
         iunit=nextunit()
         open(unit=iunit,file=infile,status='old',iostat=ierr)
         if(ierr.ne.0)then
           write(6,435) trim(afile)
           go to 470
         end if
         if(headerspec.eq.'yes')then
           read(iunit,*,err=9000,end=9000) mcol,mrow
           if((mrow.ne.nrow).or.(mcol.ne.ncol))then
             call addquote(gridspec%specfile,bfile)
             write(amessage,191) trim(afile),trim(bfile)
             go to 9890
           end if
         end if
         do irow=1,nrow
           read(iunit,*,iostat=ierr) (fresh(icol,irow),icol=1,ncol)
           if(ierr.ne.0)then
             write(amessage,440) trim(afile)
             go to 9890
           end if
         end do
         close(unit=iunit)
         write(6,450) trim(infile)
       end if

510    write(6,520,advance='no')
520    format(' Enter sea water concentration value or formatted real array file: ')
       read(5,'(a)') afile
       if(afile.eq.' ') go to 510
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e')then
         write(6,*)
         go to 470
       end if
       anum=afile(1:20)
       call char2num(ifail,anum,rtemp)
       if(ifail.eq.0)then
         salt=rtemp ! an array
       else
         ibeg=1
         iend=len_trim(afile)
         call getfile(ifail,afile,infile,ibeg,iend)
         if(ifail.ne.0) go to 510
         call addquote(infile,afile)
         open(unit=iunit,file=infile,status='old',iostat=ierr)
         if(ierr.ne.0)then
           write(6,435) trim(afile)
           go to 510
         end if
         if(headerspec.eq.'yes')then
           read(iunit,*,err=9000,end=9000) mcol,mrow
           if((mrow.ne.nrow).or.(mcol.ne.ncol))then
             call addquote(gridspec%specfile,bfile)
             write(amessage,191) trim(afile),trim(bfile)
             go to 9890
           end if
         end if
         do irow=1,nrow
           read(iunit,*,iostat=ierr) (salt(icol,irow),icol=1,ncol)
           if(ierr.ne.0)then
             write(amessage,440) trim(afile)
             go to 9890
           end if
         end do
         close(unit=iunit)
         write(6,450) trim(afile)
       end if

! -- The handling of cells with activity of -1 is now ascertained.


530    write(6,540,advance='no')
540    format(' Assign sea water concentration to cells with negative activity?  [y/n]: ')
       read(5,'(a)')  asw
       call casetrans(asw,'lo')
       if(asw.eq.' ') go to 530
       asw=adjustl(asw)
       if(asw.eq.'e')then
         write(6,*)
         go to 510
       else if(asw.eq.'y')then
         continue
       else if(asw.eq.'n')then
         continue
       else
         go to 530
       end if

       write(6,*)
541    write(6,542,advance='no')
542    format(' Enter density of water of zero concentration (<Enter> if 1000): ')
       i=key_read(rof)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 530
       else if(i.eq.-1) then
         rof=1000.0
       else if(i.ne.0) then
         write(6,144)
         go to 541
       end if
       if(rof.le.0.0) then
         write(6,144)
         go to 541
       end if
543    write(6,544,advance='no')
544    format(' Enter concentration of water at zero head (<Enter> if 35): ')
       i=key_read(cs)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 541
       else if(i.eq.-1) then
         cs=35
       else if(i.ne.0) then
         write(6,144)
         go to 543
       end if
       if(cs.le.0.0) then
         write(6,144)
         go to 543
       end if
545    write(6,546,advance='no')
546    format(' Enter delta-density/delta-concentration slope (<Enter> if 0.7143): ')
       i=key_read(slope)
       if(escset.ne.0)then
         escset=0
         write(6,*)
         go to 543
       else if(i.eq.-1) then
         slope=0.7143
       else if(i.ne.0) then
         write(6,144)
         go to 545
       end if
       if(slope.lt.0.0) then
         write(6,144)
         go to 545
       end if


! -- Names of output files are now acquired.

       write(6,*)
550    write(6,560,advance='no')
560    format(' Enter filename base for formatted concentration output files: ')
       read(5,'(a)') afile
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         write(6,*)
         go to 545
       end if
       nbb=len_trim(afile)
       call getfile(ifail,afile,outbase,1,nbb)
       if(ifail.ne.0) go to 550
570    aprompt = ' Enter name for unformatted concentration output file: '
       call open_output_file(ifail,aprompt,uoutfile,uoutunit,file_format='unformatted')
       if(ifail.ne.0) go to 9900
       if(escset.eq.1)then
         escset=0
         write(6,*)
         go to 550
       end if

! -- The name for zero flow head output files is obtained.

       write(6,*)
700    write(6,710,advance='no')
710    format(' Enter filename base for formatted "zero flow head" output files: ')
       read(5,'(a)') afile
       atemp=afile(1:2)
       call casetrans(atemp,'lo')
       if(atemp(1:2).eq.'e ')then
         write(6,*)
         close(unit=uoutunit)
         go to 570
       end if
       nbb=len_trim(afile)
       call getfile(ifail,afile,houtbase,1,nbb)
       if(ifail.ne.0) go to 700

       concbase=' '
       if((laylo.gt.1).or.(layhi.lt.nlay))then
         write(6,*)
582      write(6,581,advance='no')
581      format(' Enter existing conc array filename base of unanalysed layers: ')
         read(5,'(a)') afile
         if(afile.eq.' ') go to 582
         atemp=afile(1:2)
         call casetrans(atemp,'lo')
         if(atemp(1:2).eq.'e ')then
           write(6,*)
           go to 700
         end if
         nbb=len_trim(afile)
         call getfile(ifail,afile,concbase,1,nbb)
         if(ifail.ne.0) go to 582
         do ilay=1,nlay
           if((ilay.ge.laylo).and.(ilay.le.layhi))cycle
           call num2char(ilay,alay)
           infile=trim(concbase)//trim(alay)//'.ref'
           call addquote(infile,afile)
           iunit=nextunit()
           open(unit=iunit,file=infile,status='old',iostat=ierr)
           if(ierr.ne.0)then
             write(amessage,588) trim(afile)
588          format(' Cannot open file ',a,'.')
             go to 9890
           end if
           if(headerspec.eq.'yes')then
             read(iunit,*,err=9000,end=9000) mcol,mrow
             if((mrow.ne.nrow).or.(mcol.ne.ncol))then
               call addquote(infile,afile)
               call addquote(gridspec%specfile,bfile)
               write(amessage,191) trim(afile),trim(bfile)
               go to 9890
             end if
           end if
           do irow=1,nrow
             read(iunit,*,err=9050,end=9050) (conc(icol,irow,ilay),icol=1,ncol)
           end do
           close(unit=iunit)
           write(6,584) trim(afile)
584        format('  - file ',a,' read ok.')
         end do
       end if

! -- We now determine the concentration in each cell in the model.

       write(6,590)
590    format(/,' Working...')
       do ilay=laylo,layhi
         do irow=1,nrow
           do icol=1,ncol
             if(active(icol,irow,ilay).eq.0) then
               conc(icol,irow,ilay)=0.0
               cycle
             else if(active(icol,irow,ilay).lt.0)then
               if(asw.eq.'y')then
                 conc(icol,irow,ilay)=salt(icol,irow)
                 cycle
               end if
             end if
             cellelev=0.5*(top(icol,irow,ilay)+top(icol,irow,ilay+1))
             diff=cellelev-interface(icol,irow)
             wid=width(icol,irow)
             if(ichoice.eq.1)then
               if(diff.gt.0.5*wid)then
                 conc(icol,irow,ilay)=fresh(icol,irow)
               else if(diff.lt.-0.5*wid)then
                 conc(icol,irow,ilay)=salt(icol,irow)
               else
                 conc(icol,irow,ilay)=salt(icol,irow)+     &
                 (fresh(icol,irow)-salt(icol,irow))*       &
                 (diff+wid*0.5)/wid
               end if
             else if(ichoice.eq.2)then
               if(diff.gt.0.5*wid)then
                 conc(icol,irow,ilay)=fresh(icol,irow)
               else if(diff.lt.-0.5*wid)then
                 conc(icol,irow,ilay)=salt(icol,irow)
               else
                 conc(icol,irow,ilay)=(sin(diff/wid*pi)+1.0)*0.5*      &
                 (fresh(icol,irow)-salt(icol,irow))+salt(icol,irow)
               end if
             else if(ichoice.eq.3)then
               a=-log(0.2)/(wid*0.5)
               rtemp=0.5*(1.0-exp(-a*abs(diff)))
               if(diff.gt.0)then
                 rtemp=0.5+rtemp
               else
                 rtemp=0.5-rtemp
               end if
               conc(icol,irow,ilay)=salt(icol,irow)+                    &
               (fresh(icol,irow)-salt(icol,irow))*rtemp
             end if
           end do
         end do
       end do

! -- Formatted output files are now written.

       write(6,600)
600    format(/,' Writing output files...')
       do ilay=laylo,layhi
         call num2char(ilay,alay)
         outfile=trim(outbase)//trim(alay)//'.ref'
         outunit=nextunit()
         open(unit=outunit,file=outfile)
         if(headerspec.eq.'yes')then
           write(outunit,601) ncol,nrow
601        format(2i10)
         end if
         do irow=1,nrow
           write(outunit,603) (conc(icol,irow,ilay),icol=1,ncol)
603        format(8(1x,1pg14.7))
         end do
         close(unit=outunit)
         write(6,620) trim(outfile)
620      format( '  - file ',a,' written ok.')
      end do

! -- The unformatted output file is written.

      ntrans=1
      kstp=1
      kper=1
      totim=1.0
      text='concentration'
      do ilay=1,nlay
        write(uoutunit) ntrans,kstp,kper,totim,text,ncol,nrow,ilay
        write(uoutunit) ((conc(icol,irow,ilay),icol=1,ncol),irow=1,nrow)
      end do
      close(unit=uoutunit)
      call addquote(uoutfile,afile)
      write(6,*)
      write(6,620) trim(afile)

! -- The concentration array is converted to a zero flow head array.

      do ilay=laylo,layhi
        do irow=1,nrow
          do icol=1,ncol
            if(active(icol,irow,ilay).ne.0)then
              rtemp=(top(icol,irow,ilay)+top(icol,irow,ilay+1))*0.5
              rtemp1=conc(icol,irow,ilay)
              conc(icol,irow,ilay)=rtemp*slope*(rtemp1-cs)/(slope*rtemp1+rof)
            else
              conc(icol,irow,ilay)=1.0e35
            end if
          end do
        end do
      end do

! -- Formatted head output files are written.

       do ilay=laylo,layhi
         call num2char(ilay,alay)
         outfile=trim(houtbase)//trim(alay)//'.ref'
         outunit=nextunit()
         open(unit=outunit,file=outfile)
         if(headerspec.eq.'yes')then
           write(outunit,601) ncol,nrow
         end if
         do irow=1,nrow
           write(outunit,603) (conc(icol,irow,ilay),icol=1,ncol)
         end do
         close(unit=outunit)
         write(6,620) trim(outfile)
      end do

      go to 9900

9000  write(amessage,9010) trim(afile)
9010  format(' Error encountered in reading array header from file ',a,'.')
      go to 9890
9050  write(amessage,9060) trim(afile)
9060  format(' Error encountered in reading array from file ',a,'.')
      go to 9890

9200  write(amessage,9210)
9210  format(' Cannot allocate sufficient memory to continue execution.')
      go to 9890

9890  continue
      call write_message(leadspace='yes',endspace='yes')

9900  continue
      deallocate(active,top,interface,width,fresh,salt,conc,stat=ierr)
      call close_files
      call free_grid_mem(gridspec)

end program elev2conc1

