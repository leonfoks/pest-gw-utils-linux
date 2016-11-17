program fem2smp

! -- Program FEM2SMP read MICROFEM "fth" and "ftq" files, creating a single bore
!    sample file from the contents of these files.

        use defn
        use inter
        implicit none

        integer              :: ifail,idate,iline,i,j,ierr,iobs,nbb,  &
                                outunit,iheader,ftqunit,fthunit,inunit,  &
                                fthnlay,fthnobs,ftqnlay,ftqnobs,maxid,maxtim, &
                                maxobs,fthline,ftqline,fthlinekp,ftqlinekp, &
                                dds,mms,yys,hhhs,mmms,ssss,nlay,nobs, &
                                fthid,ftqid,oldnobs,ifile,fthnobstot, &
                                ilay,jtime,itime,oldiobs,icount,nt,ftqnobstot
        integer, parameter   :: MAXLAY=80
        integer              :: fthnlays(MAXLAY),ftqnlays(MAXLAY),nlays(MAXLAY)
        integer, allocatable, dimension(:)  :: dd(:),mm(:),yy(:),hh(:),nn(:),ss(:)
        real                                :: day_convert,rtime
        real, allocatable, dimension(:)     :: rval                        

        character (len=5)    :: anum,aline
        character (len=10)   :: aid
        character (len=15)   :: atemp,adate,atime
        character (len=120)  :: outfile,fthfile,ftqfile,bfile,infile
        character (len=25), allocatable :: bid(:)




        write(amessage,5)
5       format(' Program FEM2SMP builds a bore sample file from ',  &
        'the contents of MicroFEM "fth" and "ftq" files.')
        call write_message(leadspace='yes',endspace='yes')

! -- Variables are initialised.

        fthnlay=0
        fthnobs=0
        ftqnlay=0
        ftqnobs=0
        maxid=0
        maxtim=0
        maxobs=0
        day_convert=1.0
        
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

10      write(6,15,advance='no')
15      format(' Enter name of "fth" file (press <Enter> if none): ')
        read(5,'(a)') fthfile
        if(fthfile.eq.' ')then
          fthunit=0
        else
          fthfile=adjustl(fthfile)
          if(index(eschar,fthfile(1:2)).ne.0) go to 9900
          nbb=len_trim(fthfile)
          call getfile(ifail,fthfile,bfile,1,nbb)
          if(ifail.ne.0) go to 10
          fthfile=bfile
          fthunit=nextunit()
          open(unit=fthunit,file=fthfile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(fthfile)
20          format(' Cannot open file ',a,' - try again.')
            call write_message()
            go to 10
          end if
        end if
22      write(6,23,advance='no')
23      format(' Enter name of "ftq" file (press <Enter> if none): ')
        read(5,'(a)') ftqfile
        if(ftqfile.eq.' ')then
          ftqunit=0
        else
          ftqfile=adjustl(ftqfile)
          if(index(eschar,ftqfile(1:2)).ne.0) then
            close(unit=fthunit)
            write(6,*)
            go to 10
          end if
          nbb=len_trim(ftqfile)
          call getfile(ifail,ftqfile,bfile,1,nbb)
          if(ifail.ne.0) go to 22
          ftqfile=bfile
          ftqunit=nextunit()
          open(unit=ftqunit,file=ftqfile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(ftqfile)
            call write_message()
            go to 22
          end if
        end if
        if((fthunit.eq.0).and.(ftqunit.eq.0))then
          write(amessage,25)
25        format(' The name of a "fth" file and/or a "ftq" file must ', &
          'be supplied.')
          go to 9800
        end if

300     write(6,*)
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
          close(unit=ftqunit)
          go to 22
        end if
        call char2date(ifail,adate,dds,mms,yys)
        if(ifail.ne.0)then
          write(6,340)
340       format(' Illegal date  - try again.')
          go to 310
        end if

360     write(6,370,advance='no')
370     format(' Enter simulation starting time [hh:mm:ss]: ')
        read(5,'(a)') atime
        if(atime.eq.' ') go to 360
        atime=adjustl(atime)
        if(index(eschar,atime(1:2)).ne.0) go to 300
        call char2time(ifail,atime,hhhs,mmms,ssss)
        if(ifail.ne.0)then
          write(6,380)
380       format(' Illegal time  - try again.')
          go to 360
        end if


140     write(6,*)
150     call open_output_file(ifail,   &
        ' Enter name for new bore sample file: ',outfile,outunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 360
        end if

! -- The number of lines in each of the fth and ftq files is now read.

        fthline=0
        if(fthunit.ne.0)then
          do
            read(fthunit,*,end=170)
            fthline=fthline+1
          end do
170       continue
          rewind(unit=fthunit)
        end if
        ftqline=0
        if(ftqunit.ne.0)then
          do
            read(ftqunit,*,end=175)
            ftqline=ftqline+1
          end do
175       continue
          rewind(unit=ftqunit)
        end if
        fthlinekp=fthline
        ftqlinekp=ftqline

! -- The first few lines of each of the files are read - until enough is read 
!    to allow dimensioning of arrays.

        fthline=0
        ftqline=0
        if(fthunit.ne.0)then
          iline=0
          infile=fthfile
          fthline=fthline+1
          iline=iline+1
          read(fthunit,*,err=9000,end=9050)
          iline=iline+1
          fthline=fthline+1
          read(fthunit,*,err=9000,end=9050)fthnlay,fthnobs
          if(fthnlay.gt.MAXLAY)then
            write(amessage,176)
176         format(' Too many model layers - increase MAXLAY and ',  &
            're-compile program.')
            go to 9800
          end if
          iline=iline+1
          fthline=fthline+1
          read(fthunit,*,err=9000,end=9050)(fthnlays(i),i=1,fthnlay)
        end if
        if(ftqunit.ne.0)then
          iline=0
          infile=ftqfile
          ftqline=ftqline+1
          iline=iline+1
          read(ftqunit,*,err=9000,end=9050)
          iline=iline+1
          ftqline=ftqline+1
          read(ftqunit,*,err=9000,end=9050)ftqnlay,ftqnobs
          if(ftqnlay.gt.MAXLAY)then
            write(amessage,176)
            go to 9800
          end if
          iline=iline+1
          ftqline=ftqline+1
          read(ftqunit,*,err=9000,end=9050)(ftqnlays(i),i=1,ftqnlay)
        end if
        fthnobstot=fthlinekp*fthnlay*fthnobs
        ftqnobstot=ftqlinekp*ftqnlay*ftqnobs
        maxobs=max(fthnobstot,ftqnobstot)
        if(maxobs.eq.0)then
          write(amessage,180) 
180       format(' No data is recorded in the "fth" and/or "ftq" files.')
          go to 9800
        end if
        allocate(rval(maxobs),stat=ierr)
        if(ierr.ne.0) go to 9200
        maxtim=max(fthlinekp,ftqlinekp)
        if(maxtim.eq.0)then
          write(amessage,180)
          go to 9800
        end if
        allocate(dd(maxtim),mm(maxtim),yy(maxtim),hh(maxtim),   &
        nn(maxtim),ss(maxtim),stat=ierr)
        if(ierr.ne.0) go to 9200
        fthid=fthnlay*fthnobs
        ftqid=ftqnlay*ftqnobs
        maxid=fthid+ftqid
        if(maxid.eq.0)then
          write(amessage,180)
          go to 9800
        end if
        allocate(bid(maxid),stat=ierr)
        if(ierr.ne.0) go to 9200
          
! -- Data is now read from each of the files and transferred to the bore
!    sample file.          

        oldnobs=0
        ifile=0
200     continue
        iline=3
        iobs=0
        ifile=ifile+1
        if(ifile.eq.1)then
          inunit=fthunit
          infile=fthfile
          nlay=fthnlay
          nobs=fthnobs
          do ilay=1,fthnlay
            nlays(ilay)=fthnlays(ilay)
          end do
        else if(ifile.eq.2)then
          inunit=ftqunit
          infile=ftqfile
          nlay=ftqnlay
          nobs=ftqnobs
          do ilay=1,ftqnlay
            nlays(ilay)=ftqnlays(ilay)
          end do          
        else
          go to 500
        end if
        if(inunit.eq.0) go to 200
        if((nlay.eq.0).or.(nobs.eq.0)) go to 200
        iline=iline+1
        read(inunit,*,err=9000,end=9050)
        read(inunit,*,err=9000,end=9050) atemp,(bid(i),i=oldnobs+1,oldnobs+nobs)
        if(nlay.gt.1)then
          do ilay=2,nlay
            call num2char(nlays(ilay),anum)
            anum=adjustl(anum)
            do i=1,nobs
              j=oldnobs+(ilay-1)*nobs+i
              bid(j)=trim(bid(oldnobs+i))//'_'//trim(anum)
            end do
          end do
          call num2char(nlays(1),anum)
          do i=1,nobs
            j=oldnobs+i
            bid(j)=trim(bid(j))//'_'//trim(anum)
          end do
        end if
        do i=1,2
          iline=iline+1
          read(inunit,*,err=9000,end=9050)
        end do
        oldiobs=0
        jtime=0
        do
          iline=iline+1
          jtime=jtime+1
          read(inunit,*,err=9000,end=305)   &
          (rtime,                            &
           (rval(oldiobs+(ilay-1)*nobs+iobs),iobs=1,nobs),ilay=1,nlay)
          oldiobs=oldiobs+nlay*nobs
          call elapsdate(rtime,day_convert,dds,mms,yys,hhhs,&
          mmms,ssss,dd(jtime),mm(jtime),yy(jtime),hh(jtime),nn(jtime),ss(jtime))
        end do
305     continue
        jtime=jtime-1
        write(amessage,315) trim(infile)
315     format(' - file ',a,' read ok.')
        call write_message()
        close(unit=inunit)
        icount=0
        do i=1,nobs*nlay
          nt=len_trim(bid(oldnobs+i))
          if(nt.gt.10) then
            bid(oldnobs+i)=bid(oldnobs+i)(nt-9:nt)
            icount=icount+1
          end if
          call casetrans(bid(oldnobs+i),'lo')
          do itime=1,jtime
            j=(itime-1)*nobs*nlay+i
            if(datespec.eq.1) then
              write(outunit,400,err=9400) trim(bid(oldnobs+i)),dd(itime),  &
              mm(itime),yy(itime),hh(itime),nn(itime),ss(itime),rval(j)
            else
              write(outunit,400,err=9400) trim(bid(oldnobs+i)),mm(itime),  &
              dd(itime),yy(itime),hh(itime),nn(itime),ss(itime),rval(j)
            end if
400         format(1x,a,t14,i2.2,'/',i2.2,'/',i4.4,t28,i2.2,':',i2.2,':',i2.2,&
            t38,1pg14.7)
          end do
        end do
        oldnobs=oldnobs+nlay*nobs        
        go to 200

500     continue

! -- The new bore sample file has now been written. Names are now checked.

        close(unit=fthunit)
        close(unit=ftqunit)
        

        imessage=0
        if((fthnlay.gt.0).or.(ftqnlay.gt.0))then
          write(amessage,510) trim(outfile)
510       format(' Warning: bore identifiers written to file ',a,  &
          ' include layer number suffixes. Measured data equivalents ', &
          'should include the same subscripts to allow matching of ', &
          'simulated and observed data.')
          call write_message(leadspace='yes',increment=1)
        end if
        if(icount.ne.0)then
          write(amessage,520) trim(outfile)
520       format(' Warning: at least some of the identifiers written ', &
          'to file ',a,' have had their names contracted from the right ', &
          'to adhere to the 10 character limit.')
          if(imessage.eq.0)then
            call write_message(leadspace='yes',increment=1)
          else
            call write_message(increment=1)
          end if
        end if
        if(oldnobs.gt.1)then
          do i=2,oldnobs
            aid=bid(i)
            do j=1,i-1
              if(bid(j).eq.aid)then
                write(amessage,530) trim(outfile)
530             format(' FATAL ERROR: bore identifiers recorded in file ',a, &
                ' are not unique. Hence this is not a valid bore sample ', &
                'file. You will need to review the names used in the ',  &
                'model to characterise observation points.')
                call write_message(leadspace='yes')
                go to 9900
              end if
            end do
          end do
        end if
        write(amessage,540) trim (outfile)
540     format(' - file ',a,' written ok.')
        call write_message(leadspace='yes')
        go to 9900
        

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim (infile)
9010    format(' Error reading line ',a,' of file ',a)
        go to 9800
9050    write(amessage,9060) trim(infile)
9060    format(' Unexpected end encountered to file ',a)
        go to 9800
9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9800
9400    write(amessage,9410) trim(outfile)
9410    format(' Cannot write to file ',a)
        go to 9800
        

9800    call write_message(leadspace='yes')
9900    call close_files
        deallocate(dd,mm,yy,hh,nn,ss,rval,bid,stat=ierr)
        write(6,*)

end program fem2smp
 