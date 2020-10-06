!     Last change:  JD   12 Apr 2003    2:49 am
program dar2smp

! -- Program DAR2SMP rewrites FEFLOW outputs as provided in a DAR file in bore sample file format.

        use defn
        use inter
        implicit none

        integer, parameter    :: MAXOBS=1000
        integer               :: ifail,nbb,iline,itemp,ierr,i,n,m
        integer               :: idate,iheader
        integer               :: day,mon,year,hour,min,sec
        integer               :: nobs,ntime,mout,nout,itime,ncol,icol,iobs,iout,iwarn
        integer               :: modunit,obsunit,outunit
        integer               :: dds,mms,yys,hhhs,mmms,ssss
        integer               :: borenum(MAXOBS)
        integer, allocatable  :: obsind(:),writetime(:)

        real                  :: day_convert
        real, allocatable     :: modval(:,:)
        real, allocatable     :: obsdays(:)

        double precision      :: time

        character (len=10)    :: aline
        character (len=10)    :: datatype
        character (len=20)    :: adate,atime
        character (len=25)    :: atemp
        character (len=200)   :: modfile,obsfile,outfile

        character (len=10)    :: boreid(MAXOBS)


	write(amessage,5)
5	format(' Program DAR2SMP extracts FEFLOW outputs from a ', &
        'DAR file and re-writes them in bore sample file format.')
	call write_message(leadspace='yes',endspace='yes')

        iwarn=0
        day_convert=1.0

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  go to 9890
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  go to 9890
	endif
	if((idate.ne.0).or.(datespec.eq.0)) then
	  write(amessage,9)
9	  format(' Cannot read date format from settings file ', &
	  'settings.fig')
	  go to 9890
	end if

! -- The DAR file is opened.

20	call open_input_file(ifail,' Enter name of FEFLOW DAR file: ',modfile,modunit)
	if(escset.ne.0) go to 9900
	if(ifail.ne.0) go to 9900

! The symbol indicating a data type for extraction is provided.

30      write(6,40,advance='no')
40      format(' Enter header for data type to extract from this file: ')
        read(5,'(a)') atemp
        if(atemp.eq.' ') go to 30
        if((atemp(1:2).eq.'e ').or.(atemp(1:2).eq.'E '))then
          close(unit=modunit)
          write(6,*)
          go to 20
        end if
        atemp=adjustl(atemp)
        nbb=len_trim(atemp)
        call getfile(ifail,atemp,datatype,1,nbb)
        if(ifail.ne.0) go to 30
        call casetrans(datatype,'hi')

! -- The observation number to boreid conversion file is opened.

        write(6,*)
50	call open_input_file(ifail,' Enter name of obs number to boreid conversion file: ',  &
        obsfile,obsunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          write(6,*)
          go to 30
        end if

! -- The starting date and time are read.

        write(6,*)
69      continue
        if(datespec.eq.1) then
          write(6,80,advance='no')
80        format(' Enter simulation starting date [dd/mm/yyyy]: ')
        else
          write(6,81,advance='no')
81        format(' Enter simulation starting date [mm/dd/yyyy]: ')
        end if
        read(5,'(a)') adate
        if(adate.eq.' ') go to 69
        adate=adjustl(adate)
        if(index(eschar,adate(1:2)).ne.0) then
          write(6,*)
          close(unit=obsunit)
          go to 50
        end if
        call char2date(ifail,adate,dds,mms,yys)
        if(ifail.ne.0)then
          write(6,90)
90        format(/,' Illegal date  - try again.',/)
          go to 69
        end if

100     write(6,110,advance='no')
110     format(' Enter simulation starting time [hh:mm:ss]: ')
        read(5,'(a)') atime
        if(atime.eq.' ') go to 100
        atime=adjustl(atime)
        if(index(eschar,atime(1:2)).ne.0) then
          write(6,*)
          go to 69
        end if
        call char2time(ifail,atime,hhhs,mmms,ssss)
        if(ifail.ne.0)then
          write(6,120)
120       format(/,' Illegal time  - try again.',/)
          go to 100
        end if

! -- The output file is now opened.

140	write(6,*)
150	call open_output_file(ifail,   &
	' Enter name for bore sample output file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 100
	end if

! -- The observation number to bore_id conversion file is now read.

        write(6,160) trim(obsfile)
160     format(/,' - reading obs number to boreid conversion file ',a,'...')
        nobs=0
        iline=0
        do
          iline=iline+1
          read(obsunit,'(a)',err=9000,end=250) cline
          if(cline.eq.' ') cycle
          cline=adjustl(cline)
          if(cline(1:1).eq.'#') cycle
          call linesplit(ifail,2)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,170) trim(aline),trim(obsfile)
170         format(' Two entries expected on line ',a,' of file ',a,'.')
            go to 9890
          end if
          itemp=char2int(ifail,1)
          if(ifail.ne.0) then
            call num2char(iline,aline)
            write(amessage,180) trim(aline),trim(obsfile)
180         format(' Error reading observation bore number from line ',a,' of file ',a,'.')
            go to 9890
          end if
          if(itemp.lt.0)then
            call num2char(iline,aline)
            write(amessage,190) trim(aline),trim(obsfile)
190         format(' Negative observation bore number not allowed at line ',a,' of file ',a,'.')
            go to 9890
          end if
          nobs=nobs+1
          if(nobs.gt.MAXOBS)then
            write(amessage,200) trim(obsfile)
200         format(' Too many observation wells cited in file ',a,'. Increase MAXOBS and re-compile program.')
            go to 9890
          end if
          borenum(nobs)=itemp
          atemp=cline(left_word(2):right_word(2))
          nbb=len_trim(atemp)
          if(nbb.gt.10)then
            call num2char(iline,aline)
            write(amessage,210) trim(aline),trim(obsfile)
210         format(' Bore identifier must not exceed 10 characters in length at line ',a,  &
            ' of file ',a,'.')
            go to 9890
          end if
          boreid(nobs)=atemp(1:10)
          call casetrans(boreid(nobs),'lo')
          if(nobs.gt.1)then
            do i=1,nobs-1
              if(boreid(i).eq.boreid(nobs))then
                call num2char(iline,aline)
                write(amessage,215) trim(aline),trim(obsfile)
215             format(' Bore identifier duplicated at line ',a,' of file ',a,'.')
                go to 9890
              end if
            end do
            do i=1,nobs-1
              if(borenum(i).eq.borenum(nobs))then
                call num2char(iline,aline)
                write(amessage,216) trim(aline),trim(obsfile)
216             format(' Observation bore number duplicated at line ',a,' of file ',a,'.')
                go to 9890
              end if
            end do
          end if
        end do
250     continue
        close(unit=obsunit)
        if(nobs.eq.0)then
          write(amessage,260) trim(obsfile)
260       format(' No observation wells are cited in file ',a,'.')
          go to 9890
        end if
        call num2char(nobs,atemp)
        write(6,270) trim(atemp),trim(obsfile)
270     format(' - data for ',a,' obs wells read from file ',a,'.')

! -- The DAR file is read a first time in order to ascertain the number of output times.
!    The maximum number of observation bores is also acquired.

        write(6,300) trim(modfile)
300     format(/,' - reading DAR file ',a,'...')
        ntime=0
        mout=0
        do
          read(modunit,'(a)',end=320) cline
          if(ntime.eq.1)mout=mout+1
          if(index(cline,'***').ne.0) ntime=ntime+1
        end do
320     ntime=ntime-1
        if(ntime.eq.0)then
          write(amessage,325) trim(modfile)
325       format(' There are no simulation times for which data is recorded in file ',a,'.')
          go to 9890
        end if
        allocate(obsdays(ntime),modval(ntime,nobs),writetime(ntime),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,330)
330       format(' Cannot allocate sufficient memory to continue execution.')
          go to 9890
        end if
        writetime=1                   ! an array
        allocate(obsind(mout),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,330)
          go to 9890
        end if

! -- The DAR file is rewound.

        rewind(unit=modunit)

! -- The first part of the DAR file is read to link observation numbers to observation names.

        iline=0
        do
          iline=iline+1
          read(modunit,'(a)',end=9100) cline
          if(index(cline,'***').ne.0) exit
        end do
        iline=iline+1
        read(modunit,'(a)',end=9100) cline
        if(index(cline,'LOCATION').eq.0) go to 9100
        iline=iline+1
        read(modunit,'(a)',end=9100) cline
        if(index(cline,'---').eq.0) go to 9100
        iline=iline+1
        read(modunit,'(a)',end=9100) cline
        call casetrans(cline,'lo')
        call linesplit(ifail,5)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,335) trim(aline),trim(modfile)
335       format(' Unexpected observation well header at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(cline(left_word(1):right_word(1)).ne.'obs')then
          call num2char(iline,aline)
          write(amessage,340) trim(aline),trim(modfile)
340       format(' "Obs" expected as first header text at line ',a,' of file ',a,'.')
          go to 9890
        end if
        iline=iline+1
        read(modunit,'(a)',end=9100) cline
        if(index(cline,'---').eq.0) go to 9100
        nout=0
        obsind=0            ! an array
        do
          iline=iline+1
          read(modunit,'(a)',end=9100) cline
          if(index(cline,'---').ne.0) go to 380
          call linesplit(ifail,1)
          if(ifail.ne.0) cycle
          itemp=char2int(ifail,1)
          if(ifail.ne.0) go to 380
          if(itemp.lt.0) go to 380
          nout=nout+1
          do i=1,nobs
            if(borenum(i).eq.itemp) then
              obsind(nout)=i
              go to 370
            end if
          end do
370       continue
        end do
380     continue
        if(all(obsind.eq.0))then
          write(amessage,384) trim(modfile),trim(obsfile)
384       format(' There are no observation well numbers listed in FEFLOW DAR file ',a,' for which ',   &
          'a corresponding bore identifier is provided in file ',a,'.')
          go to 9890
        end if

! -- We now read model outputs at observation wells from the rest of the file.

        modval=-1.1e35        ! an array
        itime=0
        do
          do
            iline=iline+1
            read(modunit,'(a)',end=500) cline
            if(index(cline,'***').ne.0) exit
          end do
          itime=itime+1
          iline=iline+1
          read(modunit,'(a)',end=500) cline
          n=index(cline,'TIME')
          if(n.eq.0) then
            call num2char(iline,aline)
            write(amessage,385) trim(aline),trim(modfile)
385         format(' "TIME =" header expected at line ',a,' of file ',a,'.')
            go to 9890
          end if
          cline(1:n+3)=' '
          cline=adjustl(cline)
          n=index(cline,'=')
          if(n.eq.0) then
            call num2char(iline,aline)
            write(amessage,385) trim(aline),trim(modfile)
            go to 9890
          end if
          cline(1:n)=' '
          call linesplit(ifail,1)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,390) trim(aline),trim(modfile)
390         format(' Cannot read simulation time from line ',a,' of file ',a,'.')
            go to 9890
          end if
          time=char2double(ifail,1)
          if(ifail.ne.0) then
            call num2char(iline,aline)
            write(amessage,390) trim(aline),trim(modfile)
            go to 9890
          end if
          if(time.lt.0.0)then
            call num2char(iline,aline)
            write(amessage,395) trim(aline),trim(modfile)
395         format(' Negative simulation time read from line ',a,' of file ',a,'.')
            go to 9890
          end if
          call linesplit(ifail,2)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,400) trim(aline),trim(modfile)
400         format(' Cannot read time units from line ',a,' of file ',a,'.')
            go to 9890
          end if
          atemp=cline(left_word(2):right_word(2))
          call casetrans(atemp,'lo')
          n=index(atemp,':')
          if(n.ne.0)then
            atemp(n:n)=' '
            atemp=adjustl(atemp)
          end if
          if(atemp.ne.'[d]')then
            call num2char(iline,aline)
            write(amessage,410) trim(aline),trim(modfile)
410         format(' Time units designator expected to be "[d]" at line ',a,' of file ',a,'.')
            go to 9890
          end if
          obsdays(itime)=time
          if(itime.eq.1)then
            iline=iline+1
            read(modunit,*,err=9200,end=9200)
            iline=iline+1
            read(modunit,'(a)',err=9200,end=9200) cline
            do
              n=index(cline,'[')
              if(n.eq.0) exit
              m=index(cline,']')
              if((m.eq.0).or.(m.lt.n))then
                call num2char(iline,aline)
                write(amessage,420) trim(aline),trim(modfile)
420             format(' Unbalanced brackets at line ',a,' of file ',a,'.')
                go to 9890
              end if
              cline(n:m)=' '
            end do
            icol=0
            do
              icol=icol+1
              call linesplit(ifail,icol)
              if(icol.eq.1)then
                if(ifail.ne.0)then
                  call num2char(iline,aline)
                  write(amessage,430) trim(aline),trim(modfile)
430               format(' First entry at line ',a,' of file ',a,' expected to be "obs".')
                  go to 9890
                else
                  atemp=cline(left_word(icol):right_word(icol))
                  call casetrans(atemp,'lo')
                  if(atemp.ne.'obs')then
                    call num2char(iline,aline)
                    write(amessage,430) trim(aline),trim(modfile)
                    go to 9890
                  end if
                end if
              else
                if(ifail.ne.0)then
                  call num2char(iline,aline)
                  write(amessage,440) trim(datatype),trim(aline),trim(modfile)
440               format(' Cannot find "',a,'" data header at line ',a,' of file ',a,'.')
                  go to 9890
                end if
                atemp=cline(left_word(icol):right_word(icol))
                call casetrans(atemp,'hi')
                if(atemp.eq.datatype) go to 450
              end if
            end do
450         ncol=icol
          else
            if(obsdays(itime).lt.obsdays(itime-1))then
              call num2char(iline,aline)
              write(amessage,451) trim(aline),trim(modfile)
451           format(' Simulation time recorded on line ',a,' of file ',a,' precedes ',  &
              'previous simulation time.')
              go to 9890
            else if (obsdays(itime).eq.obsdays(itime-1))then
              iwarn=1
              writetime(itime-1)=0
            end if
            do i=1,2
              iline=iline+1
              read(modunit,*,end=500)
            end do
          end if
          iline=iline+1
          read(modunit,*,end=500)
          do iout=1,nout
            iline=iline+1
            read(modunit,'(a)',end=500) cline
            i=obsind(iout)
            if(i.ne.0)then
              call linesplit(ifail,ncol)
              if(ifail.ne.0)then
                call num2char(iline,aline)
                write(amessage,460) trim(aline),trim(modfile)
460             format(' Insufficient entries on line ',a,' of file ',a,'.')
                go to 9890
              end if
              modval(itime,i)=char2double(ifail,ncol)
              if(ifail.ne.0)then
                call num2char(iline,aline)
                write(amessage,470) trim(datatype),trim(aline),trim(modfile)
470             format(' Error reading "',a,'" data type from line ',a,' of file ',a,'.')
                go to 9890
              end if
            end if
          end do
        end do

500     continue
        ntime=itime
        close(unit=modunit)
        call num2char(ntime,atemp)
        write(6,510) trim(atemp),trim(modfile)
510     format(' - data for ',a,' times read from file ',a,'.')

! -- The bore sample file is now written.

        write(6,540) trim(outfile)
540     format(/,' - writing bore sample file ',a,'...')
        do iobs=1,nobs
          do itime=1,ntime
            if(writetime(itime).ne.0)then
              if(modval(itime,iobs).gt.-1.0e35)then
                call elapsdate(obsdays(itime),day_convert,dds,mms,yys,hhhs,&
                mmms,ssss,day,mon,year,hour,min,sec)
                if(datespec.eq.1)then
                  write(outunit,550) trim(boreid(iobs)), &
                  day,mon,year,hour,min,sec,modval(itime,iobs)
                else
                  write(outunit,550) trim(boreid(iobs)), &
                  mon,day,year,hour,min,sec,modval(itime,iobs)
                end if
550             format(1x,a,t14,i2.2,'/',i2.2,'/',i4.4,t28,i2.2,':',i2.2,':',i2.2,t40,1pg14.7)
              end if
            end if
          end do
        end do
        close(unit=outunit)
        write(6,560) trim(outfile)
560     format(' - file ',a,' written ok.')
        if(iwarn.ne.0)then
          write(amessage,570) trim(modfile)
570       format(' Warning: at least one simulation time provided in file ',a,     &
          ' is equal to the previous simulation time provided in the same file. ', &
          ' Data for the first of two identical simulation times is ignored by DAR2SMP.')
          call write_message(leadspace='yes')
        end if


	go to 9900

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(obsfile)
9010    format(' Cannot read line ',a,' of file ',a,'.')
        go to 9890
9100    write(amessage,9110) trim(modfile)
9110    format(' Cannot read observation bore numbers from "LOCATION" section of file ',a,'.')
        go to 9890
9200    call num2char(iline,aline)
        write(amessage,9210) trim(modfile),trim(aline)
9210    format(' Error reading model output data block from file ',a,': error occurs at line ',a,'.')
        go to 9890

9890	call write_message(leadspace='yes')
9900	call close_files
        if(allocated(obsdays))   deallocate(obsdays,stat=ierr)
        if(allocated(obsind))    deallocate(obsind,stat=ierr)
        if(allocated(modval))    deallocate(modval,stat=ierr)
        if(allocated(writetime)) deallocate(writetime,stat=ierr)
	write(6,*)

end program dar2smp
