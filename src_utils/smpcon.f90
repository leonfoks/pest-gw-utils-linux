!     Last change:  JD    9 May 2003   11:50 am
program smpcon

! -- Program SMPCON converts a bore/site sample file between binary and ascii
!    format.

	use defn
	use inter
	implicit none

        logical              :: lexist
	integer              :: ifail,idate,itry,inunit,ierr, &
                                iline,cols,ndays,nsecs,idays, &
                                dd,mm,yy,hhh,mmm,sss,iheader,outunit

	double precision     :: value

        character (len=3)    :: acode
	character (len=5)    :: anum
        character (len=12)   :: atemp1
        character (len=15)   :: adate,atime,aline
        character (len=23)   :: bcode
	character (len=200)  :: infile,outfile



	write(amessage,5)
5	format(' Program SMPCON converts bore/site sample data between ', &
        'ASCII and binary format.')
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
	if((idate.ne.0).or.(datespec.eq.0)) then
	  write(amessage,9)
9	  format(' Cannot read date format from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

24      write(6,15,advance='no')
15      format(' Convert ascii-to-binary or binary-to-ascii? [a2b/b2a]: ')
        read(5,'(a)') acode
        call casetrans(acode,'lo')
	if(acode.eq.' ') go to 24
	acode=adjustl(acode)
	if(index(eschar,acode(1:2)).ne.0) go to 9900
        if((acode.ne.'a2b').and.(acode.ne.'b2a')) go to 24

21      itry=0
        if(acode.eq.'a2b')then
22        continue
          if(itry.gt.5) go to 9900
          write(6,23,advance='no')
23        format(' Enter name of existing ASCII site sample file: ')
          read(5,'(a)') infile
          if(index(eschar,infile(1:2)).ne.0) then
            write(6,*)
            go to 24
          end if
          inquire(file=infile,exist=lexist)
          if(.not.lexist)then
            write(amessage,28) trim(infile)
28          format(' File ',a,' does not exist - try again.')
            call write_message()
            itry=itry+1
            go to 22
          end if
          inunit=nextunit()
          open(unit=inunit,file=infile,form='formatted',status='old',   &
          iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,29)
29          format(' Cannot open file ',a,'. Is this file truly an ASCII file? ', &
            'Is it being used by another program?')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            go to 22
          end if
        else
32        continue
          if(itry.gt.5) go to 9900
          write(6,33,advance='no')
33        format(' Enter name of existing binary site sample file: ')
          read(5,'(a)') infile
          if(index(eschar,infile(1:2)).ne.0) then
            write(6,*)
            go to 24
          end if
          inquire(file=infile,exist=lexist)
          if(.not.lexist)then
            write(amessage,28) trim(infile)
            call write_message()
            itry=itry+1
            go to 32
          end if
          inunit=nextunit()
          open(unit=inunit,file=infile,form='binary',status='old',   &
          iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,39)
39          format(' Cannot open file ',a,'. Is this file truly a binary file? ', &
            'Is it being used by another program?')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            go to 32
          end if
          read(inunit,iostat=ierr) bcode
          if(ierr.ne.0)then
            write(amessage,42) trim(infile)
42          format(' Cannot read header to file binary site sample file ',a,  &
            ': are you sure that this is a file of the correct type? Try again.')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            close(unit=inunit)
            go to 32
          end if
          if(bcode.ne.'binary_site_sample_file')then
            write(amessage,45) trim(infile)
45          format(' Incorrect header to file ',a,': are you sure that this is a ', &
            'binary site sample file? Try again.')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            close(unit=inunit)
            go to 32
          end if
        end if

        write(6,*)
49      itry=0
        if(acode.eq.'a2b')then
50        if(itry.gt.5) go to 9900
          write(6,60,advance='no')
60        format(' Enter name for new binary output file: ')
          read(5,'(a)') outfile
          if(index(eschar,outfile(1:2)).ne.0) then
            write(6,*)
            close(unit=inunit)
            go to 21
          end if
          outunit=nextunit()
          open(unit=outunit,file=outfile,form='binary', &
          iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,70) trim(outfile)
70          format(' Cannot open file ',a,' for output. Try again.')
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            go to 50
          end if
          write(outunit) 'binary_site_sample_file'
        else
130       if(itry.gt.5) go to 9900
          write(6,140,advance='no')
140       format(' Enter name for new ASCII output file: ')
          read(5,'(a)') outfile
          if(index(eschar,outfile(1:2)).ne.0) then
            write(6,*)
            close(unit=inunit)
            go to 21
          end if
          outunit=nextunit()
          open(unit=outunit,file=outfile,iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,70) trim(outfile)
            call write_message(leadspace='yes',endspace='yes')
            itry=itry+1
            go to 130
          end if
        end if

! -- The output file is now written.

        iline=0
        if(acode.eq.'a2b')then
	  read_asc_smp_file: do
	    iline=iline+1
	    read(inunit,'(a)',end=500) cline
	    cols=5
	    call linesplit(ifail,5)
	    if(ifail.lt.0) cycle read_asc_smp_file
	    if(ifail.gt.0)then
	      cols=4
	      call linesplit(ifail,4)
	      if(ifail.ne.0) then
	        call num2char(iline,aline)
	        write(amessage,150) trim(aline),trim(infile)
150             format('Insufficient entries on line ',a,' of site sample file ',a)
	        go to 9890
	      end if
	    end if
	    atemp1=cline(left_word(1):right_word(1))
	    call casetrans(atemp1,'hi')
	    call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value, &
	    iline,infile)
	    if(ifail.ne.0) go to 9900
	    if(value.lt.-1.0e38) cycle read_asc_smp_file
            write(outunit)atemp1,ndays,nsecs,value
          end do read_asc_smp_file
        else
	  read_bin_smp_file: do
	    iline=iline+1
	    read(inunit,end=500,iostat=ierr) atemp1,ndays,nsecs,value
            if(ierr.ne.0)then
              call num2char(iline,aline)
              write(amessage,160) trim(aline),trim(infile)
160           format('Cannot read record ',a,' of site sample file ',a)
              go to 9890
            end if
	    call casetrans(atemp1,'hi')
	    call newdate(ndays,1,1,1970,dd,mm,yy)
	    hhh=nsecs/3600
	    mmm=(nsecs-hhh*3600)/60
	    sss=nsecs-hhh*3600-mmm*60
	    if(datespec.eq.1) then
	      write(outunit,170) trim(atemp1),dd,mm,yy,hhh,mmm,sss,value
170	      format(1x,a,t15,i2.2,'/',i2.2,'/',i4.4,3x,i2.2,':',i2.2,':',   &
	      i2.2,3x,1pg15.8)
	    else
	      write(outunit,170) trim(atemp1),mm,dd,yy,hhh,mmm,sss,value
	    endif
          end do read_bin_smp_file
        end if

500     write(6,*)
        write(6,510) trim(infile)
510     format(' - file ',a,' read ok.')
        write(6,520) trim(outfile)
520     format(' - file ',a,' written ok.')

        go to 9900

9890    call write_message(leadspace='yes')
9900    call close_files
        write(6,*)

end program smpcon

