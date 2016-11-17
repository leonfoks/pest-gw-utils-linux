subroutine read_feflow_elem_prop_file(ifail,aprompt)

! -- Subroutine READ_FEFLOW_ELEM_PROP_FILE reads a FEFLOW element property file.

        use defn
        use inter
        implicit none

        integer, intent(out)            :: ifail
        character (len=*), intent(in)   :: aprompt

        integer                         :: iunit,iend,ierr,itemp,ielem,i,ilay,iline,icount
        double precision                :: dtemp
        character*10                    :: atemp,aline
        character*10                    :: atemp1,atemp2

! -- Initialisation.

        ifail=0
        imessage=0

! -- The name of the FEFLOW element property file is acquired.

        icount=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
        read(5,'(a)') cline
        if(cline.eq.' ') go to 10
        cline=adjustl(cline)
        if(index(eschar,cline(1:2)).ne.0) then
          escset=1
          return
        end if
        iend=len_trim(cline)
        call getfile(ifail,cline,epfile_f,1,iend)
        if(ifail.ne.0) go to 10

! -- The FEFLOW element property file is opened.

	iunit=nextunit()
	open(unit=iunit,file=epfile_f,status='old',iostat=ierr)
	if(ierr.ne.0)then
          icount=icount+1
          if(icount.eq.5) go to 9891
	  write(6,20) trim(epfile_f)
20        format(/,' Cannot open file ',a,' - try again.',/)
	  go to 10
	end if
        imessage=0

! -- The first line of the file is checked for the correct header.

        read(iunit,'(a)',iostat=ierr) cline
        if(ierr.ne.0)then
          write(amessage,32) trim(epfile_f)
32        format(' Error reading first line of file ',a,'.')
          go to 9890
        end if
        call casetrans(cline,'hi')
        call linesplit(ifail,6)
        if(ifail.ne.0)then
          write(amessage,35) trim(epfile_f)
35        format(' Unexpected header at first line of file ',a,'.')
          go to 9890
        end if
        if((cline(left_word(1):right_word(1)).ne.'ELEM').or.   &
           (cline(left_word(2):right_word(2)).ne.'LAYER').or.  &
           (cline(left_word(3):right_word(3)).ne.'X').or.      &
           (cline(left_word(4):right_word(4)).ne.'Y').or.      &
           (cline(left_word(5):right_word(5)).ne.'Z').or.      &
           (cline(left_word(6):right_word(6)).ne.'F'))then
           write(amessage,35) trim(epfile_f)
           go to 9890
        end if

        allocate(eastelem_f(numelem_f),northelem_f(numelem_f),zoneelem_f(numelem_f), &
        stat=ierr)
        if(ierr.ne.0)then
          write(amessage,111)
111       format(' Cannot allocate sufficient memory to continue execution.')
          go to 9890
        end if

! -- Element data are now read.

        iline=1
        do ielem=1,numelem_f
109       iline=iline+1
          read(iunit,'(a)',err=115,end=116) cline
          go to 119
115       call num2char(iline,atemp)
          write(amessage,112) trim(atemp),trim(epfile_f)
112       format(' Error reading data at line ',a,' of file ',a,'.')
          go to 9890
116       write(amessage,117) trim(epfile_f)
117       format(' End encountered to file ',a,' before data read for all elements.')
          go to 9890
119       continue
          if(cline.eq.' ') go to 109
          call linesplit(ifail,6)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,130) trim(aline),trim(epfile_f)
130         format(' Insufficient items on line ',a,' of file ',a,'. Six entries ',  &
            'are expected.')
            go to 9890
          end if
          i=char2int(ifail,1)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,140) trim(aline),trim(epfile_f)
140         format(' Cannot read element number from line ',a,' of file ',a,'.')
            go to 9890
          end if
          if(i.ne.ielem)then
            call num2char(ielem,atemp1)
            call num2char(i,atemp2)
            call num2char(iline,aline)
            write(amessage,113) trim(atemp2),trim(atemp1),trim(aline),trim(epfile_f)
113         format(' Data is provided for element ',a,' when data is expected for element ',   &
            a,' at line ',a,' of file ',a,'.')
            go to 9890
          end if
          ilay=char2int(ifail,2)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,150) trim(aline),trim(epfile_f)
150         format(' Cannot read layer number from line ',a,' of file ',a,'.')
            go to 9890
          end if
          eastelem_f(ielem)=char2double(ifail,3)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,160) trim(aline),trim(epfile_f)
160         format(' Cannot read element easting from line ',a,' of file ',a,'.')
            go to 9890
          end if
          northelem_f(ielem)=char2double(ifail,4)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,170) trim(aline),trim(epfile_f)
170         format(' Cannot read element northing from line ',a,' of file ',a,'.')
            go to 9890
          end if
          dtemp=char2double(ifail,5)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,180) trim(aline),trim(epfile_f)
180         format(' Cannot read element "Z" value from line ',a,' of file ',a,'.')
            go to 9890
          end if
          dtemp=char2double(ifail,6)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,190) trim(aline),trim(epfile_f)
190         format(' Cannot read element "F" value from line ',a,' of file ',a,'.')
            go to 9890
          end if
          zoneelem_f(ielem)=nint(dtemp)
        end do

! -- Wrapping up.

        call num2char(numelem_f,atemp)
        write(amessage,220) trim(atemp),trim(epfile_f)
220     format('  - data for ',a,' elements read from element property file ',a)
        call write_message
        go to 9900

9890    call write_message(leadspace='yes')
9891    ifail=1

9900    continue
        close(unit=iunit,iostat=ierr)
        imessage=0
        return

end subroutine read_feflow_elem_prop_file




