subroutine read_gms_real_mesh_data_file(ifail,aprompt,realdata)

! -- Subroutine READ_GMS_REAL_MESH_DATA_FILE reads a GMS mesh data file in which
!    stored values are real numbers.

        use defn
        use inter
        implicit none

        integer, intent(out)             :: ifail
        character (len=*), intent(inout) :: aprompt
        real, intent(out)                :: realdata(:)

        integer                          :: iline,iunit,iend,i,ierr,itemp,iitemp,j
        real                             :: rtemp
        character*10                     :: aline
        character*5                      :: atemp
        character*120                    :: afile

! -- Initialisation.

        iline=0
        ifail=0
        imessage=0

! -- The name of the GMS data file is acquired.

10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
        read(5,'(a)') cline
        if(cline.eq.' ') go to 10
        cline=adjustl(cline)
        if(index(eschar,cline(1:2)).ne.0) then
          escset=1
          return
        end if
        iend=len_trim(cline)
        call getfile(ifail,cline,afile,1,iend)
        if(ifail.ne.0) go to 10

! -- The GMS data file is opened.

        iunit=nextunit()
        open(unit=iunit,file=afile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          if(imessage.eq.5) go to 9890
          write(amessage,20) trim(afile)
20        format('cannot open index file ',a,' - try again.')
          call write_message(error='yes',increment=1)
          go to 10
        end if
        imessage=0
        write(initial_message,30) trim(afile)
30      format(' Errors in index file ',a,' ----->')

! -- The first line of the file is checked for the correct header.

        iline=iline+1
        read(iunit,'(a)',err=9000,end=9100) cline
        call casetrans(cline,'hi')
        if(cline(1:7).ne.'DATASET')then
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          write(amessage,32) trim(afile)
32        format(' First line of file ',a,' expected to contain "DATASET" header.')
          call write_message(increment=1)
          go to 9890
        end if

! -- The "ND" keyword is searched for.

        do
          iline=iline+1
          read(iunit,'(a)',err=9000,end=35) cline
          if(cline.eq.' ') cycle
          call linesplit(ifail,1)
          atemp=cline(left_word(1):right_word(1))
          call casetrans(atemp,'hi')
          if(atemp.eq.'ND') go to 38
        end do

35      continue
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        write(amessage,36) trim(afile)
36      format(' Cannot find "ND" keyword in file ',a,'.')
        call write_message(increment=1)
        go to 9890

38      continue
        call linesplit(ifail,2)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,39) trim(aline),trim(afile)
39        format(' Insufficient entries on line ',a,' of file ',a,'.')
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          call write_message(increment=1)
          go to 9890
        end if
        itemp=char2int(ifail,2)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,42) trim(aline),trim(afile)
42        format(' Cannot read integer following "ND" keyword at line ',a,' of file ',a,'.')
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          call write_message(increment=1)
          go to 9890
        end if
        if(itemp.ne.numelem_g)then
          call num2char(iline,aline)
          write(amessage,43) trim(aline),trim(afile)
43        format(' Integer following "ND" keyword at line ',a,' of file ',a,' is not the ', &
          'same as the number of elements in the mesh.')
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          call write_message(increment=1)
          go to 9890
        end if

! -- The data is now searched for.

        do
          iline=iline+1
          read(iunit,'(a)',err=9000,end=25) cline
          if(cline.eq.' ') cycle
          call linesplit(ifail,1)
          atemp=cline(left_word(1):right_word(1))
          call casetrans(atemp,'hi')
          if(atemp.eq.'TS') go to 50
        end do

25      continue
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        write(amessage,26) trim(afile)
26      format(' Cannot find "TS" keyword in file ',a,'.')
        call write_message(increment=1)
        go to 9890

! -- The actual data is now read.

50      continue
        do i=1,numelem_g
51        continue
          iline=iline+1
          read(iunit,'(a)',err=9000,end=9100) cline
          if(cline.eq.' ') go to 51
          call linesplit(ifail,1)
          realdata(elemindex_g(i))=char2real(ifail,1)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(iline,aline)
            write(amessage,52) trim(aline),trim(afile)
52          format(' Error reading data from line ',a,' of file ',a,'.')
            call write_message(increment=1)
          end if
        end do
        if(imessage.ne.0) then
          go to 9890
        else
          aprompt=afile
          write(amessage,220) trim(afile)
220       format('  - file ',a,' read ok')
          call write_message
          go to 9900
        end if

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error encountered when reading line ',a,' of file ',a,'.')
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        call write_message(increment=1)
        go to 9890

9100    write(amessage,9110) trim(afile)
9110    format(' Unexpected end encountered to file ',a,'.')
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        call write_message(increment=1)
        go to 9890

9890    ifail=1

9900    continue
        close(unit=iunit,iostat=ierr)
        imessage=0
        return

end subroutine read_gms_real_mesh_data_file
   
