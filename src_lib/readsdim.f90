!     Last change:  JD   29 Jan 2001    6:12 pm


subroutine read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)

! -- Subroutine READ_STRUCTURE_FILE_DIM reads a structure file, obtaining the
!    size of arrays that must be dimensioned in the main program.

! -- Arguments are as follows:-
!        ifail:       returned as zero unless an error condition is encountered
!        structunit:  unit number from which structure file is read
!        numstruct:   number of geostatistical structure definitions in file
!        numvario:    number of variogram definitions in file
!        structfile:  name of structure file

        use defn
        use inter
        implicit none

        integer, intent(out)            :: ifail
        integer, intent(in)             :: structunit
        integer, intent(out)            :: numstruct,numvario
        character (len=*), intent(in)   :: structfile

        integer                         :: istructure,ivariogram,iline,lifail, &
                                           ishortflag
        character (len=5)               :: aline
        character (len=20)              :: atemp,aatemp

        ifail=0
        imessage=0
        istructure=0
        ivariogram=0
        numstruct=0
        numvario=0

	write(initial_message,5) trim(structfile)
5       format(' Errors in structure file ',a,' ----->')

        iline=0
10      iline=iline+1
        read(structunit,'(a)',err=9000,end=1000) cline
        if(cline.eq.' ') go to 10
        if(cline(1:1).eq.'#') go to 10
        ishortflag=0
        call linesplit(lifail,2)
        if(lifail.ne.0)then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  call num2char(iline,aline)
	  write(amessage,20) trim(aline)
20        format('   Line ',a,': insufficient entries.')
	  call write_message(increment=1)
	  call linesplit(lifail,1)
          ishortflag=1
        end if
        atemp=cline(left_word(1):right_word(1))
        call casetrans(atemp,'hi')
        if(atemp.eq.'STRUCTURE')then
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(iline,aline)
            write(amessage,30) trim(aline)
30          format('   Unexpected STRUCTURE keyword at line ',a)
            call write_message(increment=1)
            go to 9990
          else if(ivariogram.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(iline,aline)
            write(amessage,30) trim(aline)
            call write_message(increment=1)
            go to 9990
          else
            istructure=1
            numstruct=numstruct+1
          end if
        else if(atemp.eq.'VARIOGRAM')then
          if(istructure.eq.0)then
            if(ivariogram.eq.1)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              call num2char(iline,aline)
              write(amessage,40) trim(aline)
40            format('   Unexpected VARIOGRAM keyword at line ',a)
              call write_message(increment=1)
              go to 9990
            else
              ivariogram=1
              numvario=numvario+1
            end if
          end if
        else if(atemp.eq.'END')then
          if(istructure.eq.1)then
            if(ishortflag.ne.1)then
              aatemp=cline(left_word(2):right_word(2))
              call casetrans(aatemp,'hi')
              if(aatemp.ne.'STRUCTURE')then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call num2char(iline,aline)
                write(amessage,50) trim(aline)
50              format('   Line ',a,' should read "END STRUCTURE".')
                call write_message(increment=1)
              end if
            end if
            istructure=0
          else if(ivariogram.eq.1)then
            if(ishortflag.ne.1)then
              aatemp=cline(left_word(2):right_word(2))
              call casetrans(aatemp,'hi')
              if(aatemp.ne.'VARIOGRAM')then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call num2char(iline,aline)
                write(amessage,55) trim(aline)
55              format('   Line ',a,' should read "END VARIOGRAM".')
                call write_message(increment=1)
              end if
            end if
            ivariogram=0
          else
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(iline,aline)
            write(amessage,60) trim(aline)
60          format('   Unexpected "END" keyword encountered on line ',a)
            call write_message(increment=1)
            go to 9990
          end if
        end if
        go to 10

1000    if(istructure.eq.1)then
          write(amessage,1010) trim(structfile)
1010      format('   Unexpected end to file ',a,' encountered while reading structure ',&
          'specifications.')
          go to 9890
        end if
        if(ivariogram.eq.1)then
          write(amessage,1020) trim(structfile)
1020      format('   Unexpected end to file ',a,' encountered while reading variogram ',&
          'specifications.')
          go to 9890
        end if
        if(imessage.ne.0) go to 9990
	rewind(unit=structunit,err=9100)
        go to 9999


9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(structfile)
9010    format('Error reading line ',a,' of structure file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990
9100    write(amessage,9110) trim(structfile)
9110    format(' Cannot rewind structure file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990

9890    call write_message(leadspace='yes')
9990    ifail=1
        close(unit=structunit)
9999    continue
        imessage=0
        return


end subroutine read_structure_file_dim
