!     Last change:  J    14 Jun 2002   11:09 pm
subroutine read_real_tab_file_int(ifail,mode,tabunit,tabfile,realarray, &
intarray,gridspec,insunit,insfile,aprefix)

! -- Subroutine read_real_tab_file_int reads a real array table file, checking
!    it for errors and transferring the results to a real array. It records cells
!    represented in the real array file to an integer array.

! -- Arguments are as follows:-
!       ifail:     returned as non-zero if error condition encountered
!       mode:      1 if fill real array, 2 if just integer array, 3 to write
!                  instruction file
!       tabunit:   unit number of real array table file
!       tabfile:   real array table filename
!       realarray: real array
!       intarray:  integer array
!       gridspec:  defined type containing grid specifications
!       insunit:   unit number of instruction file
!       insfile:   an instruction file to be written
!       aprefix:   prefix for observation names


        use defn
        use inter
        implicit none

        integer, intent(out)                    :: ifail
        integer, intent(in)                     :: mode
        integer, intent(in)                     :: tabunit
        character(len=*), intent(in)            :: tabfile
        real, dimension(:,:),intent(inout)      :: realarray
        integer, dimension(:,:),intent(inout)   :: intarray
        type(modelgrid), intent(in)             :: gridspec
        integer, intent(in)                     :: insunit
        character(len=*), intent(in)            :: insfile,aprefix

        integer                                 :: icol,irow,iline
        real                                    :: rval
        character (len=10)                      :: aline
        character (len=5)                       :: acol,arow
        character (len=12)                      :: aobs

        if(mode.eq.3)then
          write(insunit,10,err=9300)
10        format('pif $')
        end if

        iline=0
        imessage=0
30      iline=iline+1
        read(tabunit,'(a)',err=9000,end=1000) cline
        call linesplit(ifail,3)
        if(ifail.lt.0) go to 30
        if(ifail.gt.0) then
          if(imessage.eq.0) call write_initial_message(leadspace='no')
          call num2char(iline,aline)
          write(amessage,50) trim(aline), trim(tabfile)
50        format(' Insufficient items on line ',a,' of file ',a,  &
          ': 3 expected.')
          call write_message(increment=1)
          go to 30
        end if

        irow=char2int(ifail,1)
        if(ifail.ne.0) then
          if(imessage.eq.0) call write_initial_message(leadspace='no')
          call num2char(iline,aline)
          write(amessage,70) trim(aline),trim(tabfile)
70        format(' Cannot read grid row number from line ',a,' of file ',a,'.')
          call write_message(increment=1)
        else
          if(irow.le.0) then
            if(imessage.eq.0) call write_initial_message(leadspace='no')
            call num2char(iline,aline)
            write(amessage,90) trim(aline),trim(tabfile)
90          format(' Grid row number less than one at line ',a,' of file ',a,'.')
            call write_message(increment=1)
          else if(irow.gt.gridspec%nrow) then
            if(imessage.eq.0) call write_initial_message(leadspace='no')
            call num2char(iline,aline)
            write(amessage,110) trim(aline),trim(tabfile)
110         format(' Grid row number exceeds grid row upper bound at line ',  &
            a,' of file ',a,'.')
            call write_message(increment=1)
          end if
        end if

        icol=char2int(ifail,2)
        if(ifail.ne.0) then
          if(imessage.eq.0) call write_initial_message(leadspace='no')
          call num2char(iline,aline)
          write(amessage,130) trim(aline),trim(tabfile)
130       format(' Cannot read grid column number at line ',a,' of file ',a,'.')
          call write_message(increment=1)
        else
          if(icol.le.0) then
            if(imessage.eq.0) call write_initial_message(leadspace='no')
            call num2char(iline,aline)
            write(amessage,150) trim(aline),trim(tabfile)
150         format(' Grid column number less than one at line ',a,' of file ',a,'.')
            call write_message(increment=1)
          else if(icol.gt.gridspec%ncol) then
            if(imessage.eq.0) call write_initial_message(leadspace='no')
            call num2char(iline,aline)
            write(amessage,160) trim(aline),trim(tabfile)
160         format(' Grid column number exceeds grid column ',&
            'upper bound at line ',a,' of file ',a,'.')
            call write_message(increment=1)
          end if
        end if

        if(mode.eq.3)then
          call num2char(irow,arow)
          call num2char(icol,acol)
          aobs=trim(aprefix)//'_'//trim(arow)//'_'//trim(acol)
          write(insunit,170) trim(aobs)
170       format('l1  [',a,']16:33')
        else
          rval=char2real(ifail,3)
          if(ifail.ne.0) then
            if(imessage.eq.0) call write_initial_message(leadspace='no')
            call num2char(iline,aline)
            write(amessage,180) trim(aline),trim(tabfile)
180         format(' Cannot read real array value from line ',a,' of file ',a,'.')
            call write_message(increment=1)
          else
            if(imessage.eq.0)then
              if(mode.eq.1)then
                realarray(icol,irow)=rval
              end if
              intarray(icol,irow)=1
            end if
          end if
        end if
        go to 30

1000    if(imessage.eq.0) then
          if(mode.ne.3)then
            iline=count(intarray.ne.0)
            call num2char(iline,aline)
            write(amessage,1030) trim(aline),trim(tabfile)
1030        format(' - ',a,' values read from real array table file ',a)
            call write_message()
          end if
          ifail=0
          return
        else
          ifail=1
          return
        end if

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(tabfile)
9010    format(' Cannot read line ',a,' of real array table file ',a)
        call write_message(leadspace='yes')
        ifail=1
        return

9300    write(amessage,9310) trim(insfile)
9310    format(' Cannot write to instruction file ',a,'.')
        ifail=1
        return

end subroutine read_real_tab_file_int

