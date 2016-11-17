!     Last change:  JD   14 Jan 2002   10:35 am


subroutine read_parameter_replacement_file(ifail,aprompt,ncol,nrow)

! -- Subroutine read_parameter_replacement_file reads a MODFLOW 2000 parameter
!    replacement file, checking it for errors.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless an error condition arises
!       prompt:    contains text used to prompt for filename
!       ncol,nrow : number of columns and rows in the finite-difference grid.

	use defn
	use inter

	integer, intent(out)            :: ifail
	character (len=*), intent(in)   :: aprompt
        integer, intent(in)             :: ncol,nrow
	integer                         :: iunit,ierr,iline,i,iend,jrep, &
                                           ireplace,irep,itemp1,itemp2,nbb, &
                                           izone,imult,inumout,icol,irow,mcol,mrow
        character (len=10)              :: anum
	character (len=20)              :: aline,akey
	character (len=20)              :: atemp1,atemp2,atemp3
        character (len=120)             :: afile,bfile

	if(associated(replace)) deallocate (replace)
	if(associated(removezone)) deallocate (removezone)
	if(associated(removemult)) deallocate (removemult)
	nullify(replace,removezone,removemult)


! -- Some variables are initialised.

        numparamreplace=0
        numzoneremove=0
        nummultremove=0
        ireplace=0

        imessage=0
        ifail=0

! -- The name of the file is obtained and the file is opened.

10	write(6,'(a)',advance='no') aprompt
	read(5,'(a)') afile
        if(afile.eq.' ') go to 10
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0) then
	  escset=1
	  return
	end if
        iend=len_trim(afile)
        call getfile(ifail,afile,repfile,1,iend)
        if(ifail.ne.0) go to 10

	iunit=nextunit()
	open(unit=iunit,file=repfile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9990
	  write(amessage,20) trim(repfile)
20        format(' Cannot open parameter replacement file ',a,' - try again.')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 10
	end if
	write(initial_message,30) trim(repfile)
30      format(' Errors in parameter replacement file ',a,' ----->')

! -- The file is now read a first time in order to dimension arrays.

        iline=0
40      iline=iline+1
        call num2char(iline,aline)
        read(iunit,'(a)',err=9000,end=55) cline
        if(cline.eq.' ') go to 40
        if(cline(1:1).eq.'#') go to 40
        call linesplit(ifail,2)
        if(ifail.ne.0)then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,50) trim(aline)
50        format('   Insufficient entries on line ',a)
	  call write_message(increment=1)
	  go to 9990
        end if
        atemp1=cline(left_word(1):right_word(1))
        call casetrans(atemp1,'hi')
        if(atemp1.eq.'REPLACE')then
          if(ireplace.eq.1)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='REPLACE'
            write(amessage,60) trim(akey),trim(aline)
60          format('   Unexpected ',a,' keyword at line ',a)
            call write_message(increment=1)
            go to 9990
          end if
          ireplace=1
          numparamreplace=numparamreplace+1
        else if(atemp1.eq.'END')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='END'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          ireplace=0
        else if(atemp1.eq.'REMOVEZONE')then
          numzoneremove=numzoneremove+1
        else if(atemp1.eq.'REMOVEMULT')then
          nummultremove=nummultremove+1
        end if
        go to 40

55      continue

        if(imessage.ne.0) go to 9990

        allocate(replace(numparamreplace),removezone(numzoneremove),  &
        removemult(nummultremove),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,80)
80        format(' Cannot allocate sufficient memory to continue execution.')
          go to 9890
        end if

        rewind(unit=iunit,iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,90) trim(repfile)
90        format(' Cannot rewind file ',a)
          go to 9890
        end if

! -- The replacement file is now read a second time and values are given to pertinent
!    variables.

! -- But first, a little more initialisation is done.

        do irep=1,numparamreplace
          replace(irep)%numout=0
          replace(irep)%numin=0
          replace(irep)%newparamprefix=' '
          replace(irep)%transformtype=-9999
          replace(irep)%replacetype=' '
          do i=1,MAX_OUT_PARAM
            replace(irep)%outparamname(i)=' '
            replace(irep)%outparamlay1(i)=-9999
            replace(irep)%outparamlay2(i)=-9999
          end do
          do i=1,MAX_IN_PARAM
            replace(irep)%newparamlay1(i)=-9999
            replace(irep)%newparamlay2(i)=-9999
          end do
          replace(irep)%minival=-1.1e35
          replace(irep)%maxival=-1.1e35
          replace(irep)%minifile=' '
          replace(irep)%maxifile=' '
          replace(irep)%minpval=-1.1e35
          replace(irep)%maxpval=-1.1e35
        end do

! -- Now the file is re-read.

        irep=0
        izone=0
        imult=0
        iline=0
100     iline=iline+1
        call num2char(iline,aline)
        read(iunit,'(a)',err=9000,end=245) cline
        if(cline.eq.' ') go to 100
        if(cline(1:1).eq.'#') go to 100
        call linesplit(ifail,2)
        atemp1=cline(left_word(1):right_word(1))
        atemp2=cline(left_word(2):right_word(2))
        call casetrans(atemp1,'hi')
        if(atemp1.eq.'REPLACE')then
          ireplace=1
          irep=irep+1
          call casetrans(atemp2,'hi')
          if((atemp2.ne.'HK').and.(atemp2.ne.'HANI').and.(atemp2.ne.'VK').and. &
             (atemp2.ne.'VANI').and.(atemp2.ne.'SS').and.(atemp2.ne.'SY').and. &
             (atemp2.ne.'VKCB').and.(atemp2.ne.'RCH'))then
	     if(imessage.eq.0) call write_initial_message(leadspace='yes')
             write(amessage,110) trim(aline)
110          format('   Parameter type must be "HK", "HANI", "VK", "VANI", "SS", "SY", ', &
             '"VKCB" or "RCH" at line ',a)
             call write_message(increment=1)
           else
             replace(irep)%replacetype=atemp2
           end if
        else if(atemp1.eq.'OUT')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='OUT'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          replace(irep)%numout=replace(irep)%numout+1
          if(replace(irep)%numout.gt.MAX_OUT_PARAM)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,115)
115         format('   Too many parameters cited as "OUT" - contact programmer.')
            call write_message(increment=1)
            go to 9990
          end if
          call casetrans(atemp2,'hi')
          inumout=replace(irep)%numout
          replace(irep)%outparamname(inumout)=atemp2
!          call linesplit(ifail,4)
!          if(ifail.ne.0)then
!	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
!            write(amessage,120) trim(aline)
!120         format('   4 entries expected on line ',a)
!            call write_message(increment=1)
!            go to 100
!          end if
!          itemp1=char2int(ifail,3)
!          if(ifail.ne.0)then
!	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
!            write(amessage,130) trim(aline)
130         format('   Cannot read first layer number at line ',a)
!            call write_message(increment=1)
!          end if
!          if(imessage.eq.0)then
!            itemp2=replace(irep)%numout
!            replace(irep)%outparamlay1(itemp2)=itemp1
!          end if
!          itemp1= char2int(ifail,4)
!          if(ifail.ne.0)then
!	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
!            write(amessage,140) trim(aline)
140         format('   Cannot read second layer number at line ',a)
!            call write_message(increment=1)
!          end if
!          if(imessage.eq.0)then
!            itemp2=replace(irep)%numout
!            replace(irep)%outparamlay2(itemp2)=itemp1
!            if(itemp1.lt.replace(irep)%outparamlay1(itemp2))then
!               if(imessage.eq.0) call write_initial_message(leadspace='yes')
!               write(amessage,155) trim(aline)
!               call write_message(increment=1)
!             end if
!             if((itemp1.le.0).or. &
!                (replace(irep)%outparamlay1(itemp2).le.0))then
!                if(imessage.eq.0) call write_initial_message(leadspace='yes')
!                write(amessage,141) trim(aline)
!141             format('   Invalid layer number at line ',a)
!                call write_message(increment=1)
!             end if
!          end if

        else if((atemp1(1:11).eq.'NEWPARAMLAY').or.   &
                (atemp1(1:11).eq.'NEWPARAMPER')) then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey=atemp1(1:11)//'S'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          if(replace(irep)%replacetype.eq.'RCH')then
            if(atemp1(1:11).eq.'NEWPARAMLAY')then
 	      if(imessage.eq.0) call write_initial_message(leadspace='yes')
              akey='NEWPARAMLAYS'
              write(amessage,60) trim(akey),trim(aline)
              call write_message(increment=1)
              go to 100
            end if
          else
            if(atemp1(1:11).eq.'NEWPARAMPER')then
 	      if(imessage.eq.0) call write_initial_message(leadspace='yes')
              akey='NEWPARAMPERS'
              write(amessage,60) trim(akey),trim(aline)
              call write_message(increment=1)
              go to 100
            end if
          end if
          replace(irep)%numin=replace(irep)%numin+1
          if(replace(irep)%numin.gt.MAX_IN_PARAM)then
            write(amessage,145)
145         format('   Too many "NEWPARAMLAY" keywords - contact programmer.')
            go to 9890
          end if
          call linesplit(ifail,3)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(aline)
150         format('   3 entries expected on line ',a)
            call write_message(increment=1)
            go to 100
          end if
          itemp1=char2int(ifail,2)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,130) trim(aline)
            call write_message(increment=1)
          else
            itemp2=replace(irep)%numin
            replace(irep)%newparamlay1(itemp2)=itemp1
          end if
          itemp1= char2int(ifail,3)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,140) trim(aline)
            call write_message(increment=1)
          else
            itemp2=replace(irep)%numin
            replace(irep)%newparamlay2(itemp2)=itemp1
          end if
          if(imessage.eq.0)then
            if(replace(irep)%newparamlay2(itemp2).lt.  &
               replace(irep)%newparamlay1(itemp2))then
                 if(imessage.eq.0) call write_initial_message(leadspace='yes')
                 if(replace(irep)%replacetype.eq.'RCH')then
                   write(amessage,154) trim(aline)
154                format(' Second stress period less than first at line ',a)
                 else
                   write(amessage,155) trim(aline)
155                format(' Second layer less than first layer at line ',a)
                 endif
                 call write_message(increment=1)
            end if
          end if
        else if(atemp1.eq.'NEWPARAMPREFIX')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='NEWPARAMPREFIX'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          nbb=len_trim(cline)
          call getfile(ifail,cline,atemp2,left_word(2),nbb)
          if(ifail.ne.0)then
            if((cline(left_word(2):left_word(2)+1).eq.'""').or.  &
               (cline(left_word(2):left_word(2)+1).eq.''''''))then
               replace(irep)%newparamprefix=' '
             else
	       if(imessage.eq.0) call write_initial_message(leadspace='yes')
               write(amessage,170) trim(aline)
170            format('   Cannot read parameter name prefix on line ',a)
              call write_message(increment=1)
            end if
          else
            if(len_trim(atemp2).gt.2)then
	      if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,180) trim(aline)
180           format('   Parameter name prefix must be 2 characters or less at line ',a)
              call write_message(increment=1)
            else
              call casetrans(atemp2,'hi')
              replace(irep)%newparamprefix=atemp2
            end if
          end if
        else if(atemp1.eq.'TRANSFORMTYPE')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='TRANSFORMTYPE'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          call casetrans(atemp2,'lo')
          nbb=len_trim(atemp2)
          call getfile(ifail,atemp2,atemp3,1,nbb)
          if(atemp3.eq.'log')then
            replace(irep)%transformtype=1
          else if(atemp3.eq.'none')then
            replace(irep)%transformtype=0
          else
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,200) trim(aline)
200         format('   TRANSFORMTYPE must be "log" or "none" at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MININTERPVAL')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MININTERPVAL'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          replace(irep)%minival=char2real(ifail,2)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,201) trim(aline)
201         format('   Cannot read value for MININTERPVAL at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MAXINTERPVAL')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MAXINTERPVAL'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          replace(irep)%maxival=char2real(ifail,2)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,202) trim(aline)
202         format('   Cannot read value for MAXINTERPVAL at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MININTERPFILE')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MININTERPFILE'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          bfile=cline(left_word(2):)
          nbb=len_trim(bfile)
          call getfile(ifail,bfile,replace(irep)%minifile,1,nbb)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,206) trim(aline)
206         format('   Cannot read filename at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MAXINTERPFILE')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MAXINTERPFILE'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          bfile=cline(left_word(2):)
          nbb=len_trim(bfile)
          call getfile(ifail,bfile,replace(irep)%maxifile,1,nbb)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,206) trim(aline)
            call write_message(increment=1)
          end if

        else if(atemp1.eq.'MINPARVAL')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MINPARVAL'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          replace(irep)%minpval=char2real(ifail,2)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,203) trim(aline)
203         format('   Cannot read value for MINPARVAL at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MAXPARVAL')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='MAXPARVAL'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 100
          end if
          replace(irep)%maxpval=char2real(ifail,2)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,204) trim(aline)
204         format('   Cannot read value for MAXPARVAL at line ',a)
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'END')then
          if(ireplace.eq.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            akey='END'
            write(amessage,60) trim(akey),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          call casetrans(atemp2,'hi')
          if(atemp2.ne.'REPLACE')then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,210) trim(aline)
210         format('   Line ',a,' should be "END REPLACE"')
            call write_message(increment=1)
          end if
          ireplace=0
        else if(atemp1.eq.'REMOVEZONE')then
          if(ireplace.eq.1)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,215) trim(aline)
215         format('   "END REPLACE" expected before "REMOVEZONE" at line ',a)
            call write_message(increment=1)
            go to 100
          end if
          izone=izone+1
          nbb=len_trim(atemp2)
          call getfile(ifail,atemp2,atemp3,1,nbb)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,220) trim(aline)
220         format('   Cannot read zone array name at line ',a)
            call write_message(increment=1)
          else
            call casetrans(atemp3,'hi')
            removezone(izone)=atemp3
          end if
        else if(atemp1.eq.'REMOVEMULT')then
          if(ireplace.eq.1)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,225) trim(aline)
225         format('   "END REPLACE" expected before "REMOVEMULT" at line ',a)
            call write_message(increment=1)
            go to 100
          end if
          imult=imult+1
          nbb=len_trim(atemp2)
          call getfile(ifail,atemp2,atemp3,1,nbb)
          if(ifail.ne.0)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,230) trim(aline)
230         format('   Cannot read multiplier array name at line ',a)
            call write_message(increment=1)
          else
            call casetrans(atemp3,'hi')
!            if(atemp3.eq.'FUNCTION')then
!	      if(imessage.eq.0) call write_initial_message(leadspace='yes')
!              write(amessage,232) trim(aline)
!232           format('   Illegal multiplier array name at line ',a)
!              call write_message(increment=1)
!            end if
            removemult(imult)=atemp3
          end if
        else
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
          write(amessage,240) trim(atemp1),trim(aline)
240       format('   Unrecognised keyword - "',a,'" at line ',a)
          call write_message(increment=1)
        end if
        go to 100

245     continue
        if(imessage.ne.0) go to 9990

! -- Now that the file has been read, some further checks are made.

        do irep=1,numparamreplace
          if(replace(irep)%newparamlay1(1).eq.-9999)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            if(replace(irep)%replacetype.eq.'RCH')then
              write(amessage,269) trim(anum)
269           format('   No NEWPARAMPERS values supplied for replacement set no. ',a)
            else
              write(amessage,270) trim(anum)
270           format('   No NEWPARAMLAYS values supplied for replacement set no. ',a)
            end if
            call write_message(increment=1)
          end if
          if(replace(irep)%transformtype.eq.-9999)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,295) trim(anum)
295         format('   No TRANSFORMTYPE value supplied for replacement set no. ',a)
            call write_message(increment=1)
          end if
          if(replace(irep)%minpval.lt.-1.0e35)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,296) trim(anum)
296         format('   No MINPARVAL value supplied for replacement set no. ',a)
            call write_message(increment=1)
          end if
          if(replace(irep)%maxpval.lt.-1.0e35)then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,297) trim(anum)
297         format('   No MAXPARVAL value supplied for replacement set no. ',a)
            call write_message(increment=1)
          end if
          if((replace(irep)%minival.gt.-1.0e35).and.(replace(irep)%minifile.ne.' '))then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,320) trim(anum)
320         format('   Values are supplied for both MININTERPVAL and MININTERPFILE for ', &
            'replacement set no. ',a,': a value for only one of these must be supplied.')
            call write_message(increment=1)
          end if
          if((replace(irep)%maxival.gt.-1.0e35).and.(replace(irep)%maxifile.ne.' '))then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,321) trim(anum)
321         format('   Values are supplied for both MAXINTERPVAL and MAXINTERPFILE for ', &
            'replacement set no. ',a,': a value for only one of these must be suppled.')
            call write_message(increment=1)
          end if
          if((replace(irep)%minival.lt.-1.0e35).and.(replace(irep)%minifile.eq.' '))then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,322) trim(anum)
322         format('   A value must be supplied for either MININTERPVAL or ',&
            'MININTERPFILE for replacement set no. ',a,'.')
            call write_message(increment=1)
          end if
          if((replace(irep)%maxival.lt.-1.0e35).and.(replace(irep)%maxifile.eq.' '))then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
            call num2char(irep,anum)
            write(amessage,323) trim(anum)
323         format('   A value must be supplied for either MAXINTERPVAL or ',&
            'MAXINTERPFILE for replacement set no. ',a,'.')
            call write_message(increment=1)
          end if
          if(irep.gt.1)then
            do jrep=1,irep-1
              if(replace(jrep)%newparamprefix.eq. &
                 replace(irep)%newparamprefix)then
	         if(imessage.eq.0) call write_initial_message(leadspace='yes')
                   write(amessage,326)
326                format('   There are at least two replacement sets with ', &
                   'the same NEWPARAMPREFIX value.')
                   call write_message(increment=1)
                   go to 329
               end if
            end do
          end if
329       continue
        end do
        if(imessage.ne.0) go to 9990


! -- If necessary the array files are read. First memory is allocated for arrays.

        do irep=1,numparamreplace
          allocate(replace(irep)%miniarray(ncol,nrow), &
                   replace(irep)%maxiarray(ncol,nrow),stat=ierr)
          if(ierr.ne.0)then
            write(amessage,350)
350         format(' Cannot allocate sufficient memory to continue execution.')
            go to 9890
          end if
          if(replace(irep)%minifile.eq.' ')then
            replace(irep)%miniarray=replace(irep)%minival
          else
            iunit=nextunit()
            afile=replace(irep)%minifile
            open(unit=iunit,file=afile,status='old',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,360) trim(afile),trim(repfile)
360           format(' Cannot open real array file ',a,' cited in file ',a)
              go to 9890
            end if
            if(headerspec.eq.'yes')then
              read(iunit,*,err=9050,end=9100)mcol,mrow
              if((mcol.ne.ncol).or.(mrow.ne.nrow))then
                write(amessage,370) trim(afile)
370             format(' Number of rows and columns cited in header of file ',a, &
                ' does not match expectations.')
                go to 9890
              end if
            end if
            do irow=1,nrow
              read(iunit,*,err=9050,end=9100) &
              (replace(irep)%miniarray(icol,irow),icol=1,ncol)
            end do
            close(unit=iunit)
          end if
          if(replace(irep)%maxifile.eq.' ')then
            replace(irep)%maxiarray=replace(irep)%maxival
          else
            iunit=nextunit()
            afile=replace(irep)%maxifile
            open(unit=iunit,file=afile,status='old',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,360) trim(afile),trim(repfile)
              go to 9890
            end if
            if(headerspec.eq.'yes')then
              read(iunit,*,err=9050,end=9100)mcol,mrow
              if((mcol.ne.ncol).or.(mrow.ne.nrow))then
                write(amessage,370) trim(afile)
                go to 9890
              end if
            end if
            do irow=1,nrow
              read(iunit,*,err=9050,end=9100) &
              (replace(irep)%maxiarray(icol,irow),icol=1,ncol)
            end do
            close(unit=iunit)
          end if
          do irow=1,nrow
            do icol=1,ncol
              if(replace(irep)%miniarray(icol,irow).gt. &
                 replace(irep)%maxiarray(icol,irow))then
                   call num2char(irep,anum)
                   write(amessage,390) trim(anum)
390                format(' The minimum interpolation value exceeds the maximum ', &
                   'interpolation value in at least one part of the model domain ', &
                   'for replacement set ',a,': check MININTERPVAL, MAXINTERPVAL, ', &
                   'MININTERPFILE and/or MAXINTERPFILE for this replacement set.')
                   go to 9890
               end if
            end do
          end do
        end do

        if(imessage.ne.0) go to 9990
        ifail=0
        go to 9999


9000    write(amessage,9010) trim(aline),trim(repfile)
9010    format(' Error encountered reading line ',a,' of file ',a)
        go to 9890

9050    write(amessage,9060) trim(afile)
9060    format(' Error reading data from real array file ',a)
        go to 9890
9100    write(amessage,9110) trim(afile)
9110    format(' Unexpected end encountered to real array file ',a)
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9990    ifail=1
9999    close(unit=iunit,iostat=ierr)
	imessage=0
	return

end subroutine read_parameter_replacement_file


