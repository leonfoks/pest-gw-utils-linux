!     Last change:  JD    9 May 2003   11:49 am
program reparray

! -- Program REPARRAY replaces an array in a MODFLOW/MT3D input file by a 
!    user-supplied array.


	use defn
	use inter

	implicit none

	logical				    :: lexist
	integer				    :: ifail,ierr,ncol,nrow,irow,icol
	integer				    :: modunit1,modunit2,n,iline,    &
					       locline,irep,iprn,iunit,idel, &
                                               iout,i96,iicol,idate,iheader,nbb,ifail1
	real, allocatable, dimension(:,:)   :: rarray,darray
	character (len=15)		    :: anum
	character (len=24)		    :: aname
	character (len=40)		    :: aheader,bheader
	character (len=80)		    :: aprompt,atext,btext,ctext
	character (len=80)		    :: modfile1,modfile2
	character (len=3000)		    :: ccline,capline
	type (modelgrid) 		    :: gridspec
	character (len=1)		    :: aa,bb


	aheader='                   1(7F14.0)            '
	bheader='INTERNAL 1.0 (FREE) '
	idel=0

	write(amessage,5)
5	format(' Program REPARRAY replaces an array in a MODFLOW/MT3D input ', &
	'file with a user-supplied array.')
	call write_message(leadspace='yes',endspace='yes')

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9999
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message
	  go to 9999
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9999
	end if

	call readfig(gridspec%specfile)
10	call spec_open(ifail,gridspec)
	if((ifail.ne.0).or.(escset.eq.1)) then
	  write(*,*)
	  go to 9998
	end if
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) then
	  write(*,*)
	  go to 9998
	end if
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(rarray(ncol,nrow),darray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run REPARRAY.')
	  go to 9997
	end if

! -- First the existing real array is read.

55	write(*,*)
	aprompt=' Enter name of real array file: '
	call read_real_array(ifail,aprompt,rarray,pm_header=headerspec, &
        rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9998
	if(escset.eq.1) then
	  escset=0
	  write(*,*)
	  deallocate(rarray,darray,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9997
	  end if
!	  call free_grid_mem(gridspec)
	  go to 10
	end if
	write(*,*)

! -- Next the MODFLOW/MT3D input file is identified.

70	call open_input_file(ifail, &
	' Enter name of MODFLOW/MT3D input file: ',modfile1,modunit1)
	if(ifail.ne.0) go to 9998
	if(escset.ne.0) then
	  escset=0
	  go to 55
	end if
75	write(6,80,advance='no')
80	format(' Is this a MODFLOW or MT3D input file? [f/t]: ')
	read(5,'(a)') bb
	if(bb.eq.' ') go to 75
	if(index(eschar,bb).ne.0)then
	  close(unit=modunit1,iostat=ierr)
	  if(ierr.ne.0) go to 9600
	  write(6,*)
	  go to 70
	end if
	call casetrans(bb,'lo')
	if((bb.ne.'f').and.(bb.ne.'t')) go to 75

90	write(6,100,advance='no')
100	format(' Locate array using text or line number?  [t/l]: ')
	read(5,'(a)') aa
	if(aa.eq.' ') go to 90
	if(index(eschar,aa).ne.0)then
	  write(*,*)
	  go to 75
	end if
	call casetrans(aa,'lo')
	if((aa.ne.'t').and.(aa.ne.'l')) go to 90

	if(aa.eq.'t')then
	  atext=' '
120	  write(6,130,advance='no')
130	  format(' Enter text: ')
	  read(5,'(a)') atext
	  if(atext.eq.' ') go to 120
	  if(index(eschar,atext(1:2)).ne.0)then
	    write(*,*)
	    go to 90
	  end if
          nbb=len_trim(atext)
          call getfile(ifail1,atext,ctext,1,nbb)
          if(ifail1.ne.0) go to 120
          atext=ctext
	  btext=atext
	  call casetrans(atext,'hi')
	else
150	  write(6,160,advance='no')
160	  format(' Enter line number: ')
	  read(5,'(a)') anum
	  if(anum.eq.' ') go to 150
	  anum=adjustl(anum)
	  if(index(eschar,anum(1:2)).ne.0) then
	    write(*,*)
	    go to 90
	  end if
	  call char2num(ifail,anum,locline)
	  if(ifail.ne.0) go to 150
	  if(locline.le.0) go to 150
	end if

	write(*,*)
	aprompt=' Enter name for new MODFLOW/MT3D input file: '
	call open_output_file(ifail,aprompt,modfile2,modunit2)
	if(ifail.ne.0) go to 9998
	if(escset.ne.0) then
	  escset=0
	  write(*,*)
	  go to 90
	end if

! -- Next we rewrite data from the old MODFLOW/MT3D input file to the new one.

	write(*,*)
	idel=1
	irep=0
	if(aa.eq.'t') then
	  do
	    read(modunit1,'(a)',end=500) ccline
	    n=len_trim(ccline)
	    if(n.gt.3000) go to 9350
	    capline=ccline
	    call casetrans(capline,'hi')
	    if(index(capline,trim(atext)).eq.0) then
	      write(modunit2,'(a)',err=9300) ccline(1:n)
	      cycle
	    end if
	    irep=irep+1
	    write(amessage,170)
170	    format(' Replacing array......')
	    call write_message(endspace='yes')
	    call u2drel(ifail,darray,aname,nrow,ncol,1,modunit1,iout,ccline, &
	    i96,iprn,iunit,iicol,modunit1,bb)
	    if((ifail.eq.1).or.(ifail.eq.2)) then
	      write(amessage,240) trim(modfile1)
240	      format(' Improper or disallowed header to selected array in ',  &
              'MODFLOW/MT3D input file ',a,'.')
	      call write_message
	      write(amessage,210)
210	      format(' Array header follows:-')
	      call write_message
	      write(*,'(a)') ccline(1:n)
	      go to 9998
	    end if
	    if(ifail.eq.3)then
	      if(bb.eq.'f')then
	        write(amessage,260) trim(modfile1)
260	        format(' LOCAT value of array header in file ',a,          &
	        ' must not be zero.')
	      else
	        write(amessage,261) trim(modfile1)
261	        format(' IREAD value of array header in file ',a,          &
	        ' must not be zero.')
	      end if
	      go to 9997
	    end if
	    if(bb.eq.'t')then
	      if((iunit.ne.100).and.(iunit.ne.103))then
	        write(amessage,265) trim(modfile1)
265	        format(' IREAD value of array header in MT3D input ',      &
                'file ',a,' must be 100 or 103.')
	        go to 9997
	      end if
	    end if
	    if(ifail.eq.4)then
	      write(amessage,268) trim(modfile1)
268	      format(' Error reading array from file ',a)
	      go to 9997
	    end if
	    if(i96.eq.0)then
	      write(aheader(1:10),'(i10)') iunit
	      write(modunit2,275,err=9300) aheader,ccline(41:n)
275	      format(a40,a)
	    else
	      if(iunit.eq.0) iprn=0
	      if(iicol.le.len(ccline))then
	        write(modunit2,280) trim(bheader),iprn,ccline(iicol:n)
280	        format(1x,a,i8,2x,a)
	      else
	        write(modunit2,281) trim(bheader),iprn
281	        format(1x,a,i8)
	      end if
	    end if
	    do irow=1,nrow
	      write(modunit2,310,err=9300) (rarray(icol,irow),icol=1,ncol)
310	      format(7(1pe14.6))  		!check carriage control
	    end do
	  end do
	else
	  do iline=1,locline-1
	    read(modunit1,'(a)',end=500) ccline
	    n=len_trim(ccline)
	    if(n.gt.3000) go to 9350
	    write(modunit2,'(a)',err=9300) ccline(1:n)
	  end do
	  write(amessage,170)
	  call write_message(endspace='yes')
	  read(modunit1,'(a)',end=500) ccline
	  n=len_trim(ccline)
	  if(n.gt.3000) go to 9350
	  call u2drel(ifail,darray,aname,nrow,ncol,1,modunit1,iout,	&
          ccline,i96,iprn,iunit,iicol,modunit1,bb)
	  if((ifail.eq.1).or.(ifail.eq.2)) then
	    write(amessage,240) trim(modfile1)
	    call write_message
	    write(amessage,210)
	    call write_message
	    write(*,'(a)') ccline(1:n)
	    go to 9998
	  end if
	  if(ifail.eq.3)then
	    if(bb.eq.'f')then
	      write(amessage,260) trim(modfile1)
	    else
	      write(amessage,261) trim(modfile1)
	    end if
	    go to 9997
	  end if
	  if(bb.eq.'t')then
	    if((iunit.ne.100).and.(iunit.ne.103))then
	      write(amessage,265) trim(modfile1)
	      go to 9997
	    end if
	  end if
	  if(ifail.eq.4)then
	    write(amessage,268) trim(modfile1)
	    go to 9997
	  end if
	  if(i96.eq.0)then
	    write(aheader(1:10),'(i10)') iunit
	    write(modunit2,275,err=9300) aheader,ccline(41:n)
	  else
	    if(iunit.eq.0) iprn=0
	    if(iicol.le.len(ccline))then
	      write(modunit2,280) trim(bheader),iprn,ccline(iicol:n)
	    else
	      write(modunit2,281) trim(bheader),iprn
	    end if
	  end if
	  do irow=1,nrow
	    write(modunit2,310,err=9300) (rarray(icol,irow),icol=1,ncol)
	  end do
	  irep=1
	  do
	    read(modunit1,'(a)',end=500) ccline
	    n=len_trim(ccline)
	    if(n.gt.3000) go to 9350
	    write(modunit2,'(a)',err=9300) ccline(1:n)
	  end do
	end if

500	if(irep.eq.0) then
	  if(aa.eq.'l') then
	    call num2char(locline,anum)
	    write(amessage,550) trim(modfile1),trim(anum)
550	    format(' No array replacement: end of file ',a,             &
            ' encountered before line ',a,' was reached.')
	    go to 9997
	  else
	    write(amessage,570) trim(btext),trim(modfile1)
570	    format(' No array replacement: locator text "',a,           &
            '" not found in file ',a,'.')
	    go to 9997
	  end if
	else
	  if(irep.eq.1) then
	    write(amessage,600) trim(modfile1)
600	    format('  - 1 array replaced in MODFLOW/MT3D input file ',a,'.')
	  else
	    call num2char(irep,anum)
	    write(amessage,620) trim(anum),trim(modfile1)
620	    format('  - ',a,' arrays replaced in MODFLOW/MT3D input ',     &
            'file ',a,'.')
	  endif
	  call write_message
	  write(amessage,640) trim(modfile2)
640	  format('  - file ',a,' is new MODFLOW/MT3D input file.')
	  call write_message
	  go to 9999
	end if


9300	write(amessage,9310) trim(modfile2)
9310	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	go to 9997
9350    write(amessage,9360) trim(modfile1)
9360	format(' File ',a,' is too wide: maximum allowed file ',            &
        'width is 3000 characters.')
	go to 9997
9600	write(amessage,9610)
9610	format(' Cannot close file.')
	go to 9997

9997	call write_message(leadspace='yes',endspace='yes')
9998    call close_files
	call free_grid_mem(gridspec)
	deallocate(rarray,darray,stat=ierr)

	if(idel.eq.1)then 
	  inquire(file=modfile2,exist=lexist)
	  if(lexist)call system('del "'//trim(modfile2)//'"')
	end if
	call exit(100)

9999	call close_files
	call free_grid_mem(gridspec)
	deallocate(rarray,darray,stat=ierr)

end program reparray


      SUBROUTINE U2DREL(ifail,A,ANAME,II,JJ,K,IN,IOUT,cntrl,i96,iprn,	&
      locat,icol,modunit1,bb)
!
!
!-----VERSION 1539 22JUNE1993 U2DREL
!     ******************************************************************
!     ROUTINE TO INPUT 2-D REAL DATA MATRICES
!       A IS ARRAY TO INPUT
!       ANAME IS 24 CHARACTER DESCRIPTION OF A
!       II IS NO. OF ROWS
!       JJ IS NO. OF COLS
!       K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --)
!              IF K=0, NO LAYER IS PRINTED
!              IF K<0, CROSS SECTION IS PRINTED)
!       IN IS INPUT UNIT
!       IOUT IS OUTPUT UNIT
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------

	integer ifail,jj,in,i96,iprn,locat,icol,nunopn,iclose,ifree,  &
                istart,istop,n,i,j,ii,k,iout,modunit1,iii
	real a,cnstnt,r


      character*(*) bb						!wnc
      CHARACTER*24 ANAME
      DIMENSION A(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*(*) CNTRL					!wnc
      CHARACTER*16 TEXT
      CHARACTER*80 FNAME
      DATA NUNOPN/99/
!     ------------------------------------------------------------------
!
!1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
!      READ(IN,'(A)') CNTRL					!wnc
!
!2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      ifail=0							!wnc
      i96=0							!wnc
      CALL URWORD(ifail,CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      if(ifail.eq.1) return
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
         i96=1							!wnc
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
         i96=1							!wnc
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         ifail=1						!wnc
         return							!wnc
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         ifail=1						!wnc
         return							!wnc
      ELSE
!
!2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=500) LOCAT,CNSTNT,FMTIN,IPRN
    1    FORMAT(I10,F10.0,A20,I10)
         if(locat.eq.0)then					!wnc
           ifail=3						!wnc
           return						!wnc
         end if							!wnc
      END IF
!
!3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(ifail,CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN)
	 if(ifail.eq.1) return
         IF(LOCAT.NE.0) THEN
            CALL URWORD(ifail,CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
	    if(ifail.eq.1) return
            FMTIN=CNTRL(ISTART:ISTOP)
            IF(ICLOSE.NE.0) THEN
               IF(FMTIN.EQ.'(BINARY)') THEN
                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM='BINARY',    &
                  iostat=iii)
	          if(iii.ne.0)then
	            OPEN(UNIT=LOCAT,FILE=FNAME,FORM='UNFORMATTED')
	          end if
               ELSE
                  OPEN(UNIT=LOCAT,FILE=FNAME)
               END IF
            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)') LOCAT=-LOCAT
            CALL URWORD(ifail,CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
	    if(ifail.eq.1) return
         END IF
      END IF
!
!4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT) 200,50,90
!
!4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
   50 DO 80 I=1,II
      DO 80 J=1,JJ
   80 A(J,I)=CNSTNT
      RETURN
!
!4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
90    CONTINUE

      if((locat.eq.103).and.(bb.eq.'t')) FMTIN='(FREE)'
      DO 100 I=1,II
      IF(FMTIN.EQ.'(FREE)') THEN
         READ(modunit1,*,err=600,end=600) (A(J,I),J=1,JJ)
      ELSE
         READ(modunit1,FMTIN,err=600,end=600) (A(J,I),J=1,JJ)
      END IF
  100 CONTINUE
      GO TO 300
!
!4C-----LOCAT<0; READ UNFORMATTED ARRAY VALUES.
  200 ifail=1
      return
!
!5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
  300 IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
!      ZERO=0.
!      IF(CNSTNT.EQ.ZERO) GO TO 320
!      DO 310 I=1,II
!      DO 310 J=1,JJ
!      A(J,I)=A(J,I)*CNSTNT
!  310 CONTINUE

      RETURN
!
!8------CONTROL RECORD ERROR.
  500 ifail=2
      return
600   ifail=4							!wnc
      return							!wnc
      END





      SUBROUTINE URWORD(ifail,LINE,ICOL,ISTART,ISTOP,NCODE,N,R,IOUT,IN)
!
!
!-----VERSION 1003 05AUG1992 URWORD
!     ******************************************************************
!     ROUTINE TO EXTRACT A WORD FROM A LINE OF TEXT, AND OPTIONALLY
!     CONVERT THE WORD TO A NUMBER.
!        ISTART AND ISTOP WILL BE RETURNED WITH THE STARTING AND
!          ENDING CHARACTER POSITIONS OF THE WORD.
!        THE LAST CHARACTER IN THE LINE IS SET TO BLANK SO THAT IF ANY
!          PROBLEMS OCCUR WITH FINDING A WORD, ISTART AND ISTOP WILL
!          POINT TO THIS BLANK CHARACTER.  THUS, A WORD WILL ALWAYS BE
!          RETURNED UNLESS THERE IS A NUMERIC CONVERSION ERROR.  BE SURE
!          THAT THE LAST CHARACTER IN LINE IS NOT AN IMPORTANT CHARACTER
!          BECAUSE IT WILL ALWAYS BE SET TO BLANK.
!        A WORD STARTS WITH THE FIRST CHARACTER THAT IS NOT A SPACE OR
!          COMMA, AND ENDS WHEN A SUBSEQUENT CHARACTER THAT IS A SPACE
!          OR COMMA.  NOTE THAT THESE PARSING RULES DO NOT TREAT TWO
!          COMMAS SEPARATED BY ONE OR MORE SPACES AS A NULL WORD.
!        FOR A WORD THAT BEGINS WITH "'", THE WORD STARTS WITH THE
!          CHARACTER AFTER THE QUOTE AND ENDS WITH THE CHARACTER
!          PRECEDING A SUBSEQUENT QUOTE.  THUS, A QUOTED WORD CAN
!          INCLUDE SPACES AND COMMAS.  THE QUOTED WORD CANNOT CONTAIN
!          A QUOTE CHARACTER.
!        IF NCODE IS 1, THE WORD IS CONVERTED TO UPPER CASE.
!        IF NCODE IS 2, THE WORD IS CONVERTED TO AN INTEGER.
!        IF NCODE IS 3, THE WORD IS CONVERTED TO A REAL NUMBER.
!        NUMBER CONVERSION ERROR IS WRITTEN TO UNIT IOUT IF IOUT IS
!          POSITIVE; ERROR IS WRITTEN TO DEFAULT OUTPUT IF IOUT IS 0;
!          NO ERROR MESSAGE IS WRITTEN IF IOUT IS NEGATIVE.
!     ******************************************************************
!
!        SPECIFICATIONS:
!     ------------------------------------------------------------------

	integer icol,istart,istop,ncode,n,linlen,i,j,idiff,k,l,iout,	&
                in,ifail
	real r


      CHARACTER*(*) LINE
      CHARACTER*20 RW,STRING
!     ------------------------------------------------------------------

	ifail=0

!
!1------Set last char in LINE to blank and set ISTART and ISTOP to point
!1------to this blank as a default situation when no word is found.  If
!1------starting location in LINE is out of bounds, do not look for a
!1------word.
      LINLEN=LEN(LINE)
      LINE(LINLEN:LINLEN)=' '
      ISTART=LINLEN
      ISTOP=LINLEN
      LINLEN=LINLEN-1
      IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
!
!2------Find start of word, which is indicated by first character that
!2------is not a blank and not a comma.
      DO 10 I=ICOL,LINLEN
      IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.',') GO TO 20
10    CONTINUE
      ICOL=LINLEN+1
      GO TO 100
!
!3------Found start of word.  Look for end.
!3A-----When word is quoted, only a quote can terminate it.
20    IF(LINE(I:I).EQ.'''') THEN
         I=I+1
         IF(I.LE.LINLEN) THEN
            DO 25 J=I,LINLEN
            IF(LINE(J:J).EQ.'''') GO TO 40
25          CONTINUE
         END IF
!
!3B-----When word is not quoted, space or comma will terminate.
      ELSE
         DO 30 J=I,LINLEN
         IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.',') GO TO 40
30       CONTINUE
      END IF
!
!3C-----End of line without finding end of word; set end of word to
!3C-----end of line.
      J=LINLEN+1
!
!4------Found end of word; set J to point to last character in WORD and
!-------set ICOL to point to location for scanning for another word.
40    ICOL=J+1
      J=J-1
      IF(J.LT.I) GO TO 100
      ISTART=I
      ISTOP=J
!
!5------Convert word to upper case and RETURN if NCODE is 1.
      IF(NCODE.EQ.1) THEN
         IDIFF=ICHAR('a')-ICHAR('A')
         DO 50 K=ISTART,ISTOP
            IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z')			&
                   LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50       CONTINUE
         RETURN
      END IF
!
!6------Convert word to a number if requested.
100   IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
         RW=' '
         L=20-ISTOP+ISTART
         IF(L.LT.1) GO TO 200
         RW(L:20)=LINE(ISTART:ISTOP)
         IF(NCODE.EQ.2) READ(RW,'(I20)',ERR=200) N
         IF(NCODE.EQ.3) READ(RW,'(F20.0)',ERR=200) R
      END IF
      RETURN
!

!7------Number conversion error.

200   ifail=1
      return


      END
 