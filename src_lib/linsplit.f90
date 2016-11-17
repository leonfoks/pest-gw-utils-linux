!*****************************************************************************
! subprograms for reading and parsing data (mainly from files) ------->
!*****************************************************************************

subroutine linesplit(ifail,num)

! -- subroutine linesplit splits a line into whitespace-delimited words

! -- Arguments are as follows:-
!       ifail:   returned as -1 if line is blank
!                returned as  1 if less than num segments
!       num:     number of words to be extracted

!    Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	integer, intent(out)            :: ifail
	integer, intent(in)             :: num
	integer                         :: nblc,j,i,nw
	character (len=3)               :: aspace

	ifail=0; nw=0; j=0
	aspace=' ,'//achar(9)   
	if(num.gt.NUM_WORD_DIM) call sub_error('LINESPLIT')
	nblc=len_trim(cline)
	if(nblc.eq.0) then
	  ifail=-1
	  return
	end if

5       if(nw.eq.num) return
	do i=j+1,nblc
	  if(index(aspace,cline(i:i)).eq.0) go to 20
	end do
	ifail=1
	return
20      nw=nw+1
	left_word(nw)=i
	do i=left_word(nw)+1,nblc
	  if(index(aspace,cline(i:i)).ne.0) go to 40
	end do
	right_word(nw)=nblc
	if(nw.lt.num) ifail=1
	return
40      right_word(nw)=i-1
	j=right_word(nw)
	go to 5

end subroutine linesplit


integer function char2int(ifail,num)

! -- Function char2int extracts an integer from a word demarcated by subroutine
!    linesplit.

! -- Arguments are as follows:-
!       ifail:    returned as zero unless an error condition arises
!       num:      the number of the word previously extracted by linesplit
!       returns   value of integer read from word

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn

	integer, intent(in)             :: num
	integer, intent(out)            :: ifail
	character (len=8)               :: afmt

	if(num.gt.NUM_WORD_DIM) call sub_error('CHAR2INT')
	if((right_word(num).lt.left_word(num)).or. &
	  (left_word(num).le.0)) call sub_error('CHAR2INT')

	ifail=0
	afmt='(i   )'
	write(afmt(3:5),'(i3)') right_word(num)-left_word(num)+1
	read(cline(left_word(num):right_word(num)),afmt,err=100) char2int
	return

100     ifail=1
	return

end function char2int


real function char2real(ifail,num)

! -- Function char2real extracts a real number from a word demarcated by
!    subroutine linesplit.

! -- Arguments are as follows:-
!       ifail:    returned as zero unless an error condition arises
!       num:      the number of the word previously extracted by linesplit
!       returns   value of real number read from word

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn

	integer, intent(in)             :: num
	integer, intent(out)            :: ifail
	integer                         :: ierr
	character (len=10)              :: afmt

	if(num.gt.NUM_WORD_DIM) call sub_error('CHAR2REAL')
	if((right_word(num).lt.left_word(num)).or. &
	  (left_word(num).le.0)) call sub_error('CHAR2REAL')

	ifail=0
	afmt='(f   .0)'
	write(afmt(3:5),'(i3)') right_word(num)-left_word(num)+1
	read(cline(left_word(num):right_word(num)),afmt, iostat=ierr) char2real
	if(ierr.ne.0) go to 110
	return

110     ifail=1
	return

end function char2real


double precision function char2double(ifail,num)

! -- Function char2double extracts a double precision number from a word
!    demarcated by subroutine linesplit.

! -- Arguments are as follows:-
!       ifail:    returned as zero unless an error condition arises
!       num:      the number of the word previously extracted by linesplit
!       returns   value of double precision number read from word

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn

	integer, intent(in)             :: num
	integer, intent(out)            :: ifail
	integer                         :: ierr
	character (len=10)              :: afmt

	if(num.gt.NUM_WORD_DIM) call sub_error('CHAR2DOUBLE')
	if((right_word(num).lt.left_word(num)).or. &
	  (left_word(num).le.0)) call sub_error('CHAR2DOUBLE')

	ifail=0
	afmt='(f   .0)'
	write(afmt(3:5),'(i3)') right_word(num)-left_word(num)+1
	read(cline(left_word(num):right_word(num)),afmt, iostat=ierr) char2double
	if(ierr.ne.0) go to 110
	return

110     ifail=1
	return

end function char2double
 