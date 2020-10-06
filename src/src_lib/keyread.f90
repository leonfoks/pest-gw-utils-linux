!*****************************************************************************
! functions comprising the generic function KEY_READ
!*****************************************************************************

! -- Function key_read reads and checks a number from the keyboard.

! -- Arguments are as follows:-
!       value:   value of the number read
!       returns  0 unless an error condition arises

! -- Revision history:-
!       June-November, 1995: version 1.

integer function int_key_read(value)

	use defn

	integer, intent(out)    :: value
	integer                 :: ifail
	character (len=50)      :: atemp

	int_key_read=0
	read(5,'(a)') atemp
	if(atemp.ne.' ') atemp=adjustl(atemp)   ! lf90 bug
	if(atemp.eq.' ') then
	  int_key_read=-1
	else if(index(eschar,atemp(1:2)).ne.0) then
	  escset=1
	else
	  call a2i(ifail,atemp,value)
	  if(ifail.ne.0) int_key_read=1
	end if
	return

end function int_key_read


integer function real_key_read(value)

	use defn

	real, intent(out)       :: value
	integer                 :: ifail
	character (len=50)      :: atemp

	real_key_read=0
	read(5,'(a)') atemp
	if(atemp.ne.' ') atemp=adjustl(atemp)   ! lf90 bug
	if(atemp.eq.' ') then
	  real_key_read=-1
	else if(index(eschar,atemp(1:2)).ne.0) then
	  escset=1
	else
	  call a2r(ifail,atemp,value)
	  if(ifail.ne.0) real_key_read=1
	end if
	return

end function real_key_read


integer function double_key_read(value)

	use defn

	double precision, intent(out)   :: value
	integer                         :: ifail
	character (len=50)              :: atemp

	double_key_read=0
	read(5,'(a)') atemp
	if(atemp.ne.' ') atemp=adjustl(atemp)   ! lf90 bug
	if(atemp.eq.' ') then
	  double_key_read=-1
	else if(index(eschar,atemp(1:2)).ne.0) then
	  escset=1
	else
	  call a2d(ifail,atemp,value)
	  if(ifail.ne.0) double_key_read=1
	end if
	return

end function double_key_read

 