!*****************************************************************************
! functions comprising the generic function NNEG_TEST
!*****************************************************************************

! -- Function nnegtest tests that a number is not negative.

! -- Arguments are as follows:-
!       value:   the number to be tested
!       string:  part of displayed error message if number is -ve.
!       returns  0 if not negative; 1 otherwise

! -- Revision history:-
!       June-November, 1995: version 1.

integer function nneg_i_test(value,string)

	use defn

	integer, intent(in)             :: value
	character (len=*), intent(in)   :: string

	nneg_i_test=0
	if(value.lt.0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must not be negative  - try again.')
	  call write_message(0,6,'no','no','no')
	  nneg_i_test=1
	end if
	return

end function nneg_i_test


integer function nneg_r_test(value,string)

	use defn

	real, intent(in)                :: value
	character (len=*), intent(in)   :: string

	nneg_r_test=0
	if(value.lt.0.0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must not be negative  - try again.')
	  call write_message(0,6,'no','no','no')
	  nneg_r_test=1
	end if
	return

end function nneg_r_test


integer function nneg_d_test(value,string)

	use defn

	double precision, intent(in)    :: value
	character (len=*), intent(in)   :: string

	nneg_d_test=0
	if(value.lt.0.0d0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must not be negative  - try again.')
	  call write_message(0,6,'no','no','no')
	  nneg_d_test=1
	end if
	return

end function nneg_d_test
 