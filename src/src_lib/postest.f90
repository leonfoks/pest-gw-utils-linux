!*****************************************************************************
! functions comprising the generic function POS_TEST
!*****************************************************************************

! -- Function pos_test tests whether a number is positive.

! -- Arguments are as follows:-
!       value:    the number to be tested
!       string:   part of the displayed error message if string is not +ve
!       returns   0 if positive; 1 otherwise

! -- Revision history:-
!       June-November, 1995: version 1.

integer function pos_i_test(value,string)

	use defn

	integer, intent(in)             :: value
	character (len=*), intent(in)   :: string

	pos_i_test=0
	if(value.le.0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must be positive  - try again.')
	  call write_message(0,6,'no','no','no')
	  pos_i_test=1
	end if
	return

end function pos_i_test


integer function pos_r_test(value,string)

	use defn

	real, intent(in)                :: value
	character (len=*), intent(in)   :: string

	pos_r_test=0
	if(value.le.0.0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must be positive  - try again.')
	  call write_message(0,6,'no','no','no')
	  pos_r_test=1
	end if
	return

end function pos_r_test


integer function pos_d_test(value,string)

	use defn

	double precision, intent(in)    :: value
	character (len=*), intent(in)   :: string

	pos_d_test=0
	if(value.le.0.0d0) then
	  write(amessage,10) trim(string)
10        format(' Error: ',a,' must be positive  - try again.')
	  call write_message(0,6,'no','no','no')
	  pos_d_test=1
	end if
	return

end function pos_d_test
 