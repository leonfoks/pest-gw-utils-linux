subroutine get_ids_and_interval(ifail,iunit,afile,nid,aid,ndays1,nsecs1, &
                                ndays2,nsecs2,ignore_x)

! -- Subroutine reads boreid's from a bore sample file. Normally
!    subroutine get_num_ids will have been run first in order to determine
!    how many different id's feature in the list.
!    This subroutine also obtains the time of the first and last sample
!    pertaining to any bore.

!    x-affected lines can be ignored if desired.

! -- Arguments are as follows:-
!       ifail:      returned as non-zero if error condition arises
!	iunit:      unit number from which information is read
!	afile:      file from which information is read
!	nid:        number of different id's cited in file
!       aid:        an array containing the identifiers
!       ndays1:     number of days since 1/1/1970 until 1st sample date
!       ndays2:     number of days since 1/1/1970 until 2nd sample date
!       nsecs1:     number of seconds since midnight until 1st sample time
!       nsecs2:     number of seconds since midnight until 2nd sample time
!       ignore_x: 'yes' if x-afected lines are ignored, 'no' otherwise


! -- Revision history:-
!	April, 1997: version 1.

	use defn
	use inter

	implicit none

	integer, intent(out)             :: ifail
	integer, intent(in)              :: iunit
	character (len=*), intent(in)    :: afile
	integer, intent(in)              :: nid
	character (len=*), intent(out)   :: aid(nid)
        integer, intent(out)             :: ndays1(nid),nsecs1(nid), &
                                            ndays2(nid),nsecs2(nid)
        character (len=*), intent(in), optional  :: ignore_x


	integer                          :: iid,iline,ndays,nsecs,cols, &
                                            odays,osecs,ifail1,ix
	double precision                 :: value
	character (len=10)               :: abore,aboreold
	character (len=15)               :: aline
	character (len=3)                :: atemp


        ix=0
        if(present(ignore_x))then
          atemp=ignore_x
          call casetrans(atemp,'lo')
          if(atemp.eq.'yes') then
            ix=1
          else if(atemp.eq.'no') then
            ix=0
          else
            call sub_error('GET_IDS_AND_INTERVAL')
          end if
        end if


	ifail=0
	iline=0
	aboreold=' '
	iid=0

	read_line: do
	  iline=iline+1
	  read(iunit,'(a)',end=9000,err=9500) cline
          cols=5
          call linesplit(ifail1,5)
          if(ifail1.lt.0) cycle
	  if(ifail1.eq.0) then
	    call casetrans(cline(left_word(5):right_word(5)),'lo')
	    if(cline(left_word(5):right_word(5)).eq.'x') cycle read_line
          else
            cols=4
            call linesplit(ifail1,4)
            if(ifail1.ne.0)then
              call num2char(iline,aline)
              write(amessage,100) trim(aline),trim(afile)
100           format('insufficient items on line ',a,' of bore sample file ',a)
              call write_message(error='yes',leadspace='yes')
              go to 9990
            end if
          end if

          if(right_word(1)-left_word(1).gt.9) then
            call num2char(iline,aline)
            write(amessage,130) trim(aline),trim(afile)
130         format('bore identifier greater than 10 characters in length ',&
            'on line ',a,' of bore sample file ',a)
            call write_message(error='yes',leadspace='yes')
            go to 9990
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
          call read_rest_of_sample_line(ifail1,cols,ndays, &
          nsecs,value,iline,afile)
          if(ifail1.ne.0) go to 9990
          if(abore.ne.aboreold) then
            iid=iid+1
	    aid(iid)=adjustl(abore)
	    aboreold=abore
	    if(iid.ne.1) then
	      ndays2(iid-1)=odays
	      nsecs2(iid-1)=osecs
	    end if
	    ndays1(iid)=ndays
	    nsecs1(iid)=nsecs
          end if
	  odays=ndays
	  osecs=nsecs
	end do read_line

9000	iid=iid+1
	if(iid.ne.1) then
	  ndays2(iid-1)=odays
	  nsecs2(iid-1)=osecs
	end if
	rewind(unit=iunit)
	go to 9999

9500	write(amessage,9510) trim(afile)
9510	format('cannot re-read file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 9990

9990	ifail=1
9999	return

end subroutine get_ids_and_interval
 