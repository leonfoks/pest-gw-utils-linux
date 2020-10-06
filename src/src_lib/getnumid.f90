subroutine get_num_ids(ifail,iunit,afile,numid,maxsamp,ignore_x)

! -- Subroutine get_num_ids reads any file in which the first column is
!    comprised of borehole identifiers. It establishes the number of different
!    identifiers present in the file and the maximum number of occurrences of
!    any identifier.

!    It is assumed that all incidences of a particular id are juxtaposed as
!    in a bore sample file. This subroutine only detects changes in the
!    identifier between lines.

!    If the file is a bore sample file, x-affected lines can be ignored if 
!    desired.
   

! -- Arguments are as follows:-
!	ifail:    returned as non-zero if error encountered in reading file
!	iunit:	  the unit number from which to read the file
!	afile:    the name of the file being read
!	numid:    number of separate identifiers (see above asumption)
!	maxsamp:  maximum number of juxtaposed occurrences of any particular id
!	ignore_x: 'yes' if x-afected lines are ignored, 'no' otherwise

! -- Revision history:-
!	April, 1997: version 1.

	use defn
	use inter

	implicit none

	integer, intent(out)                     :: ifail
	integer, intent(in)                      :: iunit
	character (len=*), intent(in)            :: afile
	integer, intent(out)                     :: numid,maxsamp
	character (len=*), optional, intent(in)  :: ignore_x

	integer                        :: iline,isamp,ix,ifail1
	character (len=10)             :: abore,oldbore
	character (len=15)             :: aline
	character (len=3)              :: atemp

	ix=0
	if(present(ignore_x))then
	  atemp=ignore_x
	  call casetrans(atemp,'lo')
	  if(atemp.eq.'yes') then
	    ix=1
	  else if(atemp.eq.'no') then
	    ix=0
	  else
	    call sub_error('GET_NUM_IDS')
	  end if
	end if

	maxsamp=0
	isamp=0
	iline=0
	ifail=0
	oldbore=' '
	numid=0

	read_line: do
	  iline=iline+1
	  read(iunit,'(a)',err=9500,end=9000) cline
	  if(cline.eq.' ') cycle
	  if(ix.eq.1) then
	    call linesplit(ifail1,5)
	    if(ifail1.eq.0) then
	      call casetrans(cline(left_word(5):right_word(5)),'lo')
	      if(cline(left_word(5):right_word(5)).eq.'x') cycle read_line
	    end if
	  else
	    call linesplit(ifail1,1)
	  end if
	  if(right_word(1)-left_word(1).gt.9) then
            call num2char(iline,aline)
            write(amessage,30) trim(aline),trim(afile)
30          format('bore identifier greater than 10 characters in length ',&
            'on line ',a,' of file ',a)
            call write_message(error='yes',leadspace='yes')
            go to 9990
          end if
          abore=cline(left_word(1):right_word(1))
          call casetrans(abore,'hi')
	  if(abore.eq.oldbore) then
	    isamp=isamp+1
	  else
	    if(isamp.gt.maxsamp) maxsamp=isamp
	    numid=numid+1
	    isamp=1
	    oldbore=abore
	  end if
	end do read_line

9000	if(isamp.gt.maxsamp) maxsamp=isamp
	rewind(unit=iunit)
	go to 9999

9500	call num2char(iline,aline)
	write(amessage,9510) trim(aline),trim(afile)
9510	format('cannot read line ',a,' of file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 9990

9990	ifail=1
9999	return

end subroutine get_num_ids
 