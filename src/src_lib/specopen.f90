!     Last change:  JD   26 Dec 2000    1:12 pm
subroutine spec_open(ifail,gridspec)

! -- Subroutine spec_open opens a grid specification file.

! -- Arguments are as follows:-
!       ifail:     returned as zero unless error condition arises.
!       gridspec:  defined variable containing grid specifications.

! -- Revision history:-
!       June-November, 1995: version 1.

	use defn
	use inter

	type (modelgrid)        :: gridspec
	integer, intent(out)    :: ifail
	integer                 :: ierr,jstep,nbb
	character (len=120)      :: atemp,aatemp,bbtemp

	ifail=0
	atemp=gridspec%specfile
	imessage=0
5       jstep=0
	do
	  jstep=jstep+1
	  if(jstep.ge.5) then
	    ifail=1
	    return
	  end if
10        if(atemp.ne.' ') then
	    write(6,20,advance='no') trim(atemp)
20          format(' Enter name of grid specification file [',a,']: ')
	  else
	    write(6,30,advance='no')
30          format(' Enter name of grid specification file: ')
	  endif
	  read(5,'(a)',err=10) aatemp
          nbb=len_trim(aatemp)
          if(aatemp.ne.' ')then
            call getfile(ifail,aatemp,bbtemp,1,nbb)
            if(ifail.ne.0) go to 10
          else
            bbtemp=' '
          end if
          gridspec%specfile=bbtemp
	  if(index(eschar,gridspec%specfile(1:2)).ne.0) then
	    escset=1
	    return
	  end if
	  if(gridspec%specfile.ne.' ') exit
	  if(atemp.ne.' ')then
	    gridspec%specfile=atemp
	    exit
	  end if
	end do

	gridspec%specunit=nextunit()
	open(unit=gridspec%specunit,file=gridspec%specfile,status='old',&
	iostat=ierr)
	if(ierr.ne.0) then
	  if(imessage.eq.5) then
	    ifail=1
	    return
	  end if
	  write(amessage,50) trim(gridspec%specfile)
50        format('cannot open file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  if(gridspec%specfile.eq.atemp) atemp=' '
	  go to 5
	end if
	ifail=0
	return

end subroutine spec_open
 