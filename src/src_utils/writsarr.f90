!     Last change:  JD   28 Dec 2000    9:20 pm
subroutine write_surf_array(ifail,aprompt,array,thresh,gridspec)

!  -- Subroutine write_surf_array writes a real array in SURFER GRID format.
!     It is assumed that the finite-difference grid has uniform DELR and DELC
!     elements.

! -- Arguments are as follows:-
!       ifail:      returned as zero unless error condition arises
!       aprompt:    prompt text for name of output file
!       array:      real array to write in SURFER grid format
!       thresh:     blanking threshold
!       gridspec:   defined type holding grid specifications

	use defn
	use inter

	integer, intent(out)                    :: ifail
	character (len=*), intent(in)           :: aprompt
	real, intent(inout),dimension(:,:)      :: array
	real, intent(in)                        :: thresh
	type (modelgrid), intent(in)            :: gridspec
	integer                                 :: nrow,ncol,i,k,iunit,nbb,ifail1
	character (len=120)                     :: afile,bfile
	character (len=1)                       :: aa
	logical                                 :: lexist
	real                                    :: emin,emax,nmin,nmax,&
						   zmin,zmax,rtemp,delta

	ifail=0

5       imessage=0
10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
	read(5,'(a)')afile
	if(afile.eq.' ') go to 10
	afile=adjustl(afile)
	if(index(eschar,afile(1:2)).ne.0)then
	  escset=1
	  return
	end if
        nbb=len_trim(afile)
        call getfile(ifail1,afile,bfile,1,nbb)
        if(ifail1.ne.0) go to 10
        afile=bfile
!	inquire(file=afile,exist=lexist)
!	if(lexist) then
!55        write(6,56,advance='no')
!56        format(' File already exists: overwrite it?  [y/n]: ')
!	  read(5,'(a)') aa
!	  if(aa.eq.' ') go to 55
!	  call casetrans(aa,'lo')
!	  if((index(eschar,aa).ne.0).or.(aa.eq.'n')) then
!	    write(6,*)
!	    go to 10
!	  end if
!	  if(aa.ne.'y') go to 55
!	end if

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	rtemp=gridspec%rotation
	delta=spacing(rtemp)
	emin=gridspec%delr(1)/2.0d0
	emax=emin+gridspec%delr(1)*(ncol-1)
	nmax=-gridspec%delc(1)/2.0d0
	nmin=nmax-gridspec%delc(1)*(nrow-1)
	if((rtemp-delta.le.0.0).and.(rtemp+delta.ge.0.0))then
	  emin=emin+gridspec%east_corner
	  emax=emax+gridspec%east_corner
	  nmin=nmin+gridspec%north_corner
	  nmax=nmax+gridspec%north_corner
	end if
	k=count(abs(array).lt.thresh)
	if(k.eq.0)then
	  zmin=0.0
	  zmax=0.0
	else
	  zmin=minval(array,mask=(abs(array).lt.thresh))
	  zmax=maxval(array,mask=(abs(array).lt.thresh))
	end if
	where (abs(array).ge.thresh) array=1.70141e38

	iunit=nextunit()
	open(unit=iunit,file=afile,err=100)
	write(iunit,'(a)',err=700)'DSAA'
	write(iunit,'(1x,i5,1x,i5)',err=700) ncol,nrow
	write(iunit,40,err=700) emin,emax
40      format(1x,1pe18.10,2x,1pe18.10)
	write(iunit,40,err=700) nmin,nmax
	write(iunit,40,err=700) zmin,zmax
	do i=nrow,1,-1
	  write(iunit,65,err=700) array(:,i)
65        format(7(1pe14.6))
	  write(iunit,*,err=700)
	end do
	write(amessage,70) trim(afile)
70      format('  - SURFER "grid" file ',a,' written ok.')
	call write_message
	if(k.eq.0)then
	  write(amessage,80)
80        format('  - warning: no abs(array elements) are below blanking threshold.')
	  call write_message
	end if
	go to 999

100     write(amessage,110) trim(afile)
110     format('cannot open file ',a,' for output.')
	call write_message(error='yes',leadspace='yes')
	go to 995
700     write(amessage,710) trim(afile)
710     format(' Cannot write data to file ',a, &
	': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 995

995     ifail=1
999     close(unit=iunit,iostat=i)
	return

end subroutine write_surf_array
 