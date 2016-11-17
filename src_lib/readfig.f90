!     Last change:  JD   21 Dec 2000    2:33 pm
subroutine readfig(specfile,coordfile,sampfile,pumpfile,pilotfile)

! -- Subroutine readfig reads a filename file.

! -- Arguments are as follows:-
!       specfile:  name of grid specification file
!       coordfile: name of bore coordinates file
!       sampfile:  name of bore sample file
!       pumpfile:  name of bore pumping file


	use defn
	use inter

	character (len=*),intent(out)           :: specfile
	character (len=*),intent(out),optional  :: coordfile,sampfile,pumpfile,pilotfile
	integer                                 :: iunit,ierr,iequals,nb,ifail
        character (len=100)                     :: filename

	specfile=' '
	if(present(coordfile)) coordfile=' '
	if(present(sampfile)) sampfile=' '
	if(present(pumpfile)) pumpfile=' '
        if(present(pilotfile)) pilotfile=' '

	iunit=nextunit()
	open(unit=iunit,file='files.fig',status='old',err=200)
	do
	  read(iunit,'(a)',err=100,end=200) cline
          cline=adjustl(cline)
	  iequals=index(cline,'=')
	  if(iequals.le.1) cycle
          nb=len_trim(cline)
          if(nb.eq.iequals) return
          call getfile(ifail,cline,filename,iequals+1,nb)
          if(ifail.ne.0) cycle
	  call casetrans(cline(1:iequals),'lo')
	  if(cline(1:10).eq.'grid_speci') then
	    specfile=filename
	  else if(cline(1:10).eq.'bore_coord') then
	    if(present(coordfile)) then
	      coordfile=filename
	    end if
	  else if(cline(1:10).eq.'bore_sampl') then
	    if(present(sampfile)) then
	      sampfile=filename
	    end if
	  else if(cline(1:10).eq.'bore_pumpi') then
	    if(present(pumpfile)) then
	      pumpfile=filename
	    end if
          else if(cline(1:12).eq.'pilot_points')then
            if(present(pilotfile))then
              pilotfile=filename
            end if
	  end if
100       continue
	end do

200     close(unit=iunit,iostat=ierr)
	return

end subroutine readfig

 