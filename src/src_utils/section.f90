!     Last change:  JD   14 Feb 2001    9:26 am
program section

! -- Program SECTION interpolates multi-layered real array data to a 
!    user-defined transect through the finite-difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer					:: ifail,ierr,outunit,i, &
						   numarray,ncol,nrow,icol,irow
	integer					:: j,ilast,numact,icellno,&
						   jcellno,numhead,itemp1, &
						   itemp2,idate,iheader
	real					:: delta,thresh,gt_thresh, &
						   fac1,fac2,fac3,fac4
	real, dimension(:,:,:), allocatable, target	:: realarray
	real, dimension(:,:), pointer		:: parray
	real, dimension(:), allocatable		:: bhead
	double precision			:: e1,n1,e2,n2,distline, &
  						   cosline,sinline,distance, &
						   east,north
	character (len=80), dimension(:),allocatable	:: realfile
	character (len=80)			:: aprompt,outfile,atempfile,bfile
	character (len=20)			:: atemp
	character (len=9), dimension(:),allocatable	:: ahead
	character (len=12)			:: aeast,anorth
	character (len=9)			:: adistance
	character (len=1)			:: imethod


	write(amessage,5)
5	format(' Program SECTION interpolates multi-layered real array data ',&
        'to a user-defined transect through the finite-difference grid.')
	call write_message(leadspace='yes',endspace='yes')

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
	  call write_message
	  go to 9900
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

	nullify (parray)
	call readfig(gridspec%specfile)
10	call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow
50	write(6,*)
	write(6,60)
60	format(' Enter end points of transect ----->')
70	write(6,80,advance='no')
80	format('    east  coordinate of transect  start: ')
	if(key_read(e1).ne.0) go to 70
	if(escset.eq.1) then
	  write(6,*)
	  escset=0
	  go to 10
	end if
90	write(6,100,advance='no')
100	format('    north coordinate of transect  start: ')
	if(key_read(n1).ne.0) go to 90
	if(escset.eq.1)then
	  escset=0
	  go to 50
	end if
110	write(6,120,advance='no')
120	format('    east  coordinate of transect finish: ')
	if(key_read(e2).ne.0) go to 110
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 90
	end if
130	write(6,140,advance='no')
140	format('    north coordinate of transect finish: ')
	if(key_read(n2).ne.0) go to 130
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
	  go to 110
	end if
	distline=sqrt((e2-e1)*(e2-e1)+(n2-n1)*(n2-n1))
	if(distline.le.0.0d0) then
	  write(amessage,150)
150	  format('transect line has zero length - try again.')
	  call write_message(error='yes',leadspace='yes')
	  go to 50
	end if
	cosline=(e2-e1)/distline
	sinline=(n2-n1)/distline
160	write(6,170,advance='no')
170	format(' Enter sampling interval along transect line: ')
	if(key_read(delta).ne.0) go to 160
	if(escset.ne.0) then
	  write(6,*)
	  escset=0
	  go to 130
	end if
	if(pos_test(delta,'sampling interval').ne.0) go to 160
	if(delta.ge.distline)then
	  write(amessage,220)
220	  format('sampling interval exceeds line length.')
	  call write_message(error='yes')
	  go to 160
	end if

230	write(6,*)
240	write(6,250,advance='no')
250	format(' How many arrays to be interpolated to section line? ')
	if(key_read(numarray).ne.0) go to 240
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 160
	end if
	if(pos_test(numarray,'number of arrays').ne.0) go to 240

	allocate(realarray(0:ncol+1,0:nrow+1,numarray),realfile(numarray), &
        bhead(numarray),ahead(numarray),stat=ierr)
	if(ierr.ne.0) go to 9000

290	write(6,300)
300	format(/,' Enter the names of the files containing real arrays ----->')
	i=0
310	i=i+1
	write(6,*)
	call num2char(i,atemp)
	write(aprompt,320) trim(atemp)
320	format('   real array number ',a,':')
	parray=>realarray(1:ncol,1:nrow,i)
	call read_real_array(ifail,aprompt,parray,pm_header=headerspec, &
        rows=gridspec%nrow,columns=gridspec%ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  if(i.eq.1) then
	    deallocate(realarray,realfile,bhead,ahead,stat=ierr)
	    nullify (parray)
	    go to 230
	  else if (i.eq.2) then
	    go to 290
	  else
	    i=i-2
	  end if
	else
	  realfile(i)=aprompt
	end if
	if(i.eq.numarray) go to 400
	go to 310
400	continue
	nullify (parray)

269	write(6,*)
270	write(6,280,advance='no')
280	format(' Enter blanking threshold value: ')
	if(key_read(thresh).ne.0) go to 270
	if(escset.eq.1) then
	  escset=0
	  i=numarray-1
	  go to 310
	endif
	if(pos_test(thresh,'blanking threshold').ne.0) go to 270
	gt_thresh=thresh+2*spacing(thresh)
	realarray(0,:,:)=gt_thresh
	realarray(ncol+1,:,:)=gt_thresh
	realarray(:,0,:)=gt_thresh
	realarray(:,nrow+1,:)=gt_thresh
385	write(6,390,advance='no')
390	format(' Interpolate to full grid or outer active cell centres? [f/c] ')
	read(5,'(a)') imethod
	if(imethod.eq.' ') go to 385
	if(index(eschar,imethod).ne.0) go to 269
	call casetrans(imethod,'lo')
	if((imethod.ne.'f').and.(imethod.ne.'c')) go to 385

	write(6,*)
	aprompt=' Enter name for output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 385
	end if
	if(ifail.ne.0) go to 9900	

	write(6,200)
200	format(/,' Working.....',/)

! -- The following code is a little unweildly because it incorporates an LF90
!    bug workaround.

	do j=1,numarray
	  atempfile=realfile(j)
	  itemp1=len_trim(atempfile)
	  itemp2=max(itemp1-9,1)
          atempfile=atempfile(itemp2:itemp1)
          call addquote(atempfile,bfile)
	  if(j.eq.numarray)then
	    realfile(j)=bfile
	  else
	    realfile(j)=trim(bfile)//','
	  end if
	end do
	write(outunit,420) (realfile(j),j=1,numarray)
420	format(1x,'distance,',1x,'easting,    ',1x,'northing,   ',1x,100(a13,1x))

	ilast=0
	numact=0
        i=-1
        travel: do
          i=i+1
	  distance=i*delta
	  east=e1+distance*cosline
	  north=n1+distance*sinline
	  if(distance.gt.distline)then
	    if(ilast.eq.0) then
	      east=e2
	      north=n2
	      distance=distline
	    else
	      exit travel
	    end if
	  end if
	  if(abs((distance-distline)/delta).lt.0.01) ilast=1
	  call factor(gridspec,east,north,fac1,fac2,fac3,fac4,icellno,jcellno)
	  if(icellno.eq.-999) cycle travel
	  numhead=0
	  do j=1,numarray
	    if(imethod.eq.'c')then
	      call point_interp(ncol,nrow,thresh,fac1,fac2,fac3,fac4,icellno,&
 	      jcellno,bhead(j),realarray(:,:,j),'inside')
	    else
	      call point_interp(ncol,nrow,thresh,fac1,fac2,fac3,fac4,icellno,&
 	      jcellno,bhead(j),realarray(:,:,j))
	    end if
	    if(bhead(j).gt.1.0e37) then
	      ahead(j)=' '
	    else
	      call num2char(bhead(j),ahead(j),8)
	      numhead=numhead+1
	    end if
	  end do
	  if(numhead.eq.0) cycle travel
	  if(numarray.ne.1)then

! -- The following code is a little unweildly because it incorporates an LF90
!    bug workaround.

	    do j=1,numarray-1
	      atemp=ahead(j)
	      atemp(9:9)=','
	      ahead(j)=atemp
	    end do
	  end if
	  call num2char(distance,adistance,8)
	  adistance(9:9)=','
	  call num2char(east,aeast,11)
	  aeast(12:12)=','
	  call num2char(north,anorth,11)
	  anorth(12:12)=','
	  write(outunit,510,err=9100) adistance,aeast,anorth, &
          (ahead(j),j=1,numarray)
510	  format(1x,a9,1x,a12,1x,a12,1x,100(a9,5x))
	  numact=numact+1

	end do travel

	if(numact.eq.0) then
	  write(amessage,610) trim(outfile)
610	  format(' No transect data written to file ',a,': transect ', &
          'line does not intersect any active grid cells.')
	else
	  call num2char(numact,atemp)
	  write(amessage,620) trim(atemp),trim(outfile)
620	  format('  - ',a,' lines of data written to file ',a,'.')
	end if
	call write_message(endspace='yes')

	go to 9900
	
9000	write(amessage,9010)
9010	format(' Cannot allocate sufficient memory to run program SECTION.')
	call write_message(leadspace='yes')
	go to 9900
9100	write(amessage,9110) trim(outfile)
9110	format(' Cannot write to output file ',a,': file inaccessible ', &
        'or disk full.')
	call write_message
	go to 9900
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(realarray,realfile,bhead,ahead,stat=ierr)

end program section

 