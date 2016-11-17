!     Last change:  JD   13 Feb 2001    6:21 pm
program pt2array

! -- Program PT2ARRAY builds a real array on the basis of data contained
!    in a bore information file.

	use defn
	use inter

	implicit none

	integer                                 :: ifail,i,ierr,outunit,irow,&
						   icol,icellno,jcellno,imsg,&
						   iunit,iline,nbore,numbore,&
						   numdat,numcol,nocol,idate,iheader, &
                                                   ifail1,nbb
	real                                    :: fac1,fac2,fac3,fac4,rval,fac
	type (modelgrid) 			:: gridspec
	real, allocatable, dimension(:,:)       :: realarray
	character (len=1)			:: adivide
	character (len=80)			:: boredatfile,aprompt,bfile
	character (len=10)			:: anum1,anum2,boreid
	character (len=100)			:: amsg
	character (len=20)			:: atemp


	write(amessage,5)
5       format(' Program PT2ARRAY builds a real array on the basis of data ', &        
	'contained in a bore information file.')
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

	call readfig(gridspec%specfile,bore_coord_file)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	allocate(realarray(gridspec%ncol,gridspec%nrow),stat=ierr)
	if(ierr.ne.0)then
	  write(amessage,20)
20	  format(' Cannot allocate sufficient memory to run program PT2ARRAY.')
	  call write_message
	  go to 9900
	end if
	realarray=0.0

30      call read_bore_coord_file(ifail, &
       ' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  deallocate(realarray)
	  go to 10
	end if

150	write(6,*)
	imessage=0
160     write(6,170,advance='no')
170	format(' Enter name of bore information file: ')
	read(5,'(a)') boredatfile
	if(boredatfile.eq.' ') go to 160
	boredatfile=adjustl(boredatfile)
	if(index(eschar,boredatfile(1:2)).ne.0) then
	  write(6,*)
	  go to 30
	end if
        bfile=boredatfile
        nbb=len_trim(bfile)
        call getfile(ifail1,bfile,boredatfile,1,nbb)
        if(ifail1.ne.0) go to 160
	iunit=nextunit()
	open(unit=iunit,file=boredatfile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9900
	  write(amessage,190) trim(boredatfile)
190       format('cannot open bore information file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  go to 160
	end if

240	imessage=0
	write(6,250,advance='no')
250	format(' Use which column of bore information file to generate real array: ')
	read(5,'(a)') anum1
	if(anum1.eq.' ') go to 240
	if(index(eschar,anum1(1:2)).ne.0) go to 150
	call char2num(ifail,anum1,numcol)
	if(ifail.ne.0) go to 240
	if(pos_test(numcol,'column number').ne.0) go to 240
	if(numcol.gt.NUM_WORD_DIM)then
	  write(amessage,260)
260	  format(' Column number out of range (too large) - try again.')
	  call write_message
	  go to 240
	else if(numcol.eq.1) then
	  write(amessage,280)
280	  format(' Column 1 contains bore identifiers - try again.')
	  call write_message
	  go to 240
	end if

510	write(6,520,advance='no')
520	format(' Enter multiplication factor for data in this column  [1.00]: ')
	read(5,'(a)') atemp

	if(atemp.eq.' ')then
	  fac=1.00
	else
	  atemp=adjustl(atemp)
	  if(index(eschar,atemp(1:2)).ne.0) then
	    write(6,*)
	    go to 240
	  else
	    call char2num(ifail,atemp,fac)
	    if(ifail.ne.0) go to 510
	  end if
	end if
550	write(6,560,advance='no')
560	format(' Divide by cell area?  [y/n] ')
	read(5,'(a)') adivide
	if(adivide.eq.' ') go to 550
	if(index(eschar,adivide).ne.0) then
	  write(6,*)
	  go to 510
	end if
	call casetrans(adivide,'lo')
	if((adivide.ne.'y').and.(adivide.ne.'n')) go to 550

	write(initial_message,300) trim(boredatfile)
300      format(' Errors in bore information file ',a,' ----->')

	imessage=0
	iline=0
	nbore=0
	numbore=0
	numdat=0
	do
	  nocol=0
	  iline=iline+1
	  if(iline.eq.1000)write(6,305)
305       format(' Working.....')
	  read(iunit,'(a)',err=9000,end=200) cline
	  call linesplit(ifail,numcol)
	  if(ifail.lt.0) cycle
	  nbore=nbore+1
          if(ifail.gt.0) then
	    nocol=1
	  end if
	  call num2char(iline,anum1)
	  write(amsg,290) trim(anum1)
290       format('  Line ',a,':')
	  imsg=len_trim(amsg)+1
	  if(right_word(1)-left_word(1)+1.gt.10) then
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,310) amsg(1:imsg)
310         format(1x,a,'bore identifier greater than 10 characters in length.')
	    call write_message(increment=1)
	    boreid='*'
	  else
	    boreid=cline(left_word(1):right_word(1))
	    call casetrans(boreid,'hi')
	    do i=1,num_bore_coord
	      if(bore_coord_id(i).eq.boreid) go to 367
	    end do
	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
	    write(amessage,366) amsg(1:imsg),trim(boreid),trim(bore_coord_file)
366         format(1x,a,'no coordinates for bore ',a, &
	    ' found in bore coordinates file ',a,'.')
	    call write_message(increment=1)
	    boreid='*'
367         continue
	  end if
	  if(nocol.eq.1) cycle
	  rval=char2real(ifail,numcol)
!	  if(ifail.ne.0) then
!	    if(imessage.eq.0) call write_initial_message(leadspace='yes')
!	    call num2char(numcol,anum1)
!	    write(amessage,390) amsg(1:imsg),trim(anum1)
!390	    format(a,'cannot read number from column ',a,'.')
!	    call write_message(increment=1)
!	  else
	  if(ifail.eq.0) then
	    if(boreid.ne.'*') then
	      numbore=numbore+1
	      call factor(gridspec,bore_coord_east(i),bore_coord_north(i),&
              fac1,fac2,fac3,fac4,icellno,jcellno)
	      if(icellno.gt.0)then
	        numdat=numdat+1
	        call cell2rc(icellno,irow,icol,gridspec)
		rval=rval*fac
		if(adivide.eq.'y') rval=rval/gridspec%delr(icol)/gridspec%delc(irow)
	        realarray(icol,irow)=realarray(icol,irow)+rval
	      end if
	    end if
	  end if
	end do

200     close(unit=iunit,iostat=ierr)
	if(imessage.ne.0) go to 9900
	if(nbore.eq.0) then
	  if(imessage.eq.0) call write_initial_message(leadspace='yes')
	  write(amessage,210)
210       format('   No bores read from file.')
	  call write_message
	  go to 9900
	else
	  call num2char(nbore,anum1)
	  write(amessage,430) trim(anum1),trim(boredatfile)
430	  format('  - information for ',a,' bores found in file ',a)
	  call write_message
	  call num2char(numbore,anum1)
	  call num2char(numcol,anum2)
	  write(amessage,440) trim(anum1),trim(anum2)
440	  format('  - ',a,' instances of useable information in column ',a,&
          ' of this file')
	  call write_message
	  if(numbore.eq.0) go to 9900
	  call num2char(numdat,anum1)
	  write(amessage,450) trim(anum1)
450	  format('  - of these, ', a,' bores are within the bounds of the ',&
          'finite-difference grid')
	  call write_message
	  if(numdat.eq.0) go to 9900
	end if
	if(imessage.ne.0) go to 9900

	write(6,*)
	aprompt=' Enter name for real array output file: '
	call write_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=gridspec%nrow, columns=gridspec%ncol)
	if(escset.ne.0)then
	  escset=0
	  go to 150
	end if

	go to 9900

9000    call num2char(iline,anum1)
	write(amessage,9010) trim(anum1),trim(boredatfile)
9010    format('Error reading line ',a,' of bore information file ',a)
	call write_message(leadspace='yes')
	go to 9900

9900	call close_files
	call free_grid_mem(gridspec)
	call free_bore_mem
	deallocate(realarray,stat=ierr)

end program pt2array
 