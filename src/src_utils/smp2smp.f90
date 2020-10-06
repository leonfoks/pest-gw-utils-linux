!     Last change:  JD   12 Apr 2003    2:49 am
program smp2smp

! -- Program SMP2SMP builds a bore sample file by temporally interpolating
!    samples contained in one bore sample file to the dates and times of 
!    samples recorded in another bore sample file.

	use defn
	use inter
	implicit none

	logical              :: active
	integer, parameter   :: MAXID=8000
	integer              :: ifail,idate,obsunit,modunit,iline,nobsid,i,  &
	                        l,cols,maxmod,ierr,ind,iout,iobs,intdays,  &
	                        intsecs,ii,dd,mm,yy,hhh,mmm,sss,j,iflag,     &
	                        outunit,iheader
	integer, dimension (MAXID)                        :: nmod,loc,nout
	integer, allocatable, dimension (:)               :: ndays,nsecs

	real                 :: elim
	double precision     :: intvalue
	double precision, allocatable, dimension (:)      :: value

	character (len=5)    :: anum,aline
	character (len=10)   :: atemp1,atemp1old
	character (len=15)   :: atemp,atempold
	character (len=80)   :: obsfle,modfle,outfle,specfle
	character (len=10), dimension(MAXID)              :: obsid



	write(amessage,5)
5	format(' Program SMP2SMP builds a bore sample file by ',  &
	'interpolating samples contained in one bore sample file to the ',   &
	'dates and times recorded in another bore sample file.')
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
	if((idate.ne.0).or.(datespec.eq.0)) then
	  write(amessage,9)
9	  format(' Cannot read date format from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

	call readfig(specfle,sampfile=obsfle)

50	call open_named_input_file(ifail, &
	' Enter name of observation bore sample file: ',obsfle,obsunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

	write(6,*)
70	call open_input_file(ifail,   &
	' Enter name of model-generated bore sample file: ',modfle,modunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  close(unit=obsunit)
	  go to 50
	end if

100	write(6,110,advance='no')
110	format(' Enter extrapolation threshold in days ',                  &
	'(fractional if necessary): ')
	read(5,'(a)') anum
	if(anum.eq.' ') go to 100
	anum=adjustl(anum)
	if(index(eschar,anum(1:2)).ne.0) then
	  close(unit=modunit)
	  write(6,*)
	  go to 70
	end if
	call char2num(ifail,anum,elim)
	if(ifail.ne.0) go to 100
	if(elim.lt.0.0) go to 100
	
140	write(6,*)
150	call open_output_file(ifail,   &
	' Enter name for new bore sample file: ',outfle,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 100
	end if

! -- First the observation bore sample file is read to obtain all bore_ids.

	atempold=' '
	iline=0
	nobsid=0
	do
	  iline=iline+1
	  read(obsunit,*,end=210) atemp
	  atemp=adjustl(atemp)
	  call casetrans(atemp,'hi')
	  l=len_trim(atemp)
	  if(l.gt.10)then
	    call num2char(iline,aline)
	    write(amessage,160) trim(aline),trim(obsfle)
160	    format(' Identifier greater than 10 characters at ',  &
	    'line ',a,' of file ',a)
	    go to 9890
	  end if
	  if(nobsid.eq.0)then
	    nobsid=nobsid+1
	    if(nobsid.gt.MAXID)then
	      write(amessage,180) trim(obsfle)
180	      format(' Too many bores cited in file ',a,                    &
	      ': increase MAXID and re-compile program.')
	      go to 9890
	    end if
	    obsid(nobsid)=atemp(1:10)
	    atempold=atemp
	  else if(atemp.eq.atempold) then
	    continue
	  else
	    call casetrans(atemp,'hi')
	    do i=1,nobsid
	      if(atemp.eq.obsid(i))then
	        call num2char(iline,aline)
	        write(amessage,200) trim(aline),trim(obsfle)
200	        format(' Identifier used previously at line ',a,          &
	        ' of file ',a)
	        go to 9890
	      end if
	    end do
	    nobsid=nobsid+1
	    if(nobsid.gt.MAXID)then
	      write(amessage,180) trim(obsfle)
	      go to 9890
	    end if
	    obsid(nobsid)=atemp(1:10)
	    atempold=atemp
	  end if
	end do
210	rewind(unit=obsunit,iostat=ierr)
	if(ierr.ne.0)then
	  write(amessage,245) trim(obsfle)
245	  format(' Cannot rewind file ',a)
	  go to 9890
	end if


! -- Next the model-generated bore sample file is read to get some 
!    specifications.

250	continue
	nmod=0				!nmod is an array
	iline=0
	atemp1old=' '
	read_model_sample_file: do
	  iline=iline+1
	  read(modunit,'(a)',end=300) cline
	  call linesplit(ifail,1)
	  if(ifail.lt.0) cycle read_model_sample_file
	  if(right_word(1)-left_word(1).gt.9) then
	    call num2char(iline,aline)
	    write(amessage,280) trim(aline),trim(modfle)
280	    format(' Identifier greater than 10 characters in length at ',  &
	    'line ',a,' of file ',a)
	    go to 9890
	  end if
	  atemp1=cline(left_word(1):right_word(1))
	  atemp1=adjustl(atemp1)
	  call casetrans(atemp1,'hi')
	  if(atemp1.eq.atemp1old)then
	    if(active) then
	      nmod(iobs)=nmod(iobs)+1
	    else
	      cycle read_model_sample_file
	    end if
	  else
	    atemp1old=atemp1
	    do i=1,nobsid
	      if(atemp1.eq.obsid(i))then
	        active=.true.
	        iobs=i
	        if(nmod(iobs).ne.0)then
	          call num2char(iline,aline)
	          write(amessage,290) trim(aline),trim(modfle)
290	          format(' Identifiers with same name not in juxtaposition ',&
	          'at line ',a,' of file ',a)
	          go to 9890
	        end if
		nmod(iobs)=1
	        cycle read_model_sample_file
	      end if
	    end do
	    active=.false.
	  end if
	end do read_model_sample_file

! -- Space is allocated for the storage of all pertinent values from the 
!    model-generated bore sample file.
	
300	maxmod=0
	do i=1,nobsid
	  maxmod=maxmod+nmod(i)
	end do
	if(maxmod.eq.0)then
	  write(amessage,310) trim(modfle),trim(obsfle)
310	  format(' No identifiers cited in file ',a,' correspond to ',      &
	  'any of the identifiers in file ',a)
	  go to 9890
	end if
	allocate(ndays(maxmod),nsecs(maxmod),value(maxmod),stat=ierr)
	if(ierr.ne.0)then
	  write(amessage,330)
330	  format(' Cannot allocate sufficient memory to continue execution.')
	  go to 9890
	end if
	rewind(unit=modunit,iostat=ierr)
	if(ierr.ne.0)then
	  write(6,340) trim(modfle)
340	  format(' Cannot rewind file ',a)
	  go to 9890
	end if

!	do i=1,nobsid			!debug
!	  write(6,*) obsid(i),nmod(i)
!	enddo

! -- The model-generated bore sample file is next re-read and sample values
!    stored.

	iline=0
	loc=0				!loc is an array
	ind=0
	atemp1old=' '
	read_model_sample_file_1: do
	  iline=iline+1
	  read(modunit,'(a)',end=390) cline
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_model_sample_file_1
	  if(ifail.gt.0)then
	    cols=4
	    call linesplit(ifail,4)
	    if(ifail.gt.0)then
	      call num2char(iline,aline)
	      write(amessage,360) trim(aline),trim(modfle)
360	      format(' Insufficient entries at line ',a,' of file ',a)
	      go to 9890
	    end if
	  end if
	  atemp1=cline(left_word(1):right_word(1))
	  call casetrans(atemp1,'hi')
	  atemp1=adjustl(atemp1)
	  if(atemp1.eq.atemp1old)then
	    if(.not.active)then
	      cycle read_model_sample_file_1
	    end if
	    ind=ind+1
	  else
	    atemp1old=atemp1
	    do i=1,nobsid
	      if(atemp1.eq.obsid(i))then
	        ind=ind+1
	        loc(i)=ind
	        active=.true.
	        go to 370
	      end if
	    end do
	    active=.false.
	    cycle read_model_sample_file_1
	  end if
370	  call read_rest_of_sample_line(ifail,cols,ndays(ind),           &
	  nsecs(ind),value(ind),iline,modfle)
	  if(ifail.ne.0) go to 9900
	  if(value(ind).lt.-1.0e38)then
	    ind=ind-1
	    nmod(i)=nmod(i)-1
	  endif
	end do read_model_sample_file_1
390	close(unit=modunit)

! --  Next the observation bore sample file is read line by line and entries
!     for the output bore sample file written where appropriate.

	iobs=0
	iout=0
	iline=0
	nout=0					! nout is an array
	atemp1old=' '
	read_obs_file: do
	  iline=iline+1
!	  write(6,*) iline			!debug
	  read(obsunit,'(a)',end=500) cline
	  cols=5
	  call linesplit(ifail,5)
	  if(ifail.lt.0) cycle read_obs_file
	  if(ifail.gt.0)then
	    cols=4
	    call linesplit(ifail,4)
	    if(ifail.ne.0) then
	      call num2char(iline,aline)
	      write(amessage,360) trim(aline),trim(obsfle)
	      go to 9890
	    end if
	  end if
	  atemp1=cline(left_word(1):right_word(1))
	  call casetrans(atemp1,'hi')
	  atemp1=adjustl(atemp1)
	  if(atemp1.ne.atemp1old)then
	    iobs=iobs+1
	    atemp1old=atemp1
	  end if
	  if(nmod(iobs).eq.0) cycle read_obs_file
	  call read_rest_of_sample_line(ifail,cols,intdays,intsecs,intvalue, &
	  iline,obsfle)
	  if(ifail.ne.0) go to 9900
	  if(intvalue.lt.-1.0e38) cycle read_obs_file
	  ii=loc(iobs)
	  call time_interp(ifail,nmod(iobs),ndays(ii),nsecs(ii),value(ii),   &
	  intdays,intsecs,1.0e30,elim,intvalue)
	  if(ifail.ne.0)then
	    write(amessage,380) trim(modfle)
380	    format(' Problem with sample file ',a,'; run SMPCHEK for ',      &
	    'more information.')
	    go to 9890
	  end if
	  if(intvalue.gt.-1.0e30)then
	    iout=iout+1
	    call newdate(intdays,1,1,1970,dd,mm,yy)
	    hhh=intsecs/3600
	    mmm=(intsecs-hhh*3600)/60
	    sss=intsecs-hhh*3600-mmm*60
	    if(datespec.eq.1) then
	      write(outunit,400) trim(atemp1),dd,mm,yy,hhh,mmm,sss,intvalue
400	      format(1x,a,t15,i2.2,'/',i2.2,'/',i4.4,3x,i2.2,':',i2.2,':',   &
	      i2.2,3x,1pg15.8)
	    else
	      write(outunit,400) trim(atemp1),mm,dd,yy,hhh,mmm,sss,intvalue
	    endif
	    nout(iobs)=nout(iobs)+1
	  endif
	end do read_obs_file

500	call num2char(iout,anum)
	write(amessage,510) trim(anum),trim(outfle)
510	format(' - ',a,' data lines written to new bore sample file ',a)
	call write_message
	close(unit=obsunit)
	close(unit=modunit)

	iflag=0
	do i=1,nobsid
	  if(nmod(i).eq.0) then
	    iflag=1
	    exit
	  end if
	end do
	if(iflag.ne.0)then
	  write(amessage,540) trim(obsfle),trim(modfle)
540	  format(' The following identifiers cited in observaton sample ',   &
	  'file ',a,' are either uncited in model observation file ',a,      &
	  ' or are all x-affected in that file.')
	  call write_message(leadspace='yes')
	  amessage=' '
	  j=4
	  do i=1,nobsid
	    if(nmod(i).eq.0)then
	      write(amessage(j:),'(a)') trim(obsid(i))
	      j=j+11
	      if(j.gt.69)then
	        call write_message
	        amessage=' '
	        j=4
	      end if
	    end if
	  end do
	  if(j.ne.4) call write_message
	end if

	iflag=0
	do i=1,nobsid
	  if((nmod(i).ne.0).and.(nout(i).eq.0)) then
	    iflag=1
	    exit
	  end if
	end do
	if(iflag.ne.0)then
	  write(amessage,580) trim(obsfle),trim(modfle)
580	  format(' The following identifiers cited in observaton sample ',   &
	  'file ',a,' are cited in model sample file ',a,'. However no ',    &
	  'observation times are within model simulation time frame.')
	  call write_message(leadspace='yes')
	  amessage=' '
	  j=4
	  do i=1,nobsid
	    if((nmod(i).ne.0).and.(nout(i).eq.0))then
	      write(amessage(j:),'(a)') trim(obsid(i))
	      j=j+11
	      if(j.gt.69)then
	        call write_message
	        amessage=' '
	        j=4
	      end if
	    end if
	  end do
	  if(j.ne.4) call write_message
	end if
	

	go to 9900



9890	call write_message(leadspace='yes')
9900	call close_files
	deallocate(ndays,nsecs,value,stat=ierr)
	write(6,*)

end program smp2smp
 