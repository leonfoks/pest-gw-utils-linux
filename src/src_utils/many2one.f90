!     Last change:  JD    9 May 2003    8:56 am
program many2one

! -- Program MANY2ONE allows a user to selectively store individual arrays
!    contained in a multiple-array MODFLOW or MT3D unformatted output file.


	use defn
	use inter

	implicit none

	logical					:: lexist
	integer					:: iunit,ierr,iarray,kstp,kper,&
                                                   ncol,nrow,ilay,ntrans,&
						   nrowold,ncolold,ifail,idate,iheader, &
                                                   ifail1,nbb
	real 					:: pertim,totim
	real, allocatable, dimension(:,:)	:: realarray
	character (len=120)			:: mfile,aprompt,bfile
	character (len=5)			:: aformat
	character (len=7)			:: outext,atype
	character (len=16)			:: text,textl
	character (len=20)			:: anum1
	character (len=1)			:: aa,ayn


	write(amessage,5)
5       format(' Program MANY2ONE allows a user to selectively store individual arrays ',&
	'contained in a multiple-array MODFLOW or MT3D unformatted output file.')
	call write_message(leadspace='yes',endspace='yes')

        include 'unformat.inc'
	
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

	ncolold=0
	nrowold=0
9	imessage=0
10	if(imessage.gt.5) go to 9900
	write(6,20,advance='no')
20	format(' Enter name of MODFLOW/MT3D unformatted output file: ')
	read(5,'(a)') mfile
	if(mfile.eq.' ') go to 10
	mfile=adjustl(mfile)
	if(index(eschar,mfile(1:2)).ne.0) go to 9900
        nbb=len_trim(mfile)
        call getfile(ifail1,mfile,bfile,1,nbb)
        if(ifail1.ne.0) go to 10
        mfile=bfile
	inquire(exist=lexist,file=mfile)
	if(.not.lexist)then
	  write(amessage,50) trim(mfile)
50	  format('cannot open file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  go to 10
	end if
	iunit=nextunit()
	open(unit=iunit,file=mfile,status='old',form='binary',iostat=ierr)
	if(ierr.ne.0)then
	  open(unit=iunit,file=mfile,status='old',form='unformatted', &
	  iostat=ierr)
	  if(ierr.eq.0) go to 90
	  inquire(file=mfile,unformatted=aformat,iostat=ierr)
	  aformat=adjustl(aformat)
	  if(ierr.ne.0) then
	    write(amessage,50) trim(mfile)
	    call write_message(error='yes',increment=1)
	    go to 10
	  end if
	  call casetrans(aformat,'lo')
	  if(aformat.ne.'yes') then
	    write(amessage,80) trim(mfile)
80	    format(' File ',a,' does not appear to be an unformatted file ',&
            ' - try again.')
	    call write_message(increment=1)
	    go to 10
	  else
	    write(amessage,50) trim(mfile)
	    call write_message(error='yes',increment=1)
	    go to 10
	  end if
	end if
	
90	write(6,100,advance='no')
100	format(' Is this a MODFLOW or MT3D output file?  [f/t]: ')
	read(5,'(a)') aa
	if(aa.eq.' ') go to 90
	if(index(eschar,aa).ne.0) then
	  close(unit=iunit,iostat=ierr)
	  write(6,*)
	  go to 9
	end if
	call casetrans(aa,'lo')
	if((aa.ne.'f').and.(aa.ne.'t')) go to 90
	if(aa.eq.'f')then
	  outext='MODFLOW'
	  atype='modflow'
	else
	  outext='MT3D'
	  atype='mt3d'
	end if
	iarray=0

120	iarray=iarray+1
	if(aa.eq.'f')then
	  read(iunit,err=9000,end=8000) kstp,kper,pertim,totim,text,ncol,&
          nrow,ilay
	else
	  read(iunit,err=9000,end=8000) ntrans,kstp,kper,totim,text,ncol,&
          nrow,ilay
	end if
	if((nrow.le.0).or.(ncol.le.0))then
	  call num2char(iarray,anum1)
	  write(amessage,125) trim(anum1),trim(outext),trim(mfile)
125	  format('grid row or column dimension supplied as zero or negative ',&
	  'in header to array number ',a,' in ',a,' output file ',a)
	  call write_message(error='yes',leadspace='yes',endspace='yes')
	  go to 9900
	end if
	if((ncol.ne.ncolold).or.(nrow.ne.nrowold))then
	  if(allocated(realarray)) deallocate(realarray,stat=ierr)
	  if(ierr.ne.0)then
	    write(amessage,130)
130	    format(' MANY2ONE memory management error.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	  allocate(realarray(ncol,nrow),stat=ierr)
	  if(ierr.ne.0)then
	    write(amessage,140)
140	    format(' Cannot allocate sufficient memory to run MANY2ONE.')
	    call write_message(leadspace='yes',endspace='yes')
	    go to 9900
	  end if
	end if
	read(iunit,err=9100,end=9100) realarray
	write(6,*)
	if(aa.eq.'f') then
	  call num2char(ilay,anum1)	
	  textl=adjustl(text)
	  call casetrans(textl,'lo')
	  write(6,270) trim(outext),trim(textl),trim(anum1)
270	  format(1x,a,1x,a,' array for layer ',a,' ----->')
	  call num2char(kper,anum1)
	  write(6,271) trim(anum1)
271	  format(' Stress period',t45,'= ',a)
	  call num2char(kstp,anum1)
	  write(6,272) trim(anum1)
272	  format(' Time step',t45,'= ',a)
	  call num2char(pertim,anum1,8)
	  write(6,240) trim(anum1)
240	  format(' Elapsed time since start of stress period',t45,'= ',a)
	  call num2char(totim,anum1,8)
	  write(6,250) trim(anum1)
250	  format(' Elapsed time since start of simulation',t45,'= ',a)
	else
	  textl=adjustl(text)
	  call casetrans(textl,'lo')
	  call num2char(ilay,anum1)
	  write(6,270) trim(outext),trim(textl),trim(anum1)
	  call num2char(kper,anum1)
	  write(6,271) trim(anum1)
	  call num2char(kstp,anum1)
	  write(6,272) trim(anum1)
	  call num2char(ntrans,anum1)
	  write(6,275) trim(anum1)
275	  format(' Transport step',t45,'= ',a)
	  call num2char(totim,anum1,8)
	  write(6,250) trim(anum1)
	end if
300	write(6,310,advance='no')
310	format(' Store array in separate file?  [y/n]: ')
	read(5,'(a)') ayn
	if(ayn.eq.' ') go to 300	
	if(index(eschar,ayn).ne.0) then
	  close(unit=iunit)
	  write(6,*)
	  go to 9
	end if
	call casetrans(ayn,'lo')
	if((ayn.ne.'y').and.(ayn.ne.'n')) go to 300
	if(ayn.eq.'y') then
	  aprompt=' Enter name for real array file: '
	  call write_real_array(ifail,aprompt,realarray,pm_header=headerspec,&
	  rows=nrow,columns=ncol,binary_header='yes',atype=atype,ntrans=ntrans,&
	  kstp=kstp,kper=kper,pertim=pertim,totim=totim,text=text,ncol=ncol,&
	  nrow=nrow,ilay=ilay)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    write(6,*)
	    go to 300
	  end if
	end if
	go to 120


8000	if(iarray.eq.1)then
	  write(amessage,8010) trim(outext),trim(mfile)
8010	  format('no arrays found in ',a,' output file ',a)
	  call write_message(error='yes',leadspace='yes',endspace='yes')
	  go to 9900
	end if
	call num2char(iarray-1,anum1)	
	write(amessage,8050) trim(anum1),trim(mfile)
8050	format(' No more arrays: ',a,' arrays read from file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 9900

9000	call num2char(iarray,anum1)
	write(amessage,9010) trim(anum1),trim(outext),trim(mfile)
9010	format(' Error encountered while reading header to array number ',a, &
        ' from unformatted ',a,' output file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 9900
9100	call num2char(iarray,anum1)
	write(amessage,9110) trim(anum1),trim(outext),trim(mfile)
9110	format(' Error encountered while reading array number ',a, &
        ' from unformatted ',a,' output file ',a)
	call write_message(leadspace='yes',endspace='yes')
	go to 9900

9900    call close_files
	deallocate(realarray,stat=ierr)

end program many2one
 