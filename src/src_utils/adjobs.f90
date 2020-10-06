!     Last change:  JD   17 Jul 2001   12:44 pm
program adjobs

! -- Program ADJOBS adjusts observation weights and assigns observation 
!    groups in an existing PEST control file.


	use defn
	use inter
	implicit none

	integer, parameter      :: MAXGROUP=200
	integer, parameter      :: MAXSUBGROUP=20
	integer                 :: ifail,punit1,i,iline,nobs,nobsgp,ierr,  &
	                           j,punit2,itemp,n,k,nprior,npargp,npar,m,l, &
	                           iobs,nnn

	integer, dimension(MAXGROUP)                 :: nsubgp,nn,mobsgp

	integer, allocatable,dimension(:)            :: iobsgp
	real                                         :: rtemp
	real, allocatable, dimension(:)              :: oval,weight
	real, dimension(MAXGROUP,MAXSUBGROUP)        :: w1,w2,w3,w4,w5,w6


	character (len=5)        :: aline,anum
	character (len=200)      :: pestfle1,pestfle2
	character (len=20)       :: atemp
        character (len=20)       :: atemp1

	character (len=1), dimension(MAXGROUP)        :: aa,ag
	character (len=12), dimension(MAXGROUP)       :: aobsgp
	character (len=1), dimension(MAXGROUP,MAXSUBGROUP)  :: aw
	character (len=12), dimension(MAXGROUP,MAXSUBGROUP) :: aasubgp
        character (len=20), dimension(MAXGROUP,MAXSUBGROUP) :: asubgp
	character (len=20), allocatable, dimension(:) :: aobs


	nsubgp=0				!an array
	ag='n'					!an array



	write(amessage,5)
5	format(' Program ADJOBS adjusts observation weights and re-assigns ',&
	'observation groups in an existing PEST control file.')
	call write_message(leadspace='yes')
        write(amessage,6)
6       format(' Note that weights assigned to prior information equations ', &
        'are NOT adjusted even if they belong to the same groups as observations.')
	call write_message(leadspace='yes',endspace='yes')

50	call open_input_file(ifail,   &
	' Enter name of existing PEST control file: ',pestfle1,punit1)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

! -- The existing PEST control file is read, and relevant parts of it 
!    stored.

	iline=1
	read(punit1,'(a)',end=9000) cline
	call casetrans(cline,'lo')
	if(cline(1:3).ne.'pcf')then
	  write(amessage,60) trim(pestfle1)
60	  format(' File ',a,' does not appear to be a PEST control file ',  &
	  '- try again.')
	  call write_message(leadspace='yes',endspace='yes')
	  close(unit=punit1)
	  go to 50
	end if
	do i=1,3
70	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  if(cline.eq.' ') go to 70
	end do
	call linesplit(ifail,5)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,80) trim(aline),trim(pestfle1)
80	  format(' Five entries are expected on line ',a,' of file ',a,     &
	  '. Either the file is in error or it uses an old PEST file format.')
	  go to 9890
	end if
	npar=char2int(ifail,1)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,85) trim(aline),trim(pestfle1)
85	  format(' Error reading number of parameters from line ',a         &
	  ,' of PEST control file ',a)
	  go to 9890
	end if
	nobs=char2int(ifail,2)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,90) trim(aline),trim(pestfle1)
90	  format(' Error reading number of observations from line ',a       &
	  ,' of PEST control file ',a)
	  go to 9890
	end if
	npargp=char2int(ifail,3)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,92) trim(aline),trim(pestfle1)
92	  format(' Error reading number of parameter groups from line ',a    &
	  ,' of PEST control file ',a)
	  go to 9890
	end if
	nprior=char2int(ifail,4)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,95) trim(aline),trim(pestfle1)
95	  format(' Error reading number of prior information items',         &
          ' from line ',a,' of PEST control file ',a)
	  go to 9890
	end if
	nobsgp=char2int(ifail,5)
	if(ifail.ne.0)then
	  call num2char(iline,aline)
	  write(amessage,100) trim(aline),trim(pestfle1)
100	  format(' Error reading number of observation groups from line ',a  &
	  ,' of PEST control file ',a)
	  go to 9890
	end if
	if(nobsgp.gt.MAXGROUP) then
	  write(amessage,102) trim(pestfle1)
102	  format(' Too many observation groups in file ',a,                  &
	  ':  increase MAXGROUP and re-compile program.')
	  go to 9890
	end if
	if(nobsgp.le.0)then
	  call num2char(iline,aline)
	  write(amessage,103) trim(aline),trim(pestfle1)
103	  format(' Number of observation groups cannot be less than 1 on ',  &
	  'line ',a,' of PEST control file ',a)
	  go to 9890
	end if
	allocate(aobs(nobs),iobsgp(nobs),weight(nobs),oval(nobs),            &
	stat=ierr)
	if(ierr.ne.0) go to 9100

	write(amessage,105) trim(pestfle1)
105	format(' Reading observation data from PEST control file ',a,'......')
	call write_message(leadspace='yes')

	do
	  iline=iline+1
	  read(punit1,'(a)',end=9050) cline
	  if(cline(1:1).ne.'*') cycle
	  call casetrans(cline,'lo')
	  if(index(cline,'* observation gr').eq.0) cycle
	  go to 150
	end do
150	do i=1,nobsgp
160	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  call linesplit(ifail,1)
	  if(ifail.lt.0) goto 160
	  if(right_word(1)-left_word(1).gt.11)then
	    call num2char(iline,aline)
	    write(amessage,180) trim(aline),trim(pestfle1)
180	    format(' Observation group name longer than 12 characters at ',   &
	    'line ',a,' of file ',a)
	    go to 9890
	  end if
	  aobsgp(i)=cline(left_word(1):right_word(1))
	  aobsgp(i)=adjustl(aobsgp(i))
	  call casetrans(aobsgp(i),'lo')
	end do

200	iline=iline+1
	read(punit1,'(a)',end=9000) cline
	if(cline.eq.' ') go to 200
	call casetrans(cline,'lo')
	cline=adjustl(cline)
	if(index(cline,'* observation da').eq.0)go to 9150

	do i=1,nobs
220	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  call linesplit(ifail,4)
	  if(ifail.lt.0) go to 220
	  if(ifail.ne.0)then
	    call num2char(iline,aline)
	    write(amessage,240) trim(aline),trim(pestfle1)
240	    format(' Insufficient entries on line ',a,' of file ',a)
	    go to 9890
	  end if
	  if(right_word(1)-left_word(1).gt.19)then
	    call num2char(iline,aline)
	    write(amessage,250) trim(aline),trim(pestfle1)
250	    format(' Observation name greater than 20 characters at line ',a  &
	    ,' of file ',a)
	    go to 9890
	  end if
	  aobs(i)=cline(left_word(1):right_word(1))
	  oval(i)=char2real(ifail,2)
	  if(ifail.ne.0)then
	    call num2char(iline,aline)
	    write(amessage,260) trim(aline),trim(pestfle1)
260	    format(' Error reading observation value from line ',a,          &
	    ' of file ',a)
	    go to 9890
	  end if
	  weight(i)=char2real(ifail,3)
	  if(ifail.ne.0)then
	    call num2char(iline,aline)
	    write(amessage,270) trim(aline),trim(pestfle1)
270	    format(' Error reading observation weight from line ',a          &
	    ,' of file ',a)
	    go to 9890
	  end if
	  if(right_word(4)-left_word(4).gt.11)then
	    call num2char(iline,aline)
	    write(amessage,280) trim(aline),trim(pestfle1)
280	    format(' Observation group name greater than 12 characters at ',  &
	    'line ',a,' of file ',a)
	    go to 9890
	  end if
	  atemp=cline(left_word(4):right_word(4))
	  atemp=adjustl(atemp)
	  call casetrans(atemp,'lo')
	  do j=1,nobsgp
	    if(atemp(1:12).eq.aobsgp(j))then
	      iobsgp(i)=j
	      go to 281
	    end if
	  end do
	  call num2char(iline,aline)
	  write(amessage,284) trim(aline),trim(pestfle1)
284	  format(' Unidentified observation group at line ',a,               &
	  ' of file ',a)
	  go to 9890
281	  continue
	end do
282	rewind(unit=punit1,iostat=ierr)
	if(ierr.ne.0)then
	  write(amessage,286) trim(pestfle1)
286	  format(' Cannot rewind file ',a)
	  go to 9890
	end if
	do i=1,nobsgp
	  mobsgp(i)=0
	  do j=1,nobs
	    if(iobsgp(j).eq.i) mobsgp(i)=mobsgp(i)+1
	  end do
	end do


! -- The user is now queried as to his/her intentions.

	i=0
285	i=i+1
	if(i.gt.nobsgp) go to 1000
290	write(6,*)
	call num2char(mobsgp(i),anum)
	write(6,300) trim(aobsgp(i)),trim(anum)
300	format(' Observation group "',a,'" (',a,' observations belong ',     &
	'to this group) ---->')
	if(mobsgp(i).eq.0)then
	  write(amessage,315)
315	  format('   No adjustments warranted as this group has zero observation members.')
	  call write_message
	  write(6,316,advance='no')
316	  format('   Press <Enter> to continue (Remember: "e" + <Enter> ',   &
	  ' for escape).')
	  read(5,'(a)') atemp
	  if(index(eschar,atemp(1:2)).ne.0)then
	    if(i.eq.1)then
	      close(unit=punit1)
	      deallocate(aobs,iobsgp,weight,oval,stat=ierr)
	      if(ierr.ne.0) go to 9200
	      write(6,*)
	      go to 50
	    else
	      i=i-1
	      go to 290
	    end if
	  end if
	  go to 285
	end if
310	write(6,320,advance='no')
320	format('   Do you wish to make any adjustments  [y/n]?: ')
	read(5,'(a)') aa(i)
	if(aa(i).eq.' ') go to 310
	if(index(eschar,aa(i)).ne.0)then
	  if(i.eq.1)then
	    close(unit=punit1)
	    deallocate(aobs,iobsgp,weight,oval,stat=ierr)
	    if(ierr.ne.0) go to 9200
	    write(6,*)
	    go to 50
	  else
	    i=i-1
	    go to 290
	  end if
	end if
	if((aa(i).eq.'y').or.(aa(i).eq.'Y')) then
	  aa(i)='y'
	else if((aa(i).eq.'n').or.(aa(i).eq.'N'))then
	  aa(i)='n'
	  go to 285
	else
	  go to 310
	end if
340	write(6,360,advance='no')
360	format('   Divide group into subgoups  [y/n]? ')
	read(5,'(a)') ag(i)
	if(ag(i).eq.' ') go to 340
	if(index(eschar,ag(i)).ne.0) then
	  go to 290
	end if
	if((ag(i).eq.'y').or.(ag(i).eq.'Y'))then
	  ag(i)='y'
600	  write(6,610)
610	  format('   Use first n characters of observation name for ',       &
	  'group definition.')
620	  write(6,630,advance='no')
630	  format('   Enter n: ')
	  itemp=key_read(n)
	  if(itemp.ne.0) go to 620
	  if(escset.ne.0)then
	    write(6,*)
	    escset=0
	    go to 340
	  end if
	  if(n.le.0) go to 620
	  if(n.gt.19) go to 620
	  nn(i)=n
	  nsubgp(i)=0
	  do j=1,nobs
	    if(iobsgp(j).ne.i) cycle
	    atemp=aobs(j)(1:n)
	    atemp=adjustl(atemp)
	    call casetrans(atemp,'lo')
	    if(nsubgp(i).eq.0) then
	      nsubgp(i)=nsubgp(i)+1
	      asubgp(i,1)=atemp
	      asubgp(i,1)=adjustl(asubgp(i,1))
	    else
	      do k=1,nsubgp(i)
	        if(atemp(1:n).eq.asubgp(i,k)(1:n))go to 680
	      end do
	      nsubgp(i)=nsubgp(i)+1
	      if(nsubgp(i).gt.MAXSUBGROUP)then
	        write(amessage,650) trim(aobsgp(i))
650	        format(' Too many subgroups in observation group "',a,       &
	        '" - try again, or increase MAXSUBGROUP and re-compile program.')
	        call write_message(leadspace='yes',endspace='yes')
	        go to 600
	      end if
	      asubgp(i,nsubgp(i))=atemp
	      asubgp(i,nsubgp(i))=adjustl(asubgp(i,nsubgp(i)))
680	      continue
	    end if
	  end do
700	  call num2char(nsubgp(i),anum)  
	  write(6,710) trim(anum)
710	  format('   ',a,' subgroups have been defined.')
	  j=0
715	  j=j+1
	  if(j.gt.nsubgp(i)) go to 285
	  write(6,*)
720	  write(6,730) trim(aobsgp(i)),trim(asubgp(i,j))
730	  format('     Observations in group "',a,'" beginning with "',   &
	  a,'" --->')
750	  write(6,760,advance='no')
760	  format('     Enter observation group name for these observations: ')
	  read(5,'(a)') atemp1
          atemp=atemp1
	  if(atemp.eq.' ') go to 750
	  if(index(eschar,atemp(1:2)).ne.0) then
	    write(6,*)
	    if(j.eq.1)then
	      go to 600
	    else
	      j=j-1
	      go to 720
	    end if
	  end if
          atemp1=adjustl(atemp1)
          if(len_trim(atemp1).gt.12)then
            write(amessage,761)
761         format(' Observation group name greater than 12 characters - try again.')
            call write_message
            go to 750
          end if
	  atemp=adjustl(atemp)
	  call casetrans(atemp,'lo')
	  do k=1,nobsgp
	    if(atemp(1:12).eq.aobsgp(k)) then
	      write(amessage,765) trim(aobsgp(k))
765	      format(' There is already an observation group named "',a,   &
	      '" - try again.')
	      call write_message
	      go to 750
	    end if
	  end do
	  do l=1,i
	    if(ag(l).eq.'y')then
	      if(l.eq.i)then
	        nnn=j-1
	      else
	        nnn=nsubgp(l)
	      end if
	      do m=1,nnn
	        if(atemp(1:12).eq.aasubgp(l,m))then
	          write(amessage,766) 
766	          format(' You have already assigned this name '         &
	          'to another subgroup - try again.')
	          call write_message
	          go to 750
	        end if
	      end do
	    end if
	  end do
	  aasubgp(i,j)=atemp(1:12)
770	  write(6,780,advance='no')
780	  format('     Adjust weights for this observation subgroup  [y/n]? ')
	  read(5,'(a)') aw(i,j)
	  if(aw(i,j).eq.' ') go to 770
	  if(index(eschar,aw(i,j)).ne.0) then
	    write(6,*)
	    go to 720
	  end if
	  if((aw(i,j).eq.'y').or.(aw(i,j).eq.'Y'))then
	    aw(i,j)='y'
790	    write(6,800)
800	    format('       Weights are calculated as ',                    &
	    'w = a*(abs[observation_value]**b)+c')
810	    write(6,820,advance='no')
820	    format('       Enter a: ')
	    itemp=key_read(w1(i,j))
	    if(itemp.ne.0) go to 810
	    if(escset.eq.1)then
	      escset=0
	      write(6,*)
	      go to 770
	    end if
850	    write(6,860,advance='no')
860	    format('       Enter b: ')
	    itemp=key_read(w3(i,j))
	    if(itemp.ne.0) go to 850
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 810
	    end if
870	    write(6,880,advance='no')
880	    format('       Enter c: ')
	    itemp=key_read(w4(i,j))
	    if(itemp.ne.0) go to 870
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 850
	    end if
890	    write(6,900,advance='no')
900	    format('       Enter maximum allowable weight: ')
	    itemp=key_read(w5(i,j))
	    if(itemp.ne.0) go to 890
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 870
	    end if
	    if(w5(i,j).lt.0.0) go to 890
910	    write(6,920,advance='no')
920	    format('       Enter minimum allowable weight: ')
	    itemp=key_read(w6(i,j))
	    if(itemp.ne.0) go to 910
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 890
	    end if
	    if(w6(i,j).lt.0.0) go to 910
	    if(w6(i,j).gt.w5(i,j))then
	      write(amessage,525)
	      call write_message(leadspace='yes',endspace='yes')
	      go to 910
	    end if
	  else if((aw(i,j).eq.'n').or.(aw(i,j).eq.'N'))then
	    aw(i,j)='n'
	  else
	    go to 770
	  end if
	  go to 715
	else if((ag(i).eq.'N').or.(ag(i).eq.'n'))then
	  ag(i)='n'
370	  write(6,380,advance='no')
380	  format('   Adjust weights for this observation group  [y/n]? ')
	  read(5,'(a)') aw(i,1)
	  if(aw(i,1).eq.' ') go to 370
	  if(index(eschar,aw(i,1)).ne.0) go to 290
	  if((aw(i,1).eq.'y').or.(aw(i,1).eq.'Y'))then
	    aw(i,1)='y'
390	    write(6,400)
400	    format('     Weights are calculated as ',                      &
	    'w = a*(abs[observation_value]**b)+c')
410	    write(6,420,advance='no')
420	    format('     Enter a: ')
	    itemp=key_read(w1(i,1))
	    if(itemp.ne.0) go to 410
	    if(escset.eq.1)then
	      escset=0
	      write(6,*)
	      go to 370
	    end if
450	    write(6,460,advance='no')
460	    format('     Enter b: ')
	    itemp=key_read(w3(i,1))
	    if(itemp.ne.0) go to 450
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 410
	    end if
470	    write(6,480,advance='no')
480	    format('     Enter c: ')
	    itemp=key_read(w4(i,1))
	    if(itemp.ne.0) go to 470
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 450
	    end if
490	    write(6,500,advance='no')
500	    format('     Enter maximum allowable weight: ')
	    itemp=key_read(w5(i,1))
	    if(itemp.ne.0) go to 490
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 470
	    end if
	    if(w5(i,1).lt.0.0) go to 490
510	    write(6,520,advance='no')
520	    format('     Enter minimum allowable weight: ')
	    itemp=key_read(w6(i,1))
	    if(itemp.ne.0) go to 510
	    if(escset.eq.1) then
	      escset=0
	      write(6,*)
	      go to 490
	    end if
	    if(w6(i,1).lt.0.0) go to 510
	    if(w6(i,1).gt.w5(i,1))then
	      write(amessage,525)
525	      format(' Cannot exceed maximum allowable weight.')
	      call write_message(leadspace='yes',endspace='yes')
	      go to 510
	    end if
	  else if((aw(i,1).eq.'n').or.(aw(i,1).eq.'N'))then
	    aw(i,1)='n'
	  else
	    goto 370
	  end if
	  go to 285
	else
	  go to 340
	end if

! -- The new PEST control file is now written.

1000	write(6,*)
1050	call open_output_file(ifail,  &
	' Enter name for new PEST control file: ',pestfle2,punit2)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  i=nobsgp-1
	  go to 285
	end if

	iline=1
	read(punit1,'(a)',end=9000) cline
	write(punit2,'(a)') trim(cline)
	do i=1,2
1070	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  if(cline.eq.' ') go to 1070
	  write(punit2,'(a)') trim(cline)
	end do
1080	iline=iline+1
	read(punit1,'(a)',end=9000) cline
	if(cline.eq.' ') go to 1080
	n=0
	do i=1,nobsgp
	  if(ag(i).eq.'n')then
	    n=n+1
	  else
	    n=n+nsubgp(i)
	  end if
	end do
	write(punit2,'(5i10)') npar,nobs,npargp,nprior,n
	do
	  iline=iline+1
	  read(punit1,'(a)') cline
	  if(cline.eq.' ') cycle
	  if(index(cline,'* observation gr').ne.0) exit
	  write(punit2,'(a)') trim(cline)
	end do
	write(punit2,1100)
1100	format('* observation groups')
	do
	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  if(index(cline,'* observation da').ne.0) exit
	end do
	do i=1,nobsgp
	  if(ag(i).eq.'n')then
	    write(punit2,'(a)') trim(aobsgp(i))
	  else
	    do j=1,nsubgp(i)
	      write(punit2,'(a)') aasubgp(i,j)
	    end do
	  end if
	end do
	write(punit2,1130)
1130	format('* observation data')
	do iobs=1,nobs
1140	  iline=iline+1
	  read(punit1,'(a)',end=9000) cline
	  if(cline.eq.' ') go to 1140
	end do
	do iobs=1,nobs
	  i=iobsgp(iobs)
	  if(aa(i).eq.'n')then
	    write(punit2,1160) trim(aobs(iobs)),oval(iobs),weight(iobs),    &
	    aobsgp(iobsgp(iobs))
1160	    format(1x,a,t23,2x,g14.7,2x,g14.7,2x,a)
	  else
	    if(ag(i).eq.'n') then
	      if(aw(i,1).eq.'y')then
	        if((w3(i,1).lt.0.0).and.(oval(iobs).eq.0.0)) then
	          rtemp=1.0e36
	        else
	          rtemp=abs(oval(iobs))**w3(i,1)
	        end if
	        weight(iobs)=w1(i,1)*rtemp+w4(i,1)
	        if(weight(iobs).gt.w5(i,1)) weight(iobs)=w5(i,1)
	        if(weight(iobs).lt.w6(i,1)) weight(iobs)=w6(i,1)
	      end if
	      write(punit2,1160) trim(aobs(iobs)),oval(iobs),weight(iobs),  &
	      aobsgp(iobsgp(iobs))
	    else
	      i=iobsgp(iobs)
	      atemp=aobs(iobs)(1:nn(i))
	      atemp=adjustl(atemp)
	      call casetrans(atemp,'lo')
	      do j=1,nsubgp(i)
	        if(atemp.eq.asubgp(i,j))then
	          if(aw(i,j).eq.'y')then
	            if((w3(i,j).lt.0.0).and.(oval(iobs).eq.0.0)) then
	              rtemp=1.0e36
	            else
	              rtemp=abs(oval(iobs))**w3(i,j)
	            end if
	            weight(iobs)=w1(i,j)*rtemp+w4(i,j)
	            if(weight(iobs).gt.w5(i,j)) weight(iobs)=w5(i,j)
	            if(weight(iobs).lt.w6(i,j)) weight(iobs)=w6(i,j)
	          end if	       
	          write(punit2,1160)trim(aobs(iobs)),oval(iobs),weight(iobs),&
	          aasubgp(i,j)
	          go to 1180
	        end if
	      end do
	      write(amessage,1170)
1170	      format(' Error in program - contact programmer.')
1180	      continue
	    end if
	  end if
	end do
	do
	  read(punit1,'(a)',end=2000) cline
	  write(punit2,'(a)') trim(cline)
	end do

2000	write(amessage,2050) trim(pestfle2)
2050	format(' - file ',a,' written ok.')
	call write_message
	go to 9900


9000	write(amessage,9010) trim(pestfle1)
9010	format(' Unexpected end encountered to file ',a)
	go to 9890
9050	write(amessage,9060) trim(pestfle1)
9060	format(' Cannot find "* observation groups" section of PEST control ' &
	'file ',a)
	go to 9890
9100	write(amessage,9110)
9110	format(' Insufficient memory to continue execution.')
	go to 9890
9150	write(amessage,9160)
9160	format(' Cannot find "* observation data" section of PEST control ',  &
	'file',a)
	go to 9890
9200	write(amessage,9250)
9250	format(' Memory management error; cannot continue execution.')
	go to 9890

	


9890	call write_message(leadspace='yes')
9900	call close_files
	deallocate(aobs,iobsgp,weight,oval,stat=ierr)
	write(6,*)

end program adjobs
 