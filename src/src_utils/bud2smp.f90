!     Last change:  JD   12 Apr 2003    2:16 pm
program bud2smp

! -- Program BUD2SMP writes a bore sample file of model-generated system
!    inflows/outflows within user-specified zones.

	use defn
	use inter

	implicit none

        integer, parameter              :: MAXZONE=100
	integer		:: ifail,idate,ncol,nrow,itemp,ierr,budunit,        &
	                   dds,mms,yys,hhhs,mmms,ssss,ilay,                 &
	                   nlay,numzone,irow,icol,ltemp,izone,outunit,      &
	                   mcol,mrow,mlay,kstp,kper,itype,nlist,icrl,i,     &
	                   itime,recunit,iflag,maxtim,kstpold,kperold,narr,j,iheader
	integer         :: naux,ii
	integer, allocatable		:: jzone(:)
	integer, allocatable            :: iarray(:,:,:),intarray(:,:)
	real                            :: fac,totim,timwrit,delt,pertim,q
	real                            :: rtemp
	real, allocatable		:: flowrate(:,:),tw(:)
	real, allocatable               :: budarray(:,:,:)
	character (len=1)		:: am,ap,at
	character (len=4)		:: alay,atemp
	character (len=5)               :: anum
	character (len=16)		:: aflow,adate,atime,text
	character (len=16)              :: auxtxt(5)
	character (len=10)              :: ctemp
	character (len=80)		:: budfle,aprompt,outfle,recfle
	character (len=10), allocatable :: cid(:)
	type(modelgrid)                 :: gridspec

	iflag=0

	allocate(jzone(MAXZONE),cid(MAXZONE),stat=ierr)
	if(ierr.ne.0) go to 9500
	jzone=0

	write(amessage,5)
5       format(' Program BUD2SMP writes a bore sample file of ',             &
	'MODFLOW-generated inflows/outflows within ',              &
	'user-specified zones.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'

! -- Read settings file.

        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          go to 9890
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          go to 9890
        endif
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
          'settings.fig')
          go to 9890
        end if
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  go to 9890
	end if




! -- Read specification file (if present).

	call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow

30	write(6,35,advance='no')
35	format(' How many layers in model? ')
	itemp=key_read(nlay)
	if(escset.ne.0)then
	  escset=0
	  call free_grid_mem(gridspec)
	  write(6,*)
	  go to 10
	else if(itemp.ne.0) then
	  go to 30
	else if(nlay.le.0) then
	  go to 30
	end if

	allocate(budarray(ncol,nrow,nlay),iarray(ncol,nrow,nlay),          &
	intarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) go to 9500
	budarray=0
	iarray=0
	intarray=0

50	write(6,*)
	call open_input_file(ifail,                                        &
	' Enter name of MODFLOW unformatted budget output file: ',         &
	budfle,budunit,file_format='unformatted')
	if(ifail.ne.0) then
	  write(amessage,60) trim(budfle)
60	  format(' Cannot open unformatted file ',a)
	  go to 9890
	end if
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  deallocate(budarray,iarray,intarray,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 30
	end if

65	write(6,66,advance='no')
66	format(' Is this a MODFLOW88 or MODFLOW96 budget file  [8/9]? ',$)
	read(5,'(a)') am
	if(am.eq.' ') go to 65
	if(index(eschar,am).ne.0) then
	  close(unit=budunit)
	  go to 50
	end if
	if((am.ne.'8').and.(am.ne.'9')) go to 65
	if(am.eq.'8')then
	  write(amessage,72)
72	  format(' BUD2SMP cannot continue execution as it can read ',   &
	  'a compact MODFLOW96 budget file only.')
	  go to 9890
	end if

39	write(6,40,advance='no')
40	format(' Enter maximum number of output times: ')
	itemp=key_read(maxtim)
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 65
	end if
	if(itemp.ne.0) go to 39
	if(maxtim.le.0) go to 39
	maxtim=maxtim+1
	allocate(tw(maxtim),flowrate(maxtim,MAXZONE),stat=ierr)
	flowrate=0.0
	if(ierr.ne.0) go to 9500

70	write(6,80,advance='no')
80	format(' Enter text to identify MODFLOW flow type: ')
	read(5,'(a)') aflow
	if(aflow.eq.' ') go to 70
	if(index(eschar,aflow(1:2)).ne.0)then
	  write(6,*)
	  deallocate(tw,flowrate,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 39
	end if
	call casetrans(aflow,'hi')

100	write(6,*)
110	if(datespec.eq.1) then
	  write(6,120,advance='no')
120	  format(' Enter simulation starting date [dd/mm/yyyy]: ')
	else
	  write(6,121,advance='no')
121	  format(' Enter simulation starting date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') adate
	if(adate.eq.' ') go to 110
	adate=adjustl(adate)
	if(index(eschar,adate(1:2)).ne.0) then
	  write(6,*)
	  go to 70
	end if
	call char2date(ifail,adate,dds,mms,yys)
	if(ifail.ne.0)then
	  write(6,130)
130	  format(' Illegal date  - try again.')
	  go to 110
	end if

150	write(6,160,advance='no')
160	format(' Enter simulation starting time [hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 150
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) go to 100
	call char2time(ifail,atime,hhhs,mmms,ssss)
	if(ifail.ne.0)then
	  write(6,180)
180	  format(' Illegal time  - try again.')
	  go to 150
	end if

190	write(6,192,advance='no')
192	format(' Enter time units employed by model [y/d/h/m/s]: ')
	read(5,'(a)') at
	if(at.eq.' ') go to 190
	if(index(eschar,at).ne.0)then
	  write(6,*)
	  go to 150
	end if
	if((at.eq.'Y').or.(at.eq.'y')) then
	  at='y'
	else if((at.eq.'D').or.(at.eq.'d')) then
	  at='d'
	else if ((at.eq.'H').or.(at.eq.'h')) then
	  at='h'
	else if((at.eq.'M').or.(at.eq.'m')) then
	  at='m'
	else if((at.eq.'S').or.(at.eq.'s')) then
	  at='s'
	else
	  go to 190
	end if

! -- The integer arrays comprising the 3-d zonation are next read.

200	write(6,*)
	ilay=0
220	ilay=ilay+1
	if(ilay.gt.nlay) go to 235
230	write(alay,'(i4)') ilay
	alay=adjustl(alay)
	aprompt=' Enter name of integer array file for layer '//  &
	trim(alay)//': '
	call read_integer_array(ifail,aprompt,iarray(:,:,ilay),  &
	pm_header=headerspec,rows=nrow,columns=ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1)then
	  escset=0
	  write(6,*)
	  if(ilay.eq.1)then
	    go to 190
	  else
	    ilay=ilay-1
	    go to 230
	  end if
	end if
	go to 220
235	continue

! -- The number of zones in the model domain is next ascertained.

	numzone=0
	do ilay=1,nlay
	  do irow=1,nrow
	    do icol=1,ncol
	      if(iarray(icol,irow,ilay).ne.0)then
	        if(numzone.eq.0)then
	          numzone=numzone+1
	          jzone(numzone)=iarray(icol,irow,ilay)
	          ltemp=jzone(numzone)
	        else
	          itemp=iarray(icol,irow,ilay)
	          if(itemp.eq.ltemp) go to 240
	          ltemp=itemp
	          do izone=1,numzone
	            if(itemp.eq.jzone(izone)) go to 240
	          end do
	          numzone=numzone+1
	          if(numzone.gt.MAXZONE)then
	            write(amessage,250)
250	            format(' Increase MAXZONE and re-compile program.')
	            go to 9890
	          end if
	          jzone(numzone)=itemp
240	          continue
	        end if
	      end if
	    end do
	  end do
	end do
	write(atemp,'(i4)') numzone
	atemp=adjustl(atemp)
	write(amessage,260) trim(atemp)
260	format(' A total of ',a,' different non-zero zones were ',           &
	'identified in integer arrays.')
	call write_message(leadspace='yes')
	write(amessage,261)
261	format(' An identifier must now be provided for each zone to ',      &
	'appear in the bore sample output file:-')
	call write_message

	i=0
268	i=i+1
	if(i.gt.numzone) go to 290
269	call num2char(jzone(i),anum)
	write(6,270,advance='no') trim(anum)
270	format('   Enter identifier for flows in zone ', a,                  &
	' (10 characters or less): ')
	read(5,'(a)') ctemp
	if(ctemp.eq.' ') go to 269
	if(index(eschar,ctemp(1:2)).ne.0) then
	  if(i.eq.1)then
	    write(6,*)
	    ilay=nlay
	    go to 230
	  else
	    write(6,*)
	    i=i-1
	    go to 269
	  end if
	end if
	cid(i)=ctemp
	cid(i)=adjustl(cid(i))
	call casetrans(cid(i),'hi')
	if(i.ne.1)then
	  do j=1,i-1
	    if(cid(i).eq.cid(j))then
	      write(amessage,287)
287	      format(' This identifier has already been used - try again.')
	      call write_message
	      go to 269
	    end if
	  end do
	end if
	go to 268
290	continue

300	write(6,*)
310	call open_output_file(ifail, &
	' Enter name for bore sample output file: ',outfle,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  i=numzone
	  go to 269
	end if

350	write(6,360,advance='no')
360	format(' Enter flow rate factor: ')
	i=key_read(fac)
	if(escset.ne.0)then
	  escset=0
	  close(unit=outunit)
	  go to 300
	else if(i.eq.-1) then
	  go to 350
	else if(i.ne.0) then
	  go to 350
	end if

370	write(6,380,advance='no')
380	format(' Assign flows to beginning, middle ',                    &
        'or finish of time step?  [b/m/f]: ')
	read(5,'(a)') ap
	if(ap.eq.' ') go to 370
	if(index(eschar,ap).ne.0) then
	  write(6,*)
	  go to 350
	end if
	if((ap.eq.'B').or.(ap.eq.'b'))then
	  ap='b'
	else if((ap.eq.'M').or.(ap.eq.'m'))then
	  ap='m'
	else if((ap.eq.'F').or.(ap.eq.'f'))then
	  ap='f'
	else
	  go to 370
	end if

390	write(6,*)
395	call open_output_file(ifail, &
	' Enter name for run record file: ',recfle,recunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 370
	end if
	write(recunit,397)
397	format(' Stress_period  Time_step  Elapsed_time    Flow_type     ',   &
	'     Flow_processed_by_BUD2SMP')

! -- A MODFLOW96 budget file is read.

800	itime=1
	kperold=0
	kstpold=0
	read_an_array: do
	  budarray=0
	  read(budunit,end=900,err=9000) kstp,kper,text,mcol,mrow,mlay
!          write(*,*) kstp,kper,text,mcol,mrow,ilay   !debug
	  if((kstp.ne.kstpold).or.(kper.ne.kperold))then
	    narr=0
	    kstpold=kstp
	    kperold=kper
	  end if
	  if(mlay.gt.0)then
!	    write(6,*) kstp,kper,text,mcol,mrow,mlay     !debug
	    write(amessage,810) trim(budfle)
810	    format(' Cell-by-cell flow terms in MODFLOW budget file ',     &
	    a,' are not stored in compressed form.')
	    call write_message(leadspace='yes')
	    write(amessage,820)
820	    format(' Consult the MODFLOW96 manual for settings required '  &
	    'to ensure compressed flow term storage.')
	    go to 9890
	  end if
	  mlay=-mlay
!          write(*,*) mcol,mrow,ncol,nrow  !debug
	  if((mcol.ne.ncol).or.(mrow.ne.nrow)) go to 9200
	  if(mlay.ne.nlay) go to 9250
	  read(budunit,end=9100,err=9000) itype,delt,pertim,totim
!          write(*,*) itype,delt,pertim,totim   !debug
	  if(itype.eq.1)then
	    call casetrans(text,'hi')
	    if(index(text,trim(aflow)).eq.0)then
	      read(budunit,err=9050,end=9100) budarray
	      write(recunit,830) kper,kstp,totim,adjustl(text),'no'
830	      format(3x,i10,5x,i5,4x,1pg13.5,3x,a,t70,a)
	      cycle read_an_array
	    end if
	    if(narr.eq.1)then
	      write(amessage,815) trim(aflow),trim(budfle)
815	      format(' Text string "',a,'" does not identify a unique ',    &
	      'flow type recorded in file ',a)
	      call write_message
	      go to 9890
	    end if
	    narr=narr+1
	    read(budunit,err=9050,end=9100) budarray
	    do izone=1,numzone
	      flowrate(itime,izone)=0.0
	      itemp=jzone(izone)
	      do ilay=1,nlay
	        do irow=1,nrow
	          do icol=1,ncol
	            if(iarray(icol,irow,ilay).eq.itemp)then
	              flowrate(itime,izone)=flowrate(itime,izone)+          &
	              budarray(icol,irow,ilay)
	            end if
	          end do
	        end do
	      end do
	    end do
	  else if (itype.eq.2)then
	    read(budunit,err=9000,end=9100) nlist
	    call casetrans(text,'hi')
	    if(index(text,trim(aflow)).eq.0)then
	      do i=1,nlist
	        read(budunit,err=9050,end=9100) icrl,q
	      end do
	      write(recunit,830) kper,kstp,totim,adjustl(text),'no'
	      cycle read_an_array
	    end if
	    if(narr.eq.1)then
	      write(amessage,815) trim(aflow),trim(budfle)
	      call write_message
	      go to 9890
	    end if
	    narr=narr+1
	    do i=1,nlist
	      read(budunit,err=9050,end=9100) icrl,q
!	      write(*,*) icrl
	      ilay=(icrl-1)/(ncol*nrow)+1
	      icrl=icrl-(ilay-1)*nrow*ncol
	      irow=(icrl-1)/ncol+1
	      icrl=icrl-(irow-1)*ncol
	      icol=icrl
!	      write(6,*) icol,irow,ilay
	      itemp=iarray(icol,irow,ilay)
	      do izone=1,numzone
	        if(jzone(izone).eq.itemp) then
	          flowrate(itime,izone)=flowrate(itime,izone)+q
	          go to 780
	        end if
	      end do
780	      continue
	    end do
	  else if(itype.eq.4) then
	    call casetrans(text,'hi')
	    if(index(text,trim(aflow)).eq.0)then
	      read(budunit,err=9050,end=9100)                         &
              ((budarray(icol,irow,1),icol=1,ncol),irow=1,nrow)
	      write(recunit,830) kper,kstp,totim,adjustl(text),'no'
	      cycle read_an_array
	    end if
	    if(narr.eq.1)then
	      write(amessage,815) trim(aflow),trim(budfle)
	      call write_message
	      go to 9890
	    end if
	    narr=narr+1
	    read(budunit,err=9050,end=9100)                           &
            ((budarray(icol,irow,1),icol=1,ncol),irow=1,nrow)
	    do izone=1,numzone
	      flowrate(itime,izone)=0.0
	      itemp=jzone(izone)
	      do irow=1,nrow
	        do icol=1,ncol
	          if(iarray(icol,irow,1).eq.itemp)then
	            flowrate(itime,izone)=flowrate(itime,izone)+      &
	            budarray(icol,irow,1)
	          end if
	        end do
	      end do
	    end do
	  else if(itype.eq.3)then
	    call casetrans(text,'hi')
	    if(index(text,trim(aflow)).eq.0)then
	      read(budunit,err=9050,end=9100)                         &
              ((intarray(icol,irow),icol=1,ncol),irow=1,nrow)
	      read(budunit,err=9050,end=9100)                         &
              ((budarray(icol,irow,1),icol=1,ncol),irow=1,nrow)
	      write(recunit,830) kper,kstp,totim,adjustl(text),'no'
	      cycle read_an_array
	    end if
	    if(narr.eq.1)then
	      write(amessage,815) trim(aflow),trim(budfle)
	      call write_message
	      go to 9890
	    end if
	    narr=narr+1
	    read(budunit,err=9050,end=9100)                             &
	    ((intarray(icol,irow),icol=1,ncol),irow=1,nrow)
	    read(budunit,err=9050,end=9100)                             &
	    ((budarray(icol,irow,1),icol=1,ncol),irow=1,nrow)
	    do irow=1,nrow
	      do icol=1,ncol
	        itemp=intarray(icol,irow)
	        if((itemp.lt.1).or.(itemp.gt.nlay))then
	          iflag=1
	        else
	          if(itemp.ne.1)then
	            budarray(icol,irow,itemp)=budarray(icol,irow,1)
	            budarray(icol,irow,1)=0.0
	          end if
	        endif
	      end do
	    end do

	    do izone=1,numzone
	      flowrate(itime,izone)=0.0
	      itemp=jzone(izone)
	      do ilay=1,nlay
	        do irow=1,nrow
	          do icol=1,ncol
	            if(iarray(icol,irow,ilay).eq.itemp)then
	              flowrate(itime,izone)=flowrate(itime,izone)+         &
	              budarray(icol,irow,ilay)
	            end if
	          end do
	        end do
	      end do
	    end do
          else if (itype.eq.5) then
            read(budunit,err=9000,end=9100) naux
            naux=naux-1
!            write(6,*) ' naux = ',naux   !debug
            if(naux.gt.0)then
              read(budunit,err=9000,end=9100) (auxtxt(i),i=1,naux)
            end if
 	    read(budunit,err=9000,end=9100) nlist
	    call casetrans(text,'hi')
	    if(index(text,trim(aflow)).eq.0)then
	      if(naux.eq.0)then
                do i=1,nlist
                  read(budunit,err=9050,end=9100)  icrl,q
                end do
              else
                do i=1,nlist
                  read(budunit,err=9050,end=9100)  icrl,q,(rtemp,ii=1,naux)
                end do
              end if
	      write(recunit,830) kper,kstp,totim,adjustl(text),'no'
	      cycle read_an_array
	    end if
	    if(narr.eq.1)then
	      write(amessage,815) trim(aflow),trim(budfle)
	      call write_message
	      go to 9890
	    end if
	    narr=narr+1
	    do i=1,nlist
	      if(naux.eq.0)then
                read(budunit,err=9050,end=9100) icrl,q
              else
                read(budunit,err=9050,end=9100) icrl,q,(rtemp,ii=1,naux)
              end if
!	      write(*,*) icrl                   !debug
	      ilay=(icrl-1)/(ncol*nrow)+1
	      icrl=icrl-(ilay-1)*nrow*ncol
	      irow=(icrl-1)/ncol+1
	      icrl=icrl-(irow-1)*ncol
	      icol=icrl
!	      write(6,*) icol,irow,ilay         !debug
	      itemp=iarray(icol,irow,ilay)
	      do izone=1,numzone
	        if(jzone(izone).eq.itemp) then
	          flowrate(itime,izone)=flowrate(itime,izone)+q
	          go to 777
	        end if
	      end do
777	      continue
	    end do
          else
            text=adjustl(text)
            write(amessage,782) trim(text),itype
782         format(' Unknown compact storage method: data type = ',a, &
            '; storage method index =',i3)
            go to 9890
	  end if
	  if(ap.eq.'b')then
	    timwrit=totim-delt
	  else if(ap.eq.'m')then
	    timwrit=totim-delt/2.0
	  else
	    timwrit=totim
	  end if
	  tw(itime)=timwrit
	  call datestring(dds,mms,yys,hhhs,mmms,ssss,timwrit,               &
          at,adate,atime)
	  text=adjustl(text)
	  write(recunit,830) kper,kstp,totim,adjustl(text),'yes'
	  do i=1,len_trim(text)
	    if(text(i:i).eq.' ') text(i:i)='_'
	  end do
	  itime=itime+1
	  if(itime.gt.maxtim) then
	    call num2char(maxtim-1,anum)
	    write(amessage,440) trim(anum),trim(budfle)
440	    format(' There are greater than ',a,' output times represented ', &
	    'in file ',a,'. Run BUD2SMP again and supply a larger value ',    &
	    'for the maximum number of output times.')
	    go to 9890
	  end if
	end do read_an_array
	go to 900

! -- Write the bore sample file.

900	continue
	itime=itime-1
	do izone=1,numzone
	  do i=1,itime
	    call datestring(dds,mms,yys,hhhs,mmms,ssss,tw(i),               &
            at,adate,atime)
	    write(outunit,905) trim(cid(izone)),trim(adate),trim(atime),    &
	    flowrate(i,izone)*fac
905	    format(1x,a,t14,a,t27,a,t40,1pg14.7)
	  end do
	end do

! -- Tidy up

	call num2char(itime,atime)
	write(amessage,910) trim(atime),trim(outfle)
910	format(' - data for ',a,' model output arrays written ',              &
	'to file ',a)
	call write_message(leadspace='yes')
	write(amessage,920) trim(recfle),trim(budfle)
920	format(' - see file ',a,' for a record of arrays found in file ',a)
	call write_message

	if(iflag.ne.0)then
	  write(amessage,940)
940	  format(' Warning: there is at least one recharge or evt ',         &
	  'indicator array in the MODFLOW input dataset for which ',         &
	  'the layer reference is out of bounds.')
	  call write_message(leadspace='yes')
	end if

	go to 9900



9000	write(amessage,9010) trim(budfle)
9010	format(' Error reading header to budget array in file ',a)
	go to 9890
9050	write(amessage,9060) trim(budfle)
9060	format(' Error reading an array from budget file ',a)
	go to 9890
9100	write(amessage,9110) trim(budfle)
9110	format(' Unexpected end encountered to budget file ',a)
	go to 9890
9200	write(amessage,9210) trim(budfle)
9210	format(' Dimensions of grid as provided in file ',a,           &
	' do not agree with those in grid specification file.')
	go to 9890
9250	write(amessage,9260)trim(budfle)
9260	format(' Number of model layers as provided in file ',a,     &
	' does not agree with that supplied by user.')
	go to 9890
9300	write(amessage,9310)
9310	format(' Memory management error: cannot continue execution.')
	go to 9890
9500	write(amessage,9510)
9510	format(' Insufficient memory to continue execution.')
	go to 9890

9890	call write_message(leadspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(budarray,iarray,flowrate,tw,cid,intarray,stat=ierr)
	write(6,*)

end program bud2smp
