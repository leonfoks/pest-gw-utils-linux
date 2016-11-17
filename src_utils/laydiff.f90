!     Last change:  JD   27 Sep 2003    6:11 am
program laydiff

! -- Program LAYDIFF writes a table of inter-layer head differences.

	use defn
	use inter

	implicit none

	integer      :: i,ifail,j,ierr,sampunit,dds,mms,yys,hhhs,mmms,ssss, &
                        ndayref,nsecref,outunit,iactive,iline,idone,cols,ndays,nsecs, &
                        ndays1,nsecs1,ndays2,nsecs2,ilist,iheader,idate,laymin,laymax, &
                        icount,ilay,jlow,jlist,jcount,kcount
	integer, allocatable, dimension(:)         :: layer
	real, allocatable, dimension(:)            :: head
        double precision                           :: value,value1,value2,dayfactor, &
                                                      exdist,dd1,dd2,mindist,dist,exday
	double precision, allocatable, dimension(:):: east,north
        character (len=1)   :: aa
	character (len=10)  :: abore, aboreold
        character (len=22)  :: aobs
	character (len=200) :: sampfile,outfile
	character (len=15)  :: adate,atime,atemp
	type (modelgrid) gridspec


!Initialisation

        dayfactor=1.0d0/86400.0d0

	write(amessage,5)
5       format(' Program LAYDIFF calculates layer data differences based on ', &
        'data contained in a bore sample file.')
	call write_message(leadspace='yes',endspace='yes')

	call readfig(gridspec%specfile,bore_coord_file,sampfile)
        call read_settings(ifail,idate,iheader)
        if(ifail.eq.1) then
          write(amessage,7)
7         format(' A settings file (settings.fig) was not found in the ', &
          'current directory.')
          call write_message
          go to 9900
        else if(ifail.eq.2) then
          write(amessage,8)
8         format(' Error encountered while reading settings file settings.fig')
          call write_message
          go to 9900
        endif
        if((idate.ne.0).or.(datespec.eq.0)) then
          write(amessage,9)
9         format(' Cannot read date format from settings file ', &
          'settings.fig')
          call write_message
          go to 9900
        end if

30      call read_bore_coord_file(ifail, &
	' Enter name of bore coordinates file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 9900
	end if
	do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) go to 80
	end do
	go to 100
80	write(amessage,90) trim(bore_coord_file)
90	format(' The following bores were not provided with layer numbers in ',&
	'bore coordinates file ',a,' ----->')
	call write_message(leadspace='yes')
	j=2
	imessage=0
	write_bores: do i=1,num_bore_coord
	  if(bore_coord_layer(i).eq.-999) then
	    write(amessage(j:),'(a)') trim(bore_coord_id(i))
	    j=j+11
	    if(j.ge.71) then
	      call write_message(increment=1)
	      if(imessage.gt.12) goto 9900
	      j=2
	    end if
	  end if
	end do write_bores
	if(j.ne.2) call write_message
	go to 9900

100     call read_bore_list_file(ifail, &
       ' Enter name of bore listing file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 30
	end if

	if(num_bore_list.gt.1) then
	  do i=1,num_bore_list-1
	    do j=i+1,num_bore_list
	      if(bore_list_id(i).eq.bore_list_id(j))then
	        write(amessage,110) trim(bore_list_file)
110		format(' Execution of program LAYDIFF cannot proceed as ',&
		'there are multiple occurrences of the same bore in bore ',&
		'listing file ',a)
		go to 9890
	      end if
	    end do
	  end do
	end if

	allocate(east(num_bore_list), north(num_bore_list), &
	layer(num_bore_list),head(num_bore_list),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,130)
130	  format(' Insufficient memory available to continue execution.')
	  go to 9890
	end if
        head=-1.0e35             ! head is an array

	imessage=0
	list: do i=1,num_bore_list
	  do j=1,num_bore_coord
	    if(bore_coord_id(j).eq.bore_list_id(i)) then
	      east(i)=bore_coord_east(j)
	      north(i)=bore_coord_north(j)
	      layer(i)=bore_coord_layer(j)
	      cycle list
	    end if
	  end do
	  write(amessage,140) trim(bore_list_id(i)),trim(bore_list_file),&
	  trim(bore_coord_file)
140	  format(' No coordinates for bore ',a,' from bore listing ',&
	  'file ',a,' are provided in bore coordinates file ',a)
	  if(imessage.eq.0) write(6,*)
	  call write_message(increment=1)
	end do list
	if(imessage.ne.0) go to 9900

! -- We now find the minimum and maximum layers of listed bores.

        laymin=1000000
        laymax=-1000000
        do ilist=1,num_bore_list
          if(layer(ilist).lt.laymin)laymin=layer(ilist)
          if(layer(ilist).gt.laymax)laymax=layer(ilist)
        end do

600	call open_named_input_file(ifail,          &
	' Enter name of bore sample file: ',sampfile,sampunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  deallocate(east,north,layer,head,stat=ierr)
	  if(ierr.ne.0) go to 9300
	  go to 100
	end if

200	write(6,*)
210	if(datespec.eq.1) then
	  write(6,220,advance='no')
220	  format(' Enter reference date [dd/mm/yyyy]: ')
	else
	  write(6,221,advance='no')
221	  format(' Enter reference date [mm/dd/yyyy]: ')
	end if
	read(5,'(a)') adate
	if(adate.eq.' ') go to 210
	adate=adjustl(adate)
	if(index(eschar,adate(1:2)).ne.0) then
	  write(6,*)
          close(unit=sampunit)
	  go to 600
	end if
	call char2date(ifail,adate,dds,mms,yys)
	if(ifail.ne.0)then
	  write(6,240)
240	  format(/,' Illegal date  - try again.',/)
	  go to 210
	end if

260	write(6,270,advance='no')
270	format(' Enter reference time [hh:mm:ss]: ')
	read(5,'(a)') atime
	if(atime.eq.' ') go to 260
	atime=adjustl(atime)
	if(index(eschar,atime(1:2)).ne.0) go to 200
	call char2time(ifail,atime,hhhs,mmms,ssss)
	if(ifail.ne.0)then
	  write(6,280)
280	  format(/,' Illegal time  - try again.',/)
	  go to 260
	end if
	ndayref=numdays(1,1,1970,dds,mms,yys)
	nsecref=numsecs(0,0,0,hhhs,mmms,ssss)

290	write(6,*)
295	write(6,296,advance='no')
296	format(' Enter max days to reference date (fractional if necessary): ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 295
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  write(6,*)
	  go to 260
	end if
	call char2num(ifail,atemp,exday)
	if(ifail.ne.0) go to 295
	if(exday.le.0.0) go to 295

310	continue
320	write(6,330,advance='no')
330	format(' Enter exclusion distance: ')
	read(5,'(a)') atemp
	if(atemp.eq.' ') go to 320
	atemp=adjustl(atemp)
	if(index(eschar,atemp(1:2)).ne.0) then
	  go to 290
	end if
	call char2num(ifail,atemp,exdist)
	if(ifail.ne.0) go to 310
	if(exdist.le.0.0) go to 310
        exdist=exdist*exdist

400	write(6,*)
410	call open_output_file(ifail, &
	' Enter name for output data file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
          write(6,*)
	  go to 310
	end if

! -- We now assign a value to each well in the listing file based on data in the
!    bore sample file (if a date is close enough).

       idone=0
       iline=0
       aboreold=' '
       iactive=0
       do
         iline=iline+1
         read(sampunit,'(a)',end=500) cline
         if(cline.eq.' ') cycle
         call linesplit(ifail,5)
         if(ifail.eq.0)then
           cols=5
         else
           call linesplit(ifail,4)
           if(ifail.ne.0) go to 9150
           cols=4
         end if
         abore=cline(left_word(1):right_word(1))
         call casetrans(abore,'hi')
435      continue
         if(abore.eq.aboreold)then
           if(iactive.eq.0) cycle
           if(idone.eq.1)cycle
           call read_rest_of_sample_line(ifail,cols,ndays,nsecs,value,iline,sampfile)
           if(ifail.ne.0) go to 9900
           if(value.gt.-1.0d30)then
             if((ndays.lt.ndayref).or.    &
                ((ndays.eq.ndayref).and.(nsecs.le.nsecref)))then
                  ndays1=ndays
                  nsecs1=nsecs
                  value1=value
                  idone=-1
             else
               ndays2=ndays
               nsecs2=nsecs
               value2=value
               idone=1
               dd2=ndays2-ndayref+float(nsecs2-nsecref)*dayfactor
               dd1=ndayref-ndays1+float(nsecref-nsecs1)*dayfactor
               if((dd1.gt.exday).and.(dd2.gt.exday))then
                 continue
               else if((value1.lt.-1.0e30).and.(value2.lt.-1.0e30))then
                 continue
               else if(value1.lt.-1.0e30)then
                 if(dd2.le.exday)head(ilist)=value2
               else if(value2.lt.-1.0e30)then
                 if(dd1.le.exday)head(ilist)=value1
               else
                 head(ilist)=value1+dd1/(dd1+dd2)*(value2-value1)
               end if
             end if
           end if
         else
           if(idone.eq.-1)then
             dd1=ndayref-ndays1+float(nsecref-nsecs1)*dayfactor
             if(dd1.le.exday)then
               head(ilist)=value1
             end if
           end if
           idone=0
           aboreold=abore
           do ilist=1,num_bore_list
             if(abore.eq.bore_list_id(ilist)) go to 450
           end do
           iactive=0
           cycle
450        iactive=1
           value1=-1.0d35
           value2=-1.0d35
           go to 435
         end if
       end do
500    continue
       if(idone.eq.-1)then
         dd1=ndayref-ndays1+float(nsecref-nsecs1)*dayfactor
         if(dd1.le.exday)then
           head(ilist)=value1
         end if
       end if
       close(unit=sampunit)
       write(*,520) trim(sampfile)
520    format(' - file ',a,' read ok.')

! -- Now for each well in the listing file - where possible - we evaluate
!    head differences.

       icount=0 
       write(outunit,530)
530    format('Well_1',t12,'Well_2',t24,'lay_1',t31,'lay_2',    &
       t40,'Easting_1',t55,'Northing_1',t70,'Easting_2',t85,'Northing_2',t98,  &
       'Value_1',t113,'Value_2',t128,'Difference')
       do ilist=1,num_bore_list
         if(head(ilist).lt.-1.0e30)cycle
         do ilay=laymin,laymax
           if(ilay.le.layer(ilist))cycle   ! We are only linking to layers below
           mindist=1.0d300
           do jlist=1,num_bore_list
             if(jlist.eq.ilist)cycle
             if(layer(jlist).ne.ilay)cycle
             if(head(jlist).lt.-1.0e30) cycle
             dist=(east(ilist)-east(jlist))*(east(ilist)-east(jlist))+   &
                  (north(ilist)-north(jlist))*(north(ilist)-north(jlist))
             if(dist.lt.mindist)then
               mindist=dist
               jlow=jlist
             end if
           end do
           if(mindist.le.exdist)then
             write(outunit,590) trim(bore_list_id(ilist)),trim(bore_list_id(jlow)),  &
             layer(ilist),layer(jlow),east(ilist),north(ilist),east(jlow),north(jlow), &
             head(ilist),head(jlow),head(ilist)-head(jlow)
590          format(a,t12,a,t24,i3,t31,i2,t38,f12.3,t53,f12.3,t68,f12.3,t83,f12.3, &
             t98,1pg13.6,t113,1pg13.6,t128,1pg13.6)
             icount=icount+1
           end if
         end do
       end do
       if(icount.gt.0)then
         write(6,620) trim(outfile)
620      format(' - file ',a,' written ok.')
       else
         write(amessage,630)
630      format(' Unable to write inter-layer differences to output file because ', &
         'no differences could be calculated on basis of provided sample time and ', &
         'distance windows.')
         call write_message(leadspace='yes',endspace='yes')
       end if
       close(unit=outunit)

       if(icount.eq.0)go to 9900

! -- An instruction file is now generated to read the LAYDIFF output file.

       write(6,*)
650    write(6,660,advance='no')
660    format(' Generate an instruction file to read output file?  [y/n]: ')
       read(5,'(a)') aa
       if(aa.eq.' ') go to 650
       call casetrans(aa,'lo')
       if((aa.ne.'y').and.(aa.ne.'n')) go to 650
       if(aa.eq.'n') go to 9900
680    call open_output_file(ifail, &
       ' Enter name for this instruction file: ',outfile,outunit)
       if(ifail.ne.0) go to 9900
       if(escset.ne.0) then
         escset=0
         write(6,*)
         go to 650
       end if

! -- The following code is a repeat of the above. Any changes to the above must be made
!    below as well.

       write(outunit,685)
685    format('pif $')
       write(outunit,684)
684    format('l1')
       jcount=0
       kcount=0
       do ilist=1,num_bore_list
         if(head(ilist).lt.-1.0e30)cycle
         do ilay=laymin,laymax
           if(ilay.le.layer(ilist))cycle   ! We are only linking to layers below
           mindist=1.0d300
           do jlist=1,num_bore_list
             if(jlist.eq.ilist)cycle
             if(layer(jlist).ne.ilay)cycle
             if(head(jlist).lt.-1.0e30) cycle
             dist=(east(ilist)-east(jlist))*(east(ilist)-east(jlist))+   &
                  (north(ilist)-north(jlist))*(north(ilist)-north(jlist))
             if(dist.lt.mindist)then
               mindist=dist
               jlow=jlist
             end if
           end do
           if(mindist.le.exdist)then
             aobs=trim(bore_list_id(ilist))//'-'//trim(bore_list_id(jlow))
             call casetrans(aobs,'lo')
             if(len_trim(aobs).gt.12)kcount=kcount+1
             write(outunit,686) trim(aobs)
686          format('l1  [',a,']128:141')
             jcount=jcount+1
           end if
         end do
       end do

       close(unit=outunit)
       if(jcount.ne.icount)then
         write(*,687)
687      format(/,' *** Programming error - contact WNC ***',/)
         stop
       end if
       write(6,688) trim(outfile)
688    format(' - instruction file ',a,' written ok.')
       if(kcount.ne.0)then
         write(amessage,689) trim(outfile)
689      format(' Note: at least one of the observation names recorded in file ',a,  &
         ' is greater than 12 characters in length. This will need to be altered ', &
         'manually.')
         call write_message(leadspace='yes')
       end if
       go to 9900

9150    call num2char(iline,atemp)
        write(amessage,9160) trim(atemp),trim(sampfile)
9160    format(' Error on line ',a,' of bore sample file ',a,': ',&
        'insufficient entries.')
        go to 9890
9300	write(amessage,9310)
9310	format(' Memory management error: cannot continue execution.')
	go to 9890
9500    call num2char(iline,atemp)
        write(amessage,9510) trim(atemp),trim(sampfile)
9510    format(' Error on line ',a,' of bore sample file ',a,': bore ',&
        'identifier greater than 10 characters in length.')
        go to 9890
9600    call num2char(iline,atemp)
        write(amessage,9610) trim(atemp),trim(sampfile)
9610    format(' Error reading line ',a,' of bore sample file ',a)
        go to 9890


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_bore_mem
	deallocate(layer,head,east,north,stat=ierr)
	write(6,*)

end program laydiff


