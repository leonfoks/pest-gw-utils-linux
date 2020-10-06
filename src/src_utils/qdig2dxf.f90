!     Last change:  JD   14 Feb 2001    9:00 am
program qdig2dxf

! -- Program QDIG2DXF produces a DXF file from a QDIGIT "line", "points" 
!    or "well" file.


	use defn
	use inter

	implicit none

	integer                                 :: ifail,ierr,iend, &
						   outunit,inunit,iline,icoord,ifail1, &
                                                   nbb
	real					:: single_num
	double precision			:: clength,easting,northing,&
						   eastmax,eastmin,northmax,&
						   northmin,delta,oldeast,&
						   oldnorth,lasteast,lastnorth
	character (len=80)			:: qfile,aprompt,outfile,bfile
	character (len=1)			:: atyp
	character (len=20)			:: anum1


	write(amessage,5)
5       format(' Program QDIG2DXF produces a DXF file from a QDIGIT ',&
	'"line", "points" or "well" file.')
	call write_message(leadspace='yes')

9	write(6,*)
10	imessage=0
20	write(6,30,advance='no')
30	format(' Enter name of QDIGIT-generated file: ')
	read(5,'(a)') qfile
	if(qfile.eq.' ') go to 20
	qfile=adjustl(qfile)
	if(index(eschar,qfile(1:2)).ne.0) go to 9900
        bfile=qfile
        nbb=len_trim(bfile)
        call getfile(ifail1,bfile,qfile,1,nbb)
        if(ifail1.ne.0) go to 20
	inunit=nextunit()
	open(unit=inunit,file=qfile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.gt.5) go to 9900
	  write(amessage,50) trim(qfile)
50	  format(' Cannot open file ',a,'  - try again.')
	  call write_message(increment=1)
	  go to 20
	end if

80	write(6,90,advance='no')
90	format(' Is this file a line, points or well file? [l/p/w]: ')
	read(5,'(a)') atyp
	if(atyp.eq.' ') go to 80
	if(index(eschar,atyp).ne.0) then
	  close(unit=inunit)
	  go to 9
	end if
	call casetrans(atyp,'lo')
	if((atyp.ne.'l').and.(atyp.ne.'p').and.(atyp.ne.'w')) go to 80
95	if(atyp.ne.'l')then
100	  write(6,110,advance='no')
110	  format(' Represent points by a cross of what length? ')
	  read(5,'(a)') anum1
	  if(anum1.eq.' ') go to 100
	  anum1=adjustl(anum1)
	  if(index(eschar,anum1(1:2)).ne.0) then
	    write(6,*)
	    go to 80
	  end if
	  call char2num(ifail,anum1,clength)
	  if(ifail.ne.0) go to 100
	  if(pos_test(clength, 'length of cross').ne.0) go to 100
	  clength=clength/2.0
	else
113	  write(6,115,advance='no')
115	  format(' Enter trimming distance: ')
	  read(5,'(a)') anum1
	  if(anum1.eq.' ') go to 113
	  anum1=adjustl(anum1)
	  if(index(eschar,anum1(1:2)).ne.0) then
	    write(6,*)
	    go to 80
	  end if
	  call char2num(ifail,anum1,clength)
	  if(ifail.ne.0) go to 113
	  if(pos_test(clength, 'trimming distance').ne.0) go to 113
	  clength=clength/2.0
	end if  

120	write(6,*)
	aprompt=' Enter name for DXF output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 95
	end if

! -- The QDIGIT file is first read in its entirety to determine geographic
!    bounds.

	iline=0
	icoord=0
	do
	  iline=iline+1
	  if(iline.eq.2000) then
	    write(6,160)
160	    format(/,' Working.....',/)
	  end if
	  read(inunit,'(a)',err=9000,end=200) cline
	  call linesplit(ifail,2)
	  if(ifail.lt.0) cycle
	  if(ifail.gt.0)then
	    if(atyp.eq.'l')then
	      if(cline(left_word(1):left_word(1)+6).ne.'ENDLINE') go to 9100
	      cycle
	    else
	      go to 9100
	    end if
	  end if
	  easting=char2double(ifail,1)
	  if(ifail.ne.0) go to 9100
	  northing=char2double(ifail,2)
	  if(ifail.ne.0) go to 9100
	  icoord=icoord+1
	  if(icoord.eq.1)then
	    eastmax=easting
	    eastmin=easting
	    northmax=northing
	    northmin=northing
	  else
	    eastmin=min(eastmin,easting)
	    eastmax=max(eastmax,easting)
	    northmin=min(northmin,northing)
	    northmax=max(northmax,northing)
	  end if
	end do
200	if(atyp.eq.'l')then
	  single_num=eastmin
	  delta=spacing(single_num)
	  eastmin=eastmin-delta
	  single_num=eastmax
	  delta=spacing(single_num)
	  eastmax=eastmax+delta
	  single_num=northmin
	  delta=spacing(single_num)
	  northmin=northmin-delta
	  single_num=northmax
	  delta=spacing(single_num)
	  northmax=northmax+delta
	else
	  eastmin=eastmin-clength
	  eastmax=eastmax+clength
	  northmin=northmin-clength
	  northmax=northmax+clength
	end if
	rewind(unit=inunit,err=9150)

! -- The DXF output file is next written.

 	call write_dxf_header(outunit,eastmin,eastmax,northmin,northmax)
	if(atyp.eq.'l')then
220	  oldeast=-1.0d300
	  oldnorth=-1.0d300
	  icoord=0
	  iend=0
	  polyline: do
	    read(inunit,*,err=300,end=280) easting,northing
	    icoord=icoord+1
	    if(icoord.eq.1)call write_dxf_polyhead(outunit)
	    if((abs(easting-oldeast).gt.clength).or. &
	       (abs(northing-oldnorth).gt.clength)) then
	       call write_dxf_vertex(outunit,easting,northing)
	       oldeast=easting
	       oldnorth=northing
	    end if
	    lasteast=easting
	    lastnorth=northing
	  end do polyline
280	  iend=1
300	  if(icoord.ne.0)then
	    if((oldeast.ne.lasteast).and.(oldnorth.ne.lastnorth))&
	      call write_dxf_vertex(outunit,lasteast,lastnorth)
	    call write_dxf_polyfin(outunit)
	  end if
	  if(iend.eq.1) go to 500
	  go to 220
	else
	  points: do
	    read(inunit,*,end=500) easting,northing
	    call write_dxf_polyhead(outunit)
	    call write_dxf_vertex(outunit,easting-clength,northing)
	    call write_dxf_vertex(outunit,easting+clength,northing)
	    call write_dxf_polyfin(outunit)
	    call write_dxf_polyhead(outunit)
	    call write_dxf_vertex(outunit,easting,northing-clength)
	    call write_dxf_vertex(outunit,easting,northing+clength)
	    call write_dxf_polyfin(outunit)
	  end do points
	end if

500	call write_dxf_finish(outunit)
	write(amessage,510) trim(outfile)
510	format('  - file ',a,' written ok.')
	go to 9900

9000	call num2char(iline,anum1)
	write(amessage,9010) trim(anum1),trim(qfile)
9010	format(' Error reading line ',a,' of QDIGIT file ',a)
	call write_message(leadspace='yes')
	go to 9900
9100	call num2char(iline,anum1)
	write(amessage,9110) trim(anum1),trim(qfile)
9110	format('unexpected word or illegal data item encountered at line ',&
        a,' of QDIGIT file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 9900
9150	write(amessage,9160) trim(qfile)
9160	format('cannot rewind QDIGIT file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 9900

9900    call close_files

end program qdig2dxf
 