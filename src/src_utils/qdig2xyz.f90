!     Last change:  JD   14 Feb 2001    9:04 am
program qdig2xyz

! -- Program QDIG2XYZ produces a "xyz" data file from a QDIGIT-generated
!    contour file.

	use defn
	use inter

	implicit none

	integer                  :: ifail,ierr,iend,outunit,inunit,iline,nbb,ifail1
	double precision	 :: clength,easting,northing,oldeast,&
				    oldnorth,lasteast,lastnorth,val
	character (len=80)  	 :: qfile,aprompt,outfile,bfile
	character (len=20)	 :: anum1


	write(amessage,5)
5       format(' Program QDIG2XYZ produces a "xyz" data file from a ', &
	'QDIGIT-generated contour file.')
	call write_message(leadspace='yes')

9	write(6,*)
10	imessage=0
20	write(6,30,advance='no')
30	format(' Enter name of QDIGIT-generated contour file: ')
	read(5,'(a)') qfile
	if(qfile.eq.' ') go to 20
	qfile=adjustl(qfile)
	if(index(eschar,qfile(1:2)).ne.0) go to 9900
        bfile=qfile
        nbb=len_trim(bfile)
        call getfile(ifail1,bfile,qfile,1,nbb)
	inunit=nextunit()
	open(unit=inunit,file=qfile,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.gt.5) go to 9900
	  write(amessage,50) trim(qfile)
50	  format(' Cannot open file ',a,'  - try again.')
	  call write_message(increment=1)
	  go to 20
	end if

	write(6,*)
113	write(6,115,advance='no')
115	format(' Enter trimming distance: ')
	read(5,'(a)') anum1
	if(anum1.eq.' ') go to 113
	anum1=adjustl(anum1)
	if(index(eschar,anum1(1:2)).ne.0) then
	  write(6,*)
	  close(unit=inunit)
	  go to 10
	end if
	call char2num(ifail,anum1,clength)
	if(ifail.ne.0) go to 113
	if(pos_test(clength, 'trimming distance').ne.0) go to 113
	clength=clength/1.4

120	write(6,*)
	aprompt=' Enter name for output "xyz" file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 113
	end if

	oldeast=-1.0d300
	oldnorth=-1.0d300
	iline=0
	iend=1
	do
	  iline=iline+1
	  if(iline.eq.2000) then
	    write(6,160)
160	    format(/,' Working.....',/)
	  end if
	  read(inunit,'(a)',err=9000,end=500) cline
	  call linesplit(ifail,3)
	  if(ifail.lt.0) cycle
	  if(ifail.gt.0)then
	    if(cline(left_word(1):left_word(1)+6).ne.'ENDLINE') go to 9100
	    if((oldeast.ne.lasteast).and.(oldnorth.ne.lastnorth).and. &
	    (iend.ne.1))write(outunit,170,err=9600) lasteast,lastnorth,val
170	    format(1x,f15.3,2x,f15.3,2x,1pg13.6)
	    iend=1
	    cycle
	  end if
	  easting=char2double(ifail,1)
	  if(ifail.ne.0) go to 9100
	  northing=char2double(ifail,2)
	  if(ifail.ne.0) go to 9100
	  val=char2double(ifail,3)
	  if(ifail.ne.0) go to 9100
	  if(iend.eq.1) then
	    write(outunit,170,err=9600) easting,northing,val
	    oldeast=easting
	    oldnorth=northing
	    iend=0
	  else if((abs(easting-oldeast).gt.clength).or. &
	          (abs(northing-oldnorth).gt.clength)) then
	    write(outunit,170,err=9600) easting,northing,val
	    oldeast=easting
	    oldnorth=northing
	  endif
	  lasteast=easting
	  lastnorth=northing
	end do

500	if(iend.ne.1) then
	  if((oldeast.ne.lasteast).and.(oldnorth.ne.lastnorth)) &
	  write(outunit,170,err=9600) lasteast,lastnorth,val
	  write(amessage,520) trim(qfile)
520	  format(' Warning: file ',a,' does not finish with "ENDLINE" ', &
          'statement.')
	end if
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
9110	format('unexpected word,illegal data or insufficient items ',&
       'encountered at line ',a,' of QDIGIT file ',a)
	call write_message(error='yes',leadspace='yes')
	go to 9900
9600	write(amessage,9610) trim(outfile)
9610	format(' Cannot write to file ',a,': file inaccessible or disk full.')
	call write_message(leadspace='yes')
	go to 9900

9900    call close_files

end program qdig2xyz
 