program getmularr

! -- Program GETMULARR extracts multiple arrays from a MODFLOW or MT3D unformatted
!    output file.

	use defn
	use inter

	implicit none

        integer             :: ifail,idate,iheader,ierr,iline
        integer             :: ilay,ilayfind,oldilayfind
        integer             :: ncol,nrow,mrow,mcol,irow,icol
        integer             :: ibeg,iend
        integer             :: kstp,kper,ntrans
        integer             :: modunit,extunit,outunit
        real                :: oldtotimfind,totimfind,totim,pertim
        character (len=1)   :: af
        character (len=10)  :: aline
        character (len=16)  :: text
        character (len=200) :: modfile,extfile,outfile
        real, allocatable   :: rarray(:,:)
	type (modelgrid)    :: gridspec


	write(amessage,5)
5       format(' Program GETMULARR extracts multiple arrays from a MODFLOW or ', &
        'MT3DMS unformatted output file and writes them to separate files.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'
	
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
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
	  call write_message
	  go to 9900
	end if

	call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow

! -- The unformatted MODFLOW/MT3D output file is opened.

150	write(6,*)
	call open_input_file(ifail, &
	' Enter name of unformatted model-generated file: ',modfile,modunit, &
	file_format='unformatted')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          go to 10
        end if

170	write(6,180,advance='no')
180	format(' Is this a MODFLOW or MT3D file?  [f/t]: ')
	read(5,'(a)') af
	if(af.eq.' ') go to 170
	call casetrans(af,'lo')
	if(index(eschar,af).ne.0) then
	  close(unit=modunit)
	  go to 150
	end if
	if((af.ne.'f').and.(af.ne.'t')) go to 170

! -- The array extraction file is now opened.

        write(6,*)
200     continue
        call open_input_file(ifail, &
        ' Enter name of array extraction file: ',extfile,extunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 170
        end if

! -- Memory is allocated for one real array.

        allocate (rarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,210)
210       format(' Cannot allocate sufficient memory to continue execution.')
          go to 9890
        end if

! -- The array extraction file is read line by line and the corresponding array
!    extracted from the MODFLOW output file.

        oldtotimfind=-1.0e35
        oldilayfind=-99999999
        iline=0
300     continue

! -- First we identify a time and layer number to read an array for.

        iline=iline+1
        call num2char(iline,aline)
        read(extunit,'(a)',end=1000) cline
        if(cline.eq.' ') go to 300
        cline=adjustl(cline)
        if(cline(1:1).eq.'#') go to 300
        call linesplit(ifail,3)
        if(ifail.ne.0) then
          write(amessage,310) trim(aline),trim(extfile)
310       format(' Three entries are expected on line ',a,' of file ',a,'.')
          go to 9890
        end if
        totimfind=char2real(ifail,1)
        if(ifail.ne.0)then
          write(amessage,315) trim(aline),trim(extfile)
315       format(' Cannot read simulation time at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(totimfind.lt.0.0)then
          write(amessage,317) trim(aline),trim(extfile)
317       format(' Negative simulation time supplied at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(totimfind.lt.oldtotimfind)then
          write(amessage,320) trim(aline),trim(extfile)
320       format(' Simulation time supplied at ',   &
          'line ',a,' of file ',a,' is less than previous simulation time.')
          go to 9890
        end if
        ilayfind=char2int(ifail,2)
        if(ifail.ne.0)then
          write(amessage,322) trim(aline),trim(extfile)
322       format(' Cannot read layer number at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(ilayfind.lt.0)then
          write(amessage,323) trim(aline),trim(extfile)
323       format(' Negative layer number supplied at line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(equals(oldtotimfind,totimfind))then
          if(ilayfind.le.oldilayfind)then
            write(amessage,330) trim(aline),trim(extfile)
330         format(' Model layer numbers are not supplied in increasing order at ', &
            'line ',a,' of file ',a,'.')
            go to 9890
          end if
        end if

! -- The name of the output file to which to write the array is also identified.

        ibeg=left_word(3)
        iend=len_trim(cline)
        call getfile(ifail,cline,outfile,ibeg,iend)
        if(ifail.ne.0)then
          write(amessage,340) trim(aline),trim(extfile)
340       format(' Cannot read filename from line ',a,' of file ',a,'.')
          go to 9890
        end if
        outunit=nextunit()
        open(unit=outunit,file=outfile,action='write',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,350) trim(outfile),trim(aline),trim(extfile)
350       format(' Cannot write to file ',a,' cited on line ',a,' of file ',a,'.')
          go to 9890
        end if

        oldtotimfind=totimfind
        oldilayfind=ilayfind

! -- Now the array is located in the MODFLOW/MT3D output file.

360     continue
        do
	  if(af.eq.'f')then
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=9300) kstp,kper,pertim,totim, &
	    text,mcol,mrow,ilay
	  else
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=9300) ntrans,kstp,kper,&
	    totim,text,mcol,mrow,ilay
	  end if
	  if((mcol.le.0).or.(mrow.le.0).or.(ilay.le.0)) go to 9100
	  if((mrow.ne.nrow).or.(mcol.ne.ncol)) then
	    write(amessage,450) trim(modfile),trim(gridspec%specfile)
450	    format(' Number of rows and columns read from header to array ',&
	    'from model output file ',a,' does not agree with ',&
	    'grid specifications as read from grid specification file ',a,   &
            '. Alternatively, unformatted file protocol might be a problem - ', &
            'try using the alternative version of this program.')
	    go to 9890
	  end if
          read(modunit,err=9200,end=9250) ((rarray(icol,irow),icol=1,ncol),irow=1,nrow)
          if(equals(totim,totimfind))then
            if(equals(ilay,ilayfind))then
              if(headerspec.eq.'yes')then
                write(outunit,460) ncol,nrow
460             format(2i10)
              end if
              do irow=1,nrow
                write(outunit,470) (rarray(icol,irow),icol=1,ncol)
470             format(8(1x,1pg14.7))
              end do
              close(unit=outunit)
              write(6,471) trim(outfile)
471           format(' - file ',a,' written ok.')
              go to 300
            end if
          end if
          if(totim.gt.totimfind) go to 9300
        end do

1000    continue
        go to 9900

9100	continue
	write(amessage,9110) trim(modfile)
9110	format(' Error reading header to array from model-generated file ',a)
	go to 9890
9200	continue
	write(amessage,9210) trim(modfile)
9210	format(' Error reading array from model output file ',a)
	go to 9890
9250	continue
	write(amessage,9260) trim(modfile)
9260	format(' Unexpected end of file encountered while reading array ',&
	' from model output file ',a)
	go to 9890
9300    write(amessage,9310) trim(aline),trim(extfile)
9310    format(' Cannot find array corresponding to simulation time and layer ',  &
        'number cited at line ',a,' of file ',a,'.')
        go to 9890


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(rarray,stat=ierr)
	write(6,*)

end program getmularr
