program getmularr1

! -- Program GETMULARR1 extracts multiple arrays pertaining to a single time from a MODFLOW or MT3D unformatted
!    output file. It re-writes these as a another pseudo MODFLOW or MT3D unformatted output file.

	use defn
	use inter

	implicit none

        integer             :: ifail,ierr,iline
        integer             :: ilay,iicount1,iicount2,i
        integer             :: ncol,nrow,mrow,mcol,irow,icol
        integer             :: kstp,kper,ntrans
        integer             :: modunit,outunit
        real                :: totimfind,totim,pertim
        character (len=1)   :: af
        character (len=10)  :: anum
        character (len=16)  :: text
        character (len=100) :: aprompt
        character (len=200) :: modfile,outfile,afile
        real, allocatable   :: rarray(:,:)
	type (modelgrid)    :: gridspec


	write(amessage,5)
5       format(' Program GETMULARR1 extracts arrays corresponding to a single ',  &
        'simulation time from a MODFLOW or MT3DMS unformatted output file and writes ',  &
        'them to a single pseudo MODFLOW or MT3DMS unformatted output file.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'

        iicount1=0
        iicount2=0

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

        write(6,*)
190     write(6,195,advance='no')
195     format(' Enter simulation time for which to extract arrays: ')
	i=key_read(totimfind)
	if(escset.ne.0)then
	  escset=0
	  write(6,*)
	  go to 170
	else if(i.eq.-1) then
	  go to 190
	else if(i.ne.0) then
	  write(6,196)
196	  format(' Illegal input  - try again.')
	  go to 190
	end if
	if(totimfind.lt.0.0) then
	  write(amessage,200)
200	  format(' Simulation time must be positive  - try again.')
	  call write_message
	  go to 190
	end if

! -- The name of the unformatted output file is acquired.

       write(6,*)
310    aprompt = ' Enter name for unformatted output file: '
       call open_output_file(ifail,aprompt,outfile,outunit,file_format='unformatted')
       if(ifail.ne.0) go to 9900
       if(escset.eq.1)then
         escset=0
         write(6,*)
         go to 190
       end if

! -- Memory is allocated for one real array.

        allocate (rarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,210)
210       format(' Cannot allocate sufficient memory to continue execution.')
          go to 9890
        end if

! -- Now the array is located in the MODFLOW/MT3D output file.

        do
	  if(af.eq.'f')then
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=1000) kstp,kper,pertim,totim, &
	    text,mcol,mrow,ilay
	  else
	    mrow=0
	    mcol=0
	    read(modunit,err=9100,end=1000) ntrans,kstp,kper,&
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
          iicount1=iicount1+1
          if(equals(totim,totimfind))then
            iicount2=iicount2+1
            if(af.eq.'f')then
              write(outunit) kstp,kper,pertim,totim,text,mcol,mrow,ilay
            else
              write(outunit) ntrans,kstp,kper,totim,text,mcol,mrow,ilay
            end if
            write(outunit) ((rarray(icol,irow),icol=1,ncol),irow=1,nrow)
          end if
        end do

1000    continue
        close(unit=modunit)
        close(unit=outunit)

        write(6,*)
        call num2char(iicount1,anum)
        call addquote(modfile,afile)
        write(6,1010) trim(anum),trim(afile)
1010    format(' - ',a,' arrays read from file ',a,'.')
        call num2char(iicount2,anum)
        call addquote(outfile,afile)
        write(6,1020) trim(anum),trim(afile)
1020    format(' - ',a,' arrays written to file ',a,'.')


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


9890	call write_message(leadspace='yes')
9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(rarray,stat=ierr)
	write(6,*)

end program getmularr1
