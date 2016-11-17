!     Last change:  J     2 May 2002    8:49 pm
program arrdet

! -- Program ARRDET reads a MODFLOW/MT3D unformatted file. It lists the details
!    of arrays found in that file.

	use defn
	use inter

	implicit none

        integer             :: ifail,modunit,outunit,iarray,ierr
        integer             :: ncol,nrow,kstp,kper,ntrans,mcol,mrow,ilay, &
                               irow,icol
        real                :: totim,pertim
        character (len=1)   :: af
        character (len=10)  :: aarray
        character (len=16)  :: text
        character (len=200) :: outfile,modfile
	type (modelgrid)    :: gridspec

        real, allocatable   :: rarray(:,:)


	write(amessage,5)
5       format(' Program ARRDET lists the contents of an unformatted ',  &
        'MODFLOW or MT3D heads/drawdowns/concentrations file.')
	call write_message(leadspace='yes',endspace='yes')

	include 'unformat.inc'

	call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

	ncol=gridspec%ncol
	nrow=gridspec%nrow

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

400	write(6,*)
410	call open_output_file(ifail, &
	' Enter name for output file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  go to 170
	end if

        if(af.eq.'f')then
          write(outunit,420)
420       format(' Time_step',t12,'Stress_period',t30,'Period_time',  &
          t45,'Total_time',t63,'Layer',t75,'Text')
        else
          write(outunit,430)
430       format(' Transport_step',t21,'Time_step',t32,'Stress_period',  &
          t50,'Total_time',t68,'Layer',t78,'Text')
        end if
	iarray=0
	do
	  if(af.eq.'f')then
	    read(modunit,err=9100,end=1000) kstp,kper,pertim,totim, &
	    text,mcol,mrow,ilay
	  else
	    read(modunit,err=9100,end=1000) ntrans,kstp,kper,&
	    totim,text,mcol,mrow,ilay
	  end if
	  iarray=iarray+1
	  if((mcol.le.0).or.(mrow.le.0).or.(ilay.le.0)) go to 9100
	  if((mrow.ne.nrow).or.(mcol.ne.ncol)) then
	    call num2char(iarray,aarray)
	    write(amessage,450) trim(aarray),trim(modfile),trim(gridspec%specfile)
450	    format(' Number of rows and columns read from header to array ',&
	    'number ',a,' from model output file ',a,' does not agree with ',&
	    'grid specifications as read from grid specification file ',a,   &
            '. Alternatively, unformatted file protocol might be a problem - ', &
            'try using the alternative version of this program.')
	    go to 9890
	  end if
          if(iarray.eq.1)then
            allocate(rarray(ncol,nrow),stat=ierr)
            if(ierr.ne.0) then
              write(amessage,470)
470           format(' Cannot allocate sufficient memory to continue ',  &
              'execution.')
              go to 9890
            end if
          end if
	  read(modunit,err=9200,end=9250)   &
          ((rarray(icol,irow),icol=1,ncol),irow=1,nrow)
          text=adjustl(text)
          if(af.eq.'f')then
            write(outunit,520) kstp,kper,pertim,totim,ilay,trim(text)
520         format(t2,i8,t12,i8,t30,1pg14.7,t45,1pg14.7,t60,i8,t75,a)
          else
            write(outunit,530) ntrans,kstp,kper,totim,ilay,trim(text)
530         format(t2,i8,t22,i8,t33,i8,t50,1pg14.7,t66,i8,t78,a)
          end if
        end do

1000    continue
        write(6,*)
        close(unit=modunit)
        call num2char(iarray,aarray)
        write(6,1010) trim(aarray),trim(modfile)
1010    format(' - ',a,' arrays read from file ',a,'.')
        close(unit=outunit)
        write(6,1020) trim(outfile)
1020    format(' - file ',a,' written ok.')

	go to 9900

9100	call num2char(iarray,aarray)
	write(amessage,9110) trim(aarray),trim(modfile)
9110	format(' Error reading header to array number ',a,&
	' from model-generated file ',a)
	go to 9890
9200	call num2char(iarray,aarray)
	write(amessage,9210) trim(aarray),trim(modfile)
9210	format(' Error reading array number ',a,' from model ',&
	'output file ',a)
	go to 9890
9250	call num2char(iarray,aarray)
	write(amessage,9260) trim(aarray),trim(modfile)
9260	format(' Unexpected end of file encountered while reading array ',&
	'number ',a,' from model output file ',a)
	go to 9890

9890    call write_message(leadspace='yes')
9900    call close_files
        deallocate(rarray,stat=ierr)
        write(6,*)

end program arrdet


