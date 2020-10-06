program zonmdef

! -- Program ZONMDEF creates a large number of zone-based parameters for use with the adjoint state process
!    of MODFLOW 2005.

	use defn
	use inter

	implicit none

        integer, parameter :: MAXR=1000
        integer, parameter :: MINR=-1000
        integer, parameter :: MAXC=1000
        integer, parameter :: MINC=-1000
        integer  :: ifail,idate,iheader,ierr,n,icount,itemp,izone,nzone,iflag,ndim,mpar
        integer  :: ncol,nrow,nlay,icol,irow,ilay,mrow,mcol,irowmax,irowmin,icolmax,icolmin, &
                    startrow,startcol,ir,ic,irhigh,irlow,ichigh,iclow,mlay
        integer  :: disunit,pestunit,outunit

        double precision :: pval,lbound,ubound

        character*1   :: ayn
        character*2   :: atrans
        character*8   :: parbase
        character*10  :: alay,atemp
        character*12  :: apar,apartrans,achange
        character*100 :: aprompt
        character*200 :: disfile,pestfile,intfile,ircfile,realfile,outfile,afile,bfile

        integer              :: rowindex(MINR:MAXR),colindex(MINC:MAXC)
        integer, allocatable :: iarray(:,:,:),numact(:),idum(:,:)

	type (modelgrid) gridspec


        write(amessage,5)
5       format(' Program ZONMDEF creates zone-based parameters for use with the adjoint state process ',  &
        'of MODFLOW 2005.')
        call write_message(leadspace='yes',endspace='yes')

! -- The settings file is read.

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
        if((iheader.ne.0).or.(headerspec.eq.' ')) then
          write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
          go to 9890
        end if

! -- Grid specifications are obtained.

        call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')
        ncol=gridspec%ncol
        nrow=gridspec%nrow

50      write(6,60,advance='no')
60      format(' How many layers in the model? ')
        if(key_read(nlay).ne.0) go to 50
        if(escset.ne.0)then
          escset=0
          call free_grid_mem(gridspec)
          write(6,*)
          go to 10
        end if
        if(nlay.le.0) go to 50

        allocate(iarray(ncol,nrow,nlay),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(numact(nlay),stat=ierr)
        if(ierr.ne.0) go to 9200
        iarray=0                    ! an array
        numact=-9999                ! an array

! -- The parameter basename is acquired.

        write(6,*)
100     write(6,110,advance='no')
110     format(' Enter basename for new parameters: ')
        read(5,*,err=100) parbase
        call casetrans(parbase,'lo')
        if(parbase(1:2).eq.'e ')then
          deallocate(iarray,numact,stat=ierr)
          if(ierr.ne.0)then
            write(amessage,120)
120         format(' Cannot deallocate memory in order to backtrack.')
            go to 9890
          end if
          go to 50
        end if
        n=len_trim(parbase)
        if(n.gt.4)then
          write(6,130)
130       format(/,' *** Parameter basename must be 4 characters or less - try again ***',/)
          go to 100
        end if
140     write(6,150,advance='no')
150     format(' Enter parameter transformation status ([l/n] for "log" or "none"): ')
        read(5,*) atrans
        call casetrans(atrans,'lo')
        if(atrans(1:2).eq.'e ')then
          write(6,*)
          go to 100
        end if
        if((atrans.ne.'l').and.(atrans.ne.'n')) go to 140
160     write(6,170,advance='no')
170     format(' Enter initial value for new parameters: ')
        if(key_read(pval).ne.0) go to 160
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 140
        end if
        if(atrans.eq.'l')then
          if(pval.le.0.0d0)then
            write(6,180)
180         format(/,' *** Must be greater than zero for log-transformed parameters ***',/)
            go to 160
          end if
        end if
190     write(6,200,advance='no')
200     format(' Enter lower bound for new parameters: ')
        if(key_read(lbound).ne.0) go to 190
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 160
        end if
        if(atrans.eq.'l')then
          if(lbound.le.0.0d0)then
            write(6,180)
            go to 190
          end if
        end if
        if(lbound.gt.pval)then
          write(6,201)
201       format(/,' *** Must not be greater than initial value - try again ***',/)
          go to 190
        end if
210     write(6,220,advance='no')
220     format(' Enter upper bound for new parameters: ')
        if(key_read(ubound).ne.0) go to 210
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 190
        end if
        if(ubound.le.lbound)then
          write(6,230)
230       format(/,' *** Must be greater than lower bound - try again ***',/)
          go to 210
        end if
        if(ubound.lt.pval)then
          write(6,231)
231       format(/,' *** Must not be less than initial value - try again ***',/)
          go to 210
        end if

! -- Parameter activity arrays are now read.

        write(6,*)
        ilay=0
235     continue
        ilay=ilay+1                           ! Start of loop
        if(ilay.gt.nlay) go to 300
          call num2char(ilay,alay)
236       continue
          write(6,240,advance='no') trim(alay)
240       format(' Does the parameter exist in layer ',a,'?  [y/n]: ')
          read(5,'(a)') ayn
          if(ayn.eq.' ') go to 236
          call casetrans(ayn,'lo')
          if(ayn.eq.'e')then
            write(6,*)
            if(ilay.eq.1)then
              go to 210
            else
              ilay=ilay-2
              go to 235
            end if
          else if(ayn.eq.'n')then
            numact(ilay)=0
          else
            if(ayn.ne.'y') go to 236
285         continue
            aprompt=' Enter integer parameter activity array file for layer '//trim(alay)//': '
            call read_integer_array(ifail,aprompt,iarray(:,:,ilay),  &
            pm_header=headerspec,rows=nrow,columns=ncol)
            if(ifail.ne.0) go to 9900
            if(escset.ne.0)then
              escset=0
	      write(6,*)
              go to 236
            end if
            icount=0
            do irow=1,nrow
              do icol=1,ncol
                if(iarray(icol,irow,ilay).ne.0) icount=icount+1
              end do
            end do
            if(icount.eq.0)then
              write(6,290)
290           format(/,' *** This array has no non-zero integers - try again ***',/)
              go to 285
            end if
            numact(ilay)=icount
          end if
        go to 235            ! end of loop
300     continue
        icount=0
        do ilay=1,nlay
          icount=icount+numact(ilay)
        end do
        if(icount.eq.0)then
          write(amessage,302)
302       format(' The parameter must be defined in at least one model layer.')
          go to 9890
        end if

! -- Parameter zone dimensions are now acquired.

        write(6,*)
420     write(6,430,advance='no')
430     format(' Enter number of cells in row direction for each parameter zone: ')
        if(key_read(mcol).ne.0) go to 420
        if(escset.ne.0)then
          escset=0
          write(6,*)
          ilay=nlay-1
          go to 235
        end if
        if(mcol.le.0) go to 420

400     write(6,410,advance='no')
410     format(' Enter number of cells in column direction for each parameter zone: ')
        if(key_read(mrow).ne.0) go to 400
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 420
        end if
        if(mrow.le.0) go to 400

! -- The names of output files are now acquired.

        write(6,*)
320     call open_output_file(ifail, &
        ' Enter name for distributed-to-PEST-parameter file: ',disfile,disunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 420
        end if

340     call open_output_file(ifail, &
        ' Enter name for PEST building block file: ',pestfile,pestunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          close(unit=disunit)
          write(6,*)
          go to 320
        end if

        write(6,*)
        write(6,350)
350     format(' - defining zonation for new parameter...')

! -- First we find the entire span of active cells.

        allocate(idum(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200
        idum=0             ! an array
        do irow=1,nrow
          do icol=1,ncol
            do ilay=1,nlay
              if(iarray(icol,irow,ilay).ne.0)then
                idum(icol,irow)=1
                go to 351
              end if
            end do
351         continue
          end do
        end do

        irowmax=0
        irowmin=99999999
        icolmax=0
        icolmin=99999999
        do irow=1,nrow
          do icol=1,ncol
            if(idum(icol,irow).ne.0)then
              itemp=icol
              if(icol.lt.icolmin) icolmin=icol
              go to 356
            end if
          end do
          go to 357
356       continue
          do icol=ncol,1,-1
            if(idum(icol,irow).ne.0)then
              itemp=icol
              if(icol.gt.icolmax) icolmax=icol
              go to 357
            end if
          end do
357       continue
        end do

        do icol=1,ncol
          do irow=1,nrow
            if(idum(icol,irow).ne.0)then
              itemp=irow
              if(itemp.lt.irowmin) irowmin=itemp
              go to 361
            end if
          end do
          go to 362
361       do irow=nrow,1,-1
            if(idum(icol,irow).ne.0)then
              itemp=irow
              if(itemp.gt.irowmax) irowmax=itemp
              go to 362
            end if
          end do
362       continue
        end do

! -- We now build a zonation index array.

        startrow=(irowmin+irowmax)/2
        if(startrow.lt.1) startrow=1
        startcol=(icolmin+icolmax)/2
        if(startcol.lt.1) startcol=1

        irhigh=0
        irlow=0
        ichigh=0
        iclow=0

        ir=-1
        irow=startrow-mrow
        do
          irow=irow+mrow
          if(irow.le.irowmax)then
            ir=ir+1
            if(ir.gt.MAXR)then
              write(amessage,395)
395           format(' Increase MAXR and re-compile program.')
              go to 9890
            end if
            rowindex(ir)=irow
          else
            go to 398
          end if
        end do
398     continue
        irhigh=ir
        if(rowindex(irhigh).lt.irowmax)then
          if(irowmax-rowindex(irhigh)+1.lt.mrow/2)then
            rowindex(irhigh)=irowmax+1
          else
            irhigh=irhigh+1
            if(irhigh.gt.MAXR)then
              write(amessage,395)
              go to 9890
            end if
            rowindex(irhigh)=irowmax+1
          end if
        else
          if(mrow.le.2)then
            irhigh=irhigh+1
            if(irhigh.gt.MAXR)then
              write(amessage,395)
              go to 9890
            end if
            rowindex(irhigh)=irowmax+1
          else
            rowindex(irhigh)=irowmax+1
          end if
        end if
        ir=0
        irow=startrow
        do
          irow=irow-mrow
          if(irow.ge.irowmin)then
            ir=ir-1
            if(ir.lt.MINR)then
              write(amessage,409)
409           format(' Decrease MINR and re-compile program.')
              go to 9890
            end if
            rowindex(ir)=irow
          else
            go to 419
          end if
        end do
419     continue
        irlow=ir
        if(rowindex(irlow).gt.irowmin)then
          if(rowindex(irlow)-irowmin+1.lt.mrow/2)then
            rowindex(irlow)=irowmin
          else
            irlow=irlow-1
            if(irlow.lt.MINR)then
              write(amessage,409)
              go to 9890
            end if
            rowindex(irlow)=irowmin
          end if
        end if

        ic=-1
        icol=startcol-mcol
        do
          icol=icol+mcol
          if(icol.le.icolmax)then
            ic=ic+1
            if(ic.gt.MAXC)then
              write(amessage,431)
431           format(' Increase MAXC and re-compile program.')
              go to 9890
            end if
            colindex(ic)=icol
          else
            go to 450
          end if
        end do
450     continue
        ichigh=ic
        if(colindex(ichigh).lt.icolmax)then
          if(icolmax-colindex(ichigh)+1.lt.mcol/2)then
            colindex(ichigh)=icolmax+1
          else
            ichigh=ichigh+1
            if(ichigh.gt.MAXC)then
              write(amessage,431)
              go to 9890
            end if
            colindex(ichigh)=icolmax+1
          end if
        else
          if(mcol.le.2)then
            ichigh=ichigh+1
            if(ichigh.gt.MAXC)then
              write(amessage,431)
              go to 9890
            end if
            colindex(ichigh)=icolmax+1
          else
            colindex(ichigh)=icolmax+1
          end if
        end if
        ic=0
        icol=startcol
        do
          icol=icol-mcol
          if(icol.ge.icolmin)then
            ic=ic-1
            if(ic.lt.MINC)then
              write(amessage,460)
460           format(' Decrease MINC and re-compile program.')
              go to 9890
            end if
            colindex(ic)=icol
          else
            go to 480
          end if
        end do
480     continue
        iclow=ic
        if(colindex(iclow).gt.icolmin)then
          if(colindex(iclow)-icolmin+1.lt.mcol/2)then
            colindex(iclow)=icolmin
          else
            iclow=iclow-1
            if(iclow.lt.MINC)then
              write(amessage,460)
              go to 9890
            end if
            colindex(iclow)=icolmin
          end if
        end if

        izone=0
        do ir=irlow,irhigh-1
          do ic=iclow,ichigh-1
            iflag=0
            do irow=rowindex(ir),rowindex(ir+1)-1
              do icol=colindex(ic),colindex(ic+1)-1
                if(idum(icol,irow).ne.0) then
                  if(iflag.eq.0)then
                    iflag=1
                    izone=izone+1
                  end if
                  idum(icol,irow)=izone
                end if
              end do
            end do
          end do
        end do
        nzone=izone

! -- Arrays in all layers are now given values.

        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            do irow=1,nrow
              do icol=1,ncol
                if(iarray(icol,irow,ilay).ne.0) iarray(icol,irow,ilay)=idum(icol,irow)
              end do
            end do
          end if
        end do
        write(6,360)
360     format(' - zonation defined ok.')

! -- Integer zonation arrays are now written for each layer.

        write(6,*)
        write(6,515)
515     format(' - writing integer zonation arrays...')
        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            outunit=nextunit()
            call num2char(ilay,alay)
            outfile=trim(parbase)//'_lay'//trim(alay)//'.inf'
            open(unit=outunit,file=outfile)
            if(headerspec.eq.'yes') then
              write(outunit,520) ncol,nrow
520           format(2i10)
            end if
            do irow=1,nrow
              write(outunit,530) (iarray(icol,irow,ilay),icol=1,ncol)
530           format(20i6)
            end do
            close(unit=outunit)
            write(6,540) trim(outfile)
540         format(' - file ',a,' written ok.')
          end if
        end do

! -- Template files for the IRC files are now written.

        write(6,*)
        write(6,550)
550     format(' - writing parameter template files...')
        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            outunit=nextunit()
            call num2char(ilay,alay)
            outfile=trim(parbase)//'_lay'//trim(alay)//'.tpl'
            open(unit=outunit,file=outfile)
            write(outunit,560)
560         format('ptf $')
            do irow=1,nrow
              do icol=1,ncol
                if(iarray(icol,irow,ilay).eq.0) then
                  write(outunit,570)
570               format(1x,'0',t15,'0.0')
                  go to 561
                end if
              end do
            end do
561         continue
            do izone=1,nzone
              do irow=1,nrow
                do icol=1,ncol
                  if(iarray(icol,irow,ilay).eq.izone) go to 572
                end do
              end do
              go to 573
572           continue
              call num2char(izone,atemp)
              apar=trim(parbase)//'_'//trim(atemp)
              write(outunit,580) trim(atemp),trim(apar)
580           format(1x,a,t15,'$',a,t30,'$')
573           continue
            end do
            close(unit=outunit)
            write(6,590) trim(outfile)
590         format(' - file ',a,' written ok.')
          end if
        end do

! -- INT2REAL keyboard input files are now written.

        write(6,*)
        write(6,610)
610     format(' - writing INT2REAL keyboard input files...')
        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            call num2char(ilay,alay)
            outfile='int2real_'//trim(parbase)//'_lay'//trim(alay)//'.in'
            outunit=nextunit()
            open(unit=outunit,file=outfile)
            call addquote(gridspec%specfile,afile)
            write(outunit,'(a)') trim(afile)
            intfile=trim(parbase)//'_lay'//trim(alay)//'.inf'
            write(outunit,'(a)') trim(intfile)
            write(outunit,620)
620         format('c')
            write(outunit,630)
630         format('f')
            ircfile=trim(parbase)//'_lay'//trim(alay)//'.irc'
            write(outunit,640) trim(ircfile)
640         format(a)
            realfile=trim(parbase)//'_lay'//trim(alay)//'.ref'
            write(outunit,640) trim(realfile)
            close(unit=outunit)
            write(6,590) trim(outfile)
          end if
        end do

! -- The distributed-to-PEST-parameter file is now written.

        write(6,*)
        write(6,645)
645     format(' - writing distributed-to-PEST-parameter file...')
        ndim=0
        do ilay=1,nlay
          ndim=ndim+numact(ilay)
        end do
        write(disunit,650) ndim,0
650     format(2i10)
        mlay=0
        do ilay=1,nlay
          if(numact(ilay).ne.0)mlay=mlay+1
        end do
        write(disunit,650) mlay
        do ilay=1,nlay
          if(numact(ilay).ne.0) write(disunit,650) ilay
        end do
        mpar=nzone
        write(disunit,650) mpar
        do izone=1,nzone
          call num2char(izone,atemp)
          apar=trim(parbase)//'_'//trim(atemp)
          write(disunit,'(a)') trim(apar)
        end do
        write(disunit,660) 1
660     format(i10)
        do irow=1,nrow
          do icol=1,ncol
            itemp=idum(icol,irow)
            if(itemp.ne.0)then
              do ilay=1,nlay
                if(iarray(icol,irow,ilay).ne.0)then
                  write(disunit,670) icol,irow,ilay,1,itemp
670               format(5i10,'    1.0')
                end if
              end do
            end if
          end do
        end do
        close(unit=disunit)
        write(6,590) trim(disfile)

! -- The PEST building block file is now written.

        write(6,*)
        write(6,675)
675     format(' - writing PEST building block file...')
        write(pestunit,700)
700     format('* parameter groups')
        write(pestunit,710) trim(parbase)
710     format(a,'  relative 0.01 0.0 switch 2.0 parabolic')
        write(pestunit,720)
720     format('* parameter data')
        if(atrans.eq.'l')then
          apartrans='log'
          achange='factor'
        else
          apartrans='none'
          achange='relative'
        end if
        do izone=1,nzone
          call num2char(izone,atemp)
          apar=trim(parbase)//'_'//trim(atemp)
          write(pestunit,730) trim(apar),trim(apartrans),trim(achange),   &
          pval,lbound,ubound,trim(parbase)
730       format(a,t15,a,t25,a,t35,3(1pg14.7,1x),2x,a,'  1.0  0.0   0')
        end do
        write(pestunit,735)
735     format('* model input/output')
        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            call num2char(ilay,alay)
            outfile=trim(parbase)//'_lay'//trim(alay)//'.tpl'
            call addquote(outfile,afile)
            ircfile=trim(parbase)//'_lay'//trim(alay)//'.irc'
            call addquote(ircfile,bfile)
            write(pestunit,736) trim(afile),trim(bfile)
736         format(a,3x,a)
          end if
        end do
        close(unit=pestunit)
        write(6,590) trim(pestfile)

! -- Some instructions are written to the screen.

        write(6,*)
        write(6,750)
750     format(' The following command(s) must be added to the model batch files:-')
        write(6,*)
        do ilay=1,nlay
          if(numact(ilay).gt.0)then
            call num2char(ilay,alay)
            outfile='int2real_'//trim(parbase)//'_lay'//trim(alay)//'.in'
            write(6,760) trim(outfile)
760         format('     int2real < ',a)
          end if
        end do
        write(amessage,770)
770     format(' Note that the IRC file(s) must be built from the template files before ',  &
        'these commands can be run.')
        call write_message(leadspace='yes')
        write(6,780)
780     format(/,' MODFLOW multiplier files are named as follows:-')
        do ilay=1,nlay
          if(numact(ilay).ne.0)then
            call num2char(ilay,alay)
            realfile=trim(parbase)//'_lay'//trim(alay)//'.ref'
            write(6,790) trim(realfile)
790         format('     ',a)
          end if
        end do

        go to 9900

9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890


9890    call write_message(leadspace='yes')
9900    call close_files

        if(allocated(iarray)) deallocate(iarray)
        if(allocated(numact)) deallocate(numact)
        if(allocated(idum)) deallocate(idum)

end program zonmdef


