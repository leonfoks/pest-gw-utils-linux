!     Last change:  J    14 Jun 2002   11:06 pm

program arrayobs

! -- Program ARRAYOBS adds observations to a PEST control file based on the values
!    contained in a real array table file.

        use defn
        use inter

        implicit none

        integer                :: ifail,iline,npar,nobs,npargp,nprior,nobsgp,ierr, &
                                  ncol,nrow,icol,irow,iobsgp,itemp, &
                                  pestunit,pestoutunit,obstabunit,modtabunit,iobs, &
                                  newregflag,mode,newobs,newgroupflag,i,iheader, &
                                  idate,iflag,insunit,ntpfle,ninsfle
        integer, allocatable   :: icovfile(:),obsintarray(:,:),modintarray(:,:)
        real                   :: phimlim,fracphim,weight
        real, allocatable      :: obsarray(:,:)
        character (len=2)      :: aprefix,areg
        character (len=5)      :: acol,arow
        character (len=12)     :: aobs
        character (len=15)     :: aline,pestmode,atemp,anewgroup
        character (len=120)    :: aprompt
        character (len=120)    :: pestfile,pestoutfile,modtabfile,obstabfile,insfile, &
                                  afile,bfile
        character (len=400)    :: dline
        character (len=12), allocatable :: aobsgp(:)
	type (modelgrid)       :: gridspec



        write(amessage,5)
5       format(' Program ARRAYOBS adds observations to a PEST ',&
        'control file based on the values contained in a real array table file.')
        call write_message(leadspace='yes',endspace='yes')

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
6         format(' Cannot read array header specification from settings file ', &
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
        allocate(obsarray(ncol,nrow),obsintarray(ncol,nrow),   &
        modintarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9150
        obsintarray=0        ! an array
        modintarray=0        ! an array

! -- The existing PEST control file is partially read to obtain some information.
!    We need to know the number of observations, number of observation groups,
!    names of observation groups and whether there is a regularisation section in
!    that file.

25      aprompt=' Enter name of PEST control file: '
        call open_input_file(ifail,aprompt,pestfile,pestunit)
        if(ifail.ne.0) go to 9900

        if(escset.eq.1) then
          escset=0
          write(6,*)
          deallocate(obsarray,obsintarray,modintarray,stat=ierr)
          if(ierr.ne.0) go to 9200
          call free_grid_mem(gridspec)
          go to 10
        end if

        iline=1
        read(pestunit,700,err=9000,end=9050) cline
        cline=adjustl(cline)
        call casetrans(cline,'lo')
        if(cline(1:4).ne.'pcf ')then
          write(amessage,30) trim(pestfile)
30        format(' File ',a,' is not a PEST control file - try again.')
          call write_message(leadspace='yes',endspace='yes')
          close(unit=pestunit)
          go to 25
        end if
        iline=iline+1
        read(pestunit,700,err=9000,end=9050) cline
        iline=iline+1
        read(pestunit,700,err=9000,end=9050) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        pestmode=cline(left_word(2):right_word(2))
        call casetrans(pestmode,'lo')
        iline=iline+1
        read(pestunit,'(a)',err=9000,end=9070) cline
        call linesplit(ifail,5)
        if(ifail.ne.0) go to 9100
        npar=char2int(ifail,1)
        if(ifail.ne.0) go to 9000
        nobs=char2int(ifail,2)
        if(ifail.ne.0) go to 9000
        npargp=char2int(ifail,3)
        if(ifail.ne.0) go to 9000
        nprior=char2int(ifail,4)
        if(ifail.ne.0) go to 9000
        nobsgp=char2int(ifail,5)
        if(ifail.ne.0) go to 9000
        allocate(aobsgp(nobsgp),icovfile(nobsgp),stat=ierr)
        if(ierr.ne.0) go to 9150
        iline=iline+1
        read(pestunit,'(a)',err=9000,end=9070) cline
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        ntpfle=char2int(ifail,1)
        if(ifail.ne.0) go to 9000
        ninsfle=char2int(ifail,2)
        if(ifail.ne.0) go to 9000
        do
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9070) cline
          call casetrans(cline,'lo')
          if(index(cline,'* observation gr').ne.0) go to 70
        end do
70      do iobsgp=1,nobsgp
          iline=iline+1
          read(pestunit,'(a)',err=9000,end=9050) cline
          cline=adjustl(cline)
          if(index(cline,'* observation da').ne.0)then
            write(amessage,74) trim(pestfile)
74          format(' Too few observation group names cited in ', &
            '"* observation groups" section of file ',a,'.')
            go to 9890
          end if
          call linesplit(ifail,1)
          if(ifail.ne.0) go to 9100
          atemp=cline(left_word(1):right_word(1))
          if(len_trim(atemp).gt.12)then
            call num2char(iline,aline)
            write(amessage,72) trim(aline),trim(pestfile)
72          format(' Observation group name greater than 12 characters at ', &
            'line ',a,' of file ',a,'.')
            go to 9890
          end if
          aobsgp(iobsgp)=atemp
          call linesplit(ifail,2)
          if(ifail.eq.0)then
            icovfile(iobsgp)=1
          else
            icovfile(iobsgp)=0
          end if
        end do
        iline=iline+1
        read(pestunit,'(a)',err=9000,end=9050) cline
        call casetrans(cline,'lo')
        cline=adjustl(cline)
        if(index(cline,'* observation data').eq.0)then
          call num2char(iline,aline)
          write(amessage,76) trim(aline),trim(pestfile)
76        format(' "* observation data" section header expected at line ',a, &
          ' of file ',a,'.')
          go to 9890
        end if
        rewind(unit=pestunit,iostat=ierr)
        if(ierr.ne.0) then
          write(amessage,605) trim(pestfile)
605       format(' Error rewinding file ',a,'.')
          go to 9890
        end if

! -- More information is now sought from the user.

        write(6,*)
100     continue
        aprompt=' Enter name of observation real array table file: '
        call open_input_file(ifail,aprompt,obstabfile,obstabunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          close(unit=pestunit,err=9500)
          deallocate(aobsgp,icovfile,stat=ierr)
          if(ierr.ne.0) go to 9200
          write(6,*)
          go to 25
        end if
150     continue
        aprompt=' Enter name of model real array table file: '
        call open_input_file(ifail,aprompt,modtabfile,modtabunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          close(unit=obstabunit,err=9500)
          write(6,*)
          go to 100
        end if

160     write(6,170,advance='no')
170     format(' Enter prefix for new observation names (two characters or less): ')
        read(5,'(a)') aprefix
        if(aprefix.eq.' ') go to 160
        aprefix=adjustl(aprefix)
        call casetrans(aprefix,'lo')
        if(aprefix.eq.'e ')then
          close(unit=modtabunit,err=9500)
          write(6,*)
          go to 150
        end if
180     write(6,190,advance='no')
190     format(' Enter weight to assign to new observations: ')
        itemp=key_read(weight)
        if(itemp.lt.0) go to 180
        if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 160
        end if
        if(itemp.gt.0)then
          write(6,195)
195       format(' Data input error  - try again.')
          go to 180
        end if
        if(weight.lt.0.0)then
          write(6,196)
196       format(' The weight cannot be negative - try again.')
          go to 180
        end if
200     write(6,210,advance='no')
210     format(' Enter group name for new observations: ')
        read(5,'(a)') anewgroup
        if(anewgroup.eq.' ') go to 200
        anewgroup=adjustl(anewgroup)
        call casetrans(anewgroup,'lo')
        if(anewgroup(1:2).eq.'e ')then
          write(6,*)
          go to 180
        end if
        if(len_trim(anewgroup).gt.12)then
          write(6,215)
215       format(' Group name must be 12 characters of less in length - try again.')
          go to 200
        end if
        if(index(trim(anewgroup),' ').ne.0)then
          write(6,216)
216       format(' Group name must not contain space character - try again.')
          go to 200
        end if

        write(6,*)
209     aprompt=' Enter name for new PEST control file: '
        call open_output_file(ifail,aprompt,pestoutfile,pestoutunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 200
        end if
        areg='n '
        newregflag=0
        if((pestmode(1:5).ne.'regul').and.(anewgroup.eq.'regul'))then
229       write(6,230,advance='no')
230       format(' Use regularisation mode for new PEST control file?  [y/n]: ')
          read(5,'(a)') areg
          if(areg.eq.' ') go to 229
          areg=adjustl(areg)
          call casetrans(areg,'lo')
          if(areg.eq.'e ')then
            close(unit=pestoutunit,err=9500)
            write(6,*)
            go to 209
          end if
          if((areg.ne.'y ').and.(areg.ne.'n ')) go to 229
          if(areg.eq.'y ')then
            newregflag=1
240         write(6,241,advance='no')
241         format(' Enter value for PHIMLIM: ')
            itemp=key_read(phimlim)
            if(itemp.lt.0) go to 240
            if(itemp.gt.0)then
              write(6,195)
              go to 240
            end if
            if(escset.ne.0) then
              escset=0
              write(6,*)
              go to 229
            end if
            if(phimlim.le.0.0)then
              write(6,250)
250           format(' PHIMLIM must be greater than zero - try again.')
              go to 240
            end if
260         write(6,261,advance='no')
261         format(' Enter value for FRACPHIM: ')
               itemp=key_read(fracphim)
            if(itemp.lt.0) go to 260
            if(itemp.gt.0)then
              write(6,195)
              go to 260
            end if
            if(escset.ne.0) then
              escset=0
              write(6,*)
              go to 240
            end if
            if(fracphim.lt.0.0)then
              write(6,270)
270           format(' FRACPHIM must not be negative - try again.')
              go to 260
            else if(fracphim.gt.0.35)then
              write(6,271)
271           format(' FRACPHIM must be less than 0.35 - try again.')
              go to 260
            end if
          end if
        end if

285     aprompt=' Enter name for instruction file: '
        call open_output_file(ifail,aprompt,insfile,insunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          close(unit=pestoutunit,err=9500)
          write(6,*)
          go to 209
        end if

! -- The "observation" real array table file is read. The values in this file will
!    be added to the "* observation data" section of the new PEST control file.

        mode=1
        call read_real_tab_file_int(ifail,mode,obstabunit,obstabfile,obsarray, &
        obsintarray,gridspec,insunit,insfile,aprefix)
        if(ifail.ne.0) go to 9900
        close(unit=obstabunit)
        newobs=count(obsintarray.ne.0)
        if(newobs.eq.0)then
          write(amessage,275) trim(obstabfile)
275       format(' No elements are cited in real array table file ',a)
          go to 9890
        end if

! -- The "model" real array table file is read. It is simply checked that this is
!    consistent with the observation real array table file.

        mode=2
        call read_real_tab_file_int(ifail,mode,modtabunit,modtabfile,obsarray, &
        modintarray,gridspec,insunit,insfile,aprefix)
        if(ifail.ne.0) go to 9900
        rewind(unit=modtabunit,iostat=ierr)
        if(ierr.ne.0) then
          write(amessage,605) trim(modtabfile)
          go to 9890
        end if

        do irow=1,nrow
          do icol=1,ncol
            if(modintarray(icol,irow).ne.obsintarray(icol,irow))then
              write(amessage,300)
300           format(' The model-generated real array table file is not ', &
              'consistent with the observation real array table file. There is ', &
              'at least one cell featured in one which is not featured in the other.')
              go to 9890
            end if
          end do
        end do

! -- We next establish whether the desired observation group name already exists in the
!    original PEST control file.

        newgroupflag=1
        do iobsgp=1,nobsgp
          if(anewgroup.eq.aobsgp(iobsgp)) then
            if(icovfile(iobsgp).eq.1)then
              write(amessage,305) trim(anewgroup), trim(pestfile)
305           format(' You have said that new observations are to be assigned ', &
              'to the observation group "',a,'". A covariance matrix file ', &
              'has been assigned to this group in file ',a,'; this is not ', &
              'permitted.')
              go to 9890
            else
              newgroupflag=0
              go to 350
            end if
          end if
        end do
350     continue

! -- The instruction file is now written.

        mode=3
        call read_real_tab_file_int(ifail,mode,modtabunit,modtabfile,obsarray, &
        modintarray,gridspec,insunit,insfile,aprefix)
        if(ifail.ne.0) go to 9900
        close(unit=modtabunit)
        close(unit=insunit)
        write(6,880) trim(insfile)

! -- The new PEST control file is now written.

        do i=1,2
          read(pestunit,700,err=9000,end=9050) cline
700       format(a)
          write(pestoutunit,701,err=9300) trim(cline)
701       format(a)
        end do
        read(pestunit,700,err=9000,end=9050) cline
        if(newregflag.eq.0)then
          write(pestoutunit,701,err=9300) trim(cline)
        else
          call linesplit(ifail,1)
          write(pestoutunit,710,err=9300) cline(left_word(1):right_word(1))
710       format(a,'   regularisation')
        end if
        read(pestunit,700,err=9000,end=9050) cline
        write(pestoutunit,715) npar,nobs+newobs,npargp,nprior, &
        nobsgp+newgroupflag
715     format(5i7)
        read(pestunit,700,err=9000,end=9050) cline
        call linesplit(ifail,2)
        write(pestoutunit,716) ntpfle,ninsfle+1,trim(cline(right_word(2)+1:))
716     format(i5,i5,2x,a)
        do
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
          call casetrans(cline,'lo')
          if(index(cline,'* observation gr').ne.0) exit
        end do
        do i=1,nobsgp
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
        end do
        if(newgroupflag.eq.1)then
          write(pestoutunit,701) trim(anewgroup)
        end if
        read(pestunit,700,err=9000,end=9050) cline
        write(pestoutunit,701,err=9300) trim(cline)
        do iobs=1,nobs
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
        end do
        do irow=1,nrow
          do icol=1,ncol
            if(obsintarray(icol,irow).ne.0)then
              call num2char(irow,arow)
              call num2char(icol,acol)
              aobs=trim(aprefix)//'_'//trim(arow)//'_'//trim(acol)
              write(pestoutunit,717) trim(aobs),obsarray(icol,irow),weight, &
              trim(anewgroup)
717           format(a,t14,1pg14.7,2x,1pg14.7,2x,a)
            end if
          end do
        end do

! -- The new observations have now been added to the PEST control file.
!    The name of the extra instruction file and model output file must be added.

        do
          read(pestunit,700,err=9000,end=9600) cline
          write(pestoutunit,701,err=9300) trim(cline)
          call casetrans(cline,'lo')
          if(index(cline,'* model input/out').ne.0) exit
        end do
        do i=1,ntpfle
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
        end do
        do i=1,ninsfle
          read(pestunit,700,err=9000,end=9050) cline
          write(pestoutunit,701,err=9300) trim(cline)
        end do
        call addquote(insfile,afile)
        call addquote(modtabfile,bfile)
        write(pestoutunit,718,err=9300) trim(afile),trim(bfile)
718     format(a,2x,a)

!    The remainder of the PEST control file is now written.

        iflag=0
        if(newregflag.eq.0)then
          do
            read(pestunit,700,err=9000,end=1000) cline
            write(pestoutunit,701,err=9300) trim(cline)
          end do
        else
          do
            read(pestunit,700,err=9000,end=800) cline
            if(cline(1:1).eq.'*')then
              dline=cline
              call casetrans(dline,'lo')
              if(index(dline,'predict').ne.0) go to 800
              if(index(dline,'regul').ne.0)then
                iflag=1
                go to 800
              end if
            end if
            if(cline.ne.' ')write(pestoutunit,701,err=9300) trim(cline)
          end do
800       continue
          write(pestoutunit,810,err=9300)
810       format('* regularisation')
          write(pestoutunit,820) phimlim,phimlim*1.05,fracphim
820       format(3(1x,1pg14.7))
          write(pestoutunit,830)
830       format(' 1.0  1.0e-8  1.0e8')
          write(pestoutunit,840)
840       format(' 1.3  1.0e-2')
        end if
1000    close(unit=pestunit)
        close(unit=pestoutunit)
        write(6,880) trim(pestoutfile)
880     format(' - file ',a,' written ok.')
        if(iflag.ne.0)then
          write(amessage,850) trim(pestfile), trim(pestoutfile)
850       format(' Warning: regularisation section in PEST control file ',a, &
          ' has been overwritten in new PEST control file ',a)
          call write_message
        end if
        go to 9992


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(pestfile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(pestfile)
9060    format(' Unexpected end encountered to file ',a,'.')
        go to 9890
9070    continue
        write(amessage,9080) trim(pestfile)
9080    format(' Unexpected end to PEST control file ',a,' - looking for ', &
        '"* observation groups" section.')
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(pestfile)
9110    format(' Insufficient entries on line ',a,' of file ',a,'.')
        go to 9890
9150    write(amessage,9160)
9160    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Memory management error - cannot de-allocate memory.')
        go to 9890
9300    continue
        write(amessage,9310) trim(pestoutfile)
9310    format(' Error writing to new PEST control file ',a,'.')
        go to 9890
9500    write(amessage,9510)
9510    format(' Cannot close file - execution cannot continue.')
        go to 9890
9600    continue
        write(amessage,9610) trim(pestfile)
9610    format(' Unexpected end to PEST control file ',a,' - looking for ', &
        '"* model input/output" section.')
        go to 9890

9890    call write_message(leadspace='yes')
9900    continue
        write(6,9991)
9991    format(/,' Execution terminated.',/)
9992    continue
        call close_files
        deallocate(obsarray,obsintarray,modintarray,aobsgp,stat=ierr)


end program arrayobs

