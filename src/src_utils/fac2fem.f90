!     Last change:  JD   28 Dec 2002   10:39 pm

program fac2fem

! -- Program FAC2FEM writes mesh data to a MICROFEM input file based on the
!    interpolation factor file produced by PPK2FACF.

	use defn
	use inter

	implicit none

        integer               :: ifail,facunit,ierr,nerr,nbb, &
                                 icellno,na,i,itrans,nnppt,colrep,propunit1, &
                                 propunit2,nummesh,iline,ifem,lt
        integer, allocatable, dimension(:)  :: ipt
        real                                :: rtemp,sum,rtemp1,rlo,rhi,rtemp2,backval
        real, allocatable, dimension(:)     :: wt
        real, allocatable, dimension(:)     :: realarray,minarray,maxarray
        character (len=1)                   :: facformat,arealformat,alimit
        character (len=20)                  :: atemp,aline
        character (len=200)                 :: aprompt,facfile,afile,arealfile, &
                                               propfile1,propfile2,atemp1
        character (len=12), allocatable, dimension(:)    :: wpoints


	write(amessage,5)
5	format(' Program FAC2FEM carries out spatial interpolation to a MICROFEM ',  &
        'input file based on interpolation factors calculated by PPK2FACF and ', &
        'pilot point values contained in a pilot points file.')
	call write_message(leadspace='yes',endspace='yes')


! -- The first two lines of the interpolation factor file are read.

10      aprompt=' Enter name of interpolation factor file: '
        call open_input_file(ifail,aprompt,facfile,facunit,form_prompt='yes', &
        fformat=facformat)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) afile
        else
          read(facunit,err=9000,end=9100) afile
        end if
        afile=adjustl(afile)
        nbb=len_trim(afile)
        call getfile (ifail,afile,pilot_points_file,1,nbb)
        if(ifail.ne.0)pilot_points_file=' '
        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) nummesh
        else
          read(facunit,err=9000,end=9100) nummesh
        end if
        allocate(realarray(nummesh),minarray(nummesh),maxarray(nummesh),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run FAC2FEM.')
	  go to 9890
	end if

        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) nnppt
        else
          read(facunit,err=9000,end=9100) nnppt
        end if
        allocate(wpoints(nnppt),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if
        if(facformat.eq.'f')then
          do i=1,nnppt
            read(facunit,'(a)',err=9000,end=9100) wpoints(i)
          end do
        else
          do i=1,nnppt
            read(facunit,err=9000,end=9100) wpoints(i)
          end do
        end if
        do i=1,nnppt
          call casetrans(wpoints(i),'hi')
        end do

! -- The pilot points file is read.

        write(6,*)
80      call read_pilot_points_file(ifail, &
	' Enter name of pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(realarray,minarray,maxarray,wpoints,stat=ierr)
	  if(ierr.ne.0) then
	    write(amessage,57)
57	    format(' Memory management error: cannot continue execution.')
	    go to 9890
	  end if
          close(unit=facunit)
	  write(6,*)
	  go to 10
	end if
        if(nnppt.ne.num_pilot_points)then
          write(amessage,59)
          go to 9890
        else
          do i=1,num_pilot_points
            if(pilot_point_id(i).ne.wpoints(i))then
              write(amessage,59)
59            format(' The pilot points in the pilot points file are not the same, ', &
              'or are not arranged in the same order, as the pilot points in the ', &
              'pilot points file read by PPK2FACF when it calculated the ', &
              'factors contained in the factor file.')
              go to 9890
            end if
          end do
        end if

! -- Interpolation upper and lower limits are supplied.

270       write(6,280,advance='no')
280       format(' Enter lower interpolation limit: ')
          if(key_read(rlo).ne.0) go to 270
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 80
          end if
          minarray=rlo

310       write(6,320,advance='no')
320       format(' Enter upper interpolation limit: ')
          if(key_read(rhi).ne.0) go to 310
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 270
          end if
          maxarray=rhi

!        if(any(minarray.gt.maxarray))then
!          write(amessage,330)
!330       format(' Based on the interpolation upper and lower limits ', &
!          'that you have supplied, there is at least one place within the model ', &
!          'domain where the lower limit exceeds the upper limit.')
!          call write_message(leadspace='yes',endspace='yes')
!          write(6,*)
!          go to 270
!        end if
        if(rlo.ge.rhi)then
          write(amessage,330)
330       format(' Upper interpolation limit does not exceed lower interpolation ', &
          'limit.')
          call write_message()
          go to 270
        end if

! -- The mesh property files are named.

        write(6,*)
1070    aprompt=' Enter name of existing mesh property file: '
        call open_input_file(ifail,aprompt,propfile1,propunit1)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 310
	end if
1080    write(6,1090,advance='no')
1090    format(' Enter property number for replacement: ')
        if(key_read(colrep).ne.0) go to 1080
        if(escset.eq.1) then
          write(6,*)
          escset=0
          close(unit=propunit1)
          go to 1070
        end if
1095    aprompt = ' Enter name for new mesh property file: '
        call open_output_file(ifail,aprompt,propfile2,propunit2)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1)then
          escset=0
          go to 1080
        end if
        atemp1=propfile2
        call casetrans(atemp1,'lo')
        ifem=0
        lt=len_trim(atemp1)
        if (lt.gt.5) then
          if(atemp1(lt-3:lt).eq.'.fem')then
            ifem=1
          end if
        end if
110     write(6,120,advance='no')
120     format(' Enter value for elements to which no interpolation takes place: ')
        if(key_read(backval).ne.0) go to 110
        if(escset.eq.1) then
          write(6,*)
          escset=0
          close(unit=propunit2)
          go to 1095
        end if
        realarray=backval

! -- The interpolation factor file is now read line by line and the factors
!    are used to undertake spatial interpolation.

        allocate(ipt(num_pilot_points),wt(num_pilot_points),stat=ierr)
        if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

        do
          if(facformat.eq.'f')then
            read(facunit,*,err=9000,end=200) icellno,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          else
            read(facunit,err=9000,end=200)   icellno,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          end if
          sum=rtemp
          do i=1,na
            if(itrans.eq.0)then
              sum=sum+pilot_point_val(ipt(i))*wt(i)
            else
              rtemp2=pilot_point_val(ipt(i))
              if(rtemp2.le.0.0)then
                write(amessage,125)
125             format(' The interpolation factor file specifies that spatial ', &
                'interpolation in at least one zone takes place on the basis of the ', &
                'logarithms of pilot point values. However at least one of the ', &
                'pertinent pilot point values is negative.')
                go to 9890
              end if
              sum=sum+log10(rtemp2)*wt(i)
            end if
          end do
          if(itrans.eq.0)then
            rtemp1=sum
          else
            rtemp1=10**sum
          end if
          if(rtemp1.gt.maxarray(icellno))then
            realarray(icellno)=maxarray(icellno)
          else if(rtemp1.lt.minarray(icellno))then
            realarray(icellno)=minarray(icellno)
          else
            realarray(icellno)=rtemp1
          end if
        end do

! -- The input mesh property file is read and the output mesh property file is written.

200     continue
        if(ifem.eq.0)then
          iline=1
          read(propunit1,'(a)',err=9300,end=9350) cline
          write(propunit2,'(a)') trim(cline)
          do i=1,nummesh
            iline=iline+1
            read(propunit1,'(a)',err=9300,end=9350) cline
            call linesplit(ifail,colrep)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,1120) trim(aline),trim(propfile1)
1120          format(' There are insufficient entries on line ',a,' of file ',a)
              go to 9890
            end if
            write(atemp,'(1pg14.7)') realarray(i)
            if(colrep.gt.1)then
              write(propunit2,1140) cline(1:right_word(colrep-1)),  &
              trim(atemp),trim(cline(right_word(colrep)+1:))
1140          format(a,2x,a,2x,a)
            else
              write(propunit2,1150) trim(atemp),trim(cline(right_word(colrep)+1:))
1150          format(a,2x,a)
            end if
          end do
          do
            read(propunit1,'(a)',err=9300,end=1500) cline
            write(propunit2,'(a)') trim(cline)
          end do
        else
          iline=0
          do i=1,19
            iline=iline+1
            read(propunit1,'(a)',err=9300,end=9350) cline
            write(propunit2,'(a)') trim(cline)
          end do
          do i=1,nummesh
            iline=iline+1
            read(propunit1,'(a)',err=9300,end=9350) cline
            write(propunit2,'(a)') trim(cline)
            iline=iline+1
            read(propunit1,'(a)',err=9300,end=9350) cline
            call linesplit(ifail,colrep)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,1120) trim(aline),trim(propfile1)
              go to 9890
            end if
            write(atemp,'(1pg14.7)') realarray(i)
            if(colrep.gt.1)then
              write(propunit2,1140) cline(1:right_word(colrep-1)),  &
              trim(atemp),trim(cline(right_word(colrep)+1:))
            else
              write(propunit2,1150) trim(atemp),trim(cline(right_word(colrep)+1:))
            end if
            iline=iline+1
            read(propunit1,'(a)',err=9300,end=9350) cline
            write(propunit2,'(a)') trim(cline)
          end do
          do
            read(propunit1,'(a)',err=9300,end=1500) cline
            write(propunit2,'(a)') trim(cline)
          end do
        end if

1500    write(6,1510) trim(propfile2)
1510    format(' - file ',a,' written ok.')

        go to 9900


9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a)
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to file ',a)
        go to 9890
9300    call num2char(iline,aline)
        write(amessage,9310) trim(propfile1)
9310    format(' Error reading line ',a,' of file ',a)
        go to 9890
9350    write(amessage,9360) trim(propfile1)
9360    format(' Unexepected end encountered to file ',a)
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(realarray,minarray,maxarray,ipt,wt,wpoints,stat=ierr)

end program fac2fem


