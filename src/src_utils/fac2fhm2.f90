!     Last change:  JD   28 Dec 2002   10:39 pm

program fac2fhm2

! -- Program FAC2FHM2 writes mesh data to a set of near-FEHM-compatible files
!    based on spatial interpolation from pilot points using factors calculated
!    by PPKFACM.

	use defn
	use inter

	implicit none

        integer, parameter    :: MAXNODE=2000
        integer               :: ifail,facunit,ierr,nerr,nbb, &
                                 icellno,na,i,itrans,nnppt,colrep, &
                                 nummesh,iline,ifem,lt
        integer               :: imat,outunit,femunit,numnode,nnode
        integer                             :: node(MAXNODE)
        integer, allocatable, dimension(:)  :: ipt
        real                                :: rtemp,sum,rtemp1,rlo,rhi,rtemp2,backval,sum1, &
                                               rlo1,rhi1
        real                                :: minarray,maxarray,minarray1,maxarray1
        real, allocatable, dimension(:)     :: wt
        real, allocatable, dimension(:)     :: realarray,realarray1
        double precision                    :: e,n
        character (len=1)                   :: facformat,arealformat,alimit
        character (len=20)                  :: atemp,aline
        character (len=120)                 :: aprompt,facfile,afile,arealfile, &
                                               atemp1,outfile
        character (len=120)                    femfile
        character (len=10), allocatable, dimension(:)    :: wpoints


	write(amessage,5)
5	format(' Program FAC2FHM2 carries out spatial interpolation for FEHM ',  &
        'input based on interpolation factors calculated by PPK2FACM and ', &
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
          read(facunit,*,err=9000,end=9100) nummesh,imat
        else
          read(facunit,err=9000,end=9100) nummesh,imat
        end if
        allocate(realarray(nummesh),realarray1(nummesh),stat=ierr)
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
        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) femfile
        else
          read(facunit,err=9000,end=9100) femfile
        end if

! -- The pilot points file is read.

        write(6,*)
80      call read_dual_pilot_points_file(ifail, &
	' Enter name of pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(realarray,wpoints,realarray1,stat=ierr)
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
              'pilot points file read by PPK2FACM when it calculated the ', &
              'factors contained in the factor file.')
              go to 9890
            end if
          end do
        end if

! -- Interpolation upper and lower limits are supplied.

270       write(6,280,advance='no')
280       format(' Enter lower interpolation limit for first data type: ')
          if(key_read(rlo).ne.0) go to 270
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 80
          end if
          minarray=rlo

310       write(6,320,advance='no')
320       format(' Enter upper interpolation limit for first data type: ')
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

271       write(6,281,advance='no')
281       format(' Enter lower interpolation limit for second data type: ')
          if(key_read(rlo1).ne.0) go to 271
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 310
          end if
          minarray1=rlo1

311       write(6,321,advance='no')
321       format(' Enter upper interpolation limit for second data type: ')
          if(key_read(rhi1).ne.0) go to 311
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 271
          end if
          maxarray1=rhi1

!        if(any(minarray.gt.maxarray))then
!          write(amessage,330)
!330       format(' Based on the interpolation upper and lower limits ', &
!          'that you have supplied, there is at least one place within the model ', &
!          'domain where the lower limit exceeds the upper limit.')
!          call write_message(leadspace='yes',endspace='yes')
!          write(6,*)
!          go to 270
!        end if
        if(rlo1.ge.rhi1)then
          write(amessage,331)
331       format(' Upper interpolation limit does not exceed lower interpolation ', &
          'limit.')
          call write_message()
          go to 271
        end if

! -- The output file is named.

        write(6,*)
1070    aprompt=' Enter name for FEHM permeability output file: '
        call open_output_file(ifail,aprompt,outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  go to 311
	end if

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
          sum1=rtemp
          do i=1,na
            if(itrans.eq.0)then
              sum=sum+pilot_point_val(ipt(i))*wt(i)
              sum1=sum1+pilot_point_val1(ipt(i))*wt(i)
            else
              rtemp2=pilot_point_val(ipt(i))
              if(rtemp2.le.0.0)then
                write(amessage,125)
125             format(' The interpolation factor file specifies that spatial ', &
                'interpolation takes place on the basis of the ', &
                'logarithms of pilot point values. However at least one of the ', &
                'pertinent pilot point values is negative.')
                go to 9890
              end if
              sum=sum+log10(rtemp2)*wt(i)
              rtemp2=pilot_point_val1(ipt(i))
              if(rtemp2.le.0.0)then
                write(amessage,125)
                go to 9890
              end if
              sum1=sum1+log10(rtemp2)*wt(i)
            end if
          end do
          if(itrans.eq.0)then
            rtemp1=sum
          else
            rtemp1=10**sum
          end if
          if(rtemp1.gt.maxarray)then
            realarray(icellno)=maxarray
          else if(rtemp1.lt.minarray)then
            realarray(icellno)=minarray
          else
            realarray(icellno)=rtemp1
          end if
          if(itrans.eq.0)then
            rtemp1=sum1
          else
            rtemp1=10**sum1
          end if
          if(rtemp1.gt.maxarray1)then
            realarray1(icellno)=maxarray1
          else if(rtemp1.lt.minarray1)then
            realarray1(icellno)=minarray1
          else
            realarray1(icellno)=rtemp1
          end if
        end do
200     continue
        close(unit=facunit)

! -- Here we must read the xyfile (femfile) and write the output file.

        femunit=nextunit()
        open(unit=femunit,file=femfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,450) trim(femfile)
450       format(' Cannot open MATXY output file ',a,'.')
          go to 9890
        end if
        read(femunit,*,err=9400,end=9450) imat,numnode
        write(outunit,452)
452     format('perm')
        icellno=0
        do
          icellno=icellno+1
          read(femunit,*,err=9400,end=600) e,n,nnode,(node(i),i=1,min(nnode,MAXNODE))
          if(nnode.gt.MAXNODE)then
            write(amessage,451)
451         format(' Increase MAXNODE and re-compile program.')
            go to 9890
          end if
          do i=1,nnode
            write(outunit,480) node(i),node(i),1,realarray(icellno),   &
            realarray(icellno),realarray1(icellno)
480         format(1x,i10,2x,i10,2x,i2,2x,1pg14.7,2x,1pg14.7,2x,1pg14.7)
!jd            write(outunit,*) node(i),e,n,realarray(icellno)
          end do
        end do
600     continue
        write(outunit,601) ' '
601     format(a)

        close(unit=outunit)
1500    write(6,1510) trim(outfile)
1510    format(' - file ',a,' written ok.')

        go to 9900


9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a)
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to file ',a)
        go to 9890
9400    write(amessage,9410) trim(femfile)
9410    format(' Error reading data from MATXY output file ',a,'.')
        go to 9890
9450    write(amessage,9460) trim(femfile)
9460    format(' Unexpected end encountered to MATXY output file ',a,'.')
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(realarray,realarray1,ipt,wt,wpoints,stat=ierr)

end program fac2fhm2


