program fac2real3d

! -- Program FAC2REAL3D writes a series of real arrays based on the interpolation factor
!    file produced by PPK2FAC3D.

	use defn
	use inter

	implicit none

        integer               :: ifail,idate,iheader,facunit,ierr,ncol,nrow,nerr,nbb, &
                                 icellno,na,icol,irow,i,itrans,nnppt
        integer               :: nlay,nrc,ilay,iunit
        integer, allocatable, dimension(:)  :: ipt
        real                                :: backval,rtemp,sum,rtemp1,rlo,rhi,rtemp2
        real, allocatable, dimension(:)     :: wt
        real, allocatable, dimension(:,:,:) :: realarray
        character (len=1)                   :: facformat,arealformat,alimit
        character (len=1)                   :: ao
        character (len=12)                  :: alay
        character (len=200)                 :: aprompt,facfile,atempf
        character (len=200)                 :: intbasename,outbasename,outfile
        character (len=12), allocatable, dimension(:)    :: wpoints


	write(amessage,5)
5	format(' Program FAC2REAL3D carries out 3D spatial interpolation based on ', &
        'interpolation factors calculated by PPK2FAC3D and pilot point values contained ', &
        'in a 3D pilot points file.')
	call write_message(leadspace='yes',endspace='yes')

! -- The settings file is read.

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  datespec=1
	  headerspec='no'
	  go to 1
!	  write(amessage,7)
!7	  format(' A settings file (settings.fig) was not found in the ', &
!	  'current directory.')
!	  call write_message
!	  go to 9900
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
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
1       continue

! -- The first two lines of the interpolation factor file are read.

10      aprompt=' Enter name of interpolation factor file: '
        call open_input_file(ifail,aprompt,facfile,facunit,form_prompt='yes', &
        fformat=facformat)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) atempf
        else
          read(facunit,err=9000,end=9100) atempf
        end if
        atempf=adjustl(atempf)
        nbb=len_trim(atempf)
        call getfile (ifail,atempf,pilot_points_file,1,nbb)
        if(ifail.ne.0)pilot_points_file=' '

        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) atempf
        else
          read(facunit,err=9000,end=9100) atempf
        end if
        atempf=adjustl(atempf)
        nbb=len_trim(atempf)
        call getfile (ifail,atempf,intbasename,1,nbb)
        if(ifail.ne.0)intbasename=' '

        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) ncol,nrow,nlay
        else
          read(facunit,err=9000,end=9100) ncol,nrow,nlay
        end if
        allocate(realarray(ncol,nrow,nlay),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run FAC2REAL3D.')
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
80      call read_3d_pilot_points_file(ifail, &
	' Enter name of 3D pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  deallocate(realarray,wpoints,stat=ierr)
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
              'pilot points file read by PPK2FAC3D when it calculated the ', &
              'factors contained in the factor file.')
              go to 9890
            end if
          end do
        end if

! -- Interpolation upper and lower limits are supplied.

        write(6,*)
270     write(6,280,advance='no')
280     format(' Enter lower interpolation limit: ')
        if(key_read(rlo).ne.0) go to 270
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 80
        end if

310     write(6,320,advance='no')
320     format(' Enter upper interpolation limit: ')
        if(key_read(rhi).ne.0) go to 310
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 270
        end if
        if(rlo.ge.rhi)then
          write(6,330)
330       format(/,' *** Must be greater than lower interpolation limit - try again ***',/)
          go to 310
        end if

! -- The output real array file is named.

        write(6,*)
91      write(6,92,advance='no')
92      format(' Write outputs to single 3D table or multiple 2D real array files?  [s/m]: ')
        read(5,'(a)') ao
        call casetrans(ao,'lo')
        if(ao.eq.'e')then
          write(6,*)
          go to 310
        end if
        if((ao.ne.'s').and.(ao.ne.'m')) go to 91
100     continue
        if(ao.eq.'m')then
          write(6,105,advance='no')
105       format(' Enter filename base for output real array files: ')
        else
          write(6,106,advance='no')
106       format(' Enter name for output file: ')
        end if
        read(5,'(a)') atempf
        if(atempf.eq.' ') go to 100
        atempf=adjustl(atempf)
        if((atempf(1:2).eq.'E ').or.(atempf(1:2).eq.'e '))then
          write(6,*)
          go to 91
        end if
        nbb=len_trim(atempf)
        call getfile(ifail,atempf,outbasename,1,nbb)
        if(ifail.ne.0) go to 100

110     write(6,120,advance='no')
120     format(' Enter value for elements to which no interpolation takes place: ')
        if(key_read(backval).ne.0) go to 110
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 100
        end if
        realarray=backval            ! An array

! -- The interpolation factor file is now read line by line and the factors
!    are used to undertake spatial interpolation.

        allocate(ipt(num_pilot_points),wt(num_pilot_points),stat=ierr)
        if(ierr.ne.0) then
	  write(amessage,50)
	  go to 9890
	end if

        nrc=ncol*nrow
        do
          if(facformat.eq.'f')then
            read(facunit,*,err=9000,end=200) icellno,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          else
            read(facunit,err=9000,end=200)   icellno,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          end if
          ilay=(icellno-1)/nrc+1
          icellno=icellno-(ilay-1)*nrc
          irow=(icellno-1)/ncol+1
          icol=icellno-((irow-1)*ncol)
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
              sum=sum+log(rtemp2)*wt(i)
            end if
          end do
          if(itrans.eq.0)then
            rtemp1=sum
          else
            rtemp1=exp(sum)
          end if
          if(rtemp1.gt.rhi)then
            realarray(icol,irow,ilay)=rhi
          else if(rtemp1.lt.rlo)then
            realarray(icol,irow,ilay)=rlo
          else
            realarray(icol,irow,ilay)=rtemp1
          end if
        end do

200     continue

! -- The output file(s) are now written.

        if(ao.eq.'m')then
          do ilay=1,nlay
            call num2char(ilay,alay)
            outfile=trim(outbasename)//trim(alay)//'.ref'
            iunit=nextunit()
            open(unit=iunit,file=outfile)
            if(headerspec.eq.'yes') then
              write(iunit,220) ncol,nrow
220           format(2i10)
            end if
            do irow=1,nrow
              write(iunit,240) (realarray(icol,irow,ilay),icol=1,ncol)
240           format(7(1x,1pg14.7))
            end do
            close(unit=iunit)
            write(6,250) trim(outfile)
250         format(' - file ',a,' written ok.')
          end do
        else
          outfile=outbasename
          iunit=nextunit()
          open(unit=iunit,file=outfile)
          write(iunit,251) ncol,nrow,nlay
251       format(3i10)
          do ilay=1,nlay
            do irow=1,nrow
              do icol=1,ncol
                write(iunit,216) realarray(icol,irow,ilay)
216             format(1pg14.7)
              end do
            end do
          end do
          close(unit=iunit)
          write(6,250) trim(outfile)
        end if

        go to 9900

9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a)
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to file ',a)
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(realarray,ipt,wt,wpoints,stat=ierr)

end program fac2real3d


