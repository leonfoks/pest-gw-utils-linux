!     Last change:  JD   12 Feb 2001    2:26 pm

program fac2real

! -- Program FAC2REAL writes a real array based on the interpolation factor
!    file produced by PPK2FAC.

	use defn
	use inter

	implicit none

        integer               :: ifail,idate,iheader,facunit,ierr,ncol,nrow,nerr,nbb, &
                                 icellno,na,icol,irow,i,itrans,nnppt
        integer, allocatable, dimension(:)  :: ipt
        real                                :: backval,rtemp,sum,rtemp1,rlo,rhi,rtemp2
        real, allocatable, dimension(:)     :: wt
        real, allocatable, dimension(:,:)   :: realarray,minarray,maxarray
        character (len=1)                   :: facformat,arealformat,alimit
        character (len=200)                 :: aprompt,facfile,afile,arealfile,intfile
        character (len=12), allocatable, dimension(:)    :: wpoints


	write(amessage,5)
5	format(' Program FAC2REAL carries out spatial interpolation based on ', &
        'interpolation factors calculated by PPK2FAC and pilot point values contained ', &
        'in a pilot points file.')
	call write_message(leadspace='yes',endspace='yes')

! -- The settings file is read.

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
	  call write_message
	  go to 9900
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
          read(facunit,'(a)',err=9000,end=9100) afile
        else
          read(facunit,err=9000,end=9100) afile
        end if
        afile=adjustl(afile)
        nbb=len_trim(afile)
        call getfile (ifail,afile,intfile,1,nbb)
        if(ifail.ne.0)intfile=' '

        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) ncol,nrow
        else
          read(facunit,err=9000,end=9100) ncol,nrow
        end if
        allocate(realarray(ncol,nrow),minarray(ncol,nrow),maxarray(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run FAC2REAL.')
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
              'pilot points file read by PPK2FAC when it calculated the ', &
              'factors contained in the factor file.')
              go to 9890
            end if
          end do
        end if

! -- Interpolation upper and lower limits are supplied.

        write(6,*)
255     write(6,260,advance='no')
260     format(' Supply lower interpolation limit as an array or single value? [a/s]: ')
        read(5,'(a)') alimit
        if((alimit.eq.'e').or.(alimit.eq.'E')) then
          write(6,*)
          go to 80
        end if
        call casetrans(alimit,'lo')
        if((alimit.ne.'a').and.(alimit.ne.'s')) go to 255
        if(alimit.eq.'s')then
270       write(6,280,advance='no')
280       format(' Enter lower interpolation limit: ')
          if(key_read(rlo).ne.0) go to 270
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 255
          end if
          minarray=rlo
        else
 	  aprompt=' Enter name of lower interpolation limit array file: '
	  call read_real_array(ifail,aprompt,minarray,pm_header=headerspec, &
          rows=nrow,columns=ncol)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    write(6,*)
	    go to 255
	  end if
        end if

        write(6,*)
290     write(6,300,advance='no')
300     format(' Supply upper interpolation limit as an array or single value? [a/s]: ')
        read(5,'(a)') alimit
        if((alimit.eq.'e').or.(alimit.eq.'E'))then
          write(6,*)
          go to 255
        end if
        call casetrans(alimit,'lo')
        if((alimit.ne.'a').and.(alimit.ne.'s')) go to 290
        if(alimit.eq.'s')then
310       write(6,320,advance='no')
320       format(' Enter upper interpolation limit: ')
          if(key_read(rhi).ne.0) go to 310
          if(escset.eq.1) then
            write(6,*)
            escset=0
            go to 290
          end if
          maxarray=rhi
        else
 	  aprompt=' Enter name of upper interpolation limit array file: '
	  call read_real_array(ifail,aprompt,maxarray,pm_header=headerspec, &
          rows=nrow,columns=ncol)
	  if(ifail.ne.0) go to 9900
	  if(escset.eq.1) then
	    escset=0
	    write(6,*)
	    go to 290
	  end if
        end if

        if(any(minarray.gt.maxarray))then
          write(amessage,330)
330       format(' Based on the interpolation upper and lower limits or arrays ', &
          'that you have supplied, there is at least one place within the model ', &
          'domain where the lower limit exceeds the upper limit.')
          call write_message(leadspace='yes',endspace='yes')
          write(6,*)
          go to 255
        end if

! -- The output real array file is named.

        write(6,*)
100     aprompt=' Enter name for output real array file: '
        call write_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow,columns=ncol,istage=1,realfile=arealfile,aaformat=arealformat)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 290
        end if
110     write(6,120,advance='no')
120     format(' Enter value for elements to which no interpolation takes place: ')
        if(key_read(backval).ne.0) go to 110
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 100
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
              sum=sum+log10(rtemp2)*wt(i)
            end if
          end do
          if(itrans.eq.0)then
            rtemp1=sum
          else
            rtemp1=10**sum
          end if
          if(rtemp1.gt.maxarray(icol,irow))then
            realarray(icol,irow)=maxarray(icol,irow)
          else if(rtemp1.lt.minarray(icol,irow))then
            realarray(icol,irow)=minarray(icol,irow)
          else
            realarray(icol,irow)=rtemp1
          end if
        end do


200     continue
        call write_real_array(ifail,aprompt,realarray,pm_header=headerspec, &
        rows=nrow,columns=ncol,istage=2,realfile=arealfile,aaformat=arealformat)


        go to 9900


9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a)
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to file ',a)
        go to 9890

9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(realarray,minarray,maxarray,ipt,wt,wpoints,stat=ierr)

end program fac2real


