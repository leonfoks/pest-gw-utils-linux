!     Last change:  JD   12 Feb 2001    2:26 pm

program fac2rsm

! -- Program FAC2RSM writes an RSM input data file for a two dimensional mesh based on interpolation
!    factors computed by PPK2FACR.

        use defn
        use inter

        implicit none

        integer               :: ifail,facunit,ierr,ncol,nrow,nerr,nbb, &
                                 icellno,na,icol,irow,i,itrans,nnppt,outunit,ielem
        integer, allocatable, dimension(:)  :: ipt
        real                                :: backval,rtemp,sum,rtemp1,rlo,rhi,rtemp2
        real, allocatable, dimension(:)     :: wt
        real, allocatable, dimension(:)     :: realarray
        character (len=1)                   :: facformat,arealformat,alimit
        character (len=10)                  :: atemp
        character (len=20)                  :: datname
        character (len=120)                 :: aprompt,facfile,afile,intfile,outfile
        character (len=12), allocatable, dimension(:)    :: wpoints


        write(amessage,5)
5        format(' Program FAC2RSM carries out spatial interpolation based on ', &
        'interpolation factors calculated by PPK2FACR and pilot point values contained ', &
        'in a pilot points file.')
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
          read(facunit,'(a)',err=9000,end=9100) afile
        else
          read(facunit,err=9000,end=9100) afile
        end if
        afile=adjustl(afile)
        if(afile.eq.' ')then
          intfile=' '
        else
          nbb=len_trim(afile)
          call getfile (ifail,afile,intfile,1,nbb)
          if(ifail.ne.0)intfile=' '
        end if

        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) numelem_g
        else
          read(facunit,err=9000,end=9100) numelem_g
        end if
        allocate(realarray(numelem_g),stat=ierr)
        if(ierr.ne.0) then
          write(amessage,50)
50        format(' Cannot allocate sufficient memory to run FAC2RSM.')
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
          wpoints(i)=adjustl(wpoints(i))
          call casetrans(wpoints(i),'hi')
        end do

! -- The pilot points file is read.

        write(6,*)
80      call read_pilot_points_file(ifail, &
        ' Enter name of pilot points file: ')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          deallocate(realarray,wpoints,stat=ierr)
          if(ierr.ne.0) then
            write(amessage,57)
57          format(' Memory management error: cannot continue execution.')
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
              'pilot points file read by PPK2FACR when it calculated the ', &
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
330       format(/,' Upper limit must be greater than lower limit - try again.',/)
          go to 310
        end if

! -- The output real array file is named.

        write(6,*)
350     call open_output_file(ifail, &
        ' Enter name for mesh property file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 310
        end if
370     write(6,380,advance='no')
380     format(' Enter NAME of data type: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 370
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          close(unit=outunit)
          write(6,*)
          go to 350
        end if
        afile=adjustl(afile)
        nbb=len_trim(afile)
        call getfile (ifail,afile,datname,1,nbb)
        if(ifail.ne.0) go to 370
410     write(6,420,advance='no')
420     format(' Enter value for elements to which no interpolation takes place: ')
        if(key_read(backval).ne.0) go to 410
        if(escset.eq.1) then
          write(6,*)
          escset=0
          go to 370
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
            read(facunit,*,err=9000,end=600) icellno,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          else
            read(facunit,err=9000,end=600)   icellno,itrans, &
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
          if(rtemp1.gt.rhi)then
            rtemp1=rhi
          else if(rtemp1.lt.rlo)then
            rtemp1=rlo
          end if
          realarray(icellno)=rtemp1
        end do
600     continue

! -- The output file is now written.

        write(outunit,510)
510     format('DATASET')
        write(outunit,520)
520     format('OBJTYPE "network"')
        write(outunit,530)
530     format('BEGSCL')
        call num2char(numelem_g,atemp)
        write(outunit,540) trim(atemp)
540     format('ND ',a)
        write(outunit,550) trim(datname)
550     format('NAME "',a,'"')
        write(outunit,560)
560     format('TS 0 0.0')
        do ielem=1,numelem_g
          write(outunit,570) realarray(ielem)
570       format(1x,1pg14.7)
        end do
        close(unit=outunit)
        write(6,580) trim(outfile)
580     format(' - file ',a,' written ok.')


        go to 9900


9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a)
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to file ',a)
        go to 9890

9890    call write_message(leadspace='yes',endspace='yes')
9900    call close_files
        deallocate(realarray,ipt,wt,wpoints,stat=ierr)

end program fac2rsm



