program rdat2tab

! -- Program RDAT2TAB reads a mesh data file. It assigns coordinates to data elements
!    and writes these, together with the data, in tabular format.

        use defn
        use inter

        implicit none

        integer               :: ifail,ierr,outunit,i
        real, allocatable     :: rdata(:)
        character (len=10)    :: atemp
        character (len=20)    :: atemp1
        character (len=200)   :: aprompt,outfile


        !open(unit=*,action='read',carriagecontrol='list')

        write(amessage,5)
5       format(' Program RDAT2TAB assigns coordinates to RSM element data and ',   &
        'writes this data in tabular form.')

        call write_message(leadspace='yes',endspace='yes')

! -- The GMS 2D mesh file is read.

10      continue
        call read_gms_2d_mesh_file(ifail,' Enter name of GMS two-dimensional mesh file: ','yes')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900
        allocate(rdata(numelem_g),stat=ierr)
        if(ierr.ne.0) go to 9200

        write(6,*)
50      aprompt=' Enter name of mesh element data file: '
        call read_gms_real_mesh_data_file(ifail,aprompt,rdata)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) then
          escset=0
          write(6,*)
          deallocate(eastnode_g,northnode_g,nodenum_g,elemnum_g,elemindex_g,elemnode_g,  &
                     eastelem_g,northelem_g,rdata,stat=ierr)
          if(ierr.ne.0) then
            write(amessage,131)
131         format(' Memory management error: cannot continue execution.')
            go to 9890
          end if
          go to 10
        end if

        write(6,*)
160     call open_output_file(ifail, &
        ' Enter name for tabular data output file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 50
        end if

        write(outunit,170)
170     format(' Element_id',t16,'Easting',t31,'Northing',t45,'Data_value',t60,'Log_data_value')
        do i=1,numelem_g
          call num2char(elemnum_g(i),atemp)
          if(rdata(i).le.0.0)then
            atemp1='   - '
          else
            write(atemp1,'(1pg14.7)') log10(rdata(i))
          end if
          write(outunit,180) trim(atemp),eastelem_g(i),northelem_g(i),rdata(i),trim(atemp1)
180       format(2x,a,t15,f13.3,t30,f13.3,t45,1pg14.7,t60,a)
        end do
        close(unit=outunit)
        write(6,200) trim(outfile)
200     format('  - file ',a,' written ok.')

        go to 9900

9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890

9890    call write_message(leadspace='yes',endspace='yes')
9900    call close_files

        deallocate(nodenum_g,eastnode_g,northnode_g,stat=ierr)
        deallocate(elemnum_g,elemnode_g,eastelem_g,northelem_g,elemindex_g,stat=ierr)

end program rdat2tab



subroutine whichone_i(ifail,npar,ipar,ival,jval)

! -- Subroutine whichone_i locates an integer value in an array of integers.

        integer npar,ipar,i
        integer ifail
        integer jval
        integer ival(npar)

        ifail=0
        if((ipar.lt.1).or.(ipar.gt.npar)) ipar=1
        if(jval.eq.ival(ipar)) return
        if(ipar.ne.npar)then
          do 20 i=ipar+1,npar
          if(jval.eq.ival(i))then
            ipar=i
            return
          end if
20        continue
        end if
        if(ipar.ne.1)then
          do 40 i=ipar-1,1,-1
          if(jval.eq.ival(i)) then
            ipar=i
            return
          end if
40        continue
        end if
        ifail=1
        return

 end subroutine whichone_i


