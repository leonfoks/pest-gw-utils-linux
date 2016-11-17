subroutine read_gms_2d_mesh_file(ifail,aprompt,centroid)

! -- Subroutine READ_GMS_2D_MESH_FILE reads a GMS mesh file.

        use defn
        use inter
        implicit none

        integer, intent(out)            :: ifail
        character (len=*), intent(in)   :: aprompt
        character (len=*), intent(in)   :: centroid

        integer                         :: iline,iunit,iend,ierr,inode,icentroid,itemp,ielem,i, &
                                           jtemp,kelem,jelem
        double precision                :: dtemp1,dtemp2
        character*10                    :: aline
        character*5                     :: atemp
        character*6                     :: atemp1,atemp2

! -- The CENTROID argument is checked.

        atemp=centroid
        call casetrans(atemp,'lo')
        if(atemp.eq.'yes')then
          icentroid=1
        else if(atemp.eq.'no')then
          icentroid=0
        else
          write(amessage,5)
5         format(' CENTROID argument in subroutine READ_GMS_2D_MESH_FILE must be "yes" or "no".')
	  call write_message(increment=1,leadspace='yes')
          ifail=1
          return
        end if

! -- Initialisation.

        iline=0
        ifail=0
        imessage=0
        numnode_g=0
        numelem_g=0

! -- The name of the GMS mesh file is acquired.

10      write(6,'(a)',advance='no') aprompt(1:len_trim(aprompt)+1)
        read(5,'(a)') cline
        if(cline.eq.' ') go to 10
        cline=adjustl(cline)
        if(index(eschar,cline(1:2)).ne.0) then
          escset=1
          return
        end if
        iend=len_trim(cline)
        call getfile(ifail,cline,meshfile_g,1,iend)
        if(ifail.ne.0) go to 10

! -- The GMS mesh file is opened.

	iunit=nextunit()
	open(unit=iunit,file=meshfile_g,status='old',iostat=ierr)
	if(ierr.ne.0)then
	  if(imessage.eq.5) go to 9890
	  write(amessage,20) trim(meshfile_g)
20        format('cannot open GMS mesh file ',a,' - try again.')
	  call write_message(error='yes',increment=1)
	  go to 10
	end if
	write(initial_message,30) trim(meshfile_g)
30      format(' Errors in mesh file ',a,' ----->')
        imessage=0

! -- The first line of the file is checked for the correct header.

        read(iunit,'(a)',err=9000,end=9100) cline
        cline=adjustl(cline)
        call casetrans(cline,'hi')
        if(cline(1:6).ne.'MESH2D')then
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          write(amessage,35) trim(meshfile_g)
35        format(' First line of file ',a,' expected to contain "MESH2D" header.')
          call write_message(increment=1)
          go to 9890
        end if

! -- Nodes and elements are counted.

        do
          iline=iline+1
          read(iunit,'(a)',err=9000,end=100) cline
          if(cline.eq.' ') cycle
          cline=adjustl(cline)
          call linesplit(ifail,1)
          atemp=cline(left_word(1):right_word(1))
          call casetrans(atemp,'hi')
          if(atemp.eq.'ND') then
            numnode_g=numnode_g+1
          else if(atemp.eq.'E3T')then
            numelem_g=numelem_g+1
          end if
        end do

100     continue
        if(numnode_g.eq.0)then
          write(amessage,110) trim(meshfile_g)
110       format(' No nodes were found in file ',a,'.')
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          call write_message(increment=1)
        end if
        if(numelem_g.eq.0)then
          write(amessage,112) trim(meshfile_g)
112       format(' No elements were found in file ',a,'.')
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          call write_message(increment=1)
        end if
        if(imessage.ne.0) go to 9890

! -- The file is rewound.

        rewind(unit=iunit)

! -- Memory is allocated.

        allocate(nodenum_g(numnode_g),eastnode_g(numnode_g),northnode_g(numnode_g),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,111)
111       format(' Cannot allocate sufficient memory to continue execution.')
          call write_message(increment=1,leadspace='yes')
          go to 9890
        end if

        allocate(elemnum_g(numelem_g),elemnode_g(numelem_g,3),elemindex_g(numelem_g),stat=ierr)
        if(ierr.ne.0)then
          write(amessage,111)
          call write_message(increment=1,leadspace='yes')
          go to 9890
        end if

! -- Node and element data are now read.

        inode=0
        ielem=0
        iline=0
        do
          iline=iline+1
          call num2char(iline,aline)
          read(iunit,'(a)',err=9000,end=200) cline
          call linesplit(ifail,1)
          atemp=cline(left_word(1):right_word(1))
          call casetrans(atemp,'hi')
          if(atemp.eq.'ND') then
            inode=inode+1
            call linesplit(ifail,4)
            if(ifail.ne.0)then
              write(amessage,140) trim(aline),trim(meshfile_g)
140           format(' Insufficent items on line ',a,' of file ',a,'.')
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              call write_message(increment=1)
            else
              nodenum_g(inode)=char2int(ifail,2)
              if(ifail.ne.0)then
                write(amessage,150) trim(aline),trim(meshfile_g)
150             format(' Cannot read node number from second data column at line ',a,   &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
              eastnode_g(inode)=char2double(ifail,3)
              if(ifail.ne.0)then
                write(amessage,160) trim(aline),trim(meshfile_g)
160             format(' Cannot read easting from third data column at line ',a,    &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
              northnode_g(inode)=char2double(ifail,4)
              if(ifail.ne.0)then
                write(amessage,170) trim(aline),trim(meshfile_g)
170             format(' Cannot read northing from fourth data column at line ',a,    &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
            end if
          else if(atemp.eq.'E3T') then
            ielem=ielem+1
            call linesplit(ifail,5)
            if(ifail.ne.0)then
              write(amessage,140) trim(aline),trim(meshfile_g)
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              call write_message(increment=1)
            else
              elemnum_g(ielem)=char2int(ifail,2)
              if(ifail.ne.0)then
                write(amessage,250) trim(aline),trim(meshfile_g)
250             format(' Cannot read element number from second data column at line ',a,   &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
              elemnode_g(ielem,1)=char2int(ifail,3)
              if(ifail.ne.0)then
                write(amessage,260) trim(aline),trim(meshfile_g)
260             format(' Cannot read node number from third data column at line ',a,    &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
              elemnode_g(ielem,2)=char2int(ifail,4)
              if(ifail.ne.0)then
                write(amessage,270) trim(aline),trim(meshfile_g)
270             format(' Cannot read node number from fourth data column at line ',a,    &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
              elemnode_g(ielem,3)=char2int(ifail,5)
              if(ifail.ne.0)then
                write(amessage,280) trim(aline),trim(meshfile_g)
280             format(' Cannot read node number from fifth data column at line ',a,    &
                ' of file ',a,'.')
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                call write_message(increment=1)
              end if
            end if
          end if
        end do
200     continue

! -- Mesh element centroid coordinates are now evaluated.

        if(icentroid.eq.1)then
          allocate(eastelem_g(numelem_g),northelem_g(numelem_g),stat=ierr)
          if(ierr.ne.0)then
            write(amessage,111)
            call write_message(increment=1,leadspace='yes')
            go to 9890
          end if

          itemp=0
          do ielem=1,numelem_g
            dtemp1=0.0d0
            dtemp2=0.0d0
            do i=1,3
              call whichone_i(ifail,numnode_g,itemp,nodenum_g,elemnode_g(ielem,i))
              if(ifail.ne.0)then
                call num2char(elemnum_g(ielem),atemp1)
                call num2char(elemnode_g(ielem,i),atemp2)
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,420) trim(atemp2),trim(atemp1),trim(meshfile_g)
420             format(' Cannot find node number ',a,' cited as belonging to cell number ',  &
                a,' in file ',a,'.')
                call write_message(increment=1)
              end if
              dtemp1=dtemp1+eastnode_g(itemp)
              dtemp2=dtemp2+northnode_g(itemp)
            end do
            eastelem_g(ielem)=dtemp1/3.0d0
            northelem_g(ielem)=dtemp2/3.0d0
          end do
        end if

! -- A sorting index array of element node numbers is now built.

        jtemp=-huge(jtemp)
        do jelem=1,numelem_g
          itemp=huge(itemp)
          do ielem=1,numelem_g
            if(elemnum_g(ielem).lt.itemp)then
              if(elemnum_g(ielem).gt.jtemp)then
                itemp=elemnum_g(ielem)
                kelem=ielem
              end if
            end if
          end do
          elemindex_g(jelem)=kelem
          jtemp=itemp
        end do

! -- Wrapping up.

        if(imessage.ne.0)then
          go to 9890
        else
          call num2char(numnode_g,atemp1)
          call num2char(numelem_g,atemp2)
          write(amessage,220) trim(atemp2),trim(atemp1),trim(meshfile_g)
220       format('  - ',a,' elements and ',a,' nodes cited in file ',a)
          call write_message
          go to 9900
        end if

9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(meshfile_g)
9010    format(' Error encountered when reading line ',a,' of file ',a,'.')
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        call write_message(increment=1)
        go to 9890

9100    write(amessage,9110) trim(meshfile_g)
9110    format(' Unexpected end encountered to file ',a,'.')
        if(imessage.eq.0) call write_initial_message(leadspace='yes')
        call write_message(increment=1)
        go to 9890

9890    ifail=1

9900    continue
        close(unit=iunit,iostat=ierr)
        imessage=0
        return

end subroutine read_gms_2d_mesh_file




