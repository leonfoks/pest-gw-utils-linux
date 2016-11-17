program rsm2srf

! -- Program RSM2SRF writes RSM mesh data in SURFER-compatible form.

        use defn
        use inter

        implicit none

        integer, parameter     :: MAXCOORD=50000
        integer                :: ifail,ierr,nbb,i,ii,in,ielem,jelem,id1,id2,ic,j,   &
                                  mcoord,inext,ibeg,lastnext
        integer                :: inn(4)

        double precision       :: ecoord(MAXCOORD),ncoord(MAXCOORD)
        character (len=200)    :: afile,outfile


        write(amessage,5)
5       format(' Program RSM2SRF writes RSM mesh data is SURFER-compatible form.')

        call write_message(leadspace='yes',endspace='yes')

! -- The GMS 2D mesh file is read.

10      continue
        call read_gms_2d_mesh_file(ifail,' Enter name of GMS two-dimensional mesh file: ','yes')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

! -- First the nodes coordinates file is written.

        write(6,*)
100     write(6,120,advance='no')
120     format(' Enter name for node coordinates file (<Enter> if none): ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 150
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          deallocate(eastnode_g,northnode_g,nodenum_g,elemnum_g,elemindex_g,elemnode_g,  &
                     eastelem_g,northelem_g,stat=ierr)
          if(ierr.ne.0) then
            write(amessage,131)
131         format(' Memory management error: cannot continue execution.')
            go to 9890
          end if
          write(6,*)
          go to 10
        end if
        nbb=len_trim(afile)
        call getfile(ifail,afile,outfile,1,nbb)
        if(ifail.ne.0) go to 100
        open(unit=20,file=outfile)
        do i=1,numnode_g
          write(20,130) eastnode_g(i),northnode_g(i),nodenum_g(i)
130       format(1x,f12.3,4x,f12.3,4x,i6)
        end do
        close(unit=20)
        write(6,140) trim(outfile)
140     format(' - file ',a,' written ok.')

! -- Next the element centroid coordinates file is written.

150     continue
        write(6,*)
160     write(6,170,advance='no')
170     format(' Enter name for centroid coordinates file (<Enter> if none): ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 200
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          write(6,*)
          go to 100
        end if
        nbb=len_trim(afile)
        call getfile(ifail,afile,outfile,1,nbb)
        if(ifail.ne.0) go to 160
        open(unit=20,file=outfile)
        do i=1,numelem_g
          write(20,130) eastelem_g(i),northelem_g(i),elemnum_g(i)
        end do
        close(unit=20)
        write(6,140) trim(outfile)

! -- Now the mesh bln file is written.

200     continue
        write(6,*)
210     write(6,220,advance='no')
220     format(' Enter name for mesh BLN file (<Enter> if none): ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 300
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          write(6,*)
          go to 160
        end if
        nbb=len_trim(afile)
        call getfile(ifail,afile,outfile,1,nbb)
        if(ifail.ne.0) go to 210
        open(unit=20,file=outfile)
        ii=0
        do ielem=1,numelem_g
          do i=1,3
            in=elemnode_g(ielem,i)
            call whichone_i(ifail,numnode_g,ii,nodenum_g,in)
            ecoord(i)=eastnode_g(ii)
            ncoord(i)=northnode_g(ii)
          end do
          write(20,240) 4,0
240       format(2i6)
          do i=1,3
            write(20,260) ecoord(i),ncoord(i)
260         format(1x,f12.3,2x,f12.3)
          end do
          write(20,260) ecoord(1),ncoord(1)
        end do
        close(unit=20)
        write(6,140) trim(outfile)

! -- Now the mesh boundary bln file is written.

300     continue
        write(6,*)
310     write(6,320,advance='no')
320     format(' Enter name for mesh boundary BLN file (<Enter> if none): ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 600
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          write(6,*)
          go to 210
        end if
        nbb=len_trim(afile)
        call getfile(ifail,afile,outfile,1,nbb)
        if(ifail.ne.0) go to 310
        open(unit=20,file=outfile)

! -- First we find a line that is only used by one element.

        do ielem=1,numelem_g
          inn(1)=elemnode_g(ielem,1)
          inn(2)=elemnode_g(ielem,2)
          inn(3)=elemnode_g(ielem,3)
          inn(4)=inn(1)
          do i=1,3
            id1=inn(i)
            id2=inn(i+1)
            do jelem=1,numelem_g
              if(jelem.eq.ielem) cycle
              ic=0
              do j=1,3
                if(elemnode_g(jelem,j).eq.id1) ic=ic+1
                if(elemnode_g(jelem,j).eq.id2) ic=ic+1
              end do
              if(ic.eq.2) go to 365
            end do
            go to 400
365         continue
          end do
        end do
        write(6,380)
380     format(/,' *** Programming error type 1 ***',/)
        stop

400     continue
        in=id1
        call whichone_i(ifail,numnode_g,ii,nodenum_g,in)
        ecoord(1)=eastnode_g(ii)
        ncoord(1)=northnode_g(ii)
        in=id2
        call whichone_i(ifail,numnode_g,ii,nodenum_g,in)
        ecoord(2)=eastnode_g(ii)
        ncoord(2)=northnode_g(ii)
        mcoord=2

! -- We now have the first two coordinates of our line. Lets find the next.

! -- First we find an element that references the second node and that is part of
!    a line that is only in one element

        inext=id2
        ibeg=id1
        lastnext=ibeg
401     continue
        do ielem=1,numelem_g
          inn(1)=elemnode_g(ielem,1)
          inn(2)=elemnode_g(ielem,2)
          inn(3)=elemnode_g(ielem,3)
          inn(4)=inn(1)
          if((inn(1).ne.inext).and.(inn(2).ne.inext).and.(inn(3).ne.inext)) cycle
          do i=1,3
            id1=inn(i)
            id2=inn(i+1)
            if((id1.ne.inext).and.(id2.ne.inext)) cycle
            if((id1.eq.lastnext).or.(id2.eq.lastnext))cycle
            do jelem=1,numelem_g
              if(jelem.eq.ielem) cycle
              ic=0
              do j=1,3
                if(elemnode_g(jelem,j).eq.id1) ic=ic+1
                if(elemnode_g(jelem,j).eq.id2) ic=ic+1
              end do
              if(ic.eq.2) go to 366
            end do
            go to 470
366         continue
          end do
        end do
        write(6,367)
367     format(/,' *** Programming error type 2 ***',/)
        stop

470     continue
        lastnext=inext
        if(id1.eq.inext)then
          inext=id2
        else
          inext=id1
        end if
        in=inext
        call whichone_i(ifail,numnode_g,ii,nodenum_g,in)
        mcoord=mcoord+1
        if(mcoord.gt.MAXCOORD)then
          write(amessage,480)
480       format(' Mesh is too big - increase MAXCOORD and re-compile program.')
          go to 9890
        end if
        ecoord(mcoord)=eastnode_g(ii)
        ncoord(mcoord)=northnode_g(ii)
        if(inext.eq.ibeg) go to 500
        go to 401

500     continue
        write(20,490) mcoord+1,0
490     format(2i7)
        do i=1,mcoord
          write(20,510) ecoord(i),ncoord(i)
510       format(1x,f12.3,1x,f12.3)
        end do
        write(20,510) ecoord(1),ncoord(1)
        close(unit=20)
        write(6,140) trim(outfile)

600     continue
        go to 9900

9890    call write_message(leadspace='yes',endspace='yes')
9900    call close_files

        deallocate(nodenum_g,eastnode_g,northnode_g,stat=ierr)
        deallocate(elemnum_g,elemnode_g,eastelem_g,northelem_g,elemindex_g,stat=ierr)

end program rsm2srf



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



