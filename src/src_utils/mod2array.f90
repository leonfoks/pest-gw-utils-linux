!     Last change:  JD    9 May 2003    8:56 am
program mod2array

! -- Program MOD2ARRAY extracts arrays from a MODFLOW/MT3D input file and writes them
!    to separate files.

	use defn
	use inter

	implicit none

        integer, parameter :: MAXZONES=1000
        integer, parameter :: MAXLAY=1000
        integer, parameter :: MAXUNIT=1000
        integer          :: ifail,ierr,nbb,jarray,n,ilay,nd,itemp
        integer          :: modunit,outunit,arrayunit
        integer          :: nameunit,iunit,nunit,newnameunit,mlocat
        integer          :: locat,iconst,irow,icol,nrow,ncol
        integer          :: nblock,iblock,i1,i2,j1,j2,nzones,izones
        integer          :: iline,iheader,idate
        integer          :: lay(0:MAXLAY)
        integer          :: munit(MAXUNIT)
        real             :: rconst,rtemp
        real             :: rzones(MAXZONES)

        character        :: aft,air,aod,acr
        character*4      :: aext
        character*10     :: anum,alay
        character*15     :: aline
        character*20     :: fmt,atemp
        character*42     :: textid,atextid
        character*200    :: modfile,amodfile,basename,outfile,arrayfile,afile,aoutfile
        character*200    :: namefile,newnamefile,bfile,cfile
        character*200    :: mfile(MAXUNIT)
        character*100    :: eline,fline
        character*3000   :: dline
        type (modelgrid) :: gridspec

        integer, allocatable :: iarray(:,:)
        real, allocatable    :: rarray(:,:)


        mlocat=30                 ! Initial unit number for new MT3D files.

	write(amessage,5)
5       format(' Program MOD2ARRAY extracts arrays from a MODFLOW/MT3D input file and writes ', &
        'them to separate files.')
        call write_message(leadspace='yes',endspace='yes')

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

	call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')
        ncol=gridspec%ncol
        nrow=gridspec%nrow

! -- The following file is not opened with the normal file opening subroutine as we may want to open it
!    with an increased record length.

        write(6,*)
19      write(6,20,advance='no')
20      format(' Enter name of MODFLOW/MT3D input file: ')
        read(5,'(a)') modfile
        if(modfile.eq.' ') go to 19
        modfile=adjustl(modfile)
        if((modfile(1:2).eq.'e ').or.(modfile(1:2).eq.'E ')) then
          write(6,*)
          go to 10
        end if
        nbb=len_trim(modfile)
        call getfile(ifail,modfile,amodfile,1,nbb)
        if(ifail.ne.0) go to 19
        modfile=amodfile
        call addquote(modfile,amodfile)
        modunit=nextunit()
        open(unit=modunit,file=modfile,status='old',recl=3000,iostat=ierr)
        if(ierr.ne.0)then
          write(6,30) trim(amodfile)
30        format(' *** Cannot open file ',a,' - try again ***')
          go to 19
        end if
40      write(6,50,advance='no')
50      format(' Does this use MODFLOW or MT3D array header convention?  [f/t]: ')
        read(5,'(a)') aft
        if(aft.eq.' ') go to 40
        call casetrans(aft,'lo')
        if(aft.eq.'e')then
          close(unit=modunit)
          write(6,*)
          go to 19
        end if
        if((aft.ne.'f').and.(aft.ne.'t')) go to 40

100     write(6,110,advance='no')
110     format(' Enter text identifier for arrays to be extracted from this file: ')
        read(5,'(a)') textid
        if(textid.eq.' ') go to 100
        call casetrans(textid,'lo')
        atemp=textid
        atemp=adjustl(atemp)
        if(atemp(1:2).eq.'e ')then
          write(6,*)
          go to 40
        end if
        nbb=len_trim(textid)
        call getfile_front_space(ifail,textid,atextid,1,nbb)
        textid=atextid
120     write(6,130,advance='no')
130     format(' Are this/these integer array(s) or real array(s) [i/r]: ')
        read(5,'(a)') air
        if(air.eq.' ') go to 120
        call casetrans(air,'lo')
        if(air.eq.'e')then
          write(6,*)
          go to 100
        end if
        if((air.ne.'i').and.(air.ne.'r')) go to 120

! -- The MODFLOW/MT3D name file is read.

799     continue
        write(6,*)
800     call open_input_file(ifail, &
        ' Enter name of MODFLOW/MT3D name file: ',namefile,nameunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 120
        end if
        iline=0
        iunit=0
        do
          iline=iline+1
          read(nameunit,'(a)',end=850) cline
          if(cline.eq.' ') cycle
          cline=adjustl(cline)
          if(cline(1:1).eq.'#') cycle
          call linesplit(ifail,3)
          if(ifail.ne.0) go to 9300
          iunit=iunit+1
          if(iunit.gt.MAXUNIT)then
            write(amessage,810)
810         format(' Increase MAXUNIT and re-compile program.')
            go to 9890
          end if
          munit(iunit)=char2int(ifail,2)
          if(ifail.ne.0) go to 9300
          nbb=len(trim(cline))
          call getfile(ifail,cline,mfile(iunit),left_word(3),nbb)
          if(ifail.ne.0) go to 9400
        end do
850     rewind(unit=nameunit)
        nunit=iunit
        bfile=modfile
        call casetrans(bfile,'lo')
        do iunit=1,nunit
          afile=mfile(iunit)
          call casetrans(afile,'lo')
          if(afile.eq.bfile) go to 870
        end do
        call addquote(namefile,afile)
        write(amessage,860) trim(amodfile),trim(afile)
860     format(' MODFLOW/MT3D input file ',a,' not cited in name file ',a,'.')
        go to 9890
870     continue


        write(6,*)
150     write(6,160,advance='no')
160     format(' Enter filename base for array output files: ')
        read(5,'(a)') basename
        if(basename.eq.' ') go to 150
        basename=adjustl(basename)
        if((basename(1:2).eq.'e ').or.(basename(1:2).eq.'E '))then
          close(unit=nameunit)
          go to 799
        end if
        nbb=len_trim(basename)
        call getfile(ifail,basename,afile,1,nbb)
        basename=afile
161     continue
        acr='n'
        if(headerspec.eq.'yes')then
163       write(6,162,advance='no')
162       format(' Include NCOL/NROW header in these files? [y/n]: ')
          read(5,'(a)') acr
          if(acr.eq.' ') go to 163
          call casetrans(acr,'lo')
          if(acr.eq.'e')then
            write(6,*)
            go to 150
          end if
          if((acr.ne.'y').and.(acr.ne.'n')) go to 163
        end if

180     call open_output_file(ifail, &
	' Enter name for altered MODFLOW/MT3D input file: ',outfile,outunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          write(6,*)
	  escset=0
          if(headerspec.eq.'no')then
            go to 150
          else
            go to 161
          end if
	end if

181     write(6,182,advance='no')
182     format(' Use OPEN/CLOSE or DATA convention for array headers in this file  [o/d]: ')
        read(5,'(a)') aod
        if(aod.eq.' ') go to 181
        call casetrans(aod,'lo')
        if((aod.eq.'e').or.(aod.eq.'E')) then
          write(6,*)
          close(unit=outunit)
          go to 180
        end if
        if((aod.ne.'o').and.(aod.ne.'d')) go to 181


190     call open_output_file(ifail, &
        ' Enter name for altered MODFLOW/MT3D name file: ',newnamefile,newnameunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          write(6,*)
          escset=0
          go to 181
        end if
        afile=outfile
        call casetrans(afile,'lo')
        cfile=modfile
        call casetrans(cfile,'lo')
        do
          read(nameunit,'(a)',end=196) cline
          if(cline.eq.' ') go to 855
          call linesplit(ifail,3)
          if(ifail.ne.0) go to 855
          if(cline(left_word(1):left_word(1)).eq.'#') go to 855
          nbb=len_trim(cline)
          call getfile(ifail,cline,bfile,left_word(3),nbb)
          call casetrans(bfile,'lo')
          if(bfile.eq.cfile)then
            call addquote(outfile,aoutfile)
            write(newnameunit,845) cline(left_word(1):right_word(2)),trim(aoutfile)
845         format(a,5x,a)
            go to 846
          end if
855       write(newnameunit,'(a)') trim(cline)
846       continue
        end do
196     continue
        close(unit=nameunit)

! -- The MODFLOW input file is now read, looking for arrays.

       write(6,*)
       lay=0                    ! an array
       if(air.eq.'i')then
         allocate(iarray(ncol,nrow),stat=ierr)
       else
         allocate(rarray(ncol,nrow),stat=ierr)
       end if
       if(ierr.ne.0) go to 9200

       dline=' '
       jarray=0
       nd=len(dline)
       do
         read(modunit,'(a)',end=1000) dline
         if(dline(nd:nd).ne.' ')then
           write(amessage,193) trim(amodfile)
193        format(' Width of file ',a,' is too large. Increase length of DLINE variable '  &
           'and file opening record length and re-compile program.')
           go to 9890
         end if
         cline=dline
         call casetrans(cline,'lo')
         if(index(cline,trim(textid)).eq.0) then
           write(outunit,'(a)') trim(dline)
           cycle
         end if
         jarray=jarray+1
         call num2char(jarray,anum)
         n=index(cline,'layer')
         if(n.eq.0)then
           ilay=0
           if(lay(0).ne.0)then
             write(amessage,199) trim(textid)
199          format(' More than one array has been found without layer identification ', &
             'corresponding to the string "',a,'" this will lead to duplication of array ', &
             'filenames.')
             go to 9890
           end if
           lay(0)=lay(0)+1
           eline=cline(1:100)
         else
           eline=cline(1:100)
           cline=cline(n+5:)
           call linesplit(ifail,1)
           if(ifail.ne.0)then
             write(amessage,200) trim(anum),trim(textid),trim(amodfile)
200          format(' Cannot read layer number associated with occurence number ',a,   &
             ' of "',a,'" array identier string in file ',a,'.')
             go to 9890
           end if
           ilay=char2int(ifail,1)
           if(ifail.ne.0)then
             write(amessage,200) trim(anum),trim(textid),trim(amodfile)
             go to 9890
           end if
           if(ilay.le.0)then
             write(amessage,210) trim(anum),trim(textid),trim(modfile)
210          format(' Zero or negative layer number associated with occurence number ',a, &
             ' of "',a,'" array identifier string in file ',a,'.')
             go to 9890
           end if
           if(lay(ilay).ne.0)then
             call num2char(ilay,alay)
             write(amessage,211) trim(alay),trim(textid)
211          format(' More than one array has been found for layer ',a,' with text ',  &
             'identifier string of "',a,'": this will lead to duplication of array filenames.')
             go to 9890
           end if
           lay(ilay)=lay(ilay)+1
         end if

! -- The MODFLOW part of the array header is now read. First on the assumption that it is formatted.

         cline=eline
         if(air.eq.'i') then
           read(cline,220,iostat=ierr) locat,iconst,fmt
220        format(i10,i10,a20)
         else
           read(cline,230,iostat=ierr) locat,rconst,fmt
230        format(i10,f10.0,a20)
         end if
         if(ierr.eq.0)then
           if(locat.eq.0)then
             if(air.eq.'i')then
               iarray=iconst
             else
               rarray=rconst
             end if
           else if(locat.lt.0)then
             write(amessage,240) trim(anum),trim(textid),trim(modfile)
240          format(' Negative LOCAT value presently not allowed by MOD2ARRAY. See occurrence number ',a,  &
             ' of "',a,'" array identifier string in file ',a,'.')
             go to 9890
           else
             if((aft.eq.'t').and.(locat.lt.100))then
               write(amessage,241) trim(anum),trim(textid),trim(modfile)
241            format(' MOD2ARRAY does not allow a LOCAT of less than 100 in array header if MT3D array ',  &
               'header convention is employed, for this implies that the array is already in an '  &
               'external file. Error occurs at occurrence number ',a,  &
               ' of "',a,'" array identifier string in file ',a,'.')
               go to 9890
             end if
             if((aft.eq.'f').or.((aft.eq.'t').and.(locat.le.100)))then
               if(air.eq.'i')then
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (iarray(icol,irow),icol=1,ncol)
                 end do
               else
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (rarray(icol,irow),icol=1,ncol)
                 end do
               end if
             else if(locat.eq.101)then
               if(air.eq.'i')then
                 iarray=0
               else
                 rarray=0.0
               end if
               read(modunit,*,err=9000)nblock
               do iblock=1,nblock
                 if(air.eq.'i')then
                   read(modunit,*,err=9000) i1,i2,j1,j2,itemp
                   do irow=i1,i2
                     do icol=j1,j2
                       iarray(icol,irow)=itemp
                     end do
                   end do
                 else
                   read(modunit,*,err=9000) i1,i2,j1,j2,rtemp
                   do irow=i1,i2
                     do icol=j1,j2
                       rarray(icol,irow)=rtemp
                     end do
                   end do
                 end if
               end do
             else if(locat.eq.102)then
               read(modunit,*,err=9000) nzones
               if(nzones.gt.MAXZONES)then
                 write(amessage,250)
250              format(' Increase MAXZONES and re-compile program.')
                 go to 9890
               end if
               read(modunit,*,err=9000) (rzones(izones),izones=1,nzones)
               if(air.eq.'i')then
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (iarray(icol,irow),icol=1,ncol)
                 end do
                 do irow=1,nrow
                   do icol=1,ncol
                     itemp=iarray(icol,irow)
                     if((itemp.lt.0).or.(itemp.gt.nzones)) go to 9000
                     if(itemp.gt.0) iarray(icol,irow)=nint(rzones(itemp))
                   end do
                 end do
               else
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (rarray(icol,irow),icol=1,ncol)
                 end do
                 do irow=1,nrow
                   do icol=1,ncol
                     itemp=nint(rarray(icol,irow))
                     if((itemp.lt.0).or.(itemp.gt.nzones)) go to 9000
                     if(itemp.eq.0)then
                       rarray(icol,irow)=0.0
                     else
                       rarray(icol,irow)=rzones(itemp)
                     end if
                   end do
                 end do
               end if
             else if (locat.eq.103)then
               if(air.eq.'i')then
                 do irow=1,nrow
                   read(modunit,*,err=9000) (iarray(icol,irow),icol=1,ncol)
                 end do
               else
                 do irow=1,nrow
                   read(modunit,*,err=9000) (rarray(icol,irow),icol=1,ncol)
                 end do
               end if
             end if
           end if
           if(locat.ne.0)then
             if(air.eq.'i')then
               if((iconst.ne.0).and.(iconst.ne.1))then
                 iarray=iarray*iconst
               end if
             else
               if((rconst.ne.0.0).and.(rconst.ne.1.0))then
                 rarray=rarray*rconst
               end if
             end if
           end if
           fline=dline(51:)
         else
           call linesplit(ifail,1)
           if(ifail.ne.0) go to 9000
           atemp=cline(left_word(1):right_word(1))
           call casetrans(atemp,'lo')
           if(atemp.eq.'constant')then
             call linesplit(ifail,2)
             if(ifail.ne.0) go to 9000
             if(air.eq.'i')then
               itemp=char2int(ifail,2)
               if(ifail.ne.0) go to 9000
               iarray=itemp
             else
               rtemp=char2real(ifail,2)
               if(ifail.ne.0) go to 9000
               rarray=rtemp
             end if
             fline=dline(right_word(2)+1:)
           else if(atemp.eq.'internal')then
             call linesplit(ifail,4)
             if(ifail.ne.0) go to 9000
             if(air.eq.'i')then
               iconst=char2int(ifail,2)
               if(ifail.ne.0) go to 9000
             else
               rconst=char2real(ifail,2)
               if(ifail.ne.0) go to 9000
             end if
             fmt=cline(left_word(3):right_word(3))
             if(air.eq.'i')then
               if(index(fmt,'free').ne.0)then
                 do irow=1,nrow
                   read(modunit,*,err=9000) (iarray(icol,irow),icol=1,ncol)
                 end do
               else
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (iarray(icol,irow),icol=1,ncol)
                 end do
               end if
               if((iconst.ne.0).and.(iconst.ne.1))then
                 iarray=iarray*iconst
               end if
             else
               if(index(fmt,'free').ne.0)then
                 do irow=1,nrow
                   read(modunit,*,err=9000) (rarray(icol,irow),icol=1,ncol)
                 end do
               else
                 do irow=1,nrow
                   read(modunit,fmt,err=9000) (rarray(icol,irow),icol=1,ncol)
                 end do
               end if
               if((rconst.ne.0.0).and.(rconst.ne.1.0))then
                 rarray=rarray*rconst
               end if
             end if
             fline=dline(right_word(4)+1:)
           else if(atemp.eq.'external')then
             write(amessage,280) trim(anum),trim(textid),trim(modfile)
280          format(' "EXTERNAL" keyword for extractable array presently not allowed by MOD2ARRAY. ',  &
             'See occurence number ',a,  &
             ' of "',a,'" array identifer string in file ',a,'.')
             go to 9890
           else if(atemp.eq.'open/close')then
             write(amessage,290) trim(anum),trim(textid),trim(modfile)
290          format(' "OPEN/CLOSE" keyword for extractable array presently not allowed by MOD2ARRAY. ',  &
             'See occurence number ',a, &
             ' of "',a,'" array identifer string in file ',a,'.')
             go to 9890
           else
             go to 9000
           end if
         end if
         if(air.eq.'i')then
           aext='.inf'
         else
           aext='.ref'
         end if
         if(ilay.eq.0)then
           arrayfile=trim(basename)//trim(aext)
         else
           call num2char(ilay,alay)
           arrayfile=trim(basename)//trim(alay)//trim(aext)
         end if
         arrayunit=nextunit()
         call addquote(arrayfile,afile)
         open(unit=arrayunit,file=arrayfile,action='write',iostat=ierr)
         if(ierr.ne.0)then
           write(amessage,310) trim(afile)
310        format(' Cannot open file ',a,' for output.')
           go to 9890
         end if
         if(acr.eq.'y')then
           write(arrayunit,311) ncol,nrow
311        format(2i10)
         end if
         if(air.eq.'i')then
           do irow=1,nrow
             write(arrayunit,320) (iarray(icol,irow),icol=1,ncol)
320          format(20i6)
           end do
         else
           do irow=1,nrow
             write(arrayunit,330) (rarray(icol,irow),icol=1,ncol)
330          format(7(1pg14.6))
           end do
         end if
         close(unit=arrayunit)
         write(6,340) trim(afile)
340      format(' - file ',a,' written ok.')
         if(aod.eq.'o')then
           if(air.eq.'r')then
             write(outunit,350) trim(arrayfile),trim(fline)
350          format(' OPEN/CLOSE ''',a,'''  1.0  ''(FREE)''  -1            ',a)
           else
             write(outunit,349) trim(arrayfile),trim(fline)
349          format(' OPEN/CLOSE ''',a,'''  1    ''(FREE)''  -1            ',a)
           end if
         else
351        mlocat=mlocat+1
           if((mlocat.ge.100).and.(mlocat.le.103)) go to 351
           do iunit=1,nunit
             if(munit(iunit).eq.mlocat) go to 351
           end do
           if(air.eq.'r')then
             write(outunit,352) mlocat,1.0,' (7f14.0)',-1,trim(fline)
352          format(i10,f10.2,a20,i10,6x,a)
           else
             write(outunit,354) mlocat,1,' (20i6)',-1,trim(fline)
354          format(i10,i10,a20,i10,6x,a)
           end if
           write(newnameunit,353) mlocat,trim(afile)
353        format('DATA   ',i4,5x,a)
         end if
       end do

1000   if(jarray.eq.0)then
         write(amessage,1010) trim(textid),trim(amodfile)
1010     format(' Text string "',a,'" not found in file ',a,'.')
         go to 9890
       end if

       write(6,*)
       call addquote(outfile,aoutfile)
       write(6,1030) trim(aoutfile)
1030   format(' - file ',a,' written ok.')
       call addquote(newnamefile,afile)
       write(6,1030) trim(afile)

       go to 9900

9000   write(amessage,9010) trim(anum),trim(textid),trim(modfile)
9010   format(' Error reading array or array header associated with occurrence number ',a,  &
       ' of "',a,'" array identifier string in file ',a,'.')
       call write_message(leadspace='yes')
       write(6,9020)
9020   format(/,' String header follows: ')
       write(6,'(a)') trim(cline)
       go to 9900
9200   write(amessage,9210)
9210   format(' Cannot allocate sufficient memory to continue execution.')
       go to 9890
9300   call num2char(iline,aline)
       call addquote(namefile,afile)
       write(amessage,9310) trim(aline),trim(afile)
9310   format(' Cannot read unit number from line ',a,' of name file ',a,'.')
       go to 9890
9400   call num2char(iline,aline)
       call addquote(namefile,afile)
       write(amessage,9410) trim(aline),trim(afile)
9410   format(' Cannot read MODFLOW/MT3D input filename from line ',a,' of name file ',a,'.')
       go to 9890

9890   call write_message(leadspace='yes')

9900   call close_files
       if(allocated(iarray)) deallocate(iarray,stat=ierr)
       if(allocated(rarray)) deallocate(rarray,stat=ierr)

end program mod2array



subroutine getfile_front_space(ifail,cline,filename,ibeg,iend)

! Subroutine getfile extracts a filename from a string.

! -- Arguments are as follows:-
!       ifail: returned as zero if filename successfully read
!       cline: a character string containing the file name
!       filename: the name of the file read from the string
!       ibeg: character position at which to begin search for filename
!       iend: on input  - character position at which to end search for filename
!             on output - character postion at which filename ends


        integer, intent(out)               :: ifail
        integer, intent(in)                :: ibeg
        integer, intent(inout)             :: iend
        character (len=*), intent(in)      :: cline
        character (len=*), intent(out)     :: filename

        integer                            :: i,j,k
        character (len=1)                  :: aa

        ifail=0
        do i=ibeg,iend
          aa=cline(i:i)
          if((aa.ne.' ').and.(aa.ne.',').and.(aa.ne.char(9)))go to 50
        end do
        ifail=1
        return

50      if((aa.eq.'"').or.(aa.eq.''''))then
          do j=i+1,iend
            if(cline(j:j).eq.aa) go to 60
          end do
          ifail=1
          return
60        iend=j
          if(i+1.gt.j-1)then
            ifail=1
            return
          else
            filename=cline(i+1:j-1)
          end if
        else
          do j=i+1,iend
            if((cline(j:j).eq.' ').or.(cline(j:j).eq.',').or.(cline(j:j).eq.char(9)))then
              k=j-1
              go to 100
            end if
          end do
          k=iend
100       filename=cline(i:k)
          if(cline(k:k).eq.'"')then
            ifail=1
            return
          else if(cline(k:k).eq.'''')then
            ifail=1
            return
          end if

          iend=k
          filename=adjustl(filename)
        end if
        return

end subroutine getfile_front_space

