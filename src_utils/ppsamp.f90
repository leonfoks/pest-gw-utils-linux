program ppsamp

! -- Program PPSAMP samples stochastic real arrays at pilot point locations
!    and generates real arrays of residuals between random and pilot-point
!    interpolated fields.


	use defn
	use inter

	implicit none

        logical          :: lexist
        integer          :: ifail,iheader,idate,i,itemp,n,ierr,iatot,ient,j
        integer          :: ibeg,iend,nnppt,natot,nent,na,info,inunit,iline,icount1,icount2,  &
                            ifile,outunit
        integer          :: npar,ipar,itranskeep,itrans,ies,jes,npar_file,jpar
        integer          :: ncol,nrow,mrow,mcol,icellno_indiv,irow,icol
        integer          :: facunit
        integer          :: ialloc_prop,ialloc_present,ialloc_rarray2
        real             :: thresh,dumval,rtemp
        real             :: pval_temp,scale_temp,offset_temp
        double precision :: rcond,sum,rtemp2
        character*1      :: facformat
        character*1      :: asamp,amodpar
        character*4      :: aext
        character*12     :: apoint
        character*10     :: aline,anumfile
        character*12     :: aprefix,aapar
        character*20     :: atemp
        character*120    :: aprompt,firstline,atempf
        character*200    :: afile,randbase1,randbase2,facfile,parfile,parbase,tempfile,infile,outfile

        type (modelgrid) :: gridspec

        integer, allocatable          :: ipt(:),icellno(:),jcellno(:)
        integer, allocatable          :: ipt_inter(:),iloc_inter(:),icellno_inter(:)
        integer, allocatable          :: iarray(:,:)
        integer, allocatable          :: ipresent(:)
        real, allocatable             :: wt(:),prop(:)
        real, allocatable             :: rmean_inter(:),wt_inter(:)
        real, allocatable             :: fac1(:),fac2(:),fac3(:),fac4(:)
        real, allocatable             :: pval_file(:),scale_file(:),offset_file(:)
        real, allocatable             :: rarray(:,:),rarray2(:,:)
        double precision, allocatable :: work(:)
        double precision, allocatable :: ltl(:,:),pval(:)
        character*12, allocatable     :: apar(:),apar_file(:)


	write(amessage,5)
5       format(' Program PPSAMP samples stochastic real arrays at pilot point locations '  &
        'and then generates real arrays of residuals between random and pilot-point-interpolated '  &
        'parameter fields.')
	call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        asamp=' '
        amodpar=' '
        thresh=1.1e35
        tempfile='t###.###'
        ialloc_prop=0
        ialloc_present=0
        ialloc_rarray2=0

! -- The settings file is read.

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
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig.')
	  call write_message
	  go to 9900
	end if

! -- The grid specification file is read.

        call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) go to 9900
        call read_spec_dim(ifail,gridspec)
        if(ifail.ne.0) go to 9900
        call read_spec_data(ifail,gridspec)
        if(ifail.ne.0) go to 9900
        call close_spec_file(gridspec,ok='yes')

        ncol=gridspec%ncol
        nrow=gridspec%nrow

! -- The filename base of existing stochastic field array files is read.

        write(6,*)
12      write(6,13,advance='no')
13      format(' Enter filename base of random field arrays: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 12
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          write(6,*)
          call free_grid_mem(gridspec)
          go to 10
        end if
        ibeg=1
        iend=len_trim(afile)
        call getfile(ifail,afile,randbase1,ibeg,iend)
        if(ifail.ne.0) go to 12

! -- The pilot points file is read.

        write(6,*)
30      call read_pilot_points_file(ifail, &
	' Enter name of pilot points file: ')
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
	  escset=0
	  write(6,*)
	  go to 12
	end if

! -- The parameter name prefix is obtained.

        write(6,*)
40      write(6,50,advance='no')
50      format(' Enter parameter prefix: ')
        read(5,'(a)') aprefix
        call casetrans(aprefix,'lo')
        aprefix=adjustl(aprefix)
        if(aprefix(1:2).eq.'e ')then
          write(6,*)
          go to 30
        end if
        if(index(aprefix,'''').ne.0) go to 40
        if(index(aprefix,'"').ne.0) go to 40

! -- The interpolation factor file is opened and the first few lines read.

        write(6,*)
60      aprompt=' Enter name of interpolation factor file: '
        call open_input_file(ifail,aprompt,facfile,facunit,form_prompt='yes', &
        fformat=facformat)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          escset=0
          write(6,*)
          go to 40
        end if
        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) atempf
        else
          read(facunit,err=9000,end=9100) atempf
        end if
        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9100) atempf
        else
          read(facunit,err=9000,end=9100) atempf
        end if
        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) mcol,mrow
        else
          read(facunit,err=9000,end=9100) mcol,mrow
        end if
        if((mcol.ne.ncol).or.(mrow.ne.nrow))then
          write(amessage,70) trim(facfile),trim(gridspec%specfile)
70        format(' Dimensions of finite difference grid as read from file ',a,  &
          ' are different from those obtained from grid specification file ',a,'.')
          go to 9890
        end if
        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9100) nnppt
        else
          read(facunit,err=9000,end=9100) nnppt
        end if
        if(nnppt.ne.num_pilot_points) go to 9050
        do i=1,nnppt
          if(facformat.eq.'f')then
            read(facunit,'(a)',err=9000,end=9100) apoint
          else
            read(facunit,err=9000,end=9100) apoint
          end if
          call casetrans(apoint,'hi')
          apoint=adjustl(apoint)
          if(apoint.ne.pilot_point_id(i)) go to 9050
        end do

! -- Some specs are now read.

        write(6,*)
71      write(6,80,advance='no')
80      format(' Employ direct or least-squares sampling of stochastic fields?  [d/l]: ')
        read(5,'(a)') asamp
        if(asamp.eq.' ') go to 71
        call casetrans(asamp,'lo')
        asamp=adjustl(asamp)
        if(asamp.eq.'e')then
          close(unit=facunit)
          write(6,*)
          go to 60
        else if((asamp.eq.'D').or.(asamp.eq.'d'))then
          asamp='d'
        else if((asamp.eq.'L').or.(asamp.eq.'l'))then
          asamp='l'
        else
          go to 71
        end if

! -- Parameter file data is now acquired.

        write(6,*)
130     write(6,140,advance='no')
140     format(' Write new parameter value files or modify existing ones?  [n/m]: ')
        read(5,'(a)') amodpar
        if(amodpar.eq.' ') go to 130
        amodpar=adjustl(amodpar)
        call casetrans(amodpar,'lo')
        if(amodpar.eq.'e')then
          write(6,*)
          go to 71
        else if((amodpar.eq.'N').or.(amodpar.eq.'n'))then
          amodpar='n'
        else if((amodpar.eq.'M').or.(amodpar.eq.'m'))then
          amodpar='m'
        else
          go to 130
        end if
149     continue
        if(amodpar.eq.'n')then
150       write(6,160,advance='no')
160       format(' Enter name of an existing parameter value file: ')
          read(5,'(a)') afile
          if(afile.eq.' ') go to 150
          afile=adjustl(afile)
          if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
            write(6,*)
            go to 130
          end if
          ibeg=1
          iend=len_trim(afile)
          call getfile(ifail,afile,parfile,ibeg,iend)
          if(ifail.ne.0) go to 150
          inquire(file=parfile,exist=lexist)
          if(.not.lexist)then
            write(6,170)
170         format(/,' *** File does not exist - try again ***',/)
            go to 150
          end if
        end if
169     continue
        if(amodpar.eq.'n')then
171       write(6,172,advance='no')
172       format(' Enter filename base for new parameter value files: ')
          read(5,'(a)') afile
          if(afile.eq.' ') go to 171
          afile=adjustl(afile)
          if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
            write(6,*)
            go to 149
          end if
          ibeg=1
          iend=len_trim(afile)
          call getfile(ifail,afile,parbase,ibeg,iend)
          if(ifail.ne.0) go to 171
        end if
179     continue
        if(amodpar.eq.'m')then
180       write(6,190,advance='no')
190       format(' Enter filename base of existing parameter value files: ')
          read(5,'(a)') afile
          if(afile.eq.' ') go to 180
          afile=adjustl(afile)
          if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
            write(6,*)
            go to 130
          end if
          ibeg=1
          iend=len_trim(afile)
          call getfile(ifail,afile,parbase,ibeg,iend)
          if(ifail.ne.0) go to 180
        end if

! -- Finally the filename base of random field residual arrays is requested.

        write(6,*)
200     write(6,210,advance='no')
210     format(' Enter filename base for random field difference arrays: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 200
        afile=adjustl(afile)
        if((afile(1:2).eq.'e ').or.(afile(1:2).eq.'E '))then
          write(6,*)
          if(amodpar.eq.'m')then
            go to 179
          else
            go to 169
          end if
        end if
        ibeg=1
        iend=len_trim(afile)
        call getfile(ifail,afile,randbase2,ibeg,iend)
        if(ifail.ne.0) go to 200
211     write(6,212,advance='no')
212     format(' Enter dummy value for inactive cells for these arrays: ')
        itemp=key_read(dumval)
        if(itemp.ne.0) go to 211
        if(escset.eq.1) then
          escset=0
          write(6,*)
          go to 200
        end if

! -- Pilot point names are all translated to lower case.

        do i=1,num_pilot_points
          call casetrans(pilot_point_id(i),'lo')
        end do

! -- Parameter names corresponding to pilot point names are formed.

        npar=num_pilot_points
        allocate(apar(npar),pval(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        do ipar=1,npar
          if(aprefix.ne.' ')then
            atemp=trim(aprefix)//trim(pilot_point_id(ipar))
          else
            atemp=trim(pilot_point_id(ipar))
          end if
          n=len_trim(atemp)
          if(n.gt.12)then
            write(amessage,220) trim(aprefix),trim(pilot_point_id(ipar))
220         format(' If the prefix "',a,'" is affixed to the pilot point name "',a,  &
            ' the resulting parameter name is greater than 12 characters in length. ', &
            'This is not allowed.')
            go to 9890
          end if
          apar(ipar)=atemp
        end do

! -- The factor file is read a first time in order to determine the active part of the grid.
!    Information is also gathered for later storage of pilot point interpolation factors.

        write(6,230) trim(facfile)
230     format(/,' - reading interpolation factors from file ',a,'....')
        allocate(iarray(ncol,nrow),stat=ierr)
        if(ierr.ne.0) go to 9200
        iarray=-999   ! an array
        allocate(ipt(num_pilot_points),wt(num_pilot_points),stat=ierr)
        if(ierr.ne.0) go to 9200
        natot=0
        nent=0
        itranskeep=-999
        do
          if(facformat.eq.'f')then
            read(facunit,*,err=9000,end=250) icellno_indiv,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          else
            read(facunit,err=9000,end=250)   icellno_indiv,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          end if
          if(itranskeep.eq.-999)then
            itranskeep=itrans
          else
            if(itrans.ne.itranskeep)then
              write(amessage,235) trim(facfile)
235           format(' Data contained in file ',a,' specifies that for some cells interpolation ', &
              'takes place on the basis of native pilot point values while for others it takes place ', &
              'on the basis of the logs of pilot point values. PPSAMP does not allow mixing of native/log ', &
              'interpolation.')
              go to 9890
            end if
          end if
          nent=nent+1
          natot=natot+na
          irow=(icellno_indiv-1)/ncol+1
          icol=icellno_indiv-((irow-1)*ncol)
          iarray(icol,irow)=itrans
        end do
250     rewind(unit=facunit)
        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9000) atempf
        else
          read(facunit,err=9000,end=9000) atempf
        end if
        if(facformat.eq.'f')then
          read(facunit,'(a)',err=9000,end=9000) atempf
        else
          read(facunit,err=9000,end=9000) atempf
        end if
        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9000) mcol,mrow
        else
          read(facunit,err=9000,end=9000) mcol,mrow
        end if
        if(facformat.eq.'f')then
          read(facunit,*,err=9000,end=9000) nnppt
        else
          read(facunit,err=9000,end=9000) nnppt
        end if
        do i=1,nnppt
          if(facformat.eq.'f')then
            read(facunit,'(a)',err=9000,end=9000) apoint
          else
            read(facunit,err=9000,end=9000) apoint
          end if
        end do

! -- If sampling of pilot point locations is direct, bilinear interpolation factors are now evaluated.

        if(asamp.eq.'d')then
          allocate(fac1(npar),fac2(npar),fac3(npar),fac4(npar),icellno(npar),jcellno(npar), &
          stat=ierr)
          if(ierr.ne.0) go to 9200
          do ipar=1,npar
            call factor(gridspec,pilot_point_east(ipar),pilot_point_north(ipar),            &
            fac1(ipar),fac2(ipar),fac3(ipar),fac4(ipar),icellno(ipar),jcellno(ipar))
            if(icellno(ipar).eq.-999)then
              write(amessage,260) trim(pilot_point_id(ipar))
260           format(' Pilot point "',a,'" does not lie within the finite difference grid. All ',  &
              'pilot points must lie within the grid if direct sampling of stochastic real arrays ', &
              'is to take place.')
              go to 9890
            end if
            irow=(icellno(ipar)-1)/ncol+1
            icol=icellno(ipar)-((irow-1)*ncol)
            if(iarray(icol,irow).eq.-999)then
              write(amessage,270) trim(facfile),trim(pilot_point_id(ipar))
270           format(' No interpolation factors are provided in file ',a,' for the cell in which pilot ',  &
              'point "',a,'" lies. Thus it is assumed to be in an inactive cell. Hence sampling of ',  &
              'this pilot point value from stochastic fields cannot be done by the direct method.')
              go to 9890
            end if
          end do
        end if

! -- Interpolation factors are stored efficiently. We also take the opportunity to build the ltl
!    array if interpolation is least squares.

        allocate(wt_inter(natot),ipt_inter(natot),iloc_inter(nent+1),icellno_inter(nent),    &
        rmean_inter(nent),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(asamp.eq.'l')then
          allocate(ltl(npar,npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          ltl=0.0d0    ! an array
        end if
        iatot=0
        do ient=1,nent
          if(facformat.eq.'f')then
            read(facunit,*,err=9000,end=9000) icellno_indiv,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          else
            read(facunit,err=9000,end=9000)   icellno_indiv,itrans, &
            na,rtemp,((ipt(i),wt(i)),i=1,na)
          end if
          icellno_inter(ient)=icellno_indiv
          rmean_inter(ient)=rtemp
          iloc_inter(ient)=iatot+1
          do i=1,na
            iatot=iatot+1
            wt_inter(iatot)=wt(i)
            ipt_inter(iatot)=ipt(i)
            if(asamp.eq.'l')then
              ies=ipt(i)
              do j=1,i
                jes=ipt(j)
                ltl(ies,jes)=ltl(ies,jes)+wt(i)*wt(j)
                if(ies.ne.jes) ltl(jes,ies)=ltl(ies,jes)
              end do
            end if
          end do
        end do
        iloc_inter(nent+1)=iatot+1

! -- If interpolation is least squares, then the ltl matrix can now be inverted.

        if(asamp.eq.'l')then
          write(6,279)
279       format(' - inverting matrix required for least squares interpolation...')
          allocate(work(npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          call dpoco(ltl,npar,npar,rcond,work,info)
          if(info.ne.0)then
            write(amessage,280)
280         format(' Cannot undertake least squares interpolation: matrix for this ',  &
            'problem cannot be inverted.')
            go to 9890
          end if
        end if

! -- If only one parameter value file is to be read, it is read now.

        if(amodpar.eq.'n')then
          write(6,281) trim(parfile)
281       format(' - reading parameter value file ',a,'...')
          inunit=nextunit()
          open(unit=inunit,file=parfile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,282) trim(parfile)
282         format(' Cannot open parameter value file ',a,'.')
            go to 9890
          end if
          iline=1
          read(inunit,'(a)',err=9400,end=9450) cline
          call linesplit(ifail,2)
          if(ifail.ne.0) then
            write(amessage,284) trim(parfile)
284         format(' First line of parameter value file ',a,   &
            ' should read "single" or "double" followed by ',  &
            '"point" or "nopoint".')
            go to 9890
          end if
          firstline=cline
          call casetrans(firstline,'lo')
          if((firstline(left_word(1):right_word(1)).ne.'single').and.   &
             (firstline(left_word(1):right_word(1)).ne.'double'))then
             write(amessage,284) trim(parfile)
             go to 9890
          end if
          if((firstline(left_word(2):right_word(2)).ne.'point').and.    &
            (firstline(left_word(2):right_word(2)).ne.'nopoint'))then
            write(amessage,284) trim(parfile)
            go to 9890
          end if

! -- The file is read a first time to ascertain the number of parameters it holds.

          npar_file=0
          do
            continue
            iline=iline+1
            read(inunit,'(a)',err=9400,end=300) cline
            if(cline.eq.' ') cycle
            npar_file=npar_file+1
          end do
300       rewind(unit=inunit)
          read(inunit,*,err=9400,end=9400)
          if(npar_file.eq.0)then
            write(amessage,310) trim(parfile)
310         format(' No parameters are cited in parameter value file ',a,'.')
            go to 9890
          end if
          allocate(apar_file(npar_file),pval_file(npar_file),scale_file(npar_file),   &
          offset_file(npar_file),stat=ierr)
          if(ierr.ne.0) go to 9200

! -- The file is now read a second time and information contained therein is stored.

          iline=1
          do ipar=1,npar_file
320         continue
            iline=iline+1
            read(inunit,'(a)',err=9400,end=9400) cline
            if(cline.eq.' ') go to 320
            call linesplit(ifail,4)
            if(ifail.ne.0) go to 9500
            apar_file(ipar)=cline(left_word(1):right_word(1))
            call casetrans(apar_file(ipar),'lo')
            pval_file(ipar)=char2real(ifail,2)
            if(ifail.ne.0) go to 9400
            scale_file(ipar)=char2real(ifail,3)
            if(ifail.ne.0) go to 9400
            offset_file(ipar)=char2real(ifail,4)
            if(ifail.ne.0) go to 9400
          end do
          close(unit=inunit)

! -- A check is made that all pilot point parameters are represented in this file.

          do ipar=1,npar
            aapar=apar(ipar)
            do jpar=1,npar_file
              if(apar_file(jpar).eq.aapar) go to 330
            end do
            write(amessage,325) trim(aapar),trim(parfile)
325         format(' Pilot point parameter "',a,'" is not represented in parameter ',  &
            'value file ',a,'.')
            go to 9890
330         continue
          end do

        end if

! -- This is the beginning of the loop where we read stochastic real arrays.

        allocate(rarray(0:ncol+1,0:nrow+1),stat=ierr)
        if(ierr.ne.0) go to 9200
        aext='.ref'
        icount1=0
        icount2=0
        do ifile=1,100000

! - A new real array file containing a stochastic array is read in.

          call num2char(ifile,anumfile)
          infile=trim(randbase1)//trim(anumfile)//trim(aext)
          inunit=nextunit()
          open(unit=inunit,file=infile,status='old',iostat=ierr)
          if(ierr.ne.0)then
            if(icount1.eq.0)then
              if(ifile.gt.200)then
                write(amessage,321) trim(randbase1),trim(aext)
321             format(' Cannot find any real array files ',    &
                'named "',a,'*',a,'".')
                go to 9890
              else
                cycle
              end if
             else if(icount2.gt.100) then
               go to 1000
             else
               if(icount1.ne.0) icount2=icount2+1
              cycle
            end if
          end if
          write(6,*)
          write(6,331) trim(infile)
331       format(' - reading real array file ',a,'...')
          icount2=0
          icount1=icount1+1
          if(headerspec.eq.'yes') then
            read(inunit,'(a)',end=9250) cline
            call linesplit(ifail,2)
            if(ifail.ne.0) go to 9300
	    mcol=char2int(ifail,1)
            if(ifail.ne.0) go to 9300
            mrow=char2int(ifail,2)
            if(ifail.ne.0) go to 9300
	    if((mcol.ne.ncol).or.(mrow.ne.nrow)) then
              write(amessage,340) trim(infile),trim(gridspec%specfile)
340           format(' Dimensions of grid as set out in NCOL/NROW header to ', &
              'real array file ',a,' differ from those provided in grid ',  &
              'specification file ',a,'.')
              go to 9890
            end if
          end if
          rarray=thresh*2.0
          do irow=1,nrow
            read(inunit,*,err=9350,end=9250) (rarray(icol,irow),icol=1,ncol)
          end do
          close(unit=inunit)

! -- Inactive parts of the array are assigned values above the inactive threshold.
! -- Also, parts of the real array corresponding to log transformation of parameters are
!    assigned the log of the array value.

          do irow=1,nrow
            do icol=1,ncol
              if(iarray(icol,irow).eq.-999) then
                rarray(icol,irow)=2.0*thresh
              else if(iarray(icol,irow).eq.1)then
                if(rarray(icol,irow).le.0.0)then
                  write(amessage,342) trim(facfile),trim(infile)
342               format(' It is specified in the interpolation factor file ',a,  &
                  ' that interpolation of pilot point values takes place ', &
                  'on the basis of the logs of those point values. Hence array values in ', &
                  'file ',a,' cannot be zero or negative as this ', &
                  'may result in negative values being sampled for at least some pilot points.')
                  go to 9890
                end if
                rarray(icol,irow)=log10(rarray(icol,irow))
              end if
            end do
          end do

! -- Sampling of the real array now takes place.

          if(asamp.eq.'d')then

! -- First if it is by direct interpolation.

            do ipar=1,npar
              call point_interp(ncol,nrow,thresh,fac1(ipar),fac2(ipar),       &
              fac3(ipar),fac4(ipar),icellno(ipar),jcellno(ipar),rtemp,   &
              rarray)
              pval(ipar)=rtemp
            end do

          else

! -- And then if it is by least squares.
! -- First a 1-d property array is formed, indexed by cell number.

            if(ialloc_prop.eq.0)then
              allocate(prop(nent),stat=ierr)
              if(ierr.ne.0) go to 9200
              ialloc_prop=1
            end if
            ient=0
            do irow=1,nrow
              do icol=1,ncol
                if(iarray(icol,irow).ne.-999)then
                  ient=ient+1
                  prop(ient)=rarray(icol,irow)
                end if
              end do
            end do
            pval=0.0d0               ! an array
            do ient=1,nent
              do iatot=iloc_inter(ient),iloc_inter(ient+1)-1
                ies=ipt_inter(iatot)
                pval(ies)=pval(ies)+wt_inter(iatot)*prop(ient)
              end do
            end do
            call dposl(ltl,npar,npar,pval)

          end if

! -- The parameter value file is now written.
! -- Note that at this stage parameter values of log-transformed parameters are in fact the logs
!    of these values.

          if(amodpar.eq.'n')then
            outunit=nextunit()
            outfile=trim(parbase)//trim(anumfile)//'.par'
            write(6,405) trim(outfile)
405         format(' - writing parameter value file ',a,'...')
            open(unit=outunit,file=outfile,action='write',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,410) trim(outfile)
410           format(' Cannot open file ',a,' for output.')
              go to 9890
            end if
            write(outunit,'(a)') trim(firstline)
            do ipar=1,npar
              aapar=apar(ipar)
              do jpar=1,npar_file
                if(apar_file(jpar).eq.aapar) then
                  if(itranskeep.eq.0)then
                    pval_file(jpar)=pval(ipar)
                  else
                    pval_file(jpar)=10**pval(ipar)
                  end if
                  go to 420
                end if
              end do
420           continue
            end do
            do ipar=1,npar_file
              write(outunit,430) trim(apar_file(ipar)),pval_file(ipar),scale_file(ipar),offset_file(ipar)
430           format(1x,a,t20,3(1pg14.7,2x))
            end do
            close(unit=outunit)
          else

! -- First we must read the parameter value file corresponding to the random real array file.

            parfile=trim(parbase)//trim(anumfile)//'.par'
            outfile=tempfile
            write(6,439) trim(parfile)
439         format(' - re-writing param value file ',a,'...')
            inunit=nextunit()
            open(unit=inunit,file=parfile,status='old',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,440) trim(parfile)
440           format(' Cannot open parameter value file ',a,'.')
              go to 9890
            end if
            outunit=nextunit()
            open(unit=outunit,file=outfile,action='write',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,410) trim(outfile)
              go to 9890
            end if
            iline=1
            read(inunit,'(a)',err=9400,end=9450) cline
            call linesplit(ifail,2)
            if(ifail.ne.0) then
              write(amessage,284) trim(parfile)
              go to 9890
            end if
            firstline=cline
            call casetrans(firstline,'lo')
            if((firstline(left_word(1):right_word(1)).ne.'single').and.   &
               (firstline(left_word(1):right_word(1)).ne.'double'))then
               write(amessage,284) trim(parfile)
               go to 9890
            end if
            if((firstline(left_word(2):right_word(2)).ne.'point').and.    &
              (firstline(left_word(2):right_word(2)).ne.'nopoint'))then
              write(amessage,284) trim(parfile)
              go to 9890
            end if
            write(outunit,'(a)') trim(firstline)

            if(ialloc_present.eq.0)then
              allocate(ipresent(npar),stat=ierr)
              if(ierr.ne.0) go to 9200
              ialloc_present=1
            end if

            iline=1
            ipresent=0
            do
              iline=iline+1
              read(inunit,'(a)',err=9400,end=480) cline
              if(cline.eq.' ') cycle
              call linesplit(ifail,4)
              if(ifail.ne.0) go to 9500
              aapar=cline(left_word(1):right_word(1))
              call casetrans(aapar,'lo')
              pval_temp=char2real(ifail,2)
              if(ifail.ne.0) go to 9400
              scale_temp=char2real(ifail,3)
              if(ifail.ne.0) go to 9400
              offset_temp=char2real(ifail,4)
              if(ifail.ne.0) go to 9400
              do ipar=1,npar
                if(apar(ipar).eq.aapar)then
                  if(itranskeep.eq.0)then
                    write(outunit,430) trim(aapar),pval(ipar),scale_temp,offset_temp
                  else
                    write(outunit,430) trim(aapar),10**pval(ipar),scale_temp,offset_temp
                  end if
                  ipresent(ipar)=1
                  go to 470
                end if
              end do
              write(outunit,430) trim(aapar),pval_temp,scale_temp,offset_temp
470           continue
            end do
480         continue
            close(unit=inunit)
            close(unit=outunit)
            do ipar=1,npar
              if(ipresent(ipar).eq.0)then
                write(amessage,482) trim(apar(ipar)),trim(parfile)
482             format(' Pilot point parameter "',a,'" is not cited in parameter value ',  &
                'file ',a,'.')
                go to 9890
              end if
            end do
            inunit=nextunit()
            open(unit=inunit,file=tempfile,status='old')
            outunit=nextunit()
            open(unit=outunit,file=parfile)
            do
              read(inunit,'(a)',end=500)cline
              write(outunit,'(a)') trim(cline)
            end do
500         continue
            close(unit=inunit,status='delete')
            close(unit=outunit)

          end if

! -- Interpolation is now undertaken from the pilot points to a new array.

          if(ialloc_rarray2.eq.0)then
            allocate(rarray2(ncol,nrow),stat=ierr)
            if(ierr.ne.0) go to 9200
            ialloc_rarray2=1
          end if

          rarray2=dumval       ! an array
          do ient=1,nent
            irow=(icellno_inter(ient)-1)/ncol+1
            icol=icellno_inter(ient)-((irow-1)*ncol)
            sum=rmean_inter(ient)
            do i=iloc_inter(ient),iloc_inter(ient+1)-1
!              if(itranskeep.eq.0)then
                sum=sum+pval(ipt_inter(i))*wt_inter(i)
!              else
!                rtemp2=pval(ipt_inter(i))
!                if(rtemp2.le.0.0)then
!                  write(amessage,510)
!510               format(' The interpolation factor file specifies that spatial ', &
!                  'interpolation takes place on the basis of the ', &
!                  'logarithms of pilot point values. However at least one ', &
!                  'pilot point value sampled from a stochastic real array ',  &
!                  'file is negative.')
!                  go to 9890
!                end if
!                sum=sum+log10(rtemp2)*wt_inter(i)
!              end if
            end do
            if(itranskeep.eq.0)then
              rarray2(icol,irow)=rarray(icol,irow)-sum
            else
              rarray2(icol,irow)=10**(rarray(icol,irow)-sum)
            end if
          end do

! -- Finally the residuals real array is written.

          outfile=trim(randbase2)//trim(anumfile)//trim(aext)
          write(6,525) trim(outfile)
525       format(' - writing real array file ',a,'...')
          outunit=nextunit()
          open(unit=outunit,file=outfile,action='write',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,530) trim(outfile)
530         format(' Cannot open file ',a,' for output.')
            go to 9890
          end if
          if(headerspec.eq.'yes') then
            write(outunit,540) ncol,nrow
540         format(2i10)
          end if
          do irow=1,nrow
            write(outunit,550) (rarray2(icol,irow),icol=1,ncol)
550         format(8(1x,1pg14.7))
          end do
          close(unit=outunit)

        end do

1000    continue

	go to 9900

9000    write(amessage,9010) trim(facfile)
9010    format(' Error encountered in reading interpolation factor file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(facfile),trim(pilot_points_file)
9060    format(' Pilot point names cited in file ',a,' are not the same, or are not ',  &
        'arranged in the same order, as those provided in pilot points file ',a,'.')
        go to 9890
9100    write(amessage,9110) trim(facfile)
9110    format(' Unexpected end encountered to interpolation factor file ',a,'.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
	go to 9890
9250    write(amessage,9260) trim(infile)
9260    format(' Unexpected end encountered to real array file ',a,'.')
        go to 9890
9300    write(amessage,9310) trim(infile)
9310    format(' Error reading NCOL/NROW header in real array file ',a,'.')
        go to 9890
9350    write(amessage,9360) trim(infile)
9360    format(' Error encountered in reading real array from file ',a,'.')
        go to 9890
9400    call num2char(iline,aline)
        write(amessage,9410) trim(aline),trim(parfile)
9410    format(' Error encountered while reading line ',a,' of parameter ',  &
        'value file ',a,'.')
        go to 9890
9450    write(amessage,9460) trim(parfile)
9460    format(' Unexpected end encountered to parameter value file ',a,'.')
        go to 9890
9500    call num2char(iline,aline)
        write(amessage,9510) trim(aline),trim(parfile)
9510    format(' Insufficient entries found on line ',a,' of parameter value file ',a,'.')
        go to 9890


9890	call write_message(leadspace='yes',endspace='yes')
9900    call close_files

        if(allocated(ipt)) deallocate(ipt,stat=ierr)
        if(allocated(icellno)) deallocate(icellno,stat=ierr)
        if(allocated(jcellno)) deallocate(jcellno,stat=ierr)
        if(allocated(ipt_inter)) deallocate(ipt_inter,stat=ierr)
        if(allocated(iloc_inter)) deallocate(iloc_inter,stat=ierr)
        if(allocated(icellno_inter)) deallocate(icellno_inter,stat=ierr)
        if(allocated(iarray)) deallocate(iarray,stat=ierr)
        if(allocated(ipresent)) deallocate(ipresent,stat=ierr)
        if(allocated(wt)) deallocate(wt,stat=ierr)
        if(allocated(prop)) deallocate(prop,stat=ierr)
        if(allocated(rmean_inter)) deallocate(rmean_inter,stat=ierr)
        if(allocated(wt_inter)) deallocate(wt_inter,stat=ierr)
        if(allocated(fac1)) deallocate(fac1,fac2,fac3,fac4,stat=ierr)
        if(allocated(pval_file)) deallocate(pval_file,scale_file,offset_file,stat=ierr)
        if(allocated(rarray)) deallocate(rarray,stat=ierr)
        if(allocated(rarray2)) deallocate(rarray2,stat=ierr)
        if(allocated(work)) deallocate(work,stat=ierr)
        if(allocated(ltl)) deallocate(ltl,stat=ierr)
        if(allocated(pval)) deallocate(pval,stat=ierr)
        if(allocated(apar)) deallocate(apar,stat=ierr)
        if(allocated(apar_file)) deallocate(apar_file,stat=ierr)

end program ppsamp

