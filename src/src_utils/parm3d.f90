!     Last change:  J    22 Oct 2004    9:39 pm
program parm3d

! -- Program PARM3D assists in manipulation of MODFLOW/MT3D-compatible real
!    arrays in assigning hydraulic properties to a three-dimensional
!    MODFLOW/MT3D model.

	use defn
	use inter

	implicit none

        integer, parameter          :: MAXFILE=5
        integer                     :: idate,iheader,ifail,parmunit,iline,ibeg,iend,  &
                                       nlay,elevflag,ilay,lt,ierr,ncol,nrow,iz,iaction, &
                                       nfile,ifile,iwrite,icount,j,irow,icol,n1,n2,iunit, &
                                       idefault,uzero
        integer                     :: kstp,kper
        integer                     :: iaction6flag,iaction3flag,iaction7flag
        integer                     :: unassigned_int
        real                        :: unassigned_real
        real                        :: rtemp,rnum,rtop,rbot,gridelev,rdefault
        real                        :: pertim,totim
        character (len=1)           :: atype,aformat
        character (len=5)           :: aext,acol,arow
        character (len=10)          :: aline,anum
        character (len=16)          :: text
        character (len=30)          :: aatemp
        character (len=200)         :: aprompt
        character (len=200)         :: parmfile,assignfile,interpfile
        character (len=200)         :: afile,bfile
        character (len=500)         :: dline
	type (modelgrid)            :: gridspec

        integer, allocatable            :: iarray(:,:,:),intarraytop(:,:), &
                                           intarraybot(:,:),itemparray(:,:)
        real, allocatable               :: rarray(:,:,:),assignarray(:,:),  &
                                           interparray(:,:,:),relev(:), &
                                           bot(:,:,:)

        character (len=200),allocatable :: intfile(:),realfile(:),elevfile(:)

! -- Functions

        real new_value,interp_value

! -- Print header message.


	write(amessage,5)
5	format(' Program PARM3D assists in manipulation of MODFLOW/MT3D compatible ', &
        'real arrays in the parameterisation of a complex three-dimensional ', &
        'model.')
	call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        iaction6flag=0
        iaction3flag=0
        iaction7flag=0
        idefault=0
        kstp=0
        kper=0
        text=' '
        pertim=0.0
        totim=0.0

	call read_settings(ifail,idate,iheader)
	if(ifail.eq.1) then
	  write(amessage,7)
7	  format(' A settings file (settings.fig) was not found in the ', &
	  'current directory.')
          go to 9890
	else if(ifail.eq.2) then
	  write(amessage,8)
8	  format(' Error encountered while reading settings file settings.fig')
          go to 9890
	endif
	if((iheader.ne.0).or.(headerspec.eq.' ')) then
	  write(amessage,6)
6	  format(' Cannot read array header specification from settings file ', &
	  'settings.fig')
          go to 9890
	end if

! -- Initialisation.

        unassigned_int=-99999999
        unassigned_real=-3.71e35

25	aprompt=' Enter name of PARM3D control file: '
	call open_input_file(ifail,aprompt,parmfile,parmunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) go to 9900

        call addquote(parmfile,afile)

! -- The control data section of the PARM3D input file is read.

        iline=0
        call get_next_line(ifail,parmunit,iline,cline)
        if(ifail.eq.1) then
          go to 9000
        else if(ifail.eq.2) then
          go to 9050
        end if
        call casetrans(cline,'lo')
        cline=adjustl(cline)
        if(cline(1:1).ne.'*')then
          write(amessage,30) trim(afile)
30        format(' File ',a,' must begin with "* control data" section header ',  &
          'line.')
          go to 9890
        end if
        if(index(cline,'control').eq.0) then
          write(amessage,30) trim(afile)
          go to 9890
        end if

        call get_next_line(ifail,parmunit,iline,cline)
        if(ifail.eq.1) then
          go to 9000
        else if(ifail.eq.2) then
          go to 9050
        end if
        cline=adjustl(cline)
        ibeg=1
        iend=len_trim(cline)
        call getfile(ifail,cline,gridspec%specfile,ibeg,iend)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,40) trim(aline),trim(afile)
40        format(' Cannot read grid specification filename from line ',a,  &
          ' of file ',a)
          go to 9890
        end if

        call get_next_line(ifail,parmunit,iline,cline)
        if(ifail.eq.1) then
          go to 9000
        else if(ifail.eq.2) then
          go to 9050
        end if
        call linesplit(ifail,2)
        if(ifail.ne.0) go to 9100
        nlay=char2int(ifail,1)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,50) trim(aline),trim(afile)
50        format(' Cannot read number of layers from line ',a,' of file ',a,'.')
          go to 9890
        end if
        if(nlay.le.0)then
          call num2char(iline,aline)
          write(amessage,52) trim(aline),trim(afile)
52        format(' Number of model layers out of range on line ',a,' of file ',a,'.')
          go to 9890
        end if
        elevflag=char2int(ifail,2)
        if(ifail.ne.0)then
          call num2char(iline,aline)
          write(amessage,60) trim(aline),trim(afile)
60        format(' Cannot read ELEVFLAG from line ',a,' of file ',a,'.')
          go to 9890
        end if

! -- The grid specification file is now read.

        call casetrans(gridspec%specfile,'lo')      ! must not do in unix
        call addquote(gridspec%specfile,bfile)
        gridspec%specunit=nextunit()
        open(unit=gridspec%specunit,file=gridspec%specfile,status='old',iostat=ierr)
        if(ierr.ne.0)then
          write(amessage,215) trim(bfile)
215       format(' Cannot open grid specification file ',a,'.')
          go to 9890
        end if
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	ncol=gridspec%ncol
	nrow=gridspec%nrow
        if((ncol.le.0).or.(nrow.le.0))then
          write(amessage,217) trim(bfile)
217       format(' Illegal grid dimensions supplied in file ',a,'.')
          go to 9890
        end if
	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
        close(unit=gridspec%specunit)
        write(6,210) trim(bfile)
210     format(' - file ',a,' read ok.')

! -- Some memory is allocated.

        allocate(iarray(ncol,nrow,nlay),rarray(ncol,nrow,nlay),stat=ierr)
        if(ierr.ne.0) go to 9300
        iarray=unassigned_int
        rtemp=1.0
        rarray=unassigned_real*(1.0+2.0*epsilon(rtemp))
        allocate(intfile(nlay),realfile(nlay),stat=ierr)
        if(ierr.ne.0) go to 9300
        if(elevflag.ne.0)then
          allocate(elevfile(0:nlay),bot(ncol,nrow,0:nlay),stat=ierr)
          if(ierr.ne.0) go to 9300
        end if

! -- The "Layer Files" section of the PARM3D input control file is read.

        call get_next_line(ifail,parmunit,iline,cline)
        if(ifail.eq.1) then
          go to 9000
        else if(ifail.eq.2) then
          go to 9050
        end if
        call casetrans(cline,'lo')
        cline=adjustl(cline)
        if(cline(1:1).ne.'*')then
          call num2char(iline,aline)
          write(amessage,70) trim(aline),trim(afile)
70        format(' "* layer files" section header expected on line ',a' of file ',a,'.')
          go to 9890
        end if
        if(index(cline,'layer').eq.0) then
          call num2char(iline,aline)
          write(amessage,70) trim(aline),trim(afile)
          go to 9890
        end if

        do ilay=1,nlay
          call get_next_line(ifail,parmunit,iline,cline)
          if(ifail.eq.1) then
            go to 9000
          else if(ifail.eq.2) then
            go to 9050
          end if
          ibeg=1
          lt=len_trim(cline)
          iend=lt
          call getfile(ifail,cline,intfile(ilay),ibeg,iend)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,80) trim(aline),trim(afile)
80          format(' Cannot read integer array filename from line ',a,  &
            ' of file ',a,'.')
            go to 9890
          end if
          ibeg=iend+1
          iend=lt
          call getfile(ifail,cline,realfile(ilay),ibeg,iend)
          if(ifail.ne.0)then
            call num2char(iline,aline)
            write(amessage,90) trim(aline),trim(afile)
90          format(' Cannot read real array filename from line ',a,  &
            ' of file ',a,'.')
            go to 9890
          end if
          call casetrans(intfile(ilay),'lo')    !not in unix
          call casetrans(realfile(ilay),'lo')   !not in unix
        end do

! -- If necessary the "elevation files" section of the PARM3D control file is read.

        if(elevflag.ne.0)then
          call get_next_line(ifail,parmunit,iline,cline)
          if(ifail.eq.1) then
            go to 9000
          else if(ifail.eq.2) then
            go to 9050
          end if
          call casetrans(cline,'lo')
          cline=adjustl(cline)
          if(cline(1:1).ne.'*')then
            call num2char(iline,aline)
            write(amessage,110) trim(aline),trim(afile)
110         format(' "* elevation files" section header expected on line ',a,  &
            ' of file ',a,'.')
            go to 9890
          end if
          if(index(cline,'elevation').eq.0) then
            call num2char(iline,aline)
            write(amessage,110) trim(aline),trim(afile)
            go to 9890
          end if

          do ilay=0,nlay
            call get_next_line(ifail,parmunit,iline,cline)
            if(ifail.eq.1) then
              go to 9000
            else if(ifail.eq.2) then
              go to 9050
            end if
            ibeg=1
            iend=len_trim(cline)
            call getfile(ifail,cline,elevfile(ilay),ibeg,iend)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,120) trim(aline),trim(afile)
120           format(' Cannot read elevation array filename from line ',a,  &
              ' of file ',a,'.')
              go to 9890
            end if
            call casetrans(elevfile(ilay),'lo')    !not in unix
          end do
        end if

! -- Integer arrays are read as necessary.

        do ilay=1,nlay
          bfile=intfile(ilay)
          call casetrans(bfile,'lo')
          bfile=adjustl(bfile)
          if(bfile(1:5).eq.'none ') cycle
          call get_array(ifail,1,ncol,nrow,iarray(1,1,ilay),rarray(1,1,1),  &
          intfile(ilay),1)
          if(ifail.ne.0) go to 9890
        end do

! -- Layer elevation arrays are read if necessary.

        if(elevflag.ne.0)then
          do ilay=0,nlay
            bfile=elevfile(ilay)
            call casetrans(bfile,'lo')
            bfile=adjustl(bfile)
            if(bfile(1:5).eq.'none ') then
              write(amessage,270)
270           format(' The name of an elevation array file must not be supplied as "none" ', &
              'in "elevation files" section of PARM3D input file.')
              go to 9890
            end if
            call get_array(ifail,2,ncol,nrow,iarray(1,1,1),bot(1,1,ilay),     &
            elevfile(ilay),1)
            if(ifail.ne.0) go to 9890
          end do
        end if

! -- All input data has been read. The "parameter value assignment" is now read and processed.

        call get_next_line(ifail,parmunit,iline,cline)
        if(ifail.eq.1) then
          go to 9000
        else if(ifail.eq.2) then
          go to 9050
        end if
        call casetrans(cline,'lo')
        cline=adjustl(cline)
        if(cline(1:1).ne.'*')then
          call num2char(iline,aline)
          write(amessage,310) trim(aline),trim(afile)
310       format(' "* parameter value assignment" section header expected on line ',  &
          a,' of file ',a,'.')
          go to 9890
        end if
        if(elevflag.eq.0)then
          if(index(cline,'elevation').ne.0)then
            call num2char(iline,aline)
            write(amessage,311) trim(aline),trim(afile)
311         format(' An "elevation files" section header was found at line ',a,  &
            ' of file ',a,'; this section is not allowed unless ELEVFLAG ',    &
            'is set to 1 in "control data" section of PARM3D input file.')
            go to 9890
          end if
        end if
        if(index(cline,'parameter').eq.0) then
          call num2char(iline,aline)
          write(amessage,310) trim(aline),trim(afile)
          go to 9890
        end if

        do
          call get_next_line(ifail,parmunit,iline,cline)
          if(ifail.eq.1) then
            go to 9000
          else if(ifail.eq.2) then
            go to 1000
          end if
          lt=len_trim(cline)
          call num2char(iline,aline)
          call linesplit(ifail,2)
          if(ifail.ne.0)then
            write(amessage,320) trim(aline),trim(afile)
320         format(' Insufficient items on line ',a,' of file ',a,'.')
            go to 9890
          end if
          aatemp=cline(left_word(1):right_word(1))
          call casetrans(aatemp,'lo')
          if(aatemp.eq.'layer')then
            atype='l'
          else if(aatemp.eq.'zone')then
            atype='z'
          else if(aatemp.eq.'default')then
            if(idefault.eq.1)then
              write(amessage,325) trim(afile)
325           format(' More than one DEFAULT value is provided in ',  &
              '"parameter value assignment" section of file ',a,'.')
              go to 9890
            end if
            rdefault=char2real(ifail,2)
            if(ifail.ne.0)then
              write(amessage,326) trim(aline),trim(afile)
326           format(' Cannot read default value from line ',a,' of file ',a,'.')
              go to 9890
            end if
            idefault=1
            go to 552
          else
            write(amessage,330) trim(aline),trim(afile)
330         format(' First entry on line ',a,' of file ',a,' must be "zone", ', &
            '"layer" or "default".')
            go to 9890
          end if
          iz=char2int(ifail,2)
          if(ifail.ne.0)then
            write(amessage,340) trim(aline),trim(afile)
340         format(' Second entry on line ',a,' of file ',a,' should be layer or zone ', &
            'number; this must be an integer.')
            go to 9890
          end if
          if(atype.eq.'l')then
            if((iz.le.0).or.(iz.gt.nlay))then
              write(amessage,350) trim(aline),trim(afile)
350           format(' Layer number out of range at line ',a,' of file ',a,'.')
              go to 9890
            end if
          end if
          ibeg=right_word(2)+1
          iend=lt
          call getfile(ifail,cline,assignfile,ibeg,iend)
          call casetrans(assignfile,'lo')     ! not in unix
          aatemp=assignfile(1:15)
          call casetrans(aatemp,'lo')
          if(aatemp.eq.'interp_arith')then
            iaction=1
          else if(aatemp.eq.'interp_geom')then
            iaction=2
          else if(aatemp.eq.'grad_arith')then
            iaction=3
          else if(aatemp.eq.'grad_geom')then
            iaction=4
          else if(aatemp.eq.'rezone')then
            if(atype.eq.'z')then
              write(amessage,352) trim(aline),trim(afile)
352           format(' First entry on line ',a,' of file ',a,' must be LAYER if ', &
              'REZONE instruction is used.')
              go to 9890
            end if
            iaction=7
          else
            call char2num(ifail,aatemp,rnum)
            if(ifail.eq.0)then
              iaction=5
            else
              if(iaction6flag.eq.0)then
                allocate(assignarray(ncol,nrow),stat=ierr)
                if(ierr.ne.0) go to 9300
                iaction6flag=1
              end if
              dline=cline
              call get_array(ifail,2,ncol,nrow,iarray(1,1,1),assignarray,assignfile,1)
              if(ifail.ne.0) go to 9890
              iaction=6
              cline=dline
            end if
          end if
          if((iaction.eq.3).or.(iaction.eq.4))then
            if(iaction3flag.eq.0)then
              allocate(interparray(ncol,nrow,MAXFILE),stat=ierr)
              if(ierr.ne.0) go to 9300
              allocate(intarraytop(ncol,nrow),intarraybot(ncol,nrow),stat=ierr)
              if(ierr.ne.0) go to 9300
              allocate(relev(MAXFILE),stat=ierr)
              if(ierr.ne.0) go to 9300
              iaction3flag=1
            end if
            ibeg=iend+1
            iend=lt
            if(iend.lt.ibeg)then
              write(amessage,370) trim(aline),trim(afile)
              go to 9890
            end if
            call getfile(ifail,cline,aatemp,ibeg,iend)
            if(ifail.ne.0)then
              write(amessage,370) trim(aline),trim(afile)
              go to 9890
            end if
            call char2num(ifail,aatemp,nfile)
            if(ifail.ne.0)then
              write(amessage,370) trim(aline),trim(afile)
370           format(' Integer expected as fourth entry on line ',a,' of file ',a,'.')
              go to 9890
            end if
            if(nfile.le.1)then
              write(amessage,380) trim(aline),trim(afile)
380           format(' Number greater than 1 expected as fourth entry on line ',a,   &
              ' of file ',a,'.')
              go to 9890
            end if
            if(nfile.gt.MAXFILE)then
              write(amessage,390) trim(aline),trim(afile)
390           format(' Fourth entry on line ',a,' of file ',a,' too large - increase ', &
              'MAXFILE and re-compile program.')
              go to 9890
            end if
            dline=cline
            do ifile=1,nfile
              ibeg=iend+1
              iend=lt
              if(iend.lt.ibeg)then
                write(amessage,400) trim(aline),trim(afile)
                go to 9890
              end if
              call getfile(ifail,dline,interpfile,ibeg,iend)
              if(ifail.ne.0)then
                write(amessage,400) trim(aline),trim(afile)
400             format(' Error reading filename from line ',a,' of file ',a,'.')
                go to 9890
              end if
              call get_array(ifail,2,ncol,nrow,iarray(1,1,1),interparray(1,1,ifile),     &
              interpfile,1)
              if(ifail.ne.0) go to 9890
            end do
            cline=dline
          else if(iaction.eq.7)then
            ibeg=iend+1
            iend=lt
            if(iend.lt.ibeg)then
              write(amessage,401) trim(aline),trim(afile)
401           format(' Integer array filename expected as fourth entry on line ',a,  &
              ' of file ',a,'.')
              go to 9890
            end if
            call getfile(ifail,cline,intfile(iz),ibeg,iend)
            if(ifail.ne.0)then
              write(amessage,402) trim(aline),trim(afile)
402           format(' Cannot read integer array filename from line ',a,' of file ',a,'.')
              go to 9890
            end if
            ibeg=iend+1
            iend=lt
            if(iend.lt.ibeg)then
              write(amessage,403) trim(aline),trim(afile)
403           format(' Last entry on line ',a,' of file ',a,' should be "use_zero" ', &
              'or "ignore_zero".')
              go to 9890
            end if
            call getfile(ifail,cline,aatemp,ibeg,iend)
            call casetrans(aatemp,'lo')
            if(aatemp.eq.'use_zero')then
              uzero=1
            else if(aatemp.eq.'ignore_zero')then
              uzero=0
            else
              write(amessage,403) trim(aline),trim(afile)
              go to 9890
            end if
          end if
          if(iaction.ne.7)then
            ibeg=iend+1
            iend=lt
            if(iend.lt.ibeg)then
              write(amessage,410) trim(aline),trim(afile)
              go to 9890
            end if
            call getfile(ifail,cline,aatemp,ibeg,iend)
            call casetrans(aatemp,'lo')
            if(aatemp.eq.'overwrite')then
              iwrite=1
            else if(aatemp.eq.'arithav')then
              iwrite=2
            else if(aatemp.eq.'geomav')then
              iwrite=3
            else if(aatemp.eq.'max')then
              iwrite=4
            else if(aatemp.eq.'min')then
              iwrite=5
            else
              write(amessage,410) trim(aline),trim(afile)
410           format(' Last entry on line ',a,' of file ',a,' should be "overwrite", ', &
              '"arithav", "geomav", "max" or "min".')
              go to 9890
            end if
          end if

! -- Direct assignment (iaction = 5 or 6).

          if((iaction.eq.5).or.(iaction.eq.6))then
            if(atype.eq.'l')then
              if(iwrite.eq.1)then
                if(iaction.eq.5)then
                  do irow=1,nrow
                    do icol=1,ncol
                      rarray(icol,irow,iz)=rnum
                    end do
                  end do
                else
                  do irow=1,nrow
                    do icol=1,ncol
                      rarray(icol,irow,iz)=assignarray(icol,irow)
                    end do
                  end do
                end if
              else
                do irow=1,nrow
                  do icol=1,ncol
                    rtemp=rarray(icol,irow,iz)
                    if(rtemp.lt.unassigned_real)then
                      if(iaction.eq.5)then
                        rarray(icol,irow,iz)=rnum
                      else
                        rarray(icol,irow,iz)=assignarray(icol,irow)
                      end if
                    else
                      if(iaction.eq.6) rnum=assignarray(icol,irow)
                      rarray(icol,irow,iz)=new_value(ifail,iwrite,rnum,rtemp,aline,afile)
                      if(ifail.ne.0) go to 9890
                    end if
                  end do
                end do
              end if
            else if(atype.eq.'z')then
              if(iwrite.eq.1)then
                do ilay=1,nlay
                  do irow=1,nrow
                    do icol=1,ncol
                      if(iarray(icol,irow,ilay).eq.iz) then
                        if(iaction.eq.5)then
                          rarray(icol,irow,ilay)=rnum
                        else
                          rarray(icol,irow,ilay)=assignarray(icol,irow)
                        end if
                      end if
                    end do
                  end do
                end do
              else
                do ilay=1,nlay
                  do irow=1,nrow
                    do icol=1,ncol
                      if(iarray(icol,irow,ilay).eq.iz)then
                        rtemp=rarray(icol,irow,ilay)
                        if(iaction.eq.6)rnum=assignarray(icol,irow)
                        if(rtemp.lt.unassigned_real)then
                          rarray(icol,irow,ilay)=rnum
                        else
                          rarray(icol,irow,ilay)=new_value(ifail,iwrite,rnum,rtemp,  &
                          aline,afile)
                          if(ifail.ne.0) go to 9890
                        end if
                      end if
                    end do
                  end do
                end do
              end if
            end if

! -- Re-zone (iaction=7).

          else if(iaction.eq.7)then
            if(iaction7flag.eq.0)then
              allocate(itemparray(ncol,nrow),stat=ierr)
              if(ierr.ne.0) go to 9300
              iaction7flag=1
            end if
            dline=cline
            call get_array(ifail,1,ncol,nrow,itemparray,rarray(1,1,1),intfile(iz),1)
            if(ifail.ne.0) go to 9890
            cline=dline
            if(uzero.eq.1)then
              do irow=1,nrow
                do icol=1,ncol
                  iarray(icol,irow,iz)=itemparray(icol,irow)
                end do
              end do
            else
              do irow=1,nrow
                do icol=1,ncol
                  if(itemparray(icol,irow).ne.0) iarray(icol,irow,iz)=itemparray(icol,irow)
                end do
              end do
            end if

! -- Vertical interpolation between already-assigned values.
!    Note: if an array element that is designated for interpolation has already been
!    assigned a value it is ignored from the point of view of interpolation; that is
!    interpolation still takes place from existing overlying or underlying elements.
!    However the "overwrite", "max" "mine", "arithav", "geomav" stuff still holds for
!    assigning it a value.

          else if((iaction.eq.1).or.(iaction.eq.2))then
            if(elevflag.eq.0)then
              write(amessage,420) trim(aline),trim(afile)
420           format(' An INTERP_ARITH or INTERP_GEOM function is not permitted ', &
              'at line ',a,  &
              ' of file ',a,' as ELEVFLAG has been set to zero and no layer ', &
              'bottom elevations have thus been read.')
              go to 9890
            end if
            if(atype.eq.'l')then
              do irow=1,nrow
                do icol=1,ncol
                  rnum=interp_value(ifail,iaction,ncol,nrow,nlay,icol,irow,iz,rarray,bot, &
                  unassigned_real,aline,afile)
                  if(ifail.ne.0) go to 9890
                  if(iwrite.eq.1)then
                    rarray(icol,irow,iz)=rnum
                  else
                    rtemp=rarray(icol,irow,iz)
                    if(rtemp.lt.unassigned_real)then
                      rarray(icol,irow,iz)=rnum
                    else
                      rarray(icol,irow,iz)=new_value(ifail,iwrite,rnum,rtemp,aline,afile)
                      if(ifail.ne.0) go to 9890
                    end if
                  end if
                end do
              end do
            else if(atype.eq.'z')then
              do irow=1,nrow
                do icol=1,ncol
                  do ilay=1,nlay
                    if(iarray(icol,irow,ilay).eq.iz)then
                      rnum=interp_value(ifail,iaction,ncol,nrow,nlay,icol,irow,ilay,rarray,bot, &
                      unassigned_real,aline,afile)
                      if(ifail.ne.0) go to 9890
                      if(iwrite.eq.1)then
                        rarray(icol,irow,ilay)=rnum
                      else
                        rtemp=rarray(icol,irow,ilay)
                        if(rtemp.lt.unassigned_real)then
                          rarray(icol,irow,ilay)=rnum
                        else
                          rarray(icol,irow,ilay)=new_value(ifail,iwrite,rnum,rtemp,  &
                          aline,afile)
                          if(ifail.ne.0) go to 9890
                        end if
                      end if
                    end if
                  end do
                end do
              end do
            end if

! -- The gradational (ie. GRAD_ARITH or GRAD_GEOM) functions are now implemented.

          else if((iaction.eq.3).or.(iaction.eq.4))then
            if(elevflag.eq.0)then
              write(amessage,425) trim(aline),trim(afile)
425           format(' A GRAD_GEOM or GRAD_ARITH function is not permitted at line ',a,  &
              ' of file ',a,' as ELEVFLAG is set to zero and no layer ', &
              'bottom elevations have thus been read.')
              go to 9890
            end if
            if(atype.eq.'l')then
              write(amessage,450) trim(aline),trim(afile)
450           format(' LAYER option is not permitted with GRAD_GEOM or ',  &
                     'GRAD_ARITH function; error at line ',a,' of file ',a,'.')
              go to 9890
            end if
            intarraytop=0          ! a 2d array
            icount=0
            do irow=1,nrow
              do icol=1,ncol
                do ilay=1,nlay
                  if(iarray(icol,irow,ilay).eq.iz) then
                    intarraytop(icol,irow)=ilay
                    icount=icount+1
                    go to 460
                  end if
                end do
460             continue
                if(ilay.eq.nlay)then
                  intarraybot(icol,irow)=nlay
                else
                  intarraybot(icol,irow)=ilay
                  do j=ilay+1,nlay
                    if(iarray(icol,irow,j).eq.iz) intarraybot(icol,irow)=j
                  end do
                end if
              end do
            end do
            if(icount.eq.0) go to 551  !xxxx
            do irow=1,nrow
              do icol=1,ncol
                if(intarraytop(icol,irow).ne.0)then
                  relev(1)=bot(icol,irow,intarraytop(icol,irow)-1)
                  relev(nfile)=bot(icol,irow,intarraybot(icol,irow))
                  if(nfile.gt.2)then
                    rtop=relev(1)
                    rbot=relev(nfile)
                    do ifile=2,nfile-1
                      relev(ifile)=rtop-float(ifile-1)/float(nfile-1)*(rtop-rbot)
                    end do
                  end if
                  do ilay=1,nlay
                    if(iarray(icol,irow,ilay).eq.iz)then
                      gridelev=0.5*(bot(icol,irow,ilay-1)+bot(icol,irow,ilay))
                      do ifile=1,nfile
                        if(gridelev.ge.relev(ifile))then
                          if(ifile.eq.1)then
                            call num2char(icol,acol)
                            call num2char(irow,arow)
                            write(amessage,519) trim(arow),trim(acol)
519                         format(' There appears to be a negative or zero thickness somewhere ', &
                            'in (ROW,COL) (',a,',',a,') of the finite difference grid.')
                            go to 9890
                          end if
                          if(iaction.eq.3)then
                            rnum=interparray(icol,irow,ifile)  &
                            +(gridelev-relev(ifile))/   &
                            (relev(ifile-1)-relev(ifile))*   &
                            (interparray(icol,irow,ifile-1)-interparray(icol,irow,ifile))
!                            write(6,*) ifile,nfile,gridelev,relev(ifile-1),relev(ifile) !debug
                          else
                            rtop=interparray(icol,irow,ifile-1)
                            rbot=interparray(icol,irow,ifile)
                            if((rtop.le.0.0).or.(rbot.le.0.0))then
                              call num2char(icol,acol)
                              call num2char(irow,arow)
                              write(amessage,520) trim(arow),trim(acol),trim(aline),trim(afile)
520                           format(' Array value for (ROW,COL) (',a,',',a,') of at least ', &
                              'one of the arrays cited on line ',a,' of file ',a,' has a ',  &
                              'negative value in zonal area; this is not permitted for the ', &
                              'GEOMAV gradation option.')
                              go to 9890
                            end if
                            rtop=log(rtop)
                            rbot=log(rbot)
                            rnum=rbot+(gridelev-relev(ifile))/  &
                                 (relev(ifile-1)-relev(ifile))*(rtop-rbot)
                            rnum=exp(rnum)
                          end if
                          if(iwrite.eq.1)then
                            rarray(icol,irow,ilay)=rnum
                          else
                            rtemp=rarray(icol,irow,ilay)
                            if(rtemp.lt.unassigned_real)then
                              rarray(icol,irow,ilay)=rnum
                            else
                              rarray(icol,irow,ilay)=new_value(ifail,iwrite,rnum,rtemp,  &
                              aline,afile)
                              if(ifail.ne.0) go to 9890
                            end if
                          end if
                          go to 550
                        end if
                      end do
550                   continue
                    end if
                  end do
                end if
              end do
            end do
551         continue
          end if
552       continue
        end do

1000    continue
        close(unit=parmunit)
        if(idefault.eq.0)then
          write(amessage,1010) trim(afile)
1010      format(' No default property value has been provided in ', &
          '"parameter value assignment" section of file ',a,'.')
          go to 9890
        end if
        do irow=1,nrow
          do icol=1,ncol
            do ilay=1,nlay
              if(rarray(icol,irow,ilay).lt.unassigned_real)  &
              rarray(icol,irow,ilay)=rdefault
            end do
          end do
        end do
        write(6,210) trim(afile)
        write(6,*)

! -- All instructions have now been processed. The array files are now written.

        do ilay=1,nlay
          bfile=realfile(ilay)
          call casetrans(bfile,'lo')
          bfile=adjustl(bfile)
          if(bfile(1:5).eq.'none ') cycle
          call addquote(realfile(ilay),bfile)
          call casetrans(bfile,'lo')            ! Not for UNIX
          n2=len_trim(realfile(ilay))
          n1=n2-3
          if(n1.lt.1)n1=1
          aext=realfile(ilay)(n1:n2)
          call casetrans(aext,'lo')
          if(aext.eq.'.reu') then
            aformat='u'
          else
            aformat='f'
          end if
          iunit=nextunit()
          if(aformat.eq.'f')then
            open(unit=iunit,file=realfile(ilay),action='write',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,570) trim(bfile)
570           format(' Cannot open real array file ',a,' for writing.')
              go to 9890
            end if
          else
            open(unit=iunit,file=realfile(ilay),form='binary',  &
            action='write',iostat=ierr)
            if(ierr.ne.0)then
              write(amessage,580) trim(bfile)
580           format(' Cannot open unformatted real array file ',a,' for writing.')
              go to 9890
            end if
          end if
          if(aformat.eq.'f')then
            if(headerspec.eq.'yes')then
              write(iunit,*) ncol,nrow
            end if
            do irow=1,nrow
              write(iunit,585) (rarray(icol,irow,ilay),icol=1,ncol)
585           format(7(1pe14.6))
!585           format(7(1x,1pg14.7))
            end do
          else
!            write(iunit) kstp,kper,pertim,totim,text,ncol,nrow,ilay
            write(iunit)
            write(iunit) ((rarray(icol,irow,ilay),icol=1,ncol),irow=1,nrow)
          end if
          write(6,590) trim(bfile)
590       format(' - file ',a,' written ok.')
        end do
        go to 9900


9000    call num2char(iline,aline)
        write(amessage,9010) trim(aline),trim(afile)
9010    format(' Error reading line ',a,' of file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(afile)
9060    format(' Unexpected end encountered to file ',a,'.')
        go to 9890
9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(afile)
9110    format(' Insufficient entries on line ',a,' of file ',a,'.')
        go to 9890
9300    write(amessage,9310)
9310    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890





9890    call write_message(leadspace='yes',endspace='yes')
9900    call close_files

        deallocate(iarray,intarraytop,intarraybot,itemparray,stat=ierr)
        deallocate(rarray,assignarray,interparray,relev,bot,stat=ierr)
        deallocate(intfile,realfile,elevfile,stat=ierr)

end program parm3d



subroutine  get_next_line(ifail,iunit,iline,cline)

! -- Subroutine GET_NEXT_LINE retreives the next line from a file.

        implicit none
        integer, intent(out)               :: ifail
        integer, intent(in)                :: iunit
        integer, intent(inout)             :: iline
        character (len=*), intent(out)     :: cline

        ifail=0
        do
          iline=iline+1
          read(iunit,'(a)',err=9000,end=9050) cline
          if(cline.eq.' ') cycle
          if(cline(1:1).eq.'#') cycle
          return
        end do

9000    ifail=1
        return
9050    ifail=2
        return

end subroutine get_next_line



real function new_value(ifail,iwrite,rnum,rexist,aline,afile)

! -- Function NEW_VALUE calculates the new value for an array element.

        use defn
        implicit none

        integer, intent(out)   :: ifail
        integer, intent(in)    :: iwrite
        real, intent(in)       :: rnum,rexist
        character (len=*)      :: aline,afile

        ifail=0
        if(iwrite.eq.2)then
          new_value=0.5*(rnum+rexist)
        else if(iwrite.eq.3)then
          if((rexist.le.0.0).or.(rnum.le.0.0))then
            write(amessage,20) trim(aline),trim(afile)
20          format(' Cannot carry out GEOMAV operation on line ',a,  &
            ' of file ',a,' as at least one existing or new array element is zero or ', &
            'negative.')
            go to 9890
          end if
          new_value=sqrt(rnum*rexist)
        else if(iwrite.eq.4)then
          new_value=max(rnum,rexist)
        else if(iwrite.eq.5)then
          new_value=min(rnum,rexist)
        end if
        return

9890    ifail=1
        return

end function new_value




subroutine get_array(ifail,itype,ncol,nrow,iarray,rarray,filename,ireport)

! -- Subroutine GET_ARRAY opens an array file and reads the array contained therein.

        use defn
        use inter
        implicit none

        integer, intent(out)          :: ifail
        integer, intent(in)           :: itype
        integer, intent(in)           :: ncol,nrow
        integer, intent(inout)        :: iarray(ncol,nrow)
        real,    intent(inout)        :: rarray(ncol,nrow)
        character (len=*), intent(in) :: filename
        integer, intent(in)           :: ireport

        integer                       :: iunit,n2,n1,ierr,jfail,nncol,nnrow,irow,icol
        integer                       :: kstp,kper,mcol,mrow,mlay
        real                          :: pertim,totim
        character (len=1)             :: aformat
        character (len=4)             :: aext
        character (len=16)            :: text
        character (len=200)           :: bfile

        ifail=0
        call addquote(filename,bfile)
        call casetrans(bfile,'lo')            ! Not for UNIX
        n2=len_trim(filename)
        n1=n2-3
        if(n1.lt.1)n1=1
        aext=filename(n1:n2)
        call casetrans(aext,'lo')
        if(itype.eq.1)then
          if(aext.eq.'.inu') then
            aformat='u'
          else
            aformat='f'
          end if
        else if(itype.eq.2)then
          if(aext.eq.'.reu') then
            aformat='u'
          else
            aformat='f'
          end if
        end if
        iunit=nextunit()
        if(aformat.eq.'f') then
	  open(unit=iunit,file=filename,status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(bfile)
20          format(' Cannot open array file ',a,'.')
            go to 9890
          end if
          if(headerspec.eq.'yes') then
            read(iunit,'(a)',err=9350,end=9500) cline
            call linesplit(jfail,2)
            if(jfail.ne.0) go to 9350
            nncol=char2int(jfail,1)
            if(jfail.ne.0) go to 9350
            nnrow=char2int(jfail,2)
            if(jfail.ne.0) go to 9350
            if((nncol.ne.ncol).or.(nnrow.ne.nrow)) go to 9400
          end if
          if(itype.eq.1)then
            do irow=1,nrow
              read(iunit,*,err=9450,end=9500) (iarray(icol,irow),icol=1,ncol)
            end do
          else if(itype.eq.2)then
            do irow=1,nrow
              read(iunit,*,err=9450,end=9500) (rarray(icol,irow),icol=1,ncol)
            end do
          end if
        else
          open(unit=iunit,file=filename,status='old',form='binary',err=9550)
          if(itype.eq.1)then
            read(iunit,err=9600,end=9600)
            read(iunit,err=9600,end=9600) ((iarray(icol,irow),icol=1,ncol),irow=1,nrow)
          else
!            read(iunit,err=9600,end=9600) kstp,kper,pertim,totim,text,mcol,mrow,mlay
            read(iunit,err=9600,end=9600)
            read(iunit,err=9600,end=9600) ((rarray(icol,irow),icol=1,ncol),irow=1,nrow)
          end if
        end if
        close(unit=iunit)
        if(ireport.ne.0)then
          write(6,30) trim(bfile)
30        format(' - file ',a,' read ok.')
        end if
        return

9350    write(amessage,9360) trim(bfile)
9360    format(' Error reading NCOL NROW header in array file ',a,'.')
        go to 9890
9400    write(amessage,9410) trim(bfile)
9410    format(' NCOL NROW header in file ',a,' does not supply same number of ', &
        'columns and rows as in grid specificatoin file.')
        go to 9890
9450    write(amessage,9460) trim(bfile)
9460    format(' Error reading array from file ',a,'.')
        go to 9890
9500    write(amessage,9510) trim(bfile)
9510    format(' Unexpected end to array file ',a,'.')
        go to 9890
9550    write(amessage,9560) trim(bfile)
9560    format(' Cannot open unformatted array file ',a,'.')
        go to 9890
9600    write(amessage,9610) trim(bfile)
9610    format(' Error reading unformatted array file ',a,'.')
        go to 9890

9890    continue
        ifail=1
        return

end subroutine get_array



real function interp_value(ifail,iaction,ncol,nrow,nlay,jcol,jrow,jlay,rarray,bot,   &
                           unassigned_real,aline,afile)

! -- Function INTERP_VALUE undertakes vertical interpolation between already-assigned
!    array values on the basis of layer elevations.

       use defn
       use inter
       implicit none

       integer, intent(out)            :: ifail
       integer, intent(in)             :: iaction
       integer, intent(in)             :: ncol,nrow,nlay
       integer, intent(in)             :: jcol,jrow,jlay
       real, intent(in)                :: rarray(ncol,nrow,nlay)
       real, intent(in)                :: bot(ncol,nrow,0:nlay)
       real, intent(in)                :: unassigned_real
       character (len=*), intent(in)   :: aline,afile

       integer                         :: ilay
       real                            :: elevtop,valuetop,elevbot,elevmid,rtemp,  &
                                          valuebot
       character (len=5)               :: acol,arow,alay

       ifail=0

       if(jlay.eq.1)then
         elevtop=-1.1e36
       else
         do ilay=jlay-1,1,-1
           if(rarray(jcol,jrow,ilay).gt.unassigned_real) then
             elevtop=0.5*(bot(jcol,jrow,ilay-1)+bot(jcol,jrow,ilay))
             valuetop=rarray(jcol,jrow,ilay)
             go to 50
           end if
         end do
         elevtop=-1.1e36
50       continue
       end if
       if(jlay.eq.nlay)then
         elevbot=-1.1e36
       else
         do ilay=jlay+1,nlay
           if(rarray(jcol,jrow,ilay).gt.unassigned_real)then
             elevbot=0.5*(bot(jcol,jrow,ilay-1)+bot(jcol,jrow,ilay))
             valuebot=rarray(jcol,jrow,ilay)
             go to 60
           end if
         end do
         elevbot=-1.1e36
       end if
60     continue

       if((elevbot.lt.-1.0e36).and.(elevtop.lt.-1.0d36))then
         call num2char(jcol,acol)
         call num2char(jrow,arow)
         call num2char(jlay,alay)
         write(amessage,70) trim(aline),trim(afile),trim(arow),trim(acol),trim(alay)
70       format(' Invalid INTERP_ARITH or INTERP_GEOM instruction on line ',a,   &
         ' of file ',a,': for (COL,ROW) ',        &
         '(',a,',',a,') there is no previously-assigned grid cell element either ',&
         'above or below target cell in layer ',a,'.')
         go to 9890
       end if
       if(elevbot.lt.-1.0e36)then
         interp_value=valuetop
       else if(elevtop.lt.-1.0e36)then
         interp_value=valuebot
       else
         if(elevtop.le.elevbot)then
           call num2char(jcol,acol)
           call num2char(jrow,arow)
           write(amessage,80) trim(arow),trim(acol)
80         format(' Zero or negative grid cell thickness encountered at ',  &
           '(ROW,COL) (',a,','a,') of finite difference grid.')
           go to 9890
         end if
         elevmid=0.5*(bot(jcol,jrow,jlay-1)+bot(jcol,jrow,jlay))
         if(iaction.eq.1)then
           interp_value=valuebot+(valuetop-valuebot)*(elevmid-elevbot)/  &
           (elevtop-elevbot)
         else if(iaction.eq.2)then
           if((valuebot.le.0.0).or.(valuetop.le.0.0))then
             call num2char(jcol,acol)
             call num2char(jrow,arow)
             write(amessage,90) trim(arow),trim(acol),trim(aline),trim(afile)
90           format(' Zero or negative property array element encountered in grid ', &
             'at (ROW,COL) (',a,',',a,') preventing use of INTERP_GEOM ', &
             'function at line ',a,' of file ',a,'.')
             go to 9890
           end if
           valuebot=log(valuebot)
           valuetop=log(valuetop)
           rtemp=valuebot+(valuetop-valuebot)*(elevmid-elevbot)/  &
           (elevtop-elevbot)
           interp_value=exp(rtemp)
         end if
       end if
       return

9890   continue
       ifail=1
       return

end function interp_value


! - Notes

! 2. Test with/without NCOL/NROW header on input and output files.
! 4. Re-compile with checking tracing etc switched off.
! 6. Point out in manual that if a zone is not present, then the pertinent
!    operation is simply ignored.
