program asenproc

! -- Program ASENPROC processes a sensitivity file written by the adjoint state version of MODFLOW 2005, writing
!    a PEST external derivatives file.

	use defn
	use inter

	implicit none

        type dispar_pestpar
          integer                            :: ndim              ! Dimension of parameter translation table
          integer                            :: itrans            ! Flag for if natural values or logs are summed to form PEST parameter.
          integer                            :: ntab              ! Actual number of entries in table
          integer                            :: mlay              ! Maximum layers to which a distributed parameter pertains
          integer                            :: mcell             ! Number of lines in translation table
          integer                            :: mpar              ! Number of PEST parameters represented by distributed parameter
          integer, dimension(:), pointer     :: par_index         ! PEST parameter number corresponding to distributed parameter number
          integer, dimension(:), pointer     :: lay_index         ! Parameter layer number corresponding to model layer number
          integer, dimension(:,:,:), pointer :: start_index       ! Starting location in translation table for a distributed parameter cell
          integer, dimension(:,:,:), pointer :: num_index         ! Number of entries in translation table for a distributed parameter cell
          integer, dimension(:), pointer     :: start_cell_index  ! Starting location where storage is per cell rather than per layer
          integer, dimension(:), pointer     :: num_cell_index    ! Number of entries where storage is per cell rather than per layer
          integer, dimension(:), pointer     :: icellno           ! The cell number in the model grid of the pertinent distributed parameter cell
          integer, dimension(:), pointer     :: parnum            ! pest parameter index of translation table
          real, dimension(:), pointer        :: parcontrib        ! pest parameter contribution of translation table
        end type dispar_pestpar

! -- Static variables.

        type (modelgrid)  :: gridspec
        integer           :: ifail,i,ierr,ibeg,iend,ndim,icount,itemp,maxent,oldmaxent,itab,nument,j, &
                             istart,ifin,mcell,icell,jcellno,iline,jline
        integer           :: ncol,nrow,nlay,ilay,icol,irow,ncr,mlay,jlay
        integer           :: nobs,npar,ipar,ndispar,idispar,mpar,jpar,ndisobs,idisobs,jobs,iobs,mobs
        integer           :: pestunit,senunit,dppunit,derunit,mmfunit
        integer           :: iaform,pvalueout,itrans
        real              :: senthresh,rnul,rtemp,rval
        double precision  :: rsen

        character (len=1)   :: asf
        character (len=5)   :: apcf
        character (len=10)  :: alay,anum1,anum2,acol,arow,aline
        character (len=10)  :: adispar10,aapar10
        character (len=12)  :: aaobs12
        character (len=12)  :: aapar
        character (len=20)  :: aaobs,adisobs
        character (len=100) :: astring
        character (len=200) :: pestfile,senfile,dppfile,derfile,afile

! -- Allocatable variables.

        type (dispar_pestpar), allocatable ::  dpp(:)
        integer, allocatable            :: iparused(:),iobsused(:),parnum(:)
        real, allocatable               :: parcontrib(:),x(:,:),rowvec(:),pval(:)
        character (len=12), allocatable :: apar(:),adispar(:)
        character (len=20), allocatable :: aobs(:)


        write(amessage,5)
5       format(' Program ASENPROC processes a sensitivity file written by the adjoint state version ',  &
        'of MODFLOW 2005, writing a PEST external derivatives file on the basis of data contained in this file.')
        call write_message(leadspace='yes',endspace='yes')

! -- Initialisation

        oldmaxent=0

! -- The grid specification file is read.

        call readfig(gridspec%specfile)
10      call spec_open(ifail,gridspec)
        if(ifail.ne.0) go to 9900
        if(escset.eq.1) go to 9900
        call read_spec_dim(ifail,gridspec)
        if(ifail.ne.0) go to 9900
        call close_spec_file(gridspec,ok='yes')
        ncol=gridspec%ncol
        nrow=gridspec%nrow
        ncr=ncol*nrow

20      write(6,25,advance='no')
25      format(' How many layers in the model? ')
        i=key_read(nlay)
        if(escset.ne.0)then
          escset=0
          write(6,*)
          go to 10
        else if(i.eq.-1) then
          go to 20
        else if(i.ne.0) then
          write(6,26)
26        format(' Illegal input  - try again.')
          go to 20
        end if
        if(nlay.le.0) then
          write(6,26)
          go to 20
        end if

! -- The PEST control file is read for parameters and observations.

        write(6,*)
30      call open_input_file(ifail,          &
        ' Enter name of PEST control file: ',pestfile,pestunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) then
          write(6,*)
          go to 20
        end if
        read(pestunit,'(a)',err=9000,end=9000) apcf
        call casetrans(apcf,'lo')
        if(apcf(1:4).ne.'pcf ') go to 9000
        read(pestunit,*,err=9000,end=9000)
        read(pestunit,*,err=9000,end=9000)
        read(pestunit,*,err=9000,end=9000) npar,nobs
        if((npar.le.0).or.(nobs.le.0)) go to 9000
        allocate(apar(npar),aobs(nobs),stat=ierr)
        if(ierr.ne.0) go to 9200
        do
          read(pestunit,'(a)',err=9000,end=9000) cline
          if(cline(1:1).eq.'*')then
            call casetrans(cline,'lo')
            if(index(cline,'parameter da').ne.0) exit
          end if
        end do
        do ipar=1,npar
          read(10,*,err=9000,end=9000) apar(ipar)
          call casetrans(apar(ipar),'lo')
        end do
        do
          read(pestunit,'(a)',err=9000,end=9000) cline
          if(cline(1:1).eq.'*')then
            call casetrans(cline,'lo')
            if(index(cline,'observation da').ne.0) exit
          end if
        end do
        do iobs=1,nobs
          read(10,*,err=9000,end=9000) aobs(iobs)
          call casetrans(aobs(iobs),'lo')
        end do
        close(unit=pestunit)
        write(6,50) trim(pestfile)
50      format(' - file ',a,' read ok.')

! -- Arrays are filled in which elements are "ticked off" as parameters and observations
!    are mentioned in the the adjoint modflow sensitivity file.

        allocate(iparused(npar),iobsused(nobs),stat=ierr)
        if(ierr.ne.0) go to 9200
        iparused=0        ! an array
        iobsused=0        ! an array

! -- The MODFLOW sensitivity file is now openend and the names of the distributed parameters
!    that we will be dealing with obtained.

        write(6,*)
100     continue
        write(6,110,advance='no')
110     format(' Enter name of MODFLOW distributed parameter sensitivity file: ')
        read(5,'(a)') afile
        if(afile.eq.' ') go to 100
        if(index(eschar,afile(1:2)).ne.0) then
          deallocate(apar,aobs,iparused,iobsused,stat=ierr)
          if(ierr.ne.0) go to 9050
          write(6,*)
          go to 30
        end if
        ibeg=1
        iend=len_trim(afile)
        call getfile(ifail,afile,senfile,ibeg,iend)
        if(ifail.ne.0) go to 100
120     write(6,130,advance='no')
130     format(' Is this a formatted or unformatted file?  [f/u]: ')
        read(5,'(a)') asf
        if(asf.eq.' ') go to 120
        call casetrans(asf,'lo')
        if(asf.eq.'e') then
          write(6,*)
          go to 100
        end if
        if((asf.ne.'f').and.(asf.ne.'u')) go to 120
        senunit=nextunit()
        if(asf.eq.'f')then
          open(unit=senunit,file=senfile,status='old',iostat=ierr)
        else
          open(unit=senunit,file=senfile,status='old',form='binary',iostat=ierr)
        end if
        if(ierr.ne.0) then
          write(6,*)
          if(asf.eq.'f') then
            write(6,140) trim(senfile)
140         format(' Cannot open formatted file ',a,' -  try again.')
          else
            write(6,141) trim(senfile)
141         format(' Cannot open unformatted file ',a,' - try again.')
          end if
          write(6,*)
          go to 100
        end if
        if(asf.eq.'f')then
          read(senunit,*,iostat=ierr) iaform,pvalueout
        else
          read(senunit,iostat=ierr) iaform,pvalueout
        end if
        if(ierr.ne.0)then
          write(amessage,150) trim(senfile)
150       format(' Error reading value for formatting control variables IAFORM and PVALUEOUT from first ', &
          'line of MODFLOW distributed parameter sensitivity file ',a,'.')
          go to 9890
        end if
        if(iaform.ne.3)then
          write(amessage,160) trim(senfile)
160       format(' IAFORM for MODFLOW sensitivity output formatting must be 3 in file ',a,'.')
          go to 9890
        end if
        if(pvalueout.ne.1)then
          write(amessage,161) trim(senfile)
161       format(' PVALUEOUT for MODFLOW sensitivity output formatting must be 1 in file ',a,'.')
          go to 9890
        end if
        if(asf.eq.'f')then
          read(senunit,*,iostat=ierr) senthresh
        else
          read(senunit,iostat=ierr) senthresh
        end if
        if(ierr.ne.0)then
          write(amessage,170) trim(senfile)
170       format(' Cannot read sensitivity threshold from second record in file ',a,'.')
          go to 9890
        end if
        if(asf.eq.'f')then
          read(senunit,*,iostat=ierr) ndispar
        else
          read(senunit,iostat=ierr) ndispar
        end if
        if(ierr.ne.0)then
          write(amessage,180) trim(senfile)
180       format(' Cannot read number of distributed parameters from third record of file ',a,'.')
          go to 9890
        end if
        if(ndispar.le.0)then
          write(amessage,181) trim(senfile)
181       format(' Number of distributed parameters must be positive in third record ', &
          'of file ',a,'.')
          go to 9890
        end if
        allocate(adispar(ndispar),stat=ierr)
        if(ierr.ne.0) go to 9200
        jline=3
        do idispar=1,ndispar
          jline=jline+1
          if(asf.eq.'f')then
            read(senunit,*,iostat=ierr) adispar(idispar)
          else
            read(senunit,iostat=ierr) adispar10
            adispar(idispar)=adjustl(adispar10)
          end if
          if(ierr.ne.0)then
            write(amessage,190) trim(senfile)
190         format(' Error reading names of distirbuted parameters from file ',a,'.')
            go to 9890
          end if
        end do

! -- Distributed-to-PEST-parameter files are now read.

        allocate(dpp(ndispar),stat=ierr)
        if(ierr.ne.0) go to 9200
        idispar=0
195     continue                          ! beginning of loop
          idispar=idispar+1
          if(idispar.gt.ndispar) go to 300
200       continue
          write(6,*)
          astring=' Enter name of distributed-to-PEST-parameter file for "'//   &
          trim(adispar(idispar))//'":'
          call open_input_file(ifail,astring,dppfile,dppunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0) then
            if(idispar.eq.1) then
              deallocate(dpp,stat=ierr)
              if(ierr.ne.0) go to 9050
              deallocate(adispar,stat=ierr)
              if(ierr.ne.0) go to 9050
              close(unit=senunit)
              escset=0
              write(6,*)
              go to 100
            else
              idispar=idispar-1
              deallocate(dpp(idispar)%par_index,dpp(idispar)%parnum,dpp(idispar)%parcontrib,stat=ierr)
              if(ierr.ne.0) go to 9050
              nullify(dpp(idispar)%par_index,dpp(idispar)%parnum,dpp(idispar)%parcontrib)
              if(dpp(idispar)%mlay.gt.0)then
                deallocate(dpp(idispar)%lay_index,dpp(idispar)%start_index,  &
                           dpp(idispar)%num_index,stat=ierr)
                if(ierr.ne.0) go to 9050
                nullify(dpp(idispar)%lay_index,dpp(idispar)%start_index,  &
                           dpp(idispar)%num_index)
              else
                deallocate(dpp(idispar)%icellno,dpp(idispar)%start_cell_index,   &
                           dpp(idispar)%num_cell_index,stat=ierr)
                if(ierr.ne.0) go to 9050
                nullify(dpp(idispar)%icellno,dpp(idispar)%start_cell_index,   &
                           dpp(idispar)%num_cell_index)
              end if
              escset=0
              write(6,*)
              go to 200
            end if
          end if
          iline=1
!          read(dppunit,*,iostat=ierr) ndim
209       read(dppunit,'(a)',iostat=ierr) cline
          if(ierr.ne.0)then
            write(amessage,210) trim(dppfile)
210         format(' Cannot read translation table dimensions and/or transformation status from first line of file ',a,'.')
            go to 9890
          end if
          call linesplit(ifail,2)
          if(ifail.ne.0)then
            write(amessage,214) trim(dppfile)
214         format(' Two entries are required on first line of file ',a,'.')
            go to 9890
          end if
          ndim=char2int(ifail,1)
          if(ifail.ne.0)then
            write(amessage,216) trim(dppfile)
216         format(' Cannot read translation table dimensions from first line of file ',a,'.')
            go to 9890
          end if
          if(ndim.le.0)then
            write(amessage,211) trim(dppfile)
211         format(' Illegal NDIM variable at first line of file ',a,'.')
            go to 9890
          end if
          itrans=char2int(ifail,2)
          if(ifail.ne.0)then
            write(amessage,212) trim(dppfile)
212         format(' Cannot read ITRANS variable from first line of file ',a,'.')
            go to 9890
          end if
          if((itrans.ne.0).and.(itrans.ne.1))then
            write(amessage,213) trim(dppfile)
213         format(' ITRANS must be 0 or 1 on first line of file ',a,'.')
            go to 9890
          end if
          dpp(idispar)%ndim=ndim
          allocate(dpp(idispar)%parnum(ndim),dpp(idispar)%parcontrib(ndim),stat=ierr)
          if(ierr.ne.0) go to 9200
          dpp(idispar)%itrans=itrans
          iline=2
          read(dppunit,*,iostat=ierr) mlay
          if(ierr.ne.0)then
            write(amessage,220) trim(dppfile)
220         format(' Cannot read "number of represented layers" variable from second line of file ',a,'.')
            go to 9890
          end if
          if(mlay.eq.0)then
            write(amessage,215) trim(dppfile)
215         format(' Value for "number of represented layers" variable must not be zero at line 2 of file ',a,'.')
            go to 9890
          end if
          if(mlay.gt.0)then
            dpp(idispar)%mlay=mlay
            dpp(idispar)%mcell=0
            allocate(dpp(idispar)%lay_index(nlay),dpp(idispar)%start_index(ncol,nrow,mlay),   &
                     dpp(idispar)%num_index(ncol,nrow,mlay),stat=ierr)
            if(ierr.ne.0) go to 9200
            dpp(idispar)%start_index=0         ! an array
            dpp(idispar)%num_index=0           ! an array
            dpp(idispar)%lay_index=0            ! an array
            icount=0
            do ilay=1,mlay
              iline=iline+1
              read(dppunit,*,iostat=ierr) itemp
              if(ierr.ne.0) then
                call num2char(ilay,alay)
                call num2char(iline,aline)
                write(amessage,221) trim(alay),trim(aline),trim(dppfile)
221             format(' Error reading layer number ',a,' from line ',a,' file ',a,'.')
                go to 9890
              end if
              if((itemp.le.0).or.(itemp.gt.nlay))then
                call num2char(ilay,alay)
                call num2char(iline,aline)
                write(amessage,230) trim(alay),trim(aline),trim(dppfile)
230             format(' Layer number ',a,' is out of range on line ',a,' of file ',a,'.')
                go to 9890
              end if
              icount=icount+1
              if(dpp(idispar)%lay_index(itemp).ne.0)then
                 call num2char(iline,aline)
                 write(amessage,240) trim(aline),trim(dppfile)
240              format(' Layer number is repeated at line ',a,' of file ',a,'.')
                 go to 9890
              else
                dpp(idispar)%lay_index(itemp)=icount
              end if
            end do
          else
            dpp(idispar)%mlay=0
            mcell=-mlay
            dpp(idispar)%mcell=mcell
            allocate(dpp(idispar)%icellno(mcell),dpp(idispar)%num_cell_index(mcell),    &
            dpp(idispar)%start_cell_index(mcell),stat=ierr)
            if(ierr.ne.0) go to 9200
          end if
          iline=iline+1
          read(dppunit,*,iostat=ierr)mpar
          if(ierr.ne.0)then
            call num2char(iline,aline)
            write(amessage,250) trim(aline),trim(dppfile)
250         format(' Error reading "number of PEST parameters" variable from line ',a,' of file ',a,'.')
            go to 9890
          end if
          if(mpar.le.0)then
            call num2char(iline,aline)
            write(amessage,260) trim(aline),trim(dppfile)
260         format(' "Number of PEST parameters" variable must be greater than zero at line ',a,   &
            ' of file ',a,'.')
            go to 9890
          end if
          dpp(idispar)%mpar=mpar
          allocate(dpp(idispar)%par_index(mpar),stat=ierr)
          if(ierr.ne.0) go to 9890
          jpar=1
          do ipar=1,mpar
            iline=iline+1
            read(dppunit,*,iostat=ierr) aapar
            call casetrans(aapar,'lo')
            call whichone(ifail,npar,jpar,apar,aapar)
            if(ifail.ne.0)then
              write(amessage,270) trim(aapar),trim(dppfile),trim(pestfile)
270           format(' Parameter "',a,'" featured in file ',a,' is not a PEST parameter ',  &
              'featured in file ',a,'.')
              go to 9890
            end if
            dpp(idispar)%par_index(ipar)=jpar
            if(iparused(jpar).ne.0)then
              itemp=iparused(jpar)
              if(dpp(itemp)%itrans.ne.dpp(idispar)%itrans)then
                write(amessage,275) trim(apar(jpar))
275             format(' The same PEST parameter "',a,'" is cited in more than one ', &
                'distributed-to-PEST-parameter file. Both of these files must have the same ', &
                'ITRANS variable, signifying combinations of cell values or log cell values ', &
                'in both cases.')
                go to 9890
              end if
            end if
            iparused(jpar)=idispar
          end do
          iline=iline+1
          read(dppunit,*,iostat=ierr)maxent
          if(ierr.ne.0)then
            call num2char(iline,aline)
            write(amessage,280) trim(aline),trim(dppfile)
280         format(' Error reading "maximum number of translation table entries per line" variable ', &
            'from line ',a,' of file ',a,'.')
            go to 9890
          end if
          if(maxent.le.0)then
            call num2char(iline,aline)
            write(amessage,281) trim(aline),trim(dppfile)
281         format(' "Maximum number of translation table entries per line" variable on line ',  &
            a,' of file ',a,' must be positive.')
            go to 9890
          end if
          if(maxent.gt.oldmaxent)then
            if(oldmaxent.gt.0)then
              deallocate(parnum,parcontrib,stat=ierr)
              if(ierr.ne.0) go to 9050
            end if
            allocate(parnum(maxent),parcontrib(maxent),stat=ierr)
            if(ierr.ne.0) go to 9200
            oldmaxent=maxent
          end if
          itab=0
          if(mlay.gt.0)then
            i=0
            do
              i=i+1
              read(dppunit,*,err=9250,end=309) icol,irow,ilay,nument,(parnum(j),parcontrib(j),j=1,nument)
              j=dpp(idispar)%lay_index(ilay)
              if(j.eq.0)then
                call num2char(i,anum1)
                call num2char(ilay,anum2)
                write(amessage,282) trim(anum1),trim(dppfile),trim(anum2)
282             format(' In record number ',a,' of the translation table in file ',a,' reference ',  &
                'is made to model layer ',a,'; this layer was not previously cited in that file as ',  &
                'being included in this super parameter.')
                go to 9890
              end if
              dpp(idispar)%start_index(icol,irow,j)=itab+1
              dpp(idispar)%num_index(icol,irow,j)=nument
              if(itab+nument.gt.ndim)then
                write(amessage,290) trim(dppfile)
290             format(' The parameter translation table in file ',a,' has more total entries ',  &
                '(i.e. records times entries on each line) than ',  &
                'indicated by the NDIM variable at the top of the file. See documentation of this ', &
                'program for meaning of this variable.')
                go to 9890
              end if
              do j=1,nument
                itab=itab+1
                dpp(idispar)%parnum(itab)=dpp(idispar)%par_index(parnum(j))
                dpp(idispar)%parcontrib(itab)=parcontrib(j)
              end do
            end do
          else
            do icell=1,mcell
              read(dppunit,*,err=9250,end=9270) icol,irow,ilay,nument,(parnum(j),parcontrib(j),j=1,nument)
              dpp(idispar)%start_cell_index(icell)=itab+1
              dpp(idispar)%num_cell_index(icell)=nument
              dpp(idispar)%icellno(icell)=(ilay-1)*ncr+(irow-1)*ncol+icol
              if(itab+nument.gt.ndim)then
                write(amessage,290) trim(dppfile)
                go to 9890
              end if
              do j=1,nument
                itab=itab+1
                dpp(idispar)%parnum(itab)=dpp(idispar)%par_index(parnum(j))
                dpp(idispar)%parcontrib(itab)=parcontrib(j)
              end do
            end do
            if(mcell.gt.1)then
              do icell=2,mcell
                if(dpp(idispar)%icellno(icell).le.dpp(idispar)%icellno(icell-1))then
                  call num2char(icell,anum1)
                  write(amessage,295) trim(dppfile),trim(anum1)
295               format(' Cells comprising the distributed parameter are not provided in ',  &
                  'increasing (COL,ROW,LAY) order (with COL cycling fastest) ', &
                  'in translation table supplied in file ',a,'; error occurs at cell number ',a,  &
                  ' in this table.')
                  go to 9890
                end if
              end do
            end if
          end if
309       continue

          dpp(idispar)%ntab=itab
          close(unit=dppunit)
          write(6,310) trim(dppfile)
310       format(' - file ',a,' read ok.')
        go to 195                           ! end of loop

300     continue

! -- The name of the PEST external derivatives file is now acquired.

        write(6,*)
291     call open_output_file(ifail, &
        ' Enter name for PEST external derivatives file: ',derfile,derunit)
	if(ifail.ne.0) go to 9900
	if(escset.ne.0) then
          idispar=ndispar
          deallocate(dpp(idispar)%par_index,dpp(idispar)%parnum,dpp(idispar)%parcontrib,stat=ierr)
          if(ierr.ne.0) go to 9050
          nullify(dpp(idispar)%par_index,dpp(idispar)%parnum,dpp(idispar)%parcontrib)
          if(dpp(idispar)%mlay.gt.0)then
            deallocate(dpp(idispar)%lay_index,dpp(idispar)%start_index,  &
                       dpp(idispar)%num_index,stat=ierr)
            if(ierr.ne.0) go to 9050
            nullify(dpp(idispar)%lay_index,dpp(idispar)%start_index,  &
                    dpp(idispar)%num_index)
          else
            deallocate(dpp(idispar)%icellno,dpp(idispar)%start_cell_index,   &
                       dpp(idispar)%num_cell_index,stat=ierr)
            if(ierr.ne.0) go to 9050
            nullify(dpp(idispar)%icellno,dpp(idispar)%start_cell_index,   &
                    dpp(idispar)%num_cell_index)
          end if
          escset=0
          go to 200
        end if

! -- If any distributed-to-PEST parameter files indicate that logs of cell sensitivities must
!    be combined to compute PEST parameter sensitivity, then a PEST-to-model message file must be found
!    in the current directory.


        itrans=0
        do idispar=1,ndispar
          if(dpp(idispar)%itrans.eq.1) itrans=1
        end do
        if(itrans.eq.1)then
          write(6,501)
501       format(/,' Reading PEST-to-model message file pest.mmf...')
          allocate(pval(npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          mmfunit=nextunit()
          open(unit=mmfunit,file='pest.mmf',status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,510)
510         format(' If any distributed-to-PEST parameter files indicate logarithmic combination ', &
            'of model parameters to PEST parameters (i.e. have an ITRANS value of 1)', &
            ', then a PEST-to-model message file (always named ', &
            '"pest.mmf") must be present in the current directory.')
            go to 9890
          end if
          read(mmfunit,'(a)',err=9600,end=9600) cline
          call casetrans(cline,'lo')
          if(index(cline,'external_der').eq.0)then
            write(amessage,520)
520         format(' First line of PEST-to-model message file pest.mmf should be "external derivatives".')
            go to 9890
          end if
          read(mmfunit,*,err=9600,end=9600) itemp
          read(mmfunit,*,err=9600,end=9600) mpar,mobs
          if(mpar.ne.npar)then
            write(amessage,530) trim(pestfile)
530         format(' Number of parameters cited in PEST-to-model message file pest.mmf does not ',  &
            'agree with number cited in PEST control file ',a,'.')
            go to 9890
          end if
          do ipar=1,npar
            read(mmfunit,*,err=9600,end=9600) aapar,pval(ipar)
            call casetrans(aapar,'lo')
            if(aapar.ne.apar(ipar))then
              write(amessage,532) trim(pestfile)
532           format(' Names and/or ordering of parameters in PEST-to-model message file pest.mmf do ', &
              'not agree with those cited in PEST control file ',a,'.')
              go to 9890
            end if
          end do
          close(unit=mmfunit)
          write(6,540)
540       format(' - PEST-to-model message file pest.mmf read ok.')
        end if

! -- Observation names are now read from the MODFLOW sensitivity file.

        write(6,315) trim(senfile)
315     format(/,' Processing sensitivity data in file ',a,'...')
        jline=jline+1
        if(asf.eq.'f')then
          read(senunit,*,iostat=ierr) ndisobs
        else
          read(senunit,iostat=ierr) ndisobs
        end if
        if(ierr.ne.0)then
          call num2char(jline,aline)
          write(amessage,320) trim(aline),trim(senfile)
320       format(' Cannot read "number-of-observations" variable from record number ',a,  &
          ' of sensitivity file ',a,'.')
          go to 9890
        end if
        jobs=1
        do idisobs=1,ndisobs
          jline=jline+1
          if(asf.eq.'f')then
            read(senunit,*,iostat=ierr) aaobs
          else
            read(senunit,iostat=ierr) aaobs12
            aaobs=adjustl(aaobs12)
          end if
          if(ierr.ne.0)then
            write(amessage,330) trim(senfile)
330         format(' Error reading names of observations from file ',a,'.')
            go to 9890
          end if
          call casetrans(aaobs,'lo')
          call whichone(ifail,nobs,jobs,aobs,aaobs)
          if(ifail.ne.0)then
            call num2char(jline,aline)
            write(amessage,340) trim(aaobs),trim(aline),trim(senfile),trim(pestfile)
340         format(' Observation "',a,'" featured on line ',a,' of file ',a,' is not a PEST observation ',  &
            'featured in file ',a,'.')
            go to 9890
          end if
          if(iobsused(jobs).ne.0)then
            write(amessage,341) trim(aaobs),trim(senfile)
341         format(' The same observation, viz "',a,'" is cited more than once in file ',a,'.')
            go to 9890
          end if
          iobsused(jobs)=idisobs
        end do

! -- A PEST sub-Jacobian matrix is allocated.

        allocate(x(ndisobs,npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        x = 0.0  ! an array  We should make this double precision.
        do idisobs=1,ndisobs
          do idispar=1,ndispar
            icell=0
            if(asf.eq.'f')then
              read(senunit,*,iostat=ierr) nument,aapar,aaobs,rtemp,rtemp
            else
              read(senunit,iostat=ierr) nument,aapar10,aaobs12,rtemp,rtemp  ! make sure that tom's character lengths agree
              aapar=adjustl(aapar10)
              aaobs=adjustl(aaobs12)
            end if
            if(ierr.ne.0)then
              do i=1,nobs
                if(iobsused(i).eq.idisobs)then
                  adisobs=aobs(i)
                  go to 351
                end if
              end do
351           continue
              write(amessage,350) trim(adisobs),trim(adispar(idispar)),trim(senfile)
350           format(' Error reading sensitivity header for observation "',a,   &
              '" and distributed parameter "',a,'" from sensitivity file ',a,'.')
              go to 9890
            end if
            do i=1,nument
              if(asf.eq.'f')then
                read(senunit,*,err=9350) icol,irow,ilay,rsen,rval
              else
                read(senunit,err=9350) icol,irow,ilay,rsen,rval
              end if
              if(rsen.ne.0.0d0)then
                if(dpp(idispar)%itrans.eq.1)then
                  rsen=rval*rsen
                end if
                if(dpp(idispar)%mlay.gt.0)then
                  jlay=dpp(idispar)%lay_index(ilay)
                  istart = dpp(idispar)%start_index(icol,irow,jlay)
                  ifin = istart+dpp(idispar)%num_index(icol,irow,jlay)-1
                  if(istart.eq.0) go to 9500
                else
                  jcellno=(ilay-1)*ncr+(irow-1)*ncol+icol
                  do j=icell+1,dpp(idispar)%mcell
                    if(dpp(idispar)%icellno(j).eq.jcellno) go to 354  ! Note that increasing ncol,nrow,lnay order is assumed
                  end do
                  go to 9500
354               istart = dpp(idispar)%start_cell_index(j)
                  ifin = istart+dpp(idispar)%num_cell_index(j)-1
                  icell=j
                end if
                do j=istart,ifin
                  x(idisobs,dpp(idispar)%parnum(j))=           &
                  x(idisobs,dpp(idispar)%parnum(j))+dpp(idispar)%parcontrib(j)*rsen
                end do
              end if
            end do
          end do
        end do
        if(itrans.eq.1)then
          do ipar=1,npar
            itemp=iparused(ipar)
            if(dpp(itemp)%itrans.eq.1)then
              if(pval(ipar).le.0.0d0)then
                write(amessage,353) trim(apar(ipar))
353             format(' Parameter "',a,'" is cited in a distributed-to-PEST-parameter file ', &
                'with an ITRANS value of 1, indicating amalgamation of log sensitivities rather than ', &
                'native sensitivities. However its value is supplied as zero or negative in PEST-to-model ', &
                'message file pest.mmf.')
                go to 9890
              end if
              rtemp=1.0/pval(ipar)
              do idisobs=1,ndisobs
                x(idisobs,ipar)=x(idisobs,ipar)*rtemp
              end do
            end if
          end do
        end if
        close(unit=senunit)
        write(6,480) trim(senfile)
480     format(' - file ',a,' read ok.')

! -- The X matrix has been filled; now it must be written out.

        allocate(rowvec(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        rnul=-1.11e33
        rowvec = rnul                ! an array
        write(derunit,450) npar,nobs
450     format(2i6)
        do iobs=1,nobs
          jobs=iobsused(iobs)
          if(jobs.eq.0)then
            write(derunit,460) (rnul,ipar=1,npar)
460         format(8(1x,1pg14.7))
          else
            do ipar=1,npar
              jpar=iparused(ipar)
              if(jpar.eq.0)then
                rowvec(ipar)=rnul
              else
                rowvec(ipar)=x(jobs,ipar)
              end if
            end do
            write(derunit,460) (rowvec(ipar),ipar=1,npar)
          end if
        end do

        close(unit=derunit)
        write(6,470) trim(derfile)
470     format(' - file ',a,' written ok.')


        go to 9900



9000    write(amessage,9010) trim(pestfile)
9010    format(' Error encountered in reading PEST control file ',a,'; check this file ', &
        'with PESTCHEK.')
        go to 9890

9050    write(amessage,9060)
9060    format(' Error in memory de-allocation.')
        go to 9890

9100    write(amessage,9110) trim(dppfile)
9110    format(' Error encountered while reading distributed-to-pest-parameter file ',a,'.')
        go to 9890

9150    write(amessage,9110) trim(dppfile)
9160    format(' Premature end encountered to distributed-to-pest-parameter file ',a,'.')
        go to 9890

9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890

9250    call num2char(i,anum1)
        write(amessage,9260) trim(anum1),trim(dppfile)
9260    format(' Error reading on or about entry number ',a, ' from parameter translation table ',  &
        'in file ',a,'.')
        go to 9890

9270    write(amessage,9280) trim(dppfile)
9280    format(' Unexpected end encountered to parameter translation table in file ',a,'.')
        go to 9890

9350    continue
        do iobs=1,nobs
          if(iobsused(iobs).eq.idisobs)then
            adisobs=aobs(iobs)
            go to 9351
          end if
        end do
9351    continue
        call num2char(i,anum1)
        write(amessage,9360) trim(anum1),trim(adisobs),trim(adispar(idispar)),trim(senfile)
9360    format(' Error reading sensitivity table entry number ',a,' for observation "',a, &
        '" and distributed parameter "',a,'" in sensitivity file ',a,'.')
        go to 9890

9500    continue
        call num2char(icol,acol)
        call num2char(irow,arow)
        call num2char(ilay,alay)
        do i=1,nobs
          if(iobsused(i).eq.idisobs)then
            adisobs=aobs(i)
            go to 9551
          end if
        end do
9551    continue
        write(amessage,9510) trim(acol),trim(arow),trim(alay),trim(adispar(idispar)),trim(adisobs), &
        trim(senfile)
9510    format(' Cell (column,row,layer) (',a,',',a,',',a') for which a sensitivity is provided ', &
        'for distributed parameter ',a,' for observation ',a,' in sensitiivty file ', &
        a,' has not been associated with that super parameter in a super parameter definition file.')
        go to 9890

9600    write(amessage,9610)
9610    format(' Error or premature end encountered to PEST-to-model message file pest.mmf.')
        go to 9890

9890    call write_message(leadspace='yes')
9900    call close_files

        if(allocated(dpp))then
          do idispar=1,ndispar
!            deallocate(dpp(idispar)%par_index,stat=ierr)
!            deallocate(dpp(idispar)%lay_index,stat=ierr)
!            deallocate(dpp(idispar)%start_index,stat=ierr)
!            deallocate(dpp(idispar)%num_index,stat=ierr)
!            deallocate(dpp(idispar)%start_cell_index,stat=ierr)
!            deallocate(dpp(idispar)%num_cell_index,stat=ierr)
!            deallocate(dpp(idispar)%icellno,stat=ierr)
!            deallocate(dpp(idispar)%parnum,stat=ierr)
!            deallocate(dpp(idispar)%parcontrib,stat=ierr)
            if(associated(dpp(idispar)%par_index))        nullify(dpp(idispar)%par_index)
            if(associated(dpp(idispar)%lay_index))        nullify(dpp(idispar)%lay_index)
            if(associated(dpp(idispar)%start_index))      nullify(dpp(idispar)%start_index)
            if(associated(dpp(idispar)%num_index))        nullify(dpp(idispar)%num_index)
            if(associated(dpp(idispar)%start_cell_index)) nullify(dpp(idispar)%start_cell_index)
            if(associated(dpp(idispar)%num_cell_index))   nullify(dpp(idispar)%num_cell_index)
            if(associated(dpp(idispar)%icellno))          nullify(dpp(idispar)%icellno)
            if(associated(dpp(idispar)%parnum))           nullify(dpp(idispar)%parnum)
            if(associated(dpp(idispar)%parcontrib))       nullify(dpp(idispar)%parcontrib)
          end do
          deallocate(dpp,stat=ierr)
        end if

        if(allocated(iparused))   deallocate(iparused,stat=ierr)
        if(allocated(iobsused))   deallocate(iobsused,stat=ierr)
        if(allocated(parnum))     deallocate(parnum,stat=ierr)
        if(allocated(parcontrib)) deallocate(parcontrib,stat=ierr)
        if(allocated(rowvec))     deallocate(rowvec,stat=ierr)
        if(allocated(apar))       deallocate(apar,stat=ierr)
        if(allocated(adispar))    deallocate(adispar,stat=ierr)
        if(allocated(aobs))       deallocate(aobs,stat=ierr)
        if(allocated(pval))       deallocate(pval,stat=ierr)

	write(6,*)

end program asenproc

subroutine whichone(ifail,npar,ipar,apar,tpar)

! -- Subroutine whichone locates a string in an array. Note that both the
!    search string and the existing array of strings are assumed to be
!    in the same case.

        integer npar,ipar,i
        integer ifail
        character*(*) tpar
        character*(*) apar(npar)

        ifail=0
        if((ipar.lt.1).or.(ipar.gt.npar)) ipar=1
        if(tpar.eq.apar(ipar)) return
        if(ipar.ne.npar)then
          do 20 i=ipar+1,npar
          if(tpar.eq.apar(i))then
            ipar=i
            return
          end if
20        continue
        end if
        if(ipar.ne.1)then
          do 40 i=ipar-1,1,-1
          if(tpar.eq.apar(i)) then
            ipar=i
            return
          end if
40        continue
        end if
        ifail=1
        return

 end subroutine whichone


! -- COMENTS

! -- Check with both formatted and unformatted sensitivity files.
!    Check with both cell and layer based distirbuted parameter definitions.
