!     Last change:  JD   25 Jan 2007    7:18 am
program zone2var2

! -- Program ZONE2VAR2 carries out the same tasks as ZONE2VAR1. However it receives most of its information
!    from an input file written by ZONE2VAR1.

        use defn
        use inter

	implicit none

        integer           ifail,ierr,dflag,icode,i,is,itemp,iline
        integer           zonunit,ircunit,varunit,derunit
        integer           npar,ipar,jpar,nobs,iobs
        integer           mlay,ilay,izone,iangle
        integer           nstructzone,nangle,istruct
        integer           maxsampcol,maxsamprow,maxsampang

        double precision  ln10,rtemp,rtemp2

        character*1       at,ayn
        character*15      alay,anum,aangle,aline
        character*200     aprompt
        character*200     zonfile,d2pfile,ircfile,varfile,derfile


        integer, allocatable :: layer(:)
        integer, allocatable :: ival(:)
        integer, allocatable :: iszone(:),structzonepar(:)
        integer, allocatable :: numsampcol(:),numsamprow(:),numsampang(:,:)
        integer, allocatable :: num_col_e(:,:),num_row_e(:,:)
        integer, allocatable :: num_ang_e(:,:,:)
        integer, allocatable :: num_par_row(:,:),num_par_col(:,:)
        integer, allocatable :: num_par_ang(:,:,:)

        double precision, allocatable :: pval(:),pvalkeep(:)
        double precision, allocatable :: angle(:)
        double precision, allocatable :: h_row(:,:),h_col(:,:),h_ang(:,:,:)
        double precision, allocatable :: gamma_row_e(:,:),gamma_col_e(:,:)
        double precision, allocatable :: gamma_ang_e(:,:,:)
        double precision, allocatable :: deriv_row(:,:),deriv_col(:,:)
        double precision, allocatable :: deriv_ang(:,:,:)
        double precision, allocatable :: x_row(:)

        character*12, allocatable     :: apar(:)
        character*20, allocatable     :: bobs(:)


        write(amessage,5)
5       format(' Program ZONE2VAR2 computes a parameter variogram, as well as derivatives of this ', &
        'variogram with respect to parameters, based on data supplied by ZONE2VAR1.')
        call write_message(leadspace='yes',endspace='yes')

! -- Initialization

        ln10=log(10.0d0)

! -- The ZONE2VAR1-generated data file is opened.

50      continue
60      call open_input_file(ifail, &
        ' Enter name of ZONE2VAR1-generated data file: ',zonfile,zonunit,file_format='unformatted')
        if(ifail.ne.0) go to 9900
        if(escset.ne.0) go to 9900

! -- Some of this file is now read.

        read(zonunit,err=9000,end=9050) npar,mlay

        allocate(layer(mlay),stat=ierr)
        if(ierr.ne.0) go to 9200
        read(zonunit,err=9000,end=9050) (layer(ilay),ilay=1,mlay)

        allocate(ival(npar),pval(npar),pvalkeep(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        pval=-1.0d300   ! an array
        read(zonunit,err=9000,end=9050) (ival(ipar),ipar=1,npar)
        read(zonunit,err=9000,end=9050) d2pfile

        allocate(apar(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        do ipar=1,npar
          read(zonunit,err=9000,end=9050) apar(ipar)
        end do
        read(zonunit,err=9000,end=9050) at

        write(6,*)
        ilay=0
204     continue
        ilay=ilay+1
        if(ilay.gt.mlay) go to 300
205       continue
          call num2char(layer(ilay),alay)
          aprompt=' Enter integer-real correspondence file for layer '//trim(alay)//': '
          call open_input_file(ifail,aprompt,ircfile,ircunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0)then
            escset=0
            if(ilay.eq.1)then
              deallocate(layer,pval,pvalkeep,ival,apar)
              close(unit=zonunit)
              write(6,*)
              go to 50
            else
              write(6,*)
              ilay=ilay-1
              go to 205
            end if
          end if
          iline=0
          jpar=1
          do
            iline=iline+1
            read(ircunit,'(a)',err=9100,end=250) cline
            if(cline.eq.' ') cycle
            call linesplit(ifail,2)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,220) trim(aline),trim(ircfile)
220           format(' Insufficient entries on line ',a,' of file ',a,'.')
              go to 9890
            end if
            izone=char2int(ifail,1)
            if(ifail.ne.0) go to 9100
            rtemp=char2double(ifail,2)
            if(ifail.ne.0) go to 9100
            call whichone_i(ifail,npar,jpar,ival,izone)
            if(ifail.ne.0)then
              call num2char(iline,aline)
              write(amessage,230) trim(aline),trim(ircfile),trim(d2pfile)
230           format(' Integer cited in column 1 at line ',a,' of file ',a,    &
              ' is not linked to any parameter cited in ZONMDEF file ',a,' read by ', &
              'ZONEVAR1 on its run previous to this.')
              go to 9890
            end if
            if(pval(jpar).gt.-1.0d299)then
              if(.not.equals(rtemp,pval(jpar)))then
                write(amessage,240) trim(apar(jpar)),trim(ircfile)
240             format(' Value provided to zone corresponding to parameter "',a,'" in file ',a,' differs ', &
                'from value provided in a previously-read integer-real correspondence file.')
                go to 9890
              end if
            else
              pval(jpar)=rtemp
            end if
          end do
250       continue
          close(unit=ircunit)
          write(6,270) trim(ircfile)
270       format('  - file ',a,' read ok.')
          go to 204
300     continue

        do ipar=1,npar
          if(pval(ipar).lt.-1.0d299)then
            write(amessage,310) trim(apar(ipar))
310         format(' No value has been assigned to a zone corresponding to parameter "',a,  &
            '" in any integer-real correspondence file.')
            go to 9890
          end if
        end do
        pvalkeep=pval    ! arrays

! -- The names of model output files are now obtained.

        write(6,*)
100     aprompt = ' Enter name for experimental variogram output file: '
	call open_output_file(ifail,aprompt,varfile,varunit)
        if(ifail.ne.0) go to 9900
        if(escset.ne.0)then
          escset=0
          write(6,*)
          ilay=0
          pval=-1.0d300
          go to 204
        end if

        write(6,*)
140     write(6,150,advance='no')
150     format(' Write a parameter derivatives file?  [y/n]: ')
        read(5,*) ayn
        if(ayn.eq.' ') go to 140
        call casetrans(ayn,'lo')
        if(ayn.eq.'e')then
          write(6,*)
          close(unit=varunit)
          go to 100
        else if(ayn.eq.'y')then
570       aprompt = ' Enter name for parameter derivatives file: '
          call open_output_file(ifail,aprompt,derfile,derunit)
          if(ifail.ne.0) go to 9900
          if(escset.ne.0)then
            escset=0
            write(6,*)
            go to 140
          end if
        else if (ayn.eq.'n')then
          derfile=' '
        else
          go to 140
        end if
        dflag=0
        if(derfile.ne.' ')dflag=1

! -- More of the ZONE2VAR1-generated file is read.

        read(zonunit,err=9000,end=9050) nstructzone,nangle
        read(zonunit,err=9000,end=9050) maxsampcol,maxsamprow,maxsampang

        if(nangle.ne.0)then
          allocate(angle(nangle),stat=ierr)
          if(ierr.ne.0) go to 9200
          read(zonunit,err=9000,end=9050) (angle(iangle),iangle=1,nangle)
        end if

        allocate(numsampcol(nstructzone),numsamprow(nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(numsampang(nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if
        do istruct=1,nstructzone
          read(zonunit,err=9000,end=9050) numsampcol(istruct),numsamprow(istruct)
          if(nangle.ne.0)then
            read(zonunit,err=9000,end=9050) (numsampang(iangle,istruct),iangle=1,nangle)
          end if
        end do

        allocate(h_row(0:maxsamprow,nstructzone),h_col(0:maxsampcol,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(h_ang(0:maxsampang,nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if
        do istruct=1,nstructzone
          read(zonunit,err=9000,end=9050) (h_row(i,istruct),i=0,numsamprow(istruct))
          read(zonunit,err=9000,end=9050) (h_col(i,istruct),i=0,numsampcol(istruct))
          if(nangle.ne.0)then
            do iangle=1,nangle
              read(zonunit,err=9000,end=9050) (h_ang(i,iangle,istruct),i=0,numsampang(iangle,istruct))
            end do
          end if
        end do

        allocate(iszone(nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        read(zonunit,err=9000,end=9050) (iszone(istruct),istruct=1,nstructzone)

        allocate(structzonepar(npar),stat=ierr)
        if(ierr.ne.0) go to 9200
        read(zonunit,err=9000,end=9050) (structzonepar(ipar),ipar=1,npar)


        write(6,658)
658     format(/,'  - computing experimental semivariograms...')

! -- If the geostatistical structure transform is log, then the log is taken of all parameter
!    values.

        if(at.eq.'l')then
          do ipar=1,npar
            if(pval(ipar).le.0.0d0)then
              write(amessage,662) trim(apar(ipar))
662           format(' Parameter "',a,'" has been provided with a zero or negative value. ', &
              'Log-transformation cannot apply to such parameters.')
              go to 9890
            end if
            pval(ipar)=log10(pval(ipar))
          end do
        end if

! -- More memory is allocated.

        allocate(gamma_row_e(maxsamprow,nstructzone), gamma_col_e(maxsampcol,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        allocate(num_col_e(maxsampcol,nstructzone),num_row_e(maxsamprow,nstructzone),stat=ierr)
        if(ierr.ne.0) go to 9200
        if(nangle.ne.0)then
          allocate(gamma_ang_e(maxsampang,nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
          allocate(num_ang_e(maxsampang,nangle,nstructzone),stat=ierr)
          if(ierr.ne.0) go to 9200
        end if
        gamma_row_e=0.0d0
        gamma_col_e=0.0d0
        num_col_e=0
        num_row_e=0
        if(nangle.ne.0)then
          gamma_ang_e=0.0d0
          num_ang_e=0
        end if

        if(dflag.ne.0)then
          allocate(num_par_row(maxsamprow,npar),num_par_col(maxsampcol,npar),  &
                   deriv_row(maxsamprow,npar),deriv_col(maxsampcol,npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          num_par_row=0           ! an array
          num_par_col=0           ! an array
          deriv_row=0.0d0         ! an array
          deriv_col=0.0d0         ! an array
          if(nangle.ne.0)then
            allocate(num_par_ang(maxsampang,nangle,npar),                &
                     deriv_ang(maxsampang,nangle,npar),stat=ierr)
            if(ierr.ne.0) go to 9200
            num_par_ang=0           ! an array
            deriv_ang=0.0d0         ! an array
          end if
        end if

! -- Computation of experimental variograms is carried out.

        do
          read(zonunit,err=9000,end=9050)icode,i,istruct,ipar,jpar
          if(icode.eq.0) exit
          rtemp=pval(ipar)-pval(jpar)
          rtemp2=rtemp*rtemp
          if(icode.eq.1)then
            gamma_row_e(i,istruct)=gamma_row_e(i,istruct)+rtemp2
!            num_row_e(i,istruct)=num_row_e(i,istruct)+1
            if(dflag.ne.0)then
!              num_par_row(i,ipar)=num_par_row(i,ipar)+1
!              num_par_row(i,jpar)=num_par_row(i,jpar)+1
              deriv_row(i,ipar)=deriv_row(i,ipar)+rtemp
              deriv_row(i,jpar)=deriv_row(i,jpar)-rtemp
            end if
          else if(icode.eq.2)then
            gamma_col_e(i,istruct)=gamma_col_e(i,istruct)+rtemp2
!            num_col_e(i,istruct)=num_col_e(i,istruct)+1
            if(dflag.ne.0)then
!              num_par_col(i,ipar)=num_par_col(i,ipar)+1
!              num_par_col(i,jpar)=num_par_col(i,jpar)+1
              deriv_col(i,ipar)=deriv_col(i,ipar)+rtemp
              deriv_col(i,jpar)=deriv_col(i,jpar)-rtemp
            end if
          else if(icode.lt.0)then
            iangle=-icode
            gamma_ang_e(i,iangle,istruct)=gamma_ang_e(i,iangle,istruct)+rtemp2
!            num_ang_e(i,iangle,istruct)=num_ang_e(i,iangle,istruct)+1
            if(dflag.ne.0)then
!              num_par_ang(i,iangle,ipar)=num_par_ang(i,iangle,ipar)+1
!              num_par_ang(i,iangle,jpar)=num_par_ang(i,iangle,jpar)+1
              deriv_ang(i,iangle,ipar)=deriv_ang(i,iangle,ipar)+rtemp
              deriv_ang(i,iangle,jpar)=deriv_ang(i,iangle,jpar)-rtemp
            end if
          end if
        end do

        write(6,750)
750     format('  - experimental semivariograms computed ok.')

! - More of the ZONE2VAR1 data file is read.

        do istruct=1,nstructzone
          read(zonunit) (num_row_e(i,istruct),i=1,maxsamprow)
          read(zonunit) (num_col_e(i,istruct),i=1,maxsampcol)
          if(nangle.ne.0)then
            do iangle=1,nangle
              read(zonunit) (num_ang_e(i,iangle,istruct),i=1,maxsampang)
            end do
          end if
        end do
        if(dflag.ne.0)then
          read(zonunit) ((num_par_row(i,ipar),i=1,maxsamprow),ipar=1,npar)
          read(zonunit) ((num_par_col(i,ipar),i=1,maxsampcol),ipar=1,npar)
          read(zonunit) (((num_par_ang(i,iangle,ipar),i=1,maxsampang),iangle=1,nangle),ipar=1,npar)
        end if

! -- Observed semi-variograms are now written to a file.

        write(6,840) trim(varfile)
840     format(/,'  - writing experimental variogram file ',a,'...')
        do istruct=1,nstructzone
          call num2char(iszone(istruct),anum)
          write(varunit,850) trim(anum)
850       format(//,' EXPERIMENTAL SEMIVARIOGRAMS FOR STRUCTURE ZONE CHARACTERIZED BY INTEGER VALUE OF ',A,'.')
          write(varunit,855)
855       format(/,' Row direction.')
          write(varunit,860)
860       format(/,' Separation          Gamma         Number_of_points')
          do i=1,numsamprow(istruct)
            if(num_row_e(i,istruct).eq.0)then
              gamma_row_e(i,istruct)=0.0d0
            else
              gamma_row_e(i,istruct)=0.5*gamma_row_e(i,istruct)/num_row_e(i,istruct)
            end if
            write(varunit,870) 0.5*(h_row(i,istruct)+h_row(i-1,istruct)),gamma_row_e(i,istruct),  &
            num_row_e(i,istruct)
870         format(1x,1pg14.7,t20,1pg14.7,t40,i7)
          end do
          write(varunit,880)
880       format(/,' Column direction.')
          write(varunit,860)
          do i=1,numsampcol(istruct)
            if(num_col_e(i,istruct).eq.0)then
              gamma_col_e(i,istruct)=0.0d0
            else
              gamma_col_e(i,istruct)=0.5*gamma_col_e(i,istruct)/num_col_e(i,istruct)
            end if
            write(varunit,870) 0.5*(h_col(i,istruct)+h_col(i-1,istruct)),gamma_col_e(i,istruct), &
            num_col_e(i,istruct)
          end do
          if(nangle.gt.0)then
            do iangle=1,nangle
              write(aangle,'(f10.1)') angle(iangle)
              aangle=adjustl(aangle)
              write(varunit,1350) trim(aangle)
1350          format(/,' Angular direction ',a,' degrees.')
              write(varunit,860)
              do i=1,numsampang(iangle,istruct)
                if(num_ang_e(i,iangle,istruct).eq.0)then
                  gamma_ang_e(i,iangle,istruct)=0.0d0
                else
                  gamma_ang_e(i,iangle,istruct)=0.5*gamma_ang_e(i,iangle,istruct)/num_ang_e(i,iangle,istruct)
                end if
                write(varunit,870) 0.5*(h_ang(i,iangle,istruct)+h_ang(i-1,iangle,istruct)),  &
                gamma_ang_e(i,iangle,istruct),num_ang_e(i,iangle,istruct)
              end do
            end do
          end if
        end do
        close(unit=varunit)
        write(6,890) trim(varfile)
890     format('  - file ',a,' written ok.')

! -- If needed, observation names are read from ZONE2VAR1 output file....

        if(dflag.ne.0)then
          read(zonunit,err=9000,end=9050) nobs
          allocate(bobs(nobs),stat=ierr)
          if(ierr.ne.0) go to 9200
          do iobs=1,nobs
            read(zonunit,err=9000,end=9050) bobs(iobs)
          end do
! -- ... and derivatives are written to the matrix file.

          write(6,1010) trim(derfile)
1010      format(/,'  - writing derivatives matrix file ',a,'...')
          allocate(x_row(npar),stat=ierr)
          if(ierr.ne.0) go to 9200
          write(derunit,1020) nobs,npar,2
          do istruct=1,nstructzone
            is=iszone(istruct)
            do i=1,numsamprow(istruct)
              do ipar=1,npar
                if(structzonepar(ipar).ne.is)then
                  x_row(ipar)=0.0d0
                else
                  itemp=num_par_row(i,ipar)
                  if(itemp.eq.0)then
                    x_row(ipar)=0.0d0
                  else
                    x_row(ipar)=deriv_row(i,ipar)/num_row_e(i,istruct)
                    if(at.eq.'l')then
                      x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                    end if
                  end if
                end if
              end do
              write(derunit,1020) (x_row(ipar),ipar=1,npar)
1020          format(8(1x,1pg14.7))
            end do
            do i=1,numsampcol(istruct)
              do ipar=1,npar
                if(structzonepar(ipar).ne.is)then
                  x_row(ipar)=0.0d0
                else
                  itemp=num_par_col(i,ipar)
                  if(itemp.eq.0)then
                    x_row(ipar)=0.0d0
                  else
                    x_row(ipar)=deriv_col(i,ipar)/num_col_e(i,istruct)
                    if(at.eq.'l')then
                      x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                    end if
                  end if
                end if
              end do
              write(derunit,1020) (x_row(ipar),ipar=1,npar)
            end do
            if(nangle.ne.0)then
              do iangle=1,nangle
                do i=1,numsampang(iangle,istruct)
                  do ipar=1,npar
                    if(structzonepar(ipar).ne.is)then
                      x_row(ipar)=0.0d0
                    else
                      itemp=num_par_ang(i,iangle,ipar)
                      if(itemp.eq.0)then
                        x_row(ipar)=0.0d0
                      else
                        x_row(ipar)=deriv_ang(i,iangle,ipar)/num_ang_e(i,iangle,istruct)
                        if(at.eq.'l')then
                          x_row(ipar)=x_row(ipar)/(pvalkeep(ipar)*ln10)
                        end if
                      end if
                    end if
                  end do
                  write(derunit,1020) (x_row(ipar),ipar=1,npar)
                end do
              end do
            end if
          end do
          write(derunit,1030)
1030      format('* row names')
          do iobs=1,nobs
            write(derunit,1035) trim(bobs(iobs))
1035        format(1x,a)
          end do
          write(derunit,1040)
1040      format('* column names')
          do ipar=1,npar
            write(derunit,1035) trim(apar(ipar))
          end do
          close(unit=derunit)
          write(6,1045) trim(derfile)
1045      format('  - file ',a,' written ok.')
       end if

       go to 9900

9000    write(amessage,9010) trim(zonfile)
9010    format(' Error reading ZONE2VAR1-generated data file ',a,'.')
        go to 9890
9050    write(amessage,9060) trim(zonfile)
9060    format(' Unexpected end encountered to ZONE2VAR1-generated data file ',a,'.')
        go to 9890

9100    call num2char(iline,aline)
        write(amessage,9110) trim(aline),trim(ircfile)
9110    format(' Error encountered in reading line ',a,' of integer-real ',  &
        'correspondence file ',a,'.')
        go to 9890
9200    write(amessage,9210)
9210    format(' Cannot allocate sufficient memory to continue execution.')
        go to 9890


9890    call write_message(leadspace='yes')

9900    call close_files

        deallocate(layer,ival,iszone,structzonepar,stat=ierr)
        deallocate(numsampcol,numsamprow,stat=ierr)
        deallocate(numsampang,stat=ierr)
        deallocate(num_col_e,num_row_e,stat=ierr)
        deallocate(num_ang_e,stat=ierr)
        deallocate(num_par_row,num_par_col,stat=ierr)
        deallocate(num_par_ang,stat=ierr)

        deallocate(pval,pvalkeep,stat=ierr)
        deallocate(angle,stat=ierr)
        deallocate(h_row,h_col,h_ang,stat=ierr)
        deallocate(gamma_row_e,gamma_col_e,stat=ierr)
        deallocate(gamma_ang_e,stat=ierr)
        deallocate(deriv_row,deriv_col,stat=ierr)
        deallocate(deriv_ang,stat=ierr)
        deallocate(x_row,stat=ierr)

        deallocate(apar,stat=ierr)
        deallocate(bobs,stat=ierr)

end program zone2var2


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



