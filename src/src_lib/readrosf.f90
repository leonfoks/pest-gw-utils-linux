subroutine read_rest_of_structure_file(ifail,structunit,numstruct,numvario,structfile, &
           structure,vario)

! -- Subroutine READ_REST_OF_STRUCTURE_FILE reads the data contained within a
!    structure file.

! -- Arguments are as follows:-
!        ifail:       returned as zero unless an error condition is encountered
!        structunit:  unit number from which structure file is read
!        numstruct:   number of geostatistical structure definitions in file
!        numvario:    number of variogram definitions in file
!        structfile:  name of structure file
!        structure:   array of geostatistical structures read from structure file
!        vario:       array of variogram specs read from structure file

        use defn
        use inter
        implicit none

        integer, intent(out)            :: ifail
        integer, intent(in)             :: structunit,numstruct,numvario
        character (len=*), intent(in)   :: structfile
        type (geostructure), intent(out):: structure(numstruct)
        type (variogram), intent(out)   :: vario(numvario)

        integer                         :: istructure,ivariogram,istruct,ivario,i,j, &
                                           iline,ierr,lifail,iv,itemp,k
        integer                         :: icount_2d,icount_3d
        real                            :: rtemp
        character (len=5)               :: aline,anum
        character (len=10)              :: atemp
        character (len=30)              :: atype,atemp1,atemp2

! -- Some variables are initialised.

        ifail=0
        imessage=0
        istructure=0
        ivariogram=0
        istruct=0
        ivario=0

        do i=1,numstruct
          structure(i)%numvariogram=-9999
          structure(i)%transform=-9999
          structure(i)%nugget=0.0
          structure(i)%mean=-1.1e35
          structure(i)%maxpowercov=10000.0
          structure(i)%numcount_3d=0
          do j=1,MAX_STRUCT_VARIO
            structure(i)%variogram_name(j)=' '
          end do
        end do
        do i=1,numvario
          vario(i)%vartype=-9999
          vario(i)%angle=-1.1e35
          vario(i)%a=-1.1e35
          vario(i)%anis=-1.1e35
          vario(i)%ang1=-1.1e35
          vario(i)%ang2=-1.1e35
          vario(i)%ang3=-1.1e35
          vario(i)%a_hmax=-1.1e35
          vario(i)%a_hmin=-1.1e35
          vario(i)%a_vert=-1.1e35
        end do

	write(initial_message,5) trim(structfile)
5       format(' Errors in structure file ',a,' ----->')

! -- The structure file is now read line by line.

        iline=0
10      iline=iline+1
        call num2char(iline,aline)
        read(structunit,'(a)',err=9000,end=1000) cline
        if(cline.eq.' ') go to 10
        if(cline(1:1).eq.'#') go to 10
        call linesplit(lifail,2)
        if(lifail.ne.0)then
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          write(amessage,11) trim(aline)
11        format('   Two entries must appear at line ',a)
          call write_message(increment=1)
          go to 9990
        end if
        atemp1=cline(left_word(1):right_word(1))
        atemp2=cline(left_word(2):right_word(2))
        call casetrans(atemp1,'hi')
        if(atemp1.eq.'STRUCTURE')then
          istructure=1
          istruct=istruct+1
          iv=0
          if(len_trim(atemp2).gt.10)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,20) trim(aline)
20          format('   Structure name must be 10 characters or less at line ',a)
            call write_message(increment=1)
          end if
          call casetrans(atemp2,'lo')
          structure(istruct)%structname=atemp2(1:10)
        else if(atemp1.eq.'NUMVARIOGRAM')then
          atype='NUMVARIOGRAM'
          if(istructure.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
30          format('   Unexpected "',a,'" string at line ',a)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,itemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
40          format('   Cannot read value for ',a,' at line ',a)
            call write_message(increment=1)
            go to 9990
          else
            if(itemp.gt.MAX_STRUCT_VARIO)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              call num2char(MAX_STRUCT_VARIO,anum)
              write(amessage,50) trim(anum),trim(aline)
50            format('   NUMVARIOGRAM must not be greater than ',a,' at line ',a)
              call write_message(increment=1)
              go to 9990
            else if(itemp.le.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,60) trim(atype),trim(aline)
60            format('   ',a,' must be greater than zero at line ',a)
              call write_message(increment=1)
              go to 9990
            else
              structure(istruct)%numvariogram=itemp
            end if
          end if
        else if(atemp1.eq.'NUGGET')then
          atype='NUGGET'
          if(istructure.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            structure(istruct)%nugget=rtemp
          end if
        else if(atemp1.eq.'TRANSFORM')then
          atype='TRANSFORM'
          if(istructure.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call casetrans(atemp2,'lo')
          if(atemp2.eq.'log')then
            structure(istruct)%transform=1
          else if(atemp2.eq.'none')then
            structure(istruct)%transform=0
          else
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,70) trim(aline)
70          format('   TRANSFORM must be "log" or "none" at line ',a,'.')
            call write_message(increment=1)
          end if
        else if(atemp1.eq.'MAXPOWERVAR')then
          atype='MAXPOWERVAR'
          if(istructure.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            if(rtemp.le.0.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,170) trim(atype),trim(aline)
              call write_message(increment=1)
            else
              structure(istruct)%maxpowercov=rtemp
            end if
          end if
        else if(atemp1.eq.'MEAN')then
          atype='MEAN'
          if(istructure.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            structure(istruct)%mean=rtemp
          end if
        else if(atemp1.eq.'VARIOGRAM')then
          atype='VARIOGRAM'
          if(istructure.eq.1)then
            iv=iv+1
            if(structure(istruct)%numvariogram.eq.-9999)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,115) trim(aline)
115           format('   A variogram reference has been read for a structure ', &
              'before a value for NUMVARIOGRAM has been supplied for that structure ', &
              'at line ',a)
              call write_message(increment=1)
              go to 9990
            end if
            if(iv.gt.structure(istruct)%numvariogram)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,120) trim(structure(istruct)%structname)
120           format('   There are more variograms in structure "',a, &
              '" than indicated by the value of NUMVARIOGRAM for this structure.')
              call write_message(increment=1)
              go to 9990
            end if
            call linesplit(lifail,3)
            if(lifail.ne.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,130) trim(aline)
130           format('   Three entries are expected on line ',a)
              call write_message(increment=1)
              go to 10
            end if
            atemp1=cline(left_word(2):right_word(2))
            atemp2=cline(left_word(3):right_word(3))
            if(len_trim(atemp1).gt.10)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,140) trim(aline)
140           format('   Variogram name must be 10 characters or less at line ',a)
              call write_message(increment=1)
            end if
            call casetrans(atemp1,'lo')
            structure(istruct)%variogram_name(iv)=atemp1(1:10)
            atype='Variogram contribution'
            call char2num(ifail,atemp2,rtemp)
            if(ifail.ne.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,40) trim(atype),trim(aline)
              call write_message(increment=1)
            else
              if(rtemp.le.0.0)then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,170) trim(atype),trim(aline)
                call write_message(increment=1)
              else
                structure(istruct)%variogram_contrib(iv)=rtemp
              end if
            end if
          else
            ivariogram=1
            ivario=ivario+1
            if(len_trim(atemp2).gt.10)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,140) trim(aline)
              call write_message(increment=1)
            end if
            call casetrans(atemp2,'lo')
            vario(ivario)%varname=atemp2(1:10)
          end if
        else if(atemp1.eq.'VARTYPE')then
          atype='VARTYPE'
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(atype),trim(aline)
150         format('   Unexpected "',a,'" string at line ',a,': is there a ',&
            'missing "END STRUCTURE" string?')
            call write_message(increment=1)
            go to 9990
          end if
          if(ivariogram.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,itemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            if((itemp.lt.1).or.(itemp.gt.4))then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,160) trim(aline)
160           format('   VARTYPE must be 1,2,3 or 4 at line ',a)
              call write_message(increment=1)
            else
              vario(ivario)%vartype=itemp
            end if
          end if
        else if(atemp1.eq.'BEARING')then
          atype='BEARING'
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          if(ivariogram.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            vario(ivario)%angle=rtemp
          end if
        else if(atemp1.eq.'A')then
          atype='A'
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          if(ivariogram.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            if(rtemp.lt.0.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,171) trim(atype),trim(aline)
171           format('   ',a,' must not be less than zero at line ',a)
              call write_message(increment=1)
            else
              vario(ivario)%a=rtemp
            end if
          end if
        else if(atemp1.eq.'ANISOTROPY')then
          atype='ANISOTROPY'
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          if(ivariogram.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            if(rtemp.le.0.0)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,170) trim(atype),trim(aline)
170           format('   ',a,' must be greater than zero at line ',a)
              call write_message(increment=1)
            else
              vario(ivario)%anis=1.0/rtemp
            end if
          end if
        else if((atemp1.eq.'ANG1').or.      &
                (atemp1.eq.'ANG2').or.      &
                (atemp1.eq.'ANG3').or.      &
                (atemp1.eq.'A_HMAX').or.    &
                (atemp1.eq.'A_HMIN').or.    &
                (atemp1.eq.'A_VERT')) then
          if(atemp1.eq.'ANG1')then
                atype='ANG1'
          else if(atemp1.eq.'ANG2')then
                atype='ANG2'
          else if(atemp1.eq.'ANG3')then
                atype='ANG3'
          else if(atemp1.eq.'A_HMAX')then
                atype='A_HMAX'
          else if(atemp1.eq.'A_HMIN')then
                atype='A_HMIN'
          else if(atemp1.eq.'A_VERT')then
                atype='A_VERT'
          end if
          if(istructure.eq.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,150) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 9990
          end if
          if(ivariogram.ne.1)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,30) trim(atype),trim(aline)
            call write_message(increment=1)
            go to 10
          end if
          call char2num(ifail,atemp2,rtemp)
          if(ifail.ne.0)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,40) trim(atype),trim(aline)
            call write_message(increment=1)
          else
            if(atype(1:3).eq.'ANG')then
              if((rtemp.le.-360.0).or.(rtemp.ge.360.0))then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,410) trim(atype),trim(aline)
410             format('   ',a,' must be greater than -360 and less than 360 at line ',a)
                call write_message(increment=1)
              else
                if(atype.eq.'ANG1')then
                  vario(ivario)%ang1=rtemp
                else if(atype.eq.'ANG2')then
                  vario(ivario)%ang2=rtemp
                else if(atype.eq.'ANG3')then
                  vario(ivario)%ang3=rtemp
                end if
              end if
            else
              if(rtemp.le.0.0)then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,411) trim(atype),trim(aline)
411             format('   ',a,' must be greater than zero at line ',a)
                call write_message(increment=1)
              else
                if(atype.eq.'A_HMAX')then
                  vario(ivario)%a_hmax=rtemp
                else if(atype.eq.'A_HMIN')then
                  vario(ivario)%a_hmin=rtemp
                else
                  vario(ivario)%a_vert=rtemp
                end if
              end if
            end if
          end if
        else if(atemp1.eq.'END')then
          if(istructure.eq.1)then
            istructure=0
          else if(ivariogram.eq.1)then
            ivariogram=0
          end if
        else
          atype=cline(left_word(1):right_word(1))
          if(imessage.eq.0) call write_initial_message(leadspace='yes')
          write(amessage,180) trim(atype),trim(aline)
180       format('   Unknown keyword "',a,'" at line ',a)
          call write_message(increment=1)
        end if
        go to 10

! -- Some final checks are made on the dataset contained in the structure file.

1000    continue
        if(imessage.ne.0) go to 9990

        if(numstruct.gt.1)then
          do i=1,numstruct-1
            atemp=structure(i)%structname
            do j=i+1,numstruct
              if(atemp.eq.structure(j)%structname)then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,205)
205             format('   At least two structures have the same name.')
                call write_message(increment=1)
                go to 9990
              end if
            end do
          end do
        end if
        if(numvario.gt.1)then
          do i=1,numvario-1
            atemp=vario(i)%varname
            do j=i+1,numvario
              if(atemp.eq.vario(j)%varname)then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,206)
206             format('   At least two variograms have the same name.')
                call write_message(increment=1)
                go to 9990
              end if
            end do
          end do
        end if

        do i=1,numstruct
          itemp=structure(i)%numvariogram
          if(itemp.eq.-9999)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,200) trim(structure(i)%structname)
200         format('   No variograms supplied for structure "',a,'".')
            call write_message(increment=1)
          else
            do j=1,itemp
              if(structure(i)%variogram_name(j).eq.' ')then
                if(imessage.eq.0) call write_initial_message(leadspace='yes')
                write(amessage,210) trim(structure(i)%structname)
210             format('   Less variograms were cited for structure "',a, &
                '" than expected from the NUMVARIOGRAM variable provided for this ', &
                'structure.')
                call write_message(increment=1)
                go to 9990
              end if
            end do
          end if
250       continue
          if(structure(i)%nugget.lt.-1.0e35)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,255) trim(structure(i)%structname)
255         format('   No NUGGET value supplied for structure "',a,'".')
            call write_message(increment=1)
          end if
          if(structure(i)%transform.eq.-9999)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,256) trim(structure(i)%structname)
256         format('   No TRANSFORM value supplied for structure "',a,'".')
            call write_message(increment=1)
          end if
        end do

        do i=1,numvario
          if(vario(i)%vartype.eq.-9999)then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,260) trim(vario(i)%varname)
260         format('   No value supplied for VARTYPE for variogram "',a,'"')
            call write_message(increment=1)
          end if
          if((vario(i)%ang1.gt.-1.0e35).or.      &
             (vario(i)%ang2.gt.-1.0e35).or.      &
             (vario(i)%ang3.gt.-1.0e35).or.      &
             (vario(i)%a_hmax.gt.-1.0e35).or.    &
             (vario(i)%a_hmin.gt.-1.0e35).or.    &
             (vario(i)%a_vert.gt.-1.0e35))then
            if((vario(i)%ang1.lt.-1.0e35).or.      &
               (vario(i)%ang2.lt.-1.0e35).or.      &
               (vario(i)%ang3.lt.-1.0e35).or.      &
               (vario(i)%a_hmax.lt.-1.0e35).or.    &
               (vario(i)%a_hmin.lt.-1.0e35).or.    &
               (vario(i)%a_vert.lt.-1.0e35))then
               if(imessage.eq.0) call write_initial_message(leadspace='yes')
               write(amessage,440) trim(vario(i)%varname)
440            format('   Variogram "',a,'": if a value is supplied for one of ANG1, ANG2, ANG3, A_HMAX, ',    &
               'A_HMIN or A_VERT (implying that the variogram is a 3D variogram) ',      &
               'then it must be supplied for all of them.')
               call write_message(increment=1)
            end if
            if((vario(i)%angle.gt.-1.0e35).or.     &
               (vario(i)%a.gt.-1.0e35).or.         &
               (vario(i)%anis.gt.-1.0e35))then
               if(imessage.eq.0) call write_initial_message(leadspace='yes')
               write(amessage,450) trim(vario(i)%varname)
450            format('   Variogram "',a,'": if a value is supplied for one of ANG1, ANG2, ANG3, A_HMAX, ',    &
               'A_HMIN or A_VERT (implying that the variogram is a 3D variogram) ',      &
               'then it must not be supplied for any of A, BEARING or ANISOTROPY.')
               call write_message(increment=1)
            end if
          else
            if(vario(i)%angle.lt.-1.0e35)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,270) trim(vario(i)%varname)
270           format('   No value supplied for BEARING for variogram "',a,'"')
              call write_message(increment=1)
            end if
            if(vario(i)%a.lt.-1.0e35)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,280) trim(vario(i)%varname)
280           format('   No value supplied for "A" for variogram "',a,'"')
              call write_message(increment=1)
            end if
            if(vario(i)%anis.lt.-1.0e35)then
              if(imessage.eq.0) call write_initial_message(leadspace='yes')
              write(amessage,290) trim(vario(i)%varname)
290           format('   No value supplied for ANISOTROPY for variogram "',a,'"')
              call write_message(increment=1)
            end if
          end if
          if((vario(i)%ang1.gt.-1.0e35).and.      &
             (vario(i)%ang2.gt.-1.0e35).and.      &
             (vario(i)%ang3.gt.-1.0e35).and.      &
             (vario(i)%a_hmax.gt.-1.0e35).and.    &
             (vario(i)%a_hmin.gt.-1.0e35).and.    &
             (vario(i)%a_vert.gt.-1.0e35))then
               vario(i)%angle=vario(i)%ang1
               vario(i)%a=vario(i)%a_hmax*cos(vario(i)%ang2/180.0d0*3.14159265)
               vario(i)%anis=(vario(i)%a_hmin*cos(vario(i)%ang3/180.0d0*3.14159265))/vario(i)%a_hmax
          end if
          if((vario(i)%vartype.eq.4).and.(vario(i)%a.gt.2.0))then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,294) trim(vario(i)%varname)
294         format('   Power in variogram "',a,'" must be 2 or less.')
            call write_message(increment=1)
          end if
          if((vario(i)%vartype.lt.4).and.(vario(i)%a.eq.0.0))then
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,296) trim(vario(i)%varname)
296         format('   "A" in variogram "',a,'" must be greater than zero.')
            call write_message(increment=1)
          end if
!          if(vario(i)%vartype.eq.4)then
!            atemp=vario(i)%varname
!            do j=1,numstruct
!              if(structure(j)%maxpowercov.gt.-1e35) go to 292
!              do k=1,structure(j)%numvariogram
!                if(atemp.eq.structure(j)%variogram_name(k))then
!                  if(imessage.eq.0) call write_initial_message(leadspace='yes')
!                  write(amessage,291) trim(structure(j)%structname)
!291               format('   Structure "',a,'" cites at least one variogram which ', &
!                  'uses a power model; however no value is provided for MAXPOWERVAR.')
!                  call write_message(increment=1)
!                  go to 292
!                end if
!              end do
!292         continue
!            end do
!          end if
        end do

        if(imessage.eq.0)then
          do i=1,numvario
            if(vario(i)%ang1.gt.-1.0e35)then
              vario(i)%three_d=.true.
            else
              vario(i)%three_d=.false.
            end if
          end do
        end if

        do i=1,numstruct
          icount_2d=0
          icount_3d=0
          do j=1,structure(i)%numvariogram
            atemp1=structure(i)%variogram_name(j)
            do k=1,numvario
              if(atemp1.eq.vario(k)%varname) go to 300
            end do
            if(imessage.eq.0) call write_initial_message(leadspace='yes')
            write(amessage,295) trim(structure(i)%variogram_name(j)), &
                                trim(structure(i)%structname)
295         format('   No variogram data provided for variogram "',a, &
            '" cited in specifications for structure "',a,'".')
            call write_message(increment=1)
300         continue
            if(vario(k)%three_d)then
              icount_3d=icount_3d+1
            end if
          end do
          structure(i)%numcount_3d=icount_3d
        end do

        if(imessage.ne.0) go to 9990
        go to 9999


9000    call num2char(iline,aline)
	write(amessage,9010) trim(aline),trim(structfile)
9010    format(' Error reading line ',a,' of structure file ',a,'.')
	call write_message(leadspace='yes',endspace='yes')
	go to 9990

9890    call write_message()
9990    ifail=1
9999    continue
        imessage=0
        close(unit=structunit,iostat=ierr)
        return

end subroutine read_rest_of_structure_file
