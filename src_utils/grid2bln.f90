!     Last change:  JD   16 Dec 2000    4:30 pm
program grid2bln

! -- Program grid2bln produces a SURFER blanking file of the
!    finite-difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer					:: ifail,ierr,ncol,nrow, &
						   outunit,idate,iheader
	integer, allocatable, dimension(:,:)	:: intarray
	real, allocatable, dimension(:,:)	:: east,north
	character (len=80)			:: aprompt,outfile
	character (len=1)			:: aa
	character (len=20)			:: atemp


	write(amessage,5)
5	format(' Program GRID2BLN produces a SURFER blanking file of the ',&
        'finite-difference grid.')
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
10	call spec_open(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) go to 9900
	call read_spec_dim(ifail,gridspec)
	if(ifail.ne.0) go to 9900

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(intarray(ncol,nrow),east(ncol,nrow),north(ncol,nrow),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50	  format(' Cannot allocate sufficient memory to run GRID2BLN')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

60	continue
	write(6,*)
	aprompt=' Enter name of window integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
        rows=gridspec%nrow,columns=gridspec%ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(intarray,east,north,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if

80	write(6,*)
	aprompt=' Enter name for SURFER blanking output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  go to 60
	end if
	if(ifail.ne.0) go to 9900

	write(6,100)
100	format(/,' Working.....',/)

	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)

	call gridlines_bln(outunit,ncol,nrow,&
	gridspec%east_corner,gridspec%north_corner,east,north,&
 	intarray,gridspec)

	write(amessage,150) trim(outfile)
150	format('  - SURFER blanking file ',a,' written ok.')
	call write_message(endspace='yes')

9900	call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,east,north,stat=ierr)

end program grid2bln



      SUBROUTINE GRIDLINES_bln(iunit,NCOL,NROW,E0,N0,ECG,NCG,IARRAY,gridspec)

! -- Subroutine gridlines_bln writes grid lines bounding active cells
!    in SURFER "blanking" format.

! -- Arguments are as follows:-
!       iunit:       unit number for output
!       ncol, nrow:  number of rows and columns in finite diffence grid
!       e0, n0:      easting and northing of top left corner of finite-
!                    difference grid
!       ecg, ncg:    x and y coordinates of grid cell centres expressed in a
!                    coordinate system where the x-direction is oriented
!                    in an easterly direction and the origin is at the top
!                    left corner of the finite-difference grid.
!       iarray:      integer "window" array
!       gridspec:    defined variable holding grid specifications

! -- Revision history:-
!       June, 1993:  version 1 of subroutine gridlines.
!       September, 1995: modified for inclusion in Groundwater Data Utilities.

	use defn
	use inter

	type (modelgrid) gridspec

      INTEGER*4 NCOL,NROW,IROW,ILAST,ICOL,iunit
      INTEGER*4 IARRAY(NCOL,NROW)
      REAL*4 E1,N1,E2,N2
      REAL*4 ECG(NCOL,NROW),NCG(NCOL,NROW)
      REAL*8 E0,N0

      DO 100 IROW=1,NROW
        ILAST=0
        DO 80 ICOL=1,NCOL
          IF(IARRAY(ICOL,IROW).EQ.0)THEN
            IF(IROW.EQ.1) GO TO 50
            IF(IARRAY(ICOL,IROW-1).EQ.0) GO TO 50
          END IF
          IF((ILAST.EQ.1).AND.(ICOL.NE.NCOL)) GO TO 80
          IF(ICOL.EQ.NCOL) THEN
            IF(ILAST.EQ.0)THEN
              CALL CORNER(1,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
              WRITE(iunit,30) 2,0
              CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
            END IF
            CALL CORNER(2,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
            CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
            ILAST=0
            GO TO 80
          END IF
          ILAST=1
          CALL CORNER(1,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
          WRITE(iunit,30) 2,0
   30     FORMAT(I5,1X,I5)
          CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
          GO TO 80
   50     IF(ILAST.EQ.0) GO TO 80
          ILAST=0
          CALL CORNER(1,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
          CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
   80   CONTINUE
  100 CONTINUE
      ILAST=0
      IROW=NROW
      DO 180 ICOL=1,NCOL
        IF(IARRAY(ICOL,IROW).NE.0) THEN
          IF((ILAST.EQ.1).AND.(ICOL.NE.NCOL)) GO TO 180
          IF(ICOL.EQ.NCOL) THEN
            IF(ILAST.EQ.0)THEN
              CALL CORNER(3,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
              WRITE(iunit,30) 2,0
              CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
            END IF
            CALL CORNER(4,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
            CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
            ILAST=0
            GO TO 180
          END IF
          ILAST=1
          CALL CORNER(3,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
          WRITE(iunit,30) 2,0
          CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
        ELSE
          IF(ILAST.EQ.0) GO TO 180
          ILAST=0
          CALL CORNER(3,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
          CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
        END IF
  180 CONTINUE

      DO 300 ICOL=1,NCOL
        ILAST=0
        DO 280 IROW=1,NROW
          IF(IARRAY(ICOL,IROW).EQ.0)THEN
            IF(ICOL.EQ.1) GO TO 350
            IF(IARRAY(ICOL-1,IROW).EQ.0) GO TO 350
          END IF
          IF((ILAST.EQ.1).AND.(IROW.NE.NROW)) GO TO 280
          IF(IROW.EQ.NROW)THEN
            IF(ILAST.EQ.0)THEN
              CALL CORNER(1,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
              WRITE(iunit,30)2,0
              CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
            END IF
            CALL CORNER(3,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
            CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
            ILAST=0
            GO TO 280
          END IF
          ILAST=1
          CALL CORNER(1,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
          WRITE(iunit,30) 2,0
          CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
          GO TO 280
  350     IF(ILAST.EQ.0) GO TO 280
          ILAST=0
          CALL CORNER(1,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
          CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
  280   CONTINUE
  300 CONTINUE
      ILAST=0
      ICOL=NCOL
      DO 380 IROW=1,NROW
        IF(IARRAY(ICOL,IROW).NE.0)THEN
          IF((ILAST.EQ.1).AND.(IROW.NE.NROW)) GO TO 380
          IF(IROW.EQ.NROW)THEN
            IF(ILAST.EQ.0)THEN
              CALL CORNER(2,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
              WRITE(iunit,30) 2,0
              CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
            END IF
            CALL CORNER(4,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
            CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
            ILAST=0
            GO TO 380
          END IF
          ILAST=1
          CALL CORNER(2,E1,N1,ICOL,IROW,gridspec,ECG,NCG)
          WRITE(iunit,30) 2,0
          CALL OWRITE(iunit,DBLE(E1)+E0,DBLE(N1)+N0)
        ELSE
          IF(ILAST.EQ.0) GO TO 380
          ILAST=0
          CALL CORNER(2,E2,N2,ICOL,IROW,gridspec,ECG,NCG)
          CALL OWRITE(iunit,DBLE(E2)+E0,DBLE(N2)+N0)
        END IF
  380 CONTINUE

      RETURN
      END


      SUBROUTINE OWRITE(iunit,A,B)

! -- Subroutine owrite writes 2 numbers in a format that depends on their
!    magnitudes.

! -- Arguments are as follows:-
!     iunit: the unit number of the output file
!     a,b:   the two numbers to be written to the output file

! -- Revision history:-
!       June, 1993:  version 1.
!       September, 1995: modified for inclusion in Groundwater Data Utilities.

      INTEGER*4 NB1,NB2,NE1,NE2,iunit
      REAL*8 A,B
      CHARACTER*25 COUT1,COUT2

      IF((ABS(A).GE.10.0).AND.(ABS(A).LT.1.0E11))THEN
        WRITE(COUT1,'(F18.5)')A
      ELSE
        WRITE(COUT1,'(1P,E15.8)')A
      END IF
      IF((ABS(B).GE.10.0).AND.(ABS(B).LT.1.0E11))THEN
        WRITE(COUT2,'(F18.5)') B
      ELSE
        WRITE(COUT2,'(1P,E15.8)')B
      END IF
      CALL BEGEND(COUT1,NB1,NE1)
      CALL BEGEND(COUT2,NB2,NE2)
      WRITE(iunit,10) COUT1(NB1:NE1), COUT2(NB2:NE2)
   10 FORMAT(1X,A,1X,A)

      RETURN
      END



      SUBROUTINE BEGEND(COUT,NB,NE)

      INTEGER*4 NB,NE,I
      CHARACTER*25 COUT

      DO 10 I=1,25
      IF(COUT(I:I).NE.' ') GO TO 20
   10 CONTINUE
   20 NB=I
      DO 30 I=NB+1,25
      IF(COUT(I:I).EQ.' ') GO TO 40
   30 CONTINUE
      I=I+1
   40 NE=I-1

      RETURN
      END


 
