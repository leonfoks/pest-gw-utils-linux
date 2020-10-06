!     Last change:  JD   16 Dec 2000    4:24 pm
program zone2dxf

! -- Program zone2dxf produces a DXF file of the zonation of the finite-
!    difference grid.

	use defn
	use inter

	implicit none

	type (modelgrid) gridspec

	integer,parameter                       :: MAXVERT=5000
	integer                                 :: ifail,ierr,ncol,nrow, &
						   outunit,idate,iheader
	integer, allocatable, dimension(:,:)    :: intarray,jarray
	real, allocatable, dimension(:,:)       :: east,north
	real, allocatable, dimension(:)         :: evert,nvert
	real, dimension(4)                      :: ecorner,ncorner
	double precision, dimension(1)          :: emin,nmin,emax,nmax
	character (len=80)                      :: aprompt,outfile
	character (len=1)                       :: aa
	character (len=20)                      :: atemp


	write(amessage,5)
5       format(' Program ZONE2DXF produces a DXF file of the ',&
	'zonation of the finite-difference grid.')
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

	ncol=gridspec%ncol
	nrow=gridspec%nrow
	allocate(intarray(ncol,nrow),jarray(ncol,nrow),east(ncol,nrow), &
	north(ncol,nrow), evert(MAXVERT), nvert(MAXVERT),stat=ierr)
	if(ierr.ne.0) then
	  write(amessage,50)
50        format(' Cannot allocate sufficient memory to run ZONE2DXF')
	  call write_message(leadspace='yes',endspace='yes')
	  go to 9900
	end if

	call read_spec_data(ifail,gridspec)
	if(ifail.ne.0) go to 9900
	call close_spec_file(gridspec,ok='yes')

60      continue
	write(6,*)
	aprompt=' Enter name of zoned integer array file: '
	call read_integer_array(ifail,aprompt,intarray,pm_header=headerspec, &
	rows=gridspec%nrow, columns=gridspec%ncol)
	if(ifail.ne.0) go to 9900
	if(escset.eq.1) then
	  escset=0
	  write(6,*)
	  deallocate(intarray,jarray,east,north,evert,nvert,stat=ierr)
	  call free_grid_mem(gridspec)
	  go to 10
	end if
	
80      write(6,*)
	aprompt=' Enter name for DXF output file: '
	call open_output_file(ifail,aprompt,outfile,outunit)
	if(escset.eq.1) then
	  escset=0
	  go to 60
	end if
	if(ifail.ne.0) go to 9900

	write(6,100)
100     format(/,' Working.....',/)
	
	call rel_centre_coords(east,north,gridspec)
	call grid2earth(east,north,gridspec)

	ecorner(1)=east(1,1)
	ecorner(2)=east(ncol,1)
	ecorner(3)=east(1,nrow)
	ecorner(4)=east(ncol,nrow)
	ncorner(1)=north(1,1)
	ncorner(2)=north(ncol,1)
	ncorner(3)=north(1,nrow)
	ncorner(4)=north(ncol,nrow)
	emin=minval(ecorner)+gridspec%east_corner
	nmin=minval(ncorner)+gridspec%north_corner
	emax=maxval(ecorner)+gridspec%east_corner
	nmax=maxval(ncorner)+gridspec%north_corner
	call write_dxf_header(outunit,emin(1),emax(1),nmin(1),nmax(1))

	call zonelines_dxf(outunit,ncol,nrow,gridspec%east_corner,&
	gridspec%north_corner,east,north,intarray,jarray,evert,nvert, &
	MAXVERT,gridspec)

	call write_dxf_finish(outunit)

	write(amessage,150) trim(outfile)
150     format('  - DXF file ',a,' written ok.')
	call write_message(endspace='yes')

9900    call close_files
	call free_grid_mem(gridspec)
	deallocate(intarray,jarray,east,north,evert,nvert,stat=ierr)

end program zone2dxf




      SUBROUTINE ZONELINES_dxf(iunit,NCOL,NROW,E0,N0,ECG,NCG,IARRAY,JARRAY,E,N,&
      MAXVERT,gridspec)

! -- Subroutine zonelines_dxf writes integer array zone boundaries in DXF
!    format.

! -- Arguments are as follows:-
!       iunit:      unit number of output file
!       ncol,nrow:  number of columns and rows in finite difference grid
!       e0,n0:      coordinates of top left corner of finite-difference grid
!       ecg,ncg:    cell centre coordinates expressed in a coordinate system
!                   in which the x-axis points east and the origin is situated
!                   at the top left corner of the finite-difference grid.
!       iarray:     zonation-defining integer array
!       jarray:     used in drawing zone boundaries
!       e,n:        contain calculated coordinates for zone boundaries
!       maxvert:    maximum number of zone boundary vertices
!       gridspec:   defined type holding grid specifications

! -- Revision history:-
!       June, 1993: version 1 of subroutine zonelines.
!       September, 1995: modified for inclusion in groundwater data utilities.

	use defn
	use inter

      INTEGER*4 NCOL,NROW,NUMVERT,MAXVERT,ICOL,IROW,JCOL,JROW,NDIR,NZERO,I
	integer iunit
      INTEGER*4 IARRAY(NCOL,NROW),JARRAY(NCOL,NROW)
      REAL*4 ECG(NCOL,NROW),NCG(NCOL,NROW),E(MAXVERT),N(MAXVERT)
      REAL*8 E0,N0

	type (modelgrid) gridspec

      DO 20 IROW=1,NROW
      DO 10 ICOL=1,NCOL
	JARRAY(ICOL,IROW)=0
   10 CONTINUE
   20 CONTINUE

      DO 500 IROW=1,NROW
      DO 499 ICOL=1,NCOL

  100   IF((ICOL.EQ.1).AND.(IROW.EQ.1)) THEN
	  GO TO 150
	ELSE IF(IROW.EQ.1) THEN
	  IF(IARRAY(ICOL-1,IROW).EQ.IARRAY(ICOL,IROW)) GO TO 499
	ELSE IF(ICOL.EQ.1) THEN
	  IF(IARRAY(ICOL,IROW-1).EQ.IARRAY(ICOL,IROW)) GO TO 499
	ELSE
	  IF((IARRAY(ICOL-1,IROW).EQ.IARRAY(ICOL,IROW)).OR. &
	  (IARRAY(ICOL,IROW-1).EQ.IARRAY(ICOL,IROW))) GO TO 499
	END IF
	IF(JARRAY(ICOL,IROW).EQ.1) GO TO 499

  150   NUMVERT=1
	CALL CORNER(1,E(1),N(1),ICOL,IROW,gridspec,ECG,NCG)
	JARRAY(ICOL,IROW)=1

	NDIR=2
	JCOL=ICOL
	JROW=IROW
  200   NUMVERT=NUMVERT+1
	IF(NUMVERT.GT.MAXVERT)THEN
	  write(amessage,210)
210       format(' Internal error in subroutine ZONELINES: increase value ', &
	  'assigned to parameter MAXVERT in ZONE2DXF source code and recompile.')
	  call write_message(endspace='yes')
	  stop
	END IF
	CALL TRAVEL(NZERO,NCOL,NROW,JCOL,JROW,NDIR,&
	IARRAY,JARRAY,E(NUMVERT),N(NUMVERT),ECG,NCG,gridspec)
	IF(NZERO.EQ.1) GO TO 499
	IF(((JCOL.EQ.ICOL).AND.(JROW.EQ.IROW)).AND.(NDIR.EQ.2)) THEN
	  call write_dxf_polyhead(iunit)
	  DO I=1,NUMVERT
	  CALL write_dxf_vertex(iunit,DBLE(E(I))+E0,DBLE(N(I))+N0)
	  end do
	  call write_dxf_polyfin(iunit)
	ELSE
	  GO TO 200
	END IF
  499 CONTINUE
  500 CONTINUE

      RETURN
      END

 