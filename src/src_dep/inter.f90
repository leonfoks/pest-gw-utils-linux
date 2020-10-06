!     Last change:  J    14 Jun 2002   10:07 pm
MODULE INTER

! -- Contains interface blocks for all subprograms.


!******************************************************************************
! generic subprograms
!******************************************************************************

interface char2num

	subroutine a2i(ifail,string,num)
	  integer, intent(out)          :: ifail
	  character (len=*), intent(in) :: string
	  integer, intent(out)          :: num
	end subroutine a2i
	subroutine a2r(ifail,string,num)
	  integer, intent(out)          :: ifail
	  character (len=*), intent(in) :: string
	  real, intent(out)             :: num
	end subroutine a2r
	subroutine a2d(ifail,string,num)
	  integer, intent(out)          :: ifail
	  character (len=*), intent(in) :: string
	  double precision, intent(out) :: num
	end subroutine a2d

end interface


interface num2char

	subroutine i2a(value,string,nchar)
	  integer, intent(in)           :: value
	  character (len=*), intent(out):: string
	  integer, intent(in), optional :: nchar
	end subroutine i2a
	subroutine r2a(value,string,nchar)
	  real, intent(in)              :: value
	  character (len=*), intent(out):: string
	  integer, intent(in), optional :: nchar
	end subroutine r2a
	subroutine d2a(value,string,nchar)
	  double precision, intent(in)  :: value
	  character (len=*), intent(out):: string
	  integer, intent(in), optional :: nchar
	end subroutine d2a

end interface


interface pos_test

	integer function pos_i_test(value,string)
	  integer, intent(in)           :: value
	  character (len=*), intent(in) :: string
	end function pos_i_test
	integer function pos_r_test(value,string)
	  real, intent(in)              :: value
	  character (len=*), intent(in) :: string
	end function pos_r_test
	integer function pos_d_test(value,string)
	  double precision, intent(in)  :: value
	  character (len=*), intent(in) :: string
	end function pos_d_test

end interface


interface nneg_test

	integer function nneg_i_test(value,string)
	  integer, intent(in)           :: value
	  character (len=*), intent(in) :: string
	end function nneg_i_test
	integer function nneg_r_test(value,string)
	  real, intent(in)              :: value
	  character (len=*), intent(in) :: string
	end function nneg_r_test
	integer function nneg_d_test(value,string)
	  double precision, intent(in)  :: value
	  character (len=*), intent(in) :: string
	end function nneg_d_test

end interface


interface key_read

	integer function int_key_read(value)
	  integer,intent(out)   :: value
	end function int_key_read
	integer function real_key_read(value)
	  real,intent(out)      :: value
	end function real_key_read
	integer function double_key_read(value)
	  double precision,intent(out)  :: value
	end function double_key_read

end interface



interface equals

       logical function equals_int(r1,r2)
          integer, intent(in)    :: r1
          integer, intent(in)    :: r2
       end function equals_int

       logical function equals_real(r1,r2)
          real, intent(in)      :: r1
          real, intent(in)      :: r2
       end function equals_real

       logical function equals_dbl(r1,r2)
          real (kind (1.0d0)), intent(in)      :: r1
          real (kind (1.0d0)), intent(in)      :: r2
       end function equals_dbl

end interface



!******************************************************************************
! other subprograms
!******************************************************************************

! dxf subprograms ------->

interface
	subroutine write_dxf_header(iunit,xlo,xhi,ylo,yhi)
	  integer, intent(in)             :: iunit
	  double precision, intent(in)    :: xlo,xhi,ylo,yhi
	end subroutine write_dxf_header

	subroutine write_dxf_polyhead(iunit)
	  integer, intent(in)     :: iunit
	end subroutine write_dxf_polyhead

	subroutine write_dxf_vertex(iunit,x,y)
	  integer, intent(in)             :: iunit
	  double precision, intent(in)    :: x,y
	end subroutine write_dxf_vertex

	subroutine write_dxf_polyfin(iunit)
	  integer, intent(in)     :: iunit
	end subroutine write_dxf_polyfin

	subroutine write_dxf_finish(iunit)
	  integer,intent(in)      :: iunit
	end subroutine write_dxf_finish

end interface

! utility subprograms ------->

interface

	subroutine casetrans(string,hi_or_lo)
	  character (len=*), intent(inout)        :: string
	  character (len=*), intent(in)           :: hi_or_lo
	end subroutine casetrans

	subroutine sub_error(subname)
	  character (len=*)               ::subname
	end subroutine sub_error

	integer function nextunit()
	end function nextunit

	subroutine close_files
	end subroutine close_files

	subroutine open_input_file(ifail,aprompt,infile,inunit,file_format,form_prompt, &
          fformat)
          integer, intent(out)                     :: ifail
          character (len=*), intent(in)            :: aprompt
          character (len=*), intent(out)           :: infile
          integer, intent(out)                     :: inunit
          character (len=*), intent(in), optional  :: file_format
          character (len=*), intent(in), optional  :: form_prompt
          character (len=*), intent(out), optional :: fformat
	end subroutine open_input_file

        subroutine open_output_file(ifail,aprompt,outfile,outunit,file_format)
          integer, intent(out)                       :: ifail
          character (len=*)                          :: aprompt,outfile
          integer, intent(out)                       :: outunit
          character (len=*), intent(in), optional    :: file_format
        end subroutine open_output_file

	subroutine wrtsig(ifail,val,word,nw,precis,tval,nopnt)
	  integer               :: ifail,nw,precis,nopnt
	  double precision      :: val,tval
	  character (len=*)     :: word
	end subroutine wrtsig

	subroutine readfig(specfile,coordfile,sampfile,pumpfile,pilotfile)
	  character (len=*), intent(out)                :: specfile
	  character (len=*), intent(out), optional      :: coordfile,sampfile,&
							   pumpfile,pilotfile
	end subroutine readfig

	subroutine read_settings(ifail,idate,iheader)
	  integer, intent(out)	:: ifail,idate,iheader
	end subroutine read_settings

	subroutine char_add(astring,achar)
          character (len=*), intent(inout)        :: astring
          character (len=*), intent(in)           :: achar
	end subroutine char_add

	subroutine int2alph(inum,alph,nsig)
	  integer, intent(in)			:: inum
	  character (len=*), intent(out)	:: alph
	  integer, optional, intent(in)		:: nsig
	end subroutine int2alph

end interface

! interpolation subprograms ------->

interface

	subroutine factor(gridspec,east,north,fac1,fac2,fac3,fac4, &
	icellno,jcellno)
	  use defn
	  type(modelgrid), intent(in)           :: gridspec
	  double precision, intent(in)          :: east,north
	  real, intent(out)                     :: fac1,fac2,fac3,fac4
	  integer, intent(out)                  :: icellno,jcellno
	end subroutine factor

	subroutine point_interp(ncol,nrow,thresh,fac1,fac2,fac3,fac4, &
	icellno,jcellno,bhead,rarray,imethod)
	  integer, intent(in)                   :: ncol,nrow
	  real, intent(in)                      :: thresh
	  real, intent(in)                      :: fac1,fac2,fac3,fac4
	  integer, intent(in)                   :: icellno,jcellno
	  real, intent(out)                     :: bhead
	  real, dimension(0:ncol+1,0:nrow+1), intent(in)  :: rarray
	  character (len=*), intent(in), optional         :: imethod
	end subroutine point_interp

end interface

! reading-a-file subprograms ------->

interface

	subroutine linesplit(ifail,num)
	  integer, intent(out)            :: ifail
	  integer, intent(in)		  :: num
	end subroutine linesplit

	integer function char2int(ifail,num)
	  integer, intent(in)             :: num
	  integer, intent(out)            :: ifail
	end function char2int

	real function char2real(ifail,num)
	  integer, intent(in)             :: num
	  integer, intent(out)            :: ifail
	end function char2real

	double precision function char2double(ifail,num)
	  integer, intent(in)             :: num
	  integer, intent(out)            :: ifail
	end function char2double

        subroutine getfile(ifail,cline,filename,ibeg,iend)
          integer, intent(out)            :: ifail
          integer, intent(in)             :: ibeg
          integer, intent(inout)          :: iend
          character (len=*), intent(in)   :: cline
          character (len=*), intent(out)  :: filename
       end subroutine getfile

       subroutine addquote(afile,aqfile)
          character (len=*), intent(in)   :: afile
          character (len=*), intent(out)  :: aqfile
       end subroutine addquote

end interface

!integer array subprograms ------->

interface

	subroutine read_integer_array(ifail,aprompt,array,pm_header,rows, &
	columns,defaultfile)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(inout)        :: aprompt
	  integer, intent(out),dimension(:,:)     :: array
	  character (len=*), intent(in), optional :: pm_header
	  integer, intent(in), optional           :: rows,columns
          character (len=*), optional             :: defaultfile
	end subroutine read_integer_array

	subroutine write_integer_array(ifail,aprompt,array,pm_header,rows, &
	columns)
	  integer, intent(out)			  :: ifail
	  character (len=*), intent(inout)	  :: aprompt
	  integer, intent(in), dimension(:,:)	  :: array
	  character (len=*), intent(in), optional :: pm_header
	  integer, intent(in), optional		  :: rows,columns
	end subroutine write_integer_array

end interface

! real array subprograms ------->

interface

	subroutine read_real_array(ifail,aprompt,array,pm_header,rows,columns)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(inout)        :: aprompt
	  real, intent(out),dimension(:,:)        :: array
	  character (len=*), intent(in), optional :: pm_header
	  integer, intent(in), optional           :: rows,columns
	end subroutine read_real_array

	subroutine write_real_array(ifail,aprompt,array,pm_header,rows,columns,&
	binary_header,atype,ntrans,kstp,kper,pertim,totim,text,ncol,nrow,ilay,istage, &
        realfile,aaformat)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(inout)        :: aprompt
	  real, intent(in),dimension(:,:)         :: array
	  character (len=*), intent(in), optional :: pm_header
	  integer, intent(in), optional           :: rows,columns
          character (len=*), optional             :: binary_header
          character (len=*), optional             :: atype
          integer, optional                       :: ntrans,kstp,kper
          real, optional                          :: pertim,totim
          character (len=16),optional             :: text
          integer, optional                       :: ncol,nrow,ilay
          integer, intent(in), optional           :: istage
          character (len=*), optional             :: realfile
          character (len=*), optional             :: aaformat
        end subroutine write_real_array

	subroutine read_surf_array(ifail,aprompt,array,gridspec)
	  use defn
	  integer, intent(out)			:: ifail
	  character (len=*), intent(in)		:: aprompt
	  real, intent(out), dimension(:,:)	:: array
	  type(modelgrid), intent(in)		:: gridspec
	end subroutine read_surf_array

	subroutine write_surf_array(ifail,aprompt,array,thresh,gridspec)
	  use defn
	  integer, intent(out)			:: ifail
	  character (len=*), intent(in)		:: aprompt
	  real, intent(inout), dimension(:,:)	:: array
	  real, intent(in)			:: thresh
	  type(modelgrid), intent(in)		:: gridspec
	end subroutine write_surf_array

	subroutine write_real_table_file(ifail,outunit,outfile,intarray,realarray)
          use defn
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: outunit
          character (len=*), intent(in)           :: outfile
          integer, dimension(:,:), intent(in)     :: intarray
          real, dimension(:,:), intent(in)        :: realarray
	end subroutine write_real_table_file

end interface

! message subprograms ------->

interface

	subroutine write_initial_message(leadspace,endspace)
	  character (len=*), intent(in), optional :: leadspace,endspace
	end subroutine write_initial_message

	subroutine write_message(increment,iunit,error,leadspace,endspace)
	  integer, intent(in), optional           ::increment,iunit
	  character (len=*), intent(in), optional ::error,leadspace,endspace
	end subroutine write_message

end interface

! grid subprograms ------->

interface

	subroutine spec_open(ifail,gridspec)
	  use defn
	  integer, intent(out)	  :: ifail
	  type (modelgrid)        :: gridspec
	end subroutine spec_open

	subroutine read_spec_dim(ifail,gridspec)
	  use defn
	  type (modelgrid)                :: gridspec
	  integer, intent(out)            :: ifail
	end subroutine read_spec_dim

	subroutine read_spec_data(ifail,gridspec)
	  use defn
	  type (modelgrid)        :: gridspec
	  integer, intent(out)    :: ifail
	end subroutine read_spec_data

	subroutine close_spec_file(gridspec,ok)
	  use defn
	  type (modelgrid)                :: gridspec
	  character (len=*),optional      :: ok
	end subroutine close_spec_file

	subroutine free_grid_mem(gridspec)
	  use defn
	  type (modelgrid)        :: gridspec
	end subroutine free_grid_mem

	subroutine rel_centre_coords(east,north,gridspec)
	  use defn
	  real, intent(out), dimension(:,:)       :: east,north
	  type(modelgrid), intent(in)             :: gridspec
	end subroutine rel_centre_coords

	subroutine grid2earth(east,north,gridspec)
	  use defn
	  real, intent(inout),dimension(:,:)      :: east,north
	  type(modelgrid), intent(in)             :: gridspec
	end subroutine grid2earth

	subroutine transform_to_earth(earth_east,earth_north,grid_east, &
	grid_north,gridspec)
	  use defn
	  real, intent(out)               :: earth_east, earth_north
	  real, intent(in)                :: grid_east, grid_north
	  type(modelgrid), intent(in)     :: gridspec
	end subroutine transform_to_earth

	subroutine rc2cell(icellno,irow,icol,gridspec)
	  use defn
	  integer, intent(out)            :: icellno
	  integer, intent(in)             :: irow,icol
	  type(modelgrid), intent(in)     :: gridspec
	end subroutine rc2cell

	subroutine rc2cell3d(icellno,irow,icol,ilay,gridspec)
	  use defn
	  integer, intent(out)            :: icellno
	  integer, intent(in)             :: irow,icol,ilay
	  type(modelgrid), intent(in)     :: gridspec
	end subroutine rc2cell3d

	subroutine cell2rc(icellno,irow,icol,gridspec)
          use defn
          integer, intent(in)             :: icellno
          integer, intent(out)            :: irow,icol
          type (modelgrid), intent(in)    :: gridspec
	end subroutine cell2rc

	subroutine corner(ind,e,n,icol,irow,gridspec,ecg,ncg)
	  use defn
	  integer, intent(in)                     :: ind
	  real, intent(out)                       :: e,n
	  integer, intent(in)                     :: icol,irow
	  type (modelgrid), intent(in)            :: gridspec
	  real, dimension(:,:), intent(in)        :: ecg,ncg
	end subroutine corner


        subroutine cell_coordinates(gridspec,east,north,icellrow,icellcol,roff,coff)
          use defn
          type(modelgrid), intent(in)             :: gridspec
          double precision, intent(in)            :: east,north
          integer, intent(out)                    :: icellrow,icellcol
          real, intent(out)                       :: roff,coff
        end subroutine cell_coordinates



end interface

! bore data manipulation subprograms ------->

interface

	subroutine read_bore_coord_file(ifail,aprompt)
	  integer, intent(out)            :: ifail
	  character (len=*), intent(in)   :: aprompt
	end subroutine read_bore_coord_file

	subroutine read_bore_list_file(ifail,aprompt,coord_check)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(in)           :: aprompt
	  character (len=*), intent(in), optional :: coord_check
	end subroutine read_bore_list_file

	subroutine read_bore_diff_file(ifail,aprompt,coord_check)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(in)           :: aprompt
	  character (len=*), intent(in), optional :: coord_check
	end subroutine read_bore_diff_file

	subroutine free_bore_mem
	end subroutine free_bore_mem

	subroutine read_rest_of_sample_line(ifail,cols,ndays,nsecs,value, &
	iline,sampfile)
          integer, intent(out)            :: ifail
          integer, intent(in)             :: cols
          integer, intent(out)            :: ndays,nsecs
          double precision, intent(out)   :: value
          integer, intent(in)             :: iline
          character (len=*), intent(in)   :: sampfile
	end subroutine read_rest_of_sample_line

	subroutine read_rest_of_pump_line(ifail,ibore,ndays,nsecs,pumped, &
	iline,pmpfile)
          integer, intent(out)                    :: ifail
          integer, intent(inout)                  :: ibore
          integer, intent(out), dimension (:)     :: ndays,nsecs
          double precision, intent(inout), dimension(:)  :: pumped
          integer, intent(in)                     :: iline
          character (len=*), intent(in)           :: pmpfile
	end subroutine read_rest_of_pump_line

	subroutine time_interp(ifail,nbore,ndays,nsecs,value,intday, &
	intsec,rnear,rconst,valinterp,extrap,direction)
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: nbore
          integer, intent(in), dimension(nbore)   :: ndays,nsecs
          double precision, intent(in), dimension(nbore)   :: value
          integer, intent(in)                     :: intday,intsec
	  real, intent(in)			  :: rnear,rconst
          double precision, intent(out)           :: valinterp
	  character (len=*), intent(in),optional  :: extrap
	  character (len=*), intent(in),optional  :: direction
	end subroutine time_interp

	subroutine get_num_ids(ifail,iunit,afile,numid,maxsamp,ignore_x)
	  integer, intent(out)                    :: ifail
          integer, intent(in)                     :: iunit
          character (len=*), intent(in)           :: afile
          integer, intent(out)                    :: numid,maxsamp
	  character (len=*), intent(in), optional :: ignore_x
	end subroutine get_num_ids

	subroutine get_ids_and_interval(ifail,iunit,afile,nid,aid,ndays1, &
                                nsecs1,ndays2,nsecs2, ignore_x)
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: iunit
          character (len=*), intent(in)           :: afile
          integer, intent(in)                     :: nid
          character (len=*), intent(out)          :: aid(nid)
          integer, intent(out)                    :: ndays1(nid),nsecs1(nid), &
                                                     ndays2(nid),nsecs2(nid)
	  character (len=*), intent(in), optional :: ignore_x
	end subroutine get_ids_and_interval

end interface

! geostat subprograms ------->

interface

	subroutine read_pilot_points_file(ifail,aprompt,accept_blank)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(in)           :: aprompt
          character (len=*), intent(in), optional :: accept_blank
	end subroutine read_pilot_points_file

	subroutine read_dual_pilot_points_file(ifail,aprompt,accept_blank)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(in)           :: aprompt
          character (len=*), intent(in), optional :: accept_blank
	end subroutine read_dual_pilot_points_file

	subroutine read_3d_pilot_points_file(ifail,aprompt,accept_blank)
	  integer, intent(out)                    :: ifail
	  character (len=*), intent(in)           :: aprompt
          character (len=*), intent(in), optional :: accept_blank
	end subroutine read_3d_pilot_points_file

	subroutine free_point_mem
	end subroutine free_point_mem

        subroutine read_structure_file_dim(ifail,structunit,numstruct,numvario,structfile)
          integer, intent(out)            :: ifail
          integer, intent(in)             :: structunit
          integer, intent(out)            :: numstruct,numvario
          character (len=*), intent(in)   :: structfile
        end subroutine read_structure_file_dim

        subroutine read_rest_of_structure_file(ifail,structunit,numstruct,numvario, &
          structfile,structure,vario)
          use defn
          integer, intent(out)            :: ifail
          integer, intent(in)             :: structunit,numstruct,numvario
          character (len=*), intent(in)   :: structfile
          type (geostructure), intent(out):: structure(numstruct)
          type (variogram), intent(out)   :: vario(numvario)
        end subroutine read_rest_of_structure_file


        subroutine read_parameter_replacement_file(ifail,aprompt,ncol,nrow)
	  integer, intent(out)            :: ifail
	  character (len=*), intent(in)   :: aprompt
          integer, intent(in)             :: ncol,nrow
        end subroutine read_parameter_replacement_file


        subroutine free_replace_mem()
        end subroutine free_replace_mem

end interface


! RSM subprograms ------->

interface

        subroutine read_gms_2d_mesh_file(ifail,aprompt,centroid)
          integer, intent(out)            :: ifail
          character (len=*), intent(in)   :: aprompt
          character (len=*), intent(in)   :: centroid
        end subroutine read_gms_2d_mesh_file

        subroutine read_gms_integer_mesh_data_file(ifail,aprompt,intdata)
          integer, intent(out)             :: ifail
          character (len=*), intent(inout) :: aprompt
          integer, intent(out)             :: intdata(:)
        end  subroutine read_gms_integer_mesh_data_file

        subroutine read_gms_real_mesh_data_file(ifail,aprompt,realdata)
          integer, intent(out)             :: ifail
          character (len=*), intent(inout) :: aprompt
          real, intent(out)                :: realdata(:)
        end  subroutine read_gms_real_mesh_data_file

end interface


! FEFLOW subprograms ------->

interface

        subroutine read_feflow_elem_prop_file(ifail,aprompt)
          integer, intent(out)            :: ifail
          character (len=*), intent(in)   :: aprompt
        end subroutine read_feflow_elem_prop_file

end interface



! date manipulation subprograms ------->

interface

	subroutine char2date(ifail,adate,dd,mm,yy)
          integer, intent(out)    	:: ifail
          character (len=*), intent(in) :: adate
          integer, intent(out) 		:: dd,mm,yy
	end subroutine char2date

	subroutine datestring(dd,mm,yy,hhh,mmm,sss,time,at,adate,atime)
	  integer, intent(in)             :: dd,mm,yy,hhh,mmm,sss
	  real, intent(in)                :: time
	  character (len=1), intent(in)   :: at
	  character (len=*), intent(out)  :: adate, atime
	end subroutine datestring

	logical function leap(year)
          integer, intent(in)     :: year
	end function leap

	integer function numdays(dr,mr,yr,d,m,y)
          integer, intent(in)     :: dr,mr,yr,d,m,y
	end function numdays

	integer function numsecs(h1,m1,s1,h2,m2,s2)
	  integer, intent(in)     :: h1,m1,s1,h2,m2,s2
	end function numsecs

	subroutine char2time(ifail,adate,hh,mm,ss)
          integer, intent(out)    	:: ifail
          character (len=*), intent(in) :: adate
          integer, intent(out) 		:: hh,mm,ss
	end subroutine char2time

	subroutine time2char(ifail,hh,mm,ss,atime)
	  integer, intent(out)            :: ifail
	  integer, intent(in)             :: hh,mm,ss
	  character (len=*), intent(out)  :: atime
	end subroutine time2char

	subroutine elapsdate(eltime,dayfactor,day1,mon1,year1,hour1,min1,sec1,&
	  day2,mon2,year2,hour2,min2,sec2)
	  real, intent(in)		:: eltime,dayfactor
	  integer, intent(in)		:: day1,mon1,year1,hour1,min1,sec1
	  integer, intent(out)		:: day2,mon2,year2,hour2,min2,sec2
	end subroutine elapsdate

	subroutine newdate(ndays,day1,mon1,year1,day2,mon2,year2)
	  integer, intent(in)		:: ndays,day1,mon1,year1
	  integer, intent(out)		:: day2,mon2,year2
	end subroutine newdate

	subroutine sectime(nsecs,sec,min,hour)
	  integer, intent(in)   :: nsecs
	  integer, intent(out)  :: sec,min,hour
	end subroutine sectime

end interface


! mapinfo interface programs ------->


interface

	subroutine write_mif_file(ifail,outunit,outfile,intarray,east,north, &
		        gridspec,izone,type)

          use defn
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: outunit
          character (len=*)                       :: outfile
          integer, dimension(:,:), intent(in)     :: intarray
          real, dimension(:,:), intent(in)        :: east,north
          type (modelgrid), intent(in)            :: gridspec
          integer, intent(in)                     :: izone
          character (len=1)                       :: type
	end subroutine write_mif_file


	subroutine write_real_mid_file(ifail,outunit,outfile,intarray,realarray)
          use defn
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: outunit
          character (len=*), intent(in)           :: outfile
          integer, dimension(:,:), intent(in)     :: intarray
          real, dimension(:,:), intent(in)        :: realarray
	end subroutine write_real_mid_file


	subroutine write_integer_mid_file(ifail,outunit,outfile,intarraywin, &
	  intarraydat)
          use defn
          integer, intent(out)                    :: ifail
          integer, intent(in)                     :: outunit
          character (len=*), intent(in)           :: outfile
          integer, dimension(:,:), intent(in)     :: intarraywin,intarraydat
	end subroutine write_integer_mid_file


	subroutine read_int_tab_file(ifail,tabunit,tabfile,intarray,gridspec)
	  use defn
	  integer, intent(out)			  :: ifail
	  integer, intent(in)			  :: tabunit
	  character (len=*), intent(in)		  :: tabfile
	  integer, dimension(:,:), intent(inout)  :: intarray
	  type(modelgrid), intent(in)		  :: gridspec
	end subroutine read_int_tab_file


	subroutine read_real_tab_file(ifail,tabunit,tabfile,realarray,gridspec)
	  use defn
	  integer, intent(out)			  :: ifail
	  integer, intent(in)			  :: tabunit
	  character (len=*), intent(in)		  :: tabfile
	  real, dimension(:,:), intent(inout)     :: realarray
	  type(modelgrid), intent(in)		  :: gridspec
	end subroutine read_real_tab_file

	subroutine read_real_tab_file_int(ifail,mode,tabunit,tabfile,realarray, &
                                          intarray,gridspec,insunit,insfile,aprefix)
	  use defn
	  integer, intent(out)			  :: ifail
          integer, intent(in)                     :: mode
	  integer, intent(in)			  :: tabunit
	  character (len=*), intent(in)		  :: tabfile
	  real, dimension(:,:), intent(inout)     :: realarray
          integer, dimension(:,:), intent(inout)  :: intarray
	  type(modelgrid), intent(in)		  :: gridspec
          integer, intent(in)                     :: insunit
          character (len=*), intent(in)           :: insfile
          character (len=*), intent(in)           :: aprefix
	end subroutine read_real_tab_file_int


end interface


END MODULE INTER
 
