!     Last change:  JD   23 Jan 2001    0:37 am
MODULE DEFN

! -- File DEFN.F90 contains definitions for defined types and global variables.


!****************************************************************************
! defined types
!****************************************************************************

	type modelgrid
	  integer                         :: nrow,ncol
	  double precision                :: east_corner,north_corner,rotation
	  real                            :: cosang,sinang
	  real, dimension(:), pointer     :: delr,delc
	  integer                         :: specunit,specline
	  character (len=200)             :: specfile
	end type modelgrid


        integer, parameter                :: MAX_STRUCT_VARIO=5
        type geostructure
          character (len=10)              :: structname
          integer                         :: numcount_3d
          integer                         :: numvariogram,transform
          real                            :: nugget,mean,maxpowercov
          real                            :: variogram_contrib(MAX_STRUCT_VARIO)
          character (len=10)              :: variogram_name(MAX_STRUCT_VARIO)
        end type geostructure

        type variogram
          logical              :: three_d
          character (len=10)   :: varname
          integer              :: vartype
          real                 :: angle,a,anis
          real                 :: ang1,ang2,ang3
          real                 :: a_hmax,a_hmin,a_vert
        end type variogram


        integer, parameter                :: MAX_OUT_PARAM=100, &
                                             MAX_IN_PARAM=20
        type mf2k_param_replace
          character (len=10)              :: replacetype
          integer                         :: numout,numin
          character (len=10)              :: outparamname(MAX_OUT_PARAM)
          integer                         :: outparamlay1(MAX_OUT_PARAM), &
                                             outparamlay2(MAX_OUT_PARAM)
          integer                         :: newparamlay1(MAX_IN_PARAM), &
                                             newparamlay2(MAX_IN_PARAM)
          character (len=3)               :: newparamprefix
          integer                         :: transformtype
          real                            :: minival,maxival,minpval,maxpval
          character (len=120)             :: minifile,maxifile
          real, dimension(:,:), pointer   :: miniarray,maxiarray
        end type mf2k_param_replace


!****************************************************************************
!global variables
!****************************************************************************

!variables for reading a file ------->

	integer, parameter              	:: NUM_WORD_DIM=100
	integer, dimension(NUM_WORD_DIM)        :: left_word,right_word
	character (len=300)             	:: cline


!variables for writing a message ------->

	integer                 :: imessage=0
	character (len=500)     :: amessage= ' '
	character (len=200)     :: initial_message=' '


!escape variables ------->

	integer                 :: escset=0
	character (len=5)       :: eschar = 'E ~e '


!variables in bore data manipulation ------->

	integer                         :: num_bore_coord, num_bore_list
	character (len=200)             :: bore_coord_file, bore_list_file
	integer, dimension(:), pointer			:: bore_coord_layer
	double precision, dimension(:), pointer         :: bore_coord_east, &
							   bore_coord_north
	character (len=10), dimension(:), pointer       :: bore_coord_id, &
                                                           bore_list_id,  &
                                                           bore_diff_id

!variables in pilot point manipulation ------->

	integer                         :: num_pilot_points
	character (len=200)             :: pilot_points_file=' '
	integer, dimension(:), pointer			:: pilot_point_zone
	double precision, dimension(:), pointer         :: pilot_point_east, &
							   pilot_point_north, &
							   pilot_point_elev, &
                                                           pilot_point_val, &
                                                           pilot_point_val1
	character (len=12), dimension(:), pointer       :: pilot_point_id

!Variables used in characterisation of a gms 2D mesh

        integer                                 :: numnode_g
        integer                                 :: numelem_g
        character (len=200)                     :: meshfile_g
        integer, dimension(:), pointer          :: nodenum_g(:)
        double precision, dimension(:), pointer :: eastnode_g(:), northnode_g(:)
        integer, dimension(:), pointer          :: elemnum_g(:),elemindex_g(:)
        integer, dimension(:), pointer          :: elemnode_g(:,:)
        double precision, dimension(:), pointer :: eastelem_g(:), northelem_g(:)

!       Note: integers in elemindex indicate the order in a sorting list of each element number.
!             That is elemnode(elemindex(3)) is the element with third lowest element number.

!Variables used in characterisation of a FEFLOW mesh

        integer                                 :: numnode_f
        integer                                 :: numelem_f
        character (len=200)                     :: epfile_f
        double precision, dimension(:), pointer :: eastelem_f(:), northelem_f(:)
        double precision, dimension(:), pointer :: datelem_f(:)
        integer, dimension(:), pointer          :: zoneelem_f(:)

!variables in MODFLOW 2000 parameter replacement ------->

        integer                     :: numparamreplace, numzoneremove, nummultremove
        character (len=12), dimension(:), pointer        :: removezone, removemult
        character (len=200)                              :: repfile
        type (mf2k_param_replace), dimension(:), pointer :: replace


!variables recording data settings ------->

	integer				:: datespec
	character (len=3)               :: headerspec


END MODULE DEFN
 
