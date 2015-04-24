!> \file general_tools.f90
!! \brief general tools definition
!! @author Innocent Souopgui
!! \details
!! 1D dynamic array, 2D dynamic array. The 2D Dynamic arrays look like java 2D, it is a 1D array of 1D array.
!<

!> Module defining general tools used in many contexts
!!\todo Turn the definition of multiple integer constants into enumeration.
!!  XDMT should be turned into XMT as there is no confusion
!!  rewrite build_regular_icoord_2D on the model of build_regular_icoord_3D
!!  accounting for the fast varying dimension in fortran
!!
!<
module general_tools
    use general_constant
    use debug_tools
implicit none
  !enumeration for the type of files
  integer, parameter :: &
      CTL_DATA  = 711,&!<control vector data
      BCTL_DATA = 712,&!<background ctl
      ACTL_DATA = 713,&!<analysed ctl
      DCTL_DATA = 714,&!<Delta ctl (increment to add to the background to get the value of the ctl)
      TCTL_DATA = 715,&!<true ctl (increment to add to the background to get the value of the ctl)
      COORD_DATA= 716,& !< coordinate data (for the cordinates of the CTL in the physical space)
      OBS_DATA  = 731,& !< observation data
      AOBS_DATA = 732,& !< observation data derived from the analysis
      TOBS_DATA = 733,&!<true observation
      DMT_DATA  = 741,&!<direct model trajectory
      ADMT_DATA = 742,&!< analysed direct model trajectory
      AMT_DATA  = 742,&!< analysed model trajectory, same thing as ADMT_DATA, prefer AMT_DATA
      BDMT_DATA = 743,&!< background direct model trajectory
      BMT_DATA  = 743,&!< background model trajectory, same thing as BDMT_DATA, prefer BMT_DATA
      TDMT_DATA = 744,&!< true direct model trajectory
      TMT_DATA  = 744,&!< true model trajectory, same thing as TDMT_DATA, prefer TMT_DATA
      SMT_DATA  = 745,&!< simulation model trajectory (temporary file used during minimization)
      DMS_DATA  = 746,&!< direct model state
      AMS_DATA  = 747,&!< analysed model state
      BMS_DATA  = 748,&!< background model state
      TMS_DATA  = 749,&!< true model state
      OGAP_DATA = 750,&
      EP_DATA   = 760,&
      COST_DATA = 770,&
      GRAD_DATA = 780,&
      CTL_SIZE  = 717,&
      CTL_BOUND_DATA = 718,&
      CTL_ALL_DATA   = 719,& !<All data about the control vector including ctl, background and bounds
      PCTL_DATA = 720,& !< preconditionned delta ctl
      PCTLF_DATA = 721,& !< tangent linear of the preconditionned delta ctl
      IMS_DATA  = 791,& !< initial model state
      BIMS_DATA = 792,& !< background initial model state
      AIMS_DATA = 793,& !< analysed initial model state
      TIMS_DATA = 794,& !< analysed initial model state
      BC_DATA   = 801,& !< boundary condition
      ENS_DATA  = 820,& !< ensemble state
      ENS_MEAN_DATA = 821,& !< ensemble mean
      ENS_VAR_DATA  = 822,& !> ensemble variance
      ENS_ME_DATA      = 823,& !> ensemble model error
      ME_MEAN_DATA = 824,& !> ensemble model error mean
      ME_VAR_DATA  = 825,& !> ensemble model error var
      OTHER_DATA= 9001   !< Other data, for non defined purposes


  !enumeration for the status or usage of files
  !todo this should be turn into enumeration type
  integer, parameter :: &
      INPUT_FILE  = 381,& !< input file
      OUTPUT_FILE = 382,& !< output file
      RTIME_FILE  = 383,& !< run time file
      BASE_FILE   = 384,& !< base file name (with relative or full path included)
      USER_FILE   = 384 !< more appriopriate name for BASE_FILE

  character(len=1) ::&
      APREF = 'a'& !< single char prefix for analysis
    , BPREF = 'b'& !< single char prefix for background
    , DPREF = 'd'& !< single char prefix for direct
    , FPREF = 'f'& !< single char prefix for forecast
    , SPREF = 's'& !< single char prefix for simulation
    , TPREF = 't' !< single char prefix for truth


  interface vec2mat
    module procedure vec2mat_sp
    module procedure vec2mat_dp
    module procedure vec2mat3d_dp
  end interface vec2mat

  interface mat2vec
    module procedure mat2vec_sp
    module procedure mat2vec_dp
  end interface mat2vec

  !> Interface for printing elapse or running time
  interface printElapseTime
    module procedure printElapseTime_cpu_time !< using cpu_time date
    module procedure printElapseTime_system_clock !< using system clock data
  end interface printElapseTime

  interface make_pFileName
    module procedure make_pFileName_from_scalar
    module procedure make_pFileName_from_fName
  end interface make_pFileName

  interface norm2
    module procedure norm2_1d_dp
    module procedure norm2_2d_dp
  end interface norm2

  interface search_pos_inc
    module procedure search_inc_int
    module procedure search_inc_dp_real
  end interface search_pos_inc

  interface search_pos_dec
    !module procedure search_dec_int
    module procedure search_dec_dp_real
  end interface search_pos_dec

  !> \brief sort routines, selection
  !<
  interface select_sort
    module procedure select_sort_int
    !module procedure select_sort_sp_real
    !module procedure select_sort_dp_real
  end interface

  !> \brief switch routines
  !<
  interface switch
    module procedure switch_int
    !module procedure switch_sp_real
    !module procedure switch_dp_real
  end interface

    !> \brief integer convertion
    !<
    interface int
        module procedure logical2int
    end interface int

    !> interface for saving vectors for plot
    interface write_vector_for_plot
        module procedure write_vector_for_plot_default
        module procedure write_vector_for_plot_1D_icoord
        module procedure write_vector_for_plot_1D_rcoord
        module procedure write_vector_for_plot_ND_rcoord
        module procedure write_vector_for_plot_domain
    end interface write_vector_for_plot

    !> interface for saving data for plot
    interface write_plot_data
        ! 1D
        module procedure write_plot_data_1D_default
        module procedure write_plot_data_1D_icoord
        module procedure write_plot_data_1D_rcoord
        module procedure write_plot_data_1D_ND_rcoord
        module procedure write_plot_data_1D_domain
        !2D
        module procedure write_plot_data_2D_default
        module procedure write_plot_data_2D_rcoord
    end interface write_plot_data

    !> @brief transform a map from vector to multi-D array into
    !!  a map from multi-D array to vector
    !! @details given a 2D array that give for each index of a vector, its indexes
    !!  in a multi-D array, this subroutine constructs the multi-D array that
    !!  gives for each position in the multi-D array its position in a vector array
    !! this is usefull when all the elements of the multi-D array are not usefull
    !! and does not need to be store when vectorizing the array.
    !! 2D example: 2D data array data with masked location
    !! that need to be vectorized and restaured for some operations.
    !! 2-row vector coord that maps the vectorized
    !! data into the 2D array, and 2D loc (same shape as data) that maps
    !! the vectorized data into the 2D array.
    !! In the vectorization operation,
    !! if the unmasked (usefull) element data[i,j] is store in the vectorized
    !! data at index l: then:
    !!  - coord[1,l] = i and coord[2,l] = j
    !!  - loc[i,j] = l
    !<
    interface coord2loc
        !module procedure coord2loc_1D
        module procedure coord2loc_2D
        module procedure coord2loc_3D
    end interface coord2loc

    !> \brief used to sample obvservation from the state vector
    interface subsample
        !> sample from 1D array to 1D array
        module procedure sample_1D_direct_obs
        !> sample from 2D array to 1D array
        module procedure subsample_2D_to_1D
        !> sample from 3D array to 1D array
        module procedure subsample_3D_to_1D
    end interface subsample

    !> \brief used to upsample obvservation to the state vector
    interface upsample
        ! !> upsample from 1D array to 1D array
        !module procedure upsample_1D_direct_obs
        !> upsample from 1D array to 2D array
        module procedure upsample_1D_to_2D
        !> upsample from 1D array to 3D array
        module procedure upsample_1D_to_3D
    end interface upsample

    interface assignment(=)
        !> assign 1D array to drv
        module procedure assign_drv_vector
        !> assign 1D array to div
        module procedure assign_div_vector
        !> assign 1D array to dcv
        module procedure assign_dsv_vector
        !> assign 1D array to d2Da
        module procedure assign_2ddra_vector
        !> assign 2D array to d2Da
        module procedure assign_2DdrA_2DA
    end interface

    interface allocate_dv
        !>real
        module procedure allocate_drv_array
        !>integer
        module procedure allocate_div_array
    end interface allocate_dv

    interface deallocate_dv
        !>real
        module procedure deallocate_drv_array
        !>integer
        module procedure deallocate_div_array
    end interface deallocate_dv

    !> \brief interface for MAVAL ADJOINT
    interface MAXVALADJ
        module procedure MAXVALADJ_1d
        module procedure MAXVALADJ_2d
    !     MODULE PROCEDURE MAXVALADJ_3d
    end interface MAXVALADJ

    !> @brief builds mask for regularly spaced coordinates
    interface build_regular_icoord_mask
        module procedure build_regular_icoord_mask_1D
        module procedure build_regular_icoord_mask_2D
        module procedure build_regular_icoord_mask_3D
    end interface build_regular_icoord_mask

    !> \brief computes the number of regularly-spaced coordinates
    interface regular_icoord_count
        module procedure regular_icoord_count_1D
        module procedure regular_icoord_count_2D
        module procedure regular_icoord_count_3D
    end interface regular_icoord_count

    !> @brief computes regularly-spaced coordinates
    interface build_regular_icoord
        module procedure build_regular_icoord_1D
        module procedure build_regular_icoord_2D
        module procedure build_regular_icoord_3D
    end interface build_regular_icoord



    !>Read variable from environment
    interface read_env
        module procedure read_env_int
        module procedure read_env_str
        module procedure read_env_dbl
    end interface read_env


!      interface minval
!         module procedure d2DrA_minval!< minval in dynamic 2D real vector
!         module procedure drv_minval!< minval in dynamic 1D real vector
!         module procedure drv_array_minval!< minval of an array of dynamic 1D real vector
!      end interface minval
!
!      interface maxval
!         module procedure d2DrA_maxval!< maxval in dynamic 2D real vector
!         module procedure drv_maxval!< maxval in dynamic 1D real vector
!         module procedure drv_array_maxval!< minval of an array of dynamic 1D real vector
!      end interface maxval

     !interface nullify
     !   module procedure drv_nullify!< nullify dynamic real vector
     !   module procedure dov_nullify!< nullify dynamic observation vector
     !end interface nullify

    integer, parameter :: OBS_FIELD = 0!< constant for data field of dynamic observation vector
    integer, parameter :: OBS_CMAT  = 1!< constant for covariance matrix field of dynamic observation vector
    integer, parameter :: OBS_XIDX  = 2!< constant for H field of dynamic observation vector

  !> \brief User defined type for environmental exchanged parameters between programs
  !! this data structure defines a set of file names and directory used for data assimilation purpose
  !! To keep coherence, each independant program should have only one copy of this data structure.
  !! @todo load and save ens_xx and eme_xx from namelist
  !<
  type environmental_exchange_param
    character(len=ip_snl) ::& !, PRIVATE (removed for compatibility with old compilers)
      default_ext= "nc"   ,& !<default file extension
      ep_bName   = "ep"   ,& !<exchange parameter file base name
      cost_bName = "cost" ,& !<cost file base name
      grad_bName = "grad" ,& !<grad file base name
      ctl_bName  = "ctl"  ,& !<control vector file base name
      bctl_bName = "b_ctl",& !<background ctl file base name
      obs_bName  = "obs"  ,& !< observation file base name
      bc_bName   = "bc"   ,& !<boundary conditions file base name
      dmt_bName  = "dmt"  ,& !<direct model traj file base name
      ims_bName  = "ims"  ,& !<initial model state file base name
      dms_bName  = "dms"  ,& !<direct model state file base name
      ens_bName      = "ensemble",& !<ensemble file base name
      ens_mean_bName = "ens_mean",& !<ens mean file base name
      ens_var_bName  = "ens_var",& !< ens var file base name
      ens_me_bName   = "ens_me",& !<ens model error file base name
      me_mean_bName  = "me_mean",& !<model error mean file base name
      me_var_bName   = "me_var",& !<model error var file base name
      ens_mean_sDir  = "",& !< ens mean subdirectory
      ens_var_sDir   = "",& !< ens var subdirectory
      ogap_bName     = "obsgap",& !< obsgap file base name
      ctl_all_bName  ="ctl_all"  ,& !<all ctl friends file base name
      ctl_bound_bName="ctl_bound",& !<ctl bounds file base name
      rt_dir     = "rt",& !<runtime directory, temporary runtime data are saved in this directory
      input_dir  = "." ,& !<input directory, input data are read from this directory
      output_dir = "."   !<output directory, output data are saved in this directory
  end type environmental_exchange_param

!   !>interface of functions used to build data assimilation file names
!   INTERFACE make_fileName
!     MODULE PROCEDURE make_fileName_old
!     MODULE PROCEDURE make_fileName_new
!   END INTERFACE make_fileName
!
!   !>interface of functions used to delete runtime files
!   INTERFACE delete_rtFile
!     MODULE PROCEDURE delete_rtFile_old
!     MODULE PROCEDURE delete_rtFile_new
!   END INTERFACE delete_rtFile
!
!   !>interface of functions used to set single file name in exchange parameter
!   INTERFACE set_fileName
!     MODULE PROCEDURE set_fileName_old
!     MODULE PROCEDURE set_fileName_new
!   END INTERFACE set_fileName
!
!   !>interface of functions used to set the name of the data files
!   INTERFACE set_data_fileNames
!     MODULE PROCEDURE set_data_fileNames_old
!     MODULE PROCEDURE set_data_fileNames_new
!   END INTERFACE set_data_fileNames
!
!   !>interface of functions used to set the name of the data files
!   INTERFACE set_data_dir
!     MODULE PROCEDURE set_data_dir_old
!     MODULE PROCEDURE set_data_dir_new
!   END INTERFACE set_data_dir

  TYPE(environmental_exchange_param), SAVE :: eep

    !> \brief dynamic real vector
    TYPE dyn_rVector
        REAL(KIND=cp), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_rVector

    !> \brief dynamic integer vector
    TYPE dyn_iVector
        INTEGER, DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_iVector

    !> \brief character vector
    TYPE dyn_sVector
        CHARACTER(LEN=ip_snl), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_sVector

    !> \brief 2D dynamic real array (array of arrays), java style 2D array
    TYPE dyn_2DrArray
        TYPE(dyn_rVector), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_2DrArray

    !> \brief 2D dynamic integer array (array of arrays), java style 2D array
    TYPE dyn_2DiArray
        TYPE(dyn_iVector), DIMENSION(:), POINTER :: dat => NULL()
    END TYPE dyn_2DiArray

    !> \brief observation vector at a given date
    !!
    !<
    TYPE dyn_oVector
        INTEGER :: ts !< time step
        LOGICAL :: H_allocated = .FALSE. !< is memory space allocated for idx!, PRIVATE
        LOGICAL :: Rm1_allocated = .FALSE. !< is memory space allocated for Rm1!, PRIVATE
        REAL(KIND=cp), DIMENSION(:), POINTER :: dat  => NULL() !< observation data
        INTEGER      , DIMENSION(:), POINTER :: H => NULL()!< observation operator, indexes of obs data in the spatial grid
        !> inverse of the observation diagonal covariance matrix, can be extended to full matrix with covariates components
        !<
        REAL(KIND=cp), DIMENSION(:), POINTER :: Rm1    => NULL()
    END TYPE dyn_oVector

CONTAINS

  !
  !**********************General routines*********************
  !

    !> @brief read an integer variable from the environment
    !! @param [in] vname variable name
    !! @param [out] var contains the read value
    !! @param [in] die flag saying to stop or not if the variable is not present
    !! @param [out] stat returns zero if succesful
    subroutine read_env_int(vname, var, die, stat)
        character(len=*), intent(in) :: vname
        integer, intent(out) :: var
        logical, optional, intent(in) :: die
        integer, optional, intent(out) :: stat
        !local variable
        character(len=ip_fnl) :: tmp
        integer :: env_stat
        logical :: ll_die

        if(present(die) )then
            ll_die = die
        else
            ll_die = .true.
        end if

        call get_environment_variable(vname, tmp, status=env_stat)
        if(env_stat==0)then
            read(tmp, *)var
        else if(ll_die)then
            call stop_program( "No env variable named "//trim(vname) )
        end if
        if( present(stat) ) stat=env_stat
    end subroutine read_env_int

    !> @brief read a double variable from the environment
    !! @param [in] vname variable name
    !! @param [out] var contains the read value
    !! @param [in] die flag saying to stop or not if the variable is not present
    !! @param [out] stat returns zero if succesful
    subroutine read_env_dbl(vname, var, die, stat)
        character(len=*), intent(in) :: vname
        real(dp), intent(out) :: var
        logical, optional, intent(in) :: die
        integer, optional, intent(out) :: stat
        !local variable
        character(len=ip_fnl) :: tmp
        integer :: env_stat
        logical :: ll_die

        if(present(die) )then
            ll_die = die
        else
            ll_die = .true.
        end if

        call get_environment_variable(vname, tmp, status=env_stat)
        if(env_stat==0)then
            read(tmp, *)var
        else if(ll_die)then
            call stop_program( "No env variable named "//trim(vname) )
        end if
        if( present(stat) ) stat=env_stat
    end subroutine read_env_dbl

    !> @brief read a character string variable from the environment
    !! @param [in] vname variable name
    !! @param [out] var contains the read value
    !! @param [in] die flag saying to stop or not if the variable is not present
    !! @param [out] stat returns zero if succesful
    subroutine read_env_str(vname, var, die, stat)
        character(len=*), intent(in) :: vname
        character(len=*), intent(in out) :: var
        logical, optional, intent(in) :: die
        integer, optional, intent(out) :: stat
        !local variable
        !character(len=ip_fnl) :: tmp
        integer :: env_stat
        logical :: ll_die

        if(present(die) )then
            ll_die = die
        else
            ll_die = .true.
        end if

        call get_environment_variable(vname, var, status=env_stat)
        if( (env_stat/=0).and.(ll_die) )then
            call stop_program( "No env variable named "//trim(vname) )
        end if
        if( present(stat) ) stat=env_stat
    end subroutine read_env_str

  !> \brief Sends a command to the shell as if it had been typed at the command line.
  !! PARAM[in] ada_command command to run
  !! PARAM[in] id_status optional argument to ruturn the running status
  !! If the id_status is present, it contains the running status of the command. Otherwise, the program aborts if the running status is different from 0
  !<
  SUBROUTINE system_run( ada_command, id_status )
    !USE IFPORT
    CHARACTER(LEN=*), INTENT(IN)   :: ada_command
    INTEGER, INTENT(OUT), OPTIONAL :: id_status
    !local var
    INTEGER il_status!, il_errnum

    CALL debug('Entering system_run ...................', tag=dTRACE)
    CALL debug(ada_command, '  command', tag=dTRACE)
    CALL SYSTEM( TRIM(ada_command) )!, il_status
    il_status = 0
    SELECT CASE (il_status)
      CASE (0) !success, nothing to do
      CASE (-1) !
        CALL debug(ada_command, 'system_run;  error running the command', tag=dALLWAYS)
        !il_errnum = ierrno( )
        !CALL debug(il_errnum, '  the error number is: ', tag=dALLWAYS)
      CASE DEFAULT
        CALL debug(ada_command, 'system_run;  abnormal things hapenned while running the command', tag=dALLWAYS)
    END SELECT

    IF( PRESENT(id_status) ) THEN
      id_status = il_status
    ELSE
      IF( il_status .NE. 0 ) THEN
        CALL debug(il_status, '  the return status from <'//TRIM(ada_command)//'> is', tag=dALLWAYS)
        CALL ABORT()
      END IF
    END IF
    CALL debug( 'Exiting system_run ...................', tag=dTRACE )

  END SUBROUTINE system_run

    function norm2_1d_dp(vec)result(norm)
        real(dp), dimension(:), intent(in) :: vec
        real(dp) :: norm
        norm = sqrt( dot_product(vec,vec) )
    end function norm2_1d_dp

    function norm2_2d_dp(mat)result(norm)
        real(dp), dimension(:, :), intent(in) :: mat
        real(dp) :: norm
        norm = sqrt( sum(mat**2) )
    end function norm2_2d_dp

    !> Check is a file exists
    !! @param [in] fName name of the file to inquire
    !<
    function fileExist(fName) result(ll_exist)
        CHARACTER(len=*), INTENT(IN) :: fName
        !Local variables
        logical ll_exist

        inquire(file=fName, exist=ll_exist)
    end function fileExist

    !> @brief print environmental exchange parameters
    !<
    subroutine print_eep()
        call debug('', 'environmental exchange parameters data structure--------', tag=dALLWAYS)
        call debug(eep%ep_bName,       'ep_bName       ', tag=dALLWAYS)
        call debug(eep%cost_bName,     'cost_bName     ', tag=dALLWAYS)
        call debug(eep%grad_bName,     'grad_bName     ', tag=dALLWAYS)
        call debug(eep%ctl_bName,      'ctl_bName      ', tag=dALLWAYS)
        call debug(eep%bctl_bName,     'bctl_bName     ', tag=dALLWAYS)
        call debug(eep%obs_bName,      'obs_bName      ', tag=dALLWAYS)
        call debug(eep%bc_bName,       'bc_bName       ', tag=dALLWAYS)
        call debug(eep%dmt_bName,      'dmt_bName      ', tag=dALLWAYS)
        call debug(eep%ims_bName,      'ims_bName      ', tag=dALLWAYS)
        call debug(eep%dms_bName,      'dms_bName      ', tag=dALLWAYS)
        call debug(eep%ens_bName,      'ens_bName      ', tag=dALLWAYS)
        call debug(eep%ens_mean_bName, 'ens_mean_bName ', tag=dALLWAYS)
        call debug(eep%ens_var_bName,  'ens_var_bName  ', tag=dALLWAYS)
        call debug(eep%ens_me_bName,   'ens_me_bName   ', tag=dALLWAYS)
        call debug(eep%me_mean_bName,  'me_mean_bName  ', tag=dALLWAYS)
        call debug(eep%me_var_bName,   'me_var_bName   ', tag=dALLWAYS)
        call debug(eep%ogap_bName,     'ogap_bName     ', tag=dALLWAYS)
        call debug(eep%ctl_all_bName,  'ctl_all_bName  ', tag=dALLWAYS)
        call debug(eep%ctl_bound_bName,'ctl_bound_bName', tag=dALLWAYS)
        call debug(eep%rt_dir,         'rt_dir         ', tag=dALLWAYS)
        call debug(eep%input_dir,      'input_dir      ', tag=dALLWAYS)
        call debug(eep%output_dir,     'output_dir     ', tag=dALLWAYS)
        call debug('', '--------------------------------------------------------', tag=dALLWAYS)
    end subroutine print_eep

    !> @brief Returns the data type associated to the statistical mean of data of type dType
    !! @param [in] dType type of data(ENS_DATA, ENS_ME_DATA, CTL_DATA, etc.)
    !! @details the program terminates if dType does not correspond to ensemble data type (ENS_DATA, ENS_ME_DATA)
    function mean_dType(dType)result(mean_type)
        INTEGER, INTENT(IN) :: dType
        INTEGER :: mean_type

        select case(dType)
        case(ENS_DATA)
            mean_type = ENS_MEAN_DATA
        case(ENS_ME_DATA)
            mean_type = ME_MEAN_DATA
        case default
            mean_type = ENS_MEAN_DATA
            call stop_program("mean_dType supports only ENS_DATA and ENS_ME_DATA")
        end select
    end function mean_dType

    !> @brief Joins two non empty strings with a given separator
    !! @param [in] sep separator
    !! @param [in] s1, s2 strings to be joined
    !! @details trailing spaces in \a s1 and \a s2 are not considered.
    !! No preprocessing is done on the separator. The user should make
    !! sure that it contains only the chosen character(s) and is of the right len.
    !! If one string is empty, the other one is returned
    !!
    !<
    function join2strings(sep, s1, s2)result(s)
        character(len=*), intent(in) :: sep, s1, s2
        character(len=ip_fnl) :: s
        integer :: l, l1, l2

        l1 = len_trim(s1)
        l2 = len_trim(s2)
        l = l1+l2
        if(l==0)then
            s = ""
        else if(l==l1)then
            s = trim(s1)
        else if(l==l2)then
            s = trim(s2)
        else
            s = trim(s1)//sep//adjustl( trim(s2) )
        end if
    end function join2strings

    !> @brief Joins up to 10 non empty strings with a given separator
    !! @param [in] sep separator
    !! @param [in] s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 strings to be joined
    !! @details trailing spaces in string to be joined are not considered.
    !! No preprocessing is done on the separator. The user should make sure
    !! that it contains only the chosen character(s) and is of the right len.
    !! Apart from the separator, trailing spaces are removed and
    !! empty strings are ignored while joining
    !<
    function join(sep, s0, s1, s2, s3, s4, s5, s6, s7, s8, s9)result(s)
        character(len=*), intent(in) :: sep
        character(len=*), optional, intent(in) :: s0, s1, s2, s3, s4&
            , s5, s6, s7, s8, s9
        character(len=ip_fnl) :: s

        if( present(s0) ) s = trim(s0)
        if(present(s1)) s = join2strings(sep, s, s1)
        if(present(s2)) s = join2strings(sep, s, s2)
        if(present(s3)) s = join2strings(sep, s, s3)
        if(present(s4)) s = join2strings(sep, s, s4)
        if(present(s5)) s = join2strings(sep, s, s5)
        if(present(s6)) s = join2strings(sep, s, s6)
        if(present(s7)) s = join2strings(sep, s, s7)
        if(present(s8)) s = join2strings(sep, s, s8)
        if(present(s9)) s = join2strings(sep, s, s9)
    end function join

    !> @brief Joins up to 10 parts of a path into one single path. Written for
    !!   unix-like system with the default separator set to '/'
    !! @param [in] p0, p1, p2, p3, p4, p5, p6, p7, p8, p9 paths to be joined
    !! @details trailing spaces in the subpath are removed.
    !! and empty subpaths are ignored
    !<
    function pJoin(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9)result(p)
        character(len=*), optional, intent(in) :: p0, p1, p2, p3, p4&
            , p5, p6, p7, p8, p9
        character(len=ip_fnl) :: p

        p = join('/', p0, p1, p2, p3, p4, p5, p6, p7, p8, p9)
    end function pJoin

    !> @brief Joins up to 10 strings into one single string with the underscore
    !!  character as a separator
    !! @param [in] s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 strings to be joined
    !! Apart from the separator, trailing spaces are removed and
    !! empty strings are ignored while joining
    !<
    function uJoin(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9)result(s)
    character(len=*), optional, intent(in) :: s0, s1, s2, s3, s4&
        , s5, s6, s7, s8, s9
        character(len=ip_fnl) :: s

        s = join('_', s0, s1, s2, s3, s4, s5, s6, s7, s8, s9)
    end function uJoin

    !> @brief Returns the file name prefix in the variable pref
    !! @param [in] dType type of data(ENS_DATA, ENS_ME_DATA, CTL_DATA, etc.)
    !! @param [in] pref file name prefix
    !<
    subroutine fName_prefix(dType, pref)
        INTEGER, INTENT(IN) :: dType
        CHARACTER(len=*), OPTIONAL, INTENT(out) :: pref

        select case(dType)
            case(AMT_DATA, ACTL_DATA, AOBS_DATA, AIMS_DATA, AMS_DATA)
                pref = APREF
            case(BMT_DATA, BCTL_DATA, BIMS_DATA, BMS_DATA)
                pref = BPREF
            case(DCTL_DATA)
                pref = DPREF
            case(SMT_DATA)
                pref = SPREF
            case(TMT_DATA, TCTL_DATA, TIMS_DATA, TMS_DATA)
                pref = TPREF
            case default
                pref = ""
        end select
    end subroutine fName_prefix

    !>short name for make_fileName
    function make_fName(id_dType, id_status, prefix, basename, num, ext, dir, subdir)result(ala_fName)
        INTEGER, INTENT(IN) :: id_dType, id_status
        CHARACTER(len=*), OPTIONAL, INTENT(IN) :: prefix, basename, ext, dir, subdir
        integer, intent(in), optional :: num
        !local variables
        CHARACTER(len=ip_fnl) :: ala_fName

        ala_fName = make_fileName(id_dType, id_status, prefix, basename, num, ext, dir, subdir)
    end function make_fName

   !> \brief make file name according to the usage
   !! @param [in] id_dType type of content of the file (CTL_DATA, etc.)
   !! @param [in] id_status file status (input, output, runtime or basename)
   !! @param [in] prefix (optional) prefix of the filename
   !! @param [in] basename (optional) basename of the filename, used when the
   !!   content type is OTHER_DATA
   !!
   !! @param [in] num optional number to be added as suffix to the filename.
   !!   \a num should be provided if the is a sequences of numbered files sharing
   !!   the same prefix name.
   !! @param [in] ext optional extension to be used as file extension
   !! @param [in] dir optional directory to be used if \a id_status is
   !!   \a BASE_FILE (shoulb \a USER_FILE)
   !! @param [in] subdir optional subdirectory to be used
   !!
   !! \details this function builds the appropriate file name of the form:
   !!  [dir/]prefix_basename[_num][.ext], given the data type to be stored in the file
   !!  and its usage. If prefix is given, it is used as prefix for the filename.
   !!  The prefix comes after the single letter prefix (a, t, d, etc.).
   !!  \a basename, \a num and \a ext are used only if the file type is
   !!  \a OTHER_DATA. Whenever \a num is provided, \a basename should not contain
   !!  the extension, and the extension should be provided by \a ext.
   !! @todo update old programs to account for the extension that has been separated
   !!  from the file name
   !<
   function make_fileName(id_dType, id_status, prefix, basename, num, ext, dir, subdir)result(ala_fName)
      integer, intent(in) :: id_dType, id_status
      character(len=*), optional, intent(in) :: prefix, basename, ext, dir, subdir
      integer, intent(in), optional :: num
      !local variables
      character(len=ip_fnl) :: ala_fName, ala_dir, ala_sDir, ala_baseName&
      , ala_prefix, ala_ext, ala_num, run_pref

      !processing the one-time prefix
      if( present(prefix) )then
         ala_prefix = trim( prefix )
      else
         ala_prefix = ""
      end if
      !processing the optional extension
      if( present(ext) )then
        ala_ext = trim( ext )
      else
         ala_ext = trim(eep%default_ext)
      end if
      !processing the optional num
      if( present(num) )then
         !write(ala_num,'(A,i3.3)')'_',num
         write(ala_num,'(i3.3)')num
      else
         ala_num = ''
      end if
      run_pref = ""
      ala_baseName = ""
      ala_sDir = ""
      !IF(LEN_TRIM(ala_prefix)>0) ala_prefix = TRIM(ala_prefix)//'_'

      select case(id_dType)
         case (CTL_DATA)
            ala_baseName = eep%ctl_bName
         case (BCTL_DATA)
            ala_baseName = eep%bctl_bName
         case (ACTL_DATA)
            ala_baseName = trim(eep%ctl_bName)
            run_pref = trim(APREF)
         case (DCTL_DATA)
            ala_baseName = trim(eep%ctl_bName)
            run_pref = trim(DPREF)
         case (TCTL_DATA)
            ala_baseName = trim(eep%ctl_bName)
            run_pref = trim(TPREF)
         case (CTL_BOUND_DATA)
            ala_baseName = eep%ctl_bound_bName
         case (CTL_ALL_DATA)
            ala_baseName = eep%ctl_all_bName
         case (OBS_DATA)
            ala_baseName = eep%obs_bName
         case (AOBS_DATA)
            ala_baseName = trim(eep%obs_bName)
            run_pref = trim(APREF)
         case (TOBS_DATA)
            ala_baseName = trim(eep%obs_bName)
            run_pref = trim(TPREF)
         case (BC_DATA)
            ala_baseName = eep%bc_bName
         case (DMT_DATA)
            ala_baseName = eep%dmt_bName
         case (AMT_DATA)
            ala_baseName = trim(eep%dmt_bName)
            run_pref = trim(APREF)
         case (BMT_DATA)
            ala_baseName = trim(eep%dmt_bName)
            run_pref = trim(BPREF)
         case (TMT_DATA)
            ala_baseName = trim(eep%dmt_bName)
            run_pref = trim(TPREF)
         case (SMT_DATA)
            ala_baseName = trim(eep%dmt_bName)
            run_pref = trim(SPREF)
         case (IMS_DATA)
            ala_baseName = eep%ims_bName
         case (BIMS_DATA)
            ala_baseName = trim(eep%ims_bName)
            run_pref = trim(BPREF)
         case (AIMS_DATA)
            ala_baseName = trim(eep%ims_bName)
            run_pref = trim(APREF)
         case (TIMS_DATA)
            ala_baseName = trim(eep%ims_bName)
            run_pref = trim(TPREF)
         case (DMS_DATA)
            ala_baseName = eep%dms_bName
         case (AMS_DATA)
            ala_baseName = trim(eep%dms_bName)
            run_pref = trim(APREF)
         case (BMS_DATA)
            ala_baseName = trim(eep%dms_bName)
            run_pref = trim(BPREF)
         case (TMS_DATA)
            ala_baseName = trim(eep%dms_bName)
            run_pref = trim(TPREF)
         case (ENS_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%ens_bName)
         case (ENS_MEAN_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%ens_mean_bName)
         case (ENS_VAR_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%ens_var_bName)
         case (ENS_ME_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%ens_me_bName)
         case (ME_MEAN_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%me_mean_bName)
         case (ME_VAR_DATA)!the processing of the ensemble data is different from other cases
            run_pref = trim(eep%me_var_bName)
         case (OGAP_DATA)
            ala_baseName = eep%ogap_bName
         case (EP_DATA)
            ala_baseName = eep%ep_bName
         case (COST_DATA)
            ala_baseName = eep%cost_bName
         case (GRAD_DATA)
            ala_baseName = eep%grad_bName
         case (OTHER_DATA)
            if( present(basename) )then
               ala_baseName = basename
            else
               call stop_program("In make_fileName: basename is required when DATA type is OTHER_DATA")
            end if
         case default
            call stop_program(id_dType, 'In make_fileName; bad data type: ')
      end select

      if( present(basename) )then
         ala_baseName = uJoin(run_pref, basename)
         if( .not.present(ext) ) ala_ext = ""
      else
         ala_baseName = uJoin(run_pref, ala_baseName)
      end if

      if( present(subdir) )then
        ala_sDir = trim (subdir)
      end if

      select case(id_status)
         case (INPUT_FILE)
            ala_dir = eep%input_dir
         case (OUTPUT_FILE)
            ala_dir = eep%output_dir
         case (RTIME_FILE)
            ala_dir = '.'
            ala_prefix =  uJoin( 'rt', ala_prefix)
         case (USER_FILE) !same as BASE_FILE
            if( present(dir) )then
                ala_dir = dir
            else
                ala_dir = '.'
            end if
         case default
            call stop_program(id_status, 'In make_fileName; file status: ')
      end select
      if( len_trim(ala_dir) == 0 )ala_dir = '.'

      ala_baseName = uJoin(ala_prefix, ala_baseName, ala_num)
      ala_baseName = join('.', ala_baseName, ala_ext)
      ala_fName = pJoin(ala_dir, ala_sDir, trim(ala_baseName)  )
      !ala_fName = TRIM(ala_dir)//'/'//TRIM(ala_prefix)&
      !              //TRIM(ala_baseName)//TRIM(ala_num)//TRIM(ala_ext)
   end function make_fileName

  !> @brief Builds stat (mean and variance) file names for ensemble data
  !! @param [out] mean_fName mean file name
  !! @param [out] var_fName variance file name
  !! @param [in] dType data type (ENS_DATA, ENS_ME_DATA)
  !! @param [in] status is the status of the file (INPUT_FILE, OUTPUT_FILE or USER_FILE)
  !! @param [in] pref optional prefix of file names (APREF, BPREF, FPREF)
  !! @param [in] basename (optional) basename of the filename, used when the
  !!   content type is OTHER_DATA
  !! @param [in] num optional file ordinal number
  !! @param [in] dir optional directory to be used is \a status is \a BASE_FILE
  !!   (should \a USER_FILE)
  !! @param [in] ext optional extension for the file name
  !<
  subroutine make_stat_fNames(mean_fName, var_fName, dType, status, pref, basename, num, dir, ext)
    character(len=ip_fnl), intent(out):: mean_fName, var_fName
    integer, intent(in) :: dType, status
    character(len=*), optional, intent(in) :: pref, basename, dir, ext
    integer, optional, intent(in) :: num

    select case (dType)
      case(ENS_DATA)
         mean_fName = make_fName(ENS_MEAN_DATA, status, prefix=pref, basename=basename, num=num, dir=dir, ext=ext)
         var_fName  = make_fName(ENS_VAR_DATA , status, prefix=pref, basename=basename, num=num, dir=dir, ext=ext)
      case(ENS_ME_DATA)
        mean_fName = make_fName(ME_MEAN_DATA, status, prefix=pref, basename=basename, num=num, dir=dir, ext=ext)
        var_fName  = make_fName(ME_VAR_DATA , status, prefix=pref, basename=basename, num=num, dir=dir, ext=ext)
      case default
        call stop_program("In make_stat_fNames, accept only ENS_DATA and ENS_ME_DATA")
    end select
  end subroutine make_stat_fNames

  !> \brief delete a runtime file
  !! @param [in] id_dType type of content of the file (CTL_DATA, etc.)
  !<
  SUBROUTINE delete_rtFile(id_dType)
    INTEGER, INTENT(IN) :: id_dType
    !local variables
    CHARACTER(len=ip_fnl) :: ala_fName

    ala_fName = make_fileName(id_dType, RTIME_FILE)
    CALL system( 'rm -f '//ala_fName )
  END SUBROUTINE delete_rtFile

  !> \brief set single file name in exchange parameter
  !! @param [in] ada_bName file name
  !! @param [in] id_dType data type (CTL_DATA, )
  !<
  SUBROUTINE set_fileName(id_dType, ada_bName)
    INTEGER, INTENT(IN) :: id_dType
    CHARACTER(LEN=*), INTENT(IN) :: ada_bName

    SELECT CASE(id_dType)
      CASE (CTL_DATA)
        eep%ctl_bName  = ADJUSTL(ada_bName)
      CASE (BCTL_DATA)
        eep%bctl_bName = ADJUSTL(ada_bName)
      CASE (CTL_BOUND_DATA)
        eep%ctl_bound_bName = ADJUSTL(ada_bName)
      CASE (CTL_ALL_DATA)
        eep%ctl_all_bName   = ADJUSTL(ada_bName)
      CASE (OBS_DATA)
        eep%obs_bName  = ADJUSTL(ada_bName)
      CASE (BC_DATA)
        eep%bc_bName  = ADJUSTL(ada_bName)
      CASE (DMT_DATA)
        eep%dmt_bName  = ADJUSTL(ada_bName)
      CASE (DMS_DATA)
        eep%dms_bName  = ADJUSTL(ada_bName)
      CASE (IMS_DATA)
        eep%ims_bName  = ADJUSTL(ada_bName)
      CASE (OGAP_DATA)
        eep%ogap_bName = ADJUSTL(ada_bName)
      CASE (EP_DATA)
        eep%ep_bName   = ADJUSTL(ada_bName)
      CASE (COST_DATA)
        eep%cost_bName = ADJUSTL(ada_bName)
      CASE (GRAD_DATA)
        eep%grad_bName = ADJUSTL(ada_bName)
      CASE DEFAULT
        WRITE(*,*) 'In set_fileName; bad data type or no field defined :', id_dType
        WRITE(*,*) '  Stopping the program ...'
        STOP
    END SELECT
  END SUBROUTINE set_fileName

    !> \brief set the base name of the data files
    !<
    SUBROUTINE set_data_fileNames(ada_ctl_bName, ada_bctl_bName,&
                ada_obs_bName, ada_dmt_bName, ada_ims_bName, ada_ogap_bName&
    )
        CHARACTER(LEN=*), INTENT(IN) :: ada_ctl_bName, ada_bctl_bName, ada_obs_bName,&
                                        ada_dmt_bName, ada_ims_bName , ada_ogap_bName

        CALL set_fileName(CTL_DATA , ada_ctl_bName)
        CALL set_fileName(BCTL_DATA, ada_bctl_bName)
        CALL set_fileName(OBS_DATA , ada_obs_bName)
        CALL set_fileName(DMT_DATA , ada_dmt_bName)
        CALL set_fileName(IMS_DATA , ada_ims_bName)
        CALL set_fileName(OGAP_DATA, ada_ogap_bName)
        !CALL set_fileName(BC_DATA , ada_bc_bName)
    END SUBROUTINE set_data_fileNames

  !> \brief set the data directories
  !! @param [in] ada_input_dir input directory
  !! @param [in] ada_output_dir output directory
  !! \details empty path is replaced by local directory
  !<
  SUBROUTINE set_data_dir(ada_input_dir, ada_output_dir)
    CHARACTER(LEN=*), INTENT(IN) :: ada_input_dir, ada_output_dir
    !local variables
    INTEGER :: il_sepPos, il_trimSize

    CALL debug('Entering set_data_dir ***********************')
    !input dir
    il_trimSize = LEN_TRIM(ada_input_dir)
    il_sepPos   = INDEX(ada_input_dir, '/', BACK = .TRUE.)! Last occurence
    IF ( (il_sepPos > 0).AND.(il_sepPos == il_trimSize) )THEN !Trailing separator
       eep%input_dir = ADJUSTL( ada_input_dir(:il_sepPos-1) )
    ELSE
       eep%input_dir = ADJUSTL( ada_input_dir )
    END IF
    IF (LEN_TRIM(eep%input_dir)==0) eep%input_dir = '.' !Replace empty string by working dir

    !output dir
    il_trimSize = LEN_TRIM(ada_output_dir)
    il_sepPos   = INDEX(ada_output_dir, '/', BACK = .TRUE.) !last occurence
    IF ( (il_sepPos > 0).AND.(il_sepPos == il_trimSize) )THEN !trailing separator
       eep%output_dir = ADJUSTL( ada_output_dir(:il_sepPos-1) )
    ELSE
       eep%output_dir = ADJUSTL( ada_output_dir )
    END IF
    IF (LEN_TRIM(eep%output_dir)==0)THEN !empty directory path
       eep%output_dir = '.'! replace empty string by working dir
    ELSE
       CALL SYSTEM('mkdir -p '//eep%output_dir) ! create the directory if it does not exist
    END IF
    CALL debug(eep%input_dir, 'eep%input_dir')
    CALL debug(eep%output_dir, 'eep%output_dir')
    CALL debug('Exiting set_data_dir  ***********************')
  END SUBROUTINE set_data_dir

    !> @brief set statistics subdirectory
    !! @param [in] mean_subdir name of ensemble mean subdirectory
    !! @param [in] var_subdir name of ensemble var subdirectory
    !!
    !! @details statistics are computed in Ensemble Kalman filters
    !!   - ensemble mean
    !!   - ensemble var
    !! The subdirectoris are relative to the directory where the ensemble is saved
    !<
    subroutine set_stat_subdir(mean_subdir, var_subdir)
        character(len=*), intent(in) :: mean_subdir, var_subdir
        !local variables
        integer :: il_sepPos, il_trimSize

        call debug('', 'Entering set_data_dir ***********************')
        !ensemble mean subdir
        il_trimSize = len_trim(mean_subdir)
        il_sepPos   = index(mean_subdir, '/', BACK = .true.)! Last occurence
        if ( (il_sepPos > 0).and.(il_sepPos == il_trimSize) )then !Trailing separator
            eep%ens_mean_sDir = adjustl( mean_subdir(:il_sepPos-1) )
        else
            eep%ens_mean_sDir = adjustl( trim(mean_subdir) )
        end if
        if (len_trim(eep%ens_mean_sDir)==0) eep%ens_mean_sDir = '.' !Replace empty string by current dir
        !
        !ensemble var subdir
        il_trimSize = len_trim(var_subdir)
        il_sepPos   = index(var_subdir, '/', BACK = .true.)! Last occurence
        if ( (il_sepPos > 0).and.(il_sepPos == il_trimSize) )then !Trailing separator
            eep%ens_var_sDir = adjustl( var_subdir(:il_sepPos-1) )
        else
            eep%ens_var_sDir = adjustl( trim(var_subdir) )
        end if
        if (len_trim(eep%ens_var_sDir)==0) eep%ens_var_sDir = '.' !Replace empty string by current dir

    end subroutine set_stat_subdir

    !> \brief make directory name associated to an experiment
    !! @param [in] id_num numerical identification of the experiment
    !! @param [out] prefix prefix of the directory name
    !! @param [in] dirName directory name
    !! this subroutine make string that can be used for the directory name in a multi experiment
    !<
    SUBROUTINE make_DirName(id_num, prefix, dirName)
        INTEGER, INTENT(IN) :: id_num
        CHARACTER(LEN=*), INTENT(IN) :: prefix
        CHARACTER(LEN=*), INTENT(INOUT) :: dirName
        IF(id_num<0) THEN
            WRITE(dirName, FMT='(A,I4.3)') TRIM(ADJUSTL(prefix)), id_num
        ELSE
            WRITE(dirName, FMT='(A,I4.4)') TRIM(ADJUSTL(prefix)), id_num
        END IF
    END SUBROUTINE make_DirName

    !> \brief get the data directory
    !! @param [in] id_status file status (input, output, or runtime)
    !<
    FUNCTION get_data_dir(id_status) RESULT(ala_dir)
        INTEGER, INTENT(IN)   :: id_status
        CHARACTER(LEN=ip_snl) :: ala_dir


        SELECT CASE(id_status)
            CASE (INPUT_FILE)
                ala_dir = eep%input_dir
            CASE (OUTPUT_FILE)
                ala_dir = eep%output_dir
            CASE (RTIME_FILE)
                ala_dir = './'
            CASE DEFAULT
                CALL stop_program(id_status, 'In make_fileName; file status: ')
        END SELECT
    END FUNCTION get_data_dir

    !> \brief generates history string for namelist data file
    !<
    FUNCTION make_history(ada_text) RESULT(ala_history)
        CHARACTER(LEN=*), INTENT(IN) :: ada_text
        !local variables
        CHARACTER(LEN=ip_fnl) :: ala_history
        CHARACTER(LEN=8)  :: d
        CHARACTER(LEN=10) :: t
        CHARACTER(LEN=5)  :: zone


        CALL date_and_time(d, t, zone)
        WRITE(ala_history, *)&
        SD, 'Automatically generated by ', TRIM(ada_text)&
            , ' on ', d(1:4)//'/'//d(5:6)//'/'//d(7:8)&
            , ' at ', t(1:2)//':'//t(3:4)//':'//t(5:6)//t(7:10), ' GMT', zone, SD
    END FUNCTION make_history

    !> Check if a given file name ends with '.nc'
    !! @param [in] fName name of the file to be checked
    function is_ncFile(fName)result(is_nc)
        character(len=*), intent(in) :: fName
        !local variables
        logical is_nc

        is_nc =( (index(fName, '.nc')+2)== len_trim(fName) )
    end function is_ncFile

    !>\brief make a plot file name from a data file name
    !<
    FUNCTION make_pFileName_from_fName(ada_fName) RESULT(ala_pFName)
        CHARACTER(LEN=*), PARAMETER :: apa_suff = '_plot.dat'
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        INTEGER :: il_idx
        CHARACTER(LEN=ip_fnl) :: ala_pFName

        il_idx = INDEX(ada_fName, '.dat')
        IF(il_idx <= 1) il_idx = INDEX(ada_fName, '.nc')
        IF(il_idx <= 1)THEN
            ala_pFName = apa_suff
        ELSE
            ala_pFName = TRIM( ada_fName(:il_idx-1) )//apa_suff
        END IF
    END FUNCTION make_pFileName_from_fName

    !>\brief make a plot file name from a data file description
    !! @param [in] id_dType type of content of the file (CTL_DATA, etc.)
    !! @param [in] id_status file status (input, output, or runtime)
    !! @param [in] prefix (optional) prefix of the filename
    !! @param [in] basename (optional) basename of the filename, used when the content type is OTHER_DATA
    !! @param [in] subdir (optional) subdirectory
    !<
    FUNCTION make_pFileName_from_scalar(id_dType, id_status, prefix, basename, subdir) RESULT(ala_pFName)
        INTEGER, INTENT(IN) :: id_dType, id_status
        CHARACTER(len=*), OPTIONAL, INTENT(IN) :: prefix, basename, subdir
        !local variables
        CHARACTER(LEN=ip_fnl) :: ala_pFName

ala_pFName = make_pFileName_from_fName( make_fileName(id_dType, id_status, prefix, basename, subdir=subdir) )
    END FUNCTION make_pFileName_from_scalar

    !> \brief print Elapsed Time from system clock in human readable format
    !! @param [in] prompt prompt message before the elapse time
    !! @param [in] initial initial number of clock ticks, resulting from a call to system_time(initial, rate)
    !! @param [in] final final number of clock ticks, resulting from a call to system_time(final, rate)
    !! @param [in] rate optional rate that gives the number of clock ticks per seconds
    !<
    subroutine printElapseTime_system_clock(prompt, initial, final, rate)
        character(len=*) :: prompt
        integer(kind=dp), intent(in) :: initial, final, rate
        !local variables
        real(dp) :: r

        r = real(rate, dp)
        call printElapseTime(prompt, initial/r, final/r)
    end subroutine printElapseTime_system_clock

    !> \brief print Elapsed Time from cpu_time in human readable format
    !! @param [in] prompt prompt message before the elapse time
    !! @param [in] initial initial time, resulting from a call to cpu_time(initial)
    !! @param [in] final final time, resulting from a call to cpu_time(final)
    !<
    SUBROUTINE printElapseTime_cpu_time(prompt, initial, final)
        INTEGER, PARAMETER ::&
        MINUTE = 60 &        !duration of a minute in seconds
        , HOUR = 60*MINUTE & !duration of an hour in seconds
        , DAY  = 24*HOUR     !duration of a day in seconds

        character(len=*) :: prompt
        real(cp), intent(in) :: initial, final
        real(cp) :: duration, s
        integer :: m, h, d
        character(len=20) :: mm, hh, dd, ss
        logical :: nonEmpty

        duration = final - initial
        d = floor(duration/DAY)
        duration = duration - int(d*DAY)
        h = floor(duration/HOUR)
        duration = duration - int(h*HOUR)
        m = floor(duration/MINUTE)
        duration = duration - int(m*MINUTE)
        s = duration

        nonEmpty = .FALSE.

        IF(d>0)THEN
            WRITE(dd, "("//IFORMAT//" ' day(s)')") d
            nonEmpty = .TRUE.
        ELSE
            dd = ""
        END IF
        IF(nonEmpty) dd = ADJUSTL(dd)

        IF(h>0)THEN
            WRITE(hh, "((I2)' hour(s)')") h
            nonEmpty = .TRUE.
        ELSE
            IF(nonEmpty)THEN
                hh = "0 hour"
            ELSE
                hh = ""
            END IF
        END IF
        IF(nonEmpty) hh = ADJUSTL(hh)

        IF(m>0)THEN
            WRITE(mm, "((I2)' minute(s)')") m
            nonEmpty = .TRUE.
        ELSE
            IF(nonEmpty)THEN
                mm = "0 minute"
            ELSE
                mm = ""
            END IF
        END IF
        IF(nonEmpty) mm = ADJUSTL(mm)

        IF(s>0)THEN
            WRITE(ss, "((F6.3)' second(s)')") s
        ELSE
            ss = "0 second"
        END IF
        IF(nonEmpty) ss = ADJUSTL(ss)
        !PRINT*, 'CPU time: '//TRIM(dd) // ' ' // TRIM(hh) // ' ' // TRIM(mm) // ' ' // TRIM(ss)
        PRINT*, join(' ', prompt, dd, hh, mm, ss)
    END SUBROUTINE printElapseTime_cpu_time

    !> \brief converting a string to uppercase
    !! From Dr. David G. Simpson web page
    !<
    subroutine uppercase(str)
        character(len=*), intent(in out) :: str
        !local var
        integer :: i, del

        del = iachar('a') - iachar('a')
        do i = 1, len_trim(str)
            if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
                str(i:i) = achar(iachar(str(i:i)) - del)
            end if
        end do
    end subroutine uppercase

    !> \brief converting a string to uppercase
    !<
    function toUpper(str)result(str_upper)
        character(len=*), intent(in) :: str
        !local var
        character(len=len(str)) :: str_upper
        integer :: i, del

        str_upper = str
        del = iachar('a') - iachar('a')
        do i = 1, len(str)
            if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
                str_upper(i:i) = achar(iachar(str(i:i)) - del)
            end if
        end do
    end function toUpper

    !> \brief converting a string to lowercase
    !! @param [in, out] str string to be converted
    !! From Dr. David G. Simpson web page
    !<
    subroutine lowercase(str)
        character(len=*), intent(in out) :: str
        !local variable
        integer :: i, del

        del = iachar('a') - iachar('a')

        do i = 1, len_trim(str)
            if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
                str(i:i) = achar(iachar(str(i:i)) + del)
            end if
        end do
    end subroutine lowercase

    !> \brief converting a string to lowercase
    !! @param [in, out] str string to be converted
    !<
    function toLower(str)result(str_lower)
        character(len=*), intent(in) :: str
        !local variable
        character(len=len(str)) :: str_lower
        integer :: i, del

        str_lower = str
        del = iachar('a') - iachar('a')
        do i = 1, len_trim(str)
            if (lge(str(i:i),'a') .and. lle(str(i:i),'z')) then
                str_lower(i:i) = achar(iachar(str(i:i)) + del)
            end if
        end do
    end function toLower

    !> @brief convert an logical value to integer
    !! @param [in] ld_logical logical value to convert
    !!
    !! @details The convention is that .true. corresponds to 1
    !!   and .false. corresponds to 0. The convention changes
    !!   when converting from integer to logical, in which case
    !!   any nonzero is treated as .true.
    !<
    elemental function logical2int(ld_logical) result(il_int)
        logical, intent(in) :: ld_logical
        integer :: il_int

        if (ld_logical)then
            il_int = 1
        else
            il_int = 0
        end if
    end function logical2int

    !> @brief convert an logical value to single precision real
    !! @param [in] ld_logical logical value to convert
    !!
    !! @details The convention is that .true. corresponds to 1
    !!   and .false. corresponds to 0. The convention changes
    !!   when converting from single precision to logical, in which case
    !!   any number greather than epsilon in absolute value is treated as .true.
    !!   and any number less than epsilon in absolute value is treated as .false.
    !<
    elemental function logical2sp(ld_logical) result(rd_sp)
        logical, intent(in) :: ld_logical
        real(sp) :: rd_sp

        if (ld_logical)then
            rd_sp = 1.0_sp
        else
            rd_sp = 0.0_sp
        end if
    end function logical2sp

    !> @brief convert an integer value to logical
    !! @param [in] id_int integer value to convert
    !!
    !! @details The convention is that any nonzero corresponds to .true.
    !<
    elemental function int2l(id_int) result(ll_logical)
        integer, intent(in) :: id_int
        logical :: ll_logical

        ll_logical = ( id_int/=0 )
    end function int2l

    !> @brief convert a single precision real value to logical
    !! @param [in] rd_val value to convert
    !!
    !! @details The convention is that any number greater than
    !!   epsilon corresponds to .true.
    !<
    elemental function sp2l(rd_val) result(ll_logical)
        real(sp), intent(in) :: rd_val
        logical :: ll_logical

        ll_logical = ( abs(rd_val).gt.epsilon(1.0_sp) )
    end function sp2l

    subroutine select_sort_int(ida_x)
        integer, dimension(:), intent(inout) :: ida_x
        integer :: ibi, ibj, il_min_idx, il_size

        il_size = size(ida_x)
        do ibi = 1, il_size-1
            il_min_idx = ibi
        do ibj = ibi+1, il_size
            if( ida_x(ibj)<ida_x(il_min_idx) ) il_min_idx = ibj
        end do
            if(ibi /= il_min_idx) then
                call switch( ida_x(ibi), ida_x(il_min_idx) )
            end if
        end do
    end subroutine select_sort_int

    subroutine switch_int(id_a, id_b)
        integer, intent(inout) :: id_a, id_b
        !local variables
        integer :: il_tmp

        il_tmp = id_a
        id_a = id_b
        id_b = il_tmp
    end subroutine switch_int

    !> \brief Returns the location of an integer in an ordered 1D array, or -1 if the value is not in the array
    !! @param [in] ida_vec input vector in which the search is performed
    !! @param [in] id_val value to search
    !<
    function location_inc(id_val, ida_vec)result(loc)
        integer, dimension(:), intent(in) :: ida_vec
        integer, intent(in) :: id_val
        integer :: loc

        loc = search_inc_int(id_val, ida_vec )
        if( ida_vec(loc)/=id_val) loc=-1
    end function location_inc

    !> \brief Search the position of an integer value in an ordered 1D array using dichotomial search
    !! @param [in] ida_vec input vector in which the search is performed
    !! @param [in] id_val value to search
    !! @param [out] ntry number of tries before finding the value
    !! the array is supposed to be sorted in increasing order
    !! The position of the largest element inferior or equal to val is returned; if the searched value is less than the smallest element of the array 1 is returned
    !<
    function search_inc_int(id_val, ida_vec, ntry )result(il_pos)
        integer, dimension(:), intent(in) :: ida_vec
        integer, intent(in) :: id_val
        integer, optional, intent(out) :: ntry
        integer :: il_pos, il_left, il_right, il_mid, il_ntry

        il_ntry = 0
        il_left = 1
        il_right = size(ida_vec)
        do while(il_left<il_right)
            il_mid = (il_left + il_right + 1)/2
            if( id_val < ida_vec(il_mid) ) then
                il_right = il_mid-1
            else
                il_left = il_mid
            end if
            il_ntry = il_ntry + 1
        end do
        il_pos = il_left
        if ( present(ntry) ) ntry = il_ntry
    end function search_inc_int

    !> \brief Search the position of a value in an ordered 1D array using dichotomial search (Double precision version)
    !! @param [in] rda_vec input vector in which the search is performed
    !! @param [in] rd_val value to search
    !! @param [out] ntry number of tries before finding the value
    !! the array is supposed to be sorted in increasing order
    !! The position of the largest element inferior or equal to val is returned; if the searched value is less than the smallest element of the array 1 is returned
    !<
    function search_inc_dp_real(rd_val, rda_vec, ntry )result(il_pos)
        real(kind=cp), dimension(:), intent(in) :: rda_vec
        real(kind=cp), intent(in) :: rd_val
        integer, optional, intent(out) :: ntry
        !local variables
        integer :: il_pos, il_left, il_right, il_mid, il_ntry

        il_ntry = 0
        il_left = 1
        il_right = size(rda_vec)
        do while(il_left<il_right)
            il_mid = (il_left + il_right + 1)/2
            if( rd_val < rda_vec(il_mid) ) then
                il_right = il_mid-1
            else
                il_left = il_mid
            end if
            il_ntry = il_ntry + 1
        end do
        il_pos = il_left
        if ( present(ntry) ) ntry = il_ntry
    end function search_inc_dp_real

    !> \brief Search the position of a value in an ordered 1D array using dichotomial search
    !! @param [in] rda_vec input vector in which the search is performed
    !! @param [in] rd_val value to search
    !! @param [out] ntry number of tries before finding the value
    !! the array is supposed to be sorted in decreasing order
    !! The position of the largest element superior or equal to val is returned; if the searched value is less than the smallest element of the array, SIZE(array) is returned
    !<
    function search_dec_dp_real(rd_val, rda_vec, ntry )result(il_pos)
        real(kind=cp), dimension(:), intent(in) :: rda_vec
        real(kind=cp), intent(in) :: rd_val
        integer, optional, intent(out) :: ntry
        integer :: il_pos, il_left, il_right, il_mid, il_ntry

        il_ntry = 0
        il_left = 1
        il_right = size(rda_vec)
        do while(il_left<il_right)
            il_mid = (il_left + il_right + 1)/2
            if( rd_val > rda_vec(il_mid) ) then
                il_right = il_mid-1
            else
                il_left = il_mid
            end if
            il_ntry = il_ntry + 1
        end do
        il_pos = il_left
        if ( present(ntry) ) ntry = il_ntry
    end function search_dec_dp_real

    !> @brief write 2D array data for plot in the friendly-plot format
    !! @param [in] rda_data vector to be saved
    !! @param [in] ada_fName path to the filename
    !!
    !! \details the 2D array is saved column by column
    !<
    subroutine write_plot_data_2D_default(rda_data, ada_fName)
        real(cp), dimension(:,:), intent(in) :: rda_data
        character(len=*), intent(in) :: ada_fName
        !local variables
        real(cp), dimension( size(rda_data,1) ) :: rla_x
        integer :: ibi, il_nb
        il_nb = size(rda_data,1)
        rla_x = real( (/(ibi, ibi=1,il_nb)/), cp )

        call write_plot_data_2D_rcoord(ada_fName, rda_data, rla_x)
    end subroutine write_plot_data_2D_default

    !> \brief write 2D array data in the friendly-plot format
    !! @param [in] rda_f 2D array to be saved
    !! @param [in] ada_fName file name
    !! @param [in] rda_x coordinates of points
    !! the array is supposed to contain the values of a vector function at given points. each column of \a rda_f is a component of the vector function
    !<
    SUBROUTINE write_plot_data_2D_rcoord(ada_fName, rda_f, rda_x)
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: rda_f
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        INTEGER :: i
        TYPE(dyn_rVector), DIMENSION(size(rda_f,2)+1) :: tla_save
        tla_save(1) = rda_x
        do i=2,size(tla_save)
            tla_save(i) = rda_f(:,i-1)
        end do

        CALL save_trj( tla_save, make_pFileName(ada_fName), ld_column = .TRUE. )
        CALL debug(make_pFileName(ada_fName), ' in write_vector_for_plot_1D_rcoord -->', tag=dALLWAYS)
        CALL reset_drv_array( tla_save )
    END SUBROUTINE write_plot_data_2D_rcoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_data vector to be saved
    !! @param [in] ada_fName file name
    !<
    SUBROUTINE write_vector_for_plot_default(rda_data, ada_fName)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_data
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        call write_plot_data_1D_default(rda_data, ada_fName)
    END SUBROUTINE write_vector_for_plot_default

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_data vector to be saved
    !! @param [in] ada_fName file name
    !<
    SUBROUTINE write_plot_data_1D_default(rda_data, ada_fName)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_data
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        REAL(cp), DIMENSION( SIZE(rda_data) ) :: rla_x
        INTEGER :: ibi, il_nb
        il_nb = SIZE(rda_data)
        rla_x = REAL( (/(ibi, ibi=1,il_nb)/), cp )

        CALL write_vector_for_plot_1D_rcoord(ada_fName, rda_data, rla_x)
    END SUBROUTINE write_plot_data_1D_default

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_data vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] rd_xmin min coordinate
    !! @param [in] rd_xmax max coordinate
    !! the vector is supposed to contain the values of a function at regularly spaced points between rd_xmin and rd_xmax (inclusive).
    !<
    SUBROUTINE write_vector_for_plot_domain(rda_data, ada_fName, rd_xmin, rd_xmax)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_data
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        REAL(cp), INTENT(IN) :: rd_xmin, rd_xmax
        call write_plot_data_1D_domain(rda_data, ada_fName, rd_xmin, rd_xmax)
    END SUBROUTINE write_vector_for_plot_domain

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_data vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] rd_xmin min coordinate
    !! @param [in] rd_xmax max coordinate
    !! the vector is supposed to contain the values of a function at regularly spaced points between rd_xmin and rd_xmax (inclusive).
    !<
    SUBROUTINE write_plot_data_1D_domain(rda_data, ada_fName, rd_xmin, rd_xmax)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_data
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        REAL(cp), INTENT(IN) :: rd_xmin, rd_xmax
        !local variables
        REAL(cp), DIMENSION( SIZE(rda_data) ) :: rla_x
        INTEGER :: ibi, il_nb
        il_nb = SIZE(rda_data)
        rla_x = REAL( (/(ibi, ibi=0,il_nb-1)/), cp )/(rd_xmax-rd_xmin) + rd_xmin

        CALL write_vector_for_plot_1D_rcoord(ada_fName, rda_data, rla_x)
    END SUBROUTINE write_plot_data_1D_domain

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] ida_x indexes of points
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_vector_for_plot_1D_icoord(ada_fName, rda_f, ida_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f
        INTEGER , DIMENSION(:), INTENT(IN) :: ida_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName

        call write_plot_data_1D_icoord(ada_fName, rda_f, ida_x)
    END SUBROUTINE write_vector_for_plot_1D_icoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] ida_x indexes of points
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_plot_data_1D_icoord(ada_fName, rda_f, ida_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f
        INTEGER , DIMENSION(:), INTENT(IN) :: ida_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        REAL(cp), DIMENSION( SIZE(rda_f) ) :: rla_x

        rla_x = REAL( ida_x, cp )
        CALL write_vector_for_plot_1D_rcoord(ada_fName, rda_f, rla_x)
    END SUBROUTINE write_plot_data_1D_icoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] rda_x coordinates of points
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_vector_for_plot_1D_rcoord(ada_fName, rda_f, rda_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f, rda_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName

        call write_plot_data_1D_rcoord(ada_fName, rda_f, rda_x)
    END SUBROUTINE write_vector_for_plot_1D_rcoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] ada_fName file name
    !! @param [in] rda_x coordinates of points
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_plot_data_1D_rcoord(ada_fName, rda_f, rda_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f, rda_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        INTEGER :: il_nb
        TYPE(dyn_rVector), DIMENSION(2) :: tla_save
        il_nb = SIZE(rda_f)
        tla_save(1) = rda_x
        tla_save(2) = rda_f

        CALL save_trj( tla_save, make_pFileName(ada_fName), ld_column = .TRUE. )
        CALL debug(make_pFileName(ada_fName), ' in write_vector_for_plot_1D_rcoord -->', tag=dALLWAYS)
        CALL reset_drv_array( tla_save )
    END SUBROUTINE write_plot_data_1D_rcoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] rda_x 2D array of coordinates (the first dimension gives the number of coordinates per component)
    !! @param [in] ada_fName file name
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_vector_for_plot_ND_rcoord(ada_fName, rda_f, rda_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: rda_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName

        call write_plot_data_1D_ND_rcoord(ada_fName, rda_f, rda_x)
    END SUBROUTINE write_vector_for_plot_ND_rcoord

    !> \brief write vector in the friendly-plot format
    !! @param [in] rda_f vector to be saved
    !! @param [in] rda_x 2D array of coordinates (the first dimension gives the number of coordinates per component)
    !! @param [in] ada_fName file name
    !! the vector is supposed to contain the values of a function at given points. Those points can be given as an integer vector or a real vector.
    !<
    SUBROUTINE write_plot_data_1D_ND_rcoord(ada_fName, rda_f, rda_x)
        REAL(cp), DIMENSION(:), INTENT(IN) :: rda_f
        REAL(cp), DIMENSION(:,:), INTENT(IN) :: rda_x
        CHARACTER(LEN=*), INTENT(IN) :: ada_fName
        !local variables
        INTEGER :: il_nb, ibi
        TYPE(dyn_rVector), DIMENSION(SIZE(rda_x,1)+1) :: tla_save
        il_nb = SIZE(rda_f)
        IF(SIZE(rda_x,2)/=il_nb)THEN
            CALL debug('the number of row in x should be equal to the size of f',tag=dALLWAYS)
            CALL stop_program()
        END IF
        DO ibi=1, SIZE(rda_x,1)
            tla_save(ibi) = rda_x(ibi, :)
        END DO
        tla_save(SIZE(rda_x,1)+1) = rda_f

        CALL save_trj( tla_save, make_pFileName(ada_fName), ld_column = .TRUE. )
        CALL reset_drv_array( tla_save )
    END SUBROUTINE write_plot_data_1D_ND_rcoord

    !> \brief assign observation data to a dynamic observation vector
    !! @param [in, out] td_dov dynamic observation vector (structure)
    !! @param [in] rda_data observation data
    !! @param [in] ida_H indexes of observation data in spatial grid
    !! @param [in] rda_Rm1 diagonal matrix, inverse of the observation's covariance error (can be extended to 2D)
    !! \details ida_idx and rda_Rm1 are optional since they can be identic for all observation date in time
    !<
    SUBROUTINE assign_dov_vector(td_dov, rda_data, rda_Rm1, ida_H)
        TYPE(dyn_oVector), INTENT(INOUT)        :: td_dov
        REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: rda_Rm1
        INTEGER,       DIMENSION(:), OPTIONAL, INTENT(IN) :: ida_H
        LOGICAL :: ll_idx, ll_Rm1

        IF( PRESENT(ida_H) )THEN
            ll_idx = .TRUE.
        ELSE
            ll_idx = .FALSE.
        END IF
        IF( PRESENT(rda_Rm1) )THEN
            ll_Rm1 = .TRUE.
        ELSE
            ll_Rm1 = .FALSE.
        END IF
        CALL reset_dov(td_dov, SIZE(rda_data), ll_idx, ll_Rm1)
        td_dov%dat = rda_data
        IF(ll_idx) td_dov%H = ida_H
        IF(ll_Rm1) td_dov%Rm1 = rda_Rm1
    END SUBROUTINE assign_dov_vector

    !> \brief nullify pointer component of dynamic observation vector
    !! @param [in, out] td_dov observation vector to be nullify
    !<
    SUBROUTINE dov_nullify(td_dov)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov

        td_dov%dat => NULL()
        td_dov%Rm1 => NULL()
        td_dov%Rm1_allocated = .FALSE.
        td_dov%H => NULL()
        td_dov%H_allocated = .FALSE.
    END SUBROUTINE dov_nullify

    !< \brief Sample observation from the state variable
    !! @param [in] rda_x state variable
    !! @param [in] ida_idx indices of values to extract
    !! @param [in, out] rda_yo observation variable
    !! @param [in] error additive error (same shape as rda_yo) to be added to the observation
    !>
    SUBROUTINE sample_1D_direct_obs(rda_x, ida_idx, rda_yo, error)
        REAL(KIND=cp), DIMENSION(:), INTENT(IN)           :: rda_x
        INTEGER, DIMENSION(:), INTENT(IN)                 :: ida_idx
        REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)        :: rda_yo
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: error
        INTEGER :: ibi

        DO ibi = 1, SIZE(rda_yo)
            rda_yo(ibi) = rda_x( ida_idx(ibi) )
        END DO
        IF (PRESENT(error)) rda_yo = rda_yo + error
    END SUBROUTINE sample_1D_direct_obs


    !< \brief Sample observation from the state variable
    !! @param [in] rda_state state variable
    !! @param [in] ida_coord indices of values to extract
    !! @param [in, out] rda_obs observation variable
    !! @param [in] error additive error (same shape as rda_yo) to be added to the observation
    !>
    SUBROUTINE subsample_2D_to_1D(rda_state, ida_coord, rda_obs, error)
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_state
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: rda_obs
        INTEGER , DIMENSION(:,:), INTENT(IN) :: ida_coord
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: error
        INTEGER :: ibi

        DO ibi = 1, SIZE(ida_coord, 2)
            rda_obs(ibi) = rda_state( ida_coord(1, ibi), ida_coord(2, ibi) )
        END DO
        IF (PRESENT(error)) rda_obs = rda_obs + error
    END SUBROUTINE subsample_2D_to_1D


    !< \brief Subsample a 3D array into a 1D array
    !! @param [in] rda_state state variable
    !! @param [in] ida_coord indices of values to extract
    !! @param [in, out] rda_sample sampled variable
    !! @param [in] error additive error (same shape as rda_yo) to be added to the observation
    !>
    SUBROUTINE subsample_3D_to_1D(rda_state, ida_coord, rda_sample, error)
        REAL(dp), DIMENSION(:,:,:), INTENT(IN) :: rda_state
        REAL(dp), DIMENSION(:), INTENT(INOUT) :: rda_sample
        INTEGER , DIMENSION(:,:), INTENT(IN) :: ida_coord
        REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: error
        INTEGER :: ibi

        DO ibi = 1, SIZE(ida_coord, 2)
            rda_sample(ibi) = rda_state( ida_coord(1, ibi), ida_coord(2, ibi), ida_coord(3, ibi) )
        END DO
        IF (PRESENT(error)) rda_sample = rda_sample + error
    END SUBROUTINE subsample_3D_to_1D


    !< @brief upsample a one dimensional array into a two dimensional array
    !! @param [in] rda_1D one dimensional array to upsample
    !! @param [in] ida_coord indexes of values of the 1d array into the 2d array
    !! @param [in, out] rda_2D two dimensional array to fill
    !! @param [in] fv fill value for locations that are undefined
    !>
    subroutine upSample_1D_to_2D(rda_1D, ida_coord, rda_2D, fv)
        real(dp), dimension(:), intent(in) :: rda_1D
        integer , dimension(:,:), intent(in) :: ida_coord
        real(dp), dimension(:,:), intent(out) :: rda_2D
        real(kind=sp), intent(in) :: fv
        !local variables
        integer :: ibi

        rda_2D = real( fv, dp )
        do ibi = 1, size(ida_coord, 2)
            rda_2D( ida_coord(1, ibi), ida_coord(2, ibi) ) = rda_1D(ibi)
        end do
    end subroutine upSample_1D_to_2D


    !< @brief upsample a one dimensional array into a three dimensional array
    !! @param [in] rda_1D one dimensional array to upsample
    !! @param [in] ida_coord indexes of values of the 1d array into the 2d array
    !! @param [in, out] rda_3D three dimensional array to fill
    !! @param [in] fv fill value for locations that are undefined
    !>
    subroutine upSample_1D_to_3D(rda_1D, ida_coord, rda_3D, fv)
        real(dp), dimension(:), intent(in) :: rda_1D
        integer , dimension(:,:), intent(in) :: ida_coord
        real(dp), dimension(:,:,:), intent(out) :: rda_3D
        real(kind=sp), intent(in) :: fv
        !local variables
        integer :: ibi

        rda_3D = real( fv, dp )
        do ibi = 1, size(ida_coord, 2)
        rda_3D( ida_coord(1, ibi), ida_coord(2, ibi), ida_coord(3, ibi) ) = rda_1D(ibi)
        end do
    end subroutine upSample_1D_to_3D

    !>@brief transform a map from vector to 2-D array into
    !!  a map from 2-D array to vector
    !! @param [in] coord 2-rows array of coordinate in the 2D array
    !! @param [out] vIdx 2D array giving the location in the vector
    !! @details this is usefull to extract and restaure some elements
    !!  of a given data array.
    !! given a 2D array of data with masked location, that need to be vectorized
    !! and restaured for some operations. In the vectorization operation,
    !! if the unmasked (usefull) element data[i,j] is store in the vectorized
    !! data at index l: then:
    !!  - coord[1,l] = i and coord[2,l] = j
    !!  - vIdx[i,j] = l
    !<
    subroutine coord2loc_2D(coord, vIdx )
        integer, dimension(:,:), intent(in) :: coord
        integer, dimension(:,:), intent(in out) :: vIdx
        !local variables
        integer :: l

        vIdx = -huge(l) !bad value to raise error in case of problems
        do l = 1, size(coord,2)
            vIdx( coord(1,l), coord(2,l) ) = l
        end do
    end subroutine coord2loc_2D

    !>@brief transform a map from vector to 3-D array into
    !!  a map from 3-D array to vector
    !! @param [in] coord 3-rows array of coordinate in the 3D array
    !! @param [out] vIdx 3D array giving the location in the vector
    !! @details this is usefull to extract and restaure some elements
    !!  of a given data array.
    !! Given a 3D array of data with masked location, that needs to be vectorized
    !! and restaured for some operations. In the vectorization operation,
    !! if the unmasked (usefull) element data[i,j,k] is stored in the vectorized
    !! data at index l: then:
    !!  - coord[1,l] = i, coord[2,l] = j and coord[3,l] = k
    !!  - vIdx[i,j,k] = l
    !<
    subroutine coord2loc_3D(coord, vIdx)
        integer, dimension(:,:), intent(in) :: coord
        integer, dimension(:,:,:), intent(in out) :: vIdx
        !local variables
        integer :: l

        vIdx = - huge(l) !bad value to raise error in case of problems
        do l = 1, size(coord,2)
            vIdx( coord(1,l),coord(2,l),coord(3,l) ) = l
        end do
    end subroutine coord2loc_3D

    !>resize the observation vector structure
    !! @param [in, out] td_dov observation vector to be resized
    !! @param [in] id_size new size of the observation vector
    !! @param [in] ld_Rm1 says if the field R (covariance) need to be allocated, default is .FALSE.
    !! @param [in] ld_H says if the field H (projection obs operator) need to be allocated, default is .FALSE.
    !<
    SUBROUTINE reset_dov(td_dov, id_size, ld_Rm1, ld_H)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov
        INTEGER, INTENT(IN) :: id_size
        LOGICAL, INTENT(IN), OPTIONAL :: ld_Rm1, ld_H
        LOGICAL :: ll_Rm1, ll_H

        IF ( PRESENT(ld_Rm1) ) THEN
            ll_Rm1 = ld_Rm1
        ELSE
            ll_Rm1 = .FALSE.
        END IF
        IF ( PRESENT(ld_H) ) THEN
            ll_H = ld_H
        ELSE
            ll_H = .FALSE.
        END IF

        !dat field
        CALL reset_dov_field(td_dov, id_size, OBS_FIELD)
        !R field
        IF(ll_Rm1)THEN
            CALL reset_dov_field(td_dov, id_size, OBS_CMAT)
        ELSE
            CALL reset_dov_field(td_dov, 0, OBS_CMAT)!resize to zero
        END IF
        !H field
        IF(ll_H)THEN
            CALL reset_dov_field(td_dov, id_size, OBS_XIDX)
        ELSE
            CALL reset_dov_field(td_dov, 0, OBS_XIDX)!resize to zero
        END IF
    END SUBROUTINE reset_dov

        !>resize a field of the observation vector structure
        !! @param [in, out] td_dov observation vector to be resized
    !! @param [in] id_size new size of the observation vector
    !! @param [in] id_field identificator for the field to be resized. accepted values : OBS_FIELD (for data field), OBS_CMAT (for covariance matrix), OBS_XIDX (for x indices)
    !!
    !<
    SUBROUTINE reset_dov_field(td_dov, id_size, id_field)
        TYPE(dyn_oVector), INTENT(INOUT) :: td_dov
        INTEGER, INTENT(IN) :: id_size, id_field
        INTEGER :: il_cSize

        !dat field
        IF( ASSOCIATED(td_dov%dat) )THEN
            il_cSize = SIZE(td_dov%dat)
        ELSE
            il_cSize = 0
        END IF

        SELECT CASE (id_field)
            CASE (OBS_FIELD)
                IF (il_cSize /= id_size) THEN !current size different from required size
                    IF(il_cSize > 0) THEN !current size greater than zero
                        DEALLOCATE(td_dov%dat)
                        NULLIFY(td_dov%dat)
                    END IF
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%dat(id_size) )
                        td_dov%dat = 0.0_cp
                    END IF
                END IF
            CASE (OBS_CMAT)
                IF(td_dov%Rm1_allocated)THEN!R has been allocated
                    IF (il_cSize /= id_size) THEN!current size different from required size
                        IF(il_cSize > 0) THEN!current size greater than zero
                            DEALLOCATE(td_dov%Rm1)
                            NULLIFY(td_dov%Rm1)
                            td_dov%Rm1_allocated = .FALSE.
                        END IF
                        IF(id_size>0)THEN!required size greater than zero
                            ALLOCATE( td_dov%Rm1(id_size) )
                            td_dov%Rm1 = 0.0_cp
                            td_dov%Rm1_allocated = .TRUE.
                        END IF
                    END IF
                ELSE!R has not been allocated
                    NULLIFY(td_dov%Rm1)
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%Rm1(id_size) )
                        td_dov%Rm1 = 0.0_cp
                        td_dov%Rm1_allocated = .TRUE.
                    END IF
                END IF
            CASE (OBS_XIDX)
                IF(td_dov%H_allocated)THEN!H has been allocated
                    IF (il_cSize /= id_size) THEN!current size different from required size
                        IF(il_cSize > 0) THEN!current size greater than zero
                            DEALLOCATE(td_dov%H)
                            NULLIFY(td_dov%H)
                            td_dov%H_allocated = .FALSE.
                        END IF
                        IF(id_size>0)THEN!required size greater than zero
                            ALLOCATE( td_dov%H(id_size) )
                            td_dov%H = 0
                            td_dov%H_allocated = .TRUE.
                        END IF
                    END IF
                ELSE!H has not been allocated
                    NULLIFY(td_dov%H)
                    IF(id_size>0)THEN!required size greater than zero
                        ALLOCATE( td_dov%H(id_size) )
                        td_dov%H = 0
                        td_dov%H_allocated = .TRUE.
                    END IF
                END IF
        END SELECT
    END SUBROUTINE reset_dov_field

    subroutine assign_drv_vector(td_dra, rd_data)
        type(dyn_rVector), intent(inout) :: td_dra
        real(kind=cp), dimension(:), intent(in) :: rd_data

        call reset_drv(td_dra, size(rd_data))
        td_dra%dat = rd_data
    end subroutine assign_drv_vector

    !> \brief nullify pointer component of dynamic real vector
    !! @param [in, out] td_drv dynamic real vector to be nullify
    !<
    subroutine drv_nullify(td_drv)
        type(dyn_rVector), intent(inout) :: td_drv

        nullify(td_drv%dat)
    end subroutine drv_nullify

    !> \brief Allocate space for array (pointer) of integer dynamic vectors
    !! @param tda_div array pointer to be allocated
    !! @param [in] id_size size to be allocated
    !! All the element of tda_div are nullified
    !<
    subroutine allocate_div_array(tda_div, id_size)
        type(dyn_iVector), pointer, dimension(:) :: tda_div
        integer, intent(in) :: id_size
        !local variables
        integer :: ibi

        if (id_size>0)then
            allocate( tda_div(id_size) )
            do ibi=1,id_size
                nullify( tda_div(ibi)%dat )
            end do
        else
            nullify(tda_div)
        end if
    end subroutine allocate_div_array

    !> \brief Allocate space for array (pointer) of real dynamic vectors
    !! @param tda_drv array pointer to be allocated
    !! @param [in] id_size size to be allocated
    !! All the element of tda_drv are nullified
    !<
    subroutine allocate_drv_array(tda_drv, id_size)
        type(dyn_rVector), pointer, dimension(:) :: tda_drv
        integer, intent(in) :: id_size
        !local variables
        integer :: ibi

        if (id_size>0)then
            allocate( tda_drv(id_size) )
            do ibi=1,id_size
                nullify( tda_drv(ibi)%dat )
            end do
        else
            nullify(tda_drv)
        end if
    end subroutine allocate_drv_array

    !> \brief Deallocate space used by an array (pointer) of dynamic vectors
    !! @param tda_div array pointer to be allocated
    !! if tda_div is associated, all its elements are deallocated if associated and nullified prior to the deallocation of tda_div
    !<
    subroutine deallocate_div_array(tda_div)
        type(dyn_iVector), pointer, dimension(:) :: tda_div
        !local variables
        integer :: ibi

        if(associated(tda_div))then
            do ibi=1,size(tda_div)
                call reset_div( tda_div(ibi), 0 )
            end do
            deallocate(tda_div)
        end if
        nullify(tda_div)
    end subroutine deallocate_div_array

    !> \brief Deallocate space used by an array (pointer) of dynamic vectors
    !! @param tda_drv array pointer to be allocated
    !! if tda_drv is associated, all its elements are deallocated if associated and nullified prior to the deallocation of tda_drv
    !<
    subroutine deallocate_drv_array(tda_drv)
        type(dyn_rVector), pointer, dimension(:) :: tda_drv
        !local variables
        integer :: ibi

        if(associated(tda_drv))then
            do ibi=1,size(tda_drv)
                call reset_drv( tda_drv(ibi), 0 )
            end do
            deallocate(tda_drv)
        end if
        nullify(tda_drv)
    end subroutine deallocate_drv_array

    !> @brief (re)allocate of free dynamic space
    !! @param [in, out] tda_drv vector of dynamic vectors
    !! @param [in] id_size (optional) new size of all dynamic vectors in the @a tda_drv
    !! @param [in] id_lb (optional) lower bound  of all dynamic vectors
    !!  in the @a tda_drv
    !! @details if @a id_lb is provided, it is used as
    !!   the lower bound in the allocation, otherwise, the default
    !!   lower bound of 1 is used.
    !!   if @a id_size is not present or is set to zero,
    !!   then the dynamic space is deallocated if it was previously allocated
    !! The same operation is done for all the element of @a tda_drv
    !<
    subroutine reset_drv_array(tda_drv, id_size, id_lb)
        type(dyn_rVector), dimension(:), intent(inout) :: tda_drv
        integer, intent(in), optional :: id_size
        integer, intent(in), optional :: id_lb!lower bound in the array if needed
        !local variables
        integer :: ibi

        do ibi=1,size(tda_drv)
            call reset_drv( tda_drv(ibi), id_size, id_lb )
        end do
    end subroutine reset_drv_array

    !> @brief (re)allocate of free dynamic space
    !! @param [in, out] td_dra
    !! @param [in] id_size (optional) new size
    !! @param [in] id_lb (optional) lower bound
    !! @details if @a id_lb is provided, it is used as
    !!   the lower bound in the allocation, otherwise, the default
    !!   lower bound of 1 is used.
    !!   if @a id_size is not present or is set to 0,
    !!   then the dynamic space is deallocated if it was previously allocated
    !!   if @a id_size is not present or is set to zero,
    !!   then the dynamic space is deallocated if it was previously allocated
    !<
    subroutine reset_drv(td_dra, id_size, id_lb)
        type(dyn_rVector), intent(inout) :: td_dra
        integer, intent(in), optional :: id_size
        integer, intent(in), optional :: id_lb!lower bound in the array if needed
        integer :: il_lb, il_size, il_stat
        logical :: ll_allocate


        if ( present(id_size) ) then
            il_size = id_size
        else
            il_size = 0
        end if
        if ( present(id_lb) ) then
            il_lb = id_lb
        else
            il_lb = 1
        end if

        if( associated(td_dra%dat) )then
            if( (il_size==size(td_dra%dat)).and.(il_lb == lbound(td_dra%dat, 1) ) )then
                ll_allocate = .false.
                !print*, '--------------------- 1'
            else
                deallocate(td_dra%dat)
                nullify(td_dra%dat)
                ll_allocate = .true.
                !print*, '--------------------- 2'
            end if
        else
            ll_allocate = .true.
            !print*, '--------------------- 3'
        end if
        il_stat=0
        if(ll_allocate.and.(il_size>0))allocate( td_dra%dat(il_lb:il_lb+il_size-1))!, stat = il_stat)
        if(il_stat/=0)then
            print*, 'in reset_drv, allocate fails, stat = ',il_stat
            call abort
        end if
    END SUBROUTINE reset_drv

  SUBROUTINE assign_div_vector(td_dia, id_data)
      TYPE(dyn_iVector), INTENT(INOUT) :: td_dia
      INTEGER, DIMENSION(:), INTENT(IN) :: id_data

      CALL reset_div(td_dia, SIZE(id_data))
      td_dia%dat = id_data
  END SUBROUTINE assign_div_vector

  SUBROUTINE reset_div(td_dia, id_size, id_lb)
      TYPE(dyn_iVector), INTENT(INOUT) :: td_dia
      INTEGER, INTENT(IN) :: id_size
      INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound in the array if needed
      INTEGER :: il_lb, il_stat
      LOGICAL :: ll_allocate

      IF ( PRESENT(id_lb) ) THEN
          il_lb = id_lb
      ELSE
          il_lb = 1
      END IF

      IF( ASSOCIATED(td_dia%dat) )THEN
          IF( (id_size==SIZE(td_dia%dat)).AND.(il_lb == LBOUND(td_dia%dat, 1) ) )THEN
              ll_allocate = .FALSE.
          ELSE
              DEALLOCATE(td_dia%dat)
              NULLIFY(td_dia%dat)
              ll_allocate = .TRUE.
          END IF
      ELSE
          ll_allocate = .TRUE.
      END IF
      il_stat=0
      IF(ll_allocate.AND.(id_size>0))ALLOCATE( td_dia%dat(il_lb:il_lb+id_size-1))!, STAT = il_stat)
      IF(il_stat/=0)THEN
          PRINT*, 'In reset_div, ALLOCATE fails, STAT = ',il_stat
          CALL ABORT
      END IF
  END SUBROUTINE reset_div

  SUBROUTINE assign_dsv_vector(td_dsa, ada_data)
      TYPE(dyn_sVector), INTENT(INOUT) :: td_dsa
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: ada_data
      INTEGER :: ibi, il_size

      il_size = SIZE(ada_data)
      CALL reset_dsv(td_dsa, il_size)
      DO ibi = 1, il_size
          td_dsa%dat(ibi) = ada_data(ibi)
      END DO
  END SUBROUTINE assign_dsv_vector

  FUNCTION dsv_maxLen(td_dsa) RESULT(il_max)
      TYPE(dyn_sVector), INTENT(IN) :: td_dsa
      INTEGER :: il_max, il_tmp, ibi
      il_max = LEN_TRIM( td_dsa%dat( LBOUND(td_dsa%dat,1) ) )
      DO ibi = LBOUND(td_dsa%dat,1)+1, UBOUND(td_dsa%dat,1)
          il_tmp = LEN_TRIM( td_dsa%dat(ibi) )
          IF (il_max < il_tmp) il_max = il_tmp
      END DO
  END FUNCTION dsv_maxLen

  SUBROUTINE reset_dsv(td_dsa, id_size, id_lb)
      TYPE(dyn_sVector), INTENT(INOUT) :: td_dsa
      INTEGER, INTENT(IN) :: id_size
      INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound in the array if needed
      INTEGER :: il_lb, il_stat
      LOGICAL :: ll_allocate

      IF ( PRESENT(id_lb) ) THEN
          il_lb = id_lb
      ELSE
          il_lb = 1
      END IF

      IF( ASSOCIATED(td_dsa%dat) )THEN
          IF( (id_size==SIZE(td_dsa%dat)).AND.(il_lb == LBOUND(td_dsa%dat, 1) ) )THEN
              ll_allocate = .FALSE.
          ELSE
              DEALLOCATE(td_dsa%dat)
              NULLIFY(td_dsa%dat)
              ll_allocate = .TRUE.
          END IF
      ELSE
          ll_allocate = .TRUE.
      END IF
      il_stat = 0
      IF(ll_allocate.AND.(id_size>0))ALLOCATE( td_dsa%dat(il_lb:il_lb+id_size-1))!, STAT = il_stat)
      IF(il_stat/=0)THEN
          PRINT*, 'In reset_dsv, ALLOCATE fails, STAT = ',il_stat
          CALL ABORT
      END IF
  END SUBROUTINE reset_dsv

  !< \brief assign 1D array data to a 2D dynamic vector, the second dimension of the dynamic array is resize to 1
  !! @param [in, out] td_dd 2D dynamic array to wich one wants to assign a 1D array
  !! @param [in] rd_data array data to be assigned
  !! \todo overload assignment operator
  !>
  SUBROUTINE assign_2DdrA_vector(td_dd, rda_data)
      TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
      REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data

      CALL reset_d2DrA(td_dd, 1)
      td_dd%dat(1) = rda_data
  END SUBROUTINE assign_2DdrA_vector

  !< \brief assign 2D array data to a 2D dynamic array
  !! @param [in, out] td_dd 2D dynamic array to wich one wants to assign a 2D array
  !! @param [in] rd_data array data to be assigned
  !! \todo overload assignment operator
  !! \details each element of td_dd receives a column of rda_data, recall that 2D dynamic array is a java style 2D array
  !! column major assumption
  !>
  SUBROUTINE assign_2DdrA_2DA(td_dd, rda_data)
      TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
      REAL(KIND=cp), DIMENSION(:, :), INTENT(IN) :: rda_data
      INTEGER :: ibi

      CALL reset_d2DrA(td_dd, size(rda_data, 2))
      DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
          td_dd%dat(ibi) = rda_data(:, ibi)
      END DO
  END SUBROUTINE assign_2DdrA_2DA

  !< \brief assign 1D array data to an element of a 2D dynamic array
  !! recall that, a 2D dynamic array is like java 2D array
  !! @param [in, out] td_dd 2D dynamic array to wich one wants to put a 1D vector
  !! @param [in] rd_data vector data to be assigned
  !! @param [in] id_idx index of the 1D vector in td_dd, recall that 2D dynamic array is a java style 2D array
  !! \details
  !! \remark call to this routine can causes segmentation fault
  !! if the first dimension of the td_dd is less than id_idx
  !>
  SUBROUTINE put_2DdrA_vector(td_dd, id_idx, rda_data)
      TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
      INTEGER, INTENT(IN)               :: id_idx
      REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_data

      td_dd%dat(id_idx) = rda_data
  END SUBROUTINE put_2DdrA_vector

  !< \brief get a 1D array data from a 2D dynamic array
  !! recall that, a 2D dynamic array is like java 2D array
  !! @param [in] td_dd 2D dynamic array from wich one wants to get a 1D vector
  !! @param [out] rda_data vector data to be get from td_dd
  !! @param [in] id_idx index of the 1D vector in td_dd, recall that 2D dynamic array is a java style 2D array
  !! \details
  !! \remark call to this routine can cause segmentation fault
  !! if the first dimension of the td_dd is less than id_idx
  !! also make sure that rl_data has the same size as
  !>
  SUBROUTINE get_2DdrA_vector(td_dd, id_idx, rda_data)
      TYPE(dyn_2DrArray), INTENT(IN) :: td_dd
      INTEGER, INTENT(IN)            :: id_idx
      REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_data

      rda_data = td_dd%dat(id_idx)%dat
  END SUBROUTINE get_2DdrA_vector

  SUBROUTINE reset_d2DrA(td_dd, id_size, id_lb)
      TYPE(dyn_2DrArray), INTENT(INOUT) :: td_dd
      INTEGER, INTENT(IN) :: id_size
      INTEGER, INTENT(IN), OPTIONAL :: id_lb!lower bound of the first dimension if needed
      INTEGER :: il_lb, ibi
      IF ( PRESENT(id_lb) ) THEN
          il_lb = id_lb
      ELSE
          il_lb = 1
      END IF

      IF(ASSOCIATED(td_dd%dat)) THEN
          DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
              CALL reset_drv(td_dd%dat(ibi), 0)
          END DO
          DEALLOCATE(td_dd%dat)
          NULLIFY(td_dd%dat)
      END IF
      IF(id_size>0)THEN
          ALLOCATE( td_dd%dat(il_lb:il_lb+id_size-1) )
          DO ibi = LBOUND(td_dd%dat,1), UBOUND(td_dd%dat,1)
              NULLIFY(td_dd%dat(ibi)%dat)
          END DO
      END IF
  END SUBROUTINE reset_d2DrA

  !1D dynamic array minval
  FUNCTION drv_minval(td_dd) RESULT(rl_min)
      TYPE(dyn_rVector), INTENT(IN) :: td_dd
      REAL(cp) :: rl_min

      rl_min = MINVAL(td_dd%dat)
  END FUNCTION drv_minval

  !array of 1D dynamic array minval
  FUNCTION drv_array_minval(tda_dd) RESULT(rl_min)
      TYPE(dyn_rVector), DIMENSION(:), INTENT(IN) :: tda_dd
      REAL(cp) :: rl_min, rl_temp
      INTEGER :: ibi

      rl_min = MINVAL(tda_dd(1)%dat)
      DO ibi = 2, SIZE(tda_dd)
          rl_temp = MINVAL(tda_dd(ibi)%dat)
          IF(rl_min>rl_temp) rl_min = rl_temp
      END DO
  END FUNCTION drv_array_minval

  !2D dynamic array minval
  FUNCTION d2DrA_minval(td_dd) RESULT(rl_min)
      TYPE(dyn_2DrArray), INTENT(IN) :: td_dd
      REAL(cp) :: rl_min, rl_temp
      INTEGER :: ibi, il_lb, il_ub

      il_lb = LBOUND(td_dd%dat,1)
      il_ub = UBOUND(td_dd%dat,1)
      rl_min = MINVAL(td_dd%dat( il_lb )%dat)
      DO ibi = il_lb+1, il_ub
          rl_temp = MINVAL(td_dd%dat(ibi)%dat)
          IF(rl_min>rl_temp) rl_min = rl_temp
      END DO
  END FUNCTION d2DrA_minval

    !1d dynamic array maxval
    function drv_maxval(td_dd) result(rl_max)
        type(dyn_rvector), intent(in) :: td_dd
        real(cp) :: rl_max

        rl_max = maxval(td_dd%dat)
    end function drv_maxval

    !array of 1d dynamic array maxval
    function drv_array_maxval(tda_dd) result(rl_max)
        type(dyn_rvector), dimension(:), intent(in) :: tda_dd
        real(cp) :: rl_max, rl_temp
        integer :: ibi

        rl_max = maxval(tda_dd(1)%dat)
        do ibi = 2, size(tda_dd)
            rl_temp = maxval(tda_dd(ibi)%dat)
            if(rl_max<rl_temp) rl_max = rl_temp
        end do
    end function drv_array_maxval

    !< 2d dynamic array maxval
    function d2dra_maxval(td_dd) result(rl_max)
        type(dyn_2drarray), intent(in) :: td_dd
        real(cp) :: rl_max, rl_temp
        integer :: ibi, il_lb, il_ub

        il_lb = lbound(td_dd%dat,1)
        il_ub = ubound(td_dd%dat,1)
        rl_max = maxval(td_dd%dat( il_lb )%dat)
        do ibi = il_lb+1, il_ub
            rl_temp = maxval(td_dd%dat(ibi)%dat)
            if(rl_max<rl_temp) rl_max = rl_temp
        end do
    end function d2dra_maxval

    !> \brief compute the number of regularly-spaced coordinates
    !! @param [in] ida_max maximum coordinates that can be generate in each dimension
    !! @param [in] shifts optional minimum coordinates that can be generate in each
    !!  dimension, default is 1 for each dimension. They give the lowest coordinates.
    !! @param [in] steps step between successive coordinates in each dimension
    !! @param [out] ncoords optional coordinates count in each dimension
    !! \details Asumming that dim_count is the number of dimensions of the problem,
    !! steps, ida_max [, shifts] [and ncoords]  are one-dimensional arrays of size
    !! dim_count each. This subroutine also compute the number of regularly-spaced
    !! coordinates in each dimension of the problem, see the optional argument ncoords
    !! \todo add a version for real coordinates and use the overloading, rename this
    !!  as regular_int_coord_count
    !<
    FUNCTION regular_coord_count(ida_max, shifts, steps, ncoords) RESULT(il_count)
        INTEGER, DIMENSION(:), INTENT(IN) :: ida_max
        INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)  :: shifts, steps
        INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: ncoords
        INTEGER, DIMENSION(SIZE(ida_max)) :: ila_shift, ila_step, ila_count
        INTEGER :: il_count

        IF( PRESENT(shifts) ) THEN
            ila_shift = shifts
        ELSE
            ila_shift = 1
        END IF
        IF( PRESENT(steps) ) THEN
            ila_step = steps
        ELSE
            ila_step = 1
        END IF
        ila_count = 1 + (ida_max - ila_shift)/ila_step
        il_count = PRODUCT(ila_count)
        IF( PRESENT(ncoords) ) ncoords = ila_count
    END FUNCTION regular_coord_count

  !> \brief build regularly-spaced coordinates of type integer
  !! @param [out] ida_coord regularly-spaced coordinates generated by the subroutine
  !! @param [in] ncoords optional array giving the number of coordinates in each dimension
  !! @param [in] max_coords optional array giving the maximum coordinate in each dimension
  !! @param [in] shifts optional array giving the shift from zero in each dimension. Default is 1 to conform to one-base indexation in fortran.
  !! @param [in] steps optional array giving the step between successive coordinates in each dimension. Default is 1.
  !! @details At least one of ncoords and max_coords must be provided. If ncoords
  !!   is provided, max_coords is not necessary. Asumming that coord_count is the
  !!   total number of coordinates to generate and dim_count is the number of
  !!   dimensions of the problem, ida_coord is of shape (dim_count, coord_count).
  !!   The second dimension can be larger than coord_count but not less. ida_step,
  !!   ida_max and ida_min (if provided) are one-dimensional array of size dim_count.
  !! This routine is used to build regularly-spaced coordinates for twin observation.
  !! The coordinates are 1-based indexed for fortran default indexation. The integer
  !! coordinates can easilly be converted to real-value coordinates if uniformly-
  !! spaced mesh is used. Shift to 0-based index and multiply by the space step
  !! (each dimension separately if the space step if different for each), then add
  !! (the real-coordinates of the origin if not 0.
  !! @todo add some clarification to the shift
  !<
  SUBROUTINE build_regular_int_coord(ida_coord, ncoords, max_coords, shifts, steps)
    INTEGER, DIMENSION(:,:), INTENT(INOUT) :: ida_coord
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)  :: max_coords, shifts, steps, ncoords
    INTEGER, DIMENSION(SIZE(ida_coord,1)) :: ila_shift, ila_step, ila_ncoord, ila_prodinf
    INTEGER :: il_ncoord, ibi, il_idx, il_ndim, il_val

    IF( PRESENT(shifts) ) THEN
      ila_shift = shifts
    ELSE
      ila_shift = 1
    END IF
    IF( PRESENT(steps) ) THEN
      ila_step = steps
    ELSE
      ila_step = 1
    END IF
    IF( PRESENT(ncoords) ) THEN
      ila_ncoord = ncoords
    ELSE IF ( PRESENT(max_coords) )THEN
      !PRINT*, 'Calling regular_coord_count'
      il_ncoord = regular_coord_count(max_coords, ila_shift, ila_step, ila_ncoord)
    END IF
    !PRINT*, 'checking PRESENT(ncoords).OR.PRESENT(max_coords)'
    IF( PRESENT(ncoords).OR.PRESENT(max_coords) ) THEN
      il_ndim = SIZE(ida_coord,1)
      il_ncoord = PRODUCT(ila_ncoord)
      !the index i of  variable ila_prodinf stores the product of the (i-1) first elts of il_ncoord, ila_obs actually give the shape of the non vectorized obs array
      ila_prodinf(1) = 1
      DO ibi = 2,il_ndim
        ila_prodinf(ibi) = ila_prodinf(ibi-1)*ila_ncoord(ibi-1)
      END DO
      !PRINT*, 'starting loop'
      !PRINT*, 'ila_ncoord   = ', ila_ncoord
      !PRINT*, 'ila_step     = ', ila_step
      !PRINT*, 'ila_shift    = ', ila_shift
      !PRINT*, 'ila_prodinf  = ', ila_prodinf
      DO il_idx = 1,il_ncoord
        !PRINT*, 'il_idx = ', il_idx
        il_val = il_idx-1 !compute zero-based coordinates
        DO ibi = il_ndim,2,-1

          !shift to one-based coordinated with appropriate shift
          ida_coord(ibi, il_idx) = il_val/ila_prodinf(ibi)*ila_step(ibi) + ila_shift(ibi)
          il_val = mod( il_val, ila_prodinf(ibi) )
        END DO
        ida_coord(1, il_idx) = il_val*ila_step(ibi) + ila_shift(ibi)
      END DO
    ELSE
      ida_coord = -999
    END IF
  END SUBROUTINE build_regular_int_coord

    !> \brief Builds regularly spaced indices
    !! @param [in, out] ida_idx array of builded indices
    !! @param [in] id_idxMin minimal value for idx
    !! @param [in] id_idxMax maximal value for idx
    !! @param [in] id_obsDist distance between two observations
    !! \remark: this subroutine was previously named build_regular_idx
    !>
    SUBROUTINE build_regular_idx(ida_idx, id_idxMin, id_idxMax, id_obsDist)
        INTEGER, DIMENSION(:), INTENT(INOUT) :: ida_idx
        INTEGER, INTENT(IN) :: id_idxMin, id_idxMax, id_obsDist
        INTEGER :: ibi, ibj

        ibj = 1
        DO ibi = id_idxMin, id_idxMax, id_obsDist
            ida_idx(ibj) = ibi
            ibj = ibj + 1
        END DO
    END SUBROUTINE build_regular_idx

  !=======version 2

    !> @brief builds mask for regularly spaced coordinates in 1D
    !! @param [out] coord_mask mask array for regularly spaced coordinates, the
    !!   indexes of unmasked positions give the regularly-spaced coord
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound
    !!   default is 0.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound
    !!   default is 0.
    !! @param [in] mask array marking the location of masked points
    !<
    subroutine build_regular_icoord_mask_1D(coord_mask, steps, lskips, uskips, mask)
        integer, dimension(:), intent(in out) :: coord_mask
        integer, intent(in) :: steps
        integer, optional, intent(in)    :: lskips, uskips
        integer, dimension(:), optional, intent(in) :: mask
        !local variables
        integer  :: il_lskips, il_uskips!, ila_count

        if( present(lskips) ) then
            il_lskips = lskips
        else
            il_lskips = 0
        end if

        if( present(uskips) ) then
            il_uskips = uskips
        else
            il_uskips = 0
        end if
        !consider all points masked
        coord_mask = MASKED
        !unmask points at sampling locations
        coord_mask( 1+il_lskips:size(coord_mask)-il_uskips:steps) =  UNMASKED
        !mask locations that are masked by input mask
        if( present(mask) )then
            if( count( shape(mask)/=shape(coord_mask) )>0 )then
                call debug("coord_mask and mask must have the same shape",tag=dALLWAYS)
                call stop_program("build_regular_icoord_mask 1D: Non conforming shapes")
            end if

            where( (coord_mask==UNMASKED).and.(mask/=UNMASKED) )
                coord_mask = MASKED
            end where
        end if
    end subroutine build_regular_icoord_mask_1D

    !> @brief builds mask for regularly spaced coordinates in 2D
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @param [out] coord_mask mask array for regularly spaced coordinates, the
    !!   indexes of unmasked positions give the regularly-spaced coord
    !! @details
    !!  steps, ida_max [, lskips, uskips] are vectors of 2 elements.
    !! mask and coord_mask are 2-dimentional arrays with named constants
    !! MASKED at masked locations and UNMASKED at unmasked location
    !<
    subroutine build_regular_icoord_mask_2D(coord_mask, steps, lskips, uskips, mask)
        integer, dimension(:,:), intent(in out) :: coord_mask
        integer, dimension(2), intent(in) :: steps
        integer, dimension(2), optional, intent(in)    :: lskips, uskips
        integer, dimension(:,:), optional, intent(in) :: mask
        !local variables
        integer, dimension(2)  :: ila_lskips, ila_uskips!, ila_count

        if( present(lskips) ) then
            ila_lskips = lskips
        else
            ila_lskips = 0
        end if

        if( present(uskips) ) then
            ila_uskips = uskips
        else
            ila_uskips = 0
        end if
        !consider all points masked
        coord_mask = MASKED
        !unmask points at sampling locations
        coord_mask( 1+ila_lskips(1):size(coord_mask,1)-ila_uskips(1):steps(1)&
                ,1+ila_lskips(2):size(coord_mask,2)-ila_uskips(2):steps(2)&
                ) =  UNMASKED
        !mask locations that are masked by input mask
        if( present(mask) )then
            if( count( shape(mask)/=shape(coord_mask) )>0 )then
                call debug("coord_mask and mask must have the same shape",tag=dALLWAYS)
                call stop_program("build_regular_icoord_mask_2D: Non conforming shapes")
            end if

            where( (coord_mask==UNMASKED).and.(mask/=UNMASKED) )
                coord_mask = MASKED
            end where
        end if
    end subroutine build_regular_icoord_mask_2D

    !> @brief builds mask for regularly spaced coordinates in 3D
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @param [out] coord_mask mask array for regularly spaced coordinates, the
    !!   indexes of unmasked positions give the regularly-spaced coord
    !! @details
    !!  steps, ida_max [, lskips, uskips] are vectors of 3 elements.
    !! mask and coord_mask are 3-dimentional arrays with named constants
    !! MASKED at masked locations and UNMASKED at unmasked location
    !<
    subroutine build_regular_icoord_mask_3D(coord_mask, steps, lskips, uskips, mask)
        integer, dimension(:,:,:), intent(in out) :: coord_mask
        integer, dimension(3), intent(in) :: steps
        integer, dimension(3), optional, intent(in)    :: lskips, uskips
        integer, dimension(:,:,:), optional, intent(in) :: mask
        !local variables
        integer, dimension(3)  :: ila_lskips, ila_uskips!, ila_count

        if( present(lskips) ) then
            ila_lskips = lskips
        else
            ila_lskips = 0
        end if

        if( present(uskips) ) then
            ila_uskips = uskips
        else
            ila_uskips = 0
        end if
        !consider all points masked
        coord_mask = MASKED
        !unmask points at sampling locations
        coord_mask( 1+ila_lskips(1):size(coord_mask,1)-ila_uskips(1):steps(1)&
                ,1+ila_lskips(2):size(coord_mask,2)-ila_uskips(2):steps(2)&
                ,1+ila_lskips(3):size(coord_mask,3)-ila_uskips(3):steps(3)&
                ) =  UNMASKED
        !mask locations that are masked by input mask
        if( present(mask) )then
            if( count( shape(mask)/=shape(coord_mask) )>0 )then
                call debug("coord_mask and mask must have the same shape",tag=dALLWAYS)
                call stop_program("build_regular_icoord_mask: Non conforming shapes for 3D mask")
            end if

            where( (coord_mask==UNMASKED).and.(mask/=UNMASKED) )
                coord_mask = MASKED
            end where
        end if
    end subroutine build_regular_icoord_mask_3D

    !> \brief computes the number of regularly-spaced coordinates in 1D
    !! @param [in] nx total number of points
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound
    !!   default is 0.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound
    !!   default is 0.
    !! @param [in] mask array marking the location of masked points
    !<
    function regular_icoord_count_1D(nx, steps, lskips, uskips, mask) result(il_count)
        integer, intent(in) :: nx
        integer, intent(in) :: steps
        integer, optional, intent(in)    :: lskips, uskips
        integer, dimension(nx), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx) :: coord_mask
        integer :: il_count

        call build_regular_icoord_mask( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
    end function regular_icoord_count_1D

    !> \brief computes the number of regularly-spaced coordinates in 2D
    !! @param [in] nx total number of points in the first (x) dimension
    !! @param [in] ny total number of points in the second (y) dimension
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound
    !!   in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @details
    !!  steps, ida_max [, lskips, uskips] are vectors of 2 elements.
    !! mask and coord_mask are 2-dimentional arrays with named constants
    !! MASKED at masked locations and UNMASKED at unmasked location
    !<
    function regular_icoord_count_2D(nx, ny, steps, lskips, uskips, mask) result(il_count)
        integer, intent(in) :: nx, ny
        integer, dimension(2), intent(in) :: steps
        integer, dimension(2), optional, intent(in)    :: lskips, uskips
        integer, dimension(nx, ny), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx, ny) :: coord_mask
        integer :: il_count

        call build_regular_icoord_mask( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
    end function regular_icoord_count_2D

    !> \brief computes the number of regularly-spaced coordinates in 3D
    !! @param [in] nx total number of points in the first (x) dimension
    !! @param [in] ny total number of points in the second (y) dimension
    !! @param [in] nz total number of points in the third (z) dimension
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the lower bound in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip in the upper bound in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @details
    !!  steps, ida_max [, lskips, uskips] are vectors of 3 elements.
    !! mask and coord_mask are 3-dimentional arrays with named constants
    !! MASKED at masked locations and UNMASKED at unmasked location
    !<
    function regular_icoord_count_3D(nx, ny, nz, steps, lskips, uskips, mask) result(il_count)
        integer, intent(in) :: nx, ny, nz
        integer, dimension(3), intent(in) :: steps
        integer, dimension(3), optional, intent(in)    :: lskips, uskips
        integer, dimension(nx, ny, nz), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx, ny, nz) :: coord_mask
        integer :: il_count

        call build_regular_icoord_mask( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
    end function regular_icoord_count_3D

    !> @brief computes regularly-spaced coordinates in 1D
    !! @param [in] nx total number of points in the first (x) dimension
    !! @param [out] icoords regularly-spaced coordinates generated by the subroutine
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip
    !!  in the lower bound in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip
    !!  in the upper bound in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @details Asumming that dim_count is the number of dimensions
    !!  of the problem, steps, ida_max [, lskips, uskips] are
    !!  one-dimensional arrays of size dim_count each.
    !! mask is a dim_count-dimentional array with @a MASKED at masked locations
    !! and @a UNMASKED at unmasked location one can call regular_icoord_count2D
    !! to computes the appropriate size of icoords
    !! @note
    !! This subroutine expects icoords to have the right shape
    !! @see regular_icoord_count for the determination of the number of coord.
    !!
    !<
    subroutine build_regular_icoord_1D(nx, icoords, steps, lskips, uskips, mask)
        integer, intent(in) :: nx
        integer, dimension(:,:), intent(out) :: icoords
        integer, intent(in) :: steps
        integer, optional, intent(in)    :: lskips, uskips
        integer, dimension(nx), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx) :: coord_mask
        integer :: il_count, pos, row

        call build_regular_icoord_mask_1D( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
        pos = 1
        do row = 1, nx
            if(coord_mask(row)==UNMASKED)then
            icoords(1, pos) = row
            pos = pos+1
            end if
        end do
    end subroutine build_regular_icoord_1D

    !> @brief computes regularly-spaced coordinates in 2D
    !! @param [in] nx total number of points in the first (x) dimension
    !! @param [in] ny total number of points in the second (y) dimension
    !! @param [out] icoords regularly-spaced coordinates generated by the subroutine
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip
    !!  in the lower bound in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip
    !!  in the upper bound in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @details Asumming that dim_count is the number of dimensions
    !!  of the problem, steps, ida_max [, lskips, uskips] are
    !!  one-dimensional arrays of size dim_count each.
    !! mask is a dim_count-dimentional array with @a MASKED at masked locations
    !! and @a UNMASKED at unmasked location one can call regular_icoord_count2D
    !! to computes the appropriate size of icoords
    !<
    subroutine build_regular_icoord_2D(nx, ny, icoords, steps, lskips, uskips, mask)
        integer, intent(in) :: nx, ny
        integer, dimension(:,:), intent(out) :: icoords
        integer, dimension(2), intent(in) :: steps
        integer, dimension(2), optional, intent(in)    :: lskips, uskips
        integer, dimension(nx, ny), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx, ny) :: coord_mask
        integer :: il_count, pos, row, col

        call build_regular_icoord_mask_2D( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
        pos = 1
        do col = 1, ny
            do row = 1, nx
                if(coord_mask(row, col)==UNMASKED)then
                icoords(1:2, pos) = [row, col]
                pos = pos+1
                end if
            end do
        end do
    end subroutine build_regular_icoord_2D

    !> \brief computes regularly-spaced coordinates in 3D
    !! @param [in] nx total number of points in the first (x) dimension
    !! @param [in] ny total number of points in the second (y) dimension
    !! @param [in] nz total number of points in the third (z) dimension
    !! @param [out] icoords regularly-spaced coordinates generated by
    !!   the subroutine
    !! @param [in] steps step between successive coordinates
    !! @param [in] lskips optional number of coordinates to skip in the
    !!    lower bound in each dimension, default is 0 for each dimension.
    !! @param [in] uskips optional number of coordinates to skip in the
    !!    upper bound in each dimension, default is 0 for each dimension.
    !! @param [in] mask array marking the location of masked points
    !! @details Asumming that dim_count is the number of dimensions
    !!  of the problem, steps, ida_max [, lskips, uskips] are
    !!  one-dimensional arrays of size dim_count each.
    !! mask is a dim_count-dimentional array with @a MASKED at masked locations
    !! and @a UNMASKED at unmasked location one can call regular_icoord_count2D
    !! to computes the appropriate size of icoords
    !<
    subroutine build_regular_icoord_3D(nx, ny, nz, icoords, steps, lskips, uskips, mask)
        integer, intent(in) :: nx, ny, nz
        integer, dimension(:,:), intent(out) :: icoords
        integer, dimension(3), intent(in) :: steps
        integer, dimension(3), optional, intent(in)    :: lskips, uskips
        integer, dimension(nx, ny, nz), optional, intent(in) :: mask
        !local variables
        integer, dimension(nx, ny, nz) :: coord_mask
        integer :: il_count, pos, i, j, k!, row, col

        call build_regular_icoord_mask_3D( coord_mask, steps, lskips, uskips, mask )
        !count the unmasked points
        il_count = count( coord_mask == UNMASKED)
        pos = 1
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if(coord_mask(i, j, k)==UNMASKED)then
                    icoords(1:3, pos) = [i, j, k]
                    pos = pos+1
                    end if
                end do
            end do
        end do
    end subroutine build_regular_icoord_3D

  !=======end of version 2

    !> \brief check if a value exist in a vector and return his index
    !! @param [in] ida_vec vector to consider
    !! @param [in] id_val value to search for
    !! @param [out] id_idx index of the last occurence of the element
    !!  in the vector, -1 if the element does not exist in the vector
    !!
    !<
    function find(ida_vec, id_val, id_idx) result(ll_find)
        integer, dimension(:), intent(in) :: ida_vec
        integer, intent(in)  :: id_val!index of the time step of interest
        integer, intent(out) :: id_idx! index of the particular obs relative to the numbers of obs
        logical              :: ll_find
        integer              :: ibi

        id_idx  = -1 !
        ll_find = .false.
        do ibi = size(ida_vec),1,-1
            if( id_val == ida_vec(ibi)) then
            !print*, 'in obsexist : id_val =', id_val
            id_idx = ibi
            ll_find = .true.
            !print*, 'in obsexist ', id_numobs, 'at time ', id_n
            end if
        end do
    end function find

  !> @brief Save 1D trajectory in ascii file
  !! @param [in] tda_trj the trajectory to be saved
  !! @param [in] ada_fileName name of the file in wich date should be saved
  !! @param [in] id_lb (optional) lower bound of the interresting part
  !!  of data to be saved at each time step
  !! @param [in] id_ub (optional) upper bound of the interresting part
  !! of data to be saved at each time step
  !! @param [in] interpolate (optional) says if interpolation is needed.
  !! This is the case if one need all data to be saved at the same location
  !! @param [in] ld_column (optional) says if each element of the trajectory
  !! should be considered as a column(usefull for gnuplot)
  !! @details Data are saved time step by time step : old the date
  !! from the first time step, then the second time step and so on.
  !! parameters id_lb, id_ub are useful in the case where one
  !! do not want to save everything, for exemple, the extension of
  !! the state variable for boundary condition
  !! Interpolation suppose that the variable is computed between
  !! grid points, the goal is to interpolate to grid point before saving.
  !! In this case, id_lb, id_ub are the bound of the computed variable.
  !<
  SUBROUTINE save_trj(tda_trj, ada_fileName, id_lb, id_ub, interpolate, ld_column)
    INTEGER, PARAMETER :: ip_fid = 41
    !> \brief format to print real
    !<
    CHARACTER(LEN=*), PARAMETER :: WFORMAT = "(ES18.6E2)"
    TYPE(dyn_rVector), DIMENSION(:), INTENT(IN) :: tda_trj
    CHARACTER(LEN=*), INTENT(IN) :: ada_fileName
    INTEGER, OPTIONAL, INTENT(IN) :: id_lb, id_ub
    LOGICAL, OPTIONAL, INTENT(IN) :: interpolate
    LOGICAL, OPTIONAL, INTENT(IN) :: ld_column
    !local variables
    REAL(KIND=cp), DIMENSION(:), ALLOCATABLE :: rla_tmp
    INTEGER :: ibk, ibj, il_lb, il_ub, il_ios, il_nbcol
    LOGICAL :: ll_interpolate, ll_column
    CHARACTER(LEN=80)   :: array_format

    !Checking optional argument
    IF( PRESENT(id_lb) ) THEN
        il_lb = id_lb
    ELSE
        il_lb = LBOUND( tda_trj(1)%dat, 1)
    END IF
    IF( PRESENT(id_ub) ) THEN
        il_ub = id_ub
    ELSE
        il_ub = UBOUND( tda_trj(1)%dat, 1)
    END IF
    IF( PRESENT(interpolate) ) THEN
        ll_interpolate = interpolate
    ELSE
        ll_interpolate = .FALSE.
    END IF
    IF( PRESENT(ld_column) ) THEN
        ll_column = ld_column
    ELSE
        ll_column = .FALSE.
    END IF

    OPEN( UNIT = ip_fid, FILE = trim(ada_fileName), STATUS = FILE_REPLACE, FORM = FILE_FORMATTED, IOSTAT = il_ios)
    IF(il_ios/=0) CALL stop_program('In save_uTrj :: error creating the file '//ada_fileName)

    IF (ll_column) THEN
      il_nbcol = SIZE(tda_trj)
      ALLOCATE( rla_tmp( il_nbcol ) )
      array_format = "("//TRIM( NUM2STR(il_nbcol) )//WFORMAT//")"
      !save the first row
      DO ibj =1, il_nbcol
        rla_tmp(ibj) = tda_trj(ibj)%dat(il_lb)
      END DO
      WRITE(UNIT = ip_fid, FMT=*) rla_tmp
      !save intermediate rows
      DO ibk = il_lb+1, il_ub
        DO ibj =1, il_nbcol
          IF( ll_interpolate ) THEN
            rla_tmp(ibj) = (tda_trj(ibj)%dat(ibk-1) + tda_trj(ibj)%dat(ibk) ) / 2.0_cp !interpolation
          ELSE
            rla_tmp(ibj) = tda_trj(ibj)%dat(ibk)
          END IF
        END DO
        WRITE(UNIT = ip_fid, FMT=*) rla_tmp
      END DO
      !save the last row
      IF( ll_interpolate ) THEN
        DO ibj =1, il_nbcol
          rla_tmp(ibj) = tda_trj(ibj)%dat(il_ub)
        END DO
        WRITE(UNIT = ip_fid, FMT=*) rla_tmp
      END IF
    ELSE
      IF( ll_interpolate ) THEN
        ALLOCATE( rla_tmp(il_lb:il_ub+1) )
      ELSE
        ALLOCATE( rla_tmp(il_lb:il_ub) )
      END IF
      il_nbcol = SIZE(rla_tmp)
      array_format = "("//TRIM( NUM2STR(il_nbcol) )//WFORMAT//")"

      DO ibk = 1, SIZE(tda_trj)
        IF (ll_interpolate) THEN ! interpolation
            rla_tmp(il_lb) = tda_trj(ibk)%dat(il_lb)
            rla_tmp(il_lb+1:il_ub) = ( tda_trj(ibk)%dat(il_lb:il_ub-1) + tda_trj(ibk)%dat(il_lb+1:il_ub) ) / 2.0_cp
            rla_tmp(il_ub+1) = tda_trj(ibk)%dat(il_ub)
        ELSE
            rla_tmp = tda_trj(ibk)%dat(il_lb:il_ub)
        END IF
        WRITE(UNIT = ip_fid, FMT=*) rla_tmp
      END DO
    END IF

    CLOSE(ip_fid)
    !PRINT*, 'SIZE(rla_tmp) = ', SIZE(rla_tmp)
    DEALLOCATE(rla_tmp)
  END SUBROUTINE save_trj


    !to be replaced by reshape
    subroutine bloc2vec(rda_f,id_i1, id_j1, id_isize, id_jsize, rda_x)
        !column major storage assumption
        !
        !** vectorize le bloc de la matrice rda_f(id_i1 : id_i1 + id_isize -1, id_j1 : id_j1 + id_jsize - 1) dans rda_x(1 : id_isize * id_jsize)
        !   column wise
        integer, intent(in) :: id_i1, id_j1, id_isize, id_jsize
        real(kind=dp), dimension(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1), intent(in) :: rda_f
        real(kind=dp), dimension(1 : id_isize * id_jsize), intent(out) :: rda_x
        integer ibi, ibj, il_idx

        il_idx = 0
        do ibj = id_j1, id_j1 + id_jsize - 1
            do ibi = id_i1, id_i1 + id_isize - 1
            il_idx = il_idx + 1
            rda_x(il_idx) = rda_f(ibi,ibj)
            enddo
        enddo

    end subroutine bloc2vec

  !to be replaced by RESHAPE
  SUBROUTINE vec2bloc(rda_x,id_i1, id_j1, id_isize, id_jsize, rda_f)
    !Column major storage assumption
    !
    !** Distribue le vecteur rda_x dans le bloc de la matrice
    !   rda_f(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1)
    !   /!\ PAR LIGNES
    INTEGER, INTENT(IN) :: id_i1, id_j1, id_isize, id_jsize
    REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1), INTENT(OUT) :: rda_f
    REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize), INTENT(IN) :: rda_x
    INTEGER ibi, ibj, il_idx

    il_idx = 0
    DO ibj = id_j1, id_j1 + id_jsize - 1
        DO ibi = id_i1, id_i1 + id_isize - 1
          il_idx = il_idx + 1
          rda_f(ibi,ibj) = rda_x(il_idx)
        ENDDO
    ENDDO

  END SUBROUTINE vec2bloc

  !to be replaced by RESHAPE
  SUBROUTINE bloc2vec2(rda_f,id_i1, id_j1, id_isize, id_jsize, rda_x)
    !Row major storage assumption for images processing
    !
    !** Vectorize le bloc de la matrice rda_f(id_i1 : id_i1 + id_isize -1, id_j1 : id_j1 + id_jsize - 1) dans rda_x(1 : id_isize * id_jsize)
    !   PAR LIGNES
    INTEGER, INTENT(IN) :: id_i1, id_j1, id_isize, id_jsize
    REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1), INTENT(IN) :: rda_f
    REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize), INTENT(OUT) :: rda_x
    INTEGER ibi, ibj, il_idx

    il_idx = 0
    DO ibi = id_i1, id_i1 + id_isize - 1
        DO ibj = id_j1, id_j1 + id_jsize - 1
          il_idx = il_idx + 1
          rda_x(il_idx) = rda_f(ibi,ibj)
        ENDDO
    ENDDO

  END SUBROUTINE bloc2vec2

  !to be replaced by RESHAPE
  SUBROUTINE vec2bloc2(rda_x,id_i1, id_j1, id_isize, id_jsize, rda_f)
    !Row major storage assumption for images processing
    !
    !** Distribue le vecteur rda_x dans le bloc de la matrice
    !   rda_f(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1)
    !   /!\ PAR LIGNES
    INTEGER, INTENT(IN) :: id_i1, id_j1, id_isize, id_jsize
    REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1), INTENT(OUT) :: rda_f
    REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize), INTENT(IN) :: rda_x
    INTEGER ibi, ibj, il_idx

    il_idx = 0
    DO ibi = id_i1, id_i1 + id_isize - 1
        DO ibj = id_j1, id_j1 + id_jsize - 1
          il_idx = il_idx + 1
          rda_f(ibi,ibj) = rda_x(il_idx)
        ENDDO
    ENDDO

  END SUBROUTINE vec2bloc2

   SUBROUTINE bloc2vec_3d(rda_f,id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize, rda_x)
      !Column major storage assumption for images processing
      !
      !** Vectorize le bloc de la matrice rda_f(id_i1 : id_i1 + id_isize -1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1)
      !   dans rda_x(1 : id_isize * id_jsize * id_ksize)
      !   PAR LIGNES
      INTEGER, INTENT(IN) :: id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize
      REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1),&
                                          INTENT(IN) :: rda_f
      REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize * id_ksize), INTENT(OUT) :: rda_x
      INTEGER ibi, ibj, ibk, il_idx

      il_idx = 0
            DO ibk = id_k1, id_k1 + id_ksize - 1
         DO ibj = id_j1, id_j1 + id_jsize - 1
      DO ibi = id_i1, id_i1 + id_isize - 1
               il_idx = il_idx + 1
               rda_x(il_idx) = rda_f(ibi,ibj,ibk)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE bloc2vec_3d

   SUBROUTINE vec2bloc_3d(rda_x,id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize, rda_f)
      !Column major storage assumption for images processing
      !
      !** Distribue le vecteur rda_x dans le bloc de la matrice
      !   rda_f(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1)
      !   /!\ PAR LIGNES
      INTEGER, INTENT(IN) :: id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize
      REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1), &
                          INTENT(OUT) :: rda_f
      REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize * id_ksize), INTENT(IN) :: rda_x
      INTEGER ibi, ibj, ibk, il_idx

      il_idx = 0
            DO ibk = id_k1, id_k1 + id_ksize - 1
         DO ibj = id_j1, id_j1 + id_jsize - 1
      DO ibi = id_i1, id_i1 + id_isize - 1
               il_idx = il_idx + 1
               rda_f(ibi,ibj,ibk) = rda_x(il_idx)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE vec2bloc_3d

   SUBROUTINE bloc2vec_img_3d(rda_f,id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize, rda_x)
      ! 2D equivalent of bloc2vec2
      !Row major storage assumption for images processing
      !
      !** Vectorize le bloc de la matrice rda_f(id_i1 : id_i1 + id_isize -1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1)
      !   dans rda_x(1 : id_isize * id_jsize * id_ksize)
      !   PAR LIGNES
      INTEGER, INTENT(IN) :: id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize
      REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1),&
                                          INTENT(IN) :: rda_f
      REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize * id_ksize), INTENT(OUT) :: rda_x
      INTEGER ibi, ibj, ibk, il_idx

      il_idx = 0
      DO ibi = id_i1, id_i1 + id_isize - 1
         DO ibj = id_j1, id_j1 + id_jsize - 1
            DO ibk = id_k1, id_k1 + id_ksize - 1
               il_idx = il_idx + 1
               rda_x(il_idx) = rda_f(ibi,ibj,ibk)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE bloc2vec_img_3d

   SUBROUTINE vec2bloc_img_3d(rda_x,id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize, rda_f)
      ! 2D equivalent of vec2bloc2
      !Row major storage assumption for images processing
      !
      !** Distribue le vecteur rda_x dans le bloc de la matrice
      !   rda_f(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1)
      !   /!\ PAR LIGNES
      INTEGER, INTENT(IN) :: id_i1, id_j1, id_k1, id_isize, id_jsize, id_ksize
      REAL(KIND=dp), DIMENSION(id_i1 : id_i1 + id_isize - 1, id_j1 : id_j1 + id_jsize - 1, id_k1 : id_k1 + id_ksize - 1), &
                          INTENT(OUT) :: rda_f
      REAL(KIND=dp), DIMENSION(1 : id_isize * id_jsize * id_ksize), INTENT(IN) :: rda_x
      INTEGER ibi, ibj, ibk, il_idx

      il_idx = 0
      DO ibi = id_i1, id_i1 + id_isize - 1
         DO ibj = id_j1, id_j1 + id_jsize - 1
            DO ibk = id_k1, id_k1 + id_ksize - 1
               il_idx = il_idx + 1
               rda_f(ibi,ibj,ibk) = rda_x(il_idx)
            ENDDO
         ENDDO
      ENDDO
   END SUBROUTINE vec2bloc_img_3d

  !> \brief (pseudo) Adjoint of MAXVAL function
  !! @param [in] rda_A original vector
  !! @param [in,out] rda_Aad adjoint vector
  !! @param [in, out] rd_maxad
  !<
  SUBROUTINE MAXVALADJ_1d(rda_A, rda_Aad, rd_maxad)
    !     Purpose: adjoint de MAXVAL, avec une matrice en parametre
    !     --------
    !     Interface : rda_A : matrice originale
    !                 rda_A : matrice adjointe
    !     ----------
    !     History:
    !     --------
    !      Version    Programmer      Date         Description
    !      -------    ----------      ----         -----------
    !       1.0        Innocent       04/09         Creation
    !*-----------------------------------------------------------------

    REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: rda_A
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)  :: rda_Aad
    REAL(KIND=dp)              , INTENT(INOUT)  :: rd_maxad
    INTEGER :: il_imax, ib_i

    il_imax = LBOUND(rda_A, 1)

    DO ib_i = LBOUND(rda_A, 1), UBOUND(rda_A, 1)
      IF( rda_A(il_imax) < rda_A(ib_i) ) THEN
        il_imax = ib_i
      END IF
    END DO
    rda_Aad(il_imax) = rda_Aad(il_imax) + rd_maxad
    rd_maxad = 0.0_cp
  END SUBROUTINE MAXVALADJ_1d

  SUBROUTINE MAXVALADJ_2d(rda_A, rda_Aad, rd_maxad)
    !     Purpose: adjoint de MAXVAL, avec une matrice en param�tre
    !     --------
    !     Interface : rda_A : matrice originale
    !                 rda_A : matrice adjointe
    !     ----------
    !     History:
    !     --------
    !      Version    Programmer      Date         Description
    !      -------    ----------      ----         -----------
    !       1.0        Innocent       04/09         Creation
    !*-----------------------------------------------------------------

    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: rda_A
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)  :: rda_Aad
    REAL(KIND=dp)                 , INTENT(INOUT)  :: rd_maxad
    INTEGER                               :: il_imax, il_jmax, ib_i, ib_j

    il_imax = LBOUND(rda_A, 1)
    il_jmax = LBOUND(rda_A, 2)

    DO ib_j = LBOUND(rda_A, 2), UBOUND(rda_A, 2)
      DO ib_i = LBOUND(rda_A, 1), UBOUND(rda_A, 1)
        IF( rda_A(il_imax, il_jmax) < rda_A(ib_i, ib_j) ) THEN
          il_imax = ib_i
          il_jmax = ib_j
        END IF
      END DO
    END DO
    rda_Aad(il_imax, il_jmax) = rda_Aad(il_imax, il_jmax) + rd_maxad
    rd_maxad = 0.0_cp
  END SUBROUTINE MAXVALADJ_2d
! Compute the value of Whittaker function W_{(\alpha-1)/2, 1/2}
   ! For the cases alpha = 1 or alpha = 3 only
   ! /!\ Need to be completed by a general formula for all alpha
   ! see Gradshteyn, I.S., Ryshik, I.M. 1980, Tables of Integrals, Academic Press
!    FUNCTION whittaker_olivier( r )
!       IMPLICIT NONE
!       REAL*8, INTENT( IN ) :: r
!       REAL*8               :: whittaker
!       REAL*8               :: w1, w3
!       ! (alpha = 1) non-isolated vortex (gaussian velocity profile)
!       ! Warning : w1 does not reach its min and max values at r = pprvm !
!       IF ( ABS(r) <= 1E-16 ) THEN
!          w1 = 0.0
!       ELSE
!          w1 = rg_vm * rg_rvm * (1.0 - exp( -0.5 * ( r / rg_rvm ) ** 2 ) ) / r
!       END IF
!
!       ! (alpha = 3) isolated vortex (gaussian h profile)
!       IF( ABS(r) <= 1E-16 ) THEN
!          w3 = 0.0
!       ELSE
!          w3 = rg_vm * exp( 0.5 ) * ( r / rg_rvm ) * exp( -0.5 * ( r / rg_rvm )** 2 )
!       END IF
!
!       ! Choose either w1 or w3 depending if the vortex is  isolated (w3) or non-isolated (w1)
!       whittaker = w3
!
!    END FUNCTION whittaker_olivier

  !> \brief Compute the value of Whittaker function W_{(alpha-1)/2, 1/2}.
  !! @param [in] r distance betwen the center and the point where the function is evaluated
  !! @param [in] rd_vm maximum azimuthal velocity
  !! @param [in] rd_rvm radius of maximum azimuthal velocity
  !! \details The Whittaker function as defined here is for the cases alpha = 1 or alpha = 3 only
  !! /!\ Need to be completed by a general formula for all alpha
  !! see Gradshteyn, I.S., Ryshik, I.M. 1980, Tables of Integrals, Academic Press
  !<
  FUNCTION whittaker( r, rd_vm, rd_rvm )
    !!CALL as xx = whittaker (r, tg_swp%r_vm, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: r, rd_vm, rd_rvm
    REAL(KIND=dp)               :: whittaker
    REAL(KIND=dp)               :: w1, w3
    ! (alpha = 1) non-isolated vortex (gaussian velocity profile)
    ! Warning : w1 does not reach its min and max values at r = pprvm !
    IF ( ABS(r) <= 1d-16 ) THEN
        w1 = 0.0_dp
    ELSE
        w1 = rd_vm * rd_rvm * (1.0 - exp( -0.5_cp * ( r / rd_rvm ) ** 2 ) ) / r
    END IF

    ! (alpha = 3) isolated vortex (gaussian h profile)
    IF( ABS(r) <= 1d-16 ) THEN
        w3 = 0.0_dp
    ELSE
        w3 = rd_vm * exp( 0.5_cp) * ( r / rd_rvm ) * exp( -0.5_cp * ( r / rd_rvm )** 2 )
    END IF

    ! Choose either w1 or w3 depending if the vortex is  isolated (w3) or non-isolated (w1)
    whittaker = w3

  END FUNCTION whittaker

  !> \brief Computes the value of the u component (x axis) of the initial velocity of the vortex at (x,y)
  !! @param [in] x x-coordinate of the point where the velocity is evaluated
  !! @param [in] y y-coordinate of the point where the velocity is evaluated
  !! @param [in] rd_xvc x-coordinate of the center of the vortex
  !! @param [in] rd_yvc y-coordinate of the center of the vortex
  !! @param [in] rd_vm maximum azimuthal velocity
  !! @param [in] rd_rvm radius of maximum azimuthal velocity
  !! Uses whittaker function
  !<
  FUNCTION u_init( x, y, rd_xvc, rd_yvc, rd_vm, rd_rvm)
    !!Call as xx = u_init(x, y, tg_swp%r_xvc, tg_swp%r_yvc, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: x, y, rd_xvc, rd_yvc, rd_rvm, rd_vm
    REAL(KIND=dp)               :: r
    REAL(KIND=dp)               :: u_init

    r = SQRT( ( x - rd_xvc ) ** 2 + ( y - rd_yvc ) ** 2 )

    IF( ABS(r) < 1d-16 ) THEN
        u_init = 0.0_dp
    ELSE
        u_init =  - ( y - rd_yvc ) * whittaker( r, rd_vm, rd_rvm ) / r
    END IF

    !* Simulation de l'initialisation par agitation dans un cylindre
    !* que l'on retire par la suite
    !IF( ( ABS(r) < 1E-16 ) .OR. ( ABS(r) > 5 * rd_rvm ) ) THEN
    !   u_init = 0.0
    !ELSE
    !   u_init =  - ( y - rd_yvc ) * whittaker( r ) / r
    !END IF

  END FUNCTION u_init

  !> \brief Computes the value of the v component (y axis) of the initial velocity of the vortex at (x,y)
  !! @param [in] x x-coordinate of the point where the velocity is evaluated
  !! @param [in] y y-coordinate of the point where the velocity is evaluated
  !! @param [in] rd_xvc x-coordinate of the center of the vortex
  !! @param [in] rd_yvc y-coordinate of the center of the vortex
  !! @param [in] rd_vm maximum azimuthal velocity
  !! @param [in] rd_rvm radius of maximum azimuthal velocity
  !! Uses whittaker function
  !<
  FUNCTION v_init( x, y, rd_xvc, rd_yvc, rd_vm, rd_rvm )
    !!Call as xx = v_init(x, y, tg_swp%r_xvc, tg_swp%r_yvc, tg_swp%r_rvm)
    IMPLICIT NONE
    REAL(KIND=dp), INTENT( IN ) :: x, y, rd_xvc, rd_yvc, rd_rvm, rd_vm
    REAL(KIND=dp)               :: r
    REAL(KIND=dp)               :: v_init

    r = SQRT( ( x - rd_xvc ) ** 2 + ( y - rd_yvc ) ** 2 )

    IF ( ABS(r) <= 1d-16 ) THEN
        v_init = 0.0_dp
    ELSE
        v_init =  ( x - rd_xvc ) * whittaker( r, rd_vm, rd_rvm ) / r
    END IF

    !* Simulation de l'initialisation par agitation dans un cylindre
    !* que l'on retire par la suite
    ! IF( ( ABS(r) <= 1E-16 ) .OR. ( ABS(r) > 5 * rd_rvm ) ) THEN
    !    v_init = 0.0
    ! ELSE
    !    v_init =  ( x - rd_xvc ) * whittaker( r ) / r
    ! END IF

  END FUNCTION v_init

  !!!! Integration


    !> \brief Compute the approximation of the integral using trapezoidal rule on a regularly spaced grid
    !! @param [in] rda_data 1D array of samples of the function that the integral is needed, sample are given at regularly spaced grid points.
    !! @param [in] ida_pt 1D array of the same size as rda_data. Gives for each sample number of dimensions in which the data point is external. ida_pt takes integer values between zero and ndim
    !! @param [in] delta_x gives the space step in each dimension
    !! \details for example, in dimension 2 with physical domain defined by [a,b]x[c,d]
    !! - points (a,c), (a,d), (b,c) and (b,d) are external in two dimension, the associated face is 2
    !! - points (a,y), (b,y), (x,c) and (x,d) with a<x<b and c<y<d are external in one dimension, the associated face is 1
    !! - all other points (x,y) with  a<x<b and c<y<d are external in zero dimension the associated face is 0
    !<
    function trapeze_regular(rda_data, ida_pt, delta_x)result(rl_trap)
        real(cp), dimension(:), intent(in) :: rda_data, delta_x
        integer, dimension(:), intent(in) :: ida_pt
        !local variables
        real(cp) :: rl_trap

        rl_trap = product(delta_x) * sum( rda_data/2**ida_pt )
    end function trapeze_regular

    subroutine print_matrix(rd_mat)
        real(kind=dp), dimension(:, :), intent(in) :: rd_mat
        integer :: ibi!, ibj
        character(len=50) ala_format, ala_nb

        write (ala_nb,'(I10)') SIZE(rd_mat, 2)
        ala_format = '('//trim(ala_nb)//'(X,E10.3))'
        !print*,''
        do ibi = 1, size(rd_mat, 1)
        !do ibj = 1, size(rd_mat, 2)
            write(*,ala_format) rd_mat(ibi, :)
        !end do
        end do
    end subroutine print_matrix

   !> \brief Fills a matrix with vector data following a given ordering
   !! @param [in] rda_vec input vector
   !! @param [out] rda_mat output matrix
   !! @param [in] major (optional) the order of elements in the vector
   !!   default is raw major ordering
   !!
   !! \details This procedure uses the reshape function
   !!
   !<
   subroutine vec2mat_sp( rda_vec, rda_mat, major )
      real(kind=sp), dimension(:), intent(in):: rda_vec
      real(kind=sp), dimension(:,:), intent(in out):: rda_mat
      integer, optional, intent(in) :: major
      !local variables
      integer, dimension(2) :: order = (/1,2/)
      integer:: il_major
      if ( present(major) )then
        il_major = major
      else
        il_major = ROW_WISE
      end if

      select case(il_major)
        case (ROW_WISE)
          order = (/2,1/)
        case (COL_WISE)
          order = (/1,2/)
      end select

      rda_mat = reshape( rda_vec, shape(rda_mat), ORDER=order )
   end subroutine vec2mat_sp

   !> \brief Fills a matrix with vector data following a given ordering
   !! @param [in] rda_vec input vector
   !! @param [out] rda_mat output matrix
   !! @param [in] major (optional) the order of elements in the vector
   !!
   !! \details This procedure uses the reshape function
   !!
   !<
   subroutine vec2mat_dp( rda_vec, rda_mat, major )
      real(kind=dp), dimension(:), intent(in):: rda_vec
      real(kind=dp), dimension(:,:), intent(in out):: rda_mat
      integer, optional, intent(in) :: major
      !local variables
      integer, dimension(2) :: order = (/1,2/)
      integer:: il_major
      if ( present(major) )then
        il_major = major
      else
        il_major = ROW_WISE
      end if

      select case(il_major)
        case (ROW_WISE)
          order = (/2,1/)
        case (COL_WISE)
          order = (/1,2/)
      end select

      rda_mat = reshape( rda_vec, shape(rda_mat), ORDER=order )
   end subroutine vec2mat_dp

   !> \brief Fills a 3D array with vector data following a given ordering
   !! @param [in] rda_vec input vector
   !! @param [out] rda_mat output matrix
   !! @param [in] major (optional) the order of elements in the vector
   !!
   !! \details This procedure uses the reshape function
   !!
   !<
   subroutine vec2mat3d_dp( rda_vec, rda_mat, major )
      real(kind=dp), dimension(:), intent(in):: rda_vec
      real(kind=dp), dimension(:,:,:), intent(in out):: rda_mat
      integer, optional, intent(in) :: major
      !local variables
      integer, dimension(3) :: order = (/1,2,3/)
      integer:: il_major

      if ( present(major) )then
        il_major = major
      else
        il_major = ROW_WISE
      end if

      select case(il_major)
        case (ROW_WISE)
          order = (/3,2,1/)
        case (COL_WISE)
          order = (/1,2,3/)
      end select

      rda_mat = reshape( rda_vec, shape(rda_mat), ORDER=order )
   end subroutine vec2mat3d_dp

   !> \brief Vectorises a matrix
   !! @param [in] rda_mat output matrix
   !! @param [out] rda_vec input vector
   !! @param [in] major (optional) the order of elements in the vector
   !<
   subroutine mat2vec_sp(rda_mat, rda_vec, major)
      real(kind=sp), dimension(:,:), intent(in):: rda_mat
      real(kind=sp), dimension(:), intent(in out):: rda_vec
      integer, optional, intent(in) :: major
      !local variables
      integer:: il_major
      if ( present(major) )then
        il_major = major
      else
        il_major = ROW_WISE
      end if

      select case(il_major)
        case (ROW_WISE)
          rda_vec = reshape( transpose(rda_mat), (/size(rda_vec)/) )
        case (COL_WISE)
          rda_vec = reshape( rda_mat, (/size(rda_vec)/) )
      end select
   end subroutine mat2vec_sp

   !> \brief Vectorises a matrix
   !! @param [in] rda_mat output matrix
   !! @param [out] rda_vec input vector
   !! @param [in] major (optional) the order of elements in the vector
   !<
   subroutine mat2vec_dp(rda_mat, rda_vec, major)
      real(kind=dp), dimension(:,:), intent(in):: rda_mat
      real(kind=dp), dimension(:), intent(in out):: rda_vec
      integer, optional, intent(in) :: major
      !local variables
      integer:: il_major
      if ( present(major) )then
        il_major = major
      else
        il_major = ROW_WISE
      end if

      select case(il_major)
        case (ROW_WISE)
          rda_vec = reshape( transpose(rda_mat), (/size(rda_vec)/) )
        case (COL_WISE)
          rda_vec = reshape( rda_mat, (/size(rda_vec)/) )
      end select
   end subroutine mat2vec_dp



    !> @brief computes the smallest power of ten greater than @a n
    !! @param [in] n number whose closest large power of the is needed
    !! @details for n between \f$ 10^k \f$ [excluded] and \f$ 10^(k+1) \f$ [included]
    !!  the program returns \f$ 10^(k+1) \f$
    !! @note this programs deals only with positive numbers.
    !<
    function nextpow10(n)result(next)
        integer, intent(in) :: n
        real(dp) :: rl_l10
        integer  :: il_l10, next

        rl_l10 = log10( real(n,dp) )
        il_l10 = int(rl_l10)
        if (il_l10 < rl_l10)then
            next = il_l10 + 1
        else
            next = il_l10
        end if
        next = 10**next
    end function nextpow10

  SUBROUTINE writeMatrix(fileName, rda_A)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN) :: rda_A
    CHARACTER(len=*)         , INTENT(IN) :: fileName
    INTEGER ib_i, ib_j, il_ios
    INTEGER :: ip_fid = 41

    IF(fileName==SCREEN) THEN
      ip_fid = 6!standard output
      PRINT*, 'rows count = ', size(rda_A, 1)
      PRINT*, 'columns count = ', size(rda_A, 2)
    ELSE
      ip_fid = 41
      OPEN( UNIT = ip_fid, FILE = fileName, STATUS = 'REPLACE', FORM = 'FORMATTED', IOSTAT = il_ios)
      IF (il_ios /= 0) CALL stop_program( 'In writeMatrix : Error creating file'//fileName )
    END IF
    PRINT*, "In writeMatrix", LBOUND(rda_A,1), UBOUND(rda_A,1), LBOUND(rda_A,2), UBOUND(rda_A,2)
    !stop
    WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 1)
    WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 2)
    DO ib_j = LBOUND(rda_A,2), UBOUND(rda_A,2)
      DO ib_i = LBOUND(rda_A,1), UBOUND(rda_A,1)
        WRITE(unit=ip_fid, fmt='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid end of line
      END DO
      WRITE(unit=ip_fid, fmt=*)''! just for adding end of line
    END DO
    IF(fileName/=SCREEN) THEN
      CLOSE(ip_fid )
    END IF
  END SUBROUTINE writeMatrix

  SUBROUTINE readMatrix(fileName, rda_A)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT) :: rda_A
    CHARACTER(len=*)         , INTENT(IN)  :: fileName
    INTEGER ib_j, il_nbRow, il_nbCol, il_ios!, ib_i
    INTEGER , PARAMETER :: ip_fid =615

    OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
    IF (il_ios /= 0) CALL stop_program( 'In readMatrix : Error opening file '//fileName )
!
    READ(UNIT=ip_fid, FMT=IFMT) il_nbRow
    READ(UNIT=ip_fid, FMT=IFMT) il_nbCol
    IF((SIZE(rda_A, 1)/=il_nbRow).OR.(SIZE(rda_A, 2)/=il_nbCol))THEN
      CALL debug('In readMatrix, bad shape of matrix parameter, based on the file data shape:', tag=dALLWAYS)
      CALL debug((/il_nbRow, il_nbCol/), 'expected shape :', tag=dALLWAYS)
      CALL debug(SHAPE(rda_A), 'provided shape :', tag=dALLWAYS)
      CALL stop_program()
    END IF
    !CALL debug( (/il_nbRow,il_nbCol/), "In readMatrix("//TRIM(fileName)//"), (il_nbRow, il_nbCol) = ")
    DO ib_j = 1, il_nbCol
      READ(UNIT=ip_fid, FMT=*)rda_A(:,ib_j)
!       DO ib_i = 1, il_nbRow
!         PRINT*, "In readMatrix, ib_j = ", ib_j, "; ib_i = ", ib_i
!         READ(UNIT=ip_fid, FMT='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid going to the next line
!       END DO
!       READ(unit=ip_fid, fmt=*) !just to skip the end of line
    END DO
    CLOSE(ip_fid )
  END SUBROUTINE readMatrix

  !!read the head of the file: nbRow and nbCol
  SUBROUTINE readInfo(fileName, id_nbRow, id_nbCol)
    INTEGER , PARAMETER :: ip_fid =616
    CHARACTER(len=*), INTENT(IN)  :: fileName
    INTEGER, INTENT(OUT) :: id_nbRow, id_nbCol
    INTEGER :: il_ios

    OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
    IF (il_ios /= 0) CALL stop_program( 'In readInfo : Error opening file '//fileName )
    READ(UNIT=ip_fid, FMT=IFMT) id_nbRow
    READ(UNIT=ip_fid, FMT=IFMT) id_nbCol
    CALL debug( (/id_nbRow,id_nbCol/), "In readInfo("//TRIM(fileName)//"), (il_nbRow, il_nbCol) = ")
    CLOSE(ip_fid )
  END SUBROUTINE readInfo

!   !> @brief initialize the random seed with a varying seed in order to ensure a
!   !>   different random number sequence for each invocation of the program.
!   !!
!   !! \details
!   !!   not usefull where reproducibility is a concern
!   !!   courtesy of GNU project: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!   !<
!   subroutine init_random_seed_optimal()
!     implicit none
!     integer, allocatable :: seed(:)
!     integer :: i, n, istat, dt(8), pid
!     integer(int64) :: t
!
!     integer, parameter :: un=703
!
!     call random_seed(size = n)
!     allocate(seed(n))
!     ! First try if the OS provides a random number generator
!     open(unit=un, file="/dev/urandom", access="stream", &
!         form="unformatted", action="read", status="old", iostat=istat)
!     if (istat == 0) then
!       read(un) seed
!       close(un)
!     else
!       ! Fallback to XOR:ing the current time and pid. The PID is
!       ! useful in case one launches multiple instances of the same
!       ! program in parallel.
!       !call system_clock(t)
!       if (t == 0) then
!         call date_and_time(values=dt)
!         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
!           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
!           + dt(3) * 24_int64 * 60 * 60 * 1000 &
!           + dt(5) * 60 * 60 * 1000 &
!           + dt(6) * 60 * 1000 + dt(7) * 1000 &
!           + dt(8)
!       end if
!       pid = getpid()
!       t = ieor( t, int(pid, kind(t)) )
!       do i = 1, n
!         seed(i) = lcg(t)
!       end do
!     end if
!     call random_seed(put=seed)
!   end subroutine init_random_seed_optimal
!
!
!     !> @brief initialize the random seed with a reproducible seed
!     !! @param [in] seed_base user provided number to generate the reproducible seed.
!     !! @details there is garantee that whenever the subroutine is call with the
!     !!   same seed, the sequence will be identical
!     !<
!     subroutine init_random_seed_simple( seed_base )
!         implicit none
!         integer, intent(in) :: seed_base
!         integer, allocatable :: seed(:)
!         integer :: i, n
!         integer(int64) :: t
!
!         integer, parameter :: un=703
!
!         call random_seed( size = n )
!         allocate( seed(n) )
!         !Add a large prime
!         t = int( abs(seed_base), kind(t) ) + 2038102361 ! 32416187567
!         do i = 1, n
!             seed(i) = lcg(t)
!         end do
!
!         call random_seed( put=seed )
!     end subroutine init_random_seed_simple
!
!   !> @brief This simple PRNG might not be good enough for real work, but is
!   !! sufficient for seeding a better PRNG.
!   !! \details courtesy of GNU project: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!   !<
!   function lcg(s)
!     integer :: lcg
!     integer(int64), intent(in out) :: s
!     if (s == 0) then
!       s = 104729
!     else
!       s = mod(s, 4294967296_int64)
!     end if
!     s = mod(s * 279470273_int64, 4294967291_int64)
!     lcg = int(mod(s, int(huge(0), int64)), kind(0))
!   end function lcg
!
!   !> \brief Ensures that a matrix contained in a file has a given shape or stop the program
!   !! \brief
!   !<
!   SUBROUTINE ensure_matrixShape_orDie(fileName, nrow, ncol)
!     CHARACTER*(ip_fnl), INTENT(IN) :: fileName
!     INTEGER, OPTIONAL, INTENT(IN)  :: nrows, ncols
!     !local variables
!     INTEGER :: il_nrow, il_ncol, il_f_nrow, il_f_ncol
!
!     IF
!     CALL readInfo(fileName, il_f_nrow, il_f_ncol)
!
!   END SUBROUTINE ensure_matrixSize_orDie

  !> \brief Says if a time step is (to be) saved accordind to subsampling parameters
  !! @param [in] id_n time step num
  !! @param [in] id_trjSub subsampling parameter, this parameter gives the number of time steps between two successive records in file
  !! @param [in] id_numInFile record number
  !! @param [in] tshift (optional) gives the time step from wich the saving starts.
  !!If the function returns, \a id_numInFile contains the record number in the nc file
  !!If \a shift is present, no time step before \a shift is saved
  !<
  FUNCTION savedstep(id_n, id_trjSub, id_numInFile, tshift)RESULT(ll_res)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id_n
    INTEGER, INTENT(IN) :: id_trjSub
    INTEGER, INTENT(OUT):: id_numInFile
    INTEGER, OPTIONAL, INTENT(IN) :: tshift
    !local variables
    LOGICAL             :: ll_res
    INTEGER :: il_num

    IF( PRESENT(tshift) )THEN
      il_num = id_n - tshift
    ELSE
      il_num = id_n
    END IF

    IF ( (il_num>=0).AND.(MOD(il_num, id_trjSub) == 0) ) THEN
      id_numInFile = (il_num)/id_trjSub+1
      ll_res = .TRUE.
    ELSE
      id_numInFile = -1 !
      ll_res = .FALSE.
    END IF
  END FUNCTION savedstep


    !Date time operations



   !> @brief returns the leap status of a given year, .true. for leap and .false. if not
   !! @param [in] year to check
   !<
   function isLeap(year) result(isL)
      integer, intent(in) :: year
      !local variables
      logical :: isL

      isL = (mod(year,4)==0.and..not.mod(year,100)==0).or.(mod(year,400)==0)
   end function isLeap

   !> @brief returns the number of days in a given month of a given year
   !! @param [in] month month to compute the number of days in
   !! @param [in] year from with the month is taken
   !<
   function daysInMonth(month,year)result(il_dim)
      integer,intent(in) :: month
      integer,intent(in) :: year
      !local variables
      integer,parameter,dimension(12) :: days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: il_dim

      if( (month < 1).or.(month > 12) )then
         print*, "In daysInMonth: bad month ", month
         stop
      else if( (month == 2).and.isLeap(year) ) then
         il_dim = 29
      else
         il_dim = days(month)
      end if
   endfunction daysInMonth

   !> @brief Compute the next date given the current date
   !! @param [in] currentDate date (format YYYYMMDD) from which the next
   !!   one is computed
   !! @param [in] inc optional number of days to add to the current date
   !! @details Only positive dates are supported. The input integer date must be
   !!    greater than 101 which is 0000-01-01.
   !<
   function nextDate(currentDate, inc) result(next)
      integer, intent(in) :: currentDate
      integer, optional, intent(in) :: inc
      !local variables
      integer :: year, month, day, tmp, daysInCurrentMonth, il_inc
      logical :: nextIter
      integer :: next

      if ( present(inc) )then
         il_inc = inc
      else
         il_inc = 1
      end if

      year = currentDate/10000
      tmp = mod(currentDate, 10000)
      month = tmp/100
      day = mod(tmp, 100)
      if(currentDate<101)then
         print*, "In nextDate: bad current date ", currentDate
         stop
      else
         if( (month < 1).or.(month > 12) )then
            print*, "In nextDate: bad month ", month
            stop
         else
            if( (day < 1).or.(day > daysInMonth(month, year)) )then
               print*, "In nextDate: bad day (for the cuurent month and year) ", day
               stop
            end if
         end if
      end if

      day = day+il_inc
      nextIter = .true.
      do while(nextIter)
         daysInCurrentMonth = daysInMonth(month, year)
         if( day > daysInCurrentMonth )then
            day = day-daysInCurrentMonth
            month = month+1
            if( month > 12 )then
               year = year+month/12
               month = mod(month,12)
            end if
         else if( day < 1 )then
            month = month-1
            if( month < 1 )then
               year = year+month/12-1
               month = 12+mod( month,12 )
            end if
            day = day+daysInMonth(month,year)
         else
            nextIter = .false.
         end if
      end do
      next = year*10000 + month*100 + day
   end function nextDate

    !> @brief returns .true. if a date (YYYYMMDD) is good or .false. if not
    !! @details Only positive dates are supported. The input integer date must be
    !!    greater than 101 which is 0000-01-01.
    function goodDate(currentDate)result(good)
      integer, intent(in) :: currentDate
      !local variables
      integer :: year, month, day, tmp
      logical :: good

      year = currentDate/10000
      tmp = mod(currentDate, 10000)
      month = tmp/100
      day = mod(tmp, 100)

      good = ( (year>=0).and.(month >=1).and.(month <= 12)&
           .and.(day >= 1).and.(day <= daysInMonth(month, year)) )
    end function goodDate

    !> @brief returns .true. if a time (HHmmSS) is good or .false. if not
    !<
    function goodTime(currentTime)result(good)
        integer, intent(in) :: currentTime
        !local variables
        integer :: h, m, s, tmp
        logical :: good

        h = currentTime/10000
        tmp = mod(currentTime, 10000)
        m = tmp/100
        s = mod(tmp, 100)

      good = ( (h>=0).and.(h < 24).and.(m >= 0).and.(m < 60 )&
           .and.(s >= 0).and.(s < 60 ) )
    end function goodTime

    !> @brief update date and time
    !! @param [in, out] currentDate, currentTime date and time to update, format YYYYMMDD, HHmmSS
    !! @param [in] hour, minute, second increment to add to the date/time
    !! @details if none of hour, minute, second is provide
    !!   there will be no update. For each increment that is provided,
    !!   the appropriate update is performed, staring with second up to year
    !! Only positive dates and times are supported. The input integer date must be
    !!    greater than 101 which is 0000-01-01.
    !<
    subroutine update_dt(currentDate, currentTime, hour, minute, second)
        integer, intent(in out) :: currentDate, currentTime
        integer, optional, intent(in) :: hour, minute, second
        !local variables
        integer :: tmp, inc_hour, inc_min, inc_sec, h, m, s, inc_day
        !logical :: nextIter
        !integer :: next

        inc_day = 0

        if ( present(hour) )then
            inc_hour = hour
        else
            inc_hour = 0
        end if

        if ( present(minute) )then
            inc_min = minute
        else
            inc_min = 0
        end if

        if ( present(second) )then
            inc_sec = second
        else
            inc_sec = 0
        end if

        h = currentTime/10000
        tmp = mod(currentTime, 10000)
        m = tmp/100
        s = mod(tmp, 100)

        if((currentTime<0).or.(currentTime>235959))then
            print*, "In update_dt: bad current time ", currentTime
            print*, "time must be tetween 0 and 235959 "
            stop
        else if( (m < 0).or.(m >= 60).or.(s<0).or.(s>=60) )then
            print*, "In update_dt: bad minute or second ", m, s
            print*, " minute and second must be between 0 and 59 "
            stop
        end if

        !updating the seconds
        s = s + inc_sec
        if(s>=60)then
            inc_min = inc_min + s/60
            s = mod(s, 60)
        else if(s<0)then
            !in modulo, the sign of the result is the sign of the denominator
            tmp = s
            s = modulo(s, 60)
            inc_min = inc_min + (tmp-s)/60
        end if

        !updating the minutes
        m = m + inc_min
        if(m>=60)then
            inc_hour = inc_hour + m/60
            m = mod(m, 60)
        else if(m<0)then
            !in modulo, the sign of the result is the sign of the denominator
            tmp = m
            m = modulo(m, 60)
            inc_hour = inc_hour + (tmp-m)/60
        end if

        !updating the hours
        h = h + inc_hour
        if(h>=24)then
            inc_day = h/24
            h = mod(h, 24)
        else if(h<0)then
            !in modulo, the sign of the result is the sign of the denominator
            tmp = h
            h = modulo(h, 24)
            inc_day = (tmp-h)/24
        end if

        currentTime = h*10000 + m*100 + s
        !updating the date
        if(inc_day/=0) currentDate = nextDate(currentDate, inc_day)
    end subroutine update_dt

end module general_tools
