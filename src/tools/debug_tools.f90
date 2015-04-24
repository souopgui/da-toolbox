!> \file debug_tools.f90
!! Set of functions/routines designed for debugging
!<

!> \brief Module defining interactive routines for debugging
!!
!<
module debug_tools
    use general_constant
    use, intrinsic :: iso_c_binding, only : C_INT
implicit none

    logical, parameter :: LP_DEBUG = .TRUE.
    !> Screen file name for output
    character (len=*), parameter :: SCREEN = "SCREEN"

    !understanding fortran io formating, see http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
    !> \brief format to print integer
    character(len=*), parameter :: IFORMAT = "(I5)"

    character(len=*), parameter :: IFMT = IFORMAT

    !> \brief format to print real
    character(len=*), parameter :: RFORMAT = "(ES14.5E3)"
    character(len=*), parameter :: RFMT = "(ES14.5E3)"

    !> \brief format to print ascii (characters string)
    character(len=*), parameter :: AFORMAT = "(A)"

    !> \brief format to print integer with comment, debug
    character(len=*), parameter :: DEBUG_IFORMAT = "(A,': ',"//IFORMAT//")"

    !> \brief format to print logical value with comment, debug
    character(len=*), parameter :: DEBUG_LFORMAT = "(A,': ',L2)"

    !> \brief format to print character value with comment, debug
    character(len=*), parameter :: DEBUG_AFORMAT = "(A,': ',A)"

    !> \brief format to print real with comment, debug
    character(len=*), parameter :: DEBUG_RFORMAT = "(A,': ',X,"//RFORMAT//")"

    !> @brief minmax format
    character(len=*), parameter:: DEBUG_MMFORMAT = &
        "(A,':  * min = ',"//RFMT//",'  * max = ',"//RFMT//")"

    !> Interface defining routines used to convert numbers to characters
    !<
    INTERFACE NUM2STR
        MODULE PROCEDURE itoa!> Integer to ASCII
        MODULE PROCEDURE rtoa!> Real single precision to ASCII
        MODULE PROCEDURE dtoa!> Real double precision to ASCII
    END INTERFACE

    INTERFACE stop_program
        MODULE PROCEDURE stop_program_with_int
        MODULE PROCEDURE stop_program_with_string
        MODULE PROCEDURE stop_program_with_sp_real
        MODULE PROCEDURE stop_program_with_dp_real
    END INTERFACE stop_program

    INTERFACE pretty_format
        MODULE PROCEDURE pretty_format_int_scalar
        MODULE PROCEDURE pretty_format_int_1D_array
        MODULE PROCEDURE pretty_format_int_2D_array
    END INTERFACE pretty_format

    INTERFACE pretty_dformat
        MODULE PROCEDURE pretty_dformat_int_scalar
    END INTERFACE pretty_dformat

    !> \brief Debug routines, used to print values of variables with a predefined comment
    !! The goal is to make it possible to activate or deactivate the printing with a global flag
    !<
    INTERFACE DEBUG
        MODULE PROCEDURE debug_string
        MODULE PROCEDURE debug_string_val
        MODULE PROCEDURE debug_integer
        MODULE PROCEDURE debug_sp_real
        MODULE PROCEDURE debug_dp_real
        MODULE PROCEDURE debug_logical
        !MODULE PROCEDURE debug_c_int_array1D
        MODULE PROCEDURE debug_integer_array1D
        MODULE PROCEDURE debug_integer_array2D
        MODULE PROCEDURE debug_sp_real_array1D
        MODULE PROCEDURE debug_sp_real_array2D
        MODULE PROCEDURE debug_dp_real_array1D
        MODULE PROCEDURE debug_dp_real_array2D
    END INTERFACE

    !> \brief Debug routines, used to print min and max values of an array variable with a predefined comment
    !! The goal is to make it possible to activate or deactivate the printing with a global flag
    !<
    INTERFACE DEBUG_MINMAX
        MODULE PROCEDURE debug_minmax_sp_real_array1D
        MODULE PROCEDURE debug_minmax_sp_real_array2D
        MODULE PROCEDURE debug_minmax_dp_real_array1D
        MODULE PROCEDURE debug_minmax_dp_real_array2D
        MODULE PROCEDURE debug_minmax_dp_real_array3D
    END INTERFACE DEBUG_MINMAX

    !> \brief Print value of variables with a predefined description
    !<
    INTERFACE print_string
        MODULE PROCEDURE print_string_val
        MODULE PROCEDURE print_string_vd
    END INTERFACE

    !> \brief Print value of variables with a predefined description
    !<
    INTERFACE PRINT_VAR
        MODULE PROCEDURE print_string_val
        MODULE PROCEDURE print_string_vd
        MODULE PROCEDURE print_integer
        MODULE PROCEDURE print_sp_real
        MODULE PROCEDURE print_dp_real
        MODULE PROCEDURE print_logical
        !MODULE PROCEDURE print_c_int_array1D
        MODULE PROCEDURE print_integer_array1D
        MODULE PROCEDURE print_integer_array2D
        MODULE PROCEDURE print_sp_real_array1D
        MODULE PROCEDURE print_sp_real_array2D
        MODULE PROCEDURE print_dp_real_array1D
        MODULE PROCEDURE print_dp_real_array2D
    END INTERFACE

    INTERFACE NOTHING
        MODULE PROCEDURE nothing_string
        MODULE PROCEDURE nothing_integer
        MODULE PROCEDURE nothing_sp_real
        MODULE PROCEDURE nothing_dp_real
        !MODULE PROCEDURE nothing_logical
        !MODULE PROCEDURE nothing_integer_array1D
        !MODULE PROCEDURE nothing_integer_array2D
        MODULE PROCEDURE nothing_sp_real_array1D
        MODULE PROCEDURE nothing_sp_real_array2D
        MODULE PROCEDURE nothing_dp_real_array1D
        MODULE PROCEDURE nothing_dp_real_array2D
        MODULE PROCEDURE nothing_isd_array1D!for the simulator
    END INTERFACE

    !> \brief Check for Nan values
    !! nan value can be check by this code
    !! ll_isnan = ( (val>0.0_dp).EQV.(val<=0.0_dp) )
    !! gfortran include the code to check it so no need to do it.
    !<
    INTERFACE isNaN2
        MODULE PROCEDURE isnan2_sp_real
        MODULE PROCEDURE isnan2_dp_real
        MODULE PROCEDURE isnan2_sp_real_array1D
        MODULE PROCEDURE isnan2_dp_real_array1D
    END INTERFACE

    !> \brief Count Nan values in an array
    !<
    INTERFACE NaNCount
        MODULE PROCEDURE nanCount_sp_real_array1D
        MODULE PROCEDURE nanCount_dp_real_array1D
        MODULE PROCEDURE nanCount_dp_real_array2D
    END INTERFACE

    !> \brief Check for Infinity values
    !<
    INTERFACE isInf
        MODULE PROCEDURE isinf_sp_real
        MODULE PROCEDURE isinf_dp_real
        MODULE PROCEDURE isinf_sp_real_array1D
        MODULE PROCEDURE isinf_dp_real_array1D
    END INTERFACE

    !> \brief Count Inf values in an array
    !<
    INTERFACE InfCount
        MODULE PROCEDURE infCount_sp_real_array1D
        MODULE PROCEDURE infCount_dp_real_array1D
    END INTERFACE

    INTEGER, PARAMETER, PRIVATE :: NB_DTAG = 8
    LOGICAL, DIMENSION (NB_DTAG), SAVE :: iga_dflags = (/&
        .TRUE.&  !> flag for debug message that are allways displayed
        , .TRUE.& !> overflow and undeflow flag
        , .TRUE.& !> flag for debug message with no tag
        , .FALSE.& !> io flag
        , .FALSE.& !> memory flag
        , .FALSE.& !> netcdf flag
        , .FALSE.& !> results flag
        , .FALSE.& !> trace flag
    /)
    ! global flags to control the debug messages
    INTEGER, PARAMETER :: &
        dALLWAYS  = 1& !>index of the debug massage that are allways displayed
        , dOVERFLOW = 2& !>index of the debug tag for overflow and undeflow results
        , dNOTAG    = 3& !>index of the debug with no tag
        , dIO       = 4& !>index of the debug tag for IO
        , dMEMORY   = 5& !>index of the debug tag for memory
        , dNETCDF   = 6& !>index of the debug tag for netcdf
        , dRESULT   = 7& !>index of the debug tag for result
        , dTRACE    = 8  !>index of the debug tag for trace, to control the trace at run
    INTEGER, PARAMETER :: dOVF = dOVERFLOW

CONTAINS

    !> convert integer to string
    FUNCTION itoa(id_val) RESULT(ala_val)
        INTEGER, INTENT(IN) :: id_val
        CHARACTER(LEN=ip_snl)   :: ala_val

        WRITE(ala_val, *) id_val
        ala_val = ADJUSTL(ala_val)
    END FUNCTION itoa

    !> convert single precision real to string
    FUNCTION rtoa(rd_val) RESULT(ala_val)
        REAL(sp), INTENT(IN) :: rd_val
        CHARACTER(LEN=ip_snl)   :: ala_val

        WRITE(ala_val, *) rd_val
        ala_val = ADJUSTL(ala_val)
    END FUNCTION rtoa

    !> convert double precision real to string
    FUNCTION dtoa(rd_val) RESULT(ala_val)
        REAL(dp), INTENT(IN) :: rd_val
        CHARACTER(LEN=ip_snl)   :: ala_val

        WRITE(ala_val, *) rd_val
        ala_val = ADJUSTL(ala_val)
    END FUNCTION dtoa

    !> \brief set the debug flags
    !! PARAM[in] notag flag of debug message with no tag
    !! PARAM[in] io flag of debug message tagged as IO
    !! PARAM[in] mem flag of debug message tagged as memory
    !! PARAM[in] nc flag of debug message tagged as netcdf
    !! PARAM[in] res flag of debug message tagged as result
    !! PARAM[in] tr flag of debug message tagged as trace
    !<
    SUBROUTINE set_debugFlag(notag, io, mem, nc, res, tr)
        LOGICAL, OPTIONAL, INTENT(IN) :: notag, io, mem, nc, res, tr
        IF ( PRESENT( notag) ) iga_dflags(dNOTAG) = notag
        IF ( PRESENT( io   ) ) iga_dflags(dIO) = io
        IF ( PRESENT( mem  ) ) iga_dflags(dMEMORY) = mem
        IF ( PRESENT( nc   ) ) iga_dflags(dNETCDF) = nc
        IF ( PRESENT( res  ) ) iga_dflags(dRESULT) = res
        IF ( PRESENT( tr   ) ) iga_dflags(dTRACE) = tr
    END SUBROUTINE set_debugFlag


    SUBROUTINE print_debugFlags()
        CALL debug(iga_dflags(dOVERFLOW), 'TRACE', tag=dALLWAYS)
        CALL debug(iga_dflags(dNOTAG), 'NOTAG', tag=dALLWAYS)
        CALL debug(iga_dflags(dIO), 'IO', tag=dALLWAYS)
        CALL debug(iga_dflags(dMEMORY), 'MEMORY', tag=dALLWAYS)
        CALL debug(iga_dflags(dNETCDF), 'NETCDF', tag=dALLWAYS)
        CALL debug(iga_dflags(dRESULT), 'RESULT', tag=dALLWAYS)
        CALL debug(iga_dflags(dTRACE), 'TRACE', tag=dALLWAYS)
    END SUBROUTINE print_debugFlags

    !> \brief load debug flags from namelist
    !!
    !<
    SUBROUTINE load_debugFlag(fName)
        INTEGER, PARAMETER           :: ip_numnam = 68
        CHARACTER(LEN=*), INTENT(IN) :: fName
        !local variables
        LOGICAL :: debug_default, debug_io, debug_memory, debug_netcdf, debug_result, debug_trace

        NAMELIST/NAM_DEBUG/&
        debug_default,&
        debug_io     ,&
        debug_memory ,&
        debug_netcdf ,&
        debug_result ,&
        debug_trace

        OPEN(ip_numnam, FILE=fName, FORM='FORMATTED', STATUS='OLD')
        READ(ip_numnam, NAM_DEBUG  )!reading the block NAM_DEBUG
        CLOSE(ip_numnam)
        CALL set_debugFlag(&
            notag = debug_default,&
            io    = debug_io     ,&
            mem   = debug_memory ,&
            nc    = debug_netcdf ,&
            res   = debug_result ,&
            tr    = debug_trace   &
        )
    END SUBROUTINE load_debugFlag

    !> \brief get the debug flag associated with a tag
    !!
    !<
    FUNCTION debug_flag(id_tag) RESULT(ll_flag)
        INTEGER, OPTIONAL, INTENT(IN) :: id_tag
        !local vars
        INTEGER :: il_tag
        LOGICAL :: ll_flag

        IF( PRESENT(id_tag) )THEN
            il_tag = MIN(id_tag, NB_DTAG)
        ELSE
            il_tag = dNOTAG
        END IF
        ll_flag = iga_dflags(il_tag)
    END FUNCTION debug_flag

    !>pause function
    SUBROUTINE dpause(ada_val)
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ada_val
        IF( PRESENT(ada_val) )THEN
            CALL print_var(' - press <Enter> key to resume', ada_val)
        ELSE
            CALL print_var('', 'Pause: press <Enter> key to resume')
        END IF
        READ(*,*)
    END SUBROUTINE dpause

    !> \brief Check if a simple precision value is NaN
    FUNCTION isnan2_sp_real(rd_val)RESULT(ll_isnan)
        REAL(sp), INTENT(IN) :: rd_val
        LOGICAL :: ll_isnan

        ll_isnan = isnan(rd_val)
    END FUNCTION isnan2_sp_real

    !> \brief Check if a double precision value is NaN
    FUNCTION isnan2_dp_real(rd_val)RESULT(ll_isnan)
        REAL(dp), INTENT(IN) :: rd_val
        LOGICAL :: ll_isnan

        ll_isnan = isnan(rd_val)
    END FUNCTION isnan2_dp_real

    !> \brief Check if a simple precision array is NaN
    FUNCTION isnan2_sp_real_array1D(rda_val)RESULT(lla_isnan)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        LOGICAL,  DIMENSION(SIZE(rda_val)) :: lla_isnan
        INTEGER :: ibi

        DO ibi=1,SIZE(rda_val)
            lla_isnan(ibi) = isnan(rda_val(ibi))
        END DO
    END FUNCTION isnan2_sp_real_array1D

    !> \brief Check if a double precision array is NaN
    FUNCTION isnan2_dp_real_array1D(rda_val)RESULT(lla_isnan)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        LOGICAL,  DIMENSION(SIZE(rda_val)) :: lla_isnan
        INTEGER :: ibi

        DO ibi=1,SIZE(rda_val)
            lla_isnan(ibi) = isnan(rda_val(ibi))
        END DO
    END FUNCTION isnan2_dp_real_array1D

    !> \brief Count NaN values in a simple precision array
    !!
    !<
    FUNCTION nanCount_sp_real_array1D(rda_val)RESULT(il_count)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: il_count

        il_count = count(isnan(rda_val))
    END FUNCTION nanCount_sp_real_array1D

    !> \brief Count NaN values in a double precision 1D array
    FUNCTION nanCount_dp_real_array1D(rda_val)RESULT(il_count)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: il_count

        il_count = count(isnan(rda_val))
    END FUNCTION nanCount_dp_real_array1D

    !> \brief Count NaN values in a double precision 2D array
    FUNCTION nanCount_dp_real_array2D(rda_val)RESULT(il_count)
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_val
        INTEGER :: il_count

        il_count = count(isnan(rda_val))
    END FUNCTION nanCount_dp_real_array2D

    !> \brief Check if a simple precision value is Inf
    FUNCTION isinf_sp_real(rd_val)RESULT(ll_isinf)
        REAL(sp), INTENT(IN) :: rd_val
        LOGICAL :: ll_isinf

        ll_isinf = ( (rd_val+1.0_sp)>=(rd_val+2.0_sp) )
    END FUNCTION isinf_sp_real

    !> \brief Check if a double precision value is Inf
    FUNCTION isinf_dp_real(rd_val)RESULT(ll_isinf)
        REAL(dp), INTENT(IN) :: rd_val
        LOGICAL :: ll_isinf

        ll_isinf = ( (rd_val+1.0_dp)>=(rd_val+2.0_dp) )
    END FUNCTION isinf_dp_real

    !> \brief Check if a simple precision array is Inf
    FUNCTION isinf_sp_real_array1D(rda_val)RESULT(lla_isinf)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        LOGICAL,  DIMENSION(SIZE(rda_val)) :: lla_isinf
        INTEGER :: ibi

        DO ibi=1,SIZE(rda_val)
            lla_isinf(ibi) = isinf( rda_val(ibi) )
        END DO
    END FUNCTION isinf_sp_real_array1D

    !> \brief Check if a double precision array is Inf
    FUNCTION isinf_dp_real_array1D(rda_val)RESULT(lla_isinf)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        LOGICAL,  DIMENSION(SIZE(rda_val)) :: lla_isinf
        INTEGER :: ibi

        DO ibi=1,SIZE(rda_val)
            lla_isinf(ibi) = isinf( rda_val(ibi) )
        END DO
    END FUNCTION isinf_dp_real_array1D

    !> \brief Count Inf values in a simple precision array
    FUNCTION infCount_sp_real_array1D(rda_val)RESULT(il_count)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: ibi, il_count

        il_count = 0
        DO ibi=1,SIZE(rda_val)
            IF ( isinf( rda_val(ibi) ) ) THEN
                il_count = il_count + 1
            END IF
        END DO
    END FUNCTION infCount_sp_real_array1D

    !> \brief Count Inf values in a simple precision array
    FUNCTION infCount_dp_real_array1D(rda_val)RESULT(il_count)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: ibi, il_count

        il_count = 0
        DO ibi=1,SIZE(rda_val)
            IF ( isinf( rda_val(ibi) ) ) THEN
                il_count = il_count + 1
            END IF
        END DO
    END FUNCTION infCount_dp_real_array1D

    !>\brief Write 2D array to file
    !<
    SUBROUTINE myio_write_matrix(fileName, rda_A)
        REAL(KIND=dp), DIMENSION(:, :), INTENT(IN) :: rda_A
        CHARACTER(LEN=*)         , INTENT(IN) :: fileName
        INTEGER ib_i, ib_j, il_ios
        INTEGER :: ip_fid = 41

        IF(fileName==SCREEN) THEN
            ip_fid = 6!standard output
            PRINT*, 'rows count = ', size(rda_A, 1)
            PRINT*, 'columns count = ', size(rda_A, 2)
        ELSE
            ip_fid = 41
            OPEN( UNIT = ip_fid, FILE = fileName, STATUS = 'REPLACE', FORM = 'FORMATTED', IOSTAT = il_ios)
            IF (il_ios /= 0) then
                WRITE(* , *) 'Error creating file', fileName
                STOP;
            END IF
        END IF
        PRINT*, "In writeMatrix", LBOUND(rda_A,1), UBOUND(rda_A,1), LBOUND(rda_A,2), UBOUND(rda_A,2)
        !stop
        WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 1)
        WRITE(unit=ip_fid, fmt=IFMT) size(rda_A, 2)
        DO ib_j = LBOUND(rda_A,2), UBOUND(rda_A,2)
            DO ib_i = LBOUND(rda_A,1), UBOUND(rda_A,1)
                WRITE(unit=ip_fid, fmt='(E18.8)', advance="no") rda_A(ib_i,ib_j)!"no" to avoid end of line
            END DO
            !WRITE(unit=ip_fid, fmt=*)''! just for adding end of line
        END DO
        IF(fileName/=SCREEN) THEN
            CLOSE(ip_fid )
        END IF
    END SUBROUTINE myio_write_matrix

    SUBROUTINE myio_read_matrix(fileName, rda_A)
        REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT) :: rda_A
        CHARACTER(LEN=*)         , INTENT(IN)  :: fileName
        INTEGER ib_i, ib_j, il_nbRow, il_nbCol, il_ios
        INTEGER , PARAMETER :: ip_fid =615

        OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
        IF (il_ios /= 0) then
            WRITE(* , *) 'Error opening file', fileName
            STOP;
        END IF
        !stop
        READ(UNIT=ip_fid, FMT=IFMT) il_nbRow
        READ(UNIT=ip_fid, FMT=IFMT) il_nbCol
        PRINT*, "In readMatrix, il_nbRow = ", il_nbRow, "; il_nbCol = ", il_nbCol
        DO ib_j = 1, il_nbCol
            DO ib_i = 1, il_nbRow
                READ(UNIT=ip_fid, FMT='(E18.8)', advance="no") rda_A(ib_i,ib_j)!"no" to avoid going to the next line
                PRINT*, 'row, col', ib_i, ib_j, ' : ', rda_A(ib_i,ib_j)
            END DO
            !READ(unit=ip_fid, fmt=*) !just to skip the end of line
        END DO
        CLOSE(ip_fid )
    END SUBROUTINE myio_read_matrix

    !!read the head of the file: nbRow and nbCol
    SUBROUTINE  myio_readInfo(fileName, id_nbRow, id_nbCol)
        INTEGER , PARAMETER :: ip_fid =616
        CHARACTER(LEN=*)         , INTENT(IN)  :: fileName
        INTEGER, INTENT(OUT) :: id_nbRow, id_nbCol
        INTEGER :: il_ios

        OPEN( UNIT = ip_fid, FILE = fileName, IOSTAT = il_ios)
        IF (il_ios /= 0) then
            WRITE(* , *) 'In readInfo : Error opening file ', fileName
            STOP;
        END IF
        READ(UNIT=ip_fid, FMT=IFMT) id_nbRow
        READ(UNIT=ip_fid, FMT=IFMT) id_nbCol
        PRINT*, "In readInfo, id_nbRow = ", id_nbRow, "; id_nbCol = ", id_nbCol
        CLOSE(ip_fid )
    END SUBROUTINE myio_readInfo

  !-!!!!!!!!!!!! stop_program routines

    !> \brief Prints  value and stops the program
    !! \param [in] description description of the calling poing and the cause of the stop
    !! \param [in] val value to be printed
    !<
    SUBROUTINE stop_program_with_string(val, description)
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: val
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
        !local variables
        integer :: l, k

        IF (PRESENT(description)) THEN
            l = len_trim(description)
            WRITE(*,*) description(1:l)
        ELSE
            WRITE(*,*) 'Something that the program can not manage, happened'
        END IF
        IF ( PRESENT(val) ) THEN
            k = len_trim(val)
            WRITE(*,*) 'the included message or problematic value is: ', val(1:k)
        ELSE
        END IF
        WRITE(*,*) '  Stopping ...'
        STOP
    END SUBROUTINE stop_program_with_string

    !> \brief Prints  value and stops the program
    !! \param [in] description description of the calling poing and the cause of the stop
    !! \param [in] val value to be printed
    !<
    SUBROUTINE stop_program_with_int(val, description)
        INTEGER, INTENT(IN) :: val
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

        CALL print_var(val, description)
        WRITE(*,*) '  Stopping the program ...'
        STOP
    END SUBROUTINE stop_program_with_int

    !> \brief Prints  value and stops the program
    !! \param [in] description description of the calling poing and the cause of the stop
    !! \param [in] val value to be printed
    !<
    SUBROUTINE stop_program_with_sp_real(val, description)
        REAL(sp), INTENT(IN) :: val
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

        CALL print_var(val, description)
        WRITE(*,*) '  Stopping the program ...'
        STOP
    END SUBROUTINE stop_program_with_sp_real

    !> \brief Prints  value and stops the program
    !! \param [in] description description of the calling poing and the cause of the stop
    !! \param [in] val value to be printed
    !<
    SUBROUTINE stop_program_with_dp_real(val, description)
        REAL(dp), INTENT(IN) :: val
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

        CALL print_var(val, description)
        WRITE(*,*) '  Stopping the program ...'
        STOP
    END SUBROUTINE stop_program_with_dp_real

    !print_var routines
    SUBROUTINE print_string_val(ada_val)
        CHARACTER(LEN=*), INTENT(IN) :: ada_val
        !
        WRITE(*, FMT="(A)")trim(ada_val)
    END SUBROUTINE print_string_val

    !print_var routines
    SUBROUTINE print_string_vd(ada_val, description)
        CHARACTER(LEN=*), INTENT(IN) :: ada_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        integer :: l, k

        l = len_trim(description)
        k = len_trim(ada_val)
        if( (l/=0).and.(k/=0) )then
            WRITE(*, FMT=DEBUG_AFORMAT) description(1:l), ada_val(1:k)
        else if(l/=0)then
            WRITE(*, FMT=DEBUG_AFORMAT) description(1:l), ''
        else if(k/=0)then
            WRITE(*, FMT=DEBUG_AFORMAT) '', ada_val(1:k)
        else
            WRITE(*, FMT="(A)")''
        end if
    END SUBROUTINE print_string_vd

    subroutine print_integer(val, description)
        integer, intent(in) :: val
        character(len=*), intent(in) :: description
        !local variables
        character(len=ip_fnl) :: dformat
        integer :: l

        l = len_trim(description)
        dformat = trim(pretty_dformat(val))
        if(l==0)then
            write(*, fmt=dformat) '-', val
        else
            write(*, fmt=dformat) description(1:l), val
        end if
    end subroutine print_integer

    SUBROUTINE print_sp_real(val, description)
        REAL(sp), INTENT(IN) :: val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        integer :: l

        l = len_trim(description)

        if(l==0)then
            write(*, fmt=DEBUG_RFORMAT) '-', val
        else
            WRITE(*, FMT=DEBUG_RFORMAT) description(1:l), val
        end if
        !call dpause("end of print_sp_real")
    END SUBROUTINE print_sp_real

    SUBROUTINE print_dp_real(val, description)
        REAL(dp), INTENT(IN) :: val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        integer :: l

        l = len_trim(description)

        if(l==0)then
            write(*, fmt=DEBUG_RFORMAT) '-', val
        else
            WRITE(*, FMT=DEBUG_RFORMAT) description(1:l), val
        end if
    END SUBROUTINE print_dp_real

    SUBROUTINE print_logical(val, description)
        LOGICAL, INTENT(IN) :: val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        integer :: l

        l = len_trim(description)
        if(l==0)then
            write(*, fmt=DEBUG_LFORMAT) '-', val
        else
            WRITE(*, FMT=DEBUG_LFORMAT) description(1:l), val
        end if
    END SUBROUTINE print_logical

    SUBROUTINE print_c_int_array1D(ida_val, description)
        INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        CHARACTER(LEN=80)   :: array_format
        integer :: vsize, l

        l = len_trim(description)
        vsize = SIZE(ida_val)
        if(vsize>0) then
            array_format = "(A,"//TRIM( NUM2STR(SIZE(ida_val)) )//TRIM(pretty_format( ida_val) )//")"
            WRITE(*, FMT=array_format) description(1:l), ida_val
        else
            call debug("zero size array", description, tag=dALLWAYS)
        end if
    END SUBROUTINE print_c_int_array1D

    SUBROUTINE print_integer_array1D(ida_val, description)
        INTEGER, DIMENSION(:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        CHARACTER(LEN=80)   :: array_format
        integer :: vsize, l

        l = len_trim(description)
        vsize = SIZE(ida_val)
        if(vsize>0) then
            array_format = "(A,': ',"//TRIM( NUM2STR(SIZE(ida_val)) )//TRIM(pretty_format( ida_val) )//")"
            WRITE(*, FMT=array_format) description(1:l), ida_val
        else
            call debug("zero size array", description, tag=dALLWAYS)
        end if
    END SUBROUTINE print_integer_array1D

    SUBROUTINE print_sp_real_array1D(rda_val, description)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        CHARACTER(LEN=80)   :: array_format
        integer :: vsize, l

        l = len_trim(description)
        vsize = SIZE(rda_val)
        if(vsize>0) then
            array_format = "(A,': ',"//TRIM( NUM2STR(SIZE(rda_val)) )//RFORMAT//")"
            WRITE(*, FMT=array_format) description(1:l), rda_val
        else
            call debug("zero size array", description, tag=dALLWAYS)
        end if
    END SUBROUTINE print_sp_real_array1D

    SUBROUTINE print_dp_real_array1D(rda_val, description)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        !local variables
        CHARACTER(LEN=80)   :: array_format
        integer :: vsize, l

        l = len_trim(description)
        vsize = SIZE(rda_val)
        if(vsize>0) then
            array_format = "(A,': ',"//TRIM( NUM2STR(vsize) )//RFORMAT//")"
            WRITE(*, FMT=array_format) description(1:l), rda_val
        else
            call debug("zero size array", description, tag=dALLWAYS)
        end if
    END SUBROUTINE print_dp_real_array1D

    FUNCTION pretty_format_int_scalar(id_val)RESULT(str_format)
        INTEGER, INTENT(IN) :: id_val
        CHARACTER(LEN=80)   :: str_format
        INTEGER :: il_npos, absval

        absval = abs(id_val)
        if (absval>=10)then
        il_npos = min(floor( log10(real(absval)) ), 10)
        else
        il_npos = 1
        end if
        !+2 for the sign and 1 space
        str_format = "(I"//TRIM( NUM2STR(il_npos+2) )//")"
    END FUNCTION pretty_format_int_scalar

    FUNCTION pretty_format_int_1D_array(ida_val)RESULT(str_format)
        INTEGER, DIMENSION(:), INTENT(IN) :: ida_val
        CHARACTER(LEN=80)   :: str_format

        str_format = pretty_format_int_scalar(maxval(abs(ida_val)))
    END FUNCTION pretty_format_int_1D_array

    FUNCTION pretty_format_int_2D_array(ida_val)RESULT(str_format)
        INTEGER, DIMENSION(:,:), INTENT(IN) :: ida_val
        CHARACTER(LEN=80)   :: str_format

        str_format = pretty_format_int_scalar(maxval(abs(ida_val)))
    END FUNCTION pretty_format_int_2D_array

    FUNCTION pretty_dformat_int_scalar(id_val)RESULT(str_format)
        INTEGER, INTENT(IN) :: id_val
        CHARACTER(LEN=80)   :: str_format

        str_format = "(A,': ',"//TRIM(pretty_format( id_val))//")"
    END FUNCTION pretty_dformat_int_scalar

    SUBROUTINE print_integer_array2D(ida_val, description)
        INTEGER, DIMENSION(:,:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        CHARACTER(LEN=80)   :: array_format
        INTEGER :: ibi

        array_format = "("//TRIM( NUM2STR(SIZE(ida_val,2)) )//TRIM(pretty_format( ida_val))//")"
        write(*,"(A,': ')")description
        DO ibi = 1, SIZE(ida_val,1)
        WRITE(*, FMT=array_format) ida_val(ibi, :)
        END DO
    END SUBROUTINE print_integer_array2D

    SUBROUTINE print_sp_real_array2D(rda_val, description)
        REAL(sp), DIMENSION(:,:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        CHARACTER(LEN=80)   :: array_format
        INTEGER :: ibi

        array_format = "("//TRIM( NUM2STR(SIZE(rda_val,2)) )//RFORMAT//")"
        write(*,"(A,': ')")description
        DO ibi = 1, SIZE(rda_val,1)
        WRITE(*, FMT=array_format), rda_val(ibi, :)
        END DO
    END SUBROUTINE print_sp_real_array2D

    subroutine print_dp_real_array2D(rda_val, description)
        real(dp), dimension(:,:), intent(in) :: rda_val
        character(len=*), intent(in) :: description
        character(len=80)   :: array_format
        integer :: ibi

        array_format = "("//trim( num2str(size(rda_val,2)) )//RFORMAT//")"
        write(*,"(A,': ')")description
        !print*,   "print_dp_real_array2D, format: ", array_format
        do ibi = 1, size(rda_val,1)
            write(*, fmt=array_format), rda_val(ibi, :)
        end do
    end subroutine print_dp_real_array2D

    !> @brief debug data and time
    subroutine debugdt(description)
        character(len=*), intent(in) :: description
        !local variables
        character(len=ip_snl):: ala_dt, mdate, mtime

        call date_and_time(DATE=mdate, TIME=mtime)!, ZONE, VALUES
        !DATE has form ccyymmdd. TIME has form hhmmss.sss
        ala_dt = mdate(1:4)//'/'//mdate(5:6)//'/'//mdate(7:8)//&
                 '--'//mtime(1:2)//':'//mtime(3:4)//':'//mtime(5:6)//mtime(7:10)

        call print_string(ala_dt, description)
    end subroutine debugdt

    !debug routines
    SUBROUTINE debug_string_val(ada_val, tag)
        CHARACTER(LEN=*), INTENT(IN) :: ada_val
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_string(ada_val)
    END SUBROUTINE debug_string_val

    SUBROUTINE debug_string(ada_val, description, tag)
        CHARACTER(LEN=*), INTENT(IN) :: ada_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_string(ada_val, description)
    END SUBROUTINE debug_string

    SUBROUTINE debug_integer(id_val, description, tag)
        INTEGER, INTENT(IN) :: id_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_integer(id_val, description)
    END SUBROUTINE debug_integer

    SUBROUTINE debug_sp_real(rd_val, description, tag)
        REAL(sp), INTENT(IN) :: rd_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_sp_real(rd_val, description)
    END SUBROUTINE debug_sp_real

    SUBROUTINE debug_dp_real(rd_val, description, tag)
        REAL(dp), INTENT(IN) :: rd_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_dp_real(rd_val, description)
    END SUBROUTINE debug_dp_real

    SUBROUTINE debug_logical(ld_val, description, tag)
        LOGICAL, INTENT(IN) :: ld_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_logical(ld_val, description)
    END SUBROUTINE debug_logical

    SUBROUTINE debug_c_int_array1D(ida_val, description, tag)
        INTEGER(C_INT), DIMENSION(:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_c_int_array1D(ida_val, description)
    END SUBROUTINE debug_c_int_array1D

    SUBROUTINE debug_integer_array1D(ida_val, description, tag)
        INTEGER, DIMENSION(:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_integer_array1D(ida_val, description)
    END SUBROUTINE debug_integer_array1D

    SUBROUTINE debug_sp_real_array1D(rda_val, description, tag)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_sp_real_array1D(rda_val, description)
    END SUBROUTINE debug_sp_real_array1D

    SUBROUTINE debug_dp_real_array1D(rda_val, description, tag)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_dp_real_array1D(rda_val, description)
    END SUBROUTINE debug_dp_real_array1D

    SUBROUTINE debug_integer_array2D(ida_val, description, tag)
        INTEGER, DIMENSION(:,:), INTENT(IN) :: ida_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_integer_array2D(ida_val, description)
    END SUBROUTINE debug_integer_array2D

    SUBROUTINE debug_sp_real_array2D(rda_val, description, tag)
        REAL(sp), DIMENSION(:,:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_sp_real_array2D(rda_val, description)
    END SUBROUTINE debug_sp_real_array2D

    SUBROUTINE debug_dp_real_array2D(rda_val, description, tag)
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_val
        CHARACTER(LEN=*), INTENT(IN) :: description
        INTEGER, OPTIONAL, INTENT(IN) :: tag

        IF( debug_flag(tag) ) CALL print_dp_real_array2D(rda_val, description)
    END SUBROUTINE debug_dp_real_array2D

    subroutine debug_minmax_sp_real_array1D(rda_val, description, tag)
        real(sp), dimension(:), intent(in) :: rda_val
        character(len=*), intent(in)  :: description
        integer, optional, intent(in) :: tag
        !local variables
        real(sp) :: mn, mx
        integer :: l

        if( debug_flag(tag) )then
            mn = minval(rda_val)
            mx = maxval(rda_val)
            l = len_trim(description)
            print*, 'DEBUG_MMFORMAT = ', DEBUG_MMFORMAT
            write(*, DEBUG_MMFORMAT) description(1:l), mn, mx
        end if
    end subroutine debug_minmax_sp_real_array1D

    subroutine debug_minmax_dp_real_array1D(rda_val, description, tag)
        real(dp), dimension(:), intent(in) :: rda_val
        character(len=*), intent(in) :: description
        integer, optional, intent(in) :: tag
        !local variables
        real(dp) :: mn, mx
        integer :: l

        if( debug_flag(tag) )then
            mn = minval(rda_val)
            mx = maxval(rda_val)
            l = len_trim(description)
            write(*, DEBUG_MMFORMAT) description(1:l), mn, mx
        end if
    end subroutine debug_minmax_dp_real_array1D

    subroutine debug_minmax_sp_real_array2D(rda_val, description, tag)
        real(sp), dimension(:,:), intent(in) :: rda_val
        character(len=*), intent(in)  :: description
        integer, optional, intent(in) :: tag
        !local variables
        real(sp) :: mn, mx
        integer :: l

        if( debug_flag(tag) )then
            mn = minval(rda_val)
            mx = maxval(rda_val)
            l = len_trim(description)
            write(*, DEBUG_MMFORMAT) description(1:l), mn, mx
        end if
    end subroutine debug_minmax_sp_real_array2D

    subroutine debug_minmax_dp_real_array2D(rda_val, description, tag)
        real(dp), dimension(:,:), intent(in) :: rda_val
        character(len=*), intent(in)  :: description
        integer, optional, intent(in) :: tag
        !local variables
        real(dp) :: mn, mx
        integer :: l

        if( debug_flag(tag) )then
            mn = minval(rda_val)
            mx = maxval(rda_val)
            l = len_trim(description)
            write(*, DEBUG_MMFORMAT) description(1:l), mn, mx
        end if
    end subroutine debug_minmax_dp_real_array2D

    subroutine debug_minmax_dp_real_array3D(rda_val, description, tag)
        real(dp), dimension(:,:,:), intent(in) :: rda_val
        character(len=*), intent(in)  :: description
        integer, optional, intent(in) :: tag
        !local variables
        real(dp) :: mn, mx
        integer :: l

        if( debug_flag(tag) )then
            mn = minval(rda_val)
            mx = maxval(rda_val)
            l = len_trim(description)
            write(*, DEBUG_MMFORMAT) description(1:l), mn, mx
        end if
    end subroutine debug_minmax_dp_real_array3D

    !NOTHINg interface
    SUBROUTINE nothing_string(ala_val)
        CHARACTER(LEN=*), INTENT(IN) :: ala_val
        INTEGER :: il_tmp
        il_tmp = LEN(ala_val)
    END SUBROUTINE nothing_string

    SUBROUTINE nothing_integer(id_val)
        INTEGER, INTENT(IN) :: id_val
        INTEGER :: il_tmp
        il_tmp = id_val
    END SUBROUTINE nothing_integer

    SUBROUTINE nothing_sp_real(rd_val)
        REAL(sp), INTENT(IN) :: rd_val
        REAL(dp) :: rl_tmp
        rl_tmp = REAL(rd_val, dp)
    END SUBROUTINE nothing_sp_real

    SUBROUTINE nothing_dp_real(rd_val)
        REAL(dp), INTENT(IN) :: rd_val
        REAL(dp) :: rl_tmp
        rl_tmp = REAL(rd_val, dp)
    END SUBROUTINE nothing_dp_real

    SUBROUTINE nothing_sp_real_array1D(rda_val)
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: il_tmp
        il_tmp = SIZE(rda_val, 1)
    END SUBROUTINE nothing_sp_real_array1D

    SUBROUTINE nothing_sp_real_array2D(rda_val)
        REAL(sp), DIMENSION(:,:), INTENT(IN) :: rda_val
        INTEGER :: il_tmp
        il_tmp = SIZE(rda_val, 1)
    END SUBROUTINE nothing_sp_real_array2D

    SUBROUTINE nothing_dp_real_array1D(rda_val)
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_val
        INTEGER :: il_tmp
        il_tmp = SIZE(rda_val, 1)
    END SUBROUTINE nothing_dp_real_array1D

    SUBROUTINE nothing_dp_real_array2D(rda_val)
        REAL(dp), DIMENSION(:,:), INTENT(IN) :: rda_val
        INTEGER :: il_tmp
        il_tmp = SIZE(rda_val, 1)
    END SUBROUTINE nothing_dp_real_array2D

    SUBROUTINE nothing_isd_array1D(ida_int, rda_sp, rda_dp)
        INTEGER, DIMENSION(:) :: ida_int
        REAL(dp), DIMENSION(:), INTENT(IN) :: rda_dp
        REAL(sp), DIMENSION(:), INTENT(IN) :: rda_sp
        INTEGER :: il_tmp
        il_tmp = SIZE(rda_dp, 1)
        il_tmp = SIZE(rda_sp, 1)
        il_tmp = SIZE(ida_int, 1)
    END SUBROUTINE nothing_isd_array1D
END MODULE debug_tools
