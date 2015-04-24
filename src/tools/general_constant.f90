!> @file general_constant.f90
!! @brief Global Constant Definition
!! @author Innocent Souopgui
!!
!<

!>
!!
!<
module general_constant
    use compiler_dependant
implicit none
    !> @brief Unknown integer value
    !<
    integer, parameter :: ip_unknown = -999!!unknown value
    !> @brief size of simple precision real variables
    !<
    integer, parameter :: sp = kind(1.0)
    !> @brief size of double precision real variables
    !<
    integer, parameter :: dp = kind(1d0)
    !> @brief size of the default precision for real variables, stand from computation precision
    !<
    integer, parameter :: cp = dp
    !   !> @brief size of the default precision for graphic, (dislin), stand from graphic precision
    !   !<
    !   integer, parameter :: gp = sp

    !>
    !! \enum noname
    !! @brief String Comparison Member Function status.
    !!
    !! \var ASC explanation 1
    !! \var DESC explanation 2
    !<
    enum, bind(c)
        enumerator :: ASC , DESC !< descending order
        enumerator :: ASCENDING = ASC, DESCENDING = DESC
    end enum

    enum, bind(c)
        enumerator :: MASKED = 1 !< masking mark
        enumerator :: UNMASKED = 0 !<
    end enum

    !> @brief pi constant
    real(cp), parameter :: rp_pi = 3.14159265358979323846_cp
    !> @brief small number for real values comparison, used as zero
    real(cp), parameter :: rp_epsilon = REAL(1d-15, cp)
    !> @brief Strings length, used for short file name or other short string variables
    !<
    integer, parameter :: ip_snl = 80
    !> @brief maximum length of the string representation of common constants
    !<
    integer, parameter :: ip_ccl = 20
    !> @brief Strings length, used for full file name or other short string variables
    !<
    integer, parameter :: ip_fnl = 255
    !> @brief format to print integer with comment
    !<
    character(len=*), parameter :: PRINT_IFORMAT = "(A,I5)"
    !> @brief format to print real with comment
    !<
    character(len=*), parameter :: PRINT_RFORMAT = "(A,E10.2)"
    !> @brief old status for file input/output
    !<
    character(len=*), parameter :: FILE_OLD = 'OLD'
    !> @brief replaced status for file input/output
    !<
    character(len=*), parameter :: FILE_REPLACE = 'REPLACE'
    !> @brief format for file input/output
    !<
    character(len=*), parameter :: FILE_FORMATTED = 'FORMATTED'
    !> @brief string delimiter
    !! @details This delimiter is used when saving string or characters data into namelist
    !! Usefull for special characters in namelist, for example \ and &
    !<
    character(len=*), parameter :: SD= "'"
    !> @brief constant for array vectorization, concatenation of rows(row by row)
    integer, parameter :: ROW_WISE =1341
    !> @brief constant for array vectorization, concatenation of columns(column by column)
    integer, parameter :: COL_WISE =1342

    !> @brief standard output file id
    !<
    integer, parameter :: STDOUT = 6

end module general_constant
