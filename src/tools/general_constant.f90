!> \file general_constant.f90
!! \brief Global Constant Definition
!! @author Innocent Souopgui
!!
!<
MODULE general_constant
IMPLICIT NONE
  !> \brief Unknown integer value
  !<
  INTEGER, PARAMETER :: IP_UNKNOWN = -999!!unknown value
  !> \brief size of simple precision real variables
  !<
  INTEGER, PARAMETER :: sp = KIND(1.0)
  !> \brief size of double precision real variables
  !<
  INTEGER, PARAMETER :: dp = KIND(1D0)
  !> \brief size of the default precision for real variables, stand from computation precision
  !<
  INTEGER, PARAMETER :: cp = dp
  !> \brief size of the default precision for graphic, (dislin), stand from graphic precision
  !<
  INTEGER, PARAMETER :: gp = sp

  !>
  !! \enum noname
  !! \brief String Comparison Member Function status.
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

  !> \brief pi constant
  Real(cp), Parameter :: rp_pi = 3.14159265358979323846_cp
  !> \brief small number for real values comparison, used as zero
  Real(cp), Parameter :: rp_epsilon = REAL(1d-15, cp)
  !> \brief Strings length, used for short file name or other short string variables
  !<
  INTEGER, PARAMETER :: ip_snl = 80
  !> \brief maximum length of the string representation of common constants
  !<
  INTEGER, PARAMETER :: ip_ccl = 20
  !> \brief Strings length, used for full file name or other short string variables
  !<
  INTEGER, PARAMETER :: ip_fnl = 255
  !> \brief format to print integer with comment
  !<
  CHARACTER(LEN=*), PARAMETER :: PRINT_IFORMAT = "(A,I5)"
  !> \brief format to print real with comment
  !<
  CHARACTER(LEN=*), PARAMETER :: PRINT_RFORMAT = "(A,E10.2)"
  !> \brief old status for file input/output
  !<
  CHARACTER(LEN=*), PARAMETER :: FILE_OLD = 'OLD'
  !> \brief replaced status for file input/output
  !<
  CHARACTER(LEN=*), PARAMETER :: FILE_REPLACE = 'REPLACE'
  !> \brief format for file input/output
  !<
  CHARACTER(LEN=*), PARAMETER :: FILE_FORMATTED = 'FORMATTED'
  !> \brief string delimiter
  !! \details This delimiter is used when saving string or characters data into namelist
  !! Usefull for special characters in namelist, for example \ and &
  !<
  CHARACTER(LEN=*), PARAMETER :: SD= "'"
  !> \brief constant for array vectorization, concatenation of rows(row by row)
  INTEGER, PARAMETER :: ROW_WISE =1341
  !> \brief constant for array vectorization, concatenation of columns(column by column)
  INTEGER, PARAMETER :: COL_WISE =1342

  !> \brief standard output file id
  !<
  INTEGER, PARAMETER :: STDOUT = 6

END MODULE general_constant
