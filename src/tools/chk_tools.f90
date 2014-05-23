!> \file chk_tools.f90
!! \brief checking module for the tools module
!! @author Innocent Souopgui
!! \details
!<
MODULE chk_tools
  USE general_constant
IMPLICIT NONE

CONTAINS

  !> \brief check argument for the function/subroutine regular_coord_count
  !! \see regular_coord_count
  !<
  SUBROUTINE chk_regular_coord_count(ida_step, ida_max, ida_min, ida_ncoord)
    INTEGER, DIMENSION(:), INTENT(IN) :: ida_step, ida_max
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)    :: ida_min
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: ida_ncoord
    INTEGER, DIMENSION(SIZE(ida_max)) :: ila_min, ila_count
    LOGICAL :: ll_conform_arg !> used to check the conformity of arguments

    IF( PRESENT(ida_min) ) THEN
      ll_conform_arg = ( (SIZE(ida_step)==SIZE(ida_max)).AND.(SIZE(ida_step)==SIZE(ida_min)) )
    ELSE
      ll_conform_arg = ( SIZE(ida_step)==SIZE(ida_max) )
    END IF
    IF( PRESENT(ida_ncoord) ) THEN
      ll_conform_arg = ( ll_conform_arg.AND.(SIZE(ida_step)==SIZE(ida_ncoord)) )
    END IF
    !Check size
    IF(.NOT.ll_conform_arg)THEN
      CALL stop_progam('', 'In chk_tools::regular_coord_count: different sizes for arrays ida_step, ida_max [, ida_min] [and ida_ncoord] ')
    END IF
  END SUBROUTINE chk_regular_coord_count
END MODULE chk_tools