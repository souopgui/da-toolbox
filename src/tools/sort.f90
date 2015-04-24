!> \file sort.f90
!! \brief defines sorting routines
!! @author Innocent Souopgui
!!
!<
MODULE sort_tools
IMPLICIT NONE
  !> \brief sort routines, selection
  !<
  INTERFACE select_sort
     MODULE PROCEDURE select_sort_int
     !MODULE PROCEDURE select_sort_sp_real
     !MODULE PROCEDURE select_sort_dp_real
  END INTERFACE
  !> \brief switch routines
  !<
  INTERFACE switch
     MODULE PROCEDURE switch_int
     !MODULE PROCEDURE switch_sp_real
     !MODULE PROCEDURE switch_dp_real
  END INTERFACE
CONTAINS
  SUBROUTINE select_sort_int(ida_x)
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ida_x
    INTEGER :: ibi, ibj, il_min_idx, il_size

    il_size = SIZE(ida_x)
    DO ibi = 1, il_size-1
      il_min_idx = ibi
      DO ibj = ibi+1, il_size
        IF( ida_x(ibj)<ida_x(il_min_idx) ) il_min_idx = ibj
      END DO
      IF(ibi /= il_min_idx) THEN
        CALL switch( ida_x(ibi), ida_x(il_min_idx) )
      END IF
    END DO
  END SUBROUTINE select_sort_int

  SUBROUTINE switch_int(id_a, id_b)
    INTEGER, INTENT(INOUT) :: id_a, id_b
    !local variables
    INTEGER :: il_tmp

    il_tmp = id_a
    id_a = id_b
    id_b = il_tmp
  END SUBROUTINE switch_int
END MODULE sort_tools