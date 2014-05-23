!> \file sort_test.f90
!! \brief 
!! @author Innocent Souopgui
!!
!<
PROGRAM SORT
USE cs_tools
USE debug_tools
USE general_tools
  
  INTEGER, PARAMETER :: ip_n = 100,&
                        ip_m = 10
  INTEGER, DIMENSION(ip_m) :: ida_vec
  
  CALL init_all()
  CALL random_indexes(ida_vec, ip_n)
  CALL debug(ida_vec, 'Before sorting: ')
  CALL select_sort(ida_vec)
  CALL debug(ida_vec, 'After  sorting: ')
  
  CALL finalize_all()
CONTAINS

  SUBROUTINE init_all()
    !CALL init_fftw()
    CALL init_cs_tools()
  END SUBROUTINE init_all
  
  SUBROUTINE finalize_all()
    !CALL finalize_fftw()
    CALL finalize_cs_tools()
  END SUBROUTINE finalize_all
END PROGRAM SORT