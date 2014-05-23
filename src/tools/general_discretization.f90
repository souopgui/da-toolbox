!> \file general_discretization.f90
!! Defines the general discretization module
!<

!> Module for the description of discretization
!!@author Innocent Souopgui
MODULE general_discretization
  USE general_constant
  USE debug_tools
  IMPLICIT NONE


  !> \brief User defined type for for the description of a discrete subdomain (2D)
  !! Defines parameters that are used to define or manage a discrete subdomain
  !<
  TYPE rectangular_subdomain
    LOGICAL :: l_allocated = .FALSE.
    INTEGER :: i_ndim = -1  !< dimensionality of the physical space
    INTEGER, DIMENSION(:), POINTER :: ia_shift => NULL()!<shift relative to the full domain
    INTEGER, DIMENSION(:), POINTER :: ia_size  => NULL()
  END TYPE
CONTAINS

  !> \brief print the parameters of a rectangular subdomain
  !! \param [in, out] td_rs data structure for the rectangular subdomain to be printed
  !<
  SUBROUTINE print_rs(td_rs)
    TYPE(rectangular_subdomain), INTENT(IN) :: td_rs

    CALL debug("Printing rectangular_subdomain -------------------", tag=dALLWAYS)
    CALL debug(td_rs%i_ndim  , "  problem dimension     = ", tag=dALLWAYS)
    CALL debug(td_rs%ia_shift, "  indexes shift         = ", tag=dALLWAYS)
    CALL debug(td_rs%ia_size , "  size of the subdomain = ", tag=dALLWAYS)
    CALL debug("--------------------------------------------------", tag=dALLWAYS)
  END SUBROUTINE print_rs

  !> \brief Sets the parameters of a rectangular subdomain
  !! \param [in, out] td_rs data structure for the rectangular subdomain to be set
  !! \param [in] ida_shift shift of the subdomaain data (relative to the entire physical domain)
  !! \param [in] ida_size size of the subdomain
  !<
  SUBROUTINE set_rs(td_rs, ida_shift, ida_size)
    TYPE(rectangular_subdomain), INTENT(IN OUT) :: td_rs
    INTEGER, DIMENSION(:), INTENT(IN) :: ida_shift, ida_size
    !local variables
    INTEGER il_size

    il_size = SIZE(ida_shift)
    IF( il_size==SIZE(ida_size) )THEN
      CALL resize_rs(td_rs, il_size)
      td_rs%ia_shift  = ida_shift
      td_rs%ia_size = ida_size
    ELSE
      CALL stop_program("In set_rs, incoherence, the size of ida_shift is diffrent from the size of ida_size")
    END IF
  END SUBROUTINE set_rs

  !> \brief Resizes data structure for rectangular subdomain
  !! \param [in, out] td_rs data structure for the rectangular subdomain to be set
  !! \param [in] id_ndim number of dimension of the physical space
  !<
  SUBROUTINE resize_rs(td_rs, id_ndim)
    TYPE(rectangular_subdomain), INTENT(IN OUT) :: td_rs
    INTEGER, INTENT(IN) :: id_ndim
    !local variables
    LOGICAL :: ll_allocate

    ll_allocate = .TRUE.
    IF(td_rs%l_allocated)THEN
      IF(id_ndim/=td_rs%i_ndim)THEN
        DEALLOCATE(td_rs%ia_shift, td_rs%ia_size)
        td_rs%i_ndim = 0
      ELSE
        ll_allocate = .FALSE.
      END IF
    END IF
    IF(ll_allocate.AND.(id_ndim>0))THEN
      ALLOCATE( td_rs%ia_shift(id_ndim), td_rs%ia_size(id_ndim) )
      td_rs%l_allocated = .TRUE.
    END IF
    td_rs%i_ndim = id_ndim
  END SUBROUTINE resize_rs

END MODULE general_discretization