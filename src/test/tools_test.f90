PROGRAM tools_test
  USE general_constant
  USE general_tools
  USE com_tools
  USE debug_tools
IMPLICIT NONE
  INTEGER, PARAMETER :: NDIM_MAX = 5
  INTEGER, DIMENSION(NDIM_MAX), PARAMETER :: IPA_SHAPE = (/101, 200, 100, 100, 100/)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ila_coord
  INTEGER, DIMENSION(:)  , ALLOCATABLE :: ila_step, ila_max, ila_shift, ila_ncoord
  INTEGER :: il_ndim, il_ncoord

  il_ndim = 3
  ALLOCATE( ila_step(il_ndim), ila_max(il_ndim), ila_shift(il_ndim), ila_ncoord(il_ndim) )
  ila_max = IPA_SHAPE(1:il_ndim)
  ila_step = 50
  il_ncoord = regular_coord_count(ila_max, steps =ila_step)
  ALLOCATE( ila_coord(il_ndim, il_ncoord) )
  CALL build_regular_int_coord(ila_coord, max_coords=ila_max, steps=ila_step)
  CALL debug(il_ncoord, 'il_ncoord = ')
  CALL debug(ila_coord, 'ila_coord = ')
END PROGRAM tools_test