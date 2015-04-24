!> @file cs_tools.f90
!! @brief general tools for compressive sensing
!! @author Innocent Souopgui
!! @details
!! generating random indexes and random coefficients for tests
!<
MODULE cs_tools
  USE general_constant
  USE debug_tools
  USE general_tools
IMPLICIT NONE

  PRIVATE init_random_seed

CONTAINS

  !> @brief generate random values and store at random indexes of a vector
  !! @param [in,out] rda_x array containing the generated values
  !! @param [in] rd_nz_ratio ratio of nonzero coefficients to be generated
  !! @param [in] rd_maxval (Optional) if present, specifies the maximum value to be generated
  !! @param [in] ida_idx (Optional) if present, contains the indexes of nonzero element of rda_x
  !! The routine first generates id_nval non redundant indexes in the range 1..size(rda_x), and then generates id_nval random numbers and their sign, and stores at generated indexes
  !! the system random generator is used first and the generated values are scaled with rd_maxval if present
  !<
  SUBROUTINE randval_randidx(rda_x, rd_nz_ratio, rd_maxval, ida_idx)
    REAL(KIND=cp), INTENT(in) :: rd_nz_ratio
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)   :: rda_x
    REAL(KIND=cp), OPTIONAL, INTENT(in)          :: rd_maxval
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: ida_idx

    !local variables
    REAL(KIND=cp) :: rl_maxval
    INTEGER :: il_nonzero, il_size
    INTEGER, DIMENSION( : ), ALLOCATABLE :: ila_idx, ila_sign
    REAL(KIND=cp), DIMENSION( : ), ALLOCATABLE :: rla_tmp

    IF( PRESENT(rd_maxval) )THEN
        rl_maxval = rd_maxval
    ELSE
        rl_maxval = 1.0_cp
    END IF
    rda_x = 0.0_cp
    il_size = SIZE(rda_x)
    IF( (rd_nz_ratio > 0.0).AND.(rd_nz_ratio < 1.0) )THEN
      il_nonzero = CEILING(rd_nz_ratio*il_size)
    ELSE
      CALL stop_program(rd_nz_ratio, 'In randval_randidx, bad ratio (]0, 1[) :')
    END IF
    ALLOCATE( ila_idx(il_nonzero), ila_sign(il_nonzero), rla_tmp(il_nonzero) )
    !generating indexes
    CALL random_indexes(ila_idx, il_size)
    !generating sign, sign is considered negative if >0.5, remember that rand number are in the intervalle [0-1]
    CALL RANDOM_NUMBER(rla_tmp)
    ila_sign = 1
    WHERE( rla_tmp>0.5_cp ) ila_sign = -1
    !generatimg values
    CALL RANDOM_NUMBER(rla_tmp)
    !scaling
    rda_x( ila_idx ) = rl_maxval*rla_tmp*REAL(ila_sign, cp)
    IF (PRESENT(ida_idx)) ida_idx(1:il_nonzero) = ila_idx

  END SUBROUTINE randval_randidx

  !>\brief generate a non redundant array of random integer values, each number being less or equal to id_max
  !!@param [in] id_max max value
  !!@param [in,out] ida_idx array containing the generated numbers
  !!@param [in] stat return status, 0 for good and other thing for bad
  !!\details This routine is use to generate a non redundant array of integer values.
  !!These values are intended to be uses as random indexes in compressive sensing.
  !! random numbeer generator is supposed to be initialized
  !<
  SUBROUTINE random_indexes(ida_idx, id_max, stat)
    INTEGER, INTENT(in) :: id_max
    INTEGER, DIMENSION(:), INTENT(INOUT) :: ida_idx
    INTEGER, OPTIONAL, INTENT(OUT) :: stat
    !local variables
    INTEGER :: il_next, il_total, il_rand, il_stat
    REAL(KIND=cp) :: rl_rand
    il_total = SIZE(ida_idx)

    il_stat = 0
    IF( il_total > id_max )THEN!size of the array bigger than the max idx
      IF( .NOT.PRESENT(stat) )THEN
        CALL stop_program('', 'In random_indexes; the number of indexes is greatter than the max index')
      END IF
      il_stat = -1
    ELSE
      ida_idx = -1
      il_next = 1
      DO WHILE(il_next <= il_total)
        CALL RANDOM_NUMBER(rl_rand)
        il_rand = MOD(CEILING( rl_rand*id_max ), id_max) + 1
        IF ( .NOT.( ANY(ida_idx==il_rand) ) )THEN
          ida_idx(il_next) = il_rand
          il_next = il_next + 1
        END IF
      END DO
      CALL select_sort( ida_idx )
    END IF
    IF( PRESENT(stat) ) stat = il_stat

  END SUBROUTINE random_indexes

  !>random seed initialization
  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed

  !>generates a row vector y of n points linearly spaced between and including a and b. For n < 2, linspace returns b.
  !!
  !<
  SUBROUTINE linspace(rd_a, rd_b, id_n, rda_x)
    REAL(KIND=cp), INTENT(in) :: rd_a, rd_b
    INTEGER, INTENT(in) :: id_n
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_x
    !local variables
    INTEGER :: ibi

    DO ibi = 1, id_n
      rda_x(ibi) = ( (id_n - ibi+1)*rd_a + (ibi-1)*rd_b )/id_n
    END DO
  END SUBROUTINE linspace

  !> \brief evaluate polynomial function at equispaced points in the interval [0, 1]
  !! @param [in] rda_a coefficients of the polynomial funtion, these ares the coefs of monomes
  !! @param [in,out] rda_y values of the evaluated polynomial
  !! \details the polynomial is evaluated at equispaced points in the interval [0, 1]
  !<
  SUBROUTINE poly_eval(rda_a, rda_y)
    REAL(KIND=cp), DIMENSION(:), INTENT(in)    :: rda_a
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_y
    !local variables
    INTEGER :: ibi, ibj, il_n
    REAL(KIND=cp), DIMENSION( SIZE(rda_a) ) :: rla_tmp, rla_x

    il_n = SIZE(rda_a)
    CALL linspace(-1.0_cp, 1.0_cp, il_n, rla_x)
    DO ibi = 1, il_n
      rla_tmp(ibi) = rda_a(1) + rda_a(2)*rla_x(ibi)
    END DO

    DO ibj = 3, il_n
      DO ibi = 1, il_n
        rda_y(ibi) = rla_tmp(ibi) + rda_a(ibj)*rla_x(ibi)**(ibj-1)
      END DO
      rla_tmp = rda_y
    END DO
  END SUBROUTINE poly_eval

  !> \brief adjoint of polynomial evaluation
  !<
  SUBROUTINE poly_evalb(rda_a, rda_ab, rda_y, rda_yb)
    REAL(KIND=cp), DIMENSION(:), INTENT(in)    :: rda_a
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_ab
    REAL(KIND=cp), DIMENSION(:), INTENT(in)    :: rda_y
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_yb
    !local variables
    INTEGER :: ibi, ibj, il_n
    REAL(KIND=cp), DIMENSION( SIZE(rda_a) ) :: rla_x, rla_tmpb
    rla_tmpb = 0.0_cp
    il_n = SIZE(rda_a)
    CALL nothing(rda_y)
    CALL linspace(-1.0_cp, 1.0_cp, il_n, rla_x)

    DO ibj = il_n, 3, -1
      rda_yb = rda_yb + rla_tmpb
      rla_tmpb = 0.0_cp
      DO ibi = il_n, 1, -1
        rla_tmpb(ibi) = rla_tmpb(ibi) + rda_yb(ibi)
        rda_ab(ibi) = rda_ab(ibi) + rda_yb(ibj)*rla_x(ibi)**(ibj-1)
        rda_yb(ibi) = 0.0_cp
      END DO
      rla_tmpb = rda_yb
    END DO

    DO ibi = il_n, 1, -1
      rda_ab(1) = rda_ab(1) + rla_tmpb(ibi)
      rda_ab(2) = rda_ab(2) + rla_x(ibi)*rla_tmpb(ibi)
      rla_tmpb(ibi) = 0.0_cp
    END DO
  END SUBROUTINE poly_evalb

  !>\brief Initializes environment for the generation of random values for compressive sensing
  !<
  SUBROUTINE init_cs_tools()
    CALL init_random_seed()
  END SUBROUTINE init_cs_tools

  !>\brief Terminates the compressive sensing environment
  !!\details This routine does nothing. Only the associated routine init_cs_tools is usefull.
  !!For conformance, it should be used so that the use of finalize make logical sense
  !<
  SUBROUTINE finalize_cs_tools()
  END SUBROUTINE finalize_cs_tools

  !>\brief Initializes environment for polynomial evaluation
  !<
  SUBROUTINE init_poly_tools()
    CALL init_random_seed()
  END SUBROUTINE init_poly_tools

  !>\brief Terminates the polynomial environment
  !!\details This routine does nothing. Only the associated routine init_poly_tools is usefull.
  !!For conformance, it should be used so that the use of finalize make logical sense
  !<
  SUBROUTINE finalize_poly_tools()
  END SUBROUTINE finalize_poly_tools

END MODULE cs_tools