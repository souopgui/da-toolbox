!> @file regul_tools.f90
!! @brief  routines to compute the gradient of a discretized (finite differences) function
!! @author Innocent Souopgui
!! @version 1.0
!<
MODULE regul_tools
  USE general_constant
  USE debug_tools
IMPLICIT NONE

  !> @brief  Method used to compute the gradient
  INTEGER, PARAMETER ::&
    BACKWARD = 17001,&
    FORWARD  = 17002,&
    CENTERED = 17003
  !>default value for constant grid spacing step

  !> @brief  interface for gradient approximation
  INTERFACE gradient
    MODULE PROCEDURE gradient_cs_1d!constant spacing
    MODULE PROCEDURE gradient_vs_1d!variable spacing
    MODULE PROCEDURE gradient_cs_2d!constant spacing
    !MODULE PROCEDURE gradient_vs_2d!variable spacing
    !MODULE PROCEDURE gradient_cs_3d!constant spacing
    !MODULE PROCEDURE gradient_vs_3d!variable spacing
  END INTERFACE gradient

  !> @brief  interface for gradient approximation using centered differences
  INTERFACE gradient_centered
    MODULE PROCEDURE gradient_ccs_1d!constant spacing
    !MODULE PROCEDURE gradient_cvs_1d!variable spacing
    MODULE PROCEDURE gradient_ccs_2d!constant spacing
    !MODULE PROCEDURE gradient_cvs_2d!variable spacing
    !MODULE PROCEDURE gradient_ccs_3d!constant spacing
    !MODULE PROCEDURE gradient_cvs_3d!variable spacing
  END INTERFACE gradient_centered

  !> @brief  interface for gradient approximation using forward differences
  INTERFACE gradient_forward
    MODULE PROCEDURE gradient_fcs_1d!constant spacing
    MODULE PROCEDURE gradient_fvs_1d!variable spacing
    MODULE PROCEDURE gradient_fcs_2d!constant spacing
    !MODULE PROCEDURE gradient_fvs_2d!variable spacing
    !MODULE PROCEDURE gradient_fcs_3d!constant spacing
    !MODULE PROCEDURE gradient_fvs_3d!variable spacing
  END INTERFACE gradient_forward

  !> @brief  interface for gradient approximation using backward differences
  INTERFACE gradient_backward
    MODULE PROCEDURE gradient_bcs_1d!constant spacing
    MODULE PROCEDURE gradient_bvs_1d!variable spacing
    MODULE PROCEDURE gradient_bcs_2d!constant spacing
    !MODULE PROCEDURE gradient_bvs_2d!variable spacing
    !MODULE PROCEDURE gradient_bcs_3d!constant spacing
    !MODULE PROCEDURE gradient_bvs_3d!variable spacing
  END INTERFACE gradient_backward

  !> @brief  interface for geostrophic_wind
  INTERFACE geostrophic_wind
    MODULE PROCEDURE geostrophic_wind_cgrid !>"staggered" Arakawa C-grid
  END INTERFACE geostrophic_wind

  !> @brief  interface for gradient regularization
  !! @details  Due to nonlinearities and the function implementation, the choice has been made to not accept default parameters; This made it easy to write adjoint code
  !<
  INTERFACE grad_regul
    MODULE PROCEDURE grad_regul_cs_1d!constant spacing
    MODULE PROCEDURE grad_regul_cs_1d_default!constant spacing with default value 1.0
    MODULE PROCEDURE grad_regul_vs_1d!variable spacing
    MODULE PROCEDURE grad_regul_cs_2d!constant spacing
    MODULE PROCEDURE grad_regul_cs_v2d!2d arrays are vectorized
    MODULE PROCEDURE grad_regul_cs_2d_default!constant spacing with default value 1.0
    MODULE PROCEDURE grad_regul_cs_v2d_default!2d arrays are vectorized
    !MODULE PROCEDURE grad_regul_vs_2d!variable spacing
    !MODULE PROCEDURE grad_regul_cs_3d!constant spacing
    !MODULE PROCEDURE grad_regul_cs_3d_default!constant spacing with default value 1.0
    !MODULE PROCEDURE grad_regul_vs_3d!variable spacing
  END INTERFACE grad_regul

  INTERFACE tikhonov_regul
    MODULE PROCEDURE tikhonov_regul_1d
    MODULE PROCEDURE tikhonov_regul_2d
  END INTERFACE tikhonov_regul

  INTERFACE geostrophic_regul
    !MODULE PROCEDURE tikhonov_regul_1d
    MODULE PROCEDURE geostrophic_regul_cgrid
  END INTERFACE geostrophic_regul

  !> @brief  interface for square L2 norm
  INTERFACE square_L2_norm
    MODULE PROCEDURE square_L2_norm_1d
    MODULE PROCEDURE square_L2_norm_2d
    !MODULE PROCEDURE square_L2_norm_3d
  END INTERFACE square_L2_norm

CONTAINS

  !>simple function used to set the discretization step
  FUNCTION disc_step(rd_dx) RESULT(rl_dx)
    REAL(KIND=cp), OPTIONAL, INTENT(IN) :: rd_dx !! dxspace step
    !local variables
    REAL(KIND=cp) :: rl_dx
    IF( PRESENT(rd_dx) )THEN
      rl_dx = rd_dx
    ELSE
      rl_dx = 1.0_cp
    END IF
  END FUNCTION disc_step

  !> @brief  Computes the square of the L2 norm of a discretized scalar function of one variable
  !! @param [in] rda_f function whose the L2 norm is requested
  !<
  FUNCTION square_L2_norm_1d(rda_f) RESULT(rl_l2)
    REAL(KIND=cp), DIMENSION(:),   INTENT(IN)  :: rda_f!!scalar field
    REAL(KIND=cp) :: rl_l2

    rl_l2 = 0.5_cp*SUM(rda_f**2)
  END FUNCTION square_L2_norm_1d

  !> @brief  Computes the square of the L2 norm of a discretized scalar function of two variables
  !! @param [in] rda_f function whose the L2 norm is requested
  !<
  FUNCTION square_L2_norm_2d(rda_f) RESULT(rl_l2)
    REAL(KIND=cp), DIMENSION(:, :),   INTENT(IN)  :: rda_f!!scalar field
    REAL(KIND=cp) :: rl_l2

    rl_l2 = 0.5_cp*SUM(rda_f**2)
  END FUNCTION square_L2_norm_2d

  !> @brief  Computes the gradient of a scalar function of one variable, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rd_dx optional, discretization step
  !! @param [in] method method used for computation
  !! @details  If the discretization step is not given, it is assumed to be 1.0
  !<
  SUBROUTINE gradient_cs_1d(rda_f, rda_df, method, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df!!gradient
    INTEGER,                     INTENT(IN)  :: method
    REAL(KIND=cp), OPTIONAL,     INTENT(IN)  :: rd_dx !! dxspace step
    !local variables
    REAL(KIND=cp) :: rl_dx

    rl_dx = disc_step(rd_dx)
    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bcs_1d(rda_f, rda_df, rl_dx)
      CASE(FORWARD)
        CALL gradient_fcs_1d(rda_f, rda_df, rl_dx)
      CASE(CENTERED)
        CALL gradient_ccs_1d(rda_f, rda_df, rl_dx)
    END SELECT
  END SUBROUTINE gradient_cs_1d

  !> @brief  Computes the gradient of a scalar function of two variables, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_fx gradient along the first direction
  !! @param [out] rda_fy gradient along the second direction
  !! @param [in] rd_dx discretization step in the first direction
  !! @param [in] rd_dy discretization step in the second direction
  !! @param [in] method method used for computation
  !! @details  If the discretization steps are not given, they are assumed to be 1.0
  !<
  SUBROUTINE gradient_cs_2d(rda_f, rda_fx, rda_fy, method, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN)  :: rda_f!!scalar vector field which gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(OUT) :: rda_fx, rda_fy!!x (resp. y) component of the gradient
    INTEGER,  INTENT(IN)  :: method
    REAL(KIND=cp), OPTIONAL,        INTENT(IN)  :: rd_dx, rd_dy!! x(resp. y) space step
    !local variables
    REAL(KIND=cp)                               :: rl_dx, rl_dy!! x(resp. y) space step

    rl_dx = disc_step(rd_dx)
    rl_dy = disc_step(rd_dy)
    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bcs_2d(rda_f, rda_fx, rda_fy, rl_dx, rl_dy)
      CASE(FORWARD)
        CALL gradient_fcs_2d(rda_f, rda_fx, rda_fy, rl_dx, rl_dy)
      CASE(CENTERED)
        CALL gradient_ccs_2d(rda_f, rda_fx, rda_fy, rl_dx, rl_dy)
    END SELECT
  END SUBROUTINE gradient_cs_2d

  !> @brief  Computes the gradient of a scalar function of one variable, variable spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rda_x discretization step in the first direction
  !! @param [in] method method used for computation
  !<
  SUBROUTINE gradient_vs_1d(rda_f, rda_df, method, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f, rda_x
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df
    INTEGER,                     INTENT(IN)  :: method

    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bvs_1d(rda_f, rda_df, rda_x)
      CASE(FORWARD)
        CALL gradient_fvs_1d(rda_f, rda_df, rda_x)
      CASE(CENTERED)
        CALL stop_program('In gradient_vs_1d: the centered approximation is not yet implemented for variable spacing grid')
        !CALL gradient_cvs_1d(rda_f, rda_df, rda_x)
    END SELECT
  END SUBROUTINE gradient_vs_1d

  !> @brief  Computes the gradient of a scalar function of one variable, centered differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rd_dx discretization step
  !! @details  gradient of a scalar fields, center diferences on the interior and decentered difference on the boundaries
  !<
  SUBROUTINE gradient_ccs_1d(rda_f, rda_df, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df!!gradient
    REAL(KIND=cp),               INTENT(IN)  :: rd_dx !! dxspace step
    INTEGER :: n

    n = SIZE(rda_f)
    rda_df(2:n-1) = ( rda_f(3:n) - rda_f(1:n-2) )/(2*rd_dx)!interior point, centered diff
    rda_df(1    ) = ( rda_f(2  ) - rda_f(1    ) )/rd_dx!left boundary, forward difference
    rda_df(n    ) = ( rda_f(n  ) - rda_f(n - 1) )/rd_dx!right boundary, backward difference
  END SUBROUTINE gradient_ccs_1d

  !> @brief  Computes the gradient of a scalar function of two variables, centered differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_fx gradient along the first direction
  !! @param [out] rda_fy gradient along the second direction
  !! @param [in] rd_dx discretization step in the first direction
  !! @param [in] rd_dy discretization step in the second direction
  !! @details  gradient of a scalar fields, center diferences on the interior and decentered difference on the boundaries
  !<
  SUBROUTINE gradient_ccs_2d(rda_f, rda_fx, rda_fy, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN)  :: rda_f!!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(OUT) :: rda_fx, rda_fy!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)  :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_f, 1)
    ny = SIZE(rda_f, 2)
    !x component
    rda_fx(2:nx-1, :) = ( rda_f(3:nx, :) - rda_f(1:nx-2, :) )/(2*rd_dx)!x : interior point, centered diff
    rda_fx(1     , :) = ( rda_f(2   , :) - rda_f(1     , :) )/rd_dx!left boundary, forward difference
    rda_fx(nx    , :) = ( rda_f(nx  , :) - rda_f(nx - 1, :) )/rd_dx!right boundary, backward difference
    !y component
    rda_fy(:, 2:ny-1) = ( rda_f(:, 3:ny) - rda_f(:, 1:ny-2) )/(2*rd_dy)!y : interior point, centered diff
    rda_fy(:, 1     ) = ( rda_f(:, 2   ) - rda_f(:, 1     ) )/rd_dy!bottom boundary, forward difference
    rda_fy(:, ny    ) = ( rda_f(:, ny  ) - rda_f(:, ny - 1) )/rd_dy!upper boundary, backward difference
  END SUBROUTINE gradient_ccs_2d

  !> @brief  Computes the gradient of a scalar function of one variable, forward differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rd_dx discretization step in the first direction
  !<
  SUBROUTINE gradient_fcs_1d(rda_f, rda_df, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f !!scalar function which gradient is required
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df!!gradient
    REAL(KIND=cp),               INTENT(IN)  :: rd_dx !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_f)
    rda_df(1:n-1) = ( rda_f(2:n) - rda_f(1:n-1) )/rd_dx!x : interior point, forward diff
    rda_df(n    ) = ( rda_f(n  ) - rda_f(n - 1) )/rd_dx!right boundary, backward difference
  END SUBROUTINE gradient_fcs_1d

  !> @brief  Computes the gradient of a scalar function of one variable, forward differences, variable spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rda_x discretization step in the first direction
  !<
  SUBROUTINE gradient_fvs_1d(rda_f, rda_df, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f, rda_x
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df
    INTEGER :: n

    n = SIZE(rda_f)
    rda_df(1:n-1) = ( rda_f(2:n) - rda_f(1:n-1) )/( rda_x(2:n)-rda_x(1:n-1) )!x : interior point, forward diff
    rda_df(n    ) = ( rda_f(n  ) - rda_f(n - 1) )/( rda_x(n  )-rda_x(n-1  ) )!right bound, back diff
  END SUBROUTINE gradient_fvs_1d

  !> @brief  Computes the gradient of a scalar function of two variables, forward differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_fx gradient along the first direction
  !! @param [out] rda_fy gradient along the second direction
  !! @param [in] rd_dx discretization step in the first direction
  !! @param [in] rd_dy discretization step in the second direction
  !<
  SUBROUTINE gradient_fcs_2d(rda_f, rda_fx, rda_fy, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :),   INTENT(IN)  :: rda_f!!scalar vector field which gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :),   INTENT(OUT) :: rda_fx, rda_fy!!x (resp. y) component of the gradient
    REAL(KIND=cp),                    INTENT(IN)  :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_f, 1)
    ny = SIZE(rda_f, 2)
    !x component
    rda_fx(1:nx-1, :) = ( rda_f(2:nx, :) - rda_f(1:nx-1, :) )/rd_dx!x : forward difference
    rda_fx(nx    , :) = ( rda_f(nx  , :) - rda_f(nx - 1, :) )/rd_dx!right boundary, backward difference
    !y component
    rda_fy(:, 1:ny-1) = ( rda_f(:, 2:ny) - rda_f(:, 1:ny-1) )/rd_dy!y : interior point, forward diff
    rda_fy(:, ny    ) = ( rda_f(:, ny  ) - rda_f(:, ny - 1) )/rd_dy!upper boundary, backward difference
  END SUBROUTINE gradient_fcs_2d

  !> @brief  Computes the gradient of a scalar function of one variable, backward differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rd_dx discretization step in the first direction
  !<
  SUBROUTINE gradient_bcs_1d(rda_f, rda_df, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f !!scalar function which gradient is required
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df!!gradient
    REAL(KIND=cp),               INTENT(IN)  :: rd_dx !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_f)
    rda_df(1  ) = ( rda_f(2  ) - rda_f(1    ) )/rd_dx!left boundary, forward difference
    rda_df(2:n) = ( rda_f(2:n) - rda_f(1:n-1) )/rd_dx!x : interior point, forward diff
  END SUBROUTINE gradient_bcs_1d

  !> @brief  Computes the gradient of a scalar function of one variable, backward differences, variable spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_df gradient
  !! @param [in] rda_x discretization points
  !<
  SUBROUTINE gradient_bvs_1d(rda_f, rda_df, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)  :: rda_f, rda_x
    REAL(KIND=cp), DIMENSION(:), INTENT(OUT) :: rda_df
    INTEGER :: n

    n = SIZE(rda_f)
    rda_df(1  ) = ( rda_f(2  ) - rda_f(1    ) )/( rda_x(2  )-rda_x(1    ) )!left bound, forward diff
    rda_df(2:n) = ( rda_f(2:n) - rda_f(1:n-1) )/( rda_x(2:n)-rda_x(1:n-1) )!x : interior point, forward diff
  END SUBROUTINE gradient_bvs_1d

  !> @brief  Computes the gradient of a scalar function of two variables, backward differences, constant spacing grid
  !! @param [in] rda_f function whose the gradient is required
  !! @param [out] rda_fx gradient along the first direction
  !! @param [out] rda_fy gradient along the second direction
  !! @param [in] rd_dx discretization step in the first direction
  !! @param [in] rd_dy discretization step in the second direction
  !<
  SUBROUTINE gradient_bcs_2d(rda_f, rda_fx, rda_fy, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :),   INTENT(IN)  :: rda_f!!scalar vector field which gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :),   INTENT(OUT) :: rda_fx, rda_fy!!x (resp. y) component of the gradient
    REAL(KIND=cp),                    INTENT(IN)  :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_f, 1)
    ny = SIZE(rda_f, 2)
    !x component
    rda_fx(1   , :) = ( rda_f(2   , :) - rda_f(1     , :) )/rd_dx!left boundary, forward difference
    rda_fx(2:nx, :) = ( rda_f(2:nx, :) - rda_f(1:nx-1, :) )/rd_dx!x : backward difference
    !y component
    rda_fy(:, 1   ) = ( rda_f(:, 2   ) - rda_f(:, 1     ) )/rd_dy!bottom boundary, forward difference
    rda_fy(:, 2:ny) = ( rda_f(:, 2:ny) - rda_f(:, 1:ny-1) )/rd_dy!y : interior point, backward diff
  END SUBROUTINE gradient_bcs_2d

  !geostrophic equilibrium tools

  !> @brief  Computes the theoretical wind that results from an exact balance between the Coriolis effect and the pressure gradient force
  !! @param [in] rda_h surface elevation or pressure
  !! @param [out] rda_u u component of the wind
  !! @param [out] rda_v v component of the wind
  !! @param [in] rd_dx space step in x direction
  !! @param [in] rd_dy space step in y direction
  !! @param [in] rd_g gravitational parameter
  !! @param [in] rda_fu Coriolis parameter at u location (fu varies only in the y direction)
  !! @param [in] rda_fv Coriolis parameter at v location (fv varies only in the y direction)
  !! @details
  !! The "staggered" Arakawa C-grid is assumed,
  !! it is supposeed that velocity components are not stored
  !! for faces that coincide with the boundary. So that
  !! the array h has 1 more element in the x direction (first index)
  !! than v, and 1 more element in the y direction (second index) than u.
  !! Decentered difference is used for the approximation of the gradient.
  !! Since each variable is stored at its own location,
  !! linear interpolation is used to approximate h at each point.
  !! call as CALL geostrophic_wind_cgrid(rda_h, rda_u, rda_v, td_gp%dx, td_gp%dy, td_sw%r_g, rga_fu, rga_fv)
  !<
  SUBROUTINE geostrophic_wind_cgrid(rda_h, rda_u, rda_v, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN)  :: rda_h
    REAL(KIND=cp), DIMENSION(:,:), INTENT(OUT) :: rda_u
    REAL(KIND=cp), DIMENSION(:,:), INTENT(OUT) :: rda_v
    REAL(KIND=cp), DIMENSION(:)  , INTENT(IN)  :: rda_fu, rda_fv
    REAL(KIND=cp), INTENT(IN) :: rd_dx, rd_dy, rd_g
    !Local variables
    INTEGER :: i0, im, j0, jm, ibi, ibj

    i0 = 1
    im = SIZE(rda_h, 1)
    j0 = 1
    jm = SIZE(rda_h, 2)
    DO ibj = 1, jm-1
      rda_u(:, ibj) = -rd_g/rda_fu(ibj) * (&
        ( rda_h(i0:im-1,ibj+1)+rda_h(i0+1:im,ibj+1) ) - (  rda_h(i0:im-1,ibj)+rda_h(i0+1:im,ibj) )&
      ) / (2.0_cp*rd_dy)
      !note, the factor 2 is used for the interpolation or average
    END DO
    !left differences for the to top row
    ibj = jm
    rda_u(:, ibj) = -rd_g/rda_fu(ibj) * (&
      ( rda_h(i0:im-1,ibj)+rda_h(i0+1:im,ibj) ) - (  rda_h(i0:im-1,ibj-1)+rda_h(i0+1:im,ibj-1) )&
    ) / (2.0_cp*rd_dy)

    DO ibi = 1, im-1
      rda_v(ibi,:) = rd_g/rda_fv * (&
        ( rda_h(ibi+1,j0:jm-1)+rda_h(ibi+1,j0+1:jm) )&
        - ( rda_h(ibi,j0:jm-1)+rda_h(ibi,j0+1:jm) )&
      ) / (2.0_cp*rd_dx)
    END DO
    !left differences for the to top row
    ibi = im
    rda_v(ibi,:) = rd_g/rda_fv * (&
      ( rda_h(ibi,j0:jm-1)+rda_h(ibi,j0+1:jm) ) - ( rda_h(ibi-1,j0:jm-1)+rda_h(ibi-1,j0+1:jm) )&
    ) / (2.0_cp*rd_dx)
  END SUBROUTINE geostrophic_wind_cgrid

  !> @brief  Computes the gradient regularization term associated with a scalar function of one variable
  !! @param [in] rda_f function whose the the gradient regularization term is required
  !! @param [in] method method used for computation
  !! @param [in] rd_dx discretization step
  !! @details  constant spacing grid
  !<
  FUNCTION grad_regul_cs_1d(rda_f, method, rd_dx) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    INTEGER,                     INTENT(IN) :: method
    REAL(KIND=cp),               INTENT(IN) :: rd_dx
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f) ) :: rla_fx!!gradient
    REAL(KIND=cp)                           :: rl_reg !! Regularization term

    CALL gradient(rda_f, rla_fx, method, rd_dx)
    rl_reg = square_L2_norm(rla_fx)
  END FUNCTION grad_regul_cs_1d

  !> @brief  Computes the gradient regularization term associated with a scalar function of one variable
  !! @param [in] rda_f function whose the the gradient regularization term is required
  !! @param [in] method method used for computation
  !! @details  constant spacing grid with default value 1.0
  !<
  FUNCTION grad_regul_cs_1d_default(rda_f, method) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    INTEGER,                     INTENT(IN) :: method
    !local variables
    !REAL(KIND=cp), DIMENSION( SIZE(rda_f) ) :: rla_fx!!gradient
    REAL(KIND=cp)                           :: rl_reg, rl_dx !! Regularization term

    rl_dx  = disc_step()
    rl_reg = grad_regul_cs_1d(rda_f, method, rl_dx)
  END FUNCTION grad_regul_cs_1d_default

  !> @brief  Computes the gradient regularization term associated with a scalar function of one variable
  !! @param [in] rda_f function whose the the gradient regularization term is required
  !! @param [in] method method used for computation
  !! @param [in] rda_x discretization points
  !! @details  variable spacing grid
  !<
  FUNCTION grad_regul_vs_1d(rda_f, method, rda_x) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f, rda_x
    INTEGER,                     INTENT(IN) :: method
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f) ) :: rla_fx!!gradient
    REAL(KIND=cp)                           :: rl_reg !! Regularization term

    CALL gradient(rda_f, rla_fx, method, rda_x)
    rl_reg = square_L2_norm(rla_fx)
  END FUNCTION grad_regul_vs_1d

  !> @brief  Computes the gradient regularization term associated with a scalar function of two variables
  !! @param [in] rda_f function whose the the gradient regularization term is required
  !! @param [in] method method used for computation
  !! @param [in] rd_dx discretization step in the first direction
  !! @param [in] rd_dy discretization step in the second direction
  !! @details  constant spacing grid
  !<
  FUNCTION grad_regul_cs_2d(rda_f, method, rd_dx, rd_dy) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN) :: rda_f
    INTEGER,                     INTENT(IN)    :: method
    REAL(KIND=cp),               INTENT(IN)    :: rd_dx, rd_dy
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f,1), SIZE(rda_f,2) )    :: rla_fx, rla_fy!!gradient
    REAL(KIND=cp)                              :: rl_reg, rl_regx, rl_regy !! Regularization term

    CALL gradient(rda_f, rla_fx, rla_fy, method, rd_dx, rd_dy)
    rl_regx = square_L2_norm(rla_fx)
    rl_regy = square_L2_norm(rla_fy)
    rl_reg  = rl_regx + rl_regy
  END FUNCTION grad_regul_cs_2d

  !> @brief  Computes the gradient regularization term associated with a scalar function of two variables
  !! \see grad_regul_cs_2d
  !! @details  the difference is that the arrays parameters are vectorized
  !<
  FUNCTION grad_regul_cs_v2d(rda_f, id_nx, id_ny, method, rd_dx, rd_dy) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    INTEGER,                     INTENT(IN) :: method, id_nx, id_ny
    REAL(KIND=cp),               INTENT(IN) :: rd_dx, rd_dy
    !local variables
    REAL(KIND=cp), DIMENSION(id_nx, id_ny)  :: rla_f
    REAL(KIND=cp)                           :: rl_reg !! Regularization term

    rla_f  = RESHAPE(rda_f, (/id_nx, id_ny/) )
    rl_reg = grad_regul_cs_2d(rla_f, method, rd_dx, rd_dy)
  END FUNCTION grad_regul_cs_v2d

  !> @brief  Computes the gradient regularization term associated with a scalar function of two variables
  !! @param [in] rda_f function whose the the gradient regularization term is required
  !! @param [in] method method used for computation
  !! @details  constant spacing grid with default value 1.0
  !<
  FUNCTION grad_regul_cs_2d_default(rda_f, method) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN) :: rda_f
    INTEGER,                        INTENT(IN) :: method
    !local variables
    REAL(KIND=cp)                           :: rl_dx, rl_dy
    REAL(KIND=cp)                           :: rl_reg !! Regularization term

    rl_dx  = disc_step()
    rl_dy  = disc_step()
    rl_reg = grad_regul_cs_2d(rda_f, method, rl_dx, rl_dy)
  END FUNCTION grad_regul_cs_2d_default

  !> @brief  Computes the gradient regularization term associated with a scalar function of two variables
  !! \see grad_regul_cs_2d_default
  !! @details  the difference is that the arrays parameters are vectorized
  !<
  FUNCTION grad_regul_cs_v2d_default(rda_f, id_nx, id_ny, method) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    INTEGER,                     INTENT(IN) :: method, id_nx, id_ny
    !local variables
    REAL(KIND=cp), DIMENSION(id_nx, id_ny)  :: rla_f
    REAL(KIND=cp)                  :: rl_dx, rl_dy
    REAL(KIND=cp)                  :: rl_reg !! Regularization term

    rl_dx  = disc_step()
    rl_dy  = disc_step()
    rla_f  = RESHAPE(rda_f, (/id_nx, id_ny/) )
    rl_reg = grad_regul_cs_2d(rla_f, method, rl_dx, rl_dy)
  END FUNCTION grad_regul_cs_v2d_default

  !> @brief  Compute the Tikhonov regularization term
  !! @param [in] rda_f variable on which the regularization is computed
  !! @param [in] rda_bf background estimation of rda_f
  !<
  FUNCTION tikhonov_regul_1d(rda_f, rda_bf) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: rda_bf
    !local variables
    REAL(KIND=cp)                  :: rl_reg !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_f,1)) :: rla_f

    IF( PRESENT(rda_bf) )THEN
      rla_f = rda_f - rda_bf
    ELSE
      rla_f = rda_f
    END IF
    rl_reg = square_L2_norm(rla_f)
  END FUNCTION tikhonov_regul_1d

  !> @brief  Compute the Tikhonov regularization term
  !! @param [in] rda_f variable on which the regularization is computed
  !! @param [in] rda_bf background estimation of rda_f
  !<
  FUNCTION tikhonov_regul_2d(rda_f, rda_bf) RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN) :: rda_f
    REAL(KIND=cp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: rda_bf
    !local variables
    REAL(KIND=cp)                  :: rl_reg !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_f,1),SIZE(rda_f,2)) :: rla_f

    IF( PRESENT(rda_bf) )THEN
      rla_f = rda_f - rda_bf
    ELSE
      rla_f = rda_f
    END IF
    rl_reg = square_L2_norm(rla_f)
  END FUNCTION tikhonov_regul_2d

  !> @brief  Computes the regulaisation term that measures the distance between the estimated wind and the geostrophic wind
  !! @param [in] rda_h surface elevation or pressure
  !! @param [out] rda_u u component of the wind
  !! @param [out] rda_v v component of the wind
  !! @param [in] rd_dx space step in x direction
  !! @param [in] rd_dy space step in y direction
  !! @param [in] rd_g gravitational parameter
  !! @param [in] rda_fu Coriolis parameter at u location (fu varies only in the y direction)
  !! @param [in] rda_fv Coriolis parameter at v location (fv varies only in the y direction)
  !! The "staggered" Arakawa C-grid is assumed
  !! call as geostrophic_regul(rda_h, rda_u, rda_v, td_gp%dx, td_gp%dy, td_sw%r_g, rga_fu, rga_fv)
  !<
  FUNCTION geostrophic_regul_cgrid(rda_h, rda_u, rda_v, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)RESULT(rl_reg)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN) :: rda_h, rda_u, rda_v
    REAL(KIND=cp), DIMENSION(:)  , INTENT(IN)  :: rda_fu, rda_fv
    REAL(KIND=cp), INTENT(IN) :: rd_dx, rd_dy, rd_g
    !local variables
    REAL(KIND=cp) :: rl_reg, rl_regu, rl_regv !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_u,1),SIZE(rda_u,2)) :: rla_gu, rla_du
    REAL(KIND=cp), DIMENSION(SIZE(rda_v,1),SIZE(rda_v,2)) :: rla_gv, rla_dv

    CALL geostrophic_wind_cgrid(rda_h, rla_gu, rla_gv, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)
    rla_du = rda_u - rla_gu
    rl_regu = square_L2_norm(rla_du)
    rla_dv = rda_v - rla_gv
    rl_regv = square_L2_norm(rla_dv)
    rl_reg = rl_regu + rl_regv
  END FUNCTION geostrophic_regul_cgrid

END MODULE regul_tools