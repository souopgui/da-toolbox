!> \file gradientadj.f90
!! \brief routines to compute the gradient of a discretized (finite differences) function
!! @author Innocent Souopgui
!! @version 1.0
!<
MODULE regul_toolsadj
  USE general_constant
  USE regul_tools
  USE debug_tools
IMPLICIT NONE

  !> \brief interface for gradient approximation
  INTERFACE gradientAdj
    MODULE PROCEDURE gradient_cs_1dAdj!constant spacing
    MODULE PROCEDURE gradient_vs_1dAdj!variable spacing
    MODULE PROCEDURE gradient_cs_2dAdj!constant spacing
    !MODULE PROCEDURE gradient_vs_2dAdj!variable spacing
    !MODULE PROCEDURE gradient_cs_3dAdj!constant spacing
    !MODULE PROCEDURE gradient_vs_3dAdj!variable spacing
  END INTERFACE gradientAdj

  !> \brief interface for gradient approximation using centered differences
  INTERFACE gradient_centeredAdj
    MODULE PROCEDURE gradient_ccs_1dAdj
    !MODULE PROCEDURE gradient_cvs_1dAdj
    MODULE PROCEDURE gradient_ccs_2dAdj
    !MODULE PROCEDURE gradient_cvs_2dAdj
    !MODULE PROCEDURE gradient_ccs_3dAdj
    !MODULE PROCEDURE gradient_cvs_3dAdj
  END INTERFACE gradient_centeredAdj

  !> \brief interface for gradient approximation using forward differences
  INTERFACE gradient_forwardAdj
    MODULE PROCEDURE gradient_fcs_1dAdj
    MODULE PROCEDURE gradient_fvs_1dAdj
    MODULE PROCEDURE gradient_fcs_2dAdj
    !MODULE PROCEDURE gradient_fvs_2dAdj
    !MODULE PROCEDURE gradient_fcs_3dAdj
    !MODULE PROCEDURE gradient_fvs_3dAdj
  END INTERFACE gradient_forwardAdj

  !> \brief interface for gradient approximation using backward differences
  INTERFACE gradient_backwardAdj
    MODULE PROCEDURE gradient_bcs_1dAdj
    MODULE PROCEDURE gradient_bvs_1dAdj
    MODULE PROCEDURE gradient_bcs_2dAdj
    !MODULE PROCEDURE gradient_bvs_2dAdj
    !MODULE PROCEDURE gradient_bcs_3dAdj
    !MODULE PROCEDURE gradient_bvs_3dAdj
  END INTERFACE gradient_backwardAdj

  !> \brief interface for the adjoint of geostrophic_wind
  INTERFACE geostrophic_windAdj
    MODULE PROCEDURE geostrophic_wind_cgridAdj !>"staggered" Arakawa C-grid
  END INTERFACE geostrophic_windAdj

  !> \brief interface for gradient regularization
  INTERFACE grad_regulAdj
    MODULE PROCEDURE grad_regul_cs_1dAdj!constant spacing
    MODULE PROCEDURE grad_regul_cs_1d_defaultAdj!constant spacing
    MODULE PROCEDURE grad_regul_vs_1dAdj!variable spacing
    MODULE PROCEDURE grad_regul_cs_2dAdj!constant spacing
    MODULE PROCEDURE grad_regul_cs_v2dAdj
    MODULE PROCEDURE grad_regul_cs_2d_defaultAdj!constant spacing
    MODULE PROCEDURE grad_regul_cs_v2d_defaultAdj
    !MODULE PROCEDURE grad_regul_vs_2dAdj!variable spacing
    !MODULE PROCEDURE grad_regul_cs_3dAdj!constant spacing
    !MODULE PROCEDURE grad_regul_cs_3d_defaultAdj!constant spacing
    !MODULE PROCEDURE grad_regul_vs_3dAdj!variable spacing
  END INTERFACE grad_regulAdj

  INTERFACE tikhonov_regulAdj
    MODULE PROCEDURE tikhonov_regul_1dAdj
    MODULE PROCEDURE tikhonov_regul_2dAdj
  END INTERFACE tikhonov_regulAdj

  INTERFACE geostrophic_regulAdj
    !MODULE PROCEDURE tikhonov_regul_1d
    MODULE PROCEDURE geostrophic_regul_cgridAdj
  END INTERFACE geostrophic_regulAdj

  !> \brief interface for square L2 norm
  INTERFACE square_L2_normAdj
    MODULE PROCEDURE square_L2_norm_1dAdj
    MODULE PROCEDURE square_L2_norm_2dAdj
    !MODULE PROCEDURE square_L2_norm_3dAdj
  END INTERFACE square_L2_normAdj

CONTAINS

  SUBROUTINE square_L2_norm_1dAdj(rda_f, rda_fad, rl_l2ad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_f!!scalar field
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad!!adjoint scalar field
    REAL(KIND=cp)              , INTENT(INOUT) :: rl_l2ad

    rda_fad = rda_fad + rl_l2ad * rda_f
    rl_l2ad = 0.0_dp
  END SUBROUTINE square_L2_norm_1dAdj

  SUBROUTINE square_L2_norm_2dAdj(rda_f, rda_fad, rl_l2ad)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(IN)    :: rda_f!!scalar field
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fad!!adjoint scalar field
    REAL(KIND=cp)                 , INTENT(INOUT) :: rl_l2ad

    rda_fad = rda_fad + rl_l2ad * rda_f
    rl_l2ad = 0.0_dp
  END SUBROUTINE square_L2_norm_2dAdj

  !> \brief Computes the gradient of a scalar function of one variable, constant spacing grid
  !<
  SUBROUTINE gradient_cs_1dAdj(rda_fad, rda_dfad, method, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad
    INTEGER,                     INTENT(IN)    :: method
    REAL(KIND=cp), OPTIONAL,     INTENT(IN)    :: rd_dx !! dxspace step
    !local variables
    REAL(KIND=cp) :: rl_dx

    rl_dx = disc_step(rd_dx)
    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bcs_1dAdj(rda_fad, rda_dfad, rl_dx)
      CASE(FORWARD)
        CALL gradient_fcs_1dAdj(rda_fad, rda_dfad, rl_dx)
      CASE(CENTERED)
        CALL gradient_ccs_1dAdj(rda_fad, rda_dfad, rl_dx)
    END SELECT
  END SUBROUTINE gradient_cs_1dAdj

  !> \brief Computes the gradient of a scalar function of two variables, constant spacing grid
  !<
  SUBROUTINE gradient_cs_2dAdj(rda_fad, rda_fxad, rda_fyad, method, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fad!!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fxad, rda_fyad!!x (resp. y) component of the gradient
    INTEGER,                        INTENT(IN)    :: method
    REAL(KIND=cp), OPTIONAL,        INTENT(IN)    :: rd_dx, rd_dy!! x(resp. y) space step
    !local variables
    REAL(KIND=cp)                                 :: rl_dx, rl_dy!! x(resp. y) space step

    rl_dx = disc_step(rd_dx)
    rl_dy = disc_step(rd_dy)
    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bcs_2dAdj(rda_fad, rda_fxad, rda_fyad, rl_dx, rl_dy)
      CASE(FORWARD)
        CALL gradient_fcs_2dAdj(rda_fad, rda_fxad, rda_fyad, rl_dx, rl_dy)
      CASE(CENTERED)
        CALL gradient_ccs_2dAdj(rda_fad, rda_fxad, rda_fyad, rl_dx, rl_dy)
    END SELECT
  END SUBROUTINE gradient_cs_2dAdj

  !> \brief Computes the gradient of a scalar function of one variable, variable spacing grid
  !! \param[in] rda_f function whose the gradient is required
  !! \param[out] rda_df gradient
  !! \param[in] rda_x discretization step in the first direction
  !! \param[in] method method used for computation
  !<
  SUBROUTINE gradient_vs_1dAdj(rda_fad, rda_dfad, method, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_x
    INTEGER,                     INTENT(IN)  :: method

    SELECT CASE (method)
      CASE(BACKWARD)
        CALL gradient_bvs_1dAdj(rda_fad, rda_dfad, rda_x)
      CASE(FORWARD)
        CALL gradient_fvs_1dAdj(rda_fad, rda_dfad, rda_x)
      CASE(CENTERED)
        CALL stop_program('In gradient_vs_1d: the centered approximation is not yet implemented for variable spacing grid')
        !CALL gradient_cvs_1dAdj(rda_fad, rda_dfad, rda_x)
    END SELECT
  END SUBROUTINE gradient_vs_1dAdj

  SUBROUTINE gradient_fcs_1dAdj(rda_fad, rda_dfad, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)    :: rd_dx   !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_fad)

    !right boundary
    rda_fad(n  ) = rda_fad(n  ) + rda_dfad(n)/rd_dx
    rda_fad(n-1) = rda_fad(n-1) - rda_dfad(n)/rd_dx
    rda_dfad(n ) = 0.0_dp
    !
    rda_fad (2:n  ) = rda_fad(2:n  ) + rda_dfad(1:n-1)/rd_dx
    rda_fad (1:n-1) = rda_fad(1:n-1) - rda_dfad(1:n-1)/rd_dx
    rda_dfad(1:n-1) = 0.0_dp
  END SUBROUTINE gradient_fcs_1dAdj

  SUBROUTINE gradient_fvs_1dAdj(rda_fad, rda_dfad, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad!!x (resp. y) component of the gradient
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_x   !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_fad)

    !right boundary
    rda_fad (n  ) = rda_fad(n  ) + rda_dfad(n)/( rda_x(n  )-rda_x(n-1  ) )
    rda_fad (n-1) = rda_fad(n-1) - rda_dfad(n)/( rda_x(n  )-rda_x(n-1  ) )
    rda_dfad(n  ) = 0.0_dp
    !
    rda_fad (2:n  ) = rda_fad(2:n  ) + rda_dfad(1:n-1)/( rda_x(2:n)-rda_x(1:n-1) )
    rda_fad (1:n-1) = rda_fad(1:n-1) - rda_dfad(1:n-1)/( rda_x(2:n)-rda_x(1:n-1) )
    rda_dfad(1:n-1) = 0.0_dp
  END SUBROUTINE gradient_fvs_1dAdj

  SUBROUTINE gradient_fcs_2dAdj(rda_fad, rda_fxad, rda_fyad, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fad!!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fxad, rda_fyad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)    :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_fad, 1)
    ny = SIZE(rda_fad, 2)

    !y component
    !upper boundary
    rda_fad (:, ny  ) = rda_fad(:, ny  ) + rda_fyad(:, ny)/rd_dy
    rda_fad (:, ny-1) = rda_fad(:, ny-1) - rda_fyad(:, ny)/rd_dy
    rda_fyad(:, ny  ) = 0.0_dp
    !
    rda_fad (:, 2:ny  ) = rda_fad(:, 2:ny  ) + rda_fyad(:, 1:ny-1)/rd_dy
    rda_fad (:, 1:ny-1) = rda_fad(:, 1:ny-1) - rda_fyad(:, 1:ny-1)/rd_dy
    rda_fyad(:, 1:ny-1) = 0.0_dp

    !x component
    !right boundary
    rda_fad(nx  , :) = rda_fad(nx  , :) + rda_fxad(nx, :)/rd_dx
    rda_fad(nx-1, :) = rda_fad(nx-1, :) - rda_fxad(nx, :)/rd_dx
    rda_fxad(nx , :) = 0.0_dp
    !
    rda_fad (2:nx  , :) = rda_fad(2:nx  , :) + rda_fxad(1:nx-1, :)/rd_dx
    rda_fad (1:nx-1, :) = rda_fad(1:nx-1, :) - rda_fxad(1:nx-1, :)/rd_dx
    rda_fxad(1:nx-1, :) = 0.0_dp
  END SUBROUTINE gradient_fcs_2dAdj

  SUBROUTINE gradient_ccs_1dAdj(rda_fad, rda_dfad, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)    :: rd_dx   !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_fad)

    !right boundary
    rda_fad (n  ) = rda_fad(n  ) + rda_dfad(n)/rd_dx
    rda_fad (n-1) = rda_fad(n-1) - rda_dfad(n)/rd_dx
    rda_dfad(n  ) = 0.0_dp
    !left boundary
    rda_fad (2) = rda_fad(2) + rda_dfad(1)/rd_dx
    rda_fad (1) = rda_fad(1) - rda_dfad(1)/rd_dx
    rda_dfad(1) = 0.0_dp
    !x : interior point
    rda_fad (3:n  ) = rda_fad(3:n  ) + rda_dfad(2:n-1)/(2*rd_dx)!centered diff
    rda_fad (1:n-2) = rda_fad(1:n-2) - rda_dfad(2:n-1)/(2*rd_dx)!centered diff
    rda_dfad(2:n-1) = 0.0_dp
  END SUBROUTINE gradient_ccs_1dAdj

  SUBROUTINE gradient_ccs_2dAdj(rda_fad, rda_fxad, rda_fyad, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fad!!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fxad, rda_fyad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)    :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_fad, 1)
    ny = SIZE(rda_fad, 2)

    !y component
    !upper boundary
    rda_fad (:, ny  ) = rda_fad(:, ny  ) + rda_fyad(:, ny)/rd_dy
    rda_fad (:, ny-1) = rda_fad(:, ny-1) - rda_fyad(:, ny)/rd_dy
    rda_fyad(:, ny  ) = 0.0_dp
    !bottom boundary
    rda_fad (:, 2) = rda_fad(:, 2) + rda_fyad(:, 1)/rd_dy
    rda_fad (:, 1) = rda_fad(:, 1) - rda_fyad(:, 1)/rd_dy
    rda_fyad(:, 1) = 0.0_dp
    !y : interior point
    rda_fad (:, 3:ny  ) = rda_fad(:, 3:ny  ) + rda_fyad(:, 2:ny-1)/(2*rd_dy)!centered diff
    rda_fad (:, 1:ny-2) = rda_fad(:, 1:ny-2) - rda_fyad(:, 2:ny-1)/(2*rd_dy)!centered diff
    rda_fyad(:, 2:ny-1) = 0.0_dp

    !x component
    !right boundary
    rda_fad(nx  , :) = rda_fad(nx  , :) + rda_fxad(nx, :)/rd_dx
    rda_fad(nx-1, :) = rda_fad(nx-1, :) - rda_fxad(nx, :)/rd_dx
    rda_fxad(nx , :) = 0.0_dp
    !left boundary
    rda_fad(2, :) = rda_fad(2, :) + rda_fxad(1, :)/rd_dx
    rda_fad(1, :) = rda_fad(1, :) - rda_fxad(1, :)/rd_dx
    rda_fxad(1, :)= 0.0_dp
    !x : interior point
    rda_fad (3:nx  , :) = rda_fad(3:nx  , :) + rda_fxad(2:nx-1, :)/(2*rd_dx)!centered diff
    rda_fad (1:nx-2, :) = rda_fad(1:nx-2, :) - rda_fxad(2:nx-1, :)/(2*rd_dx)!centered diff
    rda_fxad(2:nx-1, :) = 0.0_dp
  END SUBROUTINE gradient_ccs_2dAdj

  SUBROUTINE gradient_bcs_1dAdj(rda_fad, rda_dfad, rd_dx)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN) :: rd_dx   !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_fad)
    !
    rda_fad (2:n  ) = rda_fad(2:n  ) + rda_dfad(2:n)/rd_dx
    rda_fad (1:n-1) = rda_fad(1:n-1) - rda_dfad(2:n)/rd_dx
    rda_dfad(2:n  ) = 0.0_dp
    !left boundary
    rda_fad (2) = rda_fad(2) + rda_dfad(1)/rd_dx
    rda_fad (1) = rda_fad(1) - rda_dfad(1)/rd_dx
    rda_dfad(1) = 0.0_dp
  END SUBROUTINE gradient_bcs_1dAdj

  SUBROUTINE gradient_bvs_1dAdj(rda_fad, rda_dfad, rda_x)
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad !!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_dfad!!x (resp. y) component of the gradient
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_x   !! x(resp. y) space step
    INTEGER :: n

    n = SIZE(rda_fad)
    !
    rda_fad (2:n  ) = rda_fad(2:n  ) + rda_dfad(2:n)/( rda_x(2:n)-rda_x(1:n-1) )
    rda_fad (1:n-1) = rda_fad(1:n-1) - rda_dfad(2:n)/( rda_x(2:n)-rda_x(1:n-1) )
    rda_dfad(2:n  ) = 0.0_dp
    !left boundary
    rda_fad (2) = rda_fad(2) + rda_dfad(1)/( rda_x(2  )-rda_x(1    ) )
    rda_fad (1) = rda_fad(1) - rda_dfad(1)/( rda_x(2  )-rda_x(1    ) )
    rda_dfad(1) = 0.0_dp
  END SUBROUTINE gradient_bvs_1dAdj

  SUBROUTINE gradient_bcs_2dAdj(rda_fad, rda_fxad, rda_fyad, rd_dx, rd_dy)
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fad!!scalar vector fieldwhich gradient is to be computed
    REAL(KIND=cp), DIMENSION(:, :), INTENT(INOUT) :: rda_fxad, rda_fyad!!x (resp. y) component of the gradient
    REAL(KIND=cp),                  INTENT(IN)    :: rd_dx, rd_dy!! x(resp. y) space step
    INTEGER :: nx, ny

    nx = SIZE(rda_fad, 1)
    ny = SIZE(rda_fad, 2)

    !y component
    !
    rda_fad (:, 2:ny  ) = rda_fad(:, 2:ny  ) + rda_fyad(:, 2:ny)/rd_dy
    rda_fad (:, 1:ny-1) = rda_fad(:, 1:ny-1) - rda_fyad(:, 2:ny)/rd_dy
    rda_fyad(:, 2:ny  ) = 0.0_dp
    !bottom boundary
    rda_fad (:, 2) = rda_fad(:, 2) + rda_fyad(:, 1)/rd_dy
    rda_fad (:, 1) = rda_fad(:, 1) - rda_fyad(:, 1)/rd_dy
    rda_fyad(:, 1) = 0.0_dp

    !x component
    !
    rda_fad (2:nx  , :) = rda_fad(2:nx  , :) + rda_fxad(2:nx, :)/rd_dx
    rda_fad (1:nx-1, :) = rda_fad(1:nx-1, :) - rda_fxad(2:nx, :)/rd_dx
    rda_fxad(2:nx  , :) = 0.0_dp
    !left boundary
    rda_fad (2, :) = rda_fad(2, :) + rda_fxad(1, :)/rd_dx
    rda_fad (1, :) = rda_fad(1, :) - rda_fxad(1, :)/rd_dx
    rda_fxad(1 , :) = 0.0_dp
  END SUBROUTINE gradient_bcs_2dAdj

  !> \brief adjoint subroutine
  !!
  !<
  SUBROUTINE geostrophic_wind_cgridAdj(rda_hb, rda_ub, rda_vb, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)
    REAL(KIND=cp),DIMENSION(:,:), INTENT(IN OUT) :: rda_hb
    REAL(KIND=cp),DIMENSION(:,:), INTENT(IN OUT) :: rda_ub
    REAL(KIND=cp),DIMENSION(:,:), INTENT(IN OUT) :: rda_vb
    REAL(KIND=cp), DIMENSION(:)  , INTENT(IN)  :: rda_fu, rda_fv
    REAL(KIND=cp), INTENT(IN) :: rd_dx, rd_dy, rd_g
    !Local variables
    INTEGER :: i0, im, j0, jm, ibi, ibj

    i0 = 1
    im = SIZE(rda_hb, 1)
    j0 = 1
    jm = SIZE(rda_hb, 2)

    !left differences for the to top row
    ibi = im!recomputation of index
    rda_hb(ibi  ,j0:jm-1) = rda_hb(ibi  ,j0:jm-1) + (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
    rda_hb(ibi  ,j0+1:jm) = rda_hb(ibi  ,j0+1:jm) + (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
    rda_hb(ibi-1,j0:jm-1) = rda_hb(ibi-1,j0:jm-1) - (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
    rda_hb(ibi-1,j0+1:jm) = rda_hb(ibi-1,j0+1:jm) - (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
    rda_vb(ibi  ,:      ) = 0.0_cp

    !right differences
    DO ibi = im-1,1,-1
      rda_hb(ibi+1,j0:jm-1) = rda_hb(ibi+1,j0:jm-1) + (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
      rda_hb(ibi+1,j0+1:jm) = rda_hb(ibi+1,j0+1:jm) + (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
      rda_hb(ibi  ,j0:jm-1) = rda_hb(ibi  ,j0:jm-1) - (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
      rda_hb(ibi  ,j0+1:jm) = rda_hb(ibi  ,j0+1:jm) - (rd_g/rda_fv)*rda_vb(ibi,:)/(2.0_cp*rd_dx)
      rda_vb(ibi  ,:      ) = 0.0_cp
    END DO

    !left differences for the to top row
    ibj = jm
    rda_hb(i0:im-1,ibj  ) = rda_hb(i0:im-1,ibj  ) + ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
    rda_hb(i0+1:im,ibj  ) = rda_hb(i0+1:im,ibj  ) + ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
    rda_hb(i0:im-1,ibj-1) = rda_hb(i0:im-1,ibj-1) - ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
    rda_hb(i0+1:im,ibj-1) = rda_hb(i0+1:im,ibj-1) - ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
    rda_ub(:      ,ibj  ) = 0.0_cp

    !right differences
    DO ibj = jm-1,1,-1
      rda_hb(i0:im-1,ibj+1) = rda_hb(i0:im-1,ibj+1) + ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
      rda_hb(i0+1:im,ibj+1) = rda_hb(i0+1:im,ibj+1) + ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
      rda_hb(i0:im-1,ibj  ) = rda_hb(i0:im-1,ibj  ) - ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
      rda_hb(i0+1:im,ibj  ) = rda_hb(i0+1:im,ibj  ) - ( rd_g/rda_fu(ibj) )*rda_ub(:, ibj)/(2.0_cp*rd_dy)
      rda_ub(:      , ibj ) = 0.0_cp
    END DO
  END SUBROUTINE geostrophic_wind_cgridAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of one variable
  !<
  SUBROUTINE grad_regul_cs_1dAdj(rda_f, rda_fad, method, rd_dx, rl_regad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_f
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad
    INTEGER,                     INTENT(IN)    :: method
    REAL(KIND=cp),               INTENT(IN)    :: rd_dx
    REAL(KIND=cp),               INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f) )    :: rla_fx, rla_fxad!!gradient
    !REAL(KIND=cp)                              :: rl_reg !! Regularization term

    !zeroing local adjoint variables
    rla_fxad = 0.0_cp
    !end of zeroing
    !recomputing
    CALL gradient_cs_1d(rda_f, rla_fx, method, rd_dx)
    !end of recomputing
    CALL square_L2_normAdj(rla_fx, rla_fxad, rl_regad)
    CALL gradientAdj(rda_fad, rla_fxad, method, rd_dx)
  END SUBROUTINE grad_regul_cs_1dAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of one variable
  !<
  SUBROUTINE grad_regul_cs_1d_defaultAdj(rda_f, rda_fad, method, rl_regad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_f
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad
    INTEGER,                     INTENT(IN)    :: method
    REAL(KIND=cp),               INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp)                              :: rl_dx !! Regularization term

    !recomputing
    rl_dx = disc_step()
    !end of recomputing
    CALL grad_regul_cs_1dAdj(rda_f, rda_fad, method, rl_dx, rl_regad)
  END SUBROUTINE grad_regul_cs_1d_defaultAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of one variable
  !<
  SUBROUTINE grad_regul_vs_1dAdj(rda_f, rda_fad, method, rda_x, rl_regad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)    :: rda_f, rda_x
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT) :: rda_fad
    INTEGER,                     INTENT(IN)    :: method
    REAL(KIND=cp),               INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f) ) :: rla_fx!!gradient
    REAL(KIND=cp), DIMENSION( SIZE(rda_f) ) :: rla_fxad!!gradient
    !REAL(KIND=cp)                           :: rl_reg !! Regularization term
    !zeroing local adjoint variables
    rla_fxad = 0.0_cp
    !end of zeroing
    !recomputing
    CALL gradient(rda_f, rla_fx, method, rda_x)
    !end of recomputing
    CALL square_L2_normAdj(rla_fx, rla_fxad, rl_regad)
    CALL gradientAdj(rda_fad, rla_fxad, method, rda_x)
  END SUBROUTINE grad_regul_vs_1dAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of two variables
  !<
  SUBROUTINE grad_regul_cs_2dAdj(rda_f, rda_fad, method, rd_dx, rd_dy, rl_regad)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN)    :: rda_f
    REAL(KIND=cp), DIMENSION(:,:), INTENT(INOUT) :: rda_fad
    INTEGER,                       INTENT(IN)    :: method
    REAL(KIND=cp),                 INTENT(IN)    :: rd_dx, rd_dy
    REAL(KIND=cp),                 INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp), DIMENSION( SIZE(rda_f,1),SIZE(rda_f,2) ) :: rla_fx, rla_fy, rla_fxad, rla_fyad!!gradient
    REAL(KIND=cp) :: rl_regxad, rl_regyad!! Regularization term
    !zeroing local adjoint variables
    rla_fxad  = 0.0_cp
    rla_fyad  = 0.0_cp
    rl_regxad = 0.0_cp
    rl_regyad = 0.0_cp
    !end of zeroing
    !recomputing
    CALL gradient(rda_f, rla_fx, rla_fy, method, rd_dx, rd_dy)
    !end of recomputing
    rl_regyad = rl_regyad + rl_regad
    rl_regxad = rl_regxad + rl_regad
    rl_regad  = 0.0_cp
    CALL square_L2_normAdj(rla_fy, rla_fyad, rl_regyad)
    CALL square_L2_normAdj(rla_fx, rla_fxad, rl_regxad)
    CALL gradientAdj(rda_fad, rla_fxad, rla_fyad, method, rd_dx, rd_dy)
  END SUBROUTINE grad_regul_cs_2dAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of two variables
  !<
  SUBROUTINE grad_regul_cs_v2dAdj(rda_f, rda_fad, id_nx, id_ny, method, rd_dx, rd_dy, rl_regad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)      :: rda_f
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)   :: rda_fad
    INTEGER,                       INTENT(IN)    :: method, id_nx, id_ny
    REAL(KIND=cp),                 INTENT(IN)    :: rd_dx, rd_dy
    REAL(KIND=cp),                 INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp), DIMENSION(id_nx, id_ny)       :: rla_f, rla_fad
    !zeroing local adjoint variables
    rla_fad = 0.0_cp
    !end of zeroing
    !recomputing
    rla_f = RESHAPE(rda_f, (/id_nx, id_ny/) )
    !end of recomputing
    CALL grad_regul_cs_2dAdj(rla_f, rla_fad, method, rd_dx, rd_dy, rl_regad)
    rda_fad = rda_fad + RESHAPE(rla_fad, (/id_nx*id_ny/) )
    rla_fad = 0.0_cp
  END SUBROUTINE grad_regul_cs_v2dAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of two variables
  !<
  SUBROUTINE grad_regul_cs_2d_defaultAdj(rda_f, rda_fad, method, rl_regad)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN)    :: rda_f
    REAL(KIND=cp), DIMENSION(:,:), INTENT(INOUT) :: rda_fad
    INTEGER,                       INTENT(IN)    :: method
    REAL(KIND=cp),                 INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp)                                :: rl_dx, rl_dy
    !recomputing
    rl_dx = disc_step()
    rl_dy = disc_step()
    !end of recomputing
    CALL grad_regul_cs_2dAdj(rda_f, rda_fad, method, rl_dx, rl_dy, rl_regad)
  END SUBROUTINE grad_regul_cs_2d_defaultAdj

  !> \brief Computes the gradient regularization term associated with a scalar function of two variables
  !<
  SUBROUTINE grad_regul_cs_v2d_defaultAdj(rda_f, rda_fad, id_nx, id_ny, method, rl_regad)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN)      :: rda_f
    REAL(KIND=cp), DIMENSION(:), INTENT(INOUT)   :: rda_fad
    INTEGER,                       INTENT(IN)    :: method, id_nx, id_ny
    REAL(KIND=cp),                 INTENT(INOUT) :: rl_regad
    !local variables
    REAL(KIND=cp), DIMENSION(id_nx, id_ny)       :: rla_f, rla_fad
    REAL(KIND=cp)                                :: rl_dx, rl_dy
    !zeroing local adjoint variables
    rla_fad = 0.0_cp
    !end of zeroing
    !recomputing
    rl_dx = disc_step()
    rl_dy = disc_step()
    rla_f = RESHAPE(rda_f, (/id_nx, id_ny/) )
    !end of recomputing
    CALL grad_regul_cs_2dAdj(rla_f, rla_fad, method, rl_dx, rl_dy, rl_regad)
    rda_fad = rda_fad + RESHAPE(rla_fad, (/id_nx*id_ny/) )
    rla_fad = 0.0_cp
  END SUBROUTINE grad_regul_cs_v2d_defaultAdj

  !<
  SUBROUTINE tikhonov_regul_1dAdj(rda_f, rda_fad, rl_regad, rda_bf)
    REAL(KIND=cp), DIMENSION(:), INTENT(IN) :: rda_f
    REAL(KIND=cp), DIMENSION(:), INTENT(IN OUT) :: rda_fad
    REAL(KIND=cp), DIMENSION(:), OPTIONAL, INTENT(IN) :: rda_bf
    !local variables
    REAL(KIND=cp)                  :: rl_regad !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_f,1)) :: rla_f, rla_fad

    !recomputation
    IF( PRESENT(rda_bf) )THEN
      rla_f = rda_f - rda_bf
    ELSE
      rla_f = rda_f
    END IF
    rla_fad = 0.0_dp
    
    CALL square_L2_normAdj(rla_f, rla_fad, rl_regad)
    
    !this part is independent of the presence or not of the background term
    rda_fad = rda_fad + rla_fad
    rla_fad = 0.0_cp
  END SUBROUTINE tikhonov_regul_1dAdj

  !<
  SUBROUTINE tikhonov_regul_2dAdj(rda_f, rda_fad, rl_regad, rda_bf)
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN) :: rda_f
    REAL(KIND=cp), DIMENSION(:,:), OPTIONAL, INTENT(IN) :: rda_bf
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN OUT) :: rda_fad
    !local variables
    REAL(KIND=cp)                  :: rl_regad !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_f,1),SIZE(rda_f,2)) :: rla_f, rla_fad

    IF( PRESENT(rda_bf) )THEN
      rla_f = rda_f - rda_bf
    ELSE
      rla_f = rda_f
    END IF
    rla_fad = 0.0_dp
    CALL square_L2_normAdj(rla_f, rla_fad, rl_regad)
    !this part is independent of the presence or not of the background term
    rda_fad = rda_fad + rla_fad
    rla_fad = 0.0_cp
  END SUBROUTINE tikhonov_regul_2dAdj

  !> \brief adjoint subroutine
  SUBROUTINE geostrophic_regul_cgridAdj(&
                rda_h, rda_had, rda_u, rda_uad, rda_v, rda_vad, rd_dx, rd_dy, rd_g, rda_fu, rda_fv, rl_regad&
             )
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN)     :: rda_h, rda_u, rda_v
    REAL(KIND=cp), DIMENSION(:,:), INTENT(IN OUT) :: rda_had, rda_uad, rda_vad
    REAL(KIND=cp), DIMENSION(:)  , INTENT(IN)     :: rda_fu, rda_fv
    REAL(KIND=cp), INTENT(IN) :: rd_dx, rd_dy, rd_g
    REAL(KIND=cp), INTENT(IN OUT) :: rl_regad
    !local variables
    REAL(KIND=cp) :: rl_reguad, rl_regvad !! Regularization term
    REAL(KIND=cp), DIMENSION(SIZE(rda_u,1),SIZE(rda_u,2)) :: rla_gu, rla_du
    REAL(KIND=cp), DIMENSION(SIZE(rda_u,1),SIZE(rda_u,2)) :: rla_guad, rla_duad
    REAL(KIND=cp), DIMENSION(SIZE(rda_v,1),SIZE(rda_v,2)) :: rla_gv, rla_dv
    REAL(KIND=cp), DIMENSION(SIZE(rda_v,1),SIZE(rda_v,2)) :: rla_gvad, rla_dvad

    !zeroing local adjoint variables
    rl_reguad = 0.0_cp
    rl_regvad = 0.0_cp
    rla_guad  = 0.0_cp
    rla_duad  = 0.0_cp
    rla_gvad  = 0.0_cp
    rla_dvad  = 0.0_cp
    !recomputation
    CALL geostrophic_wind_cgrid(rda_h, rla_gu, rla_gv, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)
    rla_du = rda_u - rla_gu
    rla_dv = rda_v - rla_gv
    !end of recomputation
    rl_reguad = rl_reguad + rl_regad
    rl_regvad = rl_regvad + rl_regad
    rl_regad = 0.0_cp
    CALL square_L2_normAdj(rla_dv, rla_dvad, rl_regvad)
    rda_vad = rda_vad + rla_dvad
    rla_gvad = rla_gvad - rla_dvad
    rla_dvad = 0.0_cp
    CALL square_L2_normAdj(rla_du, rla_duad, rl_reguad)
    rda_uad = rda_uad + rla_duad
    rla_guad = rla_guad - rla_duad
    rla_duad = 0.0_cp
    CALL geostrophic_wind_cgridAdj(rda_had, rla_guad, rla_gvad, rd_dx, rd_dy, rd_g, rda_fu, rda_fv)
  END SUBROUTINE geostrophic_regul_cgridAdj

END MODULE regul_toolsadj