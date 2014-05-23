!> \file gd_tools.f90
!! \brief Generalised diffusion routines for BALAISE
!! @author Innocent Souopgui
!! @version 2.0
!! \details In this version, the mirror boundary condition can be replaced by and extrapolation which is the default.
!! The mirror boundary condition seems to give extra importance to internal points at the boundaries.
!! The removal of the mirror boundary condition do not ensure the theoritical result but, in practice, it produces better results
!! Instead of adding extra-cells to the grid with mirror values, the extra-cells are set to extrapolated values. In 1D, the value of the left extra-cell is given by:
!! f0 = f1 - (1/2f1' + 1/2f2'), this is the value at border cell translated by a linear combination
!! of the derivatives (decentered differences) at the two most external grid cells. Assuming the default 1-based indices in fortran, f0 is the value at the extra-cell
!! We used the default weight 1/2 for the linear combination. This value can be replaced if necessary
!<
MODULE gd_tools
  USE debug_tools
  USE checkpoint
IMPLICIT NONE
  !>
  REAL(KIND=cp), PARAMETER :: &
    epsilon_gd = 1.0d-8,&
    rm_w_bc = 0.5_cp
  INTEGER, PARAMETER :: &
    LINEAR_EXTRAPOLATION = 13001,&
    MIRROR_BOUNDARY      = 13002
  INTEGER, PARAMETER :: gd_boundary = LINEAR_EXTRAPOLATION
  INTEGER, PARAMETER :: &
    chkp_unit_1v_1d = 1,&
    chkp_unit_1v_2d = 2,&
    chkp_unit_2d_2v = 3

  !> \brief interface for generalised diffusion (quasi-)projection
  INTERFACE gd_projection
    MODULE PROCEDURE gd_projection_1v_1d!one variable on 1d
    MODULE PROCEDURE gd_projection_1v_2d!one variable on 2d
!     MODULE PROCEDURE gd_projection_2v_2d!two variables on 2d
  END INTERFACE gd_projection

  !> \brief interface for ensuring mirror boundaries
  INTERFACE gd_projection_bc
    MODULE PROCEDURE gd_projection_bc_1d!one dimension
    MODULE PROCEDURE gd_projection_bc_2d!two dimensions
!     MODULE PROCEDURE gd_projection_bc_3d!three dimensions
  END INTERFACE gd_projection_bc

CONTAINS

!SUBROUTINE gd_smooth(rda_v1, rda_v2, rda_u1, rda_u2, rda_phi1, rda_phi2, rd_dx, rd_dy, id_niter)

  !> \brief Generalised diffusion (quasi-)projection of one variable in one dimensional space
  !! \param[in] rda_v variable to be projected
  !! \param[out] rda_u result of the projection
  !! \param[in] rda_phi Trust function
  !! \param[in] rd_dx space discretization step
  !! \param[in] id_niter, number of iterations for the iterative process
  !! \param[in] c_unit id of the checkpointing unit, this must be an integer between 1 and im_nunit
  !! \details See checkpoint module for details on checkpointing units. The checkpointing unit is used to save intermediate results for adjoint calculation. If a configuration of the projection is used more than once and the adjoint calculation is necessary, the user should give different value to c_unit for each use of the configuration. If a configuration is used only once, the system can manage it. The configuration of the projection is defined by the dimensionality and the number of input to be projected. Example of configurations are : one variable in one dimension, two variables in one dimension, one variable in two dimensions, two variables in two dimension.
  !!
  !<
  SUBROUTINE gd_projection_1v_1d(rda_v, rda_u, rda_phi, rd_dx, id_niter, c_unit)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)   :: rda_v, rda_phi
    REAL(KIND=dp), DIMENSION( SIZE(rda_v) ), INTENT(OUT) :: rda_u!!result of the projection
    REAL(KIND=dp) , INTENT(IN) :: rd_dx!!space discretization step
    INTEGER, INTENT(IN) :: id_niter!!iteration count, amount of diffusion iterations
    INTEGER, OPTIONAL, INTENT(IN) :: c_unit!!id of the checkpointing unit
    !!local variables
    REAL(KIND=dp), DIMENSION( SIZE(rda_v) ) :: rla_tmp, rla_b, rla_c
    !!local extended variables used to handle mirror boundary conditions.  rla_ux, x for extended
    REAL(KIND=dp), DIMENSION( 0:SIZE(rda_v)+1 ) :: rla_ux

    REAL(KIND=dp)       :: rl_dt, rl_tmp, rl_a, rl_norm, rl_normr, rl_tmpnum, rl_tmpden, rl_quot
    INTEGER             :: ib_iter, il_ub, il_lb, chkp_unit

    IF( PRESENT(c_unit) )THEN
      chkp_unit = c_unit
    ELSE
      chkp_unit = chkp_unit_1v_1d
    END IF
    il_lb = 1
    il_ub = SIZE(rda_v)

    !To satisfy the condition of Courant-Friedrichs-Lewy, the time step is defined as follows
    rl_tmp= rd_dx**2/4.0_cp
    rl_dt = rl_tmp-MIN(epsilon_gd, rl_tmp/10.0_cp)!rd_dx**2 / ( 4.0_cp - epsilon_gd )

    rl_a  = rl_dt / ( rd_dx**2 )

    rla_b = -(rda_phi + 2.0_cp/(rd_dx**2) )*rl_dt + 1.0_cp

    rla_c = rl_dt * (rda_phi*rda_v)

    rla_ux(il_lb:il_ub) = rda_v
    CALL gd_projection_bc(rla_ux)
    !creating file for saving intermediate results used in adjoint calculations !chkp_unit
    !gd_init_checkpoint(id_chkpId, id_nVar, id_varSize, nrec)
    CALL chkp_open( chkp_unit, 1, il_ub )

    DO ib_iter = 1, id_niter
      !Saving intermediate result for adjoint calculation
      !chkp_save(unit, rda_u, id_rec)
      CALL chkp_save(chkp_unit, rla_ux(il_lb:il_ub), ib_iter)
      rla_tmp = &
        rl_a*( rla_ux(il_lb+1:il_ub+1) + rla_ux(il_lb-1:il_ub-1) )&
        + rla_ux(il_lb:il_ub)*rla_b + rla_c

      rla_ux(il_lb:il_ub) = rla_tmp
      CALL gd_projection_bc(rla_ux)
    END DO
    ! Saving intermediate result for adjoint calculation
    CALL chkp_save(chkp_unit, rla_ux(il_lb:il_ub), id_niter+1)
    CALL chkp_save(chkp_unit, rda_v, id_niter+2)
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    rl_norm = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb:il_ub)**2
    rl_normr= MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb:il_ub)
    IF ( rl_normr .GT. epsilon_gd ) THEN
      rl_tmpnum = SQRT(rl_norm)
      rl_tmpden = SQRT(rl_normr)
      rl_quot   = rl_tmpnum/rl_tmpden
      CALL debug(rl_quot, 'rl_mag/rl_magr = ')
      rda_u = rl_quot * rla_tmp
    ELSE
      rda_u = rla_tmp
    END IF
    CALL chkp_close(chkp_unit)
  END SUBROUTINE gd_projection_1v_1d

  !> \brief Generalised diffusion (quasi-)projection of one variable in two dimensional space
  !! \param[in] rda_v variable to be projected
  !! \param[out] rda_u result of the projection
  !! \param[in] rda_phi Trust function
  !! \param[in] rd_dx discretization step along the x direction
  !! \param[in] rd_dy discretization step along the y direction
  !! \param[in] id_niter, number of iterations for the iterative process
  !! \param[in] c_unit id of the checkpointing unit, this must be an integer between 1 and im_nunit
  !! \details See checkpoint module for details on checkpointing units. The checkpointing unit is used to save intermediate results for adjoint calculation. If a configuration of the projection is used more than once and the adjoint calculation is necessary, the user should give different value to c_unit for each use of the configuration. If a configuration is used only once, the system can manage it. The configuration of the projection is defined by the dimensionality and the number of input to be projected. Example of configurations are : one variable in one dimension, two variables in one dimension, one variable in two dimensions, two variables in two dimension.
  !!
  !<
  SUBROUTINE gd_projection_1v_2d(rda_v, rda_u, rda_phi, rd_dx, rd_dy, id_niter, c_unit)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: rda_v, rda_phi
    REAL(KIND=dp), DIMENSION( SIZE(rda_v,1), SIZE(rda_v,2) ), INTENT(OUT) :: rda_u!!result of the projection
    REAL(KIND=dp) , INTENT(IN) :: rd_dx, rd_dy!!space discretization step
    INTEGER, INTENT(IN) :: id_niter!!iteration count, amount of diffusion iterations
    INTEGER, OPTIONAL, INTENT(IN) :: c_unit!!id of the checkpointing unit
    !!local variables
    REAL(KIND=dp), DIMENSION( SIZE(rda_v,1), SIZE(rda_v,2) ) :: rla_tmp, rla_b, rla_c
    !!local extended variables used to handle mirror boundary conditions.  rla_ux, x for extended
    REAL(KIND=dp), DIMENSION( 0:SIZE(rda_v,1)+1, 0:SIZE(rda_v,2)+1 ) :: rla_ux

    REAL(KIND=dp)       :: rl_dt, rl_tmp, rl_a1, rl_a2, rl_norm, rl_normr, rl_tmpnum, rl_tmpden, rl_quot
    INTEGER             :: ib_iter, il_ub1, il_ub2, il_lb1, il_lb2, chkp_unit

    CALL debug('', 'In gd_projection_1v_2d       -500')
    IF( PRESENT(c_unit) )THEN
      chkp_unit = c_unit
    ELSE
      chkp_unit = chkp_unit_1v_2d
    END IF
    il_lb1 = 1
    il_lb2 = 1
    il_ub1 = SIZE(rda_v,1)
    il_ub2 = SIZE(rda_v,2)

    !To satisfy the condition of Courant-Friedrichs-Lewy, the time step is defined as follows
    rl_tmp = (rd_dx*rd_dy)**2/( 4.0_cp*(rd_dx**2 + rd_dy**2) )
    rl_dt = rl_tmp-MIN(epsilon_gd, rl_tmp/10.0_cp)

    rl_a1  = rl_dt / ( rd_dx**2 )
    rl_a2  = rl_dt / ( rd_dy**2 )

    rla_b = -(rda_phi + 2.0_cp/(rd_dx**2) + 2.0_cp/(rd_dy**2) )*rl_dt + 1.0_cp

    rla_c = rl_dt * (rda_phi*rda_v)

    rla_ux(il_lb1:il_ub1, il_lb2:il_ub2) = rda_v
    CALL debug('', 'In gd_projection_1v_2d       000')
    CALL gd_projection_bc(rla_ux)
    CALL debug('', 'In gd_projection_1v_2d       100')
    !creating file for saving intermediate results used in adjoint calculations !chkp_unit
    !gd_init_checkpoint(id_chkpId, id_nVar, id_varSize, nrec)
    CALL chkp_open( chkp_unit, 1, il_ub1*il_ub2 )
    CALL debug('', 'In gd_projection_1v_2d       200')
    DO ib_iter = 1, id_niter
      !Saving intermediate result for adjoint calculation
      !chkp_save(unit, rda_u, id_rec)
      CALL chkp_save(chkp_unit, rla_ux(il_lb1:il_ub1, il_lb2:il_ub2), ib_iter)
      rla_tmp = &
          rl_a1 * rla_ux(il_lb1+1:il_ub1+1, il_lb2:il_ub2)&
        + rl_a1 * rla_ux(il_lb1-1:il_ub1-1, il_lb2:il_ub2)&
        + rl_a2 * rla_ux(il_lb1:il_ub1, il_lb2+1:il_ub2+1)&
        + rl_a2 * rla_ux(il_lb1:il_ub1, il_lb2-1:il_ub2-1)&
        + rla_b * rla_ux(il_lb1:il_ub1, il_lb2:il_ub2)&
        + rla_c

      rla_ux(il_lb1:il_ub1, il_lb2:il_ub2) = rla_tmp
      CALL gd_projection_bc(rla_ux)
    END DO
    CALL debug('', 'In gd_projection_1v_2d       300')
    ! Saving intermediate result for adjoint calculation
    CALL chkp_save(chkp_unit, rla_ux(il_lb1:il_ub1, il_lb2:il_ub2), id_niter+1)
    CALL chkp_save(chkp_unit, rda_v, id_niter+2)
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    rl_norm = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb1:il_ub1, il_lb2:il_ub2)**2
    rl_normr= MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb1:il_ub1, il_lb2:il_ub2)
    IF ( rl_normr .GT. epsilon_gd ) THEN
      rl_tmpnum = SQRT(rl_norm)
      rl_tmpden = SQRT(rl_normr)
      rl_quot   = rl_tmpnum/rl_tmpden
      CALL debug(rl_quot, 'rl_mag/rl_magr = ')
      rda_u = rl_quot * rla_tmp
    ELSE
      rda_u = rla_tmp
    END IF
    CALL chkp_close(chkp_unit)
  END SUBROUTINE gd_projection_1v_2d

  !> \brief Ensure boundary conditions
  !! \param[out] rda_A array to ensured boundary conditions
  !<
  SUBROUTINE gd_projection_bc_1d(rda_A)
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)  :: rda_A
    INTEGER                    :: il_dm, il_fm

    il_dm = LBOUND(rda_A,1)
    il_fm = UBOUND(rda_A,1)
    SELECT CASE (gd_boundary)
      !the direction of the substraction operation is adjusted to have uniform sign at every boundaries
      CASE(LINEAR_EXTRAPOLATION)
        rda_A(il_dm) = rda_A(il_dm+1) &
          + rm_w_bc         *( rda_A(il_dm+1)-rda_A(il_dm+2) )&
          + (1.0_cp-rm_w_bc)*( rda_A(il_dm+2)-rda_A(il_dm+3) )
        rda_A(il_fm) = rda_A(il_fm-1) &
          + rm_w_bc         *( rda_A(il_fm-1)-rda_A(il_fm-2) )&
          + (1.0_cp-rm_w_bc)*( rda_A(il_fm-2)-rda_A(il_fm-3) )
      CASE(MIRROR_BOUNDARY)
        rda_A(il_dm) = rda_A(il_dm+2)
        rda_A(il_fm) = rda_A(il_fm-2)
    END SELECT
  END SUBROUTINE gd_projection_bc_1d

  SUBROUTINE gd_projection_bc_2d(rda_A)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)  :: rda_A
    INTEGER                              :: il_dm, il_fm, il_dn, il_fn
    !Make sure that the matrice A of size >3x3 has mirror boundary condition
    il_dm = LBOUND(rda_A,1)
    il_fm = UBOUND(rda_A,1)
    il_dn = LBOUND(rda_A,2)
    il_fn = UBOUND(rda_A,2)
    SELECT CASE (gd_boundary)
      !the direction of the substraction operation is adjusted to have uniform sign at every boundaries
      CASE(LINEAR_EXTRAPOLATION)
        rda_A(il_dm, :) = rda_A(il_dm+1, :)&
          + rm_w_bc         *( rda_A(il_dm+1, :)-rda_A(il_dm+2, :) )&
          + (1.0_cp-rm_w_bc)*( rda_A(il_dm+2, :)-rda_A(il_dm+3, :) )
        rda_A(il_fm, :) = rda_A(il_fm-1, :)&
          + rm_w_bc         *( rda_A(il_fm-1, :)-rda_A(il_fm-2, :) )&
          + (1.0_cp-rm_w_bc)*( rda_A(il_fm-2, :)-rda_A(il_fm-3, :) )
        rda_A(:, il_dn) = rda_A(:, il_dn+1)&
          + rm_w_bc         *( rda_A(:, il_dn+1)-rda_A(:, il_dn+2) )&
          + (1.0_cp-rm_w_bc)*( rda_A(:, il_dn+2)-rda_A(:, il_dn+3) )
        rda_A(:, il_fn) = rda_A(:, il_fn-1)&
          + rm_w_bc         *( rda_A(:, il_fn-1)-rda_A(:, il_fn-2) )&
          + (1.0_cp-rm_w_bc)*( rda_A(:, il_fn-2)-rda_A(:, il_fn-3) )
      CASE(MIRROR_BOUNDARY)
        rda_A(il_dm, :) = rda_A(il_dm+2, :)
        rda_A(il_fm, :) = rda_A(il_fm-2, :)
        rda_A(:, il_dn) = rda_A(:, il_dn+2)
        rda_A(:, il_fn) = rda_A(:, il_fn-2)
    END SELECT
  END SUBROUTINE gd_projection_bc_2d

END MODULE gd_tools