!> \file gd_toolsadj.f90
!! @brief Generalised diffusion routines for BALAISE
!! @author Innocent Souopgui
!! @version 2.0
!! @details In this version, the mirror boundary condition is replaced by and extrapolation.
!! The mirror boundary condition seems to give extra importance to internal points at the boundaries.
!! The removal of the mirror boundary condition do not ensure the theoritical result but, in practice, it produces better results
!! Instead of adding extra-cells to the grid with mirror values, the extra-cells are set to extrapolated values. In 1D, the value of the left extra-cell is given by:
!! f0 = f1 - (1/2f1' + 1/2f2'), this is the value at border cell translated by a linear combination
!! of the derivatives (decentered differences) at the two most external grid cells. Assuming the default 1-based indices in fortran, f0 is the value at the extra-cell
!! We used the default weight 1/2 for the linear combination. This value can be replaced if necessary
!<
MODULE gd_toolsadj
  USE general_tools
  USE gd_tools
  USE debug_tools
  USE checkpoint
  !USE adjoint_tools
IMPLICIT NONE

  !> @brief interface for generalised diffusion (quasi-)projection
  INTERFACE gd_projectionadj
    MODULE PROCEDURE gd_projection_1v_1dadj!one variable on 1d
    MODULE PROCEDURE gd_projection_1v_2dadj!one variable on 2d
!     MODULE PROCEDURE gd_projection_2v_2dadj!two variables on 2d
  END INTERFACE gd_projectionadj

  !> @brief interface for ensuring mirror boundaries
  INTERFACE gd_projection_bcadj
    MODULE PROCEDURE gd_projection_bc_1dadj!one dimension
    MODULE PROCEDURE gd_projection_bc_2dadj!two dimensions
!     MODULE PROCEDURE gd_projection_bc_3dadj!three dimensions
  END INTERFACE gd_projection_bcadj

CONTAINS

  !> @brief Generalised diffusion (quasi-)projection of one variable in one dimensional space
  !! @param [in] rda_v variable to be projected
  !! @param [in,out] rda_vad adjoint variable
  !! @param [in,out] rda_uad adjoint variable
  !! @param [in] rda_phi Trust function
  !! @param [in] rd_dx space discretization step
  !! @param [in] id_niter number of iterations for the iterative process
  !! @param [in] c_unit id of the checkpointing unit, this must be
  !!  an integer between 1 and im_nunit
  !! @details
  !! @todo check the extra parameter rda_ua. it seems to be there by error.
  !<
  SUBROUTINE gd_projection_1v_1dadj(rda_v, rda_vad, rda_uad, rda_phi, rd_dx, id_niter, c_unit)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)    :: rda_v, rda_phi
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT) :: rda_vad, rda_uad!
    REAL(KIND=dp) , INTENT(IN) :: rd_dx!!space discretization step
    INTEGER, INTENT(IN) :: id_niter!!iteration count, amount of diffusion iterations
    INTEGER, OPTIONAL, INTENT(IN) :: c_unit!!id of the checkpointing unit
    !!local variables
    REAL(KIND=dp), DIMENSION( SIZE(rda_v) ) :: rla_tmp, rla_b, rla_c, rla_tmpad, rla_bad, rla_cad
    !!local extended variables used to handle mirror boundary conditions.  rla_ux, x for extended
    REAL(KIND=dp), DIMENSION( 0:SIZE(rda_v)+1 ) :: rla_ux, rla_uxad

    REAL(KIND=dp) :: rl_dt, rl_tmp, rl_a, rl_norm, rl_normr, rl_tmpnum, rl_tmpden, rl_quot,&
      rl_aad, rl_normad, rl_normrad, rl_tmpnumad, rl_tmpdenad, rl_quotad
    INTEGER             :: ib_iter, il_ub, il_lb, chkp_unit

    !zeroing local adjoint variables
    rla_tmpad = 0.0_cp
    rla_bad = 0.0_cp
    rla_cad = 0.0_cp
    rla_uxad = 0.0_cp
    rl_aad = 0.0_cp
    rl_normad = 0.0_cp
    rl_normrad = 0.0_cp
    rl_tmpnumad = 0.0_cp
    rl_tmpdenad = 0.0_cp
    rl_quotad = 0.0_cp

    !recomputing local variables
    IF( PRESENT(c_unit) )THEN
      chkp_unit = c_unit
      CALL debug(c_unit, 'In gd_projection_1v_2dadj; using c_unit = ')
    ELSE
      chkp_unit = chkp_unit_1v_1d
      CALL debug(chkp_unit_1v_1d, 'In gd_projection_1v_2dadj; using chkp_unit_1v_1d = ')
    END IF
    il_lb = 1
    il_ub = SIZE(rda_v)
    !To satisfy the condition of Courant-Friedrichs-Lewy, the time step is defined as follows
    rl_tmp = rd_dx**2/4.0_cp
    rl_dt = rl_tmp-MIN(epsilon_gd, rl_tmp/10.0_cp)
    rl_a  = rl_dt / ( rd_dx**2 )
    rla_b = -(rda_phi + 2.0_cp/(rd_dx**2) )*rl_dt + 1.0_cp
    rla_c = rl_dt * (rda_phi*rda_v)

    !Adjoint computations
    CALL chkp_open( chkp_unit )!open for read
    !recomputing rl_normr and rl_norm
    CALL chkp_restaure(chkp_unit, rla_ux(il_lb:il_ub), id_niter+1)
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    rl_norm = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb:il_ub)**2
    rl_normr = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb:il_ub)
    !end of recomputing
    IF ( rl_normr .GT. epsilon_gd ) THEN
      !recomputing rl_tmpnum, rl_tmpden and rl_quot
      rl_tmpnum = SQRT(rl_norm)
      rl_tmpden = SQRT(rl_normr)
      rl_quot = rl_tmpnum/rl_tmpden
      !**End of recomputing
      rla_tmpad = rla_tmpad + rl_quot * rda_uad
      rl_quotad = rl_quotad + SUM( rla_tmp*rda_uad )
      rda_uad = 0.0_cp
      rl_tmpnumad = rl_tmpnumad + (1 / rl_tmpden) * rl_quotad
      rl_tmpdenad = rl_tmpdenad - (rl_tmpnum / rl_tmpden**2) * rl_quotad
      rl_quotad = 0.0_dp
      rl_normrad = rl_normrad + (1/ (2*SQRT(rl_normr) ) )*rl_tmpdenad
      rl_tmpdenad = 0.0_dp
      rl_normad = rl_normad + (1/(2*SQRT(rl_norm)) )*rl_tmpnumad
      rl_tmpnumad = 0.0_dp
    ELSE
      rla_tmpad = rla_tmpad + rda_uad
      rda_uad = 0.0_cp
    END IF
    rla_uxad(il_lb:il_ub) = rla_uxad(il_lb:il_ub) + rla_tmpad
    rla_tmpad = 0.0_cp
    !Recomputing rla_tmp
    rla_tmp = rla_ux(il_lb:il_ub)**2
    !End of recomputing
    CALL MAXVALADJ(rla_tmp, rla_tmpad, rl_normrad)
    rla_uxad(il_lb:il_ub) = rla_uxad(il_lb:il_ub) + 2.0_cp*rla_ux(il_lb:il_ub)*rla_tmpad
    rla_tmpad = 0.0_cp
    !Recomputing rla_tmp
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    !End of recomputing
    CALL MAXVALADJ(rla_tmp, rla_tmpad, rl_normad)
    rda_vad = rda_vad + 2.0_cp*rda_v*rla_tmpad*rda_phi/MAXVAL(rda_phi)
    rla_tmpad = 0.0_cp
    DO ib_iter = id_niter, 1, -1
      !recomputing rla_ux
      CALL chkp_restaure(chkp_unit, rla_ux(il_lb:il_ub), id_niter)
      CALL gd_projection_bc(rla_ux)
      !End of recomputing
      CALL gd_projection_bcadj(rla_uxad)
      rla_tmpad = rla_tmpad + rla_uxad(il_lb:il_ub)
      rla_uxad(il_lb:il_ub) = 0.0_cp
      rla_uxad(il_lb+1:il_ub+1) = rla_uxad(il_lb+1:il_ub+1) + rl_a  * rla_tmpad
      rla_uxad(il_lb-1:il_ub-1) = rla_uxad(il_lb-1:il_ub-1) + rl_a  * rla_tmpad
      rla_uxad(il_lb:il_ub)     = rla_uxad(il_lb:il_ub)     + rla_b * rla_tmpad
      rla_cad                   = rla_cad                   + rla_tmpad
      rla_tmpad = 0.0_cp
    END DO
    CALL chkp_close(chkp_unit)
    CALL gd_projection_bcadj(rla_uxad)
    rda_vad = rda_vad + rla_uxad(il_lb:il_ub)
    rla_uxad(il_lb:il_ub) = 0.0_cp
    rda_vad = rda_vad + rl_dt * rda_phi * rla_cad
    rla_cad = 0.0_cp
    !End of adjoint computations
  END SUBROUTINE gd_projection_1v_1dadj

  !> @brief Generalised diffusion (quasi-)projection of one variable in two dimensional space
  !! @param [in] rda_v variable to be projected
  !! @param [in,out] rda_vad adjoint variable
  !! @param [in,out] rda_uad adjoint variable
  !! @param [in] rda_phi Trust function
  !! @param [in] rd_dx discretization step along the x direction
  !! @param [in] rd_dy discretization step along the y direction
  !! @param [in] id_niter number of iterations for the iterative process
  !! @param [in] c_unit id of the checkpointing unit, this must be
  !!  an integer between 1 and im_nunit
  !! @details
  !! See checkpoint module for details on checkpointing units.
  !! The checkpointing unit is used to save intermediate results
  !! for adjoint calculation. If a configuration of the projection
  !! is used more than once and the adjoint calculation is necessary,
  !! the user should give different value to c_unit for each use of
  !! the configuration. If a configuration is used only once,
  !! the system can manage it. The configuration of the projection
  !! is defined by the dimensionality and the number of input to
  !! be projected. Example of configurations are :
  !!   - one variable in one dimension,
  !!   - two variables in one dimension,
  !!   - one variable in two dimensions,
  !!   - two variables in two dimension.
  !!  @todo check this additional parameter @ rda_uad. It seems to be there by error
  !!
  !<
  SUBROUTINE gd_projection_1v_2dadj(rda_v, rda_vad, rda_uad, rda_phi, rd_dx, rd_dy, id_niter, c_unit)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)   :: rda_v, rda_phi
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)   :: rda_vad
    REAL(KIND=dp), DIMENSION( SIZE(rda_v,1), SIZE(rda_v,2) ), INTENT(INOUT) :: rda_uad!!result of the projection
    REAL(KIND=dp) , INTENT(IN) :: rd_dx, rd_dy!!space discretization step
    INTEGER, INTENT(IN) :: id_niter!!iteration count, amount of diffusion iterations
    INTEGER, OPTIONAL, INTENT(IN) :: c_unit!!id of the checkpointing unit
    !!local variables
    REAL(KIND=dp), DIMENSION( SIZE(rda_v,1), SIZE(rda_v,2) ) :: rla_tmp, rla_b, rla_c
    REAL(KIND=dp), DIMENSION( SIZE(rda_v,1), SIZE(rda_v,2) ) :: rla_tmpad, rla_bad, rla_cad
    !!local extended variables used to handle mirror boundary conditions.  rla_ux, x for extended
    REAL(KIND=dp), DIMENSION( 0:SIZE(rda_v,1)+1, 0:SIZE(rda_v,2)+1 ) :: rla_ux
    REAL(KIND=dp), DIMENSION( 0:SIZE(rda_v,1)+1, 0:SIZE(rda_v,2)+1 ) :: rla_uxad

    REAL(KIND=dp)       :: rl_dt, rl_tmp, rl_a1, rl_a2, rl_norm, rl_normr, rl_tmpnum, rl_tmpden, rl_quot
    REAL(KIND=dp)       :: rl_a1ad, rl_a2ad, rl_normad, rl_normrad, rl_tmpnumad, rl_tmpdenad, rl_quotad
    INTEGER             :: ib_iter, il_ub1, il_ub2, il_lb1, il_lb2, chkp_unit

    CALL debug('', 'Entering gd_projection_1v_2dadj ++++++++++++++++++++++++++++++')
    !zeroing local adjoint variables
    rla_tmpad = 0.0_cp
    rla_bad = 0.0_cp
    rla_cad = 0.0_cp
    rla_uxad = 0.0_cp
    rl_a1ad = 0.0_cp
    rl_a2ad = 0.0_cp
    rl_normad = 0.0_cp
    rl_normrad = 0.0_cp
    rl_tmpnumad = 0.0_cp
    rl_tmpdenad = 0.0_cp
    rl_quotad = 0.0_cp
    !End of zeroing
    !recomputing local variables
    CALL debug(300, 'In gd_projection_1v_2dadj :')
    IF( PRESENT(c_unit) )THEN
      CALL debug(310, 'In gd_projection_1v_2dadj :')
      chkp_unit = c_unit
      CALL debug(c_unit, 'In gd_projection_1v_2dadj; using c_unit = ')
    ELSE
      CALL debug(320, 'In gd_projection_1v_2dadj :')
      chkp_unit = chkp_unit_1v_2d
      CALL debug(chkp_unit_1v_2d, 'In gd_projection_1v_2dadj; chkp_unit_1v_2d = ')
      CALL debug(chkp_unit_1v_2d, 'In gd_projection_1v_2dadj; using chkp_unit_1v_2d = ')
    END IF
    CALL debug(350, 'In gd_projection_1v_2dadj :')
    il_lb1 = 1
    il_lb2 = 1
    il_ub1 = SIZE(rda_v,1)
    il_ub2 = SIZE(rda_v,2)
    rl_tmp = (rd_dx*rd_dy)**2/( 4.0_cp*(rd_dx**2 + rd_dy**2) )
    rl_dt = rl_tmp-MIN(epsilon_gd, rl_tmp/10.0_cp)
    rl_a1  = rl_dt / ( rd_dx**2 )
    rl_a2  = rl_dt / ( rd_dy**2 )
    rla_b = -(rda_phi + 2.0_cp/(rd_dx**2) + 2.0_cp/(rd_dy**2) )*rl_dt + 1.0_cp
    rla_c = rl_dt * (rda_phi*rda_v)
    !End of recomputing
    CALL debug(400, 'In gd_projection_1v_2dadj :')
    !Adjoint computation
    CALL chkp_open( chkp_unit )!open for read
    !recomputing rl_normr and rl_norm
    CALL chkp_restaure(chkp_unit, rla_ux(il_lb1:il_ub1,il_lb2:il_ub2), id_niter+1)
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    rl_norm = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb1:il_ub1,il_lb2:il_ub2)**2
    rl_normr = MAXVAL(rla_tmp)
    rla_tmp = rla_ux(il_lb1:il_ub1,il_lb2:il_ub2)
    CALL debug(500, 'In gd_projection_1v_2dadj :')
    !end of recomputing
    IF ( rl_normr .GT. epsilon_gd ) THEN
      !recomputing rl_tmpnum, rl_tmpden and rl_quot
      rl_tmpnum = SQRT(rl_norm)
      rl_tmpden = SQRT(rl_normr)
      rl_quot = rl_tmpnum/rl_tmpden
      !**End of recomputing
      rla_tmpad = rla_tmpad + rl_quot * rda_uad
      rl_quotad = rl_quotad + SUM( rla_tmp*rda_uad )
      rda_uad = 0.0_cp
      rl_tmpnumad = rl_tmpnumad + (1 / rl_tmpden) * rl_quotad
      rl_tmpdenad = rl_tmpdenad - (rl_tmpnum / rl_tmpden**2) * rl_quotad
      rl_quotad = 0.0_dp
      rl_normrad = rl_normrad + (1/ (2*SQRT(rl_normr) ) )*rl_tmpdenad
      rl_tmpdenad = 0.0_dp
      rl_normad = rl_normad + (1/(2*SQRT(rl_norm)) )*rl_tmpnumad
      rl_tmpnumad = 0.0_dp
    ELSE
      rla_tmpad = rla_tmpad + rda_uad
      rda_uad = 0.0_cp
    END IF
    rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2) = rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2) + rla_tmpad
    rla_tmpad = 0.0_cp
    !Recomputing rla_tmp
    rla_tmp = rla_ux(il_lb1:il_ub1,il_lb2:il_ub2)**2
    !End of recomputing
    CALL MAXVALADJ(rla_tmp, rla_tmpad, rl_normrad)
    rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2) = rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2)&
      + 2.0_cp*rla_ux(il_lb1:il_ub1,il_lb2:il_ub2)*rla_tmpad
    rla_tmpad = 0.0_cp
    !Recomputing rla_tmp
    rla_tmp = rda_v**2*rda_phi/MAXVAL(rda_phi)
    !End of recomputing
    CALL MAXVALADJ(rla_tmp, rla_tmpad, rl_normad)
    rda_vad = rda_vad + 2.0_cp*rda_v*rla_tmpad*rda_phi/MAXVAL(rda_phi)
    rla_tmpad = 0.0_cp
    DO ib_iter = id_niter, 1, -1
      !recomputing rla_ux
      CALL chkp_restaure(chkp_unit, rla_ux(il_lb1:il_ub1,il_lb2:il_ub2), id_niter)
      CALL gd_projection_bc(rla_ux)
      !End of recomputing
      CALL gd_projection_bcadj(rla_uxad)
      rla_tmpad = rla_tmpad + rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2)
      rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2) = 0.0_cp
      rla_uxad(il_lb1+1:il_ub1+1, il_lb2:il_ub2) = rla_uxad(il_lb1+1:il_ub1+1, il_lb2:il_ub2) + rl_a1 * rla_tmpad
      rla_uxad(il_lb1-1:il_ub1-1, il_lb2:il_ub2) = rla_uxad(il_lb1-1:il_ub1-1, il_lb2:il_ub2) + rl_a1 * rla_tmpad
      rla_uxad(il_lb1:il_ub1, il_lb2+1:il_ub2+1) = rla_uxad(il_lb1:il_ub1, il_lb2+1:il_ub2+1) + rl_a2 * rla_tmpad
      rla_uxad(il_lb1:il_ub1, il_lb2-1:il_ub2-1) = rla_uxad(il_lb1:il_ub1, il_lb2-1:il_ub2-1) + rl_a2 * rla_tmpad
      rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2)      = rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2)      + rla_b * rla_tmpad
      rla_cad                                    = rla_cad                                    + rla_tmpad
      rla_tmpad = 0.0_cp
    END DO
    CALL debug(600, 'In gd_projection_1v_2dadj :')
    CALL chkp_close(chkp_unit)
    CALL gd_projection_bcadj(rla_uxad)
    rda_vad = rda_vad + rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2)
    rla_uxad(il_lb1:il_ub1,il_lb2:il_ub2) = 0.0_cp
    rda_vad = rda_vad + rl_dt * rda_phi * rla_cad
    rla_cad = 0.0_cp
    !End of Adjoint computation
    CALL debug('', 'Exiting gd_projection_1v_2dadj ++++++++++++++++++++++++++++++')

  END SUBROUTINE gd_projection_1v_2dadj

!   !> @brief Generalised diffusion (quasi)projection of one variables in two dimensional space
!   SUBROUTINE gd_projection_1v_2d()
!   END SUBROUTINE gd_projection_1v_2d

  !> @brief Ensure boundary conditions
  !! @param [out] rda_Aad array to ensured boundary conditions
  SUBROUTINE gd_projection_bc_1dadj(rda_Aad)
    REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)  :: rda_Aad
    INTEGER                    :: il_dm, il_fm

    !recomputing local variables
    il_dm = LBOUND(rda_Aad,1)
    il_fm = UBOUND(rda_Aad,1)
    !End of recomputing
    SELECT CASE (gd_boundary)
      CASE(LINEAR_EXTRAPOLATION)
        rda_Aad(il_fm-1) = rda_Aad(il_fm-1) + rm_w_bc         *rda_Aad(il_fm)
        rda_Aad(il_fm-2) = rda_Aad(il_fm-2) + (1.0_cp-rm_w_bc)*rda_Aad(il_fm)
        rda_Aad(il_fm-2) = rda_Aad(il_fm-2) - rm_w_bc         *rda_Aad(il_fm)
        rda_Aad(il_fm-3) = rda_Aad(il_fm-3) - (1.0_cp-rm_w_bc)*rda_Aad(il_fm)
        rda_Aad(il_fm) = 0.0_cp
        rda_Aad(il_dm+1) = rda_Aad(il_dm+1) + rm_w_bc         *rda_Aad(il_dm)
        rda_Aad(il_dm+2) = rda_Aad(il_dm+2) + (1.0_cp-rm_w_bc)*rda_Aad(il_dm)
        rda_Aad(il_dm+2) = rda_Aad(il_dm+2) - rm_w_bc         *rda_Aad(il_dm)
        rda_Aad(il_dm+3) = rda_Aad(il_dm+3) - (1.0_cp-rm_w_bc)*rda_Aad(il_dm)
        rda_Aad(il_dm) = 0.0_cp
      CASE(MIRROR_BOUNDARY)
        rda_Aad(il_fm-2) = rda_Aad(il_fm-2) + rda_Aad(il_fm)
        rda_Aad(il_fm) = 0.0_cp
        rda_Aad(il_dm+2) = rda_Aad(il_dm+2) + rda_Aad(il_dm)
        rda_Aad(il_dm) = 0.0_cp
    END SELECT
  END SUBROUTINE gd_projection_bc_1dadj

  SUBROUTINE gd_projection_bc_2dadj(rda_aad)
    REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT) :: rda_aad
    INTEGER                                :: il_dm, il_fm, il_dn, il_fn

    il_dm = LBOUND(rda_aad,1)
    il_fm = UBOUND(rda_aad,1)
    il_dn = LBOUND(rda_aad,2)
    il_fn = UBOUND(rda_aad,2)

    SELECT CASE (gd_boundary)
      !the direction of the substraction operation is adjusted to have uniform sign at every boundaries
      CASE(LINEAR_EXTRAPOLATION)
        rda_Aad(:, il_fn-1) = rda_Aad(:, il_fn-1) + rm_w_bc         *rda_Aad(:, il_fn)
        rda_Aad(:, il_fn-2) = rda_Aad(:, il_fn-2) + (1.0_cp-rm_w_bc)*rda_Aad(:, il_fn)
        rda_Aad(:, il_fn-2) = rda_Aad(:, il_fn-2) - rm_w_bc         *rda_Aad(:, il_fn)
        rda_Aad(:, il_fn-3) = rda_Aad(:, il_fn-3) - (1.0_cp-rm_w_bc)*rda_Aad(:, il_fn)
        rda_Aad(:, il_fn  ) = 0.0_cp
        rda_Aad(:, il_dn+1) = rda_Aad(:, il_dn+1) + rm_w_bc         *rda_Aad(:, il_dn)
        rda_Aad(:, il_dn+2) = rda_Aad(:, il_dn+2) + (1.0_cp-rm_w_bc)*rda_Aad(:, il_dn)
        rda_Aad(:, il_dn+2) = rda_Aad(:, il_dn+2) - rm_w_bc         *rda_Aad(:, il_dn)
        rda_Aad(:, il_dn+3) = rda_Aad(:, il_dn+3) - (1.0_cp-rm_w_bc)*rda_Aad(:, il_dn)
        rda_Aad(:, il_dn  ) = 0.0_cp

        rda_Aad(il_fm-1, :) = rda_Aad(il_fm-1, :) + rm_w_bc         *rda_Aad(il_fm, :)
        rda_Aad(il_fm-2, :) = rda_Aad(il_fm-2, :) + (1.0_cp-rm_w_bc)*rda_Aad(il_fm, :)
        rda_Aad(il_fm-2, :) = rda_Aad(il_fm-2, :) - rm_w_bc         *rda_Aad(il_fm, :)
        rda_Aad(il_fm-3, :) = rda_Aad(il_fm-3, :) - (1.0_cp-rm_w_bc)*rda_Aad(il_fm, :)
        rda_Aad(il_fm, :) = 0.0_cp
        rda_Aad(il_dm+1, :) = rda_Aad(il_dm+1, :) + rm_w_bc         *rda_Aad(il_dm, :)
        rda_Aad(il_dm+2, :) = rda_Aad(il_dm+2, :) + (1.0_cp-rm_w_bc)*rda_Aad(il_dm, :)
        rda_Aad(il_dm+2, :) = rda_Aad(il_dm+2, :) - rm_w_bc         *rda_Aad(il_dm, :)
        rda_Aad(il_dm+3, :) = rda_Aad(il_dm+3, :) - (1.0_cp-rm_w_bc)*rda_Aad(il_dm, :)
        rda_Aad(il_dm, :) = 0.0_cp
      CASE(MIRROR_BOUNDARY)
        rda_aad(:, il_fn-2) = rda_aad(:, il_fn-2) + rda_aad(:, il_fn)
        rda_aad(:, il_fn  ) = 0.0_dp
        rda_aad(:, il_dn+2) = rda_aad(:, il_dn+2) + rda_aad(:, il_dn)
        rda_aad(:, il_dn  ) = 0.0_dp
        rda_aad(il_fm-2,: ) = rda_aad(il_fm-2, :) + rda_aad(il_fm, :)
        rda_aad(il_fm,  : ) = 0.0_dp
        rda_aad(il_dm+2,: ) = rda_aad(il_dm+2, :) + rda_aad(il_dm, :)
        rda_aad(il_dm,  : ) = 0.0_dp
    END SELECT

  END SUBROUTINE gd_projection_bc_2dadj

END MODULE gd_toolsadj