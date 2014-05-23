PROGRAM gd_test
  USE randmod
  USE general_tools
  USE debug_tools
  USE checkpoint
  USE gd_tools
IMPLICIT NONE
  INTEGER, PARAMETER :: &
    ip_xsize=512,&
    ip_ysize=510
  REAL(KIND=cp), PARAMETER ::&
    rp_xmin = -2*rp_pi,&
    rp_xmax = +2*rp_pi,&
    rp_ymin = -2*rp_pi,&
    rp_ymax = +2*rp_pi

  !CALL test_1D()
  CALL test_2D()

CONTAINS

  SUBROUTINE test_1D()
    INTEGER :: ibi
    REAL(KIND=cp) :: rl_dx, rl_rand
    REAL(KIND=cp), DIMENSION(ip_xsize) :: rla_x, rla_phi, rla_f, rla_f_noise, rla_f_gd, rla_res
    TYPE(dyn_rVector), DIMENSION(4) :: tla_save

    rl_dx = (rp_xmax - rp_xmin)/REAL(ip_xsize-1)
    DO ibi = 1, ip_xsize
      rla_x(ibi) = rp_xmin+REAL(ibi-1, cp)*rl_dx
    END DO
    rla_f = 1.0*sin(rla_x)
    rla_phi = 1.0
    !initializing random numbers generator
    CALL init_normal_rand(mu=0.0_cp, sigma=0.5_cp, radius=1.5_cp, nright=256, stat=.FALSE.)
    !Adding noise
    DO ibi = 1, ip_xsize
      CALL normal_rand(rl_rand)
      rla_f_noise(ibi) = rla_f(ibi) + rl_rand
      rla_phi(ibi) = 1.0 - MIN( 1.0, ABS(rl_rand) )
    END DO

    CALL gd_projection(rla_f_noise, rla_f_gd, rla_phi, rl_dx, 50)
    rla_res = rla_f_gd-rla_f
    tla_save(1) = rla_x
    tla_save(2) = rla_f
    tla_save(3) = rla_f_noise
    tla_save(4) = rla_f_gd
    CALL debug('', 'Original noise *****************************')
    CALL debug(MAXVAL( ABS(rla_f_noise) ),      '  MAXVAL    = ')
    CALL debug(SUM( ABS(rla_f_noise) )/ip_xsize, '  MEAN      = ')
    CALL debug(SQRT( SUM(rla_f_noise**2) ),     '  L2 norm   = ')
    CALL debug('', '********************************************')
    CALL debug('', 'Residual from GD ***************************')
    CALL debug(MAXVAL( ABS(rla_res) ),      '  MAXVAL    = ')
    CALL debug(SUM( ABS(rla_res) )/ip_xsize, '  MEAN      = ')
    CALL debug(SQRT( SUM(rla_res**2) ),     '  L2 norm   = ')
    CALL debug('', '********************************************')

    CALL save_trj( tla_save(2:4), 'runtime_gd_test.dat', ld_column = .TRUE. )
    CALL reset_drv_array( tla_save )
    CALL finilize_normal_rand()
    CALL debug('', 'GD_test 1D end')
  END SUBROUTINE test_1D

  SUBROUTINE test_2D()
    INTEGER :: ibi, ibj
    REAL(KIND=cp) :: rl_dx, rl_dy, rl_rand
    REAL(KIND=cp), DIMENSION(ip_xsize, ip_ysize) :: rla_x, rla_y, rla_phi, rla_f, rla_f_noise, rla_f_gd, rla_res
    TYPE(dyn_rVector), DIMENSION(4) :: tla_save

    rl_dx = (rp_xmax - rp_xmin)/REAL(ip_xsize-1)
    rl_dy = (rp_ymax - rp_ymin)/REAL(ip_ysize-1)
    DO ibi = 1, ip_xsize
      rla_x(ibi, :) = rp_xmin+REAL(ibi-1, cp)*rl_dx
    END DO
    DO ibj = 1, ip_ysize
      rla_y(:, ibj) = rp_ymin+REAL(ibj-1, cp)*rl_dy
    END DO
    rla_f = 1.0*sin(SQRT(rla_x**2 + rla_y**2))
    rla_phi = 1.0
    !initializing random numbers generator
    CALL init_normal_rand(mu=0.0_cp, sigma=0.5_cp, radius=1.5_cp, nright=256, stat=.FALSE.)
    !Adding noise
    DO ibj = 1, ip_ysize
      DO ibi = 1, ip_xsize
        CALL normal_rand(rl_rand)
        rla_f_noise(ibi, ibj) = rla_f(ibi, ibj) + rl_rand
        rla_phi(ibi, ibj) = 1.0 - MIN( 1.0, ABS(rl_rand) )
      END DO
    END DO

    CALL debug('', 'Calling gd_projection')
    CALL gd_projection(rla_f_noise, rla_f_gd, rla_phi, rl_dx, rl_dy, 20)
    rla_res = rla_f_gd-rla_f
    CALL debug('', 'After gd_projection')
    CALL debug('', 'Original noise *****************************')
    CALL debug(MAXVAL( ABS(rla_f_noise) ),      '  MAXVAL    = ')
    CALL debug(SUM( ABS(rla_f_noise) )/(ip_xsize*ip_xsize), '  MEAN      = ')
    CALL debug(SQRT( SUM(rla_f_noise**2) ),     '  L2 norm   = ')
    CALL debug('', '********************************************')
    CALL debug('', 'Residual from GD ***************************')
    CALL debug(MAXVAL( ABS(rla_res) ),      '  MAXVAL    = ')
    CALL debug(SUM( ABS(rla_res) )/(ip_xsize*ip_xsize), '  MEAN      = ')
    CALL debug(SQRT( SUM(rla_res**2) ),     '  L2 norm   = ')
    CALL debug('', '********************************************')

!     tla_save(1) = rla_x
!     tla_save(2) = rla_f
!     tla_save(3) = rla_f_noise
!     tla_save(4) = rla_f_gd
!     CALL save_trj( tla_save(2:4), 'runtime_gd_test.dat', ld_column = .TRUE. )
!     CALL reset_drv_array( tla_save )
    CALL finilize_normal_rand()
    CALL debug('', 'GD_test 2D end')
  END SUBROUTINE test_2D

END PROGRAM gd_test