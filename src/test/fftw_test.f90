PROGRAM test_fftw
  USE general_tools
  USE cs_tools
  USE debug_tools
  USE fftw_tools
IMPLICIT NONE

  INTEGER, PARAMETER :: IP_SIZE = 512
  INTEGER :: ibi
  REAL(cp) :: rp_scale = 3.0
  REAL(cp) :: rl_nz_ratio
  real(C_DOUBLE) :: rl_dx
  type(C_PTR) :: plan_forward, plan_backward
  real(C_DOUBLE), dimension(IP_SIZE) :: rla_x, rla_y, rla_in, rla_out, rla_signal, rla_fft_coef
  TYPE(dyn_rVector), DIMENSION(5) :: tla_save

  CALL init_all()
  
  rl_dx = 2.0*rp_pi/(IP_SIZE-1)
  rla_x = (/ ((ibi-1)*rl_dx, ibi=1,IP_SIZE) /)
  tla_save(1) = rla_x
  rla_fft_coef = 0.0_cp
  rla_signal   = 0.0_cp
  rl_nz_ratio  = 0.1
  CALL RANDOM_NUMBER(rla_fft_coef)
  rla_fft_coef = rp_scale*(rla_fft_coef - 0.5_cp)
  CALL randval_randidx(rla_fft_coef, rl_nz_ratio, rp_scale)
  tla_save(3) = rla_fft_coef
  CALL ifftw_trans(rla_fft_coef, rla_signal)
  tla_save(2) = rla_signal
  CALL RANDOM_NUMBER(rla_signal)
  rla_signal = rp_scale*(rla_signal - 0.5_cp)
  CALL randval_randidx(rla_signal, rl_nz_ratio, rp_scale)
  tla_save(4) = rla_signal
  CALL fftw_trans(rla_signal, rla_fft_coef)
  tla_save(5) = rla_fft_coef
  CALL save_trj( tla_save(1:5) , TRIM('output/test.dat' ), ld_column = .TRUE. )
  
  CALL finalize_all()
CONTAINS

  SUBROUTINE init_all()
    CALL init_fftw()
    CALL init_cs_tools()
  END SUBROUTINE init_all
  
  SUBROUTINE finalize_all()
    CALL finalize_fftw()
    CALL finalize_cs_tools()
  END SUBROUTINE finalize_all
  
  SUBROUTINE test()
    !complex(C_DOUBLE_COMPLEX), dimension(IP_SIZE/2 + 1) :: cla_out
    !fftw_plan fftw_plan_r2r_1d(int n, double *in, double *out, fftw_r2r_kind kind, unsigned flags); FFTW_R2HC/FFTW_HC2R
    rl_dx = 2.0*rp_pi/(IP_SIZE-1)
    DO ibi=1, IP_SIZE
      rla_x(ibi) = (ibi-1)*rl_dx!-rp_pi
    END DO
    rla_y = cos(22*rla_x) -cos(40*rla_x)! - sin(16*rla_x)
    tla_save(1) = rla_x
    tla_save(2) = rla_y
    plan_forward  = fftw_plan_r2r_1d(IP_SIZE, rla_in , rla_out, FFTW_R2HC, FFTW_ESTIMATE)
    plan_backward = fftw_plan_r2r_1d(IP_SIZE, rla_out, rla_in , FFTW_HC2R, FFTW_ESTIMATE)
    !...
    rla_in = rla_y
    call fftw_execute_r2r(plan_forward , rla_in, rla_out)
    rla_out = rla_out/SQRT(IP_SIZE*1.0_cp)
    tla_save(3) = rla_out
    rla_in = 0.0_cp
    WHERE (ABS(rla_out)<MAXVAL(rla_out)/10.0_cp) rla_out=0.0_cp
    tla_save(4) = rla_out
    call fftw_execute_r2r(plan_backward, rla_out, rla_in)
    rla_in = rla_in/SQRT(IP_SIZE*1.0_cp)
    tla_save(5) = rla_in
    !...
    !...
    call fftw_destroy_plan(plan_forward )
    call fftw_destroy_plan(plan_backward)
    CALL save_trj( tla_save , TRIM('output/fft.dat' ), ld_column = .TRUE. )
  END SUBROUTINE test
END PROGRAM test_fftw