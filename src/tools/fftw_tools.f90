!> \file fftw_tools.f90
!! @brief general tools for Fourier transform
!! @author Innocent Souopgui
!! \details
!!
!<
MODULE fftw_tools
  USE, intrinsic :: iso_c_binding
  USE general_constant
  USE debug_tools
IMPLICIT NONE
  include 'fftw3.f03'
  INTEGER, PARAMETER :: FFTW_REAL = C_DOUBLE
  INTEGER, PARAMETER :: FFTW_CPLX = C_DOUBLE_COMPLEX
  INTERFACE init_fftw_engine
    MODULE PROCEDURE init_fftw_1d_engine!> init 1D engine
    !MODULE PROCEDURE init_fftw_2d_engine!> init 2D engine
    !MODULE PROCEDURE init_fftw_3d_engine!> init 3D engine
  END INTERFACE

  INTERFACE fftw_trans
    MODULE PROCEDURE fftw_1d_r2r!> run 1D direct transform
    !MODULE PROCEDURE fftw_2d_r2r!> run 2D direct transform
    !MODULE PROCEDURE fftw_3d_r2r!> run 3D direct transform
  END INTERFACE

  INTERFACE ifftw_trans
    MODULE PROCEDURE ifftw_1d_r2r!> run 1D inverse transform
    !MODULE PROCEDURE ifftw_2d_r2r!> run 2D inverse transform
    !MODULE PROCEDURE ifftw_3d_r2r!> run 3D inverse transform
  END INTERFACE

  !>@brief engine for 1D fftw transform
  TYPE, PRIVATE :: fftw_1d_engine
    INTEGER :: i_size = -1 !> size of the transform supported by the engine
    INTEGER(C_INT) :: i_kind
    INTEGER(C_INT) :: i_flags
    TYPE(C_PTR) :: t_plan !> fftw plan
    REAL(FFTW_REAL), DIMENSION(:), POINTER :: ra_in =>NULL()!>  real input vector for fftw transform
    REAL(FFTW_REAL), DIMENSION(:), POINTER :: ra_out=>NULL()!> real output vector for fftw transform
    !> complex input vector for fftw transform
    !<
    COMPLEX(FFTW_CPLX), DIMENSION(:), POINTER :: ca_in =>NULL()
    !> complex output vector for fftw transform
    !<
    COMPLEX(FFTW_CPLX), DIMENSION(:), POINTER :: ca_out=>NULL()
  END TYPE fftw_1d_engine

  TYPE(fftw_1d_engine), PRIVATE, SAVE :: &
    tm_r2r_1d_feng,& !> forward engine for r2r transform (R2HC)
    tm_r2r_1d_beng   !> backward engine for r2r transform (HC2R)
CONTAINS

  !>@brief Initializes environment for fftw transform
  !!@details This routine does nothing. Only the associated routine finalize_fftw is usefull. For conformance, it should be used so that the use of finalize make logical sense
  !<
  SUBROUTINE init_fftw()
  END SUBROUTINE init_fftw

  !>@brief Terminates the fftw environment
  !!@details This routine clean the environment: destroy plans and free dynamically allocated spaces for fftw transform
  !<
  SUBROUTINE finalize_fftw()
    !reset engines for 1D fftw transform
    CALL reset_fftw_1d_engine(tm_r2r_1d_feng)
    CALL reset_fftw_1d_engine(tm_r2r_1d_beng)
  END SUBROUTINE finalize_fftw

  !>@brief initialises engine for 1D fftw transform
  !!@param [in,out] td_eng Engine for fftw transform
  !!@param [in] id_size size of the transform
  !!@param [in] id_kind king of the transform
  !!@param [in] id_flags flags for fftw plan
  !<
  SUBROUTINE init_fftw_1d_engine(td_eng, id_size, id_kind, id_flags)
    TYPE(fftw_1d_engine), INTENT(INOUT) :: td_eng
    INTEGER, INTENT(IN) :: id_size
    INTEGER(C_INT), INTENT(IN) :: id_kind
    INTEGER(C_INT), INTENT(IN), OPTIONAL :: id_flags
    !local variables
    INTEGER(C_INT)::il_flags

    IF(PRESENT(id_flags) )THEN
      il_flags = id_flags
    ELSE
      il_flags = FFTW_ESTIMATE
    END IF

    CALL reset_fftw_1d_engine(td_eng)!FFTW_HC2R
    !memory allocation
    SELECT CASE(id_kind)
      CASE (FFTW_R2HC, FFTW_HC2R)
        ALLOCATE( td_eng%ra_in(id_size), td_eng%ra_out(id_size) )
      CASE DEFAULT
        CALL stop_progam(id_kind, 'In init_fftw_1d_engine; unknown fftw kind: ')
    END SELECT
    !creating plan
    SELECT CASE(id_kind)
      CASE (FFTW_R2HC)
        td_eng%t_plan  = fftw_plan_r2r_1d(id_size, td_eng%ra_in , td_eng%ra_out, id_kind, il_flags)
      CASE (FFTW_HC2R)
        td_eng%t_plan  = fftw_plan_r2r_1d(id_size, td_eng%ra_in , td_eng%ra_out, id_kind, il_flags)
      CASE DEFAULT
        CALL stop_progam(id_kind, 'In init_fftw_1d_engine; unknown fftw kind: ')
    END SELECT
    td_eng%i_size = id_size
    td_eng%i_kind = id_kind
    td_eng%i_flags= il_flags
  END SUBROUTINE init_fftw_1d_engine

  !>@brief reset a 1D fftw transform engine
  SUBROUTINE reset_fftw_1d_engine(td_eng)
    TYPE(fftw_1d_engine), INTENT(INOUT) :: td_eng
    IF(td_eng%i_size > 0)THEN
      SELECT CASE(td_eng%i_kind)
        CASE (FFTW_R2HC, FFTW_HC2R)
          CALL fftw_destroy_plan(td_eng%t_plan)
          DEALLOCATE( td_eng%ra_in, td_eng%ra_out )
        CASE DEFAULT
          CALL stop_progam(td_eng%i_kind, 'In reset_fftw_1d_engine; unknown fftw kind: ')
      END SELECT
      td_eng%i_size = -1
      td_eng%i_kind = -999
      td_eng%i_flags= -999
    END IF
  END SUBROUTINE reset_fftw_1d_engine

  !>@brief fftw direct r2r 1D transform
  !! @param [in] rda_signal signal to be transform
  !! @param [out] rda_fft_coef Fourier coefficient of the signal
  !<
  SUBROUTINE fftw_1d_r2r(rda_signal, rda_fft_coef)
    REAL(FFTW_REAL), DIMENSION(:), INTENT(IN)  :: rda_signal
    REAL(FFTW_REAL), DIMENSION(:), INTENT(OUT) :: rda_fft_coef
    !local variables
    INTEGER :: il_size

    il_size = SIZE(rda_signal)
    IF(tm_r2r_1d_feng%i_size /= il_size) THEN
      CALL init_fftw_1d_engine(tm_r2r_1d_feng, il_size, FFTW_R2HC)
    END IF
    tm_r2r_1d_feng%ra_in = rda_signal
    CALL fftw_execute_r2r(tm_r2r_1d_feng%t_plan, tm_r2r_1d_feng%ra_in, tm_r2r_1d_feng%ra_out)
    rda_fft_coef = tm_r2r_1d_feng%ra_out/SQRT(tm_r2r_1d_feng%i_size*1.0_cp)
  END SUBROUTINE fftw_1d_r2r

  !>@brief fftw inverse r2r 1D transform
  !! @param [in] rda_fft_coef Fourier coefficient of the signal
  !! @param [out] rda_signal computed signal
  !<
  SUBROUTINE ifftw_1d_r2r(rda_fft_coef, rda_signal)
    REAL(FFTW_REAL), DIMENSION(:), INTENT(IN) :: rda_fft_coef
    REAL(FFTW_REAL), DIMENSION(:), INTENT(OUT):: rda_signal
    !local variables
    INTEGER :: il_size

    il_size = SIZE(rda_signal)
    IF(tm_r2r_1d_beng%i_size /= il_size) THEN
      CALL init_fftw_1d_engine(tm_r2r_1d_beng, il_size, FFTW_HC2R)
    END IF
    tm_r2r_1d_beng%ra_in = rda_fft_coef
    CALL fftw_execute_r2r(tm_r2r_1d_beng%t_plan, tm_r2r_1d_beng%ra_in, tm_r2r_1d_beng%ra_out)
    rda_signal = tm_r2r_1d_beng%ra_out/SQRT(tm_r2r_1d_beng%i_size*1.0_cp)
  END SUBROUTINE ifftw_1d_r2r

END MODULE fftw_tools